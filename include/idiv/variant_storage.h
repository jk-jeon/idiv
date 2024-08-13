// Copyright 2024 Junekey Jeon
//
// The contents of this file may be used under the terms of
// the Apache License v2.0 with LLVM Exceptions.
//
//    (See accompanying file LICENSE-Apache or copy at
//     https://llvm.org/foundation/relicensing/LICENSE.txt)
//
// Alternatively, the contents of this file may be used under the terms of
// the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE-Boost or copy at
//     https://www.boost.org/LICENSE_1_0.txt)
//
// Unless required by applicable law or agreed to in writing, this software
// is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.

#ifndef JKJ_HEADER_IDIV_VARIANT_STORAGE
#define JKJ_HEADER_IDIV_VARIANT_STORAGE

#include "tmp.h"

namespace jkj {
    namespace util {
        // Storage for std::variant-like classes.
        // Follows double-storage strategy to guarantee strong exception safety.
        // The reasons for reinventing this are:
        //   - std::variant used to be no constexpr in C++20; fixed by the defect report P2231R1,
        //     but the patch is not widely available, and
        //   - I want to guarantee strong exception safety.
        // For reference-like proxy types, we follow the assign-through semantics because it is
        // easier to implement consistently without pessimizing the case of value types. Pure
        // reference types are disallowed for the ease of implementation.
        namespace detail {
            // Helpers for variant_storage.
            // Copied and modified from
            // https://github.com/microsoft/STL/blob/18c09c48f5666e6b1ea2a3724c5f6f9917c4c6fb/stl/inc/variant#L865.

            template <std::size_t idx, class TargetType>
            auto construct_array(TargetType (&&)[1])
                -> tmp::typelist<std::integral_constant<std::size_t, idx>, TargetType>;

            template <std::size_t idx, class TargetType, class InitializerType>
            using variant_type_resolver =
                decltype(construct_array<idx, TargetType>({std::declval<InitializerType>()}));

            template <std::size_t idx, class TargetType>
            struct variant_init_single_overload {
                template <class InitializerType>
                auto operator()(TargetType, InitializerType&&)
                    -> variant_type_resolver<idx, TargetType, InitializerType>;
            };

            template <class Indices, class... Types>
            struct variant_init_overload_set_impl;

            template <std::size_t... Indices, class... Types>
            struct variant_init_overload_set_impl<std::index_sequence<Indices...>, Types...>
                : variant_init_single_overload<Indices, Types>... {
                using variant_init_single_overload<Indices, Types>::operator()...;
            };

            template <class... Types>
            using variant_init_overload_set =
                variant_init_overload_set_impl<std::index_sequence_for<Types...>, Types...>;

            // Failure case (has no member "type").
            template <class Enable, class T, class... Types>
            struct variant_init_helper {};

            template <class T, class... Types>
            struct variant_init_helper<std::void_t<decltype(variant_init_overload_set<Types...>{}(
                                           std::declval<T>(), std::declval<T>()))>,
                                       T, Types...> {
                // Perform overload resolution to determine the unique alternative that should be
                // initialized in variant<Types...> from an argument expression with type and value
                // category T.
                using type = decltype(variant_init_overload_set<Types...>{}(std::declval<T>(),
                                                                            std::declval<T>()));
            };
        }

        // Extract the type from variant_init_helper.
        template <class T, class... Types>
        using variant_init_type =
            tmp::get_type<1, typename detail::variant_init_helper<void, T, Types...>::type>;

        // Extract the index from variant_init_helper.
        template <class T, class... Types>
        static constexpr std::size_t variant_init_index =
            tmp::get_type<0, typename detail::variant_init_helper<void, T, Types...>::type>::value;

        template <class... Types>
        class variant_storage {
            using alternative_list = tmp::typelist<Types...>;

            template <std::size_t... I>
            static constexpr std::size_t max() noexcept {
                std::size_t ret_value = 0;
                return ((ret_value = (I > ret_value ? I : ret_value)), ...);
            }

            static_assert(
                (std::is_nothrow_destructible_v<Types> && ...),
                "types with throwing destructor cannot be used as alternatives for a variant");
            static_assert((!std::is_reference_v<Types> && ...),
                          "reference types cannot be used as alternatives for a variant");

            static constexpr std::size_t alignment = max<alignof(Types)...>();
            static constexpr std::size_t size = max<sizeof(Types)...>();

            template <bool = (std::is_nothrow_move_constructible_v<Types> && ...),
                      class dummy = void>
            struct storage_wrapper;

            template <class dummy>
            struct storage_wrapper<true, dummy> {
                static constexpr bool is_using_double_storage = false;
                alignas(alignment) array<std::byte, size> storage;

                constexpr std::byte* data() noexcept { return storage.data(); }
                constexpr std::byte const* data() const noexcept { return storage.data(); }
            };

            template <class dummy>
            struct storage_wrapper<false, dummy> {
                static constexpr bool is_using_double_storage = true;
                alignas(alignment) array<std::byte, size> storage1;
                alignas(alignment) array<std::byte, size> storage2;

                std::byte* storage_pointer_ = storage1.data();
                constexpr std::byte* data() noexcept { return storage_pointer_; }
                constexpr std::byte const* data() const noexcept { return storage_pointer_; }
            };

            storage_wrapper<> wrapped_;

            template <auto discriminator_container = make_index_array<sizeof...(Types)>(),
                      std::size_t discriminator_idx = 0, class Discriminator, class Functor>
            static constexpr decltype(auto)
            visit_helper_helper(Discriminator current_alternative, Functor&& f) noexcept(
                (noexcept(static_cast<Functor&&>(f).template operator()<Types>()) && ...)) {
                if constexpr (discriminator_idx == discriminator_container.size() - 1) {
                    constexpr_assert<error_msgs::index_out_of_range>(
                        current_alternative == discriminator_container[discriminator_idx]);
                    return static_cast<Functor&&>(f)
                        .template operator()<tmp::get_type<discriminator_idx, alternative_list>>();
                }
                else {
                    if (current_alternative == discriminator_container[discriminator_idx]) {
                        return static_cast<Functor&&>(f).template
                        operator()<tmp::get_type<discriminator_idx, alternative_list>>();
                    }
                    else {
                        return visit_helper_helper<discriminator_container, discriminator_idx + 1>(
                            current_alternative, static_cast<Functor&&>(f));
                    }
                }
            }

            template <auto discriminator_container = make_index_array<sizeof...(Types)>(),
                      class Self, class Discriminator, class Functor>
            static constexpr decltype(auto)
            visit_helper(Self&& self, Discriminator current_alternative, Functor&& f) noexcept(
                (noexcept(static_cast<Functor&&>(f)(std::declval<Types>())) && ...)) {
                return visit_helper_helper<discriminator_container>(
                    current_alternative,
                    [&self, &f]<class T>() noexcept(noexcept(static_cast<Functor>(f)(
                        static_cast<Self&&>(self).template get<T>()))) -> decltype(auto) {
                        return static_cast<Functor>(f)(static_cast<Self&&>(self).template get<T>());
                    });
            }

        public:
            template <class T>
            static constexpr bool is_alternative() noexcept {
                return (std::is_same_v<T, Types> || ...);
            }

            static constexpr bool is_nothrow_destructible() noexcept {
                return (std::is_nothrow_destructible_v<Types> && ...);
            }

            variant_storage() = default;
            variant_storage(variant_storage const&) = delete;
            variant_storage(variant_storage&&) = delete;
            variant_storage& operator=(variant_storage const&) = delete;
            variant_storage& operator=(variant_storage&&) = delete;

            // Cannot do T{args...}, only does T(args...).
            template <class T, class... Args>
                requires(is_alternative<T>() && std::is_constructible_v<T, Args...>)
            constexpr T&
            construct(Args&&... args) noexcept(std::is_nothrow_constructible_v<T, Args...>) {
                return *std::construct_at<T>(reinterpret_cast<T*>(wrapped_.data()),
                                             std::forward<Args>(args)...);
            }

            template <class T>
                requires(is_alternative<T>())
            constexpr void destroy() noexcept {
                std::destroy_at<T>(reinterpret_cast<T*>(wrapped_.data()));
            }

            template <class T>
                requires(is_alternative<T>())
            constexpr T& get() & noexcept {
                return *std::launder<T>(reinterpret_cast<T*>(wrapped_.data()));
            }
            template <class T>
                requires(is_alternative<T>())
            constexpr T const& get() const& noexcept {
                return *std::launder<T const>(reinterpret_cast<T const*>(wrapped_.data()));
            }
            template <class T>
                requires(is_alternative<T>())
            constexpr T&& get() && noexcept {
                return static_cast<T&&>(*std::launder<T>(reinterpret_cast<T*>(wrapped_.data())));
            }
            template <class T>
                requires(is_alternative<T>())
            constexpr T const&& get() const&& noexcept {
                return static_cast<T const&&>(
                    *std::launder<T const>(reinterpret_cast<T const*>(wrapped_.data())));
            }

            template <auto discriminator_container = make_index_array<sizeof...(Types)>(),
                      class Discriminator, class Functor>
            decltype(auto) visit(Discriminator current_alternative, Functor&& f) & noexcept(
                noexcept(visit_helper<discriminator_container>(*this, current_alternative,
                                                               static_cast<Functor&&>(f)))) {
                return visit_helper<discriminator_container>(*this, current_alternative,
                                                             static_cast<Functor&&>(f));
            }
            template <auto discriminator_container = make_index_array<sizeof...(Types)>(),
                      class Discriminator, class Functor>
            decltype(auto) visit(Discriminator current_alternative, Functor&& f) const& noexcept(
                noexcept(visit_helper<discriminator_container>(*this, current_alternative,
                                                               static_cast<Functor&&>(f)))) {
                return visit_helper<discriminator_container>(*this, current_alternative,
                                                             static_cast<Functor&&>(f));
            }
            template <auto discriminator_container = make_index_array<sizeof...(Types)>(),
                      class Discriminator, class Functor>
            decltype(auto)
            visit(Discriminator current_alternative, Functor&& f) && noexcept(noexcept(
                visit_helper<discriminator_container>(static_cast<variant_storage&&>(*this),
                                                      current_alternative,
                                                      static_cast<Functor&&>(f)))) {
                return visit_helper<discriminator_container>(*this, current_alternative,
                                                             static_cast<Functor&&>(f));
            }
            template <auto discriminator_container = make_index_array<sizeof...(Types)>(),
                      class Discriminator, class Functor>
            decltype(auto)
            visit(Discriminator current_alternative, Functor&& f) const&& noexcept(noexcept(
                visit_helper<discriminator_container>(static_cast<variant_storage const&&>(*this),
                                                      current_alternative,
                                                      static_cast<Functor&&>(f)))) {
                return visit_helper<discriminator_container>(*this, current_alternative,
                                                             static_cast<Functor&&>(f));
            }

            template <class T, auto discriminator_container = make_index_array<sizeof...(Types)>(),
                      class Discriminator, class... Args>
                requires(is_alternative<T>())
            T& emplace(Discriminator current_alternative,
                       Args&&... args) noexcept(std::is_nothrow_constructible_v<T, Args...>) {
                // No-throw case is simple.
                if constexpr (std::is_nothrow_constructible_v<T, Args...>) {
                    visit_helper_helper<discriminator_container>(
                        current_alternative,
                        [this]<class Current>() noexcept { destroy<Current>(); });
                    return construct<T>(static_cast<Args&&>(args)...);
                }
                // Throwing case.
                else {
                    if constexpr (storage_wrapper<>::is_using_double_storage) {
                        // Construct a new object on the alternative storage.
                        std::byte* alt_storage_ptr =
                            wrapped_.storage_pointer_ == wrapped_.storage1.data()
                                ? wrapped_.storage2.data()
                                : wrapped_.storage1.data();
                        auto& ret_value = *std::construct_at<T>(
                            reinterpret_cast<T*>(alt_storage_ptr), static_cast<Args&&>(args)...);

                        // Destroy the current object.
                        visit_helper_helper<discriminator_container>(
                            current_alternative,
                            [this]<class Current>() noexcept { destroy<Current>(); });

                        // Swap the storage pointer.
                        using std::swap;
                        swap(wrapped_.storage_pointer_, alt_storage_ptr);

                        return ret_value;
                    }
                    else {
                        // Construct a temporary object.
                        auto temp = T(static_cast<Args&&>(args)...);

                        // Destroy the current object.
                        visit_helper_helper<discriminator_container>(
                            current_alternative,
                            [this]<class Current>() noexcept { destroy<Current>(); });

                        // Move the temporary object to the storage.
                        return construct<T>(std::move(temp));
                    }
                }
            }

            template <auto discriminator_container = make_index_array<sizeof...(Types)>(),
                      class Discriminator, class AssignFrom>
                requires(
                    std::is_constructible_v<variant_init_type<AssignFrom, Types...>, AssignFrom> &&
                    std::is_assignable_v<variant_init_type<AssignFrom, Types...>&, AssignFrom>)
            void assign_from(Discriminator current_alternative, AssignFrom&& other) noexcept(
                std::is_nothrow_assignable_v<variant_init_type<AssignFrom, Types...>, AssignFrom> &&
                std::is_nothrow_constructible_v<variant_init_type<AssignFrom, Types...>,
                                                AssignFrom>) {
                using init_type = variant_init_type<AssignFrom, Types...>;
                constexpr auto init_index = variant_init_index<AssignFrom, Types...>;
                if (current_alternative == discriminator_container[init_index]) {
                    get<init_type>() = static_cast<AssignFrom&&>(other);
                }
                else {
                    emplace<init_type, discriminator_container>(current_alternative,
                                                                static_cast<AssignFrom&&>(other));
                }
            }
        };
    }
}

#endif
