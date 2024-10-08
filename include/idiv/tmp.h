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

#ifndef JKJ_HEADER_IDIV_TMP
#define JKJ_HEADER_IDIV_TMP

#include "util.h"

namespace jkj {
    // Some metaprogramming utilities.
    namespace tmp {
        // Copy reference/cv-qualifiers from From to To.
        template <class From, class To>
        using forward_reference =
            std::conditional_t<std::is_lvalue_reference_v<From>, To&,
                               std::conditional_t<std::is_rvalue_reference_v<From>, To&&, To>>;
        template <class From, class To>
        using forward_cv = std::conditional_t<
            std::is_const_v<From>,
            std::conditional_t<std::is_volatile_v<From>, std::add_cv_t<To>, std::add_const_t<To>>,
            std::conditional_t<std::is_volatile_v<From>, std::add_volatile_t<To>, To>>;
        template <class From, class To>
        using forward_cvref =
            forward_reference<From, forward_cv<std::remove_reference_t<From>, To>>;

        // Check if Type is a specialization of Template.
        // Only works when Template only takes type template paramters.
        namespace detail {
            template <template <class...> class Template, class... Args>
            constexpr std::bool_constant<true> is_specialization_of_impl(Template<Args...>&&);
            constexpr std::bool_constant<false> is_specialization_of_impl(...);
        }
        template <class Type, template <class...> class Template>
        constexpr bool is_specialization_of() noexcept {
            return decltype(detail::is_specialization_of_impl<Template>(
                std::declval<Type>()))::value;
        }

        template <class... Types>
        struct typelist {
            static constexpr std::size_t size = sizeof...(Types);
        };

        template <class Type, class... Types>
        constexpr bool is_in(typelist<Types...>) noexcept {
            return (... || std::is_same_v<Type, Types>);
        }

        // Is {LhsTypes...} a subset of {RhsTypes...}?
        template <class... LhsTypes, class... RhsTypes>
        constexpr bool is_contained_in(typelist<LhsTypes...>, typelist<RhsTypes...>) noexcept {
            return (... && is_in<LhsTypes>(typelist<RhsTypes...>{}));
        }

        template <class Type, class... Types>
        constexpr std::size_t find_first_index(typelist<Types...>) noexcept {
            bool found = false;
            std::size_t count = 0;
            auto impl = [&](auto arg) {
                if (!found) {
                    if (std::is_same_v<Type, typename decltype(arg)::type>) {
                        found = true;
                        return;
                    }
                    ++count;
                }
            };
            (impl(std::type_identity<Types>{}), ...);
            return count;
        }

        namespace detail {
            template <auto...>
            struct value_placeholder {
                template <class T>
                constexpr value_placeholder(T&&) noexcept {}
            };

            template <std::size_t... I, class... Types>
            constexpr auto get_type_helper(std::index_sequence<I...>, typelist<Types...>) noexcept {
                return [](value_placeholder<I>..., auto nth, auto...) {
                    return nth;
                }(std::type_identity<Types>{}...);
            }

            template <class T, std::size_t index>
            struct indexed_type_placeholder {
                using type = T;

                template <class U>
                constexpr indexed_type_placeholder(U&&) noexcept {}
            };

            template <std::size_t... I, class... Types>
            constexpr auto back_sublist_helper(std::index_sequence<I...>,
                                               typelist<Types...>) noexcept {
                return [](value_placeholder<I>..., auto... args) {
                    return typelist<typename decltype(args)::type...>{};
                }(std::type_identity<Types>{}...);
            }
        }

        template <std::size_t N, class Typelist>
        using get_type = typename decltype(detail::get_type_helper(std::make_index_sequence<N>{},
                                                                   Typelist{}))::type;

        template <std::size_t N, class Typelist>
        using back_sublist = decltype(detail::back_sublist_helper(
            std::make_index_sequence<Typelist::size - N>{}, Typelist{}));

        namespace detail {
            template <auto prefix_sum, class... Types>
            constexpr auto prefix_sum_compaction(typelist<Types...>) noexcept {
                if constexpr (sizeof...(Types) == 0) {
                    return typelist<>{};
                }
                else {
                    using list = typelist<Types...>;
                    constexpr std::size_t size = prefix_sum[sizeof...(Types) - 1];
                    constexpr auto index_array = [] {
                        util::array<std::size_t, size> result{};
                        for (std::size_t i = sizeof...(Types); i > 0; --i) {
                            result[prefix_sum[i - 1] - 1] = i - 1;
                        }
                        return result;
                    }();

                    return [&index_array]<std::size_t... I>(std::index_sequence<I...>) {
                        return typelist<get_type<index_array[I], list>...>{};
                    }(std::make_index_sequence<size>{});
                }
            }

            template <class... Types>
            constexpr auto remove_duplicate_impl(typelist<Types...>) noexcept {
                if constexpr (sizeof...(Types) == 0) {
                    return typelist<>{};
                }
                else {
                    using list = typelist<Types...>;
                    constexpr auto prefix_sum = [] {
                        std::size_t count = 0;
                        std::size_t index = 0;
                        util::array<std::size_t, sizeof...(Types)> result{};
                        auto impl = [&](auto arg) {
                            if (find_first_index<typename decltype(arg)::type>(list{}) == index) {
                                result[index] = ++count;
                            }
                            else {
                                result[index] = count;
                            }
                            ++index;
                        };
                        (impl(std::type_identity<Types>{}), ...);
                        return result;
                    }();

                    return prefix_sum_compaction<prefix_sum>(list{});
                }
            }
        }

        // Guranteed to preserve the order and preserve the earliest item when duplicated.
        template <class Typelist>
        using remove_duplicate = decltype(detail::remove_duplicate_impl(Typelist{}));

        template <class... Types>
        constexpr bool has_duplicate(typelist<Types...>) noexcept {
            using list = typelist<Types...>;
            std::size_t index = 0;
            bool result = false;
            auto impl = [&](auto arg) {
                if (find_first_index<typename decltype(arg)::type>(list{}) != index) {
                    result = true;
                }
                ++index;
            };
            (impl(std::type_identity<Types>{}), ...);
            return result;
        }

        namespace detail {
            template <class Type, class... Types>
            constexpr typelist<Types..., Type> push_back_impl(typelist<Types...>) noexcept {
                return {};
            }

            constexpr typelist<> join_impl() noexcept { return {}; }
            template <class... Types>
            constexpr typelist<Types...> join_impl(typelist<Types...>) noexcept {
                return {};
            }
            template <class... Types1, class... Types2>
            constexpr typelist<Types1..., Types2...> join_impl(typelist<Types1...>,
                                                               typelist<Types2...>) noexcept {
                return {};
            }
            template <class... Types1, class... Types2, class... Typelists>
            constexpr auto join_impl(typelist<Types1...> first, typelist<Types2...> second,
                                     Typelists... remaining) noexcept {
                return join_impl(join_impl(first, second), join_impl(remaining...));
            }
        }

        template <class Typelist, class Type>
        using push_back = decltype(detail::push_back_impl<Type>(Typelist{}));

        template <class... Typelists>
        using join = decltype(detail::join_impl(Typelists{}...));

        namespace detail {
            template <class Predicate, class... Types>
            constexpr auto filter_impl(typelist<Types...>) noexcept {
                using list = typelist<Types...>;
                constexpr auto prefix_sum = [] {
                    std::size_t count = 0;
                    std::size_t index = 0;
                    util::array<std::size_t, sizeof...(Types)> result{};
                    auto impl = [&](auto arg) {
                        if (Predicate{}(arg)) {
                            result[index] = ++count;
                        }
                        else {
                            result[index] = count;
                        }
                        ++index;
                    };
                    (impl(std::type_identity<Types>{}), ...);
                    return result;
                }();

                return prefix_sum_compaction<prefix_sum>(list{});
            }
        }

        // Guaranteed to preserve the order.
        // The predicate is evaluated in the form Predicate{}(std::type_identity<T>{}) on type T.
        template <class Typelist, class Predicate>
        using filter = decltype(detail::filter_impl<Predicate>(Typelist{}));

        namespace detail {
            template <template <class...> class Map, class... BindParameters, class... T>
            constexpr auto map_impl(tmp::typelist<T...>) noexcept {
                return tmp::typelist<Map<T, BindParameters...>...>{};
            }
        }

        // Apply Map< . , BindParameters...> into Typelist.
        template <class Typelist, template <class...> class Map, class... BindParameters>
        using map = decltype(detail::map_impl<Map, BindParameters...>(Typelist{}));

        // Turn an array into std::integer_sequence.
        namespace detail {
            template <auto packed, std::size_t... I>
            constexpr auto unpack_array_helper(std::index_sequence<I...>) noexcept {
                return std::integer_sequence<typename decltype(packed)::value_type, packed[I]...>{};
            }
        }
        template <auto packed>
        using unpack_array = decltype(detail::unpack_array_helper<packed>(
            std::make_index_sequence<packed.size()>{}));
    }
}

#endif
