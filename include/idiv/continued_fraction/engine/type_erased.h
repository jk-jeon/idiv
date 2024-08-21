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

#ifndef JKJ_HEADER_IDIV_CONTINUED_FRACTION_ENGINE_TYPE_ERASED
#define JKJ_HEADER_IDIV_CONTINUED_FRACTION_ENGINE_TYPE_ERASED

#include "../generator.h"
#include "../../unique_ptr.h"

namespace jkj {
    namespace cntfrc {
        template <class PartialFractionType>
        struct default_type_erased_engine_traits {
            using partial_fraction_type = PartialFractionType;
            using required_mixins = tmp::typelist<>;

            struct abstract_interface {};

            template <class Base>
            struct interface_mixin : Base {};
        };

        namespace engine {
            template <class TypeErasedEngineTraits>
            class type_erased;
        }

        namespace type_erased_detail {
            struct empty_base {};

            template <class TypeErasedEngineTraits>
            struct impl {
                using partial_fraction_type =
                    typename TypeErasedEngineTraits::partial_fraction_type;
                using required_mixins = typename TypeErasedEngineTraits::required_mixins;

                template <class Mixin>
                using tracking_data_type = typename Mixin::template traits<
                    engine::type_erased<TypeErasedEngineTraits>>::tracking_data;

                template <class Mixin>
                struct proxy_mixin_wrapper_base {
                    constexpr virtual tracking_data_type<Mixin> const& tracking_data() const = 0;
                };

                template <class MixinList = required_mixins>
                class proxy_type;

                template <class... Mixins>
                class proxy_type<tmp::typelist<Mixins...>>
                    : util::noncopyable<>,
                      util::nonmovable<>,
                      public Mixins::template proxy_mixin<
                          engine::type_erased<TypeErasedEngineTraits>, proxy_type>...,
                      virtual private proxy_mixin_wrapper_base<Mixins>... {
                    template <class Mixin>
                    friend struct tracking_data_accessor;

                    template <class Mixin>
                    constexpr tracking_data_type<Mixin> const& tracking_data() const noexcept {
                        return static_cast<proxy_mixin_wrapper_base<Mixin> const&>(*this)
                            .tracking_data();
                    }

                public:
                    using engine_type = engine::type_erased<TypeErasedEngineTraits>;
                    using included_mixins = tmp::typelist<Mixins...>;

                    constexpr virtual ~proxy_type() noexcept {}
                };

                template <class Mixin, class Derived>
                struct proxy_wrapper_fragment : virtual public proxy_mixin_wrapper_base<Mixin> {
                    constexpr tracking_data_type<Mixin> const& tracking_data() const override {
                        return static_cast<Derived const&>(*this).template tracking_data<Mixin>();
                    }
                };

                template <class Proxy, class MixinList = required_mixins>
                class proxy_wrapper;

                template <class Proxy, class... Mixins>
                class proxy_wrapper<Proxy, tmp::typelist<Mixins...>>
                    : public proxy_type<tmp::typelist<Mixins...>>,
                      private proxy_wrapper_fragment<Mixins, proxy_wrapper<Proxy>>... {
                    template <class Mixin, class Derived>
                    friend struct proxy_wrapper_fragment;

                    template <class Mixin>
                    constexpr tracking_data_type<Mixin> const& tracking_data() const {
                        return tracking_data_accessor<Mixin>::template tracking_data<Proxy>(
                            static_cast<proxy_mixin_of<
                                Mixin, engine::type_erased<TypeErasedEngineTraits>, Proxy> const&>(
                                state_));
                    }

                    Proxy state_;

                public:
                    template <class Callback>
                    constexpr proxy_wrapper(Callback&& callback)
                        : state_{callback.current_state()} {}
                };

                // We do not need virtual destructors for these, because the callback object is only
                // constructed as a local variable of type_erased::with_next_partial_fraction and
                // type_erased::with_zeroth_partial_fraction.
                template <callback_type_tag_t callback_type_tag, class MixinList = required_mixins>
                class callback_type;

                template <class Mixin>
                struct callback_mixin_wrapper_base {
                    constexpr virtual tracking_data_type<Mixin> const& tracking_data() const = 0;
                };

                template <callback_type_tag_t callback_type_tag, class Callback,
                          class MixinList = required_mixins>
                class callback_wrapper;

                template <callback_type_tag_t callback_type_tag, class... Mixins>
                class callback_type<callback_type_tag, tmp::typelist<Mixins...>>
                    : util::noncopyable<>,
                      util::nonmovable<>,
                      public Mixins::template callback_mixin<
                          callback_type<callback_type_tag, tmp::typelist<Mixins...>>>...,
                      virtual private callback_mixin_wrapper_base<Mixins>... {
                    template <class Mixin>
                    friend struct tracking_data_accessor;

                    template <class Mixin>
                    constexpr tracking_data_type<Mixin> const& tracking_data() const noexcept {
                        return static_cast<callback_mixin_wrapper_base<Mixin> const&>(*this)
                            .tracking_data();
                    }

                    template <callback_type_tag_t, class, class>
                    friend class callback_wrapper;

                    util::unique_ptr<proxy_type<>> state_;

                    template <class Callback>
                    explicit constexpr callback_type(Callback&& callback)
                        : state_{
                              util::make_unique<proxy_wrapper<decltype(callback.current_state())>>(
                                  callback)} {}

                public:
                    using engine_type = engine::type_erased<TypeErasedEngineTraits>;
                    using included_mixins = tmp::typelist<Mixins...>;

                    constexpr virtual void on_next_partial_fraction(
                        partial_fraction_type const& next_partial_fraction) const = 0;

                    constexpr proxy_type<> const& current_state() const noexcept { return *state_; }
                };

                template <class Mixin, class Derived>
                struct callback_wrapper_fragment
                    : virtual public callback_mixin_wrapper_base<Mixin> {
                    constexpr tracking_data_type<Mixin> const& tracking_data() const override {
                        return static_cast<Derived const&>(*this).template tracking_data<Mixin>();
                    }
                };

                template <callback_type_tag_t callback_type_tag, class Callback, class... Mixins>
                class callback_wrapper<callback_type_tag, Callback, tmp::typelist<Mixins...>>
                    : public callback_type<callback_type_tag>,
                      private callback_wrapper_fragment<
                          Mixins, callback_wrapper<callback_type_tag, Callback,
                                                   tmp::typelist<Mixins...>>>... {
                    template <class Mixin, class Derived>
                    friend struct callback_wrapper_fragment;

                    template <class Mixin>
                    constexpr tracking_data_type<Mixin>& tracking_data() const {
                        return tracking_data_accessor<Mixin>::template tracking_data<Callback>(
                            static_cast<callback_mixin_of<
                                Mixin, engine::type_erased<TypeErasedEngineTraits>,
                                Callback> const&>(callback_));
                    }

                    Callback& callback_;

                public:
                    constexpr callback_wrapper(Callback& callback)
                        : callback_type<callback_type_tag, tmp::typelist<Mixins...>>{callback},
                          callback_{callback} {}

                    constexpr void on_next_partial_fraction(
                        partial_fraction_type const& next_partial_fraction) const override {
                        callback_.on_next_partial_fraction(next_partial_fraction);
                    }
                };

                using abstract_interface = typename TypeErasedEngineTraits::abstract_interface;

                struct engine_wrapper_base {
                    constexpr virtual abstract_interface& get_abstract_interface() = 0;
                    constexpr virtual void with_next_partial_fraction(
                        callback_type<callback_type_tag_t::normal>& callback) = 0;
                    constexpr virtual void with_zeroth_partial_fraction(
                        callback_type<callback_type_tag_t::normal>& callback) = 0;
                    constexpr virtual util::unique_ptr<engine_wrapper_base> clone() const = 0;
                    constexpr virtual ~engine_wrapper_base() = default;
                };

                template <class Engine>
                struct middle_wrapper {
                    Engine engine_;

                    constexpr middle_wrapper(Engine&& engine) : engine_{std::move(engine)} {}

                    constexpr Engine& implementation() noexcept { return engine_; }
                    constexpr Engine const& implementation() const noexcept { return engine_; }
                };

                template <class Engine>
                class engine_wrapper : public engine_wrapper_base,
                                       private TypeErasedEngineTraits::template interface_mixin<
                                           middle_wrapper<Engine>, abstract_interface> {

                    using mixin_base_type =
                        typename TypeErasedEngineTraits::template interface_mixin<
                            middle_wrapper<Engine>, abstract_interface>;

                public:
                    using mixin_base_type::mixin_base_type;

                    constexpr abstract_interface& get_abstract_interface() override {
                        return *this;
                    }

                    constexpr void with_next_partial_fraction(
                        callback_type<callback_type_tag_t::normal>& callback) override {
                        mixin_base_type::implementation().with_next_partial_fraction(callback);
                    }

                    constexpr void with_zeroth_partial_fraction(
                        callback_type<callback_type_tag_t::normal>& callback) override {
                        if constexpr (requires {
                                          mixin_base_type::implementation()
                                              .with_zeroth_partial_fraction(callback);
                                      }) {
                            mixin_base_type::implementation().with_zeroth_partial_fraction(
                                callback);
                        }
                        else {
                            mixin_base_type::implementation().with_next_partial_fraction(callback);
                        }
                    }

                    constexpr util::unique_ptr<engine_wrapper_base> clone() const override {
                        return util::make_unique<engine_wrapper>(*this);
                    }
                };

                class interface_adapter {
                    util::unique_ptr<engine_wrapper_base> impl_ptr_;

                public:
                    template <class Engine>
                    constexpr interface_adapter(Engine&& engine)
                        : impl_ptr_{util::make_unique<engine_wrapper<std::remove_cvref_t<Engine>>>(
                              static_cast<Engine&&>(engine))} {}

                    constexpr interface_adapter(interface_adapter const& other)
                        : impl_ptr_{other.impl_ptr_->clone()} {}

                    constexpr interface_adapter(interface_adapter&& other)
                        : impl_ptr_{std::move(other.impl_ptr_)} {}

                    constexpr interface_adapter& operator=(interface_adapter const& other) {
                        impl_ptr_ = other.impl_ptr_->clone();
                        return *this;
                    }

                    constexpr interface_adapter& operator=(interface_adapter&& other) {
                        std::swap(impl_ptr_, other.impl_ptr_);
                        return *this;
                    }

                    template <class Callback>
                    constexpr void with_zeroth_partial_fraction(Callback&& callback) {
                        // Do double indirection.
                        callback_wrapper<callback_type_tag_t::normal, Callback> wrapped_callback{
                            callback};
                        impl_ptr_->with_zeroth_partial_fraction(wrapped_callback);
                    }

                    template <class Callback>
                    constexpr void with_next_partial_fraction(Callback&& callback) {
                        // Do double indirection.
                        callback_wrapper<callback_type_tag_t::normal, Callback> wrapped_callback{
                            callback};
                        impl_ptr_->with_next_partial_fraction(wrapped_callback);
                    }

                private:
                    friend
                        typename TypeErasedEngineTraits::template interface_mixin<interface_adapter,
                                                                                  empty_base>;

                    constexpr abstract_interface& implementation() noexcept {
                        return impl_ptr_->get_abstract_interface();
                    }
                    constexpr abstract_interface const& implementation() const noexcept {
                        return impl_ptr_->get_abstract_interface();
                    }
                };
            };
        }

        namespace engine {
            // Type-erased wrapper for any continued fraction engine classes.
            template <class TypeErasedEngineTraits>
            class type_erased
                : public TypeErasedEngineTraits::template interface_mixin<
                      typename type_erased_detail::impl<TypeErasedEngineTraits>::interface_adapter,
                      typename type_erased_detail::empty_base> {

                using base_type = typename TypeErasedEngineTraits::template interface_mixin<
                    typename type_erased_detail::impl<TypeErasedEngineTraits>::interface_adapter,
                    typename type_erased_detail::empty_base>;

            public:
                using partial_fraction_type =
                    typename TypeErasedEngineTraits::partial_fraction_type;
                using required_mixins = typename TypeErasedEngineTraits::required_mixins;

                template <class Engine>
                    requires(!std::is_same_v<Engine, type_erased>)
                constexpr type_erased(Engine engine) : base_type{std::move(engine)} {
                    static_assert(
                        std::is_same_v<typename Engine::partial_fraction_type,
                                       partial_fraction_type>,
                        "the engine type does not define compatible partial_fraction_type");
                    static_assert(
                        tmp::is_contained_in(typename find_required_mixin_list<Engine>::type{},
                                             required_mixins{}),
                        "some required mixin is not covered in the specified list of "
                        "required mixins");
                }

            private:
                template <class Engine, class MixinList = required_mixins>
                struct find_required_mixin_list;

                template <class Engine, class... Mixins>
                struct find_required_mixin_list<Engine, tmp::typelist<Mixins...>> {
                    using type = decltype(find_sorted_required_mixin_list<Engine, Mixins...>());
                };
            };
        }
    }
}

#endif
