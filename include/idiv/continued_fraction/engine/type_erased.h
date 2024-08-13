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

#ifndef JKJ_HEADER_TYPE_ERASED_CONTINUED_FRACTION
#define JKJ_HEADER_TYPE_ERASED_CONTINUED_FRACTION

#include "continued_fraction.h"
#include <memory>

namespace jkj {
    namespace cntfrc {
        namespace impl {
            // Type-erased wrapper for any continued fraction implementation classes.
            // It does not support reference types.
            template <class PartialFractionType, class ConvergentType, class IntervalType>
            class type_erased {
            public:
                using partial_fraction_type = PartialFractionType;
                using convergent_type = ConvergentType;
                using interval_type = IntervalType;

            private:
                // We do not need a virtual destructor for this.
                struct callback_wrapper_base {
                    constexpr virtual void on_next_partial_fraction(
                        partial_fraction_type const& next_partial_fraction) const = 0;

                    constexpr virtual void
                    on_next_interval(interval_type const& next_interval) const = 0;
                };

                template <class Callback>
                class callback_wrapper : public callback_wrapper_base {
                    Callback& callback_;

                public:
                    constexpr callback_wrapper(Callback& callback) : callback_{callback} {}

                    constexpr void on_next_partial_fraction(
                        partial_fraction_type const& next_partial_fraction) const override {
                        callback_.on_next_partial_fraction(next_partial_fraction);
                    }

                    constexpr void
                    on_next_interval(interval_type const& next_interval) const override {
                        callback_.on_next_interval(next_interval);
                    }
                };

                struct base {
                    constexpr virtual partial_fraction_type initial_partial_fraction() const = 0;
                    constexpr virtual interval_type initial_interval() const = 0;
                    constexpr virtual void
                    with_next_partial_fraction(callback_wrapper_base& callback) = 0;
                    virtual std::unique_ptr<base> clone() const = 0;
                    virtual ~base() = default;
                };

                template <class Impl>
                class impl_wrapper : public base {
                    Impl impl_;

                public:
                    constexpr impl_wrapper(Impl&& impl) : impl_{std::move(impl)} {}

                    constexpr partial_fraction_type initial_partial_fraction() const override {
                        return impl_.initial_partial_fraction();
                    }
                    constexpr interval_type initial_interval() const override {
                        return impl_.initial_interval();
                    }
                    constexpr void
                    with_next_partial_fraction(callback_wrapper_base& callback) override {
                        impl_.with_next_partial_fraction(callback);
                    }
                    std::unique_ptr<base> clone() const override {
                        return std::make_unique<impl_wrapper>(*this);
                    }
                };

                std::unique_ptr<base> impl_ptr_;

            public:
                template <class Impl>
                    requires(!std::is_same_v<Impl, type_erased>)
                type_erased(Impl impl)
                    : impl_ptr_{std::make_unique<impl_wrapper<Impl>>(std::move(impl))} {
                    static_assert(
                        std::is_same_v<typename Impl::partial_fraction_type, partial_fraction_type>,
                        "the implementation type does not define compatible partial_fraction_type");
                    static_assert(
                        std::is_same_v<typename Impl::convergent_type, convergent_type>,
                        "the implementation type does not define compatible convergent_type");
                    static_assert(
                        std::is_same_v<typename Impl::interval_type, interval_type>,
                        "the implementation type does not define compatible interval_type");
                }

                constexpr type_erased(type_erased const& other)
                    : impl_ptr_{other.impl_ptr_->clone()} {}

                constexpr type_erased(type_erased&& other)
                    : impl_ptr_{std::move(other.impl_ptr_)} {}

                type_erased& operator=(type_erased const& other) {
                    impl_ptr_ = other.impl_ptr_->clone();
                    return *this;
                }

                constexpr type_erased& operator=(type_erased&& other) {
                    std::swap(impl_ptr_, other.impl_ptr_);
                    return *this;
                }

                constexpr partial_fraction_type initial_partial_fraction() const {
                    return impl_ptr_->initial_partial_fraction();
                }

                constexpr interval_type initial_interval() const {
                    return impl_ptr_->initial_interval();
                }

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback& callback) {
                    // Do double indirection.
                    callback_wrapper<Callback> wrapped_callback{callback};
                    impl_ptr_->with_next_partial_fraction(wrapped_callback);
                }
            };
        }
    }
}

#endif
