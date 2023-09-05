// Copyright 2022-2023 Junekey Jeon
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

#ifndef JKJ_HEADER_CONTINUED_FRACTION
#define JKJ_HEADER_CONTINUED_FRACTION

#include "interval.h"
#include "projective_rational.h"
#include "util.h"

namespace jkj {
    // An interface for generalized continued fraction calculator for real numbers.
    // Partial numerators and denominators can be signed integers, but the denominator of the
    // convergents are assumed to be always nonnegative.
    // Given continued fraction expansion
    //
    // b0 + a1 / (b1 + a2 / (b2 + a3 / (b3 + ... ) ) ),
    //
    // we call an/bn the nth "partial fraction", and also we call the fraction obtained by
    // truncating the continued fraction at the nth partial fraction as the nth "convergent".

    namespace cntfrc {
        // Declare needed type aliases by specializing this template.
        template <class Impl>
        struct continued_fraction_traits;

        struct mixin_initializer_token {};

        template <class Impl, template <class> class... Mixins>
        class continued_fraction_base : public continued_fraction_traits<Impl>,
                                        public Mixins<Impl>... {
            bool terminated_ = false;

        public:
            using traits_type = continued_fraction_traits<Impl>;

            template <class MixinInitializer>
            constexpr continued_fraction_base(MixinInitializer&& mixin_initializer = {})
                : Mixins<Impl>{mixin_initializer_token{}, mixin_initializer}... {}

            // Returns true if succeeded obtaining a further partial fraction.
            constexpr bool update() {
                if (!terminated_) {
                    terminated_ = !static_cast<Impl&>(*this).with_next_partial_fraction(
                        [&](auto&& next_partial_fraction) {
                            (static_cast<Mixins<Impl>&>(*this).update(static_cast<Impl&>(*this),
                                                                      next_partial_fraction),
                             ...);
                        });

                    if (terminated_) {
                        auto invoke_final_update = [&impl =
                                                        static_cast<Impl&>(*this)](auto&& mixin) {
                            if constexpr (requires { mixin.final_update(impl); }) {
                                mixin.final_update(impl);
                            }
                        };
                        (invoke_final_update(static_cast<Mixins<Impl>&>(*this)), ...);
                    }
                }
                return !terminated_;
            }

            constexpr bool terminated() const noexcept { return terminated_; }
        };

        template <class Impl>
        class index_tracker {
            int current_index_ = -1;

            template <class, template <class> class...>
            friend class continued_fraction_base;

            template <class MixinInitializer>
            constexpr index_tracker(mixin_initializer_token, MixinInitializer&&) {}

            constexpr void update(Impl&, auto&&) { ++current_index_; }

        public:
            constexpr int current_index() const noexcept { return current_index_; }
        };

        template <class Impl>
        class partial_fraction_tracker {
            using partial_fraction_type =
                typename continued_fraction_traits<Impl>::partial_fraction_type;

            partial_fraction_type current_partial_fraction_;

            template <class, template <class> class...>
            friend class continued_fraction_base;

            template <class MixinInitializer>
            constexpr partial_fraction_tracker(mixin_initializer_token,
                                               MixinInitializer&& mixin_initializer)
                : current_partial_fraction_{mixin_initializer.initial_partial_fraction()} {}

            constexpr void update(Impl&, partial_fraction_type const& next_partial_fraction) {
                current_partial_fraction_ = next_partial_fraction;
            }

        public:
            constexpr partial_fraction_type const& current_partial_fraction() const noexcept {
                return current_partial_fraction_;
            }
        };

        template <class Impl>
        class convergent_tracker {
            using convergent_type = typename continued_fraction_traits<Impl>::convergent_type;

            convergent_type current_convergent_{1, 0u};
            convergent_type previous_convergent_{0, 1u};

            template <class, template <class> class...>
            friend class continued_fraction_base;

            template <class MixinInitializer>
            constexpr convergent_tracker(mixin_initializer_token, MixinInitializer&&) {}

            template <class PartialFractionType>
            constexpr void update(Impl&, PartialFractionType const& next_partial_fraction) {
                auto next_numerator =
                    next_partial_fraction.denominator * current_convergent_numerator() +
                    next_partial_fraction.numerator * previous_convergent_numerator();
                auto next_denominator =
                    next_partial_fraction.denominator * current_convergent_denominator() +
                    next_partial_fraction.numerator * previous_convergent_denominator();

                previous_convergent_ = static_cast<convergent_type&&>(current_convergent_);
                current_convergent_ = convergent_type{
                    static_cast<decltype(next_numerator)&&>(next_numerator),
                    abs(static_cast<decltype(next_denominator)&&>(next_denominator))};
            }

        public:
            constexpr convergent_type const& current_convergent() const noexcept {
                return current_convergent_;
            }
            constexpr auto const& current_convergent_numerator() const noexcept {
                return current_convergent().numerator;
            }
            constexpr auto const& current_convergent_denominator() const noexcept {
                return current_convergent().denominator;
            }

            constexpr convergent_type const& previous_convergent() const noexcept {
                return previous_convergent_;
            }
            constexpr auto const& previous_convergent_numerator() const noexcept {
                return previous_convergent().numerator;
            }
            constexpr auto const& previous_convergent_denominator() const noexcept {
                return previous_convergent().denominator;
            }
        };

        template <class Impl>
        class interval_tracker {
            using convergent_type = typename continued_fraction_traits<Impl>::convergent_type;
            using interval_type = typename continued_fraction_traits<Impl>::interval_type;

            interval_type current_interval_;

            template <class, template <class> class...>
            friend class continued_fraction_base;

            template <class MixinInitializer>
            constexpr interval_tracker(mixin_initializer_token,
                                       MixinInitializer&& mixin_initializer)
                : current_interval_{mixin_initializer.initial_interval()} {}

            template <class PartialFractionType>
            constexpr void update(Impl& impl, PartialFractionType const&) {
                if constexpr (requires { impl.next_interval(); }) {
                    current_interval_ = impl.next_interval();
                }
                else {
                    // Use convergents. Note that this is valid only for regular continued
                    // fractions.
                    if (impl.current_index() >= 1) {
                        if (impl.current_index() % 2 == 0) {
                            current_interval_ =
                                cyclic_interval<typename interval_type::value_type,
                                                cyclic_interval_type_t::left_closed_right_open>{
                                    impl.current_convergent(), impl.previous_convergent()};
                        }
                        else {
                            current_interval_ =
                                cyclic_interval<typename interval_type::value_type,
                                                cyclic_interval_type_t::left_open_right_closed>{
                                    impl.previous_convergent(), impl.current_convergent()};
                        }
                    }
                }
            }
            constexpr void final_update(Impl& impl) {
                current_interval_ = cyclic_interval<typename interval_type::value_type,
                                                    cyclic_interval_type_t::single_point>{
                    impl.current_convergent()};
            }

        public:
            interval_type const& current_interval() const noexcept { return current_interval_; }

            // Refine the current value so that the maximum possible error is no more than the
            // given bound.
            template <class ErrorValue>
            convergent_type const& progress_until(ErrorValue const& error_bound) {
                // Translate the lower bound by error_bound and see if it belongs to the interval.
                while (current_interval().visit([&](auto&& itv) {
                    if constexpr (itv.interval_type() == cyclic_interval_type_t::empty) {
                        util::constexpr_assert(false);
                        return false;
                    }
                    else if constexpr (itv.interval_type() ==
                                       cyclic_interval_type_t::single_point) {
                        return false;
                    }
                    else if constexpr (itv.interval_type() == cyclic_interval_type_t::entire) {
                        return true;
                    }
                    else {
                        auto translated = linear_fractional_translation(
                            error_bound.numerator, error_bound.denominator)(itv.lower_bound());

                        if (!itv.contains(translated)) {
                            return false;
                        }

                        if (itv.interval_type() == cyclic_interval_type_t::closed &&
                            translated == itv.upper_bound()) {
                            return false;
                        }

                        return true;
                    }
                })) {
                    static_cast<Impl&>(*this).update();
                }
                return static_cast<Impl&>(*this).current_convergent();
            }
        };
    }
}

#endif
