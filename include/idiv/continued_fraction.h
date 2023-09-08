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
        template <class Impl, template <class, class> class... Mixins>
        class continued_fraction : public Mixins<Impl, continued_fraction<Impl, Mixins...>>... {
        public:
            using impl_type = Impl;
            using partial_fraction_type = typename Impl::partial_fraction_type;
            using convergent_type = typename Impl::convergent_type;
            using interval_type = typename Impl::interval_type;

        private:
            Impl impl_;
            bool terminated_ = false;

            struct callback_type {
            private:
                continued_fraction& generator_;

                friend continued_fraction;
                explicit constexpr callback_type(continued_fraction& generator) noexcept
                    : generator_{generator} {}

            public:
                constexpr continued_fraction const& get_generator() const noexcept {
                    return generator_;
                }

                constexpr void operator()(partial_fraction_type const& next_partial_fraction) {
                    (static_cast<Mixins<Impl, continued_fraction>&>(generator_)
                         .update(next_partial_fraction, generator_.impl_),
                     ...);
                    generator_.terminated_ = false;
                }
            };

        public:
            explicit constexpr continued_fraction(Impl impl)
                : Mixins<Impl, continued_fraction>{static_cast<Impl const&>(impl)}...,
                  impl_{static_cast<Impl&&>(impl)} {}

            // Returns true if succeeded obtaining a further partial fraction.
            constexpr bool update() {
                if (!terminated_) {
                    terminated_ = true;
                    impl_.with_next_partial_fraction(callback_type{*this});

                    if (terminated_) {
                        auto invoke_final_update = [this](auto&& mixin) {
                            if constexpr (requires { mixin.final_update(impl_); }) {
                                mixin.final_update(impl_);
                            }
                        };
                        (invoke_final_update(static_cast<Mixins<Impl, continued_fraction>&>(*this)),
                         ...);
                    }
                }
                return !terminated_;
            }

            constexpr bool terminated() const noexcept { return terminated_; }
        };

        template <template <class, class> class... Mixins, class Impl>
        constexpr auto make_continued_fraction_generator(Impl&& impl) {
            return continued_fraction<std::remove_cvref_t<Impl>, Mixins...>{
                static_cast<Impl&&>(impl)};
        }

        template <class Impl, class ContinuedFractionGenerator>
        class index_tracker {
            int current_index_ = -1;

            friend ContinuedFractionGenerator;

            explicit constexpr index_tracker(Impl const&) noexcept {}

            constexpr void update(auto&&, Impl const&) noexcept { ++current_index_; }

        public:
            constexpr int current_index() const noexcept { return current_index_; }
        };

        template <class Impl, class ContinuedFractionGenerator>
        class partial_fraction_tracker {
            using partial_fraction_type = typename Impl::partial_fraction_type;

            partial_fraction_type current_partial_fraction_;

            friend ContinuedFractionGenerator;

            explicit constexpr partial_fraction_tracker(Impl const& impl)
                : current_partial_fraction_{impl.initial_partial_fraction()} {}

            constexpr void update(partial_fraction_type const& next_partial_fraction, Impl const&) {
                current_partial_fraction_ = next_partial_fraction;
            }

        public:
            constexpr partial_fraction_type const& current_partial_fraction() const noexcept {
                return current_partial_fraction_;
            }
        };

        template <class Impl, class ContinuedFractionGenerator>
        class convergent_tracker {
            using partial_fraction_type = typename Impl::partial_fraction_type;
            using convergent_type = typename Impl::convergent_type;

            convergent_type current_convergent_{1, 0u};
            convergent_type previous_convergent_{0, 1u};

            friend ContinuedFractionGenerator;

            explicit constexpr convergent_tracker(Impl const&) noexcept {}

            constexpr void update(partial_fraction_type const& next_partial_fraction, Impl const&) {
                auto next_numerator =
                    next_partial_fraction.denominator * current_convergent_numerator() +
                    next_partial_fraction.numerator * previous_convergent_numerator();
                auto next_denominator =
                    next_partial_fraction.denominator * current_convergent_denominator() +
                    next_partial_fraction.numerator * previous_convergent_denominator();

                previous_convergent_ = static_cast<convergent_type&&>(current_convergent_);
                current_convergent_ =
                    convergent_type{std::move(next_numerator), abs(std::move(next_denominator))};
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

        template <class Impl, class ContinuedFractionGenerator>
        class interval_tracker {
            using partial_fraction_type = typename Impl::partial_fraction_type;
            using convergent_type = typename Impl::convergent_type;
            using interval_type = typename Impl::interval_type;

            interval_type current_interval_;

            friend ContinuedFractionGenerator;

            explicit constexpr interval_tracker(Impl const& impl)
                : current_interval_{impl.initial_interval()} {}

            constexpr void update(partial_fraction_type const&, Impl& impl) {
                auto const& self = static_cast<ContinuedFractionGenerator const&>(*this);
                if constexpr (requires { impl.next_interval(self); }) {
                    current_interval_ = impl.next_interval(self);
                }
                else {
                    // Use convergents. Note that this is valid only for regular continued
                    // fractions.
                    if (self.current_index() >= 1) {
                        if (self.current_index() % 2 == 0) {
                            current_interval_ =
                                cyclic_interval<typename interval_type::value_type,
                                                cyclic_interval_type_t::left_closed_right_open>{
                                    self.current_convergent(), self.previous_convergent()};
                        }
                        else {
                            current_interval_ =
                                cyclic_interval<typename interval_type::value_type,
                                                cyclic_interval_type_t::left_open_right_closed>{
                                    self.previous_convergent(), self.current_convergent()};
                        }
                    }
                }
            }
            constexpr void final_update(Impl const&) {
                current_interval_ = cyclic_interval<typename interval_type::value_type,
                                                    cyclic_interval_type_t::single_point>{
                    static_cast<ContinuedFractionGenerator const&>(*this).current_convergent()};
            }

        public:
            interval_type const& current_interval() const noexcept { return current_interval_; }

            // Refine the current value so that the maximum possible error is no more than the
            // given bound.
            template <class ErrorValue>
            convergent_type const& progress_until(ErrorValue const& error_bound) {
                // Translate the lower bound by error_bound and see if it belongs to the interval.
                while (current_interval().visit([&](auto&& itv) {
                    static_assert(itv.interval_type() != cyclic_interval_type_t::empty);
                    if constexpr (itv.interval_type() == cyclic_interval_type_t::single_point) {
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
                    static_cast<ContinuedFractionGenerator&>(*this).update();
                }
                return static_cast<ContinuedFractionGenerator const&>(*this).current_convergent();
            }
        };
    }
}

#endif
