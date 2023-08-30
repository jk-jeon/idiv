// Copyright 2023 Junekey Jeon
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

#ifndef JKJ_HEADER_GOSPER_CONTINUED_FRACTION
#define JKJ_HEADER_GOSPER_CONTINUED_FRACTION

#include "continued_fraction.h"
#include "frac.h"
#include "interval.h"
#include <cstdlib>
#include <type_traits>

namespace jkj {
    template <class Int>
    struct unary_gosper_coeff {
        // f(x) = (ax+b)/(cx+d).
        struct coeff {
            Int x_coeff;
            Int const_coeff;
        };
        coeff numerator;
        coeff denominator;
    };

    template <class Int>
    struct binary_gosper_coeff {
        struct coeff {
            Int const_coeff;
            Int x_coeff;
            Int y_coeff;
            Int xy_coeff;

            template <class PartialFraction>
            constexpr void feed_x(PartialFraction const& partial_fraction) {
                // s/t, (a, b, c, d) |-> (b, sa + tb, d, sc + td)
                const_coeff *= partial_fraction.numerator;
                const_coeff += partial_fraction.denominator * x_coeff;
                y_coeff *= partial_fraction.numerator;
                y_coeff += partial_fraction.denominator * xy_coeff;
                swap(const_coeff, x_coeff);
                swap(y_coeff, xy_coeff);
            }

            template <class PartialFraction>
            constexpr void feed_y(PartialFraction const& partial_fraction) {
                // s/t, (a, b, c, d) |-> (c, d, sa + tc, sb + td)
                const_coeff *= partial_fraction.numerator;
                const_coeff += partial_fraction.denominator * y_coeff;
                x_coeff *= partial_fraction.numerator;
                x_coeff += partial_fraction.denominator * xy_coeff;
                swap(const_coeff, y_coeff);
                swap(x_coeff, xy_coeff);
            }
        };

        coeff numerator;
        coeff denominator;
    };

    template <class ContinuedFractionImpl, class Unity = unity>
    class unary_gosper_continued_fraction {
    public:
        using int_type = decltype(ContinuedFractionImpl::convergent_type::numerator);
        using uint_type = decltype(ContinuedFractionImpl::convergent_type::denominator);
        using partial_fraction_type = frac<Unity, int_type>;
        using convergent_type = typename ContinuedFractionImpl::convergent_type;

    private:
        convergent_generator<ContinuedFractionImpl, interval_tracker> cf_;
        unary_gosper_coeff<int_type> coeff_;
        bool is_terminated_ = false;
        bool negative_determinant = true;

        constexpr void progress() {
            cf_.update();

            // If the current coefficient is the last one, we have to use a different update
            // formula from now on.
            if (cf_.is_terminated()) {
                is_terminated_ = true;
                coeff_.numerator.const_coeff =
                    coeff_.numerator.x_coeff * cf_.current_convergent().numerator +
                    coeff_.numerator.const_coeff * cf_.current_convergent().denominator;
                coeff_.denominator.const_coeff =
                    coeff_.denominator.x_coeff * cf_.current_convergent().numerator +
                    coeff_.denominator.const_coeff * cf_.current_convergent().denominator;

                if (is_strictly_negative(coeff_.denominator.const_coeff)) {
                    coeff_.numerator.const_coeff =
                        invert_sign(static_cast<int_type&&>(coeff_.numerator.const_coeff));
                    coeff_.denominator.const_coeff =
                        invert_sign(static_cast<int_type&&>(coeff_.denominator.const_coeff));
                }
            }
        }

    public:
        constexpr unary_gosper_continued_fraction(ContinuedFractionImpl cf,
                                                  unary_gosper_coeff<int_type> coeff)
            : cf_{static_cast<ContinuedFractionImpl>(cf)},
              coeff_{static_cast<unary_gosper_coeff<int_type>&&>(coeff)} {
            auto const determinant = coeff_.numerator.x_coeff * coeff_.denominator.const_coeff -
                                     coeff_.numerator.const_coeff * coeff_.denominator.x_coeff;

            if (is_strictly_positive(determinant)) {
                negative_determinant = false;
            }
            else if (is_strictly_negative(determinant)) {
                negative_determinant = true;
            }
            else {
                // Zero determinant.
                // Assumes cx + d is not zero.
                is_terminated_ = true;

                if (is_zero(coeff_.denominator.const_coeff)) {
                    util::constexpr_assert<util::error_msgs::divide_by_zero>(
                        !is_zero(coeff_.denominator.x_coeff));
                    coeff_.denominator.const_coeff =
                        static_cast<int_type&&>(coeff_.denominator.x_coeff);
                    coeff_.numerator.const_coeff =
                        static_cast<int_type&&>(coeff_.numerator.x_coeff);
                }
            }
        }

        constexpr next_partial_fraction_return<partial_fraction_type> next_partial_fraction() {
            using util::swap;

            while (true) {
                if (is_terminated_) {
                    // Proceed as in the case of usual rational continued fractions.
                    util::constexpr_assert<util::error_msgs::divide_by_zero>(
                        is_strictly_positive(coeff_.denominator.const_coeff));

                    using std::div;
                    auto div_result =
                        div(coeff_.numerator.const_coeff, abs(coeff_.denominator.const_coeff));
                    coeff_.numerator.const_coeff =
                        static_cast<int_type&&>(coeff_.denominator.const_coeff);
                    coeff_.denominator.const_coeff =
                        int_type(static_cast<uint_type&&>(div_result.rem));

                    // Terminate if all coefficients in the denominator has become zero.
                    return {partial_fraction_type{{}, static_cast<int_type&&>(div_result.quot)},
                            is_zero(coeff_.denominator.const_coeff)};
                }
                else {
                    // Read the current cyclic interval.
                    auto const domain_itv = cf_.current_interval();
                    // If the current interval is the entire real line, nothing can be done, so
                    // progress.
                    if (!is_zero(domain_itv.lower_bound().denominator) ||
                        !is_zero(domain_itv.upper_bound().denominator)) {
                        // Find the mapped cyclic interval.
                        cyclic_interval<frac<int_type, uint_type>, cyclic_interval_type_t::open>
                            itv{make_frac_from_signed(coeff_.numerator.x_coeff *
                                                              domain_itv.lower_bound().numerator +
                                                          coeff_.numerator.const_coeff *
                                                              domain_itv.lower_bound().denominator,
                                                      coeff_.denominator.x_coeff *
                                                              domain_itv.lower_bound().numerator +
                                                          coeff_.denominator.const_coeff *
                                                              domain_itv.lower_bound().denominator),
                                make_frac_from_signed(
                                    coeff_.numerator.x_coeff * domain_itv.upper_bound().numerator +
                                        coeff_.numerator.const_coeff *
                                            domain_itv.upper_bound().denominator,
                                    coeff_.denominator.x_coeff *
                                            domain_itv.upper_bound().numerator +
                                        coeff_.denominator.const_coeff *
                                            domain_itv.upper_bound().denominator)};

                        if (negative_determinant) {
                            swap(itv.lower_bound(), itv.upper_bound());
                        }

                        // Check if itv is a bounded nonempty open interval in the real line.
                        if (!is_zero(itv.lower_bound().denominator) &&
                            !is_zero(itv.upper_bound().denominator)) {
                            if (itv.lower_bound() < itv.upper_bound()) {
                                // See if the endpoints agree on the integer part.
                                auto lb = div_floor(itv.lower_bound().numerator,
                                                    itv.lower_bound().denominator);
                                auto ub = div_floor(itv.upper_bound().numerator,
                                                    itv.upper_bound().denominator);

                                if (lb == ub) {
                                    // Found the next coefficient.
                                    coeff_.numerator.x_coeff -= lb * coeff_.denominator.x_coeff;
                                    coeff_.numerator.const_coeff -=
                                        lb * coeff_.denominator.const_coeff;
                                    swap(coeff_.numerator, coeff_.denominator);
                                    negative_determinant = !negative_determinant;

                                    // Terminate if all coefficients in the denominator has become
                                    // zero.
                                    return {partial_fraction_type{{}, static_cast<int_type&&>(lb)},
                                            is_zero(coeff_.denominator.x_coeff) &&
                                                is_zero(coeff_.denominator.const_coeff)};
                                }
                            }
                        }
                    }

                    progress();
                }
            }
        }
    };

    template <class ContinuedFractionImplX, class ContinuedFractionImplY, class Unity = unity>
        requires(std::is_same_v<typename ContinuedFractionImplX::partial_fraction_type,
                                typename ContinuedFractionImplY::partial_fraction_type> &&
                 std::is_same_v<typename ContinuedFractionImplX::convergent_type,
                                typename ContinuedFractionImplY::convergent_type>)
    class binary_gosper_continued_fraction {
    public:
        using int_type = decltype(ContinuedFractionImplX::convergent_type::numerator);
        using uint_type = decltype(ContinuedFractionImplX::convergent_type::denominator);
        using partial_fraction_type = frac<Unity, int_type>;
        using convergent_type = typename ContinuedFractionImplX::convergent_type;

    private:
        ContinuedFractionImplX x_cf_;
        ContinuedFractionImplY y_cf_;
        binary_gosper_coeff<int_type> coeff_;
        bool is_x_terminated_ = false;
        bool is_y_terminated_ = false;

        constexpr void progress_x() {
            auto result = x_cf_.next_partial_fraction();
            coeff_.numerator.feed_x(result.partial_fraction);
            coeff_.denominator.feed_x(result.partial_fraction);

            // If the current coefficient is the last one, we have to use a different update
            // formula from now on.
            if (result.is_last) {
                is_x_terminated_ = true;
                using util::swap;
                swap(coeff_.numerator.const_coeff, coeff_.numerator.x_coeff);
                swap(coeff_.numerator.y_coeff, coeff_.numerator.xy_coeff);
                swap(coeff_.denominator.const_coeff, coeff_.denominator.x_coeff);
                swap(coeff_.denominator.y_coeff, coeff_.denominator.xy_coeff);
            }
        }
        constexpr void progress_y() {
            auto result = y_cf_.next_partial_fraction();
            coeff_.numerator.feed_y(result.partial_fraction);
            coeff_.denominator.feed_y(result.partial_fraction);

            // If the current coefficient is the last one, we have to use a different update
            // formula from now on.
            if (result.is_last) {
                is_y_terminated_ = true;
                using util::swap;
                swap(coeff_.numerator.const_coeff, coeff_.numerator.y_coeff);
                swap(coeff_.numerator.x_coeff, coeff_.numerator.xy_coeff);
                swap(coeff_.denominator.const_coeff, coeff_.denominator.y_coeff);
                swap(coeff_.denominator.x_coeff, coeff_.denominator.xy_coeff);
            }
        }

    public:
        constexpr binary_gosper_continued_fraction(ContinuedFractionImplX x_cf,
                                                   ContinuedFractionImplY y_cf,
                                                   binary_gosper_coeff<int_type> coeff)
            : x_cf_{static_cast<ContinuedFractionImplX>(x_cf)},
              y_cf_{static_cast<ContinuedFractionImplY>(y_cf)},
              coeff_{static_cast<binary_gosper_coeff<int_type>&&>(coeff)} {}

        constexpr next_partial_fraction_return<partial_fraction_type> next_partial_fraction() {
            while (true) {
                if (is_x_terminated_ && is_y_terminated_) {
                    // Proceed as in the case of usual rational continued fractions.
                    util::constexpr_assert<util::error_msgs::divide_by_zero>(
                        !is_zero(coeff_.denominator.const_coeff));
                    util::constexpr_assert<util::error_msgs::divide_by_zero>(
                        is_strictly_positive(coeff_.denominator.const_coeff));

                    using std::div;
                    auto div_result =
                        div(coeff_.numerator.const_coeff, abs(coeff_.denominator.const_coeff));
                    coeff_.numerator.const_coeff =
                        static_cast<int_type&&>(coeff_.denominator.const_coeff);
                    coeff_.denominator.const_coeff =
                        int_type(static_cast<uint_type&&>(div_result.rem));

                    // Terminate if all coefficients in the denominator has become zero.
                    return {partial_fraction_type{{}, static_cast<int_type&&>(div_result.quot)},
                            is_zero(coeff_.denominator.const_coeff)};
                }
                else if (is_x_terminated_) {
                    if (!is_zero(coeff_.denominator.const_coeff) &&
                        !is_zero(coeff_.denominator.y_coeff)) {
                        auto const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto y_output =
                            div_floor(coeff_.numerator.y_coeff, coeff_.denominator.y_coeff);

                        if (const_output == y_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.y_coeff -= const_output * coeff_.denominator.y_coeff;

                            using util::swap;
                            swap(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);

                            // Terminate if all coefficients in the denominator has become zero.
                            return {
                                partial_fraction_type{{}, static_cast<int_type&&>(const_output)},
                                is_zero(coeff_.denominator.const_coeff) &&
                                    is_zero(coeff_.denominator.y_coeff)};
                        }
                    }

                    // If two endpoints do not agree, then refine the region.
                    progress_y();
                }
                else if (is_y_terminated_) {
                    if (!is_zero(coeff_.denominator.const_coeff) &&
                        !is_zero(coeff_.denominator.x_coeff)) {
                        auto const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto x_output =
                            div_floor(coeff_.numerator.x_coeff, coeff_.denominator.x_coeff);

                        if (const_output == x_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.x_coeff -= const_output * coeff_.denominator.x_coeff;

                            using util::swap;
                            swap(coeff_.numerator, coeff_.denominator);

                            // Terminate if all coefficients in the denominator has become zero.
                            return {
                                partial_fraction_type{{}, static_cast<int_type&&>(const_output)},
                                is_zero(coeff_.denominator.const_coeff) &&
                                    is_zero(coeff_.denominator.x_coeff)};
                        }
                    }

                    // If two endpoints do not agree, then refine the region.
                    progress_x();
                }
                else {
                    if (!is_zero(coeff_.denominator.const_coeff) &&
                        !is_zero(coeff_.denominator.x_coeff) &&
                        !is_zero(coeff_.denominator.y_coeff) &&
                        !is_zero(coeff_.denominator.xy_coeff)) {
                        auto const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto x_output =
                            div_floor(coeff_.numerator.x_coeff, coeff_.denominator.x_coeff);
                        auto y_output =
                            div_floor(coeff_.numerator.y_coeff, coeff_.denominator.y_coeff);
                        auto xy_output =
                            div_floor(coeff_.numerator.xy_coeff, coeff_.denominator.xy_coeff);

                        if (const_output == x_output && x_output == y_output &&
                            y_output == xy_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.x_coeff -= const_output * coeff_.denominator.x_coeff;
                            coeff_.numerator.y_coeff -= const_output * coeff_.denominator.y_coeff;
                            coeff_.numerator.xy_coeff -= const_output * coeff_.denominator.xy_coeff;

                            using util::swap;
                            swap(coeff_.numerator, coeff_.denominator);

                            // Terminate if all coefficients in the denominator has become zero.
                            return {
                                partial_fraction_type{{}, static_cast<int_type&&>(const_output)},
                                is_zero(coeff_.denominator.const_coeff) &&
                                    is_zero(coeff_.denominator.x_coeff) &&
                                    is_zero(coeff_.denominator.y_coeff) &&
                                    is_zero(coeff_.denominator.xy_coeff)};
                        }
                    }

                    // If four endpoints do not agree, then refine the region.
                    // Compare |b/f - a/e| and |c/g - a/e|, and choose x if the former is larger.
                    // I.e., we choose x if |bge - afg| >= |cfe - afg|.

                    // If e & f are both zero, progressing x will never make them nonzero, so
                    // progress y in this case.
                    if (is_zero(coeff_.denominator.const_coeff) &&
                        is_zero(coeff_.denominator.x_coeff)) {
                        progress_y();
                    }
                    // If e & g are both zero, progressing y will never make them nonzero, so
                    // progress x in this case.
                    else if (is_zero(coeff_.denominator.const_coeff) &&
                             is_zero(coeff_.denominator.y_coeff)) {
                        progress_x();
                    }
                    else {
                        // afg
                        auto const base_point = coeff_.numerator.const_coeff *
                                                coeff_.denominator.x_coeff *
                                                coeff_.denominator.y_coeff;
                        // |bge - afg|
                        auto const const_to_x_length =
                            abs(coeff_.numerator.x_coeff * coeff_.denominator.y_coeff *
                                    coeff_.denominator.const_coeff -
                                base_point);
                        // |bge - afg|
                        auto const const_to_y_length =
                            abs(coeff_.numerator.y_coeff * coeff_.denominator.x_coeff *
                                    coeff_.denominator.const_coeff -
                                base_point);

                        if (const_to_x_length >= const_to_y_length) {
                            progress_x();
                        }
                        else {
                            progress_y();
                        }
                    }
                }
            }
        }
    };
}

#endif
