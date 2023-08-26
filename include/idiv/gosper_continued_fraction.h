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
#include <cstdlib>
#include <type_traits>

// Implements Gosper's algorithm.
// Gosper's algorithm computes the continued fractions of the number
// z = (a + bx + cy + dxy) / (e + fx + gy + hxy)
// in terms of the continued fractions of x and y.

namespace jkj {
    template <class Int>
    struct gosper_coeff {
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

    template <class ContinuedFractionImplX, class ContinuedFractionImplY, class Unity = unity,
              template <class> class... Mixin>
        requires(std::is_same_v<typename ContinuedFractionImplX::partial_fraction_type,
                                typename ContinuedFractionImplY::partial_fraction_type> &&
                 std::is_same_v<typename ContinuedFractionImplX::convergent_type,
                                typename ContinuedFractionImplY::convergent_type>)
    class gosper_continued_fraction {
    public:
        using int_type = decltype(ContinuedFractionImplX::convergent_type::numerator);
        using uint_type = decltype(ContinuedFractionImplX::convergent_type::denominator);
        using partial_fraction_type = frac<Unity, int_type>;
        using convergent_type = typename ContinuedFractionImplX::convergent_type;

    private:
        ContinuedFractionImplX x_cf_;
        ContinuedFractionImplY y_cf_;
        gosper_coeff<int_type> coeff_;
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
        constexpr gosper_continued_fraction(ContinuedFractionImplX x_cf,
                                            ContinuedFractionImplY y_cf,
                                            gosper_coeff<int_type> coeff)
            : x_cf_{static_cast<ContinuedFractionImplX>(x_cf)},
              y_cf_{static_cast<ContinuedFractionImplY>(y_cf)},
              coeff_{static_cast<gosper_coeff<int_type>&&>(coeff)} {}

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
                            swap(coeff_.numerator, coeff_.denominator);

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
                    // afg
                    auto const base_point = coeff_.numerator.const_coeff *
                                            coeff_.denominator.x_coeff * coeff_.denominator.y_coeff;
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
    };
}

#endif
