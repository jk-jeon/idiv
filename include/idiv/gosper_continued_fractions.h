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

#ifndef JKJ_HEADER_GOSPER_ALGORITHM
#define JKJ_HEADER_GOSPER_ALGORITHM

#include "continued_fractions.h"
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

            constexpr void feed_x(Int const& continued_fractions_coeff) {
                // (a, b, c, d) |-> (b, a + tb, d, c + td)
                const_coeff += continued_fractions_coeff * x_coeff;
                y_coeff += continued_fractions_coeff * xy_coeff;
                swap(const_coeff, x_coeff);
                swap(y_coeff, xy_coeff);
            }

            constexpr void feed_y(Int const& continued_fractions_coeff) {
                // (a, b, c, d) |-> (c, d, a + tc, b + td)
                const_coeff += continued_fractions_coeff * y_coeff;
                x_coeff += continued_fractions_coeff * xy_coeff;
                swap(const_coeff, y_coeff);
                swap(x_coeff, xy_coeff);
            }
        };

        coeff numerator;
        coeff denominator;
    };

    template <class ContinuedFractionsCalcX, class ContinuedFractionsCalcY>
        requires(std::is_same_v<typename ContinuedFractionsCalcX::int_type,
                                typename ContinuedFractionsCalcY::int_type> &&
                 std::is_same_v<typename ContinuedFractionsCalcX::uint_type,
                                typename ContinuedFractionsCalcY::uint_type>)
    class gosper_continued_fractions
        : public continued_fractions<
              gosper_continued_fractions<ContinuedFractionsCalcX, ContinuedFractionsCalcY>,
              typename ContinuedFractionsCalcX::int_type,
              typename ContinuedFractionsCalcX::uint_type> {
    public:
        using int_type = typename ContinuedFractionsCalcX::int_type;
        using uint_type = typename ContinuedFractionsCalcX::uint_type;

    private:
        using crtp_base = continued_fractions<
            gosper_continued_fractions<ContinuedFractionsCalcX, ContinuedFractionsCalcY>, int_type,
            uint_type>;
        friend crtp_base;

        ContinuedFractionsCalcX x_cf_;
        ContinuedFractionsCalcY y_cf_;
        gosper_coeff<int_type> coeff_;

        constexpr void progress_x() {
            x_cf_.update();
            coeff_.numerator.feed_x(x_cf_.current_coefficient());
            coeff_.denominator.feed_x(x_cf_.current_coefficient());

            // If the current coefficient is the last one, we have to use a different update
            // formula from now on.
            if (x_cf_.is_terminated()) {
                swap(coeff_.numerator.const_coeff, coeff_.numerator.x_coeff);
                swap(coeff_.numerator.y_coeff, coeff_.numerator.xy_coeff);
            }
        }
        constexpr void progress_y() {
            y_cf_.update();
            coeff_.numerator.feed_y(y_cf_.current_coefficient());
            coeff_.denominator.feed_y(y_cf_.current_coefficient());

            // If the current coefficient is the last one, we have to use a different update
            // formula from now on.
            if (y_cf_.is_terminated()) {
                swap(coeff_.numerator.const_coeff, coeff_.numerator.y_coeff);
                swap(coeff_.numerator.x_coeff, coeff_.numerator.xy_coeff);
            }
        }

        constexpr int_type compute_next_coefficient() {
            util::constexpr_assert<util::error_msgs::no_error_msg>(!x_cf_.is_terminated() ||
                                                                   !y_cf_.is_terminated());
            while (true) {
                if (x_cf_.is_terminated()) {
                    if (!is_zero(coeff_.denominator.const_coeff) &&
                        !is_zero(coeff_.denominator.y_coeff)) {
                        auto const const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto const y_output =
                            div_floor(coeff_.numerator.y_coeff, coeff_.denominator.y_coeff);

                        if (const_output == y_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.y_coeff -= const_output * coeff_.denominator.y_coeff;

                            swap(coeff_.numerator, coeff_.denominator);

                            // Terminate if all coefficients in the denominator has become zero.
                            if (is_zero(coeff_.denominator.const_coeff) &&
                                is_zero(coeff_.denominator.y_coeff)) {
                                crtp_base::set_terminate_flag();
                            }

                            return const_output;
                        }
                    }

                    // If two endpoints do not agree, then refine the region.
                    progress_y();
                }
                else if (y_cf_.is_terminated()) {
                    if (!is_zero(coeff_.denominator.const_coeff) &&
                        !is_zero(coeff_.denominator.x_coeff)) {
                        auto const const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto const x_output =
                            div_floor(coeff_.numerator.x_coeff, coeff_.denominator.x_coeff);

                        if (const_output == x_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.x_coeff -= const_output * coeff_.denominator.x_coeff;

                            swap(coeff_.numerator, coeff_.denominator);

                            // Terminate if all coefficients in the denominator has become zero.
                            if (is_zero(coeff_.denominator.const_coeff) &&
                                is_zero(coeff_.denominator.x_coeff)) {
                                crtp_base::set_terminate_flag();
                            }

                            return const_output;
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
                        auto const const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto const x_output =
                            div_floor(coeff_.numerator.x_coeff, coeff_.denominator.x_coeff);
                        auto const y_output =
                            div_floor(coeff_.numerator.y_coeff, coeff_.denominator.y_coeff);
                        auto const xy_output =
                            div_floor(coeff_.numerator.xy_coeff, coeff_.denominator.xy_coeff);

                        if (const_output == x_output && x_output == y_output &&
                            y_output == xy_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.x_coeff -= const_output * coeff_.denominator.x_coeff;
                            coeff_.numerator.y_coeff -= const_output * coeff_.denominator.y_coeff;
                            coeff_.numerator.xy_coeff -= const_output * coeff_.denominator.xy_coeff;

                            swap(coeff_.numerator, coeff_.denominator);

                            // Terminate if all coefficients in the denominator has become zero.
                            if (is_zero(coeff_.denominator.const_coeff) &&
                                is_zero(coeff_.denominator.x_coeff) &&
                                is_zero(coeff_.denominator.y_coeff) &&
                                is_zero(coeff_.denominator.xy_coeff)) {
                                crtp_base::set_terminate_flag();
                            }

                            return const_output;
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

    public:
        constexpr gosper_continued_fractions(ContinuedFractionsCalcX x_cf,
                                             ContinuedFractionsCalcY y_cf,
                                             gosper_coeff<int_type> coeff)
            : x_cf_{static_cast<ContinuedFractionsCalcX>(x_cf)},
              y_cf_{static_cast<ContinuedFractionsCalcY>(y_cf)},
              coeff_{static_cast<gosper_coeff<int_type>&&>(coeff)} {}
    };
}

#endif