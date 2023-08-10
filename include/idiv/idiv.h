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


#ifndef JKJ_HEADER_IDIV
#define JKJ_HEADER_IDIV

#include "big_uint.h"
#include "rational_continued_fractions.h"
#include "best_rational_approx.h"
#include <optional>

namespace jkj {
    namespace idiv {
        enum class strategy { multiply_shift, multiply_shift_large, multiply_add_shift };

        struct multiply_shift_info {
            big_uint::var multiplier;
            std::size_t shift_amount;
        };

        // Precondition: x is in its reduced form.
        constexpr inline multiply_shift_info
        convert_to_multiply_shift_effectively_rational(frac<big_uint::var, big_uint::var> const& x,
                                                       big_uint::var const& nmax) {
            util::constexpr_assert<util::error_msgs::no_error_msg>(x.denominator <= nmax);

            multiply_shift_info ret_value{};

            // Compute the modular inverse of -x.numerator.
            auto const mod_inv =
                find_best_rational_approx<rational_continued_fractions<big_uint::var>>(
                    x, x.denominator - 1)
                    .above.denominator;

            // v = floor((nmax - b) / q) * q + b.
            auto v = ((nmax - mod_inv) / x.denominator) * x.denominator;
            v += mod_inv;

            auto const reciprocal_interval_length = v * x.denominator;

            // k0 = ceil(log2(1/Delta)).
            auto k0 = bit_width(reciprocal_interval_length);
            if (reciprocal_interval_length.is_power_of_2()) {
                --k0;
            }

            ret_value.multiplier = div_ceil((x.numerator << k0), x.denominator);
            ret_value.shift_amount = k0;

            if (ret_value.multiplier.is_even()) {
                ret_value.shift_amount -= ret_value.multiplier.factor_out_power_of_2();
            }
            else {
                auto left_end_plus_1 = ret_value.multiplier + 1;

                // If the left_end_plus_1 is still in the interval, take that instead.
                if (left_end_plus_1 * v < (((v * x.numerator + 1) / x.denominator) << k0)) {
                    ret_value.multiplier = static_cast<big_uint::var&&>(left_end_plus_1);
                    ret_value.shift_amount -= ret_value.multiplier.factor_out_power_of_2();
                }
            }

            return ret_value;
        }

        struct multiply_add_shift_info {
            big_uint::var multiplier;
            big_uint::var adder;
            unsigned int shift_amount;
        };

        constexpr inline std::optional<multiply_add_shift_info>
        convert_to_multiply_add_shift_effectively_rational(
            frac<big_uint::var, big_uint::var> const& x, big_uint::var const& nmax,
            big_uint::var const& max_allowed_value) {

            util::constexpr_assert<util::error_msgs::divide_by_zero>(x.denominator != 0);
            util::constexpr_assert<util::error_msgs::no_error_msg>(x.denominator <= nmax);

            multiply_add_shift_info ret_value{};

            big_uint::var n_L0, n_U0;
            if (x.denominator != 1) {
                // Find the largest multiple of x.denominator <= nmax.
                n_L0 = (nmax / x.denominator) * x.denominator;

                // Compute the modular inverse of -x.numerator.
                auto const mod_inv =
                    find_best_rational_approx<rational_continued_fractions<big_uint::var>>(
                        x, x.denominator - 1)
                        .above.denominator;

                // v = floor((nmax - b) / q) * q + b.
                n_U0 = ((nmax - mod_inv) / x.denominator) * x.denominator;
                n_U0 += mod_inv;
            }
            else {
                n_L0 = nmax;
                n_U0 = 1;
            }

            using ufrac = frac<big_uint::var, big_uint::var>;

            ufrac zeta_max{0, 1}, zeta_Lmax{0, 1}, zeta_Umax{0, 1}, zeta_min;
            big_uint::var n_L1 = 0, n_U1 = 0;
            big_uint::var floor_n_L0_x, floor_n_U0_x_p1;

            while (zeta_max.numerator != zeta_max.denominator) {
                // Update zeta_Lmax if necessary.
                if (zeta_max == zeta_Lmax) {
                    n_L0 += n_L1;
                    floor_n_L0_x = (n_L0 * x.numerator) / x.denominator;

                    if (n_L0 == nmax) {
                        zeta_Lmax = {1, 1};
                    }
                    else {
                        auto const new_nmax = nmax - n_L0;
                        auto const best_approx =
                            find_best_rational_approx<rational_continued_fractions<big_uint::var>>(
                                x, new_nmax)
                                .below;
                        auto const largest_multiplier = new_nmax / best_approx.denominator;
                        n_L1 = largest_multiplier * best_approx.denominator;

                        zeta_Lmax = {(n_L1 * floor_n_L0_x) -
                                         n_L0 * (largest_multiplier * best_approx.numerator),
                                     n_L1};

                        // Truncate to 1 if necessary.
                        if (zeta_Lmax.numerator > zeta_Lmax.denominator) {
                            zeta_Lmax = {1, 1};
                        }
                    }
                }

                // Update zeta_Umax if necessary.
                if (zeta_max == zeta_Umax) {
                    n_U0 -= n_U1;
                    floor_n_U0_x_p1 = ((n_U0 * x.numerator) / x.denominator) + 1;

                    if (n_U0 == 1) {
                        zeta_Umax = {1, 1};
                    }
                    else {
                        auto const new_nmax = n_U0 - 1;
                        auto const best_approx =
                            find_best_rational_approx<rational_continued_fractions<big_uint::var>>(
                                x, new_nmax)
                                .below;
                        auto const largest_multiplier = new_nmax / best_approx.denominator;
                        n_U1 = largest_multiplier * best_approx.denominator;

                        zeta_Umax = {(n_U1 * floor_n_U0_x_p1) -
                                         n_U0 * (largest_multiplier * best_approx.numerator),
                                     n_U1};

                        // Truncate to 1 if necessary.
                        if (zeta_Umax.numerator > zeta_Umax.denominator) {
                            zeta_Umax = {1, 1};
                        }
                    }
                }

                zeta_min = static_cast<ufrac&&>(zeta_max);
                zeta_max = zeta_Lmax < zeta_Umax ? zeta_Lmax : zeta_Umax;

                auto const left_end =
                    ufrac{floor_n_L0_x * zeta_max.denominator - zeta_max.numerator,
                          n_L0 * zeta_max.denominator};
                auto const right_end =
                    ufrac{floor_n_U0_x_p1 * zeta_min.denominator - zeta_min.numerator,
                          n_U0 * zeta_min.denominator};

                std::size_t k0;
                {
                    auto delta = ufrac{right_end.numerator * left_end.denominator,
                                       left_end.denominator * right_end.denominator};
                    auto numerator_diff = right_end.denominator * left_end.numerator;

                    // If the interval is empty, move to the next subinterval for zeta.
                    if (delta.numerator <= numerator_diff) {
                        continue;
                    }

                    delta.numerator -= numerator_diff;
                    k0 = trunc_floor_log2_div(delta.denominator, delta.numerator);
                }

                ret_value.multiplier = ((left_end.numerator << k0) / left_end.denominator) + 1;

                // If t goes out of the interval, then increase k0.
                if (ret_value.multiplier * right_end.denominator >= (right_end.numerator << k0)) {
                    ++k0;
                    ret_value.multiplier = ((left_end.numerator << k0) / left_end.denominator) + 1;
                    ret_value.shift_amount = k0;
                }
                else {
                    ret_value.shift_amount = k0 - ret_value.multiplier.factor_out_power_of_2();
                }

                // Truncate zeta0 from 0 to avoid underflow.
                auto zeta0 = ufrac{(floor_n_L0_x << ret_value.shift_amount),
                                   big_uint::var::power_of_2(ret_value.shift_amount)};
                {
                    auto numerator_diff = n_L0 * ret_value.multiplier;
                    if (zeta0.numerator > numerator_diff) {
                        zeta0.numerator -= numerator_diff;
                    }
                    else {
                        zeta0.numerator = 0;
                    }
                }

                bool initial_k0 = true;
                while (ret_value.multiplier * nmax + zeta0.numerator <= max_allowed_value) {
                    // Loop over all odd numerators in the interval.
                    while (true) {
                        // These branches do not touch ret_value.multiplier and
                        // ret_value.shift_amount (unless it returns), but may modify
                        // ret_value.adder.
                        if (zeta0 >= zeta_min) {
                            ret_value.adder = zeta0.numerator;

                            // Check admissibility of xi with respect to zeta.
                            if (ret_value.multiplier * n_U0 + ret_value.adder <
                                (floor_n_U0_x_p1 << ret_value.shift_amount)) {
                                // Found.
                                return ret_value;
                            }
                        }
                        else {
                            auto const delta_zeta = zeta_min - zeta0;
                            auto zeta0_numerator_shifted = zeta0.numerator;
                            auto m = ret_value.multiplier;
                            auto k = ret_value.shift_amount;
                            do {
                                ret_value.adder =
                                    zeta0_numerator_shifted +
                                    div_ceil((delta_zeta.numerator << k), delta_zeta.denominator);

                                // Check zeta < zeta_max.
                                if (ret_value.adder * zeta_max.denominator <
                                    (zeta_max.numerator << k)) {
                                    // Check the max_allowed_value constraint.
                                    if (m * nmax + ret_value.adder <= max_allowed_value) {
                                        // Check admissibility of xi with respect to zeta.
                                        if (m * n_U0 + ret_value.adder < (floor_n_U0_x_p1 << k)) {
                                            // Found.
                                            ret_value.multiplier = static_cast<big_uint::var&&>(m);
                                            ret_value.shift_amount = k;
                                            return ret_value;
                                        }
                                    }
                                }

                                ++k;
                                m <<= 1;
                                zeta0_numerator_shifted <<= 1;

                                // Stop the iteration if the max_allowed_value constraint can no
                                // longer be satisfied.
                            } while (m * nmax + zeta0_numerator_shifted <= max_allowed_value);
                        }

                        // The first interval is guaranteed to have only one candidate.
                        if (initial_k0) {
                            initial_k0 = false;
                            break;
                        }

                        // Try the next odd number.
                        ret_value.multiplier += 2;
                        if (ret_value.multiplier * right_end.denominator >=
                            (right_end.numerator << k0)) {
                            break;
                        }

                        auto const numerator_diff = (n_L0 << 1);
                        if (zeta0.numerator > numerator_diff) {
                            zeta0.numerator -= numerator_diff;
                        }
                        else {
                            zeta0.numerator = 0;
                        }
                    }

                    // Increase k0 and recompute t, zeta0.
                    ++k0;
                    ret_value.multiplier = ((left_end.numerator << k0) / left_end.denominator) + 1;
                    // If there are only even integers in the interval, double the interval and try
                    // again until we find an odd integer.
                    if (ret_value.multiplier.is_even()) {
                        ++ret_value.multiplier;
                        while (ret_value.multiplier * right_end.denominator >=
                               (right_end.numerator << k0)) {
                            ++k0;
                            ret_value.multiplier =
                                ((left_end.numerator << k0) / left_end.denominator) + 1;

                            if (ret_value.multiplier.is_even()) {
                                ++ret_value.multiplier;
                            }
                            else {
                                break;
                            }
                        }
                    }

                    ret_value.shift_amount = k0;
                    zeta0 = {(floor_n_L0_x << ret_value.shift_amount),
                             big_uint::var::power_of_2(ret_value.shift_amount)};
                    {
                        auto numerator_diff = n_L0 * ret_value.multiplier;
                        if (zeta0.numerator > numerator_diff) {
                            zeta0.numerator -= numerator_diff;
                        }
                        else {
                            zeta0.numerator = 0;
                        }
                    }
                }
            }

            // FAIL.
            return {};
        }
    }
}

#endif