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

#include "bigint.h"
#include "rational_continued_fractions.h"
#include "best_rational_approx.h"

namespace jkj {
    namespace idiv {
        enum class strategy { multiply_shift, multiply_shift_large, multiply_add_shift };

        struct multiply_shift_info {
            bigint::uint_var multiplier;
            std::size_t shift_amount;
        };

        // Precondition: x is in its reduced form.
        constexpr inline multiply_shift_info convert_to_multiply_shift_effectively_rational(
            frac<bigint::uint_var, bigint::uint_var> const& x, bigint::uint_var const& nmax) {
            util::constexpr_assert<util::error_msgs::no_error_msg>(x.denominator <= nmax);

            multiply_shift_info ret_value{};

            bigint::uint_var v;
            if (x.denominator != 1u) {
                // Compute the modular inverse of -x.numerator.
                auto const mod_inv =
                    find_best_rational_approx<
                        rational_continued_fractions<bigint::uint_var, bigint::uint_var>>(
                        x, x.denominator - 1u)
                        .above.denominator;

                // v = floor((nmax - b) / q) * q + b.
                v = ((nmax - mod_inv) / x.denominator) * x.denominator;
                v += mod_inv;
            }
            else {
                v = nmax;
            }

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
                auto left_end_plus_1 = ret_value.multiplier + 1u;

                // If the left_end_plus_1 is still in the interval, take that instead.
                if (left_end_plus_1 * v < (((v * x.numerator + 1u) / x.denominator) << k0)) {
                    ret_value.multiplier = static_cast<bigint::uint_var&&>(left_end_plus_1);
                    ret_value.shift_amount -= ret_value.multiplier.factor_out_power_of_2();
                }
            }

            return ret_value;
        }

        struct multiply_add_shift_info {
            bool succeeded = false;
            bigint::uint_var multiplier = {};
            bigint::uint_var adder = {};
            unsigned int shift_amount = 0;
        };

        constexpr inline multiply_add_shift_info convert_to_multiply_add_shift_effectively_rational(
            frac<bigint::uint_var, bigint::uint_var> const& x, bigint::uint_var const& nmax,
            bigint::uint_var const& max_allowed_value) {

            util::constexpr_assert<util::error_msgs::divide_by_zero>(x.denominator != 0);
            util::constexpr_assert<util::error_msgs::no_error_msg>(x.denominator <= nmax);

            multiply_add_shift_info ret_value;

            bigint::uint_var n_L0, n_U0;
            if (x.denominator != 1u) {
                // Find the largest multiple of x.denominator <= nmax.
                n_L0 = (nmax / x.denominator) * x.denominator;

                // Compute the modular inverse of -x.numerator.
                auto const mod_inv =
                    find_best_rational_approx<
                        rational_continued_fractions<bigint::uint_var, bigint::uint_var>>(
                        x, x.denominator - 1u)
                        .above.denominator;

                // v = floor((nmax - b) / q) * q + b.
                n_U0 = ((nmax - mod_inv) / x.denominator) * x.denominator;
                n_U0 += mod_inv;
            }
            else {
                n_L0 = nmax;
                n_U0 = nmax;
            }

            using ufrac = frac<bigint::uint_var, bigint::uint_var>;

            ufrac zeta_max{0, 1}, zeta_Lmax{0, 1}, zeta_Umax{0, 1}, zeta_min;
            bigint::uint_var n_L1 = 0, n_U1 = 0;
            bigint::uint_var floor_n_L0_x, floor_n_U0_x_p1;

            while (zeta_max.numerator != zeta_max.denominator) {
                multiply_add_shift_info candidate;

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
                            find_best_rational_approx<
                                rational_continued_fractions<bigint::uint_var, bigint::uint_var>>(
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
                    floor_n_U0_x_p1 = ((n_U0 * x.numerator) / x.denominator) + 1u;

                    if (n_U0 == 1) {
                        zeta_Umax = {1, 1};
                    }
                    else {
                        auto const new_nmax = n_U0 - 1u;
                        auto const best_approx =
                            find_best_rational_approx<
                                rational_continued_fractions<bigint::uint_var, bigint::uint_var>>(
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

                {
                    auto delta = ufrac{right_end.numerator * left_end.denominator,
                                       left_end.denominator * right_end.denominator};
                    auto numerator_diff = right_end.denominator * left_end.numerator;

                    // If the interval is empty, move to the next subinterval for zeta.
                    if (delta.numerator <= numerator_diff) {
                        continue;
                    }

                    delta.numerator -= numerator_diff;
                    util::constexpr_assert<util::error_msgs::no_error_msg>(delta.denominator >=
                                                                           delta.numerator);
                    candidate.shift_amount =
                        trunc_floor_log2_div(delta.denominator, delta.numerator);
                }

                candidate.multiplier =
                    ((left_end.numerator << candidate.shift_amount) / left_end.denominator) + 1u;

                // If t goes out of the interval, then increase k0.
                if (candidate.multiplier * right_end.denominator >=
                    (right_end.numerator << candidate.shift_amount)) {
                    ++candidate.shift_amount;
                    candidate.multiplier =
                        ((left_end.numerator << candidate.shift_amount) / left_end.denominator) +
                        1u;
                }
                else {
                    candidate.shift_amount -= candidate.multiplier.factor_out_power_of_2();
                }

                // Truncate zeta0 from 0 to avoid underflow.
                auto zeta0 = ufrac{(floor_n_L0_x << candidate.shift_amount),
                                   bigint::uint_var::power_of_2(candidate.shift_amount)};
                {
                    auto numerator_diff = n_L0 * candidate.multiplier;
                    if (zeta0.numerator > numerator_diff) {
                        zeta0.numerator -= numerator_diff;
                    }
                    else {
                        zeta0.numerator = 0;
                    }
                }

                while (candidate.multiplier * nmax + zeta0.numerator <= max_allowed_value) {
                    // Loop over all numerators in the interval.
                    while (true) {
                        // These branches do not touch candidate.multiplier and
                        // candidate.shift_amount, but may modify candidate.adder.
                        if (zeta0 >= zeta_min) {
                            candidate.adder = zeta0.numerator;

                            // Check admissibility of xi with respect to zeta.
                            if (candidate.multiplier * n_U0 + candidate.adder <
                                (floor_n_U0_x_p1 << candidate.shift_amount)) {
                                // Found.
                                candidate.succeeded = true;
                                break;
                            }
                        }
                        else {
                            auto const delta_zeta = zeta_min - zeta0;
                            candidate.adder =
                                zeta0.numerator +
                                div_ceil((delta_zeta.numerator << candidate.shift_amount),
                                         delta_zeta.denominator);

                            // Check zeta < zeta_max.
                            if (candidate.adder * zeta_max.denominator <
                                (zeta_max.numerator << candidate.shift_amount)) {
                                // Check the max_allowed_value constraint.
                                if (candidate.multiplier * nmax + candidate.adder <=
                                    max_allowed_value) {
                                    // Check admissibility of xi with respect to zeta.
                                    if (candidate.multiplier * n_U0 + candidate.adder <
                                        (floor_n_U0_x_p1 << candidate.shift_amount)) {
                                        // Found.
                                        candidate.succeeded = true;
                                        break;
                                    }
                                }
                            }
                        }

                        // Try the next numerator.
                        ++candidate.multiplier;
                        if (candidate.multiplier * right_end.denominator >=
                            (right_end.numerator << candidate.shift_amount)) {
                            break;
                        }

                        if (zeta0.numerator > n_L0) {
                            zeta0.numerator -= n_L0;
                        }
                        else {
                            zeta0.numerator = 0;
                        }

                        if (candidate.multiplier * nmax + zeta0.numerator > max_allowed_value) {
                            break;
                        }
                    }

                    if (candidate.succeeded) {
                        // If this is the first success, record it and move to the next
                        // subinterval.
                        if (!ret_value.succeeded) {
                            ret_value = candidate;
                        }
                        // Otherwise, compare it with the previous best one and replace if
                        // appropriate.
                        else {
                            if (ret_value.shift_amount > candidate.shift_amount) {
                                ret_value = candidate;
                            }
                            else if (ret_value.shift_amount == candidate.shift_amount) {
                                if (ret_value.multiplier > candidate.multiplier) {
                                    ret_value = candidate;
                                }
                                else if (ret_value.multiplier == candidate.multiplier) {
                                    if (ret_value.adder > candidate.adder) {
                                        ret_value = candidate;
                                    }
                                }
                            }
                        }

                        break;
                    }

                    // Increase k0 and recompute t, zeta0.
                    ++candidate.shift_amount;
                    candidate.multiplier =
                        ((left_end.numerator << candidate.shift_amount) / left_end.denominator) +
                        1u;

                    zeta0.numerator = (floor_n_L0_x << candidate.shift_amount);
                    zeta0.denominator <<= 1;
                    {
                        auto numerator_diff = n_L0 * candidate.multiplier;
                        if (zeta0.numerator > numerator_diff) {
                            zeta0.numerator -= numerator_diff;
                        }
                        else {
                            zeta0.numerator = 0;
                        }
                    }
                }
            }

            return ret_value;
        }
    }
}

#endif