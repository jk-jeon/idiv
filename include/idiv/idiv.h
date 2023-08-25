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

#include "best_rational_approx.h"
#include "bigint.h"
#include "caching_continued_fractions.h"
#include "interval.h"
#include "rational_continued_fractions.h"

namespace jkj {
    namespace idiv {
        // For a given real number x and a positive integer nmax, find the interval
        // [max_n floor(nx)/n, min_n (floor(nx)+1)/n), where n ranges from {1, ... , nmax}.
        // The number x is given in terms of its continued fractions. The first parameter cf is the
        // continued fractions calculator for x. It must be initialized, i.e., it should start with
        // the "convergent" 1/0 when evaluated.
        template <class ContinuedFractionsCalc, class UInt>
        constexpr interval<frac<typename ContinuedFractionsCalc::int_type,
                                typename ContinuedFractionsCalc::uint_type>,
                           interval_type_t::bounded_left_closed_right_open>
        find_floor_quotient_range(ContinuedFractionsCalc& cf, UInt const& nmax) {
            util::constexpr_assert<util::error_msgs::no_error_msg>(is_strictly_positive(nmax));

            using frac_t = frac<typename ContinuedFractionsCalc::int_type,
                                typename ContinuedFractionsCalc::uint_type>;

            // First, find the last convergent and the last semiconvergent whose denominator is
            // bounded above by nmax.
            frac_t previous_previous_convergent;

            // This lambda replaces previous_previous_convergent to the last semiconvergent, and
            // return it as an rvalue reference.
            auto get_last_semiconvergent = [&]() -> decltype(auto) {
                auto semiconvergent_coeff =
                    (nmax - previous_previous_convergent.denominator) / cf.previous_denominator();

                previous_previous_convergent.numerator +=
                    semiconvergent_coeff * cf.previous_numerator();
                previous_previous_convergent.denominator +=
                    semiconvergent_coeff * cf.previous_denominator();

                return static_cast<frac_t&&>(previous_previous_convergent);
            };

            while (cf.current_denominator() <= nmax) {
                // Obtain the next convergent.
                previous_previous_convergent = cf.previous_convergent();
                cf.update();

                // If we reach to the perfect approximation, then we have to find the largest
                // positive integer v <= nmax such that vp == -1 (mod q), where x = p/q. Then the
                // lower bound is p/q, while the upper bound is ((vp+1)/q) / v.
                if (cf.is_terminated()) {
                    // To compute v, we find the modular inverse b of -p. This can be done by
                    // observing that the best rational approximation from above whose denominator
                    // is strictly less than q must be precisely ((bp+1)/q) / b. Then
                    //        v = b + floor((nmax - b)/q)q and
                    // (vp+1)/q = (bp+1)/q + floor((nmax - b)/q)p.

                    // If we ended at an even convergent, the last convergent is the best rational
                    // approximation from above. Otherwise, the last semiconvergent is the best
                    // rational approximation from above.
                    auto upper_bound = cf.current_index() % 2 == 0 ? cf.previous_convergent()
                                                                   : get_last_semiconvergent();

                    // At this point, upper_bound is ((bp+1)/q) / b, so we adjust the numerator and
                    // the denominator by floor((nmax - b)/q).
                    auto max_quotient =
                        (nmax - upper_bound.denominator) / cf.current_convergent().denominator;
                    upper_bound.numerator += max_quotient * cf.current_convergent().numerator;
                    upper_bound.denominator += max_quotient * cf.current_convergent().denominator;

                    return {cf.current_convergent(), static_cast<frac_t&&>(upper_bound)};
                }
            }

            // If there the last convergent is still not a perfect approximation, then we return
            // the last semiconvergent and the convergent as the bounds. Which one is the lower
            // bound and which one is the upper bound is determined by the parity of
            // cf.current_index(). Note that cf.current_index() is the index of the first
            // convergent with the denominator strictly larger than nmax, so if this index is
            // even, then the semiconvergent is the lower bound and the convergent is the upepr
            // bound, and if the index is odd, then the other way around.
            if (cf.current_index() % 2 == 0) {
                return {get_last_semiconvergent(), cf.previous_convergent()};
            }
            else {
                return {cf.previous_convergent(), get_last_semiconvergent()};
            }
        }


        struct multiply_shift_info {
            bigint::int_var multiplier;
            std::size_t shift_amount;
        };

        // For a given real number x and a positive integer nmax, find the smallest nonnegative
        // integer k such that there exists an integer m satisfying
        // floor(nx) = floor(nm/2^k) for all n = 1, ... , nmax.
        // The number x is given in terms of its continued fractions. The first parameter cf is the
        // continued fractions calculator for x. It must be initialized, i.e., it should start with
        // the "convergent" 1/0 when evaluated.
        template <class ContinuedFractionsCalc>
        constexpr multiply_shift_info convert_to_multiply_shift(ContinuedFractionsCalc& cf,
                                                                bigint::uint_var const& nmax) {
            multiply_shift_info ret_value{};
            auto range = find_floor_quotient_range(cf, nmax);

            // k0 = ceil(log2(1/Delta)).
            auto k0 = [&] {
                auto const delta = range.upper_bound() - range.lower_bound();
                util::constexpr_assert<util::error_msgs::no_error_msg>(
                    is_strictly_positive(delta.numerator));

                return trunc_ceil_log2_div(delta.denominator, abs(delta.numerator));
            }();

            ret_value.multiplier =
                div_ceil((range.lower_bound().numerator << k0), range.lower_bound().denominator);
            ret_value.shift_amount = k0;

            if (ret_value.multiplier.is_even()) {
                ret_value.shift_amount -= ret_value.multiplier.factor_out_power_of_2();
            }
            else {
                auto left_end_plus_1 = ret_value.multiplier + 1u;

                // If the left_end_plus_1 is still in the interval, take that instead.
                if (left_end_plus_1 * range.upper_bound().denominator <
                    (range.upper_bound().numerator << k0)) {
                    ret_value.multiplier =
                        static_cast<decltype(left_end_plus_1)&&>(left_end_plus_1);
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
            util::constexpr_assert<util::error_msgs::divide_by_zero>(!x.denominator.is_zero());
            util::constexpr_assert<util::error_msgs::no_error_msg>(x.denominator <= nmax);

            using continued_fractions_calculator_type = caching_continued_fractions<
                rational_continued_fractions<bigint::uint_var, bigint::uint_var>, std::vector>;

            multiply_add_shift_info ret_value;
            continued_fractions_calculator_type continued_fractions_calculator{x};

            bigint::uint_var n_L0, n_U0;
            if (x.denominator != 1u) {
                // Find the largest multiple of x.denominator <= nmax.
                n_L0 = (nmax / x.denominator) * x.denominator;

                // Compute the modular inverse of -x.numerator.
                auto const mod_inv =
                    find_best_rational_approx(continued_fractions_calculator, x.denominator - 1u)
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
                        continued_fractions_calculator.rewind();
                        auto const best_approx =
                            find_best_rational_approx(continued_fractions_calculator, new_nmax)
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
                        continued_fractions_calculator.rewind();
                        auto const best_approx =
                            find_best_rational_approx(continued_fractions_calculator, new_nmax)
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