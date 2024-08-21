// Copyright 2023-2024 Junekey Jeon
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

#ifndef JKJ_HEADER_IDIV_FRACTIONAL_PART_EXTREMIZER
#define JKJ_HEADER_IDIV_FRACTIONAL_PART_EXTREMIZER

#include "continued_fraction/engine/rational.h"
#include "simultaneous_floor.h"

namespace jkj {
    namespace idiv {
        // Given real numbers x, y and a range [nmin:nmax] of integers, find the smallest minimizer
        // and the largest maximizer of (nx+y) - floor(nx+y).
        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY>
        constexpr extremizers_of_fractional_part<bigint::int_var>
        find_extremizers_of_fractional_part(
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            static_assert(
                cntfrc::has_mixins<ContinuedFractionGeneratorX, cntfrc::index_tracker,
                                   cntfrc::convergent_tracker, cntfrc::interval_estimate_provider,
                                   cntfrc::rewinder>(),
                "the first continued fraction generator must include index_tracker, "
                "convergent_tracker, interval_estimate_provider and rewinder");
            static_assert(
                cntfrc::has_mixins<ContinuedFractionGeneratorY, cntfrc::interval_estimate_provider,
                                   cntfrc::rewinder>(),
                "the second continued fraction generator must include interval_estimate_provider "
                "and rewinder");

            extremizers_of_fractional_part<bigint::int_var> result{nrange.lower_bound(),
                                                                   nrange.lower_bound()};

            if (nrange.lower_bound() == nrange.upper_bound()) {
                return result;
            }

            // First, find fine enough approximations of x and y.
            auto approx_info = find_simultaneous_multiply_add_shift(
                std::forward<ContinuedFractionGeneratorX>(xcf),
                std::forward<ContinuedFractionGeneratorY>(ycf), nrange);
            auto xi_cf = cntfrc::make_rational_generator<cntfrc::convergent_tracker>(
                approx_info.multiplier, bigint::uint_var::power_of_2(approx_info.shift_amount));

            // RHS times 2^k.
            auto compute_scaled_threshold_for_maximizer = [&] {
                auto temp = result.largest_maximizer * approx_info.multiplier + approx_info.adder;
                temp =
                    (((temp >> approx_info.shift_amount) + 1) << approx_info.shift_amount) - temp;
                return util::abs(std::move(temp));
            };
            auto scaled_threshold_for_maximizer = compute_scaled_threshold_for_maximizer();
            auto compute_scaled_threshold_for_minimizer = [&] {
                auto temp = result.smallest_minimizer * approx_info.multiplier + approx_info.adder;
                temp = temp - ((temp >> approx_info.shift_amount) << approx_info.shift_amount);
                return util::abs(std::move(temp));
            };
            auto scaled_threshold_for_minimizer = compute_scaled_threshold_for_minimizer();

            bool found_minimizer = false;
            bool found_maximizer = false;
            auto itr = xi_cf.begin();

            util::constexpr_assert(itr != xi_cf.end());

            // qx - p or p - qx.
            auto previous_scaled_lhs_for_convergent =
                util::abs((itr->previous_convergent_numerator() << approx_info.shift_amount) -
                          approx_info.multiplier * itr->previous_convergent_denominator());

            while (true) {
                // Maximizer.
                // Find a new even convergent.
                auto scaled_lhs_for_convergent =
                    util::abs(approx_info.multiplier * itr->current_convergent_denominator() -
                              (itr->current_convergent_numerator() << approx_info.shift_amount));

                if (!found_maximizer) {
                    auto increment = itr->current_convergent_denominator();

                    if (itr != xi_cf.end()) {
                        while (scaled_lhs_for_convergent < scaled_threshold_for_maximizer) {
                            auto scaled_lhs = scaled_lhs_for_convergent;
                            {
                                auto s = util::div_ceil(scaled_threshold_for_maximizer -
                                                            scaled_lhs_for_convergent,
                                                        previous_scaled_lhs_for_convergent) -
                                         1u;

                                scaled_lhs += s * previous_scaled_lhs_for_convergent;
                                increment -= std::move(s) * itr->previous_convergent_denominator();
                            }

                            if (!util::is_zero(scaled_lhs)) {
                                auto new_estimate =
                                    result.largest_maximizer +
                                    (util::div_ceil(scaled_threshold_for_maximizer, scaled_lhs) -
                                     1u) *
                                        increment;
                                if (new_estimate <= nrange.upper_bound()) {
                                    result.largest_maximizer = std::move(new_estimate);
                                    scaled_threshold_for_maximizer =
                                        compute_scaled_threshold_for_maximizer();
                                    increment = itr->current_convergent_denominator();
                                    continue;
                                }
                            }

                            found_maximizer = true;
                            break;
                        } // while (scaled_lhs_for_convergent < scaled_threshold_for_maximizer)
                    }     // if (itr != rg.end())
                    else {
                        found_maximizer = true;
                    }

                    if (found_maximizer) {
                        result.largest_maximizer +=
                            util::div_floor(nrange.upper_bound() - result.largest_maximizer,
                                            increment) *
                            increment;
                    }
                } // if (!found_maximizer)
                previous_scaled_lhs_for_convergent = scaled_lhs_for_convergent;

                // Minimizer.
                // Find a new odd convergent.
                ++itr;
                scaled_lhs_for_convergent =
                    util::abs((itr->current_convergent_numerator() << approx_info.shift_amount) -
                              approx_info.multiplier * itr->current_convergent_denominator());

                if (!found_minimizer) {
                    if (itr != xi_cf.end()) {
                        while (scaled_lhs_for_convergent <= scaled_threshold_for_minimizer) {
                            auto scaled_lhs = scaled_lhs_for_convergent;
                            auto increment = itr->current_convergent_denominator();
                            {
                                auto s = util::div_floor(scaled_threshold_for_minimizer -
                                                             scaled_lhs_for_convergent,
                                                         previous_scaled_lhs_for_convergent);

                                scaled_lhs += s * previous_scaled_lhs_for_convergent;
                                increment -= std::move(s) * itr->previous_convergent_denominator();
                            }

                            if (!util::is_zero(scaled_lhs)) {
                                auto new_estimate =
                                    result.smallest_minimizer +
                                    util::div_floor(scaled_threshold_for_minimizer, scaled_lhs) *
                                        increment;
                                if (new_estimate <= nrange.upper_bound()) {
                                    result.smallest_minimizer = std::move(new_estimate);
                                    scaled_threshold_for_minimizer =
                                        compute_scaled_threshold_for_minimizer();
                                    continue;
                                }

                                result.smallest_minimizer +=
                                    util::div_floor(nrange.upper_bound() -
                                                        result.smallest_minimizer,
                                                    increment) *
                                    increment;
                            }

                            found_minimizer = true;
                            break;
                        } // while (scaled_lhs_for_convergent <= scaled_threshold_for_minimizer)
                    }     // if (itr != rg.end())
                    else {
                        found_minimizer = true;
                    }
                } // if (!found_minimizer)
                previous_scaled_lhs_for_convergent = scaled_lhs_for_convergent;

                if (found_minimizer && found_maximizer) {
                    break;
                }

                ++itr;
            } // while (true)

            return result;
        }
    }
}

#endif
