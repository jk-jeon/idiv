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

#ifndef JKJ_HEADER_IDIV_OPTIMAL_MULTIPLY_SHIFT
#define JKJ_HEADER_IDIV_OPTIMAL_MULTIPLY_SHIFT

#include "best_rational_approx.h"

namespace jkj {
    namespace idiv {
        // Given an interval of rational numbers, find the smallest nonnegative integer k such that
        // at least one number of the form m/2^k for an integer m belongs to the interval, find
        // such m with the smallest absolute value, and then return (m,k).

        template <class Int>
        struct binary_lattice_point {
            Int numerator;
            std::size_t shift_amount;
        };

        template <
            class RationalInterval,
            class Int = std::remove_cvref_t<
                decltype(std::declval<typename std::remove_cvref_t<RationalInterval>::value_type>()
                             .numerator)>>
        constexpr binary_lattice_point<Int>
        find_optimal_binary_lattice_point(RationalInterval const& vitv) {
            return vitv.visit([](auto&& itv) -> binary_lattice_point<Int> {
                using enum interval_type_t;
                using itv_type = std::remove_cvref_t<decltype(itv)>;
                static_assert(itv_type::interval_type() != empty);

                if constexpr (itv_type::interval_type() == entire) {
                    return {0, 0u};
                }
                else if constexpr (itv_type::interval_type() == bounded_below_open ||
                                   itv_type::interval_type() == bounded_below_closed) {
                    if (util::is_strictly_negative(itv.lower_bound().numerator)) {
                        return {0, 0u};
                    }

                    auto numerator = itv_type::left_boundary_type() == boundary_type_t::open
                                         ? util::div_floor(itv.lower_bound().numerator,
                                                           itv.lower_bound().denominator) +
                                               1u
                                         : util::div_ceil(itv.lower_bound().numerator,
                                                          itv.lower_bound().denominator);
                    return {std::move(numerator), 0u};
                }
                else if constexpr (itv_type::interval_type() == bounded_above_open ||
                                   itv_type::interval_type() == bounded_above_closed) {
                    if (util::is_strictly_positive(itv.upper_bound().numerator)) {
                        return {0, 0u};
                    }

                    auto numerator = itv_type::right_boundary_type() == boundary_type_t::open
                                         ? util::div_ceil(itv.upper_bound().numerator,
                                                          itv.upper_bound().denominator) -
                                               1u
                                         : util::div_floor(itv.upper_bound().numerator,
                                                           itv.upper_bound().denominator);
                    return {std::move(numerator), 0u};
                }
                else {
                    auto interval_sign = util::sign_t::positive;
                    if (util::is_zero(itv.lower_bound().numerator)) {
                        if constexpr (itv_type::left_boundary_type() == boundary_type_t::closed) {
                            return {0, 0u};
                        }
                    }
                    else if (util::is_zero(itv.upper_bound().numerator)) {
                        if constexpr (itv_type::right_boundary_type() == boundary_type_t::closed) {
                            return {0, 0u};
                        }
                        interval_sign = util::sign_t::negative;
                    }
                    else if (util::is_strictly_negative(itv.lower_bound().numerator) &&
                             util::is_strictly_positive(itv.upper_bound().numerator)) {
                        return {0, 0u};
                    }

                    // k = ceil(log2(1/Delta))       if itv is not open,
                    // k = floor(log2(1/Delta)) + 1  if itv is open.
                    auto k = [&] {
                        auto const delta = itv.upper_bound() - itv.lower_bound();
                        util::constexpr_assert(util::is_strictly_positive(delta.numerator));

                        return itv_type::interval_type() == bounded_open
                                   ? util::trunc_floor_log2_div(delta.denominator,
                                                                util::abs(delta.numerator)) +
                                         1u
                                   : util::trunc_ceil_log2_div(delta.denominator,
                                                               util::abs(delta.numerator));
                    }();

                    auto numerator = [&] {
                        if (interval_sign == util::sign_t::positive) {
                            // Take the left-most lattice point.
                            if constexpr (itv_type::left_boundary_type() == boundary_type_t::open) {
                                return util::div_floor((itv.lower_bound().numerator << k),
                                                       itv.lower_bound().denominator) +
                                       1u;
                            }
                            else {
                                return util::div_ceil((itv.lower_bound().numerator << k),
                                                      itv.lower_bound().denominator);
                            }
                        }
                        else {
                            // Take the right-most lattice point.
                            if constexpr (itv_type::right_boundary_type() ==
                                          boundary_type_t::open) {
                                return util::div_ceil((itv.upper_bound().numerator << k),
                                                      itv.upper_bound().denominator) -
                                       1u;
                            }
                            else {
                                return util::div_floor((itv.upper_bound().numerator << k),
                                                       itv.upper_bound().denominator);
                            }
                        }
                    }();

                    if (k > 0) {
                        auto factor_out_power_of_2_limited = [&] {
                            auto pow2_factors = util::factor_out_power_of_2(numerator);
                            if (pow2_factors > k) {
                                numerator <<= (pow2_factors - k);
                                k = 0;
                            }
                            else {
                                k -= pow2_factors;
                            }
                        };

                        if (util::is_even(numerator)) {
                            factor_out_power_of_2_limited();
                        }
                        else {
                            if (interval_sign == util::sign_t::positive) {
                                auto next_lattice_point = numerator + 1u;
                                if ((itv_type::right_boundary_type() == boundary_type_t::open &&
                                     next_lattice_point * itv.upper_bound().denominator <
                                         (itv.upper_bound().numerator << k)) ||
                                    (itv_type::right_boundary_type() == boundary_type_t::closed &&
                                     next_lattice_point * itv.upper_bound().denominator <=
                                         (itv.upper_bound().numerator << k))) {
                                    numerator = std::move(next_lattice_point);
                                    factor_out_power_of_2_limited();
                                }
                            }
                            else {
                                auto next_lattice_point = numerator - 1u;
                                if ((itv_type::left_boundary_type() == boundary_type_t::open &&
                                     next_lattice_point * itv.lower_bound().denominator >
                                         (itv.lower_bound().numerator << k)) ||
                                    (itv_type::left_boundary_type() == boundary_type_t::closed &&
                                     next_lattice_point * itv.lower_bound().denominator <=
                                         (itv.lower_bound().numerator << k))) {
                                    numerator = std::move(next_lattice_point);
                                    factor_out_power_of_2_limited();
                                }
                            }
                        }
                    }

                    return {std::move(numerator), k};
                }
            });
        }

        // For a given real number x and a positive integer nmax, find the smallest nonnegative
        // integer k such that there exists an integer m satisfying
        // floor(nx) = floor(nm/2^k) for all n = 1, ... , nmax.
        // The number x is specified in terms of a continued fraction generator giving its continued
        // fraction expansion. The generator needs to have index_tracker and
        // previous_previous_convergent_tracker within it, and it also needs to be at its initial
        // stage, i.e., the call to current_index() without calling
        // proceed_to_next_partial_fraction() should return -1. After the function returns, the
        // generator is terminated if x is rational and its denominator is at most nmax.

        template <class Int>
        struct multiply_shift_info {
            Int multiplier;
            std::size_t shift_amount;
        };

        template <class ContinuedFractionGenerator, class UInt,
                  class Int = std::remove_cvref_t<
                      decltype(std::declval<typename std::remove_cvref_t<
                                   ContinuedFractionGenerator>::partial_fraction_type>()
                                   .denominator)>>
        constexpr multiply_shift_info<Int>
        find_optimal_multiply_shift(ContinuedFractionGenerator&& cf, UInt const& nmax) {
            static_assert(cntfrc::has_mixins<ContinuedFractionGenerator, cntfrc::index_tracker,
                                             cntfrc::convergent_tracker>(),
                          "the passed continued fraction generator must include index_tracker and "
                          "convergent_tracker");

            auto result = find_optimal_binary_lattice_point(
                find_floor_quotient_range(std::forward<ContinuedFractionGenerator>(cf), nmax));
            return {std::move(result.numerator), result.shift_amount};
        }
    }
}

#endif
