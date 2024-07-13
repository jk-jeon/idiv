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

#ifndef JKJ_HEADER_IDIV
#define JKJ_HEADER_IDIV

#include "optimal_multiply_shift.h"
#include "xi_zeta_region.h"

namespace jkj {
    namespace idiv {
        // Given real numbers x, y, zeta and a range [nmin:nmax] of integers, find a maximizer of
        // (floor(nx+y) - zeta) / n.
        // Precondition: 0 < nmin <= nmax.
        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY,
                  class ContinuedFractionGeneratorZeta>
        constexpr bigint::int_var find_maximizer_of_floor_subtract_quotient_positive_range(
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            ContinuedFractionGeneratorZeta&& zetacf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorX>::template is_implementing_mixins<
                    cntfrc::previous_previous_convergent_tracker, cntfrc::interval_tracker>(),
                "the first continued fraction generator must implement "
                "previous_previous_convergent_tracker and "
                "interval_tracker");
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorY>::template is_implementing_mixins<
                    cntfrc::interval_tracker>(),
                "the second continued fraction generator must implement interval_tracker");
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorZeta>::
                    template is_implementing_mixins<cntfrc::index_tracker,
                                                    cntfrc::previous_previous_convergent_tracker>(),
                "the third continued fraction generator must implement index_tracker and "
                "previous_previous_convergent_tracker");

            util::constexpr_assert(util::is_strictly_positive(nrange.lower_bound()));

            // Find good enough approximations of x and y.
            auto approx_x_y_info =
                find_simultaneous_multiply_add_shift(xcf.copy(), ycf.copy(), nrange);
            auto xcf_copy = cntfrc::make_caching_generator(
                cntfrc::make_generator<cntfrc::partial_fraction_tracker,
                                       cntfrc::convergent_tracker>(cntfrc::impl::rational{
                    approx_x_y_info.multiplier,
                    bigint::uint_var::power_of_2(approx_x_y_info.shift_amount)}));

            // Find a good enough approximation of zeta.
            // Because of potential aliasing of zetacf with xcf/ycf, we make a copy here.
            // "+ 1u" is just to make the case nmin == nmax work correctly.
            auto approx_zeta_info = find_best_rational_approx(
                zetacf, util::abs(nrange.upper_bound() - nrange.lower_bound()) + 1u);

            // Find the smallest minimizer of the fractional part.
            auto const n00 =
                find_extremizers_of_fractional_part(xcf, ycf, nrange).smallest_minimizer;

            // Solve the maximization problem on the left.
            auto left_maximizer = n00;
            while (left_maximizer > nrange.lower_bound()) {
                // n1 is the largest minimizer of (floor(nx) + 1) / n, which is the largest multiple
                // of the largest maximizer of nx - floor(nx), where n is in [1:n0 - nmin].
                auto const n1 = [&] {
                    xcf_copy.rewind();
                    auto const nmax = util::abs(left_maximizer - nrange.lower_bound());
                    auto smallest_minimizer =
                        find_extremizers_of_fractional_part(xcf_copy, nmax).largest_maximizer;
                    return util::div_floor(nmax, smallest_minimizer) * smallest_minimizer;
                }();

                // Check if (floor(n1 x) + 1)/n1 > (floor(n0 x + y) - zeta)/n0, or equivalently,
                // ceil(n1 zeta) > n1 floor(n0 x + y) - n0 (floor(n1 x) + 1).
                if (util::div_ceil(n1 * approx_zeta_info.above.numerator,
                                   approx_zeta_info.above.denominator) >
                    n1 * ((left_maximizer * approx_x_y_info.multiplier + approx_x_y_info.adder) >>
                          approx_x_y_info.shift_amount) -
                        left_maximizer *
                            (((n1 * approx_x_y_info.multiplier) >> approx_x_y_info.shift_amount) +
                             1)) {
                    break;
                }
                left_maximizer -= n1;
            }

            // Solve the maximization problem on the right.
            auto right_maximizer = n00;
            while (right_maximizer < nrange.upper_bound()) {
                // n1 is the largest maximizer of floor(nx)/ n, which is the largest multiple of the
                // smallest minimizer of nx - floor(nx), where n is in [1:nmax - n0].
                auto const n1 = [&] {
                    xcf_copy.rewind();
                    auto const nmax = util::abs(nrange.upper_bound() - right_maximizer);
                    auto smallest_maximizer =
                        find_extremizers_of_fractional_part(xcf_copy, nmax).smallest_minimizer;
                    return util::div_floor(nmax, smallest_maximizer) * smallest_maximizer;
                }();

                // Check if floor(n1 x)/n1 < (floor(n0 x + y) - zeta)/n0, or equivalently,
                // floor(n1 zeta) < n1 floor(n0 x + y) - n0 floor(n1 x).
                if (util::div_floor(n1 * approx_zeta_info.below.numerator,
                                    approx_zeta_info.below.denominator) <
                    n1 * ((right_maximizer * approx_x_y_info.multiplier + approx_x_y_info.adder) >>
                          approx_x_y_info.shift_amount) -
                        right_maximizer *
                            ((n1 * approx_x_y_info.multiplier) >> approx_x_y_info.shift_amount)) {
                    break;
                }
                right_maximizer += n1;
            }

            // Compare the maximizers on the left and the right, and choose the better one.
            // (floor(nl x + y) - zeta)/nl <= (floor(nr x + y) - zeta)/nr holds if and only if
            // nr floor(nl x + y) - nl floor(nr x + y) <= floor((nr - nl) zeta).
            auto const lhs =
                right_maximizer *
                    ((left_maximizer * approx_x_y_info.multiplier + approx_x_y_info.adder) >>
                     approx_x_y_info.shift_amount) -
                left_maximizer *
                    ((right_maximizer * approx_x_y_info.multiplier + approx_x_y_info.adder) >>
                     approx_x_y_info.shift_amount);

            auto const rhs = util::div_floor((right_maximizer - left_maximizer) *
                                                 approx_zeta_info.below.numerator,
                                             approx_zeta_info.below.denominator);

            return lhs <= rhs ? right_maximizer : left_maximizer;
        }

        // Given real numbers x, y, zeta and a range [nmin:nmax] of integers, find a minimizer of
        // (floor(nx+y) - zeta) / n.
        // Precondition: 0 < nmin <= nmax.
        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY,
                  class ContinuedFractionGeneratorZeta>
        constexpr bigint::int_var find_minimizer_of_floor_subtract_quotient_positive_range(
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            ContinuedFractionGeneratorZeta&& zetacf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorX>::template is_implementing_mixins<
                    cntfrc::previous_previous_convergent_tracker, cntfrc::interval_tracker>(),
                "the first continued fraction generator must implement "
                "previous_previous_convergent_tracker and "
                "interval_tracker");
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorY>::template is_implementing_mixins<
                    cntfrc::interval_tracker>(),
                "the second continued fraction generator must implement interval_tracker");
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorZeta>::
                    template is_implementing_mixins<cntfrc::index_tracker,
                                                    cntfrc::previous_previous_convergent_tracker>(),
                "the third continued fraction generator must implement index_tracker and "
                "previous_previous_convergent_tracker");

            util::constexpr_assert(nrange.lower_bound() > 0);

            // Find good enough approximations of x and y.
            auto approx_x_y_info =
                find_simultaneous_multiply_add_shift(xcf.copy(), ycf.copy(), nrange);
            auto xcf_copy = cntfrc::make_caching_generator(
                cntfrc::make_generator<cntfrc::partial_fraction_tracker,
                                       cntfrc::convergent_tracker>(cntfrc::impl::rational{
                    approx_x_y_info.multiplier,
                    bigint::uint_var::power_of_2(approx_x_y_info.shift_amount)}));

            // Find a good enough approximation of zeta.
            // Because of potential aliasing of zetacf with xcf/ycf, we make a copy here.
            // "+ 1u" is just to make the case nmin == nmax work correctly.
            auto approx_zeta_info = find_best_rational_approx(
                zetacf, util::abs(nrange.upper_bound() - nrange.lower_bound()) + 1u);

            // Find the largest maximizer of the fractional part.
            auto const n00 =
                find_extremizers_of_fractional_part(xcf, ycf, nrange).largest_maximizer;

            // Solve the minimization problem on the left.
            auto left_minimizer = n00;
            while (left_minimizer > nrange.lower_bound()) {
                // n1 is the largest maximizer of floor(nx)/ n, which is the largest multiple of the
                // smallest minimizer of nx - floor(nx), where n is in [1:n0 - nmin].
                auto const n1 = [&] {
                    xcf_copy.rewind();
                    auto const nmax = util::abs(left_minimizer - nrange.lower_bound());
                    auto smallest_maximizer =
                        find_extremizers_of_fractional_part(xcf_copy, nmax).smallest_minimizer;
                    return util::div_floor(nmax, smallest_maximizer) * smallest_maximizer;
                }();

                // Check if floor(n1 x)/n1 < (floor(n0 x + y) - zeta)/n0, or equivalently,
                // floor(n1 zeta) < n1 floor(n0 x + y) - n0 floor(n1 x).
                if (util::div_floor(n1 * approx_zeta_info.below.numerator,
                                    approx_zeta_info.below.denominator) <
                    n1 * ((left_minimizer * approx_x_y_info.multiplier + approx_x_y_info.adder) >>
                          approx_x_y_info.shift_amount) -
                        left_minimizer *
                            ((n1 * approx_x_y_info.multiplier) >> approx_x_y_info.shift_amount)) {
                    break;
                }
                left_minimizer -= n1;
            }

            // Solve the minimization problem on the right.
            auto right_minimizer = n00;
            while (right_minimizer < nrange.upper_bound()) {
                // n1 is the largest minimizer of (floor(nx) + 1) / n, which is the largest multiple
                // of the largest maximizer of nx - floor(nx), where n is in [1:nmax - n0].
                auto const n1 = [&] {
                    xcf_copy.rewind();
                    auto const nmax = util::abs(nrange.upper_bound() - right_minimizer);
                    auto smallest_minimizer =
                        find_extremizers_of_fractional_part(xcf_copy, nmax).largest_maximizer;
                    return util::div_floor(nmax, smallest_minimizer) * smallest_minimizer;
                }();

                // Check if (floor(n1 x) + 1)/n1 > (floor(n0 x + y) - zeta)/n0, or equivalently,
                // ceil(n1 zeta) > n1 floor(n0 x + y) - n0 (floor(n1 x) + 1).
                if (util::div_ceil(n1 * approx_zeta_info.above.numerator,
                                   approx_zeta_info.above.denominator) >
                    n1 * ((right_minimizer * approx_x_y_info.multiplier + approx_x_y_info.adder) >>
                          approx_x_y_info.shift_amount) -
                        right_minimizer *
                            (((n1 * approx_x_y_info.multiplier) >> approx_x_y_info.shift_amount) +
                             1)) {
                    break;
                }
                right_minimizer += n1;
            }

            // Compare the minimizers on the left and the right, and choose the better one.
            // (floor(nl x + y) - zeta)/nl <= (floor(nr x + y) - zeta)/nr holds if and only if
            // nr floor(nl x + y) - nl floor(nr x + y) <= floor((nr - nl) zeta).
            auto const lhs =
                right_minimizer *
                    ((left_minimizer * approx_x_y_info.multiplier + approx_x_y_info.adder) >>
                     approx_x_y_info.shift_amount) -
                left_minimizer *
                    ((right_minimizer * approx_x_y_info.multiplier + approx_x_y_info.adder) >>
                     approx_x_y_info.shift_amount);

            auto const rhs = util::div_floor((right_minimizer - left_minimizer) *
                                                 approx_zeta_info.below.numerator,
                                             approx_zeta_info.below.denominator);

            return lhs <= rhs ? left_minimizer : right_minimizer;
        }
    }
}

#endif
