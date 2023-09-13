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
#include "gosper_continued_fraction.h"
#include "rational_continued_fraction.h"
#include "bigint.h"

namespace jkj {
    namespace idiv {
        struct multiply_shift_info {
            bigint::int_var multiplier;
            std::size_t shift_amount;
        };

        // Given an interval of rational numbers, find the smallest nonnegative integer k such that
        // at least one number of the form m/2^k for an integer m belongs to the interval, find
        // such m with the smallest absolute value, and then return (m,k).
        template <class RationalInterval>
        constexpr multiply_shift_info find_optimal_multiply_shift(RationalInterval const& itv) {
            return itv.visit([](auto&& itv) -> multiply_shift_info {
                using enum interval_type_t;
                constexpr auto itv_type = itv.interval_type();
                static_assert(itv_type != empty);

                if constexpr (itv_type == entire) {
                    return {0u, 0u};
                }
                else if constexpr (itv_type == bounded_below_open ||
                                   itv_type == bounded_below_closed) {
                    if (util::is_strictly_negative(itv.lower_bound().numerator)) {
                        return {0u, 0u};
                    }

                    auto multiplier = itv.left_endpoint_type() == endpoint_type_t::open
                                          ? util::div_floor(itv.lower_bound().numerator,
                                                            itv.lower_bound().denominator) +
                                                1u
                                          : util::div_ceil(itv.lower_bound().numerator,
                                                           itv.lower_bound().denominator);
                    return {std::move(multiplier), 0u};
                }
                else if constexpr (itv_type == bounded_above_open ||
                                   itv_type == bounded_above_closed) {
                    if (util::is_strictly_positive(itv.upper_bound().numerator)) {
                        return {0u, 0u};
                    }

                    auto multiplier = itv.right_endpoint_type() == endpoint_type_t::open
                                          ? util::div_ceil(itv.upper_bound().numerator,
                                                           itv.upper_bound().denominator) -
                                                1u
                                          : util::div_floor(itv.upper_bound().numerator,
                                                            itv.upper_bound().denominator);
                    return {std::move(multiplier), 0u};
                }
                else {
                    bigint::sign_t interval_sign = bigint::sign_t::positive;
                    if (util::is_zero(itv.lower_bound().numerator)) {
                        if constexpr (itv.left_endpoint_type() == endpoint_type_t::closed) {
                            return {0u, 0u};
                        }
                    }
                    else if (util::is_zero(itv.upper_bound().numerator)) {
                        if constexpr (itv.right_endpoint_type() == endpoint_type_t::closed) {
                            return {0u, 0u};
                        }
                        interval_sign = bigint::sign_t::negative;
                    }
                    else if (util::is_strictly_negative(itv.lower_bound().numerator) &&
                             util::is_strictly_positive(itv.upper_bound().numerator)) {
                        return {0u, 0u};
                    }

                    // k = ceil(log2(1/Delta))       if itv is not open,
                    // k = floor(log2(1/Delta)) + 1  if itv is open.
                    auto k = [&] {
                        auto const delta = itv.upper_bound() - itv.lower_bound();
                        util::constexpr_assert(util::is_strictly_positive(delta.numerator));

                        return itv_type == bounded_open
                                   ? trunc_floor_log2_div(delta.denominator,
                                                          util::abs(delta.numerator)) +
                                         1u
                                   : trunc_ceil_log2_div(delta.denominator,
                                                         util::abs(delta.numerator));
                    }();

                    auto multiplier = [&] {
                        if (interval_sign == bigint::sign_t::positive) {
                            // Take the left-most lattice point.
                            if constexpr (itv.left_endpoint_type() == endpoint_type_t::open) {
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
                            if constexpr (itv.right_endpoint_type() == endpoint_type_t::open) {
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

                    if (util::is_even(multiplier)) {
                        k -= factor_out_power_of_2(multiplier);
                    }
                    else {
                        if (interval_sign == bigint::sign_t::positive) {
                            auto next_lattice_point = multiplier + 1u;
                            if ((itv.right_endpoint_type() == endpoint_type_t::open &&
                                 next_lattice_point * itv.upper_bound().denominator <
                                     (itv.upper_bound().numerator << k)) ||
                                (itv.right_endpoint_type() == endpoint_type_t::closed &&
                                 next_lattice_point * itv.upper_bound().denominator <=
                                     (itv.upper_bound().numerator << k))) {
                                multiplier = std::move(next_lattice_point);
                                k -= factor_out_power_of_2(multiplier);
                            }
                        }
                        else {
                            auto next_lattice_point = multiplier - 1u;
                            if ((itv.left_endpoint_type() == endpoint_type_t::open &&
                                 next_lattice_point * itv.lower_bound().denominator >
                                     (itv.lower_bound().numerator << k)) ||
                                (itv.left_endpoint_type() == endpoint_type_t::closed &&
                                 next_lattice_point * itv.lower_bound().denominator <=
                                     (itv.lower_bound().numerator << k))) {
                                multiplier = std::move(next_lattice_point);
                                k -= factor_out_power_of_2(multiplier);
                            }
                        }
                    }

                    return {std::move(multiplier), k};
                }
            });
        }

        // For a given real number x and a positive integer nmax, find the smallest nonnegative
        // integer k such that there exists an integer m satisfying
        // floor(nx) = floor(nm/2^k) for all n = 1, ... , nmax.
        // The number x is specified in terms of a continued fraction generator giving its continued
        // fraction expansion. The generator needs to have index_tracker and
        // previous_previous_convergent_tracker within it, and it also needs to be at its initial
        // stage, i.e., the call to current_index() without calling update() should return -1.
        // After the function returns, the generator is terminated if x is rational and its
        // denominator is at most nmax.
        template <class ContinuedFractionGenerator>
        constexpr multiply_shift_info find_optimal_multiply_shift(ContinuedFractionGenerator& cf,
                                                                  bigint::uint_var const& nmax) {
            return find_optimal_multiply_shift(find_floor_quotient_range(cf, nmax));
        }

        struct multiply_add_shift_info {
            bigint::int_var multiplier;
            bigint::int_var adder;
            std::size_t shift_amount;
        };

        // Given real numbers x, y and a range [nmin:nmax] of integers, find a triple (k,m,s) of
        // integers such that
        // (1) k >= 0,
        // (2) floor(nx + y) = floor((nm + s)/2^k) holds for all n in [nmin:nmax], and
        // (3) floor(nx) = floor(nm/2^k) holds for all n in [0:nmax-nmin].
        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY>
        constexpr multiply_add_shift_info find_suboptimal_multiply_add_shift(
            ContinuedFractionGeneratorX& xcf, ContinuedFractionGeneratorY& ycf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            // TODO: deal with possible rational dependence between x and y.
            using impl_type_x = typename ContinuedFractionGeneratorX::impl_type;
            using impl_type_y = typename ContinuedFractionGeneratorY::impl_type;
            auto xcf_copy = xcf;

            util::constexpr_assert(nrange.upper_bound() > nrange.lower_bound());
            auto const& nmin = nrange.lower_bound();
            auto const nlength = util::abs(nrange.upper_bound() - nrange.lower_bound());

            // Step 1. Find the range of xi satisfying floor(nx) = floor(nxi) for all n.
            auto xi_range = find_floor_quotient_range(xcf_copy, nlength);
            auto xi_info = find_optimal_multiply_shift(xi_range);

            // Step 2. Subtract out the integer part of y.
            auto floor_y = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper<impl_type_x, impl_type_y>{
                        xcf.copy_internal_implementation(),
                        ycf.copy_internal_implementation(),
                        {// numerator
                         0, nmin, 1, 0,
                         // denominator
                         0, 0, 0, 1}});
                cf.update();
                return cf.current_partial_fraction().denominator;
            }();

            // Step 3. Determine if any of L, R is empty.
            bool is_L_empty = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper<impl_type_x, impl_type_y>{
                        xcf.copy_internal_implementation(),
                        ycf.copy_internal_implementation(),
                        {// numerator
                         0, util::to_signed(xi_range.lower_bound().denominator) + nmin, 1, 0,
                         // denominator
                         0, 0, 0, 1}});
                cf.update();

                auto floor = cf.current_partial_fraction().denominator;
                return floor > xi_range.lower_bound().numerator + floor_y;
            }();
            bool is_R_empty = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper<impl_type_x, impl_type_y>{
                        xcf.copy_internal_implementation(),
                        ycf.copy_internal_implementation(),
                        {// numerator
                         0, util::to_signed(xi_range.upper_bound().denominator) + nmin, 1, 0,
                         // denominator
                         0, 0, 0, 1}});
                cf.update();

                auto floor = cf.current_partial_fraction().denominator;
                return floor < xi_range.upper_bound().numerator + floor_y;
            }();

            bigint::int_var adder = 0;
            if (is_L_empty || is_R_empty) {
                if (is_L_empty) {
                    util::constexpr_assert(!is_R_empty);
                    adder = ((((xi_info.multiplier * xi_range.lower_bound().denominator) >>
                               xi_info.shift_amount) +
                              1)
                             << xi_info.shift_amount) -
                            xi_info.multiplier * xi_range.lower_bound().denominator;
                }
                else {
                    util::constexpr_assert(is_R_empty);
                    adder = 0;
                }
            }
            else {
                // Step 4. Find mu and nu.
                // Find floor(q_* y).
                auto floor_qstar_y = [&] {
                    auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                        cntfrc::impl::binary_gosper<impl_type_x, impl_type_y>{
                            xcf.copy_internal_implementation(),
                            ycf.copy_internal_implementation(),
                            {// numerator
                             0, nmin * xi_range.lower_bound().denominator,
                             util::to_signed(xi_range.lower_bound().denominator), 0,
                             // denominator
                             0, 0, 0, 1}});
                    cf.update();
                    return cf.current_partial_fraction().denominator;
                }();
                bigint::uint_var mu = 0u;
                bigint::uint_var nu = 0u;
                if (xcf_copy.terminated()) {
                    // When x is effectively rational.
                    // Use the relation vp == -1 (mod q).
                    mu = util::abs(((floor_qstar_y + 1u) * xi_range.upper_bound().denominator) %
                                   xi_range.lower_bound().denominator);
                    mu += ((nlength - mu) / xi_range.lower_bound().denominator) *
                          xi_range.lower_bound().denominator;

                    nu = util::abs((floor_qstar_y * xi_range.upper_bound().denominator) %
                                   xi_range.lower_bound().denominator);
                }
                else {
                    // When x is effectively irrational.
                    // Use the relation p_* q^* == -1 (mod q_*).
                    bool computed_mu = false;
                    nu = util::abs((floor_qstar_y * xi_range.upper_bound().denominator) %
                                   xi_range.lower_bound().denominator);
                    unsigned int l = 1u;
                    while (true) {
                        auto ceiling = [&] {
                            auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                                cntfrc::impl::binary_gosper<impl_type_x, impl_type_y>{
                                    xcf.copy_internal_implementation(),
                                    ycf.copy_internal_implementation(),
                                    {// numerator
                                     0, nmin * xi_range.lower_bound().denominator,
                                     util::to_signed(xi_range.lower_bound().denominator),
                                     -floor_qstar_y - l,
                                     // denominator
                                     0, util::to_signed(xi_range.lower_bound().denominator), 0,
                                     -xi_range.lower_bound().numerator}});
                            cf.update();
                            return -cf.current_partial_fraction().denominator;
                        }();

                        if (ceiling > nlength) {
                            break;
                        }

                        auto b =
                            util::abs(((floor_qstar_y + l) * xi_range.upper_bound().denominator) %
                                      xi_range.lower_bound().denominator);

                        if (!computed_mu) {
                            if (b < ceiling) {
                                mu = b;
                                mu +=
                                    util::abs(util::div_floor(ceiling - b - 1,
                                                              xi_range.lower_bound().denominator) *
                                              xi_range.lower_bound().denominator);
                                computed_mu = true;
                            }
                        }

                        if (b >= ceiling) {
                            nu = b;
                        }
                        else {
                            b += util::div_ceil(util::abs(ceiling - b),
                                                xi_range.lower_bound().denominator) *
                                 xi_range.lower_bound().denominator;
                            if (b <= nlength) {
                                nu = b;
                            }
                        }

                        ++l;
                    }

                    if (!computed_mu) {
                        mu = util::abs(((floor_qstar_y + l) * xi_range.upper_bound().denominator) %
                                       xi_range.lower_bound().denominator);
                        mu += ((nlength - mu) / xi_range.lower_bound().denominator) *
                              xi_range.lower_bound().denominator;
                    }
                }

                // If xi is precisely p_*/q_* and mu < nu, we may need to be careful.
                if (xi_info.multiplier == xi_range.lower_bound().numerator &&
                    bigint::uint_var::power_of_2(xi_info.shift_amount) ==
                        xi_range.lower_bound().denominator &&
                    mu < nu && (nu - mu).factor_out_power_of_2() >= xi_info.shift_amount) {
                    // Find xi with the next smallest k.
                    while (true) {
                        ++xi_info.shift_amount;
                        xi_info.multiplier = util::div_ceil(
                            (xi_range.lower_bound().numerator << xi_info.shift_amount),
                            xi_range.lower_bound().denominator);
                        if (util::is_even(xi_info.multiplier)) {
                            ++xi_info.multiplier;

                            if (xi_info.multiplier * xi_range.upper_bound().denominator <
                                (xi_range.upper_bound().numerator << xi_info.shift_amount)) {
                                break;
                            }
                        }
                        else {
                            break;
                        }
                    }
                }
                adder = ((((xi_info.multiplier * nu) >> xi_info.shift_amount) + 1)
                         << xi_info.shift_amount) -
                        xi_info.multiplier * nu;
            }

            adder += (floor_y <<= xi_info.shift_amount);
            adder -= nmin * xi_info.multiplier;
            return {std::move(xi_info.multiplier), std::move(adder), xi_info.shift_amount};
        }

        // Given real numbers x, y and a range [nmin:nmax] of integers, find the smallest minimizer
        // and the largest maximizer of (nx+y) - floor(nx+y).
        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY>
        constexpr extrema_of_fractional_part_output<bigint::int_var>
        find_extrema_of_fractional_part(
            ContinuedFractionGeneratorX& xcf, ContinuedFractionGeneratorY& ycf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            extrema_of_fractional_part_output<bigint::int_var> result{nrange.lower_bound(),
                                                                      nrange.lower_bound()};

            // First, find fine enough approximations of x and y.
            auto approx_info = find_suboptimal_multiply_add_shift(xcf, ycf, nrange);
            auto xi_cf =
                cntfrc::make_generator<cntfrc::index_tracker, cntfrc::partial_fraction_tracker,
                                       cntfrc::previous_previous_convergent_tracker>(
                    cntfrc::impl::rational{cntfrc::projective_rational{
                        approx_info.multiplier,
                        bigint::uint_var::power_of_2(approx_info.shift_amount)}});

            // RHS times 2^k.
            auto compute_scaled_threshold_for_maximizer = [&](auto const& n) {
                auto temp = n * approx_info.multiplier + approx_info.adder;
                temp =
                    (((temp >> approx_info.shift_amount) + 1) << approx_info.shift_amount) - temp;
                return util::abs(std::move(temp));
            };
            auto scaled_threshold_for_maximizer =
                compute_scaled_threshold_for_maximizer(nrange.lower_bound());
            auto compute_scaled_threshold_for_minimizer = [&](auto const& n) {
                auto temp = n * approx_info.multiplier + approx_info.adder;
                temp = temp - ((temp >> approx_info.shift_amount) << approx_info.shift_amount);
                return util::abs(std::move(temp));
            };
            auto scaled_threshold_for_minimizer =
                compute_scaled_threshold_for_minimizer(nrange.lower_bound());

            bool found_minimizer = false;
            bool found_maximizer = false;
            while (!found_minimizer || !found_maximizer) {
                // Maximizer.
                // Find a new even convergent.
                xi_cf.update();
                if (!found_maximizer) {
                    if (xi_cf.terminated()) {
                        // If we already have reached to the exact value, just add the multiple of
                        // its denominator as many times as allowed.
                        result.largest_maximizer +=
                            (((nrange.upper_bound() - result.largest_maximizer) >>
                              approx_info.shift_amount)
                             << approx_info.shift_amount);
                        found_maximizer = true;
                    }
                    else {
                        // Otherwise, see if the current convergent satisfies the condition.
                        auto scaled_lefthand_side_for_convergent =
                            approx_info.multiplier * xi_cf.current_convergent_denominator() -
                            (xi_cf.current_convergent_numerator() << approx_info.shift_amount);

                        auto semiconvergent = xi_cf.current_index() >= 2
                                                  ? xi_cf.previous_previous_convergent()
                                                  : xi_cf.current_convergent();
                        auto scaled_lefthand_side =
                            xi_cf.current_index() >= 2
                                ? approx_info.multiplier * semiconvergent.denominator -
                                      (semiconvergent.numerator << approx_info.shift_amount)
                                : scaled_lefthand_side_for_convergent;

                        std::size_t semiconvergent_coeff = 0;
                        while (!found_maximizer) {
                            // If the current convergent does not satisfy the condition, then move
                            // on to the next convergent.
                            if (scaled_lefthand_side_for_convergent >=
                                scaled_threshold_for_maximizer) {
                                break;
                            }

                            // Otherwise, find the first even semiconvergent still satisfying the
                            // condition.
                            if (xi_cf.current_index() >= 2) {
                                do {
                                    ++semiconvergent_coeff;
                                    semiconvergent.numerator +=
                                        xi_cf.previous_convergent_numerator();
                                    semiconvergent.denominator +=
                                        xi_cf.previous_convergent_denominator();

                                    scaled_lefthand_side -= (xi_cf.previous_convergent_numerator()
                                                             << approx_info.shift_amount);
                                    scaled_lefthand_side += approx_info.multiplier *
                                                            xi_cf.previous_convergent_denominator();
                                } while (semiconvergent_coeff <
                                             xi_cf.current_partial_fraction().denominator &&
                                         scaled_lefthand_side >= scaled_threshold_for_maximizer);
                            }

                            // Update the current estimate of the maximizer with the denominator of
                            // the found semiconvergent.
                            auto new_estimate = result.largest_maximizer +
                                                (util::div_ceil(scaled_threshold_for_maximizer,
                                                                scaled_lefthand_side) -
                                                 1) *
                                                    semiconvergent.denominator;

                            if (new_estimate >= nrange.upper_bound()) {
                                new_estimate =
                                    result.largest_maximizer +
                                    util::div_floor(nrange.upper_bound() - result.largest_maximizer,
                                                    semiconvergent.denominator) *
                                        semiconvergent.denominator;
                                found_maximizer = true;
                            }
                            result.largest_maximizer = std::move(new_estimate);
                            scaled_threshold_for_maximizer =
                                compute_scaled_threshold_for_maximizer(result.largest_maximizer);
                        }
                    }
                }

                // Minimizer.
                // Find a new odd convergent.
                xi_cf.update();
                if (!found_minimizer) {
                    if (xi_cf.terminated()) {
                        // If we already have reached to the exact value, there is nothing else to
                        // do.
                        found_minimizer = true;
                    }
                    else {
                        // Otherwise, see if the current convergent satisfies the condition.
                        // If quantity below is zero, then the current convergent is the exact
                        // value.
                        auto scaled_lefthand_side_for_convergent =
                            (xi_cf.current_convergent_numerator() << approx_info.shift_amount) -
                            approx_info.multiplier * xi_cf.current_convergent_denominator();

                        auto semiconvergent = xi_cf.previous_previous_convergent();
                        auto scaled_lefthand_side =
                            (semiconvergent.numerator << approx_info.shift_amount) -
                            approx_info.multiplier * semiconvergent.denominator;

                        std::size_t semiconvergent_coeff = 0;
                        while (!found_minimizer) {
                            // If the current convergent does not satisfy the condition, then move
                            // on to the next convergent.
                            if (scaled_lefthand_side_for_convergent >
                                scaled_threshold_for_minimizer) {
                                break;
                            }

                            // Otherwise, find the first even semiconvergent still satisfying the
                            // condition.
                            do {
                                ++semiconvergent_coeff;
                                semiconvergent.numerator += xi_cf.previous_convergent_numerator();
                                semiconvergent.denominator +=
                                    xi_cf.previous_convergent_denominator();

                                scaled_lefthand_side += (xi_cf.previous_convergent_numerator()
                                                         << approx_info.shift_amount);
                                scaled_lefthand_side -= approx_info.multiplier *
                                                        xi_cf.previous_convergent_denominator();
                            } while (semiconvergent_coeff <
                                         xi_cf.current_partial_fraction().denominator &&
                                     scaled_lefthand_side > scaled_threshold_for_minimizer);

                            // Update the current estimate of the maximizer with the denominator of
                            // the found semiconvergent.
                            if (util::is_zero(scaled_lefthand_side_for_convergent) &&
                                semiconvergent_coeff ==
                                    xi_cf.current_partial_fraction().denominator) {
                                // If the current convergent is the exact value and there is no
                                // semiconvergent satisfying the condition, then the set is empty.
                                // There is nothing else to do in that case.
                                found_minimizer = true;
                            }
                            else {
                                auto new_estimate = result.smallest_minimizer +
                                                    util::div_floor(scaled_threshold_for_minimizer,
                                                                    scaled_lefthand_side) *
                                                        semiconvergent.denominator;

                                if (new_estimate >= nrange.upper_bound()) {
                                    new_estimate = result.smallest_minimizer +
                                                   util::div_floor(nrange.upper_bound() -
                                                                       result.smallest_minimizer,
                                                                   semiconvergent.denominator) *
                                                       semiconvergent.denominator;
                                    found_minimizer = true;
                                }
                                result.smallest_minimizer = std::move(new_estimate);
                                scaled_threshold_for_minimizer =
                                    compute_scaled_threshold_for_minimizer(
                                        result.smallest_minimizer);
                            }
                        }
                    }
                }
            }

            return result;
        }
#if 0
				
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
            util::constexpr_assert(x.denominator <= nmax);

            using continued_fractions_calculator_type =
                caching_continued_fraction<rational<bigint::uint_var, bigint::uint_var>>;

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
                    util::constexpr_assert(delta.denominator >= delta.numerator);
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
#endif // 0

    }
}

#endif
