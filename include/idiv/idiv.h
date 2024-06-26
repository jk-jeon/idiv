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

#include "best_rational_approx.h"
#include "caching_generator.h"
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
                using itv_type = std::remove_cvref_t<decltype(itv)>;
                static_assert(itv_type::interval_type() != empty);

                if constexpr (itv_type::interval_type() == entire) {
                    return {0u, 0u};
                }
                else if constexpr (itv_type::interval_type() == bounded_below_open ||
                                   itv_type::interval_type() == bounded_below_closed) {
                    if (util::is_strictly_negative(itv.lower_bound().numerator)) {
                        return {0u, 0u};
                    }

                    auto multiplier = itv_type::left_endpoint_type() == endpoint_type_t::open
                                          ? util::div_floor(itv.lower_bound().numerator,
                                                            itv.lower_bound().denominator) +
                                                1u
                                          : util::div_ceil(itv.lower_bound().numerator,
                                                           itv.lower_bound().denominator);
                    return {std::move(multiplier), 0u};
                }
                else if constexpr (itv_type::interval_type() == bounded_above_open ||
                                   itv_type::interval_type() == bounded_above_closed) {
                    if (util::is_strictly_positive(itv.upper_bound().numerator)) {
                        return {0u, 0u};
                    }

                    auto multiplier = itv_type::right_endpoint_type() == endpoint_type_t::open
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
                        if constexpr (itv_type::left_endpoint_type() == endpoint_type_t::closed) {
                            return {0u, 0u};
                        }
                    }
                    else if (util::is_zero(itv.upper_bound().numerator)) {
                        if constexpr (itv_type::right_endpoint_type() == endpoint_type_t::closed) {
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

                        return itv_type::interval_type() == bounded_open
                                   ? trunc_floor_log2_div(delta.denominator,
                                                          util::abs(delta.numerator)) +
                                         1u
                                   : trunc_ceil_log2_div(delta.denominator,
                                                         util::abs(delta.numerator));
                    }();

                    auto multiplier = [&] {
                        if (interval_sign == bigint::sign_t::positive) {
                            // Take the left-most lattice point.
                            if constexpr (itv_type::left_endpoint_type() == endpoint_type_t::open) {
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
                            if constexpr (itv_type::right_endpoint_type() ==
                                          endpoint_type_t::open) {
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
                            auto pow2_factors = factor_out_power_of_2(multiplier);
                            if (pow2_factors > k) {
                                multiplier <<= (pow2_factors - k);
                                k = 0;
                            }
                            else {
                                k -= pow2_factors;
                            }
                        };

                        if (util::is_even(multiplier)) {
                            factor_out_power_of_2_limited();
                        }
                        else {
                            if (interval_sign == bigint::sign_t::positive) {
                                auto next_lattice_point = multiplier + 1u;
                                if ((itv_type::right_endpoint_type() == endpoint_type_t::open &&
                                     next_lattice_point * itv.upper_bound().denominator <
                                         (itv.upper_bound().numerator << k)) ||
                                    (itv_type::right_endpoint_type() == endpoint_type_t::closed &&
                                     next_lattice_point * itv.upper_bound().denominator <=
                                         (itv.upper_bound().numerator << k))) {
                                    multiplier = std::move(next_lattice_point);
                                    factor_out_power_of_2_limited();
                                }
                            }
                            else {
                                auto next_lattice_point = multiplier - 1u;
                                if ((itv_type::left_endpoint_type() == endpoint_type_t::open &&
                                     next_lattice_point * itv.lower_bound().denominator >
                                         (itv.lower_bound().numerator << k)) ||
                                    (itv_type::left_endpoint_type() == endpoint_type_t::closed &&
                                     next_lattice_point * itv.lower_bound().denominator <=
                                         (itv.lower_bound().numerator << k))) {
                                    multiplier = std::move(next_lattice_point);
                                    factor_out_power_of_2_limited();
                                }
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
        constexpr multiply_shift_info find_optimal_multiply_shift(ContinuedFractionGenerator&& cf,
                                                                  bigint::uint_var const& nmax) {
            static_assert(
                std::remove_cvref_t<ContinuedFractionGenerator>::template is_implementing_mixins<
                    cntfrc::index_tracker, cntfrc::previous_previous_convergent_tracker>(),
                "the passed continued fraction generator must implement index_tracker and "
                "previous_previous_convergent_tracker");

            return find_optimal_multiply_shift(
                find_floor_quotient_range(std::forward<ContinuedFractionGenerator>(cf), nmax));
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
        constexpr multiply_add_shift_info find_simultaneous_multiply_add_shift(
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            // TODO: deal with possible rational dependence between x and y.

            auto xcf_copy = cntfrc::make_generator<cntfrc::index_tracker,
                                                   cntfrc::previous_previous_convergent_tracker>(
                xcf.copy_internal_implementation());

            auto const& nmin = nrange.lower_bound();
            auto const nlength = util::abs(nrange.upper_bound() - nrange.lower_bound());

            // Step 1. Subtract out the integer part of y.
            auto floor_y = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::make_binary_gosper_from_impl(xcf.copy_internal_implementation(),
                                                         ycf.copy_internal_implementation(),
                                                         {// numerator
                                                          0, nmin, 1, 0,
                                                          // denominator
                                                          0, 0, 0, 1}));
                cf.update();
                return cf.current_partial_fraction().denominator;
            }();

            // Step 2. Find the range of xi satisfying floor(nx) = floor(nxi) for all n.
            if (util::is_zero(nlength)) {
                // When nmin = nmax, we can take xi to be a whatever number.
                // For sanity, we choose xi = floor(x).
                xcf_copy.update();
                auto xi = xcf_copy.current_convergent().numerator;
                auto zeta = floor_y - nrange.lower_bound() * xi;
                return {std::move(xi), std::move(zeta), 0};
            }
            auto xi_range = find_floor_quotient_range(xcf_copy, nlength);
            auto xi_info = find_optimal_multiply_shift(xi_range);

            // Step 3. Determine if any of L, R is empty.
            bool is_L_empty = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::make_binary_gosper_from_impl(
                        xcf.copy_internal_implementation(), ycf.copy_internal_implementation(),
                        {// numerator
                         0, util::to_signed(xi_range.lower_bound().denominator) + nmin, 1, 0,
                         // denominator
                         0, 0, 0, 1}));
                cf.update();

                auto floor = cf.current_partial_fraction().denominator;
                return floor > xi_range.lower_bound().numerator + floor_y;
            }();
            bool is_R_empty = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::make_binary_gosper_from_impl(
                        xcf.copy_internal_implementation(), ycf.copy_internal_implementation(),
                        {// numerator
                         0, util::to_signed(xi_range.upper_bound().denominator) + nmin, 1, 0,
                         // denominator
                         0, 0, 0, 1}));
                cf.update();

                auto floor = cf.current_partial_fraction().denominator;
                return floor < xi_range.upper_bound().numerator + floor_y;
            }();

            bigint::int_var adder = 0;
            if (is_L_empty || is_R_empty) {
                if (is_R_empty) {
                    adder = 0;
                }
                else {
                    util::constexpr_assert(is_L_empty);
                    adder = ((((xi_info.multiplier * xi_range.lower_bound().denominator) >>
                               xi_info.shift_amount) +
                              1)
                             << xi_info.shift_amount) -
                            xi_info.multiplier * xi_range.lower_bound().denominator;
                }
            }
            else {
                // Step 4. Find mu and nu.
                // Find floor(q_* y).
                auto floor_qstar_y = [&] {
                    auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                        cntfrc::make_binary_gosper_from_impl(
                            xcf.copy_internal_implementation(), ycf.copy_internal_implementation(),
                            {// numerator
                             0, nmin * xi_range.lower_bound().denominator,
                             util::to_signed(xi_range.lower_bound().denominator), 0,
                             // denominator
                             0, 0, 0, 1}));
                    cf.update();
                    return cf.current_partial_fraction().denominator -
                           xi_range.lower_bound().denominator * floor_y;
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
                    unsigned int l = 0u;
                    while (true) {
                        auto ceiling = [&] {
                            auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                                cntfrc::make_binary_gosper_from_impl(
                                    xcf.copy_internal_implementation(),
                                    ycf.copy_internal_implementation(),
                                    {// numerator
                                     0, nmin * xi_range.lower_bound().denominator,
                                     util::to_signed(xi_range.lower_bound().denominator),
                                     -floor_qstar_y - l -
                                         xi_range.lower_bound().denominator * floor_y,
                                     // denominator
                                     0, util::to_signed(xi_range.lower_bound().denominator), 0,
                                     -xi_range.lower_bound().numerator}));
                            cf.update();
                            return -cf.current_partial_fraction().denominator;
                        }();

                        if (ceiling > nlength) {
                            break;
                        }

                        auto b =
                            util::abs(((floor_qstar_y + l) * xi_range.upper_bound().denominator) %
                                      xi_range.lower_bound().denominator);

                        if (!computed_mu && l != 0u) {
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
                adder = ((((nu * xi_info.multiplier) >> xi_info.shift_amount) + 1)
                         << xi_info.shift_amount) -
                        nu * xi_info.multiplier;
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
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            extrema_of_fractional_part_output<bigint::int_var> result{nrange.lower_bound(),
                                                                      nrange.lower_bound()};

            if (nrange.lower_bound() == nrange.upper_bound()) {
                return result;
            }

            // First, find fine enough approximations of x and y.
            auto approx_info = find_simultaneous_multiply_add_shift(
                std::forward<ContinuedFractionGeneratorX>(xcf),
                std::forward<ContinuedFractionGeneratorY>(ycf), nrange);
            auto xi_cf =
                cntfrc::make_generator<cntfrc::index_tracker, cntfrc::partial_fraction_tracker,
                                       cntfrc::previous_previous_convergent_tracker>(
                    cntfrc::impl::rational{approx_info.multiplier,
                                           bigint::uint_var::power_of_2(approx_info.shift_amount)});

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

        // Given real numbers x, y, zeta and a range [nmin:nmax] of integers, find a maximizer of
        // (floor(nx+y) - zeta) / n.
        // Precondition: 0 < nmin <= nmax.
        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY,
                  class ContinuedFractionGeneratorZeta>
        constexpr bigint::int_var find_maximizer_of_floor_subtract_quotient_positive_range(
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            ContinuedFractionGeneratorZeta&& zetacf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
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
                cntfrc::make_generator<cntfrc::index_tracker,
                                       cntfrc::previous_previous_convergent_tracker>(
                    zetacf.copy_internal_implementation()),
                util::abs(nrange.upper_bound() - nrange.lower_bound()) + 1u);

            // Find the smallest minimizer of the fractional part.
            auto const n00 = find_extrema_of_fractional_part(xcf, ycf, nrange).smallest_minimizer;

            // Solve the maximization problem on the left.
            auto left_maximizer = n00;
            while (left_maximizer > nrange.lower_bound()) {
                // n1 is the largest minimizer of (floor(nx) + 1) / n, which is the largest multiple
                // of the largest maximizer of nx - floor(nx), where n is in [1:n0 - nmin].
                auto const n1 = [&] {
                    xcf_copy.rewind();
                    auto const nmax = util::abs(left_maximizer - nrange.lower_bound());
                    auto smallest_minimizer =
                        find_extrema_of_fractional_part(xcf_copy, nmax).largest_maximizer;
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
                        find_extrema_of_fractional_part(xcf_copy, nmax).smallest_minimizer;
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
                cntfrc::make_generator<cntfrc::index_tracker,
                                       cntfrc::previous_previous_convergent_tracker>(
                    zetacf.copy_internal_implementation()),
                util::abs(nrange.upper_bound() - nrange.lower_bound()) + 1u);

            // Find the largest maximizer of the fractional part.
            auto const n00 = find_extrema_of_fractional_part(xcf, ycf, nrange).largest_maximizer;

            // Solve the minimization problem on the left.
            auto left_minimizer = n00;
            while (left_minimizer > nrange.lower_bound()) {
                // n1 is the largest maximizer of floor(nx)/ n, which is the largest multiple of the
                // smallest minimizer of nx - floor(nx), where n is in [1:n0 - nmin].
                auto const n1 = [&] {
                    xcf_copy.rewind();
                    auto const nmax = util::abs(left_minimizer - nrange.lower_bound());
                    auto smallest_maximizer =
                        find_extrema_of_fractional_part(xcf_copy, nmax).smallest_minimizer;
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
                        find_extrema_of_fractional_part(xcf_copy, nmax).largest_maximizer;
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
