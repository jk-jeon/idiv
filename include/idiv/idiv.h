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
#include <algorithm>
#include <ranges>

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
        constexpr multiply_shift_info find_optimal_multiply_shift(RationalInterval const& vitv) {
            return vitv.visit([](auto&& itv) -> multiply_shift_info {
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

            // TODO: deal with possible rational dependence between x and y.

            auto xcf_copy = xcf.copy();

            auto const& nmin = nrange.lower_bound();
            auto const nlength = util::abs(nrange.upper_bound() - nrange.lower_bound());

            // Step 1. Subtract out the integer part of y.
            auto floor_y = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper{xcf.copy(),
                                                ycf.copy(),
                                                {// numerator
                                                 0, nmin, 1, 0,
                                                 // denominator
                                                 0, 0, 0, 1}});
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
                    cntfrc::impl::binary_gosper{
                        xcf.copy(),
                        ycf.copy(),
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
                    cntfrc::impl::binary_gosper{
                        xcf.copy(),
                        ycf.copy(),
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
                        cntfrc::impl::binary_gosper{
                            xcf.copy(),
                            ycf.copy(),
                            {// numerator
                             0, nmin * xi_range.lower_bound().denominator,
                             util::to_signed(xi_range.lower_bound().denominator), 0,
                             // denominator
                             0, 0, 0, 1}});
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
                                cntfrc::impl::binary_gosper{
                                    xcf.copy(),
                                    ycf.copy(),
                                    {// numerator
                                     0, nmin * xi_range.lower_bound().denominator,
                                     util::to_signed(xi_range.lower_bound().denominator),
                                     -floor_qstar_y - l -
                                         xi_range.lower_bound().denominator * floor_y,
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

        // Given real numbers x, y, and a range [nmin:nmax] of integers, find the set of (xi,zeta)
        // such that floor(nx + y) = floor(n xi + zeta) holds for all n in [nmin:nmax]. The set is
        // represented as the union of nonempty sets of one of the forms
        // {(xi,zeta) in R x J | (a-zeta)/b < xi < (c-zeta)/d},
        // {(xi,zeta) in R x J | (a-zeta)/b < xi <= (c-zeta)/d},
        // {(xi,zeta) in R x J | (a-zeta)/b <= xi < (c-zeta)/d} and
        // {(xi,zeta) in R x J | (a-zeta)/b <= xi <= (c-zeta)/d},
        // where J is a nonempty interval, b, d are positive integers, and a, c are integers,
        // and the projection J onto the zeta-component of each of these sets is disjoint from each
        // other. We can in fact assume that J is either open and bounded, a singleton set, or the
        // entire R.
        struct elementary_xi_zeta_region {
            variable_shape_interval<frac<bigint::int_var, bigint::uint_var>,
                                    interval_type_t::bounded_open, interval_type_t::bounded_closed,
                                    interval_type_t::entire>
                zeta_range;

            bigint::int_var xi_left_endpoint_numerator;
            bigint::int_var xi_left_endpoint_denominator;

            bigint::int_var xi_right_endpoint_numerator;
            bigint::int_var xi_right_endpoint_denominator;

            bool xi_left_endpoint_included;
            bool xi_right_endpoint_included;
        };

        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY,
                  std::ranges::range RangeOfIntegerIntervals>
        constexpr std::vector<elementary_xi_zeta_region>
        find_xi_zeta_region(ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
                            RangeOfIntegerIntervals&& range_of_nranges) {
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

            using frac_t = frac<bigint::int_var, bigint::uint_var>;
            using nrange_t = interval<bigint::int_var, interval_type_t::bounded_closed>;

            // Sanity check.
            util::constexpr_assert(!std::ranges::empty(range_of_nranges));

            // Get the reduced form of the number num/den, where num and den are rationals.
            auto get_reduced_quotient = [](auto const& num, auto const& den) {
                auto num_int = num.numerator * den.denominator;
                auto den_int = num.denominator * den.numerator;
                if (util::is_strictly_negative(den_int)) {
                    num_int = -std::move(num_int);
                    den_int = -std::move(den_int);
                }
                return find_best_rational_approx(
                           cntfrc::make_generator<cntfrc::index_tracker,
                                                  cntfrc::previous_previous_convergent_tracker>(
                               cntfrc::impl::rational{num_int, util::abs(den_int)}),
                           util::abs(den_int))
                    .below;
            };


            ////////////////////////////////////////////////////////////////////////////////////
            // Step 1 - Write the region as the intersection of half-spaces.
            ////////////////////////////////////////////////////////////////////////////////////

            enum class elementary_problem_sign : bool { positive, negative };
            struct half_space_info {
                frac_t zeta_coeff;
                frac_t eta_coeff;
                enum class boundary_type_t : bool { inclusive, exclusive } boundary_type;
            };
            std::vector<half_space_info> right_half_spaces; // Lower bounds for xi.
            std::vector<half_space_info> left_half_spaces;  // Upper bounds for xi.
            bool nrange_contains_zero = false;
            bigint::int_var floor_y; // Used in Step 4.

            {
                // First, rewrite the domain into a disjoint union of intervals.
                std::vector<nrange_t> normalized_nranges(
                    std::forward<RangeOfIntegerIntervals>(range_of_nranges));
                std::ranges::sort(normalized_nranges, [](auto&& lhs, auto&& rhs) {
                    return lhs.lower_bound() < rhs.lower_bound();
                });

                auto current_pos = normalized_nranges.begin();
                for (auto merge_target_pos = current_pos + 1;
                     merge_target_pos < normalized_nranges.end(); ++merge_target_pos) {
                    // Merge if possible.
                    if (current_pos->upper_bound() + 1 >= merge_target_pos->lower_bound()) {
                        if (current_pos->upper_bound() < merge_target_pos->upper_bound()) {
                            current_pos->upper_bound() = std::move(merge_target_pos->upper_bound());
                        }
                    }
                    else {
                        ++current_pos;
                        util::constexpr_assert(current_pos <= merge_target_pos);
                        if (current_pos != merge_target_pos) {
                            *current_pos = std::move(*merge_target_pos);
                        }
                    }
                }
                normalized_nranges.erase(++current_pos, normalized_nranges.end());

                // Compute good enough approximations of (x,y) and (-x,y) for future computations.
                auto const approx_plus_x_y_info = find_simultaneous_multiply_add_shift(
                    xcf.copy(), ycf.copy(),
                    nrange_t{util::is_strictly_negative(normalized_nranges.front().lower_bound())
                                 ? 0
                                 : normalized_nranges.front().lower_bound(),
                             util::is_strictly_negative(normalized_nranges.back().upper_bound())
                                 ? 0
                                 : normalized_nranges.back().upper_bound()});

                auto const approx_minus_x_y_info = find_simultaneous_multiply_add_shift(
                    cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                           cntfrc::interval_tracker>(
                        cntfrc::impl::unary_gosper{xcf, {-1, 0, 0, 1}}),
                    ycf,
                    nrange_t{util::is_strictly_positive(normalized_nranges.back().upper_bound())
                                 ? 0
                                 : -normalized_nranges.back().upper_bound(),
                             util::is_strictly_positive(normalized_nranges.front().lower_bound())
                                 ? 0
                                 : -normalized_nranges.front().lower_bound()});

                auto xcf_plus_side = cntfrc::make_caching_generator(
                    cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                           cntfrc::interval_tracker>(cntfrc::impl::rational{
                        approx_plus_x_y_info.multiplier,
                        bigint::uint_var::power_of_2(approx_plus_x_y_info.shift_amount)}));

                auto ycf_plus_side = cntfrc::make_caching_generator(
                    cntfrc::make_generator<cntfrc::interval_tracker>(cntfrc::impl::rational{
                        approx_plus_x_y_info.adder,
                        bigint::uint_var::power_of_2(approx_plus_x_y_info.shift_amount)}));

                auto xcf_minus_side = cntfrc::make_caching_generator(
                    cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                           cntfrc::interval_tracker>(cntfrc::impl::rational{
                        approx_minus_x_y_info.multiplier,
                        bigint::uint_var::power_of_2(approx_minus_x_y_info.shift_amount)}));

                auto ycf_minus_side = cntfrc::make_caching_generator(
                    cntfrc::make_generator<cntfrc::interval_tracker>(cntfrc::impl::rational{
                        approx_minus_x_y_info.adder,
                        bigint::uint_var::power_of_2(approx_minus_x_y_info.shift_amount)}));

                floor_y = (approx_plus_x_y_info.adder >> approx_plus_x_y_info.shift_amount);

                // For each disjoint component, decompose it further into elementary domains.
                // For each elementary problem, find all half-spaces that determine the solution
                // region.
                {
                    // For the maximization on the left, the half-space at the base point is
                    // included.
                    auto solve_maximization_on_left = [&](bigint::uint_var base_point,
                                                          bigint::uint_var max_diff,
                                                          elementary_problem_sign sign) {
                        auto const& approx_x_y_info = sign == elementary_problem_sign::positive
                                                          ? approx_plus_x_y_info
                                                          : approx_minus_x_y_info;
                        auto& approx_xcf = sign == elementary_problem_sign::positive
                                               ? xcf_plus_side
                                               : xcf_minus_side;
                        auto& half_spaces = sign == elementary_problem_sign::positive
                                                ? right_half_spaces
                                                : left_half_spaces;

                        while (!util::is_zero(max_diff)) {
                            half_spaces.push_back(
                                {frac_t{1, base_point},
                                 frac_t{-((base_point * approx_x_y_info.multiplier +
                                           approx_x_y_info.adder) >>
                                          approx_x_y_info.shift_amount),
                                        base_point},
                                 half_space_info::boundary_type_t::inclusive});

                            auto movement = find_extrema_of_fractional_part(approx_xcf, max_diff)
                                                .largest_maximizer;
                            approx_xcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point -= movement;
                            max_diff -= std::move(movement);
                        }

                        half_spaces.push_back({frac_t{1, base_point},
                                               frac_t{-((base_point * approx_x_y_info.multiplier +
                                                         approx_x_y_info.adder) >>
                                                        approx_x_y_info.shift_amount),
                                                      base_point},
                                               half_space_info::boundary_type_t::inclusive});
                    };
                    // For the maximization on the right, the half-space at the base point is not
                    // included.
                    auto solve_maximization_on_right = [&](bigint::uint_var base_point,
                                                           bigint::uint_var max_diff,
                                                           elementary_problem_sign sign) {
                        auto const& approx_x_y_info = sign == elementary_problem_sign::positive
                                                          ? approx_plus_x_y_info
                                                          : approx_minus_x_y_info;
                        auto& approx_xcf = sign == elementary_problem_sign::positive
                                               ? xcf_plus_side
                                               : xcf_minus_side;
                        auto& half_spaces = sign == elementary_problem_sign::positive
                                                ? right_half_spaces
                                                : left_half_spaces;

                        while (!util::is_zero(max_diff)) {
                            auto movement = find_extrema_of_fractional_part(approx_xcf, max_diff)
                                                .smallest_minimizer;
                            approx_xcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point += movement;
                            max_diff -= std::move(movement);

                            half_spaces.push_back(
                                {frac_t{1, base_point},
                                 frac_t{-((base_point * approx_x_y_info.multiplier +
                                           approx_x_y_info.adder) >>
                                          approx_x_y_info.shift_amount),
                                        base_point},
                                 half_space_info::boundary_type_t::inclusive});
                        }
                    };
                    // For the minimization on the left, the half-space at the base point is
                    // included.
                    auto solve_minimization_on_left = [&](bigint::uint_var base_point,
                                                          bigint::uint_var max_diff,
                                                          elementary_problem_sign sign) {
                        auto const& approx_x_y_info = sign == elementary_problem_sign::positive
                                                          ? approx_plus_x_y_info
                                                          : approx_minus_x_y_info;
                        auto& approx_xcf = sign == elementary_problem_sign::positive
                                               ? xcf_plus_side
                                               : xcf_minus_side;
                        auto& half_spaces = sign == elementary_problem_sign::positive
                                                ? left_half_spaces
                                                : right_half_spaces;

                        while (!util::is_zero(max_diff)) {
                            half_spaces.push_back(
                                {frac_t{1, base_point},
                                 frac_t{1 - ((base_point * approx_x_y_info.multiplier +
                                              approx_x_y_info.adder) >>
                                             approx_x_y_info.shift_amount),
                                        base_point},
                                 half_space_info::boundary_type_t::exclusive});

                            auto movement = find_extrema_of_fractional_part(approx_xcf, max_diff)
                                                .smallest_minimizer;
                            approx_xcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point -= movement;
                            max_diff -= std::move(movement);
                        }

                        half_spaces.push_back(
                            {frac_t{1, base_point},
                             frac_t{1 - ((base_point * approx_x_y_info.multiplier +
                                          approx_x_y_info.adder) >>
                                         approx_x_y_info.shift_amount),
                                    base_point},
                             half_space_info::boundary_type_t::exclusive});
                    };
                    // For the minimization on the right, the half-space at the base point is not
                    // included.
                    auto solve_minimization_on_right = [&](bigint::uint_var base_point,
                                                           bigint::uint_var max_diff,
                                                           elementary_problem_sign sign) {
                        auto const& approx_x_y_info = sign == elementary_problem_sign::positive
                                                          ? approx_plus_x_y_info
                                                          : approx_minus_x_y_info;
                        auto& approx_xcf = sign == elementary_problem_sign::positive
                                               ? xcf_plus_side
                                               : xcf_minus_side;
                        auto& half_spaces = sign == elementary_problem_sign::positive
                                                ? left_half_spaces
                                                : right_half_spaces;

                        while (!util::is_zero(max_diff)) {
                            auto movement = find_extrema_of_fractional_part(approx_xcf, max_diff)
                                                .largest_maximizer;
                            approx_xcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point += movement;
                            max_diff -= std::move(movement);

                            half_spaces.push_back(
                                {frac_t{1, base_point},
                                 frac_t{1 - ((base_point * approx_x_y_info.multiplier +
                                              approx_x_y_info.adder) >>
                                             approx_x_y_info.shift_amount),
                                        base_point},
                                 half_space_info::boundary_type_t::exclusive});
                        }
                    };

                    auto process_single_sign_interval = [&](nrange_t const& nrange,
                                                            elementary_problem_sign sign,
                                                            auto&& pm_xcf, auto&& pm_ycf) {
                        auto base_points = find_extrema_of_fractional_part(pm_xcf, pm_ycf, nrange);
                        pm_xcf.rewind();
                        pm_ycf.rewind();

                        solve_maximization_on_left(
                            util::abs(base_points.smallest_minimizer),
                            util::abs(base_points.smallest_minimizer - nrange.lower_bound()), sign);
                        solve_maximization_on_right(
                            util::abs(base_points.smallest_minimizer),
                            util::abs(nrange.upper_bound() - base_points.smallest_minimizer), sign);
                        solve_minimization_on_left(
                            util::abs(base_points.largest_maximizer),
                            util::abs(base_points.largest_maximizer - nrange.lower_bound()), sign);
                        solve_minimization_on_right(
                            util::abs(base_points.largest_maximizer),
                            util::abs(nrange.upper_bound() - base_points.largest_maximizer), sign);
                    };

                    for (auto const& nrange : normalized_nranges) {
                        // Negative interval.
                        if (util::is_strictly_negative(nrange.upper_bound())) {
                            process_single_sign_interval(
                                nrange_t{-nrange.upper_bound(), -nrange.lower_bound()},
                                elementary_problem_sign::negative, xcf_minus_side, ycf_minus_side);
                        }
                        // Positive interval.
                        else if (util::is_strictly_positive(nrange.lower_bound())) {
                            process_single_sign_interval(nrange, elementary_problem_sign::positive,
                                                         xcf_plus_side, ycf_plus_side);
                        }
                        else if (util::is_zero(nrange.upper_bound())) {
                            nrange_contains_zero = true;
                            process_single_sign_interval(nrange_t{1, -nrange.lower_bound()},
                                                         elementary_problem_sign::negative,
                                                         xcf_minus_side, ycf_minus_side);
                        }
                        else if (util::is_zero(nrange.lower_bound())) {
                            nrange_contains_zero = true;
                            process_single_sign_interval(nrange_t{1, nrange.upper_bound()},
                                                         elementary_problem_sign::positive,
                                                         xcf_plus_side, ycf_plus_side);
                        }
                        else {
                            // nrange.lower_bound() < 0 < nrange.upper_bound()
                            nrange_contains_zero = true;
                            process_single_sign_interval(nrange_t{1, -nrange.lower_bound()},
                                                         elementary_problem_sign::negative,
                                                         xcf_minus_side, ycf_minus_side);
                            process_single_sign_interval(nrange_t{1, nrange.upper_bound()},
                                                         elementary_problem_sign::positive,
                                                         xcf_plus_side, ycf_plus_side);
                        }
                    }
                }
            }


            ////////////////////////////////////////////////////////////////////////////////////
            // Step 2 - Find the intersection of lower bounds for xi and upper bounds for xi,
            // respectively. The intersection is described in terms of its horizontal slices.
            // To simplify further processing, we make all the vertical projections of these
            // slices to be either open or singleton.
            ////////////////////////////////////////////////////////////////////////////////////

            struct elementary_one_sided_region {
                variable_shape_interval<
                    frac_t, interval_type_t::bounded_open, interval_type_t::bounded_closed,
                    interval_type_t::bounded_below_open, interval_type_t::bounded_above_open,
                    interval_type_t::entire>
                    zeta_range;

                bigint::int_var xi_endpoint_numerator;
                bigint::int_var xi_endpoint_denominator;
                bool xi_endpoint_included;
            };

            // Find the extreme points for the lower/upper bounds by finding the convex hull of the
            // dual problem projected onto the plane xi = +-1.
            enum class bound_direction_t : bool { lower, upper };
            auto compute_intersection = [&](bound_direction_t bound_direction) {
                struct vec2d {
                    frac_t zeta_coord;
                    frac_t eta_coord;

                    constexpr frac_t dot(vec2d const& other) const {
                        return zeta_coord * other.zeta_coord + eta_coord * other.eta_coord;
                    }
                    constexpr frac_t normsq() const { return dot(*this); }
                };

                auto& half_spaces = bound_direction == bound_direction_t::lower ? right_half_spaces
                                                                                : left_half_spaces;
                std::vector<elementary_one_sided_region> result;

                util::constexpr_assert(!half_spaces.empty());

                // Start from the one with the smallest/largest zeta-coordinate, depending on the
                // bound direction. This corresponds to the half-space with the lowest boundary
                // line.
                auto first_elmt = bound_direction == bound_direction_t::lower
                                      ? std::ranges::min_element(std::as_const(half_spaces), {},
                                                                 &half_space_info::zeta_coeff)
                                      : std::ranges::max_element(std::as_const(half_spaces), {},
                                                                 &half_space_info::zeta_coeff);

                // If there is only one half-space, return immediately.
                if (half_spaces.size() == 1) {
                    result.push_back(
                        {interval<frac_t, interval_type_t::entire>{},
                         util::is_nonnegative(first_elmt->zeta_coeff.numerator)
                             ? -first_elmt->eta_coeff.numerator
                             : first_elmt->eta_coeff.numerator,
                         util::is_nonnegative(first_elmt->zeta_coeff.numerator)
                             ? -util::to_signed(first_elmt->zeta_coeff.denominator)
                             : util::to_signed(first_elmt->zeta_coeff.denominator),
                         first_elmt->boundary_type == half_space_info::boundary_type_t::inclusive});
                    return result;
                }

                // We are at the left/right-end and we want to travel counterclockwise, when viewed
                // from the positive xi-axis. To do so, we set the initial direction to be along the
                // negative/positive eta-axis, depending on the direction of the bound.
                auto prev_direction_vec =
                    vec2d{frac_t{0, 1u},
                          frac_t{bound_direction == bound_direction_t::lower ? -1 : 1, 1u}};
                auto last_elmt = first_elmt;
                frac_t prev_turning_point_zeta{0, 1u};

                struct angle_info {
                    typename std::vector<half_space_info>::const_iterator itr;
                    vec2d direction_vec;
                    frac_t cos_square;
                    bool is_cos_strictly_negative;
                    bool is_inclusive_at_turning_point;
                };
                auto compute_angle_info = [&](auto itr) {
                    auto direction_vec = vec2d{itr->zeta_coeff - last_elmt->zeta_coeff,
                                               itr->eta_coeff - last_elmt->eta_coeff};
                    auto dot_product = prev_direction_vec.dot(direction_vec);
                    bool is_cos_strictly_negative =
                        util::is_strictly_negative(dot_product.numerator);
                    return angle_info{
                        itr, direction_vec, (dot_product * dot_product) / direction_vec.normsq(),
                        is_cos_strictly_negative,
                        last_elmt->boundary_type == half_space_info::boundary_type_t::inclusive &&
                            itr->boundary_type == half_space_info::boundary_type_t::inclusive};
                };
                auto compare_angle_info = [](angle_info const& left,
                                             angle_info const& right) -> std::strong_ordering {
                    if (left.is_cos_strictly_negative) {
                        if (right.is_cos_strictly_negative) {
                            return right.cos_square <=> left.cos_square;
                        }
                        else {
                            return std::strong_ordering::less;
                        }
                    }
                    else {
                        if (right.is_cos_strictly_negative) {
                            return std::strong_ordering::greater;
                        }
                        else {
                            return left.cos_square <=> right.cos_square;
                        }
                    }
                };

                while (true) {
                    // Find the point whose direction vector is the closest in angle to the previous
                    // direction vector.
                    auto itr = half_spaces.cbegin();
                    if (itr == last_elmt) {
                        ++itr;
                    }
                    auto current_angle_info = compute_angle_info(itr);

                    for (++itr; itr != half_spaces.cend(); ++itr) {
                        if (itr == last_elmt) {
                            continue;
                        }

                        auto new_angle_info = compute_angle_info(itr);
                        auto const compare_result =
                            compare_angle_info(current_angle_info, new_angle_info);

                        // If current < new, found a better one.
                        if (compare_result < 0) {
                            current_angle_info = std::move(new_angle_info);
                        }
                        // If current == new.
                        else if (compare_result == 0) {
                            // We choose the one located further, but take account that the
                            // inclusivity of the turning point may change.
                            bool const is_inclusive_at_turning_point =
                                (current_angle_info.is_inclusive_at_turning_point &&
                                 new_angle_info.is_inclusive_at_turning_point);
                            if (current_angle_info.direction_vec.normsq() <
                                new_angle_info.direction_vec.normsq()) {
                                current_angle_info = std::move(new_angle_info);
                            }
                            current_angle_info.is_inclusive_at_turning_point =
                                is_inclusive_at_turning_point;
                        }
                    }
                    // Found one.

                    // Find a normal vector to the newly found face.
                    // For the case of lower bound, every functional should have nonnegative inner
                    // product with this normal, while for the case of upper bound, they should have
                    // nonpositive inner product.
                    auto face_normal_xi =
                        last_elmt->zeta_coeff * current_angle_info.itr->eta_coeff -
                        last_elmt->eta_coeff * current_angle_info.itr->zeta_coeff;
                    auto face_normal_zeta =
                        last_elmt->eta_coeff - current_angle_info.itr->eta_coeff;
                    auto face_normal_eta =
                        current_angle_info.itr->zeta_coeff - last_elmt->zeta_coeff;

                    // Project it down to the plane eta = 1. The resulting point must be a turning
                    // point.
                    auto turning_point_zeta =
                        get_reduced_quotient(face_normal_zeta, face_normal_eta);

                    auto xi_endpoint_numerator =
                        util::is_nonnegative(last_elmt->zeta_coeff.numerator)
                            ? -last_elmt->eta_coeff.numerator
                            : last_elmt->eta_coeff.numerator;
                    auto xi_endpoint_denominator =
                        util::is_nonnegative(last_elmt->zeta_coeff.numerator)
                            ? -util::to_signed(last_elmt->zeta_coeff.denominator)
                            : util::to_signed(last_elmt->zeta_coeff.denominator);

                    // The unbounded bottom boundary line.
                    if (last_elmt == first_elmt) {
                        util::constexpr_assert(
                            (bound_direction == bound_direction_t::lower &&
                             util::is_strictly_positive(face_normal_eta.numerator)) ||
                            (bound_direction == bound_direction_t::upper &&
                             util::is_strictly_negative(face_normal_eta.numerator)));

                        // The unbounded open region.
                        result.push_back({interval<frac_t, interval_type_t::bounded_above_open>{
                                              turning_point_zeta},
                                          xi_endpoint_numerator, xi_endpoint_denominator,
                                          last_elmt->boundary_type ==
                                              half_space_info::boundary_type_t::inclusive});

                        // The horizontal ray right above it.
                        result.push_back({interval<frac_t, interval_type_t::bounded_closed>{
                                              turning_point_zeta, turning_point_zeta},
                                          xi_endpoint_numerator, xi_endpoint_denominator,
                                          last_elmt->boundary_type ==
                                                  half_space_info::boundary_type_t::inclusive &&
                                              current_angle_info.itr->boundary_type ==
                                                  half_space_info::boundary_type_t::inclusive});
                    }
                    // The unbounded top boundary line.
                    else if (current_angle_info.itr == first_elmt) {
                        util::constexpr_assert(
                            (bound_direction == bound_direction_t::lower &&
                             util::is_strictly_negative(face_normal_eta.numerator)) ||
                            (bound_direction == bound_direction_t::upper &&
                             util::is_strictly_positive(face_normal_eta.numerator)));

                        // The unbounded open region.
                        result.push_back({interval<frac_t, interval_type_t::bounded_below_open>{
                                              prev_turning_point_zeta},
                                          xi_endpoint_numerator, xi_endpoint_denominator,
                                          last_elmt->boundary_type ==
                                              half_space_info::boundary_type_t::inclusive});

                        // End of the while (true) {...} loop.
                        break;
                    }
                    // Bounded middle boundary lines.
                    else {
                        util::constexpr_assert(
                            (bound_direction == bound_direction_t::lower &&
                             util::is_strictly_positive(face_normal_eta.numerator)) ||
                            (bound_direction == bound_direction_t::upper &&
                             util::is_strictly_negative(face_normal_eta.numerator)));

                        // The bounded open region.
                        result.push_back({interval<frac_t, interval_type_t::bounded_open>{
                                              prev_turning_point_zeta, turning_point_zeta},
                                          xi_endpoint_numerator, xi_endpoint_denominator,
                                          last_elmt->boundary_type ==
                                              half_space_info::boundary_type_t::inclusive});

                        // The horizontal ray right above it.
                        result.push_back({interval<frac_t, interval_type_t::bounded_closed>{
                                              turning_point_zeta, turning_point_zeta},
                                          xi_endpoint_numerator, xi_endpoint_denominator,
                                          last_elmt->boundary_type ==
                                                  half_space_info::boundary_type_t::inclusive &&
                                              current_angle_info.itr->boundary_type ==
                                                  half_space_info::boundary_type_t::inclusive});
                    }
                    // End of the branching on top/middle/bottom regions.

                    prev_direction_vec = current_angle_info.direction_vec;
                    last_elmt = current_angle_info.itr;
                    prev_turning_point_zeta = turning_point_zeta;
                } // while (true)

                return result;
            };


            ////////////////////////////////////////////////////////////////////////////////////
            // Step 3 - Find the intersection of the region from the lower bounds and the region
            // from the upper bounds.
            ////////////////////////////////////////////////////////////////////////////////////

            std::vector<elementary_xi_zeta_region> result;
            {
                auto lower_bound_region = compute_intersection(bound_direction_t::lower);
                auto upper_bound_region = compute_intersection(bound_direction_t::upper);

                // Sweep from below to above.
                // By the construction, these arrays should be sorted according to the
                // zeta-coordinate, from below to above.
                util::constexpr_assert(!lower_bound_region.empty() && !upper_bound_region.empty());
                auto lower_bound_itr = lower_bound_region.cbegin();
                auto upper_bound_itr = upper_bound_region.cbegin();

                // Except for the unique exceptional case of having only one n, these two bounds
                // must meet at exactly two points.
                if (lower_bound_region.size() == 1) {
                    // For the said exceptional case, we should have an infinite parallelogram.
                    util::constexpr_assert(upper_bound_region.size() == 1);
                    util::constexpr_assert(lower_bound_itr->zeta_range.interval_type() ==
                                           interval_type_t::entire);
                    util::constexpr_assert(upper_bound_itr->zeta_range.interval_type() ==
                                           interval_type_t::entire);
                    util::constexpr_assert(lower_bound_itr->xi_endpoint_denominator ==
                                           upper_bound_itr->xi_endpoint_denominator);
                    util::constexpr_assert(
                        lower_bound_itr->xi_endpoint_numerator +
                            (util::is_strictly_positive(lower_bound_itr->xi_endpoint_denominator)
                                 ? 1
                                 : -1) ==
                        upper_bound_itr->xi_endpoint_numerator);

                    result.push_back({interval<frac_t, interval_type_t::entire>{},
                                      lower_bound_itr->xi_endpoint_numerator,
                                      lower_bound_itr->xi_endpoint_denominator,
                                      upper_bound_itr->xi_endpoint_numerator,
                                      upper_bound_itr->xi_endpoint_denominator,
                                      lower_bound_itr->xi_endpoint_included,
                                      upper_bound_itr->xi_endpoint_included});
                }
                else {
                    frac_t previous_zeta_endpoint{0, 0u};

                    auto push_new_region = [&] {
                        auto xi_left_endpoint_numerator =
                            std::move(lower_bound_itr->xi_endpoint_numerator);
                        auto xi_left_endpoint_denominator =
                            std::move(lower_bound_itr->xi_endpoint_denominator);
                        auto xi_right_endpoint_numerator =
                            std::move(upper_bound_itr->xi_endpoint_numerator);
                        auto xi_right_endpoint_denominator =
                            std::move(upper_bound_itr->xi_endpoint_denominator);
                        auto xi_left_endpoint_included = lower_bound_itr->xi_endpoint_included;
                        auto xi_right_endpoint_included = upper_bound_itr->xi_endpoint_included;

                        // Compare two zeta ranges, and choose the one with the smaller
                        // right endpoint to construct the region.
                        auto zeta_range =
                            [&]() -> variable_shape_interval<frac_t, interval_type_t::bounded_open,
                                                             interval_type_t::bounded_closed> {
                            // If any of the lower bound and the upper bound is a horizontal ray.
                            if (lower_bound_itr->zeta_range.interval_type() ==
                                interval_type_t::bounded_closed) {
                                ++lower_bound_itr;
                                if (upper_bound_itr->zeta_range.interval_type() ==
                                    interval_type_t::bounded_closed) {
                                    ++upper_bound_itr;
                                }
                                return interval<frac_t, interval_type_t::bounded_closed>{
                                    previous_zeta_endpoint, previous_zeta_endpoint};
                            }
                            else if (upper_bound_itr->zeta_range.interval_type() ==
                                     interval_type_t::bounded_closed) {
                                ++upper_bound_itr;
                                return interval<frac_t, interval_type_t::bounded_closed>{
                                    previous_zeta_endpoint, previous_zeta_endpoint};
                            }

                            return lower_bound_itr->zeta_range.with_upper_bound(
                                [&](auto const& ub1) -> variable_shape_interval<
                                                         frac_t, interval_type_t::bounded_open,
                                                         interval_type_t::bounded_closed> {
                                    return upper_bound_itr->zeta_range.with_upper_bound(
                                        [&](auto const& ub2) {
                                            auto cmp_result = ub1 <=> ub2;
                                            if (cmp_result <= 0) {
                                                ++lower_bound_itr;
                                            }
                                            if (cmp_result >= 0) {
                                                ++upper_bound_itr;
                                            }
                                            auto ret_value =
                                                interval<frac_t, interval_type_t::bounded_open>{
                                                    std::move(previous_zeta_endpoint),
                                                    cmp_result <= 0 ? ub1 : ub2};
                                            previous_zeta_endpoint = ret_value.upper_bound();

                                            return ret_value;
                                        },
                                        [&] {
                                            // Upper bound has unbounded zeta interval.
                                            ++lower_bound_itr;
                                            auto ret_value =
                                                interval<frac_t, interval_type_t::bounded_open>{
                                                    std::move(previous_zeta_endpoint), ub1};
                                            previous_zeta_endpoint = ub1;
                                            return ret_value;
                                        });
                                },
                                // Lower bound has unbounded zeta interval.
                                [&]() -> variable_shape_interval<frac_t,
                                                                 interval_type_t::bounded_open,
                                                                 interval_type_t::bounded_closed> {
                                    return upper_bound_itr->zeta_range.with_upper_bound(
                                        [&](auto const& ub2)
                                            -> variable_shape_interval<
                                                frac_t, interval_type_t::bounded_open,
                                                interval_type_t::bounded_closed> {
                                            ++upper_bound_itr;
                                            auto ret_value =
                                                interval<frac_t, interval_type_t::bounded_open>{
                                                    std::move(previous_zeta_endpoint), ub2};
                                            previous_zeta_endpoint = ub2;
                                            return ret_value;
                                        },
                                        [&]() -> variable_shape_interval<
                                                  frac_t, interval_type_t::bounded_open,
                                                  interval_type_t::bounded_closed> {
                                            // Impossible to reach here.
                                            return interval<frac_t,
                                                            interval_type_t::bounded_closed>{
                                                previous_zeta_endpoint, previous_zeta_endpoint};
                                        });
                                });
                        }();

                        result.push_back({std::move(zeta_range),
                                          std::move(xi_left_endpoint_numerator),
                                          std::move(xi_left_endpoint_denominator),
                                          std::move(xi_right_endpoint_numerator),
                                          std::move(xi_right_endpoint_denominator),
                                          xi_left_endpoint_included, xi_right_endpoint_included});
                    };

                    bool found_first_intersection = false;
                    while (true) {
                        // If two lines are not parallel.
                        if (lower_bound_itr->xi_endpoint_denominator !=
                            upper_bound_itr->xi_endpoint_denominator) {
                            auto intersection_zeta = get_reduced_quotient(
                                frac{lower_bound_itr->xi_endpoint_numerator *
                                             upper_bound_itr->xi_endpoint_denominator -
                                         upper_bound_itr->xi_endpoint_numerator *
                                             lower_bound_itr->xi_endpoint_denominator,
                                     cntfrc::unity{}},
                                frac{upper_bound_itr->xi_endpoint_denominator -
                                         lower_bound_itr->xi_endpoint_denominator,
                                     cntfrc::unity{}});

                            // If zeta is in the range, we found an intersection.
                            if (lower_bound_itr->zeta_range.contains(intersection_zeta) &&
                                upper_bound_itr->zeta_range.contains(intersection_zeta)) {
                                // When this is the first intersection, start pushing regions.
                                if (!found_first_intersection) {
                                    // push_new_region() will not push the bottom vertex if both of
                                    // the zeta-intervals are open.
                                    if (lower_bound_itr->zeta_range.interval_type() !=
                                            interval_type_t::bounded_closed &&
                                        lower_bound_itr->zeta_range.interval_type() !=
                                            interval_type_t::bounded_closed) {
                                        result.push_back(
                                            {interval<frac_t, interval_type_t::bounded_closed>{
                                                 intersection_zeta, intersection_zeta},
                                             lower_bound_itr->xi_endpoint_numerator,
                                             lower_bound_itr->xi_endpoint_denominator,
                                             upper_bound_itr->xi_endpoint_numerator,
                                             upper_bound_itr->xi_endpoint_denominator,
                                             lower_bound_itr->xi_endpoint_included,
                                             upper_bound_itr->xi_endpoint_included});
                                    }

                                    previous_zeta_endpoint = intersection_zeta;
                                    push_new_region();
                                    found_first_intersection = true;
                                    continue;
                                }
                                // When this is the second intersection, we completed the process.
                                else {
                                    // When the intersection is found between open intervals, then
                                    // add both the open triangular region and the vertex.
                                    // Otherwise, just add the vertex.
                                    if (lower_bound_itr->zeta_range.interval_type() !=
                                            interval_type_t::bounded_closed &&
                                        lower_bound_itr->zeta_range.interval_type() !=
                                            interval_type_t::bounded_closed) {
                                        result.push_back(
                                            {interval<frac_t, interval_type_t::bounded_open>{
                                                 previous_zeta_endpoint, intersection_zeta},
                                             lower_bound_itr->xi_endpoint_numerator,
                                             lower_bound_itr->xi_endpoint_denominator,
                                             upper_bound_itr->xi_endpoint_numerator,
                                             upper_bound_itr->xi_endpoint_denominator,
                                             lower_bound_itr->xi_endpoint_included,
                                             upper_bound_itr->xi_endpoint_included});
                                    }

                                    result.push_back(
                                        {interval<frac_t, interval_type_t::bounded_closed>{
                                             intersection_zeta, intersection_zeta},
                                         lower_bound_itr->xi_endpoint_numerator,
                                         lower_bound_itr->xi_endpoint_denominator,
                                         upper_bound_itr->xi_endpoint_numerator,
                                         upper_bound_itr->xi_endpoint_denominator,
                                         lower_bound_itr->xi_endpoint_included,
                                         upper_bound_itr->xi_endpoint_included});

                                    // End of the while (true) {...} loop.
                                    break;
                                }
                            }
                        }

                        // Found no intersection.
                        // Add a new region if we already have found the first intersection.
                        if (found_first_intersection) {
                            push_new_region();
                        }
                        else {
                            // Choose the one with smaller upper bound and move to the next region.
                            lower_bound_itr->zeta_range.with_upper_bound(
                                [&](auto const& ub1) {
                                    return upper_bound_itr->zeta_range.with_upper_bound(
                                        [&](auto const& ub2) {
                                            auto cmp_result = ub1 <=> ub2;
                                            if (cmp_result <= 0) {
                                                ++lower_bound_itr;
                                            }
                                            if (cmp_result >= 0) {
                                                ++upper_bound_itr;
                                            }
                                        },
                                        [&] {
                                            // Upper bound has unbounded zeta interval.
                                            ++lower_bound_itr;
                                        });
                                },
                                // Lower bound has unbounded zeta interval.
                                [&] {
                                    return upper_bound_itr->zeta_range.with_upper_bound(
                                        [&](auto const& ub2) { ++upper_bound_itr; },
                                        [&] {
                                            // Impossible to reach here.
                                        });
                                });
                        }
                    } // while (true)
                }     // End of the branching on lower_bound_region.size().
            }


            ////////////////////////////////////////////////////////////////////////////////////
            // Step 4 - Horizontally cut the region if the constraint from n = 0 is present.
            ////////////////////////////////////////////////////////////////////////////////////

            if (nrange_contains_zero) {
                // zeta should satisfy the inequality floor_y <= zeta < floor_y + 1.
                auto floor_y_frac = frac_t{floor_y, 1u};
                auto floor_y_p1_frac = frac_t{floor_y + 1, 1u};
                std::vector<elementary_xi_zeta_region> trimmed;
                auto src_itr = result.begin();

                // Find floor_y.
                for (; src_itr != result.end(); ++src_itr) {
                    if (src_itr->zeta_range.contains(floor_y_frac)) {
                        bool also_contains_p1 = src_itr->zeta_range.contains(floor_y_p1_frac);

                        if (src_itr->zeta_range.interval_type() == interval_type_t::bounded_open) {
                            // Split the interval.
                            trimmed.push_back({interval<frac_t, interval_type_t::bounded_closed>{
                                                   floor_y_frac, floor_y_frac},
                                               src_itr->xi_left_endpoint_numerator,
                                               src_itr->xi_left_endpoint_denominator,
                                               src_itr->xi_right_endpoint_numerator,
                                               src_itr->xi_right_endpoint_denominator,
                                               src_itr->xi_left_endpoint_included,
                                               src_itr->xi_right_endpoint_included});

                            src_itr->zeta_range = src_itr->zeta_range.with_upper_bound(
                                [&](auto const& ub) {
                                    return interval<frac_t, interval_type_t::bounded_open>{
                                        floor_y_frac, also_contains_p1 ? floor_y_p1_frac : ub};
                                },
                                [&] {
                                    // Cannot reach here.
                                    return interval<frac_t, interval_type_t::bounded_open>{
                                        floor_y_frac, floor_y_p1_frac};
                                });
                        }
                        trimmed.push_back(std::move(*src_itr));

                        if (also_contains_p1) {
                            trimmed.push_back(
                                {interval<frac_t, interval_type_t::bounded_closed>{floor_y_p1_frac,
                                                                                   floor_y_p1_frac},
                                 trimmed.back().xi_left_endpoint_numerator,
                                 trimmed.back().xi_left_endpoint_denominator,
                                 trimmed.back().xi_right_endpoint_numerator,
                                 trimmed.back().xi_right_endpoint_denominator, false, false});

                            src_itr = result.end();
                        }
                        else {
                            ++src_itr;
                        }
                        break;
                    }
                } // End of the for loop.

                // Find floor_y + 1.
                for (; src_itr != result.end(); ++src_itr) {
                    if (src_itr->zeta_range.contains(floor_y_p1_frac)) {
                        if (src_itr->zeta_range.interval_type() == interval_type_t::bounded_open) {
                            // Split the interval.
                            src_itr->zeta_range = src_itr->zeta_range.with_lower_bound(
                                [&](auto const& lb) {
                                    return interval<frac_t, interval_type_t::bounded_open>{
                                        lb, floor_y_p1_frac};
                                },
                                [&] {
                                    // Cannot reach here.
                                    return interval<frac_t, interval_type_t::bounded_open>{
                                        floor_y_frac, floor_y_p1_frac};
                                });

                            trimmed.push_back(*src_itr);
                        }

                        trimmed.push_back(
                            {interval<frac_t, interval_type_t::bounded_closed>{floor_y_p1_frac,
                                                                               floor_y_p1_frac},
                             std::move(src_itr->xi_left_endpoint_numerator),
                             std::move(src_itr->xi_left_endpoint_denominator),
                             std::move(src_itr->xi_right_endpoint_numerator),
                             std::move(src_itr->xi_right_endpoint_denominator), false, false});

                        break;
                    }
                    trimmed.push_back(std::move(*src_itr));
                }

                result = std::move(trimmed);
            } // if (nrange_contains_zero)

            return result;
        }
    }
}

#endif
