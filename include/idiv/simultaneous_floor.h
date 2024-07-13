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

#ifndef JKJ_HEADER_SIMULTANEOUS_FLOOR
#define JKJ_HEADER_SIMULTANEOUS_FLOOR

#include "best_rational_approx.h"
#include "caching_generator.h"
#include "gosper_continued_fraction.h"
#include "bigint.h"

namespace jkj {
    namespace idiv {
        // Represents the set either of the form
        // {(xi, zeta) | xi_min <= xi < xi_max, a - b*xi <= zeta < c - d*xi}
        // or of the form
        // {(xi, zeta) | a - b*xi <= zeta < c - d*xi}.
        struct xi_zeta_trapezoid {
            variable_shape_interval<frac<bigint::int_var, bigint::uint_var>,
                                    interval_type_t::bounded_left_closed_right_open,
                                    interval_type_t::entire>
                xi_range;
            bigint::int_var zeta_left_endpoint_constant_coeff;
            bigint::int_var zeta_left_endpoint_negative_linear_coeff;
            bigint::int_var zeta_right_endpoint_constant_coeff;
            bigint::int_var zeta_right_endpoint_negative_linear_coeff;
        };

        // Given real numbers x, y and a range [nmin:nmax] of integers, find the trapezoidal region
        // of (xi, zeta) such that
        // (1) floor(nx + y) = floor(nxi + zeta) holds for all n in [nmin:nmax], and
        // (2) floor(nx) = floor(nxi) holds for all n in [0:nmax-nmin].
        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY>
        constexpr xi_zeta_trapezoid find_xi_zeta_region_simultaneous_floor(
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorX>::template is_implementing_mixins<
                    cntfrc::convergent_tracker, cntfrc::interval_tracker>(),
                "the first continued fraction generator must implement convergent_tracker and "
                "interval_tracker");
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorY>::template is_implementing_mixins<
                    cntfrc::interval_tracker>(),
                "the second continued fraction generator must implement interval_tracker");

            // TODO: deal with possible rational dependence between x and y.

            using frac_t = frac<bigint::int_var, bigint::uint_var>;
            auto caching_xcf = cntfrc::caching_generator<ContinuedFractionGeneratorX>{
                std::forward<ContinuedFractionGeneratorX>(xcf)};
            auto caching_ycf = cntfrc::caching_generator<ContinuedFractionGeneratorY>{
                std::forward<ContinuedFractionGeneratorY>(ycf)};

            using caching_xcf_ref = decltype(caching_xcf)&;
            using caching_ycf_ref = decltype(caching_ycf)&;

            auto const& nmin = nrange.lower_bound();
            auto const nlength = util::abs(nrange.upper_bound() - nrange.lower_bound());

            // floor(nmin * x + y).
            auto floor_nmin_x_plus_y = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper<caching_xcf_ref, caching_ycf_ref>{caching_xcf,
                                                                                  caching_ycf,
                                                                                  {// numerator
                                                                                   0, nmin, 1, 0,
                                                                                   // denominator
                                                                                   0, 0, 0, 1}});
                cf.update();
                caching_xcf.rewind();
                caching_ycf.rewind();
                return cf.current_partial_fraction().denominator;
            }();

            // Early return for the case nmin == nmax.
            if (util::is_zero(nlength)) {
                return xi_zeta_trapezoid{interval<frac_t, interval_type_t::entire>{},
                                         floor_nmin_x_plus_y, nmin, floor_nmin_x_plus_y + 1u, nmin};
            }

            // Find the range of xi, p_*/q_*, p^*/q^*, and the modular inverse q of -p_* with
            // respect to q_*. Recall that q = 1 if q_* = 1, and it is the denominator of the last
            // odd semiconvergent before p_*/q_*.
            struct x_info_t {
                using convergent_type =
                    typename std::remove_cvref_t<ContinuedFractionGeneratorX>::convergent_type;
                frac_t best_below;
                frac_t best_above;
                frac_t xi_upper_bound;
                bigint::uint_var mod_inv;
            } x_info = [&] {
                auto xi_range = find_floor_quotient_range(caching_xcf, nlength);

                if (xi_range.lower_bound().denominator == 1u) {
                    return x_info_t{xi_range.lower_bound(), xi_range.lower_bound(),
                                    xi_range.upper_bound(), 1u};
                }

                // If q_* = q^*.
                if (caching_xcf.terminated()) {
                    // If x is its even convergent, then the last proper odd semiconvergent is the
                    // previous convergent. If x is its odd convergent, then the denominator of the
                    // last proper odd semiconvergent is the current denominator minust the last
                    // convergent's denominator.
                    return x_info_t{xi_range.lower_bound(), xi_range.lower_bound(),
                                    xi_range.upper_bound(),
                                    caching_xcf.current_index() % 2 == 0
                                        ? caching_xcf.previous_convergent_denominator()
                                        : caching_xcf.current_convergent_denominator() -
                                              caching_xcf.previous_convergent_denominator()};
                }

                // If the last semiconvergent is an even semiconvergent, then the last odd
                // semiconvergent with the denominator < q_* is the previous convergent. If the last
                // semiconvergent is an odd semiconvergent, then the last odd semiconvergent with
                // the denominator < q_* is the previous previous convergent.
                return x_info_t{xi_range.lower_bound(), xi_range.upper_bound(),
                                xi_range.upper_bound(),
                                caching_xcf.current_index() % 2 == 0
                                    ? caching_xcf.previous_convergent_denominator()
                                    : caching_xcf.previous_previous_convergent_denominator()};
            }();
            bool xcf_terminated = caching_xcf.terminated();
            caching_xcf.rewind();

            // See if L is empty.
            {
                auto const lhs = [&] {
                    auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                        cntfrc::impl::binary_gosper<caching_xcf_ref, caching_ycf_ref>{
                            caching_xcf,
                            caching_ycf,
                            {// numerator
                             0, x_info.best_below.denominator + nmin, 1, 0,
                             // denominator
                             0, 0, 0, 1}});
                    cf.update();
                    caching_xcf.rewind();
                    caching_ycf.rewind();
                    return cf.current_partial_fraction().denominator;
                }();
                auto rhs = x_info.best_below.numerator + floor_nmin_x_plus_y;

                if (lhs > rhs) {
                    // L is empty.
                    auto mu = xcf_terminated ? x_info.mod_inv : x_info.best_above.denominator;

                    return xi_zeta_trapezoid{
                        interval<frac_t, interval_type_t::bounded_left_closed_right_open>{
                            x_info.best_below, x_info.xi_upper_bound},
                        rhs + 1u, x_info.best_below.denominator + nmin,
                        util::div_floor(mu * x_info.best_below.numerator,
                                        x_info.best_below.denominator) +
                            floor_nmin_x_plus_y + 2u,
                        mu + nmin};
                }
            }

            auto const floor_qstar_yprime = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper<caching_xcf_ref, caching_ycf_ref>{
                        caching_xcf,
                        caching_ycf,
                        {// numerator
                         0, x_info.best_below.denominator * nmin,
                         util::to_signed(x_info.best_below.denominator), 0,
                         // denominator
                         0, 0, 0, 1}});
                cf.update();
                caching_xcf.rewind();
                caching_ycf.rewind();
                return cf.current_partial_fraction().denominator -
                       x_info.best_below.denominator * floor_nmin_x_plus_y;
            }();

            // See if R is empty.
            {
                if (xcf_terminated) {
                    if (util::is_zero(floor_qstar_yprime)) {
                        // R is empty.
                        return xi_zeta_trapezoid{
                            interval<frac_t, interval_type_t::bounded_left_closed_right_open>{
                                x_info.best_below, x_info.xi_upper_bound},
                            x_info.best_below.numerator + floor_nmin_x_plus_y,
                            x_info.best_below.denominator + nmin,
                            util::div_floor(x_info.mod_inv * x_info.best_below.numerator,
                                            x_info.best_below.denominator) +
                                floor_nmin_x_plus_y + 1u,
                            x_info.mod_inv + nmin};
                    }
                }
                else {
                    auto const lhs = [&] {
                        auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                            cntfrc::impl::binary_gosper<caching_xcf_ref, caching_ycf_ref>{
                                caching_xcf,
                                caching_ycf,
                                {// numerator
                                 0, x_info.best_above.denominator + nmin, 1, 0,
                                 // denominator
                                 0, 0, 0, 1}});
                        cf.update();
                        caching_xcf.rewind();
                        caching_ycf.rewind();
                        return cf.current_partial_fraction().denominator;
                    }();
                    auto rhs = x_info.best_above.numerator + floor_nmin_x_plus_y;

                    if (lhs < rhs) {
                        // R is empty.
                        return xi_zeta_trapezoid{
                            interval<frac_t, interval_type_t::bounded_left_closed_right_open>{
                                x_info.best_below, x_info.xi_upper_bound},
                            x_info.best_below.numerator + floor_nmin_x_plus_y,
                            x_info.best_below.denominator + nmin,
                            util::div_floor(x_info.best_above.denominator *
                                                x_info.best_below.numerator,
                                            x_info.best_below.denominator) +
                                floor_nmin_x_plus_y + 1u,
                            x_info.best_above.denominator + nmin};
                    }
                }
            }

            // Both L & R are not empty but x is rational with denominator <= nlength.
            if (xcf_terminated) {
                auto mu =
                    ((floor_qstar_yprime + 1u) * x_info.mod_inv) % x_info.best_below.denominator;
                auto nu = (floor_qstar_yprime * x_info.mod_inv) % x_info.best_below.denominator;

                return xi_zeta_trapezoid{
                    interval<frac_t, interval_type_t::bounded_left_closed_right_open>{
                        x_info.best_below, x_info.xi_upper_bound},
                    util::div_floor(nu * x_info.best_below.numerator,
                                    x_info.best_below.denominator) +
                        floor_nmin_x_plus_y + 1u,
                    nu + nmin,
                    util::div_floor(mu * x_info.best_below.numerator,
                                    x_info.best_below.denominator) +
                        floor_nmin_x_plus_y + 1u,
                    mu + nmin};
            }

            // Both L & R are not empty but x is not rational with denominator <= nlength.
            // Compute common constants.
            auto const& qstar = x_info.best_below.denominator;
            auto const& pstar = x_info.best_below.numerator;
            auto const& q = x_info.mod_inv;

            using gosper_coefficients = cntfrc::bilinear_fractional_transform<
                bigint::int_var, bigint::int_var, bigint::int_var, bigint::int_var, bigint::int_var,
                bigint::int_var, bigint::int_var, bigint::int_var>;

            auto k_common = [&](auto const& r) {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper<caching_xcf_ref, caching_ycf_ref>{
                        caching_xcf, caching_ycf,
                        gosper_coefficients{
                            // numerator
                            0, qstar * (nmin - q * r), util::to_signed(qstar),
                            (q * pstar + 1) * r - qstar * (floor_nmin_x_plus_y + 1u),
                            // denominator
                            0, util::to_signed(qstar * qstar), 0, qstar * -pstar}});
                cf.update();
                caching_xcf.rewind();
                caching_ycf.rewind();
                return -cf.current_partial_fraction().denominator;
            };

            auto r_common = [&](auto const& k) {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper<caching_xcf_ref, caching_ycf_ref>{
                        caching_xcf, caching_ycf,
                        gosper_coefficients{// numerator
                                            0, qstar * (qstar * k + nmin), util::to_signed(qstar),
                                            qstar * -(pstar * k + floor_nmin_x_plus_y + 1u),
                                            // denominator
                                            0, util::to_negative(q * qstar), 0, q * pstar + 1u}});
                cf.update();
                caching_xcf.rewind();
                caching_ycf.rewind();
                return -cf.current_partial_fraction().denominator;
            };

            auto const ceil_rT = qstar - floor_qstar_yprime;
            auto const ceil_rB = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::binary_gosper<caching_xcf_ref, caching_ycf_ref>{
                        caching_xcf, caching_ycf,
                        gosper_coefficients{// numerator
                                            0, qstar * nrange.upper_bound(), util::to_signed(qstar),
                                            0,
                                            // denominator
                                            0, 0, 0, 1}});
                cf.update();
                caching_xcf.rewind();
                caching_ycf.rewind();
                return qstar * (floor_nmin_x_plus_y + 1u) + pstar * nlength -
                       cf.current_partial_fraction().denominator;
            }();

            struct k_r_pair {
                bigint::int_var k;
                bigint::int_var r;
            };

            // Compute mu.
            auto k_r_pair_for_mu = [&]() -> k_r_pair {
                auto const floor_kPR = util::div_floor(nlength + q * (ceil_rB - 1u), qstar);
                auto const ceil_kTR = k_common(ceil_rT - 1u);
                auto const floor_kTL = util::div_floor(q * (ceil_rT - 1u), qstar);
                auto const ceil_rTR = r_common(ceil_kTR);
                auto const ceil_rTL = util::div_ceil(qstar * floor_kTL, q);
                auto const ceil_kTRR = k_common(ceil_rTR - 1u);
                auto const ceil_kTLR = k_common(ceil_rTL - 1u);

                if (ceil_rT - 1u < ceil_rB) {
                    return {floor_kPR, ceil_rB - 1u};
                }
                else if (ceil_kTR - 1u > floor_kTL) {
                    return {ceil_kTR - 1u, ceil_rT - 1u};
                }

                if (ceil_rTR - 1u >= ceil_rB) {
                    if (ceil_rTL >= ceil_rTR) {
                        return {ceil_kTLR - 1u, ceil_rTL - 1u};
                    }
                    else {
                        return {ceil_kTRR - 1u, ceil_rTR - 1u};
                    }
                }
                else {
                    if (ceil_rTL - 1u >= ceil_rB) {
                        return {ceil_kTLR - 1u, ceil_rTL - 1u};
                    }
                    else {
                        return {floor_kPR, ceil_rB - 1u};
                    }
                }
            }();
            util::constexpr_assert(k_r_pair_for_mu.r >= 0 && k_r_pair_for_mu.r < qstar);
            auto mu = qstar * std::move(k_r_pair_for_mu.k) - q * std::move(k_r_pair_for_mu.r);
            util::constexpr_assert(mu >= 1 && mu <= nlength);

            // Compute nu.
            auto k_r_pair_for_nu = [&]() -> k_r_pair {
                auto const rBprime = util::is_nonnegative(ceil_rB) ? ceil_rB : bigint::int_var{0};

                auto const floor_kPL = util::div_floor(q * ceil_rT, qstar);
                auto const ceil_kBL = k_common(rBprime);
                auto const floor_kBR = util::div_floor(nlength + q * rBprime, qstar);
                auto const ceil_rBL = r_common(floor_kBR);
                auto const ceil_rBR = util::div_ceil(nlength + qstar * ceil_kBL, q);
                auto const ceil_kBLL = k_common(ceil_rBL);
                auto const ceil_kBRL = k_common(ceil_rBR);

                if (rBprime >= ceil_rT) {
                    return {floor_kPL + 1u, ceil_rT};
                }
                else if (ceil_kBL <= floor_kBR) {
                    return {ceil_kBL, rBprime};
                }

                if (ceil_rBL < ceil_rT) {
                    if (ceil_rBR <= ceil_rBL) {
                        return {ceil_kBRL, ceil_rBR};
                    }
                    else {
                        return {ceil_kBLL, ceil_rBL};
                    }
                }
                else {
                    if (ceil_rBR < ceil_rT) {
                        return {ceil_kBRL, ceil_rBR};
                    }
                    else {
                        return {floor_kPL + 1u, ceil_rT};
                    }
                }
            }();
            util::constexpr_assert(k_r_pair_for_nu.r >= 0 && k_r_pair_for_nu.r < qstar);
            auto nu = qstar * std::move(k_r_pair_for_nu.k) - q * std::move(k_r_pair_for_nu.r);
            util::constexpr_assert(nu >= 1 && nu <= nlength);

            return xi_zeta_trapezoid{
                interval<frac_t, interval_type_t::bounded_left_closed_right_open>{
                    x_info.best_below, x_info.xi_upper_bound},
                util::div_floor(nu * pstar, qstar) + floor_nmin_x_plus_y + 1u, nu + nmin,
                util::div_floor(mu * pstar, qstar) + floor_nmin_x_plus_y + 1u, mu + nmin};
        }

        // Given real numbers x, y and a range [nmin:nmax] of integers, loop over all triples
        // (k,m,s) of integers such that
        // (1) k >= 0,
        // (2) floor(nx + y) = floor((nm + s)/2^k) holds for all n in [nmin:nmax],
        // (3) floor(nx) = floor(nm/2^k) holds for all n in [0:nmax-nmin], and
        // (4) k is the smallest among all integers satisfying the above three.
        // The parameter callback is called with each such triple.
        // callback should be able to take a parameter of type multiply_add_shift_info and return
        // bool. If it returns true, the loop continues, and if it returns false, the function
        // returns immediately.

        struct multiply_add_shift_info {
            bigint::int_var multiplier;
            bigint::int_var adder;
            std::size_t shift_amount;
        };

        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY,
                  class Callback>
        constexpr void for_each_simultaneous_multiply_add_shift(
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange,
            Callback&& callback) {
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorX>::template is_implementing_mixins<
                    cntfrc::convergent_tracker, cntfrc::interval_tracker>(),
                "the first continued fraction generator must implement convergent_tracker and "
                "interval_tracker");
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorY>::template is_implementing_mixins<
                    cntfrc::interval_tracker>(),
                "the second continued fraction generator must implement interval_tracker");

            using frac_t = frac<bigint::int_var, bigint::uint_var>;
            using nrange_t = interval<bigint::int_var, interval_type_t::bounded_closed>;
            auto const trapezoid = find_xi_zeta_region_simultaneous_floor(xcf, ycf, nrange);

            trapezoid.xi_range.visit([&callback, &trapezoid](auto const& itv) {
                using itv_type = std::remove_cvref_t<decltype(itv)>;

                if constexpr (itv_type::interval_type() == interval_type_t::entire) {
                    // xi can be any real number, and for each xi, there uniquely exists an integer
                    // zeta such that floor(nmin x + y) = floor(nmin xi + zeta) holds.
                    // Thus, we simply take k = 0, and loop over all integers xi, sorted with the
                    // absolute value.
                    bigint::uint_var abs_xi = 0u;
                    auto zeta_for_positive_xi = trapezoid.zeta_left_endpoint_constant_coeff;
                    auto zeta_for_negative_xi = trapezoid.zeta_left_endpoint_constant_coeff;

                    bool should_continue = callback(
                        multiply_add_shift_info{util::to_signed(abs_xi), zeta_for_positive_xi, 0});
                    while (should_continue) {
                        ++abs_xi;
                        zeta_for_positive_xi -= trapezoid.zeta_left_endpoint_negative_linear_coeff;
                        zeta_for_negative_xi += trapezoid.zeta_left_endpoint_negative_linear_coeff;

                        should_continue = callback(multiply_add_shift_info{
                            util::to_signed(abs_xi), zeta_for_positive_xi, 0});
                        if (!should_continue) {
                            break;
                        }
                        should_continue = callback(multiply_add_shift_info{
                            util::to_negative(abs_xi), zeta_for_negative_xi, 0});
                    }
                }
                else {
                    static_assert(itv_type::interval_type() ==
                                  interval_type_t::bounded_left_closed_right_open);

                    // k = ceil(log2(1/Delta)) since itv is not open.
                    auto k = [&] {
                        auto const delta = itv.upper_bound() - itv.lower_bound();
                        util::constexpr_assert(util::is_strictly_positive(delta.numerator));

                        return trunc_ceil_log2_div(delta.denominator, util::abs(delta.numerator));
                    }();
                    auto multiplier = util::div_ceil(itv.lower_bound().numerator << k,
                                                     itv.lower_bound().denominator);

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
                            auto next_lattice_point = multiplier + 1u;
                            if (next_lattice_point * itv.upper_bound().denominator <
                                (itv.upper_bound().numerator << k)) {
                                multiplier = std::move(next_lattice_point);
                                factor_out_power_of_2_limited();
                            }
                        }
                    }

                    auto right_end = util::div_ceil(itv.upper_bound().numerator << k,
                                                    itv.upper_bound().denominator);
                    bool succeeded = false;
                    while (true) {
                        for (; multiplier < right_end; ++multiplier) {
                            auto smin =
                                (trapezoid.zeta_left_endpoint_constant_coeff << k) -
                                multiplier * trapezoid.zeta_left_endpoint_negative_linear_coeff;
                            auto smax =
                                (trapezoid.zeta_right_endpoint_constant_coeff << k) -
                                multiplier * trapezoid.zeta_right_endpoint_negative_linear_coeff;

                            if (smin < smax) {
                                succeeded = true;
                                bool should_continue = true;
                                for (auto s = smin; s < smax; ++s) {
                                    should_continue =
                                        callback(multiply_add_shift_info{multiplier, s, k});
                                    if (!should_continue) {
                                        break;
                                    }
                                }
                                if (!should_continue) {
                                    break;
                                }
                            }
                        }

                        if (succeeded) {
                            break;
                        }
                        else {
                            --multiplier;
                            auto left_gap = frac_t{multiplier * itv.lower_bound().denominator -
                                                       (itv.lower_bound().numerator << k),
                                                   itv.lower_bound().denominator};
                            auto right_gap = frac_t{(itv.upper_bound().numerator << k) -
                                                        multiplier * itv.upper_bound().denominator,
                                                    itv.upper_bound().denominator};
                            util::constexpr_assert(util::is_nonnegative(left_gap.numerator) &&
                                                   util::is_strictly_positive(right_gap.numerator));

                            if (left_gap >= right_gap) {
                                k += trunc_ceil_log2_div(left_gap.denominator,
                                                         util::abs(left_gap.numerator));
                            }
                            else {
                                k += trunc_floor_log2_div(right_gap.denominator,
                                                          util::abs(right_gap.numerator)) +
                                     1u;
                            }

                            multiplier = util::div_ceil(itv.lower_bound().numerator << k,
                                                        itv.lower_bound().denominator);
                            right_end = util::div_ceil(itv.upper_bound().numerator << k,
                                                       itv.upper_bound().denominator);
                        }
                    } // while (true)
                }
            });
        }

        // Given real numbers x, y and a range [nmin:nmax] of integers, find a triple (k,m,s) of
        // integers such that
        // (1) k >= 0,
        // (2) floor(nx + y) = floor((nm + s)/2^k) holds for all n in [nmin:nmax],
        // (3) floor(nx) = floor(nm/2^k) holds for all n in [0:nmax-nmin], and
        // (4) k is the smallest among all integers satisfying the above three.

        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY>
        constexpr multiply_add_shift_info find_simultaneous_multiply_add_shift(
            ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
            interval<bigint::int_var, interval_type_t::bounded_closed> const& nrange) {
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorX>::template is_implementing_mixins<
                    cntfrc::convergent_tracker, cntfrc::interval_tracker>(),
                "the first continued fraction generator must implement convergent_tracker and "
                "interval_tracker");
            static_assert(
                std::remove_cvref_t<ContinuedFractionGeneratorY>::template is_implementing_mixins<
                    cntfrc::interval_tracker>(),
                "the second continued fraction generator must implement interval_tracker");

            multiply_add_shift_info result;
            for_each_simultaneous_multiply_add_shift(xcf, ycf, nrange, [&](auto&& info) {
                result = std::move(info);
                return false;
            });
            return result;
        }
    }
}

#endif
