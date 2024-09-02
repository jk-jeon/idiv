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

#ifndef JKJ_HEADER_IDIV_SIMULTANEOUS_FLOOR
#define JKJ_HEADER_IDIV_SIMULTANEOUS_FLOOR

#include "best_rational_approx.h"
#include "continued_fraction/engine/caching.h"
#include "continued_fraction/engine/gosper.h"
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

            // TODO: deal with possible rational dependence between x and y.

            using frac_t = frac<bigint::int_var, bigint::uint_var>;
            using binary_gosper_t = cntfrc::engine::binary_gosper<ContinuedFractionGeneratorX&,
                                                                  ContinuedFractionGeneratorY&>;

            auto const& nmin = nrange.lower_bound();
            auto const nlength = util::abs(nrange.upper_bound() - nrange.lower_bound());

            // floor(nmin * x + y).
            auto floor_nmin_x_plus_y = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    binary_gosper_t{xcf,
                                    ycf,
                                    {// numerator
                                     0, nmin, 1, 0,
                                     // denominator
                                     0, 0, 0, 1}});
                return cf.current_partial_fraction().denominator;
            }();
            xcf.rewind();
            ycf.rewind();

            // Early return for the case nmin == nmax.
            if (util::is_zero(nlength)) {
                return xi_zeta_trapezoid{interval<frac_t, interval_type_t::entire>{},
                                         floor_nmin_x_plus_y, nmin, floor_nmin_x_plus_y + 1u, nmin};
            }

            auto xi_range = find_floor_quotient_range(xcf, nlength);
            bool xcf_terminated = xcf.terminated();
            xcf.rewind();

            // See if L is empty.
            {
                auto const lhs = [&] {
                    auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                        binary_gosper_t{xcf,
                                        ycf,
                                        {// numerator
                                         0, xi_range.lower_bound().denominator + nmin, 1, 0,
                                         // denominator
                                         0, 0, 0, 1}});
                    return cf.current_partial_fraction().denominator;
                }();
                xcf.rewind();
                ycf.rewind();
                auto rhs = xi_range.lower_bound().numerator + floor_nmin_x_plus_y;

                if (lhs > rhs) {
                    // L is empty.
                    return xi_zeta_trapezoid{xi_range, rhs + 1u,
                                             xi_range.lower_bound().denominator + nmin,
                                             floor_nmin_x_plus_y + 1u, nmin};
                }
            }

            auto const floor_qstar_yprime = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    binary_gosper_t{xcf,
                                    ycf,
                                    {// numerator
                                     0, xi_range.lower_bound().denominator * nmin,
                                     util::to_signed(xi_range.lower_bound().denominator), 0,
                                     // denominator
                                     0, 0, 0, 1}});
                return cf.current_partial_fraction().denominator -
                       xi_range.lower_bound().denominator * floor_nmin_x_plus_y;
            }();
            xcf.rewind();
            ycf.rewind();

            // See if R is empty.
            {
                if (xcf_terminated) {
                    if (util::is_zero(floor_qstar_yprime)) {
                        // R is empty.
                        return xi_zeta_trapezoid{
                            xi_range, floor_nmin_x_plus_y, nmin,
                            util::div_floor(xi_range.upper_bound().denominator *
                                                xi_range.lower_bound().numerator,
                                            xi_range.lower_bound().denominator) +
                                floor_nmin_x_plus_y + 1u,
                            xi_range.upper_bound().denominator + nmin};
                    }
                }
                else {
                    auto const lhs = [&] {
                        auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                            binary_gosper_t{xcf,
                                            ycf,
                                            {// numerator
                                             0, xi_range.upper_bound().denominator + nmin, 1, 0,
                                             // denominator
                                             0, 0, 0, 1}});
                        return cf.current_partial_fraction().denominator;
                    }();
                    xcf.rewind();
                    ycf.rewind();
                    auto rhs = xi_range.upper_bound().numerator + floor_nmin_x_plus_y;

                    if (lhs < rhs) {
                        // R is empty.
                        return xi_zeta_trapezoid{
                            xi_range, floor_nmin_x_plus_y, nmin,
                            util::div_floor(xi_range.upper_bound().denominator *
                                                xi_range.lower_bound().numerator,
                                            xi_range.lower_bound().denominator) +
                                floor_nmin_x_plus_y + 1u,
                            xi_range.upper_bound().denominator + nmin};
                    }
                }
            }

            // Both L & R are not empty but x is rational with denominator <= nlength.
            if (xcf_terminated) {
                auto mu = ((floor_qstar_yprime + 1u) * xi_range.upper_bound().denominator) %
                          xi_range.lower_bound().denominator;
                auto nu = (floor_qstar_yprime * xi_range.upper_bound().denominator) %
                          xi_range.lower_bound().denominator;

                return xi_zeta_trapezoid{xi_range,
                                         util::div_floor(nu * xi_range.lower_bound().numerator,
                                                         xi_range.lower_bound().denominator) +
                                             floor_nmin_x_plus_y + 1u,
                                         nu + nmin,
                                         util::div_floor(mu * xi_range.lower_bound().numerator,
                                                         xi_range.lower_bound().denominator) +
                                             floor_nmin_x_plus_y + 1u,
                                         mu + nmin};
            }

            // Both L & R are not empty but x is not rational with denominator <= nlength.
            // Compute common constants.
            auto const& qlstar = xi_range.lower_bound().denominator;
            auto const& plstar = xi_range.lower_bound().numerator;
            auto const& qustar = xi_range.upper_bound().denominator;

            using gosper_coefficients = cntfrc::bilinear_fractional_mapping<bigint::int_var>;

            auto k_common = [&](auto const& r) {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(binary_gosper_t{
                    xcf, ycf,
                    gosper_coefficients{
                        // numerator
                        0, qlstar * (nmin - qustar * r), util::to_signed(qlstar),
                        (qustar * plstar + 1) * r - qlstar * (floor_nmin_x_plus_y + 1u),
                        // denominator
                        0, util::to_signed(qlstar * qlstar), 0, qlstar * -plstar}});
                auto ret_value = -cf.current_partial_fraction().denominator;
                xcf.rewind();
                ycf.rewind();
                return ret_value;
            };

            auto r_common = [&](auto const& k) {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(binary_gosper_t{
                    xcf, ycf,
                    gosper_coefficients{// numerator
                                        0, qlstar * (qlstar * k + nmin), util::to_signed(qlstar),
                                        qlstar * -(plstar * k + floor_nmin_x_plus_y + 1u),
                                        // denominator
                                        0, util::to_negative(qustar * qlstar), 0,
                                        qustar * plstar + 1u}});
                auto ret_value = -cf.current_partial_fraction().denominator;
                xcf.rewind();
                ycf.rewind();
                return ret_value;
            };

            auto const ceil_rT = qlstar - floor_qstar_yprime;
            auto const ceil_rB = [&] {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    binary_gosper_t{xcf, ycf,
                                    gosper_coefficients{// numerator
                                                        0, qlstar * nrange.upper_bound(),
                                                        util::to_signed(qlstar), 0,
                                                        // denominator
                                                        0, 0, 0, 1}});
                return qlstar * (floor_nmin_x_plus_y + 1u) + plstar * nlength -
                       cf.current_partial_fraction().denominator;
            }();
            xcf.rewind();
            ycf.rewind();

            struct k_r_pair {
                bigint::int_var k;
                bigint::int_var r;
            };

            // Compute mu.
            auto k_r_pair_for_mu = [&]() -> k_r_pair {
                auto const floor_kPR = util::div_floor(nlength + qustar * (ceil_rB - 1u), qlstar);
                auto const ceil_kTR = k_common(ceil_rT - 1u);
                auto const floor_kTL = util::div_floor(qustar * (ceil_rT - 1u), qlstar);
                auto const ceil_rTR = r_common(ceil_kTR);
                auto const ceil_rTL = util::div_ceil(qlstar * floor_kTL, qustar);
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
            util::constexpr_assert(k_r_pair_for_mu.r >= 0 && k_r_pair_for_mu.r < qlstar);
            auto mu = qlstar * std::move(k_r_pair_for_mu.k) - qustar * std::move(k_r_pair_for_mu.r);
            util::constexpr_assert(mu >= 1 && mu <= nlength);

            // Compute nu.
            auto k_r_pair_for_nu = [&]() -> k_r_pair {
                auto const rBprime = util::is_nonnegative(ceil_rB) ? ceil_rB : bigint::int_var{0};

                auto const floor_kPL = util::div_floor(qustar * ceil_rT, qlstar);
                auto const ceil_kBL = k_common(rBprime);
                auto const floor_kBR = util::div_floor(nlength + qustar * rBprime, qlstar);
                auto const ceil_rBL = r_common(floor_kBR);
                auto const ceil_rBR = util::div_ceil(qlstar * ceil_kBL - nlength, qustar);
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
            util::constexpr_assert(k_r_pair_for_nu.r >= 0 && k_r_pair_for_nu.r < qlstar);
            auto nu = qlstar * std::move(k_r_pair_for_nu.k) - qustar * std::move(k_r_pair_for_nu.r);
            util::constexpr_assert(nu >= 1 && nu <= nlength);

            return xi_zeta_trapezoid{
                xi_range, util::div_floor(nu * plstar, qlstar) + floor_nmin_x_plus_y + 1u,
                nu + nmin, util::div_floor(mu * plstar, qlstar) + floor_nmin_x_plus_y + 1u,
                mu + nmin};
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

            using frac_t = frac<bigint::int_var, bigint::uint_var>;
            auto const trapezoid = find_xi_zeta_region_simultaneous_floor(
                std::forward<ContinuedFractionGeneratorX>(xcf),
                std::forward<ContinuedFractionGeneratorY>(ycf), nrange);

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

                        return util::trunc_ceil_log2_div(delta.denominator,
                                                         util::abs(delta.numerator));
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
                                k += util::trunc_ceil_log2_div(left_gap.denominator,
                                                               util::abs(left_gap.numerator));
                            }
                            else {
                                k += util::trunc_floor_log2_div(right_gap.denominator,
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

            multiply_add_shift_info result;
            for_each_simultaneous_multiply_add_shift(std::forward<ContinuedFractionGeneratorX>(xcf),
                                                     std::forward<ContinuedFractionGeneratorY>(ycf),
                                                     nrange, [&](auto&& info) {
                                                         result = std::move(info);
                                                         return false;
                                                     });
            return result;
        }
    }
}

#endif
