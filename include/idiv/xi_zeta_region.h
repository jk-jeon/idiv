// Copyright 2024 Junekey Jeon
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

#ifndef JKJ_HEADER_XI_ZETA_REGION
#define JKJ_HEADER_XI_ZETA_REGION

#include "caching_generator.h"
#include "fractional_part_extremizer.h"
#include <algorithm>
#include <ranges>
#include <variant>

namespace jkj {
    namespace idiv {
        // Given real numbers x, y, and a range [nmin:nmax] of integers, find the set of (xi,zeta)
        // such that floor(nx + y) = floor(n xi + zeta) holds for all n in [nmin:nmax]. The set is
        // represented as the union of nonempty sets of one of the forms
        // {(xi,zeta) in J x R | a - b * xi < zeta < c - d * xi},
        // {(xi,zeta) in J x R | a - b * xi < zeta <= c - d * xi},
        // {(xi,zeta) in J x R | a - b * xi <= zeta < c - d * xi} and
        // {(xi,zeta) in J x R | a - b * xi <= zeta <= c - d * xi},
        // where J is a nonempty interval, a, b, c, d are rational numbers,
        // and the projection J onto the xi-axis of each of these sets is disjoint from each
        // other. We can in fact assume that J is either open and bounded, a singleton set, or the
        // entire R.
#if 0
        namespace xi_zeta_region {
            enum class boundary_type_t : bool { exclusive, inclusive };
            enum class region_type_t { entire_plane, infinite_parallelogram, bounded_polygon };
            using frac_t = frac<bigint::int_var, bigint::uint_var>;

            class entire_plane {
            public:
                static constexpr region_type_t region_type = region_type_t::entire_plane;
            };

            // Can be an infinite line when value_gap == 0.
            class infinite_parallelogram {
            public:
                static constexpr region_type_t region_type = region_type_t::infinite_paralellogram;

                bigint::int_var xi_coeff;
                bigint::uint_var zeta_coeff;

                bigint::int_var min_value;
                bigint::uint_var value_gap;

                boundary_type_t lower_boundary_type;
                boundary_type_t upper_boundary_type;
            };

            // Can be a line segment or a single point.
            class bounded_polygon {
            public:
                static constexpr region_type_t region_type = region_type_t::bounded_polygon;

                struct horizontally_open_trapezoid {
                    interval<frac_t const&, interval_type_t::bounded_open> xi_range;

                    frac_t const& lower_boundary_constant_coeff;
                    frac_t const& lower_boundary_negative_linear_coeff;
                    frac_t const& upper_boundary_constant_coeff;
                    frac_t const& upper_boundary_negative_linear_coeff;

                    boundary_type_t const lower_boundary_type;
                    boundary_type_t const upper_boundary_type;
                };

                struct vertical_line_segment {
                    frac_t const& xi;
                    variable_shape_interval<frac_t const&, interval_type_t::bounded_open,
                                            interval_type_t::bounded_left_open_right_closed,
                                            interval_type_t::left_closed_right_open,
                                            interval_type_t::bounded_closed> const zeta_range;
                };

                struct single_point {
                    frac_t const& xi;
                    frac_t const& zeta;
                };

                template <class Functor>
                void for_each_horizontal_slice(Functor&& f) const {
                    util::constexpr_assert(!turning_points.empty());
                    util::constexpr_assert(!boundary_lines_.empty());

                    auto on_turning_point = [](auto const& turning_point) {
                        if (turning_point.lower_boundary_included) {
                            if (turning_point.upper_boundary_included) {
                                f(vertical_line_segment{
                                    turning_point.xi,
                                    interval<frac_t, interval_type_t::bounded_closed>{
                                        turning_point.lower_zeta, turning_point.upper_zeta}});
                            }
                            else {
                                f(vertical_line_segment{
                                    turning_point.xi,
                                    interval<frac_t,
                                             interval_type_t::bounded_left_closed_right_open>{
                                        turning_point.lower_zeta, turning_point.upper_zeta}});
                            }
                        }
                        else {
                            if (turning_point.upper_boundary_included) {
                                f(vertical_line_segment{
                                    turning_point.xi,
                                    interval<frac_t,
                                             interval_type_t::bounded_left_open_right_closed>{
                                        turning_point.lower_zeta, turning_point.upper_zeta}});
                            }
                            else {
                                f(vertical_line_segment{
                                    turning_point.xi,
                                    interval<frac_t, interval_type_t::bounded_open>{
                                        turning_point.lower_zeta, turning_point.upper_zeta}});
                            }
                        }
                    };

                    auto on_last_turning_point = [](auto const& turning_point) {
                        if (turning_point.lower_zeta == turning_point.upper_zeta) {
                            util::constexpr_assert(turning_point.lower_boundary_included &&
                                                   turning_point.upper_boundary_included);
                            f(single_point{turning_point.xi, turning_point.lower_zeta});
                        }
                        else {
                            on_turning_point(turning_point);
                        }
                    };

                    if (left_boundary_included()) {
                        on_last_turning_point(turning_points().front());
                    }

                    for (std::size_t turning_point_idx = 1;
                         turning_point_idx < turning_points.size(); ++turning_point_idx) {
                        f(horizontally_open_trapezoid{
                            interval<frac_t, interval_type_t::bounded_open>{
                                turning_points()[turning_point_idx - 1].xi,
                                turning_points()[turning_point_idx].xi},
                            boundary_line_pairs()[turning_point_idx - 1]
                                .lower_boundary_constant_coeff,
                            boundary_line_pairs()[turning_point_idx - 1]
                                .lower_boundary_negative_linear_coeff,
                            boundary_line_pairs()[turning_point_idx - 1]
                                .upper_boundary_constant_coeff,
                            boundary_line_pairs()[turning_point_idx - 1]
                                .upper_boundary_negative_linear_coeff,
                            boundary_line_pairs()[turning_point_idx - 1].lower_boundary_included,
                            boundary_line_pairs()[turning_point_idx - 1].upper_boundary_included});

                        if (turning_point_idx < turning_points.size() - 1) {
                            on_turning_point(turning_points()[turning_point_idx]);
                        }
                        else if (include_right_boundary()) {
                            on_last_turning_point(turning_points().back());
                        }
                    }
                }

                struct turning_point_info {
                    frac_t xi;
                    frac_t lower_zeta;
                    frac_t upper_zeta;

                    boundary_type_t lower_boundary_type;
                    boundary_type_t upper_boundary_type;
                };

                struct boundary_line_pair_info {
                    frac_t lower_boundary_constant_coeff;
                    frac_t lower_boundary_negative_linear_coeff;
                    frac_t upper_boundary_constant_coeff;
                    frac_t upper_boundary_negative_linear_coeff;

                    boundary_type_t lower_boundary_type;
                    boundary_type_t upper_boundary_type;
                };

                auto const& turning_points() const noexcept { return turning_points_; }
                auto const& boundary_line_pairs() const noexcept { return boundary_line_pairs_; }
                boundary_type_t left_boundary_type() const noexcept { return left_boundary_type_; }
                boundary_type_t right_boundary_type() const noexcept { return right_boundary_type_; }


                template <std::ranges::common_range TurningPoints,
                          std::ranges::common_range BoundaryLinePairs>
                bounded_polygon(TurningPoints&& turning_point_rg,
                                BoundaryLinePairs&& boundary_line_pair_rg,
                    boundary_type_t left_boundary_type_in, boundary_type_t right_boundary_type_in)
                    : turning_points_(turning_point_rg.cbegin(), turning_point_rg.cend()),
                      boundary_line_pairs_(boundary_line_pair_rg.cbegin(),
                                           boundary_line_pair_rg.cend()),
                    left_boundary_type_{ left_boundary_type_in },
                    right_boundary_type_{ right_boundary_type_in } {}

            private:
                std::vector<turning_point_info> turning_points_;
                std::vector<boundary_line_pair_info> boundary_line_pairs_;
                boundary_type_t left_boundary_type_;
                boundary_type_t right_boundary_type_;
            };

            using variable_shape_region =
                std::variant<entire_plane, infinite_parallelogram, bounded_polygon>;
        }

        struct floor_constraint_spec {
            using frac_t = frac<bigint::int_var, bigint::uint_var>;

            // (x, y) |-> (linear_coeff.xx * x + linear_coeff.xy * y + constant_coeff_x,
            //             linear_coeff.yx * x + linear_coeff.yy * y + constant_coeff_y).
            struct affine_transform {
                struct linear_transform {
                    frac_t xx; // x to x.
                    frac_t xy; // y to x.
                    frac_t yx; // x to y.
                    frac_t yy; // y to y.
                } linear_coeff;

                frac_t constant_coeff_x;
                frac_t constant_coeff_y;
            } affine_coeff;

            interval<bigint::int_var, interval_type_t::bounded_closed> nrange;
        };

        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY,
                  std::ranges::range RangeOfConstraintSpec>
        constexpr xi_zeta_region::variable_shape_region
        find_xi_zeta_region(ContinuedFractionGeneratorX&& xcf, ContinuedFractionGeneratorY&& ycf,
                            RangeOfConstraintSpec&& range_of_constraint_specs) {
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
                frac_t xi_coeff;
                frac_t zeta_coeff;
                frac_t eta_coeff;
                enum class boundary_type_t : bool { inclusive, exclusive } boundary_type;
            };
            std::vector<half_space_info> right_half_spaces; // Lower bounds for xi.
            std::vector<half_space_info> left_half_spaces;  // Upper bounds for xi.

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

                            auto movement =
                                find_extremizers_of_fractional_part(approx_xcf, max_diff)
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
                            auto movement =
                                find_extremizers_of_fractional_part(approx_xcf, max_diff)
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
                                {frac_t{-1, base_point},
                                 frac_t{1 + ((base_point * approx_x_y_info.multiplier +
                                              approx_x_y_info.adder) >>
                                             approx_x_y_info.shift_amount),
                                        base_point},
                                 half_space_info::boundary_type_t::exclusive});

                            auto movement =
                                find_extremizers_of_fractional_part(approx_xcf, max_diff)
                                    .smallest_minimizer;
                            approx_xcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point -= movement;
                            max_diff -= std::move(movement);
                        }

                        half_spaces.push_back(
                            {frac_t{-1, base_point},
                             frac_t{1 + ((base_point * approx_x_y_info.multiplier +
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
                            auto movement =
                                find_extremizers_of_fractional_part(approx_xcf, max_diff)
                                    .largest_maximizer;
                            approx_xcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point += movement;
                            max_diff -= std::move(movement);

                            half_spaces.push_back(
                                {frac_t{-1, base_point},
                                 frac_t{1 + ((base_point * approx_x_y_info.multiplier +
                                              approx_x_y_info.adder) >>
                                             approx_x_y_info.shift_amount),
                                        base_point},
                                 half_space_info::boundary_type_t::exclusive});
                        }
                    };

                    auto process_single_sign_interval = [&](nrange_t const& nrange,
                                                            elementary_problem_sign sign,
                                                            auto&& pm_xcf, auto&& pm_ycf) {
                        auto base_points =
                            find_extremizers_of_fractional_part(pm_xcf, pm_ycf, nrange);
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
                        // Zero.
                        if (util::is_zero(nrange.lower_bound()) &&
                            util::is_zero(nrange.upper_bound())) {
                            nrange_contains_zero = true;
                        }
                        // Negative interval.
                        else if (util::is_strictly_negative(nrange.upper_bound())) {
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

            // If n = 0 is the only constraint, then return early.
            if (right_half_spaces.empty()) {
                util::constexpr_assert(nrange_contains_zero && left_half_spaces.empty());
                // zeta should satisfy the inequality floor_y <= zeta < floor_y + 1.
                auto floor_y_frac = frac_t{floor_y, 1u};
                auto floor_y_p1_frac = frac_t{floor_y + 1, 1u};
                return {
                    elementary_xi_zeta_region{interval<frac_t, interval_type_t::bounded_closed>{
                                                  floor_y_frac, floor_y_frac},
                                              floor_y_frac.numerator, 0u, floor_y_p1_frac.numerator,
                                              0u, true, false},
                    elementary_xi_zeta_region{interval<frac_t, interval_type_t::bounded_open>{
                                                  floor_y_frac, floor_y_p1_frac},
                                              floor_y_frac.numerator, 0u, floor_y_p1_frac.numerator,
                                              0u, true, false},
                };
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
            enum class bounding_direction_t : bool { lower, upper };
            auto compute_one_sided_intersection = [&](bounding_direction_t bounding_direction) {
                struct vec2d {
                    frac_t zeta_coord;
                    frac_t eta_coord;

                    constexpr frac_t dot(vec2d const& other) const {
                        return zeta_coord * other.zeta_coord + eta_coord * other.eta_coord;
                    }
                    constexpr frac_t normsq() const { return dot(*this); }
                };

                auto const& half_spaces = bounding_direction == bounding_direction_t::lower
                                              ? right_half_spaces
                                              : left_half_spaces;
                util::constexpr_assert(!half_spaces.empty());

                auto invert_sign_wrt_bounding_direction =
                    [bounding_direction](bigint::int_var const& n) {
                        if (bounding_direction == bounding_direction_t::lower) {
                            return n;
                        }
                        else {
                            return -n;
                        }
                    };

                std::vector<elementary_one_sided_region> result;

                // Start from the one with the largest zeta-coordinate, which corresponds to the
                // half-space with the highst boundary line.
                auto first_elmt = std::ranges::max_element(std::as_const(half_spaces), {},
                                                           &half_space_info::zeta_coeff);

                // If there is only one half-space, return immediately.
                if (half_spaces.size() == 1) {
                    result.push_back(
                        {interval<frac_t, interval_type_t::entire>{},
                         util::is_nonnegative(first_elmt->zeta_coeff.numerator)
                             ? -first_elmt->eta_coeff.numerator
                             : first_elmt->eta_coeff.numerator,
                         invert_sign_wrt_bounding_direction(
                             util::is_nonnegative(first_elmt->zeta_coeff.numerator)
                                 ? util::to_signed(first_elmt->zeta_coeff.denominator)
                                 : -util::to_signed(first_elmt->zeta_coeff.denominator)),
                         first_elmt->boundary_type == half_space_info::boundary_type_t::inclusive});
                    return result;
                }

                // We are at the right-end and we want to travel clockwise, when viewed from the
                // positive xi-axis. To do so, we set the initial direction to be along the
                // negative eta-axis.
                auto prev_direction_vec = vec2d{frac_t{0, 1u}, frac_t{-1, 1u}};
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
                    // Find the point whose direction vector is the closest in angle to the
                    // previous direction vector.
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
                    // For the case of lower bound, every functional should have nonnegative
                    // inner product with this normal, while for the case of upper bound, they
                    // should have nonpositive inner product.
                    auto face_normal_zeta =
                        current_angle_info.itr->eta_coeff - last_elmt->eta_coeff;
                    auto face_normal_eta =
                        last_elmt->zeta_coeff - current_angle_info.itr->zeta_coeff;

                    // Project it down to the plane eta = 1. The resulting point must be a
                    // turning point.
                    auto turning_point_zeta =
                        get_reduced_quotient(face_normal_zeta, face_normal_eta);

                    auto xi_endpoint_numerator =
                        util::is_nonnegative(last_elmt->zeta_coeff.numerator)
                            ? -last_elmt->eta_coeff.numerator
                            : last_elmt->eta_coeff.numerator;
                    auto xi_endpoint_denominator = invert_sign_wrt_bounding_direction(
                        util::is_nonnegative(last_elmt->zeta_coeff.numerator)
                            ? util::to_signed(last_elmt->zeta_coeff.denominator)
                            : -util::to_signed(last_elmt->zeta_coeff.denominator));

                    // The unbounded bottom boundary line.
                    if (last_elmt == first_elmt) {
                        util::constexpr_assert(
                            util::is_strictly_positive(face_normal_eta.numerator));

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
                            util::is_strictly_negative(face_normal_eta.numerator));

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
                            util::is_strictly_positive(face_normal_eta.numerator));

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
                auto lower_bound_region =
                    compute_one_sided_intersection(bounding_direction_t::lower);
                auto upper_bound_region =
                    compute_one_sided_intersection(bounding_direction_t::upper);

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
                                        [&](auto const&) { ++upper_bound_itr; },
                                        [&] {
                                            // Impossible to reach here.
                                        });
                                });
                        }
                    } // while (true)
                }     // End of the branching on lower_bound_region.size().
            }

            // Remove top/bottom vertices if they do not belong to the region.
            util::constexpr_assert(result.front().zeta_range.interval_type() ==
                                   interval_type_t::bounded_closed);
            util::constexpr_assert(result.back().zeta_range.interval_type() ==
                                   interval_type_t::bounded_closed);
            if (!result.front().xi_left_endpoint_included ||
                !result.front().xi_right_endpoint_included) {
                result.erase(result.begin());
            }
            if (!result.back().xi_left_endpoint_included ||
                !result.back().xi_right_endpoint_included) {
                result.pop_back();
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

                // Find floor(y).
                bool found_both_at_the_same_time = false;
                for (; src_itr != result.end(); ++src_itr) {
                    if (src_itr->zeta_range.contains(floor_y_frac)) {
                        found_both_at_the_same_time = src_itr->zeta_range.contains(floor_y_p1_frac);

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
                                        floor_y_frac,
                                        found_both_at_the_same_time ? floor_y_p1_frac : ub};
                                },
                                [&] {
                                    // Cannot reach here.
                                    return interval<frac_t, interval_type_t::bounded_open>{
                                        floor_y_frac, floor_y_p1_frac};
                                });
                        }
                        trimmed.push_back(std::move(*src_itr));
                        break;
                    }
                } // End of the for loop.

                // If both are found in the same region, then we are done.
                if (!found_both_at_the_same_time) {
                    // If floor(y) wasn't found, go back to the beginning.
                    if (src_itr == result.end()) {
                        src_itr = result.begin();
                    }
                    // If floor(y) was found, move to the next region.
                    else {
                        ++src_itr;
                    }

                    // Find floor(y) + 1.
                    for (; src_itr != result.end(); ++src_itr) {
                        if (src_itr->zeta_range.contains(floor_y_p1_frac)) {
                            if (src_itr->zeta_range.interval_type() ==
                                interval_type_t::bounded_open) {
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
                            break;
                        }
                        trimmed.push_back(std::move(*src_itr));
                    }
                } // if (!found_both_at_the_same_time)

                result = std::move(trimmed);
            } // if (nrange_contains_zero)

            return result;
        }
#else
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

                            auto movement =
                                find_extremizers_of_fractional_part(approx_xcf, max_diff)
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
                            auto movement =
                                find_extremizers_of_fractional_part(approx_xcf, max_diff)
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
                                {frac_t{-1, base_point},
                                 frac_t{1 + ((base_point * approx_x_y_info.multiplier +
                                              approx_x_y_info.adder) >>
                                             approx_x_y_info.shift_amount),
                                        base_point},
                                 half_space_info::boundary_type_t::exclusive});

                            auto movement =
                                find_extremizers_of_fractional_part(approx_xcf, max_diff)
                                    .smallest_minimizer;
                            approx_xcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point -= movement;
                            max_diff -= std::move(movement);
                        }

                        half_spaces.push_back(
                            {frac_t{-1, base_point},
                             frac_t{1 + ((base_point * approx_x_y_info.multiplier +
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
                            auto movement =
                                find_extremizers_of_fractional_part(approx_xcf, max_diff)
                                    .largest_maximizer;
                            approx_xcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point += movement;
                            max_diff -= std::move(movement);

                            half_spaces.push_back(
                                {frac_t{-1, base_point},
                                 frac_t{1 + ((base_point * approx_x_y_info.multiplier +
                                              approx_x_y_info.adder) >>
                                             approx_x_y_info.shift_amount),
                                        base_point},
                                 half_space_info::boundary_type_t::exclusive});
                        }
                    };

                    auto process_single_sign_interval = [&](nrange_t const& nrange,
                                                            elementary_problem_sign sign,
                                                            auto&& pm_xcf, auto&& pm_ycf) {
                        auto base_points =
                            find_extremizers_of_fractional_part(pm_xcf, pm_ycf, nrange);
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
                        // Zero.
                        if (util::is_zero(nrange.lower_bound()) &&
                            util::is_zero(nrange.upper_bound())) {
                            nrange_contains_zero = true;
                        }
                        // Negative interval.
                        else if (util::is_strictly_negative(nrange.upper_bound())) {
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

            // If n = 0 is the only constraint, then return early.
            if (right_half_spaces.empty()) {
                util::constexpr_assert(nrange_contains_zero && left_half_spaces.empty());
                // zeta should satisfy the inequality floor_y <= zeta < floor_y + 1.
                auto floor_y_frac = frac_t{floor_y, 1u};
                auto floor_y_p1_frac = frac_t{floor_y + 1, 1u};
                return {
                    elementary_xi_zeta_region{interval<frac_t, interval_type_t::bounded_closed>{
                                                  floor_y_frac, floor_y_frac},
                                              floor_y_frac.numerator, 0u, floor_y_p1_frac.numerator,
                                              0u, true, false},
                    elementary_xi_zeta_region{interval<frac_t, interval_type_t::bounded_open>{
                                                  floor_y_frac, floor_y_p1_frac},
                                              floor_y_frac.numerator, 0u, floor_y_p1_frac.numerator,
                                              0u, true, false},
                };
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
            enum class bounding_direction_t : bool { lower, upper };
            auto compute_one_sided_intersection = [&](bounding_direction_t bounding_direction) {
                struct vec2d {
                    frac_t zeta_coord;
                    frac_t eta_coord;

                    constexpr frac_t dot(vec2d const& other) const {
                        return zeta_coord * other.zeta_coord + eta_coord * other.eta_coord;
                    }
                    constexpr frac_t normsq() const { return dot(*this); }
                };

                auto const& half_spaces = bounding_direction == bounding_direction_t::lower
                                              ? right_half_spaces
                                              : left_half_spaces;
                util::constexpr_assert(!half_spaces.empty());

                auto invert_sign_wrt_bounding_direction =
                    [bounding_direction](bigint::int_var const& n) {
                        if (bounding_direction == bounding_direction_t::lower) {
                            return n;
                        }
                        else {
                            return -n;
                        }
                    };

                std::vector<elementary_one_sided_region> result;

                // Start from the one with the largest zeta-coordinate, which corresponds to the
                // half-space with the highst boundary line.
                auto first_elmt = std::ranges::max_element(std::as_const(half_spaces), {},
                                                           &half_space_info::zeta_coeff);

                // If there is only one half-space, return immediately.
                if (half_spaces.size() == 1) {
                    result.push_back(
                        {interval<frac_t, interval_type_t::entire>{},
                         util::is_nonnegative(first_elmt->zeta_coeff.numerator)
                             ? -first_elmt->eta_coeff.numerator
                             : first_elmt->eta_coeff.numerator,
                         invert_sign_wrt_bounding_direction(
                             util::is_nonnegative(first_elmt->zeta_coeff.numerator)
                                 ? util::to_signed(first_elmt->zeta_coeff.denominator)
                                 : -util::to_signed(first_elmt->zeta_coeff.denominator)),
                         first_elmt->boundary_type == half_space_info::boundary_type_t::inclusive});
                    return result;
                }

                // We are at the right-end and we want to travel clockwise, when viewed from the
                // positive xi-axis. To do so, we set the initial direction to be along the
                // negative eta-axis.
                auto prev_direction_vec = vec2d{frac_t{0, 1u}, frac_t{-1, 1u}};
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
                    // Find the point whose direction vector is the closest in angle to the
                    // previous direction vector.
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
                    // For the case of lower bound, every functional should have nonnegative
                    // inner product with this normal, while for the case of upper bound, they
                    // should have nonpositive inner product.
                    auto face_normal_zeta =
                        current_angle_info.itr->eta_coeff - last_elmt->eta_coeff;
                    auto face_normal_eta =
                        last_elmt->zeta_coeff - current_angle_info.itr->zeta_coeff;

                    // Project it down to the plane eta = 1. The resulting point must be a
                    // turning point.
                    auto turning_point_zeta =
                        get_reduced_quotient(face_normal_zeta, face_normal_eta);

                    auto xi_endpoint_numerator =
                        util::is_nonnegative(last_elmt->zeta_coeff.numerator)
                            ? -last_elmt->eta_coeff.numerator
                            : last_elmt->eta_coeff.numerator;
                    auto xi_endpoint_denominator = invert_sign_wrt_bounding_direction(
                        util::is_nonnegative(last_elmt->zeta_coeff.numerator)
                            ? util::to_signed(last_elmt->zeta_coeff.denominator)
                            : -util::to_signed(last_elmt->zeta_coeff.denominator));

                    // The unbounded bottom boundary line.
                    if (last_elmt == first_elmt) {
                        util::constexpr_assert(
                            util::is_strictly_positive(face_normal_eta.numerator));

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
                            util::is_strictly_negative(face_normal_eta.numerator));

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
                            util::is_strictly_positive(face_normal_eta.numerator));

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
                auto lower_bound_region =
                    compute_one_sided_intersection(bounding_direction_t::lower);
                auto upper_bound_region =
                    compute_one_sided_intersection(bounding_direction_t::upper);

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
                                        [&](auto const&) { ++upper_bound_itr; },
                                        [&] {
                                            // Impossible to reach here.
                                        });
                                });
                        }
                    } // while (true)
                }     // End of the branching on lower_bound_region.size().
            }

            // Remove top/bottom vertices if they do not belong to the region.
            util::constexpr_assert(result.front().zeta_range.interval_type() ==
                                   interval_type_t::bounded_closed);
            util::constexpr_assert(result.back().zeta_range.interval_type() ==
                                   interval_type_t::bounded_closed);
            if (!result.front().xi_left_endpoint_included ||
                !result.front().xi_right_endpoint_included) {
                result.erase(result.begin());
            }
            if (!result.back().xi_left_endpoint_included ||
                !result.back().xi_right_endpoint_included) {
                result.pop_back();
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

                // Find floor(y).
                bool found_both_at_the_same_time = false;
                for (; src_itr != result.end(); ++src_itr) {
                    if (src_itr->zeta_range.contains(floor_y_frac)) {
                        found_both_at_the_same_time = src_itr->zeta_range.contains(floor_y_p1_frac);

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
                                        floor_y_frac,
                                        found_both_at_the_same_time ? floor_y_p1_frac : ub};
                                },
                                [&] {
                                    // Cannot reach here.
                                    return interval<frac_t, interval_type_t::bounded_open>{
                                        floor_y_frac, floor_y_p1_frac};
                                });
                        }
                        trimmed.push_back(std::move(*src_itr));
                        break;
                    }
                } // End of the for loop.

                // If both are found in the same region, then we are done.
                if (!found_both_at_the_same_time) {
                    // If floor(y) wasn't found, go back to the beginning.
                    if (src_itr == result.end()) {
                        src_itr = result.begin();
                    }
                    // If floor(y) was found, move to the next region.
                    else {
                        ++src_itr;
                    }

                    // Find floor(y) + 1.
                    for (; src_itr != result.end(); ++src_itr) {
                        if (src_itr->zeta_range.contains(floor_y_p1_frac)) {
                            if (src_itr->zeta_range.interval_type() ==
                                interval_type_t::bounded_open) {
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
                            break;
                        }
                        trimmed.push_back(std::move(*src_itr));
                    }
                } // if (!found_both_at_the_same_time)

                result = std::move(trimmed);
            } // if (nrange_contains_zero)

            return result;
        }
#endif
    }
}

#endif
