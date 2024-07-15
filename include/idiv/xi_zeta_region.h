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
        // Given real numbers x, y, and a set of "constraints", which are pairs of a rational 2D
        // affine transform and a range [nmin:nmax] of integers, find the set of (xi,zeta) such that
        // for each constraint (T,[nmin:nmax]), floor(nx' + y') = floor(n xi' + zeta') holds for all
        // n in [nmin:nmax] where (x',y') = T(x,y) and (xi',zeta') = T(xi,zeta). The returned set is
        // one of the following kinds:
        //
        // 1. Entire plane R2,
        // 2. A single point in R2,
        // 3. A line segment represented as {base_point + t * direction_vector | t in [0,1]},
        // 4. An infinite parallelogram represented as {(xi,zeta) | <normal, (xi,zeta)> in I} for
        //    some bounded interval I (which can be a single point, in which case the parallelogram
        //    degenerates into a line), or
        // 5. A bounded polygon (with nonempty interior) represented as a union of vertically
        //    parallel trapezoids (which are in fact further divided into vertical boundary lines
        //    and horizontally-open trapezoids). A vertically parallel trapezoid is represented by a
        //    pair of affine functions of xi so that zeta is in between those two.
        //

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

        namespace xi_zeta_region {
            enum class region_type_t {
                entire_plane,
                single_point,
                line_segment,
                infinite_parallelogram,
                bounded_polygon
            };
            using frac_t = frac<bigint::int_var, bigint::uint_var>;

            class entire_plane {
            public:
                static constexpr region_type_t region_type = region_type_t::entire_plane;
            };

            class single_point {
            public:
                static constexpr region_type_t region_type = region_type_t::single_point;

                frac_t xi;
                frac_t zeta;
            };

            // {base_point + t * direction_vector | t in [0,1]}.
            class line_segment {
            public:
                static constexpr region_type_t region_type = region_type_t::line_segment;

                frac_t base_point_xi;
                frac_t base_point_zeta;
                frac<bigint::uint_var, bigint::uint_var> direction_vector_xi;
                frac_t direction_vector_zeta;
                boundary_type_t left_boundary_type;
                boundary_type_t right_boundary_type;
            };

            // Can be an infinite line when value_gap == 0.
            class infinite_parallelogram {
            public:
                static constexpr region_type_t region_type = region_type_t::infinite_parallelogram;

                // Inward normal vector.
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

                enum class slice_type_t {
                    single_point,
                    vertical_line_segment,
                    horizontally_open_trapezoid
                };

                struct single_point {
                    static constexpr slice_type_t slice_type = slice_type_t::single_point;

                    frac_t const& xi;
                    frac_t const& zeta;

                    constexpr operator xi_zeta_region::single_point() const { return {xi, zeta}; }
                };

                struct vertical_line_segment {
                    static constexpr slice_type_t slice_type = slice_type_t::vertical_line_segment;

                    frac_t const& xi;
                    variable_shape_interval<frac_t, interval_type_t::bounded_open,
                                            interval_type_t::bounded_left_open_right_closed,
                                            interval_type_t::bounded_left_closed_right_open,
                                            interval_type_t::bounded_closed> const zeta_range;
                };

                struct horizontally_open_trapezoid {
                    static constexpr slice_type_t slice_type =
                        slice_type_t::horizontally_open_trapezoid;

                    interval<frac_t const&, interval_type_t::bounded_open> xi_range;

                    frac_t const& lower_boundary_linear_coeff;
                    frac_t const& lower_boundary_constant_coeff;
                    frac_t const& upper_boundary_linear_coeff;
                    frac_t const& upper_boundary_constant_coeff;

                    boundary_type_t const lower_boundary_type;
                    boundary_type_t const upper_boundary_type;
                };

                template <class Functor>
                void for_each_vertical_slice(Functor&& f) const {
                    util::constexpr_assert(!turning_points_.empty());
                    util::constexpr_assert(!boundary_line_pairs_.empty());

                    auto on_turning_point = [&f](auto const& turning_point) {
                        if (turning_point.lower_boundary_type == boundary_type_t::closed) {
                            if (turning_point.upper_boundary_type == boundary_type_t::closed) {
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
                            if (turning_point.upper_boundary_type == boundary_type_t::closed) {
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

                    auto on_last_turning_point = [&](auto const& turning_point) {
                        if (turning_point.lower_zeta == turning_point.upper_zeta) {
                            util::constexpr_assert(
                                turning_point.lower_boundary_type == boundary_type_t::closed &&
                                turning_point.upper_boundary_type == boundary_type_t::closed);
                            f(single_point{turning_point.xi, turning_point.lower_zeta});
                        }
                        else {
                            on_turning_point(turning_point);
                        }
                    };

                    if (left_boundary_type() == boundary_type_t::closed) {
                        on_last_turning_point(turning_points().front());
                    }

                    for (std::size_t turning_point_idx = 1;
                         turning_point_idx < turning_points().size(); ++turning_point_idx) {
                        f(horizontally_open_trapezoid{
                            .xi_range =
                                interval<frac_t const&, interval_type_t::bounded_open>{
                                    turning_points()[turning_point_idx - 1].xi,
                                    turning_points()[turning_point_idx].xi},
                            .lower_boundary_linear_coeff =
                                boundary_line_pairs()[turning_point_idx - 1]
                                    .lower_boundary_linear_coeff,
                            .lower_boundary_constant_coeff =
                                boundary_line_pairs()[turning_point_idx - 1]
                                    .lower_boundary_constant_coeff,
                            .upper_boundary_linear_coeff =
                                boundary_line_pairs()[turning_point_idx - 1]
                                    .upper_boundary_linear_coeff,
                            .upper_boundary_constant_coeff =
                                boundary_line_pairs()[turning_point_idx - 1]
                                    .upper_boundary_constant_coeff,
                            .lower_boundary_type =
                                boundary_line_pairs()[turning_point_idx - 1].lower_boundary_type,
                            .upper_boundary_type =
                                boundary_line_pairs()[turning_point_idx - 1].upper_boundary_type});

                        if (turning_point_idx < turning_points().size() - 1) {
                            on_turning_point(turning_points()[turning_point_idx]);
                        }
                        else if (right_boundary_type() == boundary_type_t::closed) {
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
                    frac_t lower_boundary_linear_coeff;
                    frac_t lower_boundary_constant_coeff;
                    frac_t upper_boundary_linear_coeff;
                    frac_t upper_boundary_constant_coeff;

                    boundary_type_t lower_boundary_type;
                    boundary_type_t upper_boundary_type;
                };

                auto const& turning_points() const noexcept { return turning_points_; }
                auto const& boundary_line_pairs() const noexcept { return boundary_line_pairs_; }
                boundary_type_t left_boundary_type() const noexcept { return left_boundary_type_; }
                boundary_type_t right_boundary_type() const noexcept {
                    return right_boundary_type_;
                }


                template <std::ranges::common_range TurningPoints,
                          std::ranges::common_range BoundaryLinePairs>
                bounded_polygon(TurningPoints&& turning_points,
                                BoundaryLinePairs&& boundary_line_pairs,
                                boundary_type_t left_boundary_type,
                                boundary_type_t right_boundary_type)
                    : turning_points_(turning_points.cbegin(), turning_points.cend()),
                      boundary_line_pairs_(boundary_line_pairs.cbegin(),
                                           boundary_line_pairs.cend()),
                      left_boundary_type_{left_boundary_type},
                      right_boundary_type_{right_boundary_type} {}

            private:
                std::vector<turning_point_info> turning_points_;
                std::vector<boundary_line_pair_info> boundary_line_pairs_;
                boundary_type_t left_boundary_type_;
                boundary_type_t right_boundary_type_;
            };

            using variable_shape_region = std::variant<entire_plane, single_point, line_segment,
                                                       infinite_parallelogram, bounded_polygon>;
        }

        template <class ContinuedFractionGeneratorX, class ContinuedFractionGeneratorY,
                  std::ranges::range RangeOfConstraintSpec>
        xi_zeta_region::variable_shape_region
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

            auto gcd = [](auto const& a, auto const& b) {
                auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                    cntfrc::impl::rational{a, b});
                while (!cf.terminated()) {
                    cf.update();
                }
                return cf.current_partial_fraction().denominator;
            };
            auto lcm = [&gcd](auto const& a, auto const& b) { return (a / gcd(a, b)) * b; };

            // Reduce a fraction.
            auto reduce_fraction = [](auto&& fr) {
                auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                    cntfrc::impl::rational{static_cast<decltype(fr)&&>(fr).numerator,
                                           static_cast<decltype(fr)&&>(fr).denominator});
                while (!cf.terminated()) {
                    cf.update();
                }
                return project_to_rational(cf.current_convergent());
            };

            // Get the reduced form of the number num/den, where num and den are rationals.
            auto reduced_division = [](auto const& num, auto const& den) {
                auto num_int = num.numerator * den.denominator;
                auto den_int = num.denominator * den.numerator;
                if (util::is_strictly_negative(den_int)) {
                    num_int = -std::move(num_int);
                    den_int = -std::move(den_int);
                }
                auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                    cntfrc::impl::rational{num_int, util::abs(den_int)});
                while (!cf.terminated()) {
                    cf.update();
                }
                return project_to_rational(cf.current_convergent());
            };


            ////////////////////////////////////////////////////////////////////////////////////
            // Step 1 - Write the region as the intersection of half-planes.
            ////////////////////////////////////////////////////////////////////////////////////

            struct affine_boundary_info {
                frac_t linear_coeff;
                frac_t constant_coeff;
                boundary_type_t boundary_type;
            };
            std::vector<affine_boundary_info> upper_half_planes; // Lower bounds on zeta.
            std::vector<affine_boundary_info> lower_half_planes; // Upper bounds on zeta.

            struct vertical_boundary_info {
                frac_t threshold;
                boundary_type_t boundary_type;
            };
            std::vector<vertical_boundary_info> right_half_planes; // Lower bounds on xi.
            std::vector<vertical_boundary_info> left_half_planes;  // Upper bounds on xi.

            {
                auto caching_xcf_for_xp = cntfrc::make_caching_generator(xcf.copy());
                auto caching_ycf_for_xp = cntfrc::make_caching_generator(ycf.copy());
                auto caching_xcf_for_yp =
                    cntfrc::caching_generator<ContinuedFractionGeneratorX&>{xcf};
                auto caching_ycf_for_yp =
                    cntfrc::caching_generator<ContinuedFractionGeneratorY&>{ycf};

                // Generate a continued fraction implementation for ax + by + c.
                auto affine_combination = [&](auto& xcf, auto& ycf, frac_t const& coeff_x,
                                              frac_t const& coeff_y, frac_t const& constant) {
                    auto denominator = lcm(coeff_x.denominator, coeff_y.denominator);
                    denominator *= (constant.denominator / gcd(denominator, constant.denominator));

                    return cntfrc::impl::binary_gosper<decltype(xcf), decltype(ycf)>{
                        xcf,
                        ycf,
                        {// numerator
                         0, coeff_x.numerator * (denominator / coeff_x.denominator),
                         coeff_y.numerator * (denominator / coeff_y.denominator),
                         constant.numerator * (denominator / constant.denominator),
                         // denominator
                         0, 0, 0, util::to_signed(denominator)}};
                };

                auto triage_and_push_half_plane =
                    [&](auto const& xi_coeff, auto const& minus_zeta_coeff, auto const& eta_coeff,
                        boundary_type_t bdy_type) {
                        // Right or left half-planes.
                        if (util::is_zero(minus_zeta_coeff.numerator)) {
                            // Ignore the ones with only the constant coefficient.
                            if (!util::is_zero(xi_coeff.numerator)) {
                                auto bdy_info = vertical_boundary_info{
                                    reduced_division(-eta_coeff, xi_coeff), bdy_type};

                                if (util::is_nonnegative(xi_coeff.numerator)) {
                                    right_half_planes.push_back(std::move(bdy_info));
                                }
                                else {
                                    left_half_planes.push_back(std::move(bdy_info));
                                }
                            }
                        }
                        // Upper or lower half-planes.
                        else {
                            auto bdy_info = affine_boundary_info{
                                reduced_division(xi_coeff, minus_zeta_coeff),
                                reduced_division(eta_coeff, minus_zeta_coeff), bdy_type};

                            if (util::is_strictly_negative(minus_zeta_coeff.numerator)) {
                                upper_half_planes.push_back(std::move(bdy_info));
                            }
                            else {
                                lower_half_planes.push_back(std::move(bdy_info));
                            }
                        }
                    };

                auto push_closed_half_plane =
                    [&](bigint::uint_view base_point, util::sign_t sign,
                        floor_constraint_spec::affine_transform const& affine_transform,
                        multiply_add_shift_info const& approx_xp_yp_info) {
                        auto signed_base_point = sign == util::sign_t::positive
                                                     ? util::to_signed(base_point)
                                                     : util::to_negative(base_point);

                        // Since approx_xp_yp_info is supposed to be an approximation of (+-x', y'),
                        // we do not multiply sign for computing floor(n(+-x')+y').
                        auto xi_coeff =
                            frac_t{affine_transform.linear_coeff.xx.numerator * signed_base_point,
                                   affine_transform.linear_coeff.xx.denominator} +
                            affine_transform.linear_coeff.yx;
                        auto minus_zeta_coeff =
                            frac_t{-affine_transform.linear_coeff.xy.numerator * signed_base_point,
                                   affine_transform.linear_coeff.xy.denominator} -
                            affine_transform.linear_coeff.yy;
                        auto eta_coeff =
                            frac_t{affine_transform.constant_coeff_x.numerator * signed_base_point,
                                   affine_transform.constant_coeff_x.denominator} +
                            affine_transform.constant_coeff_y -
                            frac_t{((base_point * approx_xp_yp_info.multiplier +
                                     approx_xp_yp_info.adder) >>
                                    approx_xp_yp_info.shift_amount),
                                   1u};

                        triage_and_push_half_plane(xi_coeff, minus_zeta_coeff, eta_coeff,
                                                   boundary_type_t::closed);
                    };

                auto push_open_half_plane =
                    [&](bigint::uint_view base_point, util::sign_t sign,
                        floor_constraint_spec::affine_transform const& affine_transform,
                        multiply_add_shift_info const& approx_xp_yp_info) {
                        auto signed_base_point = sign == util::sign_t::positive
                                                     ? util::to_signed(base_point)
                                                     : util::to_negative(base_point);

                        // Since approx_xp_yp_info is supposed to be an approximation of (+-x', y'),
                        // we do not multiply sign for computing floor(n(+-x')+y').
                        auto xi_coeff =
                            frac_t{-affine_transform.linear_coeff.xx.numerator * signed_base_point,
                                   affine_transform.linear_coeff.xx.denominator} -
                            affine_transform.linear_coeff.yx;
                        auto minus_zeta_coeff =
                            frac_t{affine_transform.linear_coeff.xy.numerator * signed_base_point,
                                   affine_transform.linear_coeff.xy.denominator} +
                            affine_transform.linear_coeff.yy;
                        auto eta_coeff =
                            frac_t{((base_point * approx_xp_yp_info.multiplier +
                                     approx_xp_yp_info.adder) >>
                                    approx_xp_yp_info.shift_amount) +
                                       1u,
                                   1u} -
                            frac_t{affine_transform.constant_coeff_x.numerator * signed_base_point,
                                   affine_transform.constant_coeff_x.denominator} -
                            affine_transform.constant_coeff_y;

                        triage_and_push_half_plane(xi_coeff, minus_zeta_coeff, eta_coeff,
                                                   boundary_type_t::open);
                    };

                // For the maximization on the left, the half-space at the base point is
                // included.
                auto solve_maximization_on_left =
                    [&](bigint::uint_var base_point, bigint::uint_var max_diff, util::sign_t sign,
                        floor_constraint_spec::affine_transform const& affine_transform,
                        multiply_add_shift_info const& approx_xp_yp_info, auto& approx_xpcf) {
                        while (!util::is_zero(max_diff)) {
                            push_closed_half_plane(base_point, sign, affine_transform,
                                                   approx_xp_yp_info);

                            auto movement =
                                find_extremizers_of_fractional_part(approx_xpcf, max_diff)
                                    .largest_maximizer;
                            approx_xpcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point -= movement;
                            max_diff -= std::move(movement);
                        }

                        push_closed_half_plane(base_point, sign, affine_transform,
                                               approx_xp_yp_info);
                    };

                // For the maximization on the right, the half-space at the base point is
                // not included.
                auto solve_maximization_on_right =
                    [&](bigint::uint_var base_point, bigint::uint_var max_diff, util::sign_t sign,
                        floor_constraint_spec::affine_transform const& affine_transform,
                        multiply_add_shift_info const& approx_xp_yp_info, auto& approx_xpcf) {
                        while (!util::is_zero(max_diff)) {
                            auto movement =
                                find_extremizers_of_fractional_part(approx_xpcf, max_diff)
                                    .smallest_minimizer;
                            approx_xpcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point += movement;
                            max_diff -= std::move(movement);

                            push_closed_half_plane(base_point, sign, affine_transform,
                                                   approx_xp_yp_info);
                        }
                    };

                // For the minimization on the left, the half-space at the base point is
                // included.
                auto solve_minimization_on_left =
                    [&](bigint::uint_var base_point, bigint::uint_var max_diff, util::sign_t sign,
                        floor_constraint_spec::affine_transform const& affine_transform,
                        multiply_add_shift_info const& approx_xp_yp_info, auto& approx_xpcf) {
                        while (!util::is_zero(max_diff)) {
                            push_open_half_plane(base_point, sign, affine_transform,
                                                 approx_xp_yp_info);

                            auto movement =
                                find_extremizers_of_fractional_part(approx_xpcf, max_diff)
                                    .smallest_minimizer;
                            approx_xpcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point -= movement;
                            max_diff -= std::move(movement);
                        }

                        push_open_half_plane(base_point, sign, affine_transform, approx_xp_yp_info);
                    };

                // For the minimization on the right, the half-space at the base point is
                // not included.
                auto solve_minimization_on_right =
                    [&](bigint::uint_var base_point, bigint::uint_var max_diff, util::sign_t sign,
                        floor_constraint_spec::affine_transform const& affine_transform,
                        multiply_add_shift_info const& approx_xp_yp_info, auto& approx_xpcf) {
                        while (!util::is_zero(max_diff)) {
                            auto movement =
                                find_extremizers_of_fractional_part(approx_xpcf, max_diff)
                                    .largest_maximizer;
                            approx_xpcf.rewind();
                            movement *= util::div_floor(max_diff, movement);
                            base_point += movement;
                            max_diff -= std::move(movement);

                            push_open_half_plane(base_point, sign, affine_transform,
                                                 approx_xp_yp_info);
                        }
                    };

                auto process_single_sign_interval =
                    [&](util::sign_t sign, multiply_add_shift_info const& approx_xp_yp_info,
                        nrange_t const& nrange,
                        floor_constraint_spec::affine_transform const& affine_coeff) {
                        auto approx_xpcf = cntfrc::make_caching_generator(
                            cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                                   cntfrc::interval_tracker>(cntfrc::impl::rational{
                                approx_xp_yp_info.multiplier,
                                bigint::uint_var::power_of_2(approx_xp_yp_info.shift_amount)}));

                        auto approx_ypcf =
                            cntfrc::make_generator<cntfrc::interval_tracker>(cntfrc::impl::rational{
                                approx_xp_yp_info.adder,
                                bigint::uint_var::power_of_2(approx_xp_yp_info.shift_amount)});

                        auto base_points =
                            find_extremizers_of_fractional_part(approx_xpcf, approx_ypcf, nrange);
                        approx_xpcf.rewind();

                        solve_maximization_on_left(
                            util::abs(base_points.smallest_minimizer),
                            util::abs(base_points.smallest_minimizer - nrange.lower_bound()), sign,
                            affine_coeff, approx_xp_yp_info, approx_xpcf);
                        solve_maximization_on_right(
                            util::abs(base_points.smallest_minimizer),
                            util::abs(nrange.upper_bound() - base_points.smallest_minimizer), sign,
                            affine_coeff, approx_xp_yp_info, approx_xpcf);
                        solve_minimization_on_left(
                            util::abs(base_points.largest_maximizer),
                            util::abs(base_points.largest_maximizer - nrange.lower_bound()), sign,
                            affine_coeff, approx_xp_yp_info, approx_xpcf);
                        solve_minimization_on_right(
                            util::abs(base_points.largest_maximizer),
                            util::abs(nrange.upper_bound() - base_points.largest_maximizer), sign,
                            affine_coeff, approx_xp_yp_info, approx_xpcf);
                    };

                // For each constraint.
                for (auto const& constraint_spec : range_of_constraint_specs) {
                    // Positive n's.
                    if (util::is_strictly_positive(constraint_spec.nrange.upper_bound())) {
                        auto nrange =
                            nrange_t{util::is_nonpositive(constraint_spec.nrange.lower_bound())
                                         ? 1
                                         : constraint_spec.nrange.lower_bound(),
                                     constraint_spec.nrange.upper_bound()};

                        // Read affine coefficients and initialize the transformed variables
                        // (x',y') accordingly. Compute a good enough approximation of (x',y') for
                        // future computations; we only need to care about numbers of the form
                        // floor(nx') and floor(nx' + y').
                        auto const approx_xp_yp_info = find_simultaneous_multiply_add_shift(
                            cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                                   cntfrc::interval_tracker>(
                                affine_combination(caching_xcf_for_xp, caching_ycf_for_xp,
                                                   constraint_spec.affine_coeff.linear_coeff.xx,
                                                   constraint_spec.affine_coeff.linear_coeff.xy,
                                                   constraint_spec.affine_coeff.constant_coeff_x)),
                            cntfrc::make_generator<cntfrc::interval_tracker>(
                                affine_combination(caching_xcf_for_yp, caching_ycf_for_yp,
                                                   constraint_spec.affine_coeff.linear_coeff.yx,
                                                   constraint_spec.affine_coeff.linear_coeff.yy,
                                                   constraint_spec.affine_coeff.constant_coeff_y)),
                            nrange);

                        caching_xcf_for_xp.rewind();
                        caching_ycf_for_xp.rewind();
                        caching_xcf_for_yp.rewind();
                        caching_ycf_for_yp.rewind();

                        process_single_sign_interval(util::sign_t::positive, approx_xp_yp_info,
                                                     nrange, constraint_spec.affine_coeff);
                    }

                    // Negative n's.
                    if (util::is_strictly_negative(constraint_spec.nrange.lower_bound())) {
                        auto nrange =
                            nrange_t{util::is_nonnegative(constraint_spec.nrange.upper_bound())
                                         ? 1
                                         : -constraint_spec.nrange.upper_bound(),
                                     -constraint_spec.nrange.lower_bound()};

                        // Read affine coefficients and initialize the transformed variables
                        // (x',y') accordingly. Compute a good enough approximation of (-x',y') for
                        // future computations; we only need to care about numbers of the form
                        // floor(n(-x')) and floor(n(-x') + y').
                        auto const approx_xp_yp_info = find_simultaneous_multiply_add_shift(
                            cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                                   cntfrc::interval_tracker>(
                                affine_combination(caching_xcf_for_xp, caching_ycf_for_xp,
                                                   -constraint_spec.affine_coeff.linear_coeff.xx,
                                                   -constraint_spec.affine_coeff.linear_coeff.xy,
                                                   -constraint_spec.affine_coeff.constant_coeff_x)),
                            cntfrc::make_generator<cntfrc::interval_tracker>(
                                affine_combination(caching_xcf_for_yp, caching_ycf_for_yp,
                                                   constraint_spec.affine_coeff.linear_coeff.yx,
                                                   constraint_spec.affine_coeff.linear_coeff.yy,
                                                   constraint_spec.affine_coeff.constant_coeff_y)),
                            nrange);

                        caching_xcf_for_xp.rewind();
                        caching_ycf_for_xp.rewind();
                        caching_xcf_for_yp.rewind();
                        caching_ycf_for_yp.rewind();

                        process_single_sign_interval(util::sign_t::negative, approx_xp_yp_info,
                                                     nrange, constraint_spec.affine_coeff);
                    }

                    // Zero.
                    if (constraint_spec.nrange.contains(0)) {
                        auto const floor_yp = [&] {
                            auto cf = cntfrc::make_generator<cntfrc::partial_fraction_tracker>(
                                affine_combination(caching_xcf_for_yp, caching_ycf_for_yp,
                                                   constraint_spec.affine_coeff.linear_coeff.yx,
                                                   constraint_spec.affine_coeff.linear_coeff.yy,
                                                   constraint_spec.affine_coeff.constant_coeff_y));
                            cf.update();
                            caching_xcf_for_yp.rewind();
                            caching_ycf_for_yp.rewind();
                            return cf.current_partial_fraction().denominator;
                        }();

                        triage_and_push_half_plane(constraint_spec.affine_coeff.linear_coeff.yx,
                                                   -constraint_spec.affine_coeff.linear_coeff.yy,
                                                   constraint_spec.affine_coeff.constant_coeff_y -
                                                       frac_t{floor_yp, 1u},
                                                   boundary_type_t::closed);

                        triage_and_push_half_plane(
                            -constraint_spec.affine_coeff.linear_coeff.yx,
                            constraint_spec.affine_coeff.linear_coeff.yy,
                            frac_t{floor_yp + 1u, 1u} -
                                constraint_spec.affine_coeff.constant_coeff_y,
                            boundary_type_t::open);
                    }
                } // for (auto const& constraint_spec : range_of_constraint_specs)
            }


            ////////////////////////////////////////////////////////////////////////////////////
            // Step 2 - Find the intersection of lower bounds for zeta and upper bounds for
            // zeta, respectively. The intersection is described in terms of a piecewise affine
            // convex/concave function.
            ////////////////////////////////////////////////////////////////////////////////////

            // Intersection for half-planes with vertical boundaries.
            auto xi_min_itr =
                std::ranges::max_element(right_half_planes, [](auto const& lhs, auto const& rhs) {
                    if (lhs.threshold < rhs.threshold) {
                        return true;
                    }
                    else if (lhs.threshold > rhs.threshold) {
                        return false;
                    }
                    else {
                        return lhs.boundary_type == boundary_type_t::closed &&
                               rhs.boundary_type == boundary_type_t::open;
                    }
                });
            auto xi_max_itr =
                std::ranges::min_element(left_half_planes, [](auto const& lhs, auto const& rhs) {
                    if (lhs.threshold < rhs.threshold) {
                        return true;
                    }
                    else if (lhs.threshold > rhs.threshold) {
                        return false;
                    }
                    else {
                        return lhs.boundary_type == boundary_type_t::open &&
                               rhs.boundary_type == boundary_type_t::closed;
                    }
                });

            // Early return if there are only constraints on xi.
            if (lower_half_planes.empty()) {
                util::constexpr_assert(upper_half_planes.empty());

                // No constraint: return the entire plane.
                if (right_half_planes.empty()) {
                    util::constexpr_assert(left_half_planes.empty());
                    return xi_zeta_region::entire_plane{};
                }
                util::constexpr_assert(!left_half_planes.empty());

                // Otherwise, the result is a vertically parallel infinite parallelogram.
                auto denominator =
                    lcm(xi_min_itr->threshold.denominator, xi_max_itr->threshold.denominator);
                auto min_value = xi_min_itr->threshold.numerator *
                                 (denominator / xi_min_itr->threshold.denominator);
                auto max_value = xi_max_itr->threshold.numerator *
                                 (denominator / xi_max_itr->threshold.denominator);
                util::constexpr_assert(min_value <= max_value);

                return xi_zeta_region::infinite_parallelogram{
                    .xi_coeff = util::to_signed(denominator),
                    .zeta_coeff = 0u,
                    .min_value = min_value,
                    .value_gap = util::abs(max_value - min_value),
                    .lower_boundary_type = xi_min_itr->boundary_type,
                    .upper_boundary_type = xi_max_itr->boundary_type};
            }
            util::constexpr_assert(!upper_half_planes.empty());


            struct piecewise_affine {
                std::vector<affine_boundary_info> boundary_lines;
                std::vector<frac_t> turning_points;
                std::vector<boundary_type_t> turning_point_types;
            };

            // Find the max/min of affine functions by finding the convex envelope of the convex
            // conjugate.
            enum class bounding_direction_t : bool { from_below, from_above };
            auto compute_one_sided_intersection =
                [&](bounding_direction_t bounding_direction) -> piecewise_affine {
                auto& half_planes = bounding_direction == bounding_direction_t::from_below
                                        ? upper_half_planes
                                        : lower_half_planes;

                std::ranges::sort(
                    half_planes, [bounding_direction](auto const& lhs, auto const& rhs) {
                        if (lhs.linear_coeff < rhs.linear_coeff) {
                            return bounding_direction == bounding_direction_t::from_below ? true
                                                                                          : false;
                        }
                        else if (lhs.linear_coeff > rhs.linear_coeff) {
                            return bounding_direction == bounding_direction_t::from_below ? false
                                                                                          : true;
                        }
                        else if (lhs.constant_coeff < rhs.constant_coeff) {
                            return bounding_direction == bounding_direction_t::from_below ? true
                                                                                          : false;
                        }
                        else if (lhs.constant_coeff > rhs.constant_coeff) {
                            return bounding_direction == bounding_direction_t::from_below ? false
                                                                                          : true;
                        }
                        else {
                            return lhs.boundary_type == boundary_type_t::open &&
                                   rhs.boundary_type == boundary_type_t::closed;
                        }
                    });

                piecewise_affine result;

                // Andrew's monotone chain algorithm.
                auto curr_itr = half_planes.cbegin();
                result.boundary_lines.push_back(*curr_itr);

                // Push the second point.
                while (true) {
                    if (++curr_itr == half_planes.cend()) {
                        // Only one affine segment.
                        return result;
                    }

                    if (curr_itr->linear_coeff == result.boundary_lines.back().linear_coeff) {
                        continue;
                    }
                    else {
                        // Found the second point.
                        result.boundary_lines.push_back(*curr_itr);

                        // Remember inclusiveness of the turning point.
                        result.turning_point_types.push_back(
                            result.boundary_lines.front().boundary_type ==
                                        boundary_type_t::closed &&
                                    curr_itr->boundary_type == boundary_type_t::closed
                                ? boundary_type_t::closed
                                : boundary_type_t::open);
                        break;
                    }
                }

                // Main iteration.
                while (++curr_itr != half_planes.cend()) {
                    if (curr_itr->linear_coeff == result.boundary_lines.back().linear_coeff) {
                        continue;
                    }

                    // Found a new point. Check if this makes a counterclockwise turn with the last
                    // two. Remember that the point we are looking are of the form
                    // (linear_coeff, -constant_coeff).
                    bool colinear = false;
                    while (result.boundary_lines.size() >= 2) {
                        auto const& first = result.boundary_lines[result.boundary_lines.size() - 2];
                        auto const& second =
                            result.boundary_lines[result.boundary_lines.size() - 1];
                        auto const& third = *curr_itr;
                        auto const cyclic_constant = (first.constant_coeff * second.linear_coeff +
                                                      second.constant_coeff * third.linear_coeff +
                                                      third.constant_coeff * first.linear_coeff) -
                                                     (first.linear_coeff * second.constant_coeff +
                                                      second.linear_coeff * third.constant_coeff +
                                                      third.linear_coeff * first.constant_coeff);

                        // When three points are colinear.
                        if (util::is_zero(cyclic_constant.numerator)) {
                            colinear = true;
                            break;
                        }
                        // When the turn is strictly valid.
                        else if (util::is_nonnegative(cyclic_constant.numerator)) {
                            break;
                        }
                        // When the turn is strictly invalid.
                        else {
                            // Remove the last point.
                            result.boundary_lines.pop_back();
                            result.turning_point_types.pop_back();
                        }
                    } // while (result.boundary_lines.size() >= 2)

                    if (colinear) {
                        // Replace the last with the newly found point, but reflect the
                        // inclusiveness of the last point.
                        result.boundary_lines.pop_back();
                        result.boundary_lines.push_back(*curr_itr);
                        result.turning_point_types.back() =
                            result.turning_point_types.back() == boundary_type_t::closed &&
                                    curr_itr->boundary_type == boundary_type_t::closed
                                ? boundary_type_t::closed
                                : boundary_type_t::open;
                    }
                    else {
                        // Push the newly found point.
                        result.boundary_lines.push_back(*curr_itr);

                        // Remember inclusiveness of the turning point.
                        result.turning_point_types.push_back(
                            result.boundary_lines.back().boundary_type == boundary_type_t::closed &&
                                    curr_itr->boundary_type == boundary_type_t::closed
                                ? boundary_type_t::closed
                                : boundary_type_t::open);
                    }
                } // while (++curr_itr != half_planes.cend())

                // Compute turning points and return.
                for (std::size_t idx = 1; idx < result.boundary_lines.size(); ++idx) {
                    result.turning_points.push_back(
                        reduced_division(result.boundary_lines[idx].constant_coeff -
                                             result.boundary_lines[idx - 1].constant_coeff,
                                         result.boundary_lines[idx - 1].linear_coeff -
                                             result.boundary_lines[idx].linear_coeff));
                }
                return result;
            };

            auto lower_envelope = compute_one_sided_intersection(bounding_direction_t::from_below);
            auto upper_envelope = compute_one_sided_intersection(bounding_direction_t::from_above);


            ////////////////////////////////////////////////////////////////////////////////////
            // Step 3 - Find the intersection of the region between convex/concave functions.
            ////////////////////////////////////////////////////////////////////////////////////

            // Early return for infinite parallelogram.
            if (lower_envelope.turning_points.empty() && upper_envelope.turning_points.empty() &&
                right_half_planes.empty()) {
                util::constexpr_assert(lower_envelope.boundary_lines.size() == 1);
                util::constexpr_assert(upper_envelope.boundary_lines.size() == 1);

                auto const& lower_bdy = lower_envelope.boundary_lines.front();
                auto const& upper_bdy = upper_envelope.boundary_lines.front();

                util::constexpr_assert(lower_bdy.linear_coeff == upper_bdy.linear_coeff);
                util::constexpr_assert(left_half_planes.empty());

                auto denominator =
                    lcm(lower_bdy.constant_coeff.denominator, upper_bdy.constant_coeff.denominator);
                denominator *= lower_bdy.linear_coeff.denominator /
                               gcd(denominator, lower_bdy.linear_coeff.denominator);
                auto min_value = lower_bdy.constant_coeff.numerator *
                                 (denominator / lower_bdy.constant_coeff.denominator);
                auto max_value = upper_bdy.constant_coeff.numerator *
                                 (denominator / upper_bdy.constant_coeff.denominator);

                util::constexpr_assert(min_value < max_value ||
                                       (min_value == max_value &&
                                        lower_bdy.boundary_type == boundary_type_t::closed &&
                                        upper_bdy.boundary_type == boundary_type_t::closed));

                return xi_zeta_region::infinite_parallelogram{
                    .xi_coeff = lower_bdy.linear_coeff.numerator *
                                (denominator / lower_bdy.linear_coeff.denominator),
                    .zeta_coeff = denominator,
                    .min_value = min_value,
                    .value_gap = util::abs(max_value - min_value),
                    .lower_boundary_type = lower_bdy.boundary_type,
                    .upper_boundary_type = upper_bdy.boundary_type};
            } // Early return for infinite parallelogram.

            // Early return for vertical line segments.
            if (!right_half_planes.empty()) {
                util::constexpr_assert(!left_half_planes.empty());
                if (xi_min_itr->threshold == xi_max_itr->threshold) {
                    util::constexpr_assert(xi_min_itr->boundary_type == boundary_type_t::closed &&
                                           xi_max_itr->boundary_type == boundary_type_t::closed);

                    struct endpoint_t {
                        frac_t zeta;
                        boundary_type_t boundary_type;
                    };

                    auto compute_endpoint = [&](auto const& envelope) -> endpoint_t {
                        if (envelope.turning_points.empty()) {
                            util::constexpr_assert(envelope.boundary_lines.size() == 1);
                            return {reduce_fraction(envelope.boundary_lines.front().linear_coeff *
                                                        xi_min_itr->threshold +
                                                    envelope.boundary_lines.front().constant_coeff),
                                    envelope.boundary_lines.front().boundary_type};
                        }
                        for (std::size_t idx = 0; idx < envelope.turning_points.size(); ++idx) {
                            auto const cmp_result =
                                xi_min_itr->threshold <=> envelope.turning_points[idx];
                            if (cmp_result < 0) {
                                return {
                                    reduce_fraction(envelope.boundary_lines[idx].linear_coeff *
                                                        xi_min_itr->threshold +
                                                    envelope.boundary_lines[idx].constant_coeff),
                                    envelope.boundary_lines[idx].boundary_type};
                            }
                            else if (cmp_result == 0) {
                                return {
                                    reduce_fraction(envelope.boundary_lines[idx].linear_coeff *
                                                        xi_min_itr->threshold +
                                                    envelope.boundary_lines[idx].constant_coeff),
                                    envelope.turning_point_types[idx]};
                            }
                        }
                        return {reduce_fraction(envelope.boundary_lines.back().linear_coeff *
                                                    xi_min_itr->threshold +
                                                envelope.boundary_lines.back().constant_coeff),
                                envelope.boundary_lines.back().boundary_type};
                    };

                    auto lower_end = compute_endpoint(lower_envelope);
                    auto upper_end = compute_endpoint(upper_envelope);
                    auto gap = reduce_fraction(upper_end.zeta - lower_end.zeta);

                    return xi_zeta_region::line_segment{
                        .base_point_xi = xi_min_itr->threshold,
                        .base_point_zeta = lower_end.zeta,
                        .direction_vector_xi = frac<bigint::uint_var, bigint::uint_var>{0u, 1u},
                        .direction_vector_zeta = std::move(gap),
                        .left_boundary_type = lower_end.boundary_type,
                        .right_boundary_type = upper_end.boundary_type};
                }
            } // Early return for vertical line segments.

            // For all remaining cases, the result must be a (possibly degenerate) bounded polygon.
            std::vector<xi_zeta_region::bounded_polygon::turning_point_info> turning_points;
            std::vector<xi_zeta_region::bounded_polygon::boundary_line_pair_info>
                boundary_line_pairs;

            bool const infinite_min_xi = lower_envelope.boundary_lines.front().linear_coeff >=
                                         upper_envelope.boundary_lines.front().linear_coeff;
            bool const infinite_max_xi = lower_envelope.boundary_lines.back().linear_coeff <=
                                         upper_envelope.boundary_lines.back().linear_coeff;

            {
                bool found_min_xi = infinite_min_xi;

                // The case of parallel affine functions.
                if (infinite_min_xi && infinite_max_xi) {
                    util::constexpr_assert(lower_envelope.boundary_lines.size() == 1 &&
                                           upper_envelope.boundary_lines.size() == 1 &&
                                           lower_envelope.boundary_lines.front().linear_coeff ==
                                               upper_envelope.boundary_lines.front().linear_coeff);
                    boundary_line_pairs.push_back(
                        xi_zeta_region::bounded_polygon::boundary_line_pair_info{
                            .lower_boundary_linear_coeff =
                                lower_envelope.boundary_lines.front().linear_coeff,
                            .lower_boundary_constant_coeff =
                                lower_envelope.boundary_lines.front().constant_coeff,
                            .upper_boundary_linear_coeff =
                                upper_envelope.boundary_lines.front().linear_coeff,
                            .upper_boundary_constant_coeff =
                                upper_envelope.boundary_lines.front().constant_coeff,
                            .lower_boundary_type =
                                lower_envelope.boundary_lines.front().boundary_type,
                            .upper_boundary_type =
                                upper_envelope.boundary_lines.front().boundary_type});

                    util::constexpr_assert(!right_half_planes.empty() && !left_half_planes.empty());
                    turning_points.push_back(xi_zeta_region::bounded_polygon::turning_point_info{
                        .xi = xi_min_itr->threshold,
                        .lower_zeta =
                            reduce_fraction(lower_envelope.boundary_lines.front().linear_coeff *
                                                xi_min_itr->threshold +
                                            lower_envelope.boundary_lines.front().constant_coeff),
                        .upper_zeta =
                            reduce_fraction(upper_envelope.boundary_lines.front().linear_coeff *
                                                xi_min_itr->threshold +
                                            upper_envelope.boundary_lines.front().constant_coeff),
                        .lower_boundary_type = lower_envelope.boundary_lines.front().boundary_type,
                        .upper_boundary_type =
                            upper_envelope.boundary_lines.front().boundary_type});
                    turning_points.push_back(xi_zeta_region::bounded_polygon::turning_point_info{
                        .xi = xi_max_itr->threshold,
                        .lower_zeta =
                            reduce_fraction(lower_envelope.boundary_lines.front().linear_coeff *
                                                xi_max_itr->threshold +
                                            lower_envelope.boundary_lines.front().constant_coeff),
                        .upper_zeta =
                            reduce_fraction(upper_envelope.boundary_lines.front().linear_coeff *
                                                xi_max_itr->threshold +
                                            upper_envelope.boundary_lines.front().constant_coeff),
                        .lower_boundary_type = lower_envelope.boundary_lines.front().boundary_type,
                        .upper_boundary_type =
                            upper_envelope.boundary_lines.front().boundary_type});

                    return xi_zeta_region::bounded_polygon{turning_points, boundary_line_pairs,
                                                           xi_min_itr->boundary_type,
                                                           xi_max_itr->boundary_type};
                }

                // Sweep from left to right. By the construction, these arrays should be sorted
                // according to the xi-coordinate, from left to right.
                util::constexpr_assert(!lower_envelope.boundary_lines.empty() &&
                                       !upper_envelope.boundary_lines.empty());

                enum class position_in_interval_t {
                    exterior,
                    interior,
                    left_endpoint,
                    right_endpoint
                };
                auto compute_position_in_domain = [&](frac_t const& xi, auto const& envelope,
                                                      std::size_t idx) -> position_in_interval_t {
                    if (idx != 0) {
                        auto cmp_result = envelope.turning_points[idx - 1] <=> xi;
                        if (cmp_result > 0) {
                            return position_in_interval_t::exterior;
                        }
                        else if (cmp_result == 0) {
                            return position_in_interval_t::left_endpoint;
                        }
                    }
                    if (idx != envelope.boundary_lines.size() - 1) {
                        auto cmp_result = envelope.turning_points[idx] <=> xi;
                        if (cmp_result < 0) {
                            return position_in_interval_t::exterior;
                        }
                        else if (cmp_result == 0) {
                            return position_in_interval_t::right_endpoint;
                        }
                    }
                    return position_in_interval_t::interior;
                };

                std::size_t lower_idx = 0;
                std::size_t upper_idx = 0;

                auto push_slice = [&] {
                    boundary_line_pairs.push_back(
                        xi_zeta_region::bounded_polygon::boundary_line_pair_info{
                            .lower_boundary_linear_coeff =
                                lower_envelope.boundary_lines[lower_idx].linear_coeff,
                            .lower_boundary_constant_coeff =
                                lower_envelope.boundary_lines[lower_idx].constant_coeff,
                            .upper_boundary_linear_coeff =
                                upper_envelope.boundary_lines[upper_idx].linear_coeff,
                            .upper_boundary_constant_coeff =
                                upper_envelope.boundary_lines[upper_idx].constant_coeff,
                            .lower_boundary_type =
                                lower_envelope.boundary_lines[lower_idx].boundary_type,
                            .upper_boundary_type =
                                upper_envelope.boundary_lines[upper_idx].boundary_type});

                    auto const cmp_result = lower_idx == lower_envelope.boundary_lines.size() - 1
                                                ? std::strong_ordering::greater
                                            : upper_idx == upper_envelope.boundary_lines.size() - 1
                                                ? std::strong_ordering::less
                                                : (lower_envelope.turning_points[lower_idx] <=>
                                                   upper_envelope.turning_points[upper_idx]);

                    util::constexpr_assert(lower_idx != lower_envelope.boundary_lines.size() - 1 ||
                                           upper_idx != upper_envelope.boundary_lines.size() - 1);

                    auto xi = cmp_result <= 0 ? lower_envelope.turning_points[lower_idx]
                                              : upper_envelope.turning_points[upper_idx];
                    auto lower_zeta =
                        reduce_fraction(lower_envelope.boundary_lines[lower_idx].linear_coeff * xi +
                                        lower_envelope.boundary_lines[lower_idx].constant_coeff);
                    auto upper_zeta =
                        reduce_fraction(upper_envelope.boundary_lines[upper_idx].linear_coeff * xi +
                                        upper_envelope.boundary_lines[upper_idx].constant_coeff);
                    auto lower_boundary_type =
                        lower_envelope.boundary_lines[lower_idx].boundary_type;
                    auto upper_boundary_type =
                        upper_envelope.boundary_lines[upper_idx].boundary_type;

                    util::constexpr_assert(lower_zeta < upper_zeta);

                    if (cmp_result <= 0) {
                        if (lower_envelope.turning_point_types[lower_idx] ==
                            boundary_type_t::open) {
                            lower_boundary_type = boundary_type_t::open;
                        }
                        ++lower_idx;
                    }
                    if (cmp_result >= 0) {
                        if (upper_envelope.turning_point_types[upper_idx] ==
                            boundary_type_t::open) {
                            upper_boundary_type = boundary_type_t::open;
                        }
                        ++upper_idx;
                    }

                    turning_points.push_back(xi_zeta_region::bounded_polygon::turning_point_info{
                        .xi = std::move(xi),
                        .lower_zeta = std::move(lower_zeta),
                        .upper_zeta = std::move(upper_zeta),
                        .lower_boundary_type = lower_boundary_type,
                        .upper_boundary_type = upper_boundary_type});
                };

                while (true) {
                    auto const& lower_bdy = lower_envelope.boundary_lines[lower_idx];
                    auto const& upper_bdy = upper_envelope.boundary_lines[upper_idx];

                    // If two lines are not parallel.
                    if (lower_bdy.linear_coeff != upper_bdy.linear_coeff) {
                        auto intersection_xi =
                            reduced_division(upper_bdy.constant_coeff - lower_bdy.constant_coeff,
                                             lower_bdy.linear_coeff - upper_bdy.linear_coeff);

                        auto const position_in_lower_domain =
                            compute_position_in_domain(intersection_xi, lower_envelope, lower_idx);
                        auto const position_in_upper_domain =
                            compute_position_in_domain(intersection_xi, upper_envelope, upper_idx);

                        // If the intersection happens within the domains.
                        if (position_in_lower_domain != position_in_interval_t::exterior &&
                            position_in_upper_domain != position_in_interval_t::exterior) {
                            util::constexpr_assert(position_in_lower_domain !=
                                                   position_in_interval_t::left_endpoint);
                            util::constexpr_assert(position_in_upper_domain !=
                                                   position_in_interval_t::left_endpoint);

                            // If the found intersection is the left-most point.
                            if (!found_min_xi) {
                                // We branch on whether the lower envelope will:
                                //   - go strictly above the upper envelope, or
                                //   - go together with the upper envelope, or
                                //   - go strictly below the upper envelope.
                                // To do this branch, all that matters is the right-hand Dini
                                // derivatives of the lower/upper envelopes.

                                // The current slope if not the right endpoint, the next slope
                                // otherwise.
                                auto lower_slope_from_right =
                                    position_in_lower_domain !=
                                            position_in_interval_t::right_endpoint
                                        ? lower_bdy.linear_coeff
                                        : lower_envelope.boundary_lines[lower_idx + 1].linear_coeff;
                                auto upper_slope_from_right =
                                    position_in_upper_domain !=
                                            position_in_interval_t::right_endpoint
                                        ? upper_bdy.linear_coeff
                                        : upper_envelope.boundary_lines[upper_idx + 1].linear_coeff;

                                auto lower_vs_upper =
                                    lower_slope_from_right <=> upper_slope_from_right;

                                // The lower goes strictly above the upper.
                                if (lower_vs_upper > 0) {
                                    // The intersection is a single point.
                                    // Make sure the set is not empty.
                                    util::constexpr_assert(
                                        lower_envelope.boundary_lines[lower_idx].boundary_type ==
                                            boundary_type_t::closed &&
                                        upper_envelope.boundary_lines[upper_idx].boundary_type ==
                                            boundary_type_t::closed &&
                                        (position_in_lower_domain !=
                                             position_in_interval_t::right_endpoint ||
                                         lower_envelope.turning_point_types[lower_idx] ==
                                             boundary_type_t::closed) &&
                                        (position_in_upper_domain !=
                                             position_in_interval_t::right_endpoint ||
                                         upper_envelope.turning_point_types[upper_idx] ==
                                             boundary_type_t::closed));
                                    util::constexpr_assert(
                                        (right_half_planes.empty() ||
                                         (xi_min_itr->threshold <= intersection_xi &&
                                          (xi_min_itr->boundary_type == boundary_type_t::closed ||
                                           xi_min_itr->threshold < intersection_xi))) &&
                                        (left_half_planes.empty() ||
                                         (xi_max_itr->threshold >= intersection_xi &&
                                          (xi_max_itr->boundary_type == boundary_type_t::closed ||
                                           xi_max_itr->threshold > intersection_xi))));

                                    return xi_zeta_region::single_point{
                                        intersection_xi,
                                        reduce_fraction(lower_bdy.linear_coeff * intersection_xi +
                                                        lower_bdy.constant_coeff)};
                                }
                                // The lower goes together with the upper.
                                else if (lower_vs_upper == 0) {
                                    // The intersection is a non-vertical line segment.
                                    // We need to find the left and the right endpoints.
                                    util::constexpr_assert(lower_bdy.boundary_type ==
                                                           boundary_type_t::closed);
                                    util::constexpr_assert(upper_bdy.boundary_type ==
                                                           boundary_type_t::closed);

                                    struct endpoint_t {
                                        frac_t xi;
                                        boundary_type_t boundary_type;
                                    };

                                    auto xi_min = [&]() -> endpoint_t {
                                        auto boundary_type = boundary_type_t::closed;
                                        if (position_in_lower_domain ==
                                            position_in_interval_t::right_endpoint) {
                                            if (lower_envelope.turning_point_types[lower_idx] ==
                                                boundary_type_t::open) {
                                                boundary_type = boundary_type_t::open;
                                            }
                                        }
                                        if (position_in_upper_domain ==
                                            position_in_interval_t::right_endpoint) {
                                            if (upper_envelope.turning_point_types[lower_idx] ==
                                                boundary_type_t::open) {
                                                boundary_type = boundary_type_t::open;
                                            }
                                        }

                                        if (right_half_planes.empty()) {
                                            return {std::move(intersection_xi), boundary_type};
                                        }
                                        else {
                                            auto const cmp_result =
                                                intersection_xi <=> xi_min_itr->threshold;
                                            if (cmp_result > 0) {
                                                return {std::move(intersection_xi), boundary_type};
                                            }
                                            else if (cmp_result == 0) {
                                                if (xi_min_itr->boundary_type ==
                                                    boundary_type_t::open) {
                                                    boundary_type = boundary_type_t::open;
                                                }
                                                return {std::move(intersection_xi), boundary_type};
                                            }
                                            else {
                                                return {std::move(xi_min_itr->threshold),
                                                        xi_min_itr->boundary_type};
                                            }
                                        }
                                    }();

                                    auto xi_max = [&]() -> endpoint_t {
                                        // If the current segment is the last for both the lower and
                                        // the upper envelope, then the intersection continues to
                                        // infinity.
                                        if (lower_idx == lower_envelope.boundary_lines.size() - 1 &&
                                            upper_idx == upper_envelope.boundary_lines.size() - 1) {
                                            util::constexpr_assert(!left_half_planes.empty());
                                            return {std::move(xi_max_itr->threshold),
                                                    xi_max_itr->boundary_type};
                                        }

                                        auto min_between_lower_and_upper = [&]() -> endpoint_t {
                                            if (lower_idx ==
                                                lower_envelope.boundary_lines.size() - 1) {
                                                return {
                                                    upper_envelope.turning_points[upper_idx],
                                                    upper_envelope.turning_point_types[upper_idx]};
                                            }
                                            if (upper_idx ==
                                                upper_envelope.boundary_lines.size() - 1) {
                                                return {
                                                    lower_envelope.turning_points[lower_idx],
                                                    lower_envelope.turning_point_types[lower_idx]};
                                            }

                                            auto cmp_result =
                                                lower_envelope.turning_points[lower_idx] <=>
                                                upper_envelope.turning_points[upper_idx];
                                            if (cmp_result < 0) {
                                                return {
                                                    lower_envelope.turning_points[lower_idx],
                                                    lower_envelope.turning_point_types[lower_idx]};
                                            }
                                            else if (cmp_result > 0) {
                                                return {
                                                    upper_envelope.turning_points[upper_idx],
                                                    upper_envelope.turning_point_types[upper_idx]};
                                            }
                                            else {
                                                return {
                                                    lower_envelope.turning_points[lower_idx],
                                                    lower_envelope.turning_point_types[lower_idx] ==
                                                                boundary_type_t::open ||
                                                            upper_envelope.turning_point_types
                                                                    [upper_idx] ==
                                                                boundary_type_t::open
                                                        ? boundary_type_t::open
                                                        : boundary_type_t::closed};
                                            }
                                        }();

                                        if (left_half_planes.empty()) {
                                            return min_between_lower_and_upper;
                                        }

                                        auto cmp_result = min_between_lower_and_upper.xi <=>
                                                          xi_max_itr->threshold;
                                        if (cmp_result < 0) {
                                            return min_between_lower_and_upper;
                                        }
                                        else if (cmp_result > 0) {
                                            return {xi_max_itr->threshold,
                                                    xi_max_itr->boundary_type};
                                        }
                                        else {
                                            return {min_between_lower_and_upper.xi,
                                                    min_between_lower_and_upper.boundary_type ==
                                                                boundary_type_t::open ||
                                                            xi_max_itr->boundary_type ==
                                                                boundary_type_t::open
                                                        ? boundary_type_t::open
                                                        : boundary_type_t::closed};
                                        }
                                    }();

                                    auto base_point_zeta =
                                        reduce_fraction(lower_bdy.linear_coeff * xi_min.xi +
                                                        lower_bdy.constant_coeff);

                                    // Check if the line actually degenerates into a point.
                                    {
                                        auto cmp_result = xi_min.xi <=> xi_max.xi;
                                        util::constexpr_assert(cmp_result <= 0);
                                        if (cmp_result == 0) {
                                            util::constexpr_assert(
                                                xi_min.boundary_type == boundary_type_t::closed &&
                                                xi_max.boundary_type == boundary_type_t::closed);

                                            return xi_zeta_region::single_point{xi_min.xi,
                                                                                base_point_zeta};
                                        }
                                    }

                                    auto direction_vector_xi =
                                        reduce_fraction(xi_max.xi - xi_min.xi);

                                    return xi_zeta_region::line_segment{
                                        .base_point_xi = std::move(xi_min.xi),
                                        .base_point_zeta = std::move(base_point_zeta),
                                        .direction_vector_xi =
                                            frac<bigint::uint_var, bigint::uint_var>{
                                                util::abs(direction_vector_xi.numerator),
                                                direction_vector_xi.denominator},
                                        .direction_vector_zeta = reduce_fraction(
                                            lower_bdy.linear_coeff * direction_vector_xi),
                                        .left_boundary_type = xi_min.boundary_type,
                                        .right_boundary_type = xi_max.boundary_type};
                                }

                                // The lower goes strictly below the upper.
                                // Start pushing slices.
                                found_min_xi = true;

                                util::constexpr_assert(
                                    lower_idx != lower_envelope.boundary_lines.size() - 1 ||
                                    upper_idx != upper_envelope.boundary_lines.size() - 1);

                                auto zeta = reduce_fraction(
                                    lower_envelope.boundary_lines[lower_idx].linear_coeff *
                                        intersection_xi +
                                    lower_envelope.boundary_lines[lower_idx].constant_coeff);
                                auto boundary_type =
                                    lower_envelope.boundary_lines[lower_idx].boundary_type ==
                                                boundary_type_t::closed &&
                                            upper_envelope.boundary_lines[upper_idx]
                                                    .boundary_type == boundary_type_t::closed
                                        ? boundary_type_t::closed
                                        : boundary_type_t::open;

                                if (position_in_lower_domain ==
                                    position_in_interval_t::right_endpoint) {
                                    if (lower_envelope.turning_point_types[lower_idx] ==
                                        boundary_type_t::open) {
                                        boundary_type = boundary_type_t::open;
                                    }
                                    ++lower_idx;
                                }
                                if (position_in_upper_domain ==
                                    position_in_interval_t::right_endpoint) {
                                    if (upper_envelope.turning_point_types[upper_idx] ==
                                        boundary_type_t::open) {
                                        boundary_type = boundary_type_t::open;
                                    }
                                    ++upper_idx;
                                }

                                turning_points.push_back(
                                    xi_zeta_region::bounded_polygon::turning_point_info{
                                        .xi = std::move(intersection_xi),
                                        .lower_zeta = zeta,
                                        .upper_zeta = zeta,
                                        .lower_boundary_type = boundary_type,
                                        .upper_boundary_type = boundary_type});

                                push_slice();
                                continue;
                            } // if (!found_min_xi)
                            // Otherwise, it must be the right-most point.
                            else {
                                util::constexpr_assert(!infinite_max_xi);
                                // Push the last slice and break.
                                boundary_line_pairs.push_back(
                                    xi_zeta_region::bounded_polygon::boundary_line_pair_info{
                                        .lower_boundary_linear_coeff = lower_bdy.linear_coeff,
                                        .lower_boundary_constant_coeff = lower_bdy.constant_coeff,
                                        .upper_boundary_linear_coeff = upper_bdy.linear_coeff,
                                        .upper_boundary_constant_coeff = upper_bdy.constant_coeff,
                                        .lower_boundary_type = lower_bdy.boundary_type,
                                        .upper_boundary_type = upper_bdy.boundary_type});

                                auto zeta = reduce_fraction(
                                    lower_envelope.boundary_lines[lower_idx].linear_coeff *
                                        intersection_xi +
                                    lower_envelope.boundary_lines[lower_idx].constant_coeff);
                                auto boundary_type =
                                    lower_envelope.boundary_lines[lower_idx].boundary_type ==
                                                boundary_type_t::closed &&
                                            upper_envelope.boundary_lines[upper_idx]
                                                    .boundary_type == boundary_type_t::closed
                                        ? boundary_type_t::closed
                                        : boundary_type_t::open;

                                if (position_in_lower_domain ==
                                        position_in_interval_t::right_endpoint &&
                                    lower_envelope.turning_point_types[lower_idx] ==
                                        boundary_type_t::open) {
                                    boundary_type = boundary_type_t::open;
                                }
                                if (position_in_upper_domain ==
                                        position_in_interval_t::right_endpoint &&
                                    upper_envelope.turning_point_types[upper_idx] ==
                                        boundary_type_t::open) {
                                    boundary_type = boundary_type_t::open;
                                }

                                turning_points.push_back(
                                    xi_zeta_region::bounded_polygon::turning_point_info{
                                        .xi = std::move(intersection_xi),
                                        .lower_zeta = zeta,
                                        .upper_zeta = zeta,
                                        .lower_boundary_type = boundary_type,
                                        .upper_boundary_type = boundary_type});

                                break;
                            }
                        } // If the intersection happens within the domains.
                    }     // If two lines are not parallel.
                    // If two lines coincide.
                    else if (lower_bdy.linear_coeff == upper_bdy.linear_coeff &&
                             lower_bdy.constant_coeff == upper_bdy.constant_coeff) {
                        // This is possible only for the first segments.
                        // Since the intersection must be bounded, the left-most point must be
                        // xi_min_itr->threshold.
                        util::constexpr_assert(lower_idx == 0 && upper_idx == 0);
                        util::constexpr_assert(!right_half_planes.empty() &&
                                               !left_half_planes.empty());

                        struct endpoint_t {
                            frac_t xi;
                            boundary_type_t boundary_type;
                        };

                        auto base_point_zeta =
                            reduce_fraction(lower_bdy.linear_coeff * xi_min_itr->threshold +
                                            lower_bdy.constant_coeff);

                        auto const xi_min_position_in_lower_domain =
                            compute_position_in_domain(xi_min_itr->threshold, lower_envelope, 0);
                        auto const xi_min_position_in_upper_domain =
                            compute_position_in_domain(xi_min_itr->threshold, upper_envelope, 0);
                        auto const xi_max_position_in_lower_domain =
                            compute_position_in_domain(xi_max_itr->threshold, lower_envelope, 0);
                        auto const xi_max_position_in_upper_domain =
                            compute_position_in_domain(xi_max_itr->threshold, upper_envelope, 0);

                        util::constexpr_assert(
                            xi_min_position_in_lower_domain != position_in_interval_t::exterior &&
                            xi_min_position_in_upper_domain != position_in_interval_t::exterior &&
                            xi_min_position_in_lower_domain !=
                                position_in_interval_t::left_endpoint &&
                            xi_min_position_in_upper_domain !=
                                position_in_interval_t::left_endpoint &&
                            xi_max_position_in_lower_domain !=
                                position_in_interval_t::left_endpoint &&
                            xi_max_position_in_upper_domain !=
                                position_in_interval_t::left_endpoint);

                        bool is_single_point = false;
                        if (xi_min_position_in_lower_domain ==
                            position_in_interval_t::right_endpoint) {
                            is_single_point = true;
                            util::constexpr_assert(lower_envelope.turning_point_types[0] !=
                                                   boundary_type_t::open);
                        }
                        if (xi_min_position_in_upper_domain ==
                            position_in_interval_t::right_endpoint) {
                            is_single_point = true;
                            util::constexpr_assert(upper_envelope.turning_point_types[0] !=
                                                   boundary_type_t::open);
                        }

                        if (is_single_point) {
                            return xi_zeta_region::single_point{std::move(xi_min_itr->threshold),
                                                                std::move(base_point_zeta)};
                        }

                        auto xi_max = [&]() -> endpoint_t {
                            auto const cmp_result = lower_envelope.turning_points[0] <=>
                                                    upper_envelope.turning_points[0];
                            if (cmp_result <= 0) {
                                if (xi_max_position_in_lower_domain ==
                                    position_in_interval_t::exterior) {
                                    return {lower_envelope.turning_points[0],
                                            lower_envelope.turning_point_types[0]};
                                }
                                else if (xi_max_position_in_lower_domain ==
                                         position_in_interval_t::interior) {
                                    return {std::move(xi_max_itr->threshold),
                                            xi_max_itr->boundary_type};
                                }
                                else {
                                    return {std::move(xi_max_itr->threshold),
                                            (lower_envelope.turning_point_types[0] ==
                                             boundary_type_t::closed) &&
                                                    (cmp_result < 0 ||
                                                     upper_envelope.turning_point_types[0] ==
                                                         boundary_type_t::closed)
                                                ? xi_max_itr->boundary_type
                                                : boundary_type_t::open};
                                }
                            }
                            else {
                                if (xi_max_position_in_upper_domain ==
                                    position_in_interval_t::exterior) {
                                    return {upper_envelope.turning_points[0],
                                            upper_envelope.turning_point_types[0]};
                                }
                                else if (xi_max_position_in_upper_domain ==
                                         position_in_interval_t::interior) {
                                    return {std::move(xi_max_itr->threshold),
                                            xi_max_itr->boundary_type};
                                }
                                else {
                                    return {std::move(xi_max_itr->threshold),
                                            upper_envelope.turning_point_types[0] ==
                                                    boundary_type_t::closed
                                                ? xi_max_itr->boundary_type
                                                : boundary_type_t::open};
                                }
                            }
                        }();

                        auto direction_vector_xi =
                            reduce_fraction(xi_max.xi - xi_min_itr->threshold);

                        return xi_zeta_region::line_segment{
                            .base_point_xi = std::move(xi_min_itr->threshold),
                            .base_point_zeta = std::move(base_point_zeta),
                            .direction_vector_xi =
                                frac<bigint::uint_var, bigint::uint_var>{
                                    util::abs(direction_vector_xi.numerator),
                                    direction_vector_xi.denominator},
                            .direction_vector_zeta =
                                reduce_fraction(lower_bdy.linear_coeff * direction_vector_xi),
                            .left_boundary_type = xi_min_itr->boundary_type,
                            .right_boundary_type = xi_max.boundary_type};
                    } // If two lines coincide.

                    // Found no intersection.
                    // Add a new slice if we already have found the left-most intersection point, or
                    // we know it is infinite.
                    if (found_min_xi) {
                        push_slice();

                        // If the right-most intersection point is infinite and both of the
                        // envelopes have reached the end, push the last slice and break.
                        if (infinite_max_xi &&
                            lower_idx == lower_envelope.boundary_lines.size() - 1 &&
                            upper_idx == upper_envelope.boundary_lines.size() - 1) {
                            boundary_line_pairs.push_back(
                                xi_zeta_region::bounded_polygon::boundary_line_pair_info{
                                    .lower_boundary_linear_coeff =
                                        lower_envelope.boundary_lines.back().linear_coeff,
                                    .lower_boundary_constant_coeff =
                                        lower_envelope.boundary_lines.back().constant_coeff,
                                    .upper_boundary_linear_coeff =
                                        upper_envelope.boundary_lines.back().linear_coeff,
                                    .upper_boundary_constant_coeff =
                                        upper_envelope.boundary_lines.back().constant_coeff,
                                    .lower_boundary_type =
                                        lower_envelope.boundary_lines.back().boundary_type,
                                    .upper_boundary_type =
                                        upper_envelope.boundary_lines.back().boundary_type});

                            break;
                        }
                    }
                    else {
                        // Choose the one with smaller upper bound and move to the next segment.
                        if (lower_idx < lower_envelope.boundary_lines.size() - 1) {
                            if (upper_idx < upper_envelope.boundary_lines.size() - 1) {
                                auto const cmp_result = lower_envelope.turning_points[lower_idx] <=>
                                                        upper_envelope.turning_points[upper_idx];
                                if (cmp_result <= 0) {
                                    ++lower_idx;
                                }
                                if (cmp_result >= 0) {
                                    ++upper_idx;
                                }
                            }
                            else {
                                ++lower_idx;
                            }
                        }
                        else {
                            ++upper_idx;
                        }
                    }
                } // while (true)
            }     // End of the computation of the intersection of lower/upper intersections.


            ////////////////////////////////////////////////////////////////////////////////////
            // Step 4 - Vertically cut the region if necessary.
            ////////////////////////////////////////////////////////////////////////////////////

            auto left_boundary_type = boundary_type_t::closed;
            auto right_boundary_type = boundary_type_t::closed;

            if (!right_half_planes.empty()) {
                util::constexpr_assert(!left_half_planes.empty());

                left_boundary_type = xi_min_itr->boundary_type;
                right_boundary_type = xi_max_itr->boundary_type;

                decltype(turning_points) trimmed_turning_points;
                decltype(boundary_line_pairs) trimmed_boundary_line_pairs;

                std::size_t idx = 0;
                for (; idx < turning_points.size(); ++idx) {
                    auto const cmp_result = xi_min_itr->threshold <=> turning_points[idx].xi;
                    if (cmp_result < 0) {
                        // Left boundary is outside of the region.
                        if (!infinite_min_xi && idx == 0) {
                            util::constexpr_assert(turning_points[0].lower_zeta ==
                                                   turning_points[0].upper_zeta);
                            left_boundary_type =
                                turning_points[0].lower_boundary_type == boundary_type_t::closed &&
                                        turning_points[0].upper_boundary_type ==
                                            boundary_type_t::closed
                                    ? boundary_type_t::closed
                                    : boundary_type_t::open;

                            trimmed_turning_points.push_back(turning_points[0]);
                            ++idx;
                            break;
                        }

                        auto const& boundary_line_pair = infinite_min_xi
                                                             ? boundary_line_pairs[idx]
                                                             : boundary_line_pairs[idx - 1];

                        auto lower_zeta = reduce_fraction(
                            boundary_line_pair.lower_boundary_linear_coeff * xi_min_itr->threshold +
                            boundary_line_pair.lower_boundary_constant_coeff);
                        auto upper_zeta = reduce_fraction(
                            boundary_line_pair.lower_boundary_linear_coeff * xi_min_itr->threshold +
                            boundary_line_pair.lower_boundary_constant_coeff);

                        util::constexpr_assert(lower_zeta < upper_zeta);

                        trimmed_turning_points.push_back(
                            xi_zeta_region::bounded_polygon::turning_point_info{
                                .xi = std::move(xi_min_itr->threshold),
                                .lower_zeta = std::move(lower_zeta),
                                .upper_zeta = std::move(upper_zeta),
                                .lower_boundary_type = boundary_line_pair.lower_boundary_type,
                                .upper_boundary_type = boundary_line_pair.upper_boundary_type});
                        break;
                    } // if (cmp_result < 0)
                    else if (cmp_result == 0) {
                        auto const& turning_point = turning_points[idx];

                        if (!infinite_max_xi && idx == turning_points.size() - 1) {
                            // Left boundary is at the right-end of the region, so the region is
                            // actually a point.
                            util::constexpr_assert(turning_point.lower_zeta ==
                                                   turning_point.upper_zeta);
                            util::constexpr_assert(
                                left_boundary_type == boundary_type_t::closed &&
                                turning_point.lower_boundary_type == boundary_type_t::closed &&
                                turning_point.upper_boundary_type == boundary_type_t::closed);

                            return xi_zeta_region::single_point{
                                std::move(turning_point.xi), std::move(turning_point.lower_zeta)};
                        }

                        if (turning_point.lower_zeta == turning_point.upper_zeta &&
                            (turning_point.lower_boundary_type == boundary_type_t::open ||
                             turning_point.upper_boundary_type == boundary_type_t::open)) {
                            left_boundary_type = boundary_type_t::open;
                        }

                        trimmed_turning_points.push_back(turning_point);
                        ++idx;
                        break;
                    } // else if (cmp_result == 0)
                }     // for (; idx < turning_points.size(); ++idx)

                while (true) {
                    auto& turning_point = turning_points[idx];
                    auto& boundary_line_pair =
                        infinite_min_xi ? boundary_line_pairs[idx] : boundary_line_pairs[idx - 1];

                    trimmed_boundary_line_pairs.push_back(boundary_line_pair);

                    auto const cmp_result = xi_max_itr->threshold <=> turning_point.xi;
                    if (cmp_result < 0) {
                        auto lower_zeta = reduce_fraction(
                            boundary_line_pair.lower_boundary_linear_coeff * xi_max_itr->threshold +
                            boundary_line_pair.lower_boundary_constant_coeff);
                        auto upper_zeta = reduce_fraction(
                            boundary_line_pair.lower_boundary_linear_coeff * xi_max_itr->threshold +
                            boundary_line_pair.lower_boundary_constant_coeff);

                        util::constexpr_assert(lower_zeta < upper_zeta);

                        trimmed_turning_points.push_back(
                            xi_zeta_region::bounded_polygon::turning_point_info{
                                .xi = std::move(xi_max_itr->threshold),
                                .lower_zeta = std::move(lower_zeta),
                                .upper_zeta = std::move(upper_zeta),
                                .lower_boundary_type = boundary_line_pair.lower_boundary_type,
                                .upper_boundary_type = boundary_line_pair.upper_boundary_type});
                        break;
                    } // if (cmp_result < 0)
                    else if (cmp_result == 0) {
                        if (!infinite_min_xi && idx == 0) {
                            // Right boundary is at the left-end of the region, so the region is
                            // actually a point.
                            util::constexpr_assert(turning_point.lower_zeta ==
                                                   turning_point.upper_zeta);
                            util::constexpr_assert(
                                right_boundary_type == boundary_type_t::closed &&
                                turning_point.lower_boundary_type == boundary_type_t::closed &&
                                turning_point.upper_boundary_type == boundary_type_t::closed);

                            return xi_zeta_region::single_point{
                                std::move(turning_point.xi), std::move(turning_point.lower_zeta)};
                        }

                        if (turning_point.lower_zeta == turning_point.upper_zeta &&
                            (turning_point.lower_boundary_type == boundary_type_t::open ||
                             turning_point.upper_boundary_type == boundary_type_t::open)) {
                            right_boundary_type = boundary_type_t::open;
                        }

                        trimmed_turning_points.push_back(std::move(turning_point));
                        break;
                    } // else if (cmp_result == 0)

                    trimmed_turning_points.push_back(std::move(turning_point));

                    // If xi_max_itr->threshold is strictly larger than the last turning point.
                    if (++idx == turning_points.size()) {
                        // If the region was unbounded to right, cut it.
                        // Otherwise, there is nothing further to do.
                        if (infinite_max_xi) {
                            trimmed_boundary_line_pairs.push_back(boundary_line_pairs.back());

                            auto lower_zeta = reduce_fraction(
                                boundary_line_pairs.back().lower_boundary_linear_coeff *
                                    xi_max_itr->threshold +
                                boundary_line_pairs.back().lower_boundary_constant_coeff);
                            auto upper_zeta = reduce_fraction(
                                boundary_line_pairs.back().lower_boundary_linear_coeff *
                                    xi_max_itr->threshold +
                                boundary_line_pairs.back().lower_boundary_constant_coeff);

                            util::constexpr_assert(lower_zeta < upper_zeta);

                            trimmed_turning_points.push_back(
                                xi_zeta_region::bounded_polygon::turning_point_info{
                                    .xi = std::move(xi_max_itr->threshold),
                                    .lower_zeta = std::move(lower_zeta),
                                    .upper_zeta = std::move(upper_zeta),
                                    .lower_boundary_type =
                                        boundary_line_pairs.back().lower_boundary_type,
                                    .upper_boundary_type =
                                        boundary_line_pairs.back().upper_boundary_type});
                        }

                        break;
                    }
                } // while(true)
                turning_points = std::move(trimmed_turning_points);
                boundary_line_pairs = std::move(trimmed_boundary_line_pairs);
            } // End of the vertical cutting.
            // If there is no vertical cutting, the left/right ends must be sharp.
            else {
                util::constexpr_assert(!infinite_min_xi && !infinite_max_xi);
                util::constexpr_assert(turning_points.front().lower_zeta ==
                                       turning_points.front().upper_zeta);
                util::constexpr_assert(turning_points.back().lower_zeta ==
                                       turning_points.back().upper_zeta);

                if (turning_points.front().lower_boundary_type == boundary_type_t::open ||
                    turning_points.front().upper_boundary_type == boundary_type_t::open) {
                    left_boundary_type = boundary_type_t::open;
                }
                if (turning_points.back().lower_boundary_type == boundary_type_t::open ||
                    turning_points.back().upper_boundary_type == boundary_type_t::open) {
                    right_boundary_type = boundary_type_t::open;
                }
            }

            util::constexpr_assert(turning_points.size() == boundary_line_pairs.size() + 1);
            return xi_zeta_region::bounded_polygon{std::move(turning_points),
                                                   std::move(boundary_line_pairs),
                                                   left_boundary_type, right_boundary_type};
        }
    }
}

#endif
