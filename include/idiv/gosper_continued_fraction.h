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

#ifndef JKJ_HEADER_GOSPER_CONTINUED_FRACTION
#define JKJ_HEADER_GOSPER_CONTINUED_FRACTION

#include "continued_fraction.h"
#include <cstdlib>

namespace jkj {
    namespace cntfrc {
        namespace impl {
            namespace detail {
                template <class Int>
                struct gosper_util {
                    // Find the kernel of a rank-1 linear fractional transform.
                    template <class NumNum, class DenNum, class NumDen, class DenDen>
                    static constexpr projective_rational<Int, Int> kernel_of_rank1_transform(
                        linear_fractional_transform<NumNum, DenNum, NumDen, DenDen> const&
                            transform) {
                        return util::is_zero(transform.num_to_num) &&
                                       util::is_zero(transform.den_to_num)
                                   ? projective_rational{-Int{transform.den_to_den},
                                                         Int{transform.num_to_den}}
                                   : projective_rational{-Int{transform.den_to_num},
                                                         Int{transform.num_to_num}};
                    }
                    // Find the range of a rank-1 linear fractional transform.
                    template <class NumNum, class DenNum, class NumDen, class DenDen>
                    static constexpr projective_rational<Int, Int>
                    range_of_rank1_transform(linear_fractional_transform<NumNum, DenNum, NumDen,
                                                                         DenDen> const& transform) {
                        return util::is_zero(transform.num_to_num) &&
                                       util::is_zero(transform.num_to_den)
                                   ? projective_rational{Int{transform.den_to_num},
                                                         Int{transform.den_to_den}}
                                   : projective_rational{Int{transform.num_to_num},
                                                         Int{transform.num_to_den}};
                    }
                    template <class CyclicInterval, class T>
                    static constexpr bool contains_in_closure(CyclicInterval const& itv,
                                                              T const& x) {
                        return itv.visit([&x](auto&& itv) {
                            using enum cyclic_interval_type_t;
                            using itv_type = std::remove_cvref_t<decltype(itv)>;
                            static_assert(itv_type::interval_type() != empty);

                            if constexpr (itv_type::interval_type() == entire) {
                                return true;
                            }
                            else if constexpr (itv_type::interval_type() == single_point) {
                                return itv.lower_bound() == x;
                            }
                            else {
                                return itv.lower_bound() == x || itv.upper_bound() == x ||
                                       cyclic_order(itv.lower_bound(), x, itv.upper_bound());
                            }
                        });
                    }
                };
            }

            template <class InternalContinuedFractionImpl, class Unity = unity>
            class unary_gosper {
            public:
                using int_type =
                    decltype(InternalContinuedFractionImpl::convergent_type::numerator);
                using uint_type =
                    decltype(InternalContinuedFractionImpl::convergent_type::denominator);
                using partial_fraction_type = projective_rational<Unity, int_type>;
                using convergent_type = typename InternalContinuedFractionImpl::convergent_type;
                using interval_type = variable_shape_cyclic_interval<
                    projective_rational<int_type, int_type>, cyclic_interval_type_t::single_point,
                    cyclic_interval_type_t::left_open_right_closed,
                    cyclic_interval_type_t::left_closed_right_open, cyclic_interval_type_t::entire>;
                using internal_continued_fraction_impl_type = InternalContinuedFractionImpl;

            private:
                generator<InternalContinuedFractionImpl, interval_tracker> cf_;
                linear_fractional_transform<int_type> coeff_;
                int determinant_sign_ = 0;
                bool is_first_ = true;

            public:
                static constexpr partial_fraction_type initial_partial_fraction() {
                    return {Unity{}, int_type{0}};
                }
                static constexpr interval_type initial_interval() {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }

                constexpr unary_gosper(internal_continued_fraction_impl_type cf_impl,
                                       linear_fractional_transform<int_type> coeff)
                    : cf_{std::move(cf_impl)},
                      coeff_{static_cast<linear_fractional_transform<int_type>&&>(coeff)} {
                    determinant_sign_ = coeff_.determinant_sign();
                    // Step 1. Get away from singularities.
                    if (determinant_sign_ == 0) {
                        // The degenerate case a = b = c = d = 0 is disallowed.
                        util::constexpr_assert(!util::is_zero(coeff_.num_to_num) ||
                                               !util::is_zero(coeff_.den_to_num));

                        auto singularity =
                            detail::gosper_util<int_type>::kernel_of_rank1_transform(coeff_);

                        while (true) {
                            // If the unique singularity is contained in the closure of the current
                            // interval, then refine the interval.
                            if (detail::gosper_util<int_type>::contains_in_closure(
                                    cf_.current_interval(), singularity)) {
                                // We do not allow the input to be at the singularity.
                                util::constexpr_assert(cf_.current_interval().interval_type() !=
                                                       cyclic_interval_type_t::single_point);
                                continue;
                            }
                            break;
                        }
                    }
                }

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    // Step 3. Compare the floor of two endpoints and if they are equal, return.
                    auto update_and_callback = [&](auto&& common_floor) {
                        is_first_ = false;
                        coeff_.translate(-common_floor);
                        coeff_.reflect();
                        determinant_sign_ *= -1;
                        callback(partial_fraction_type{
                            Unity{}, static_cast<decltype(common_floor)&&>(common_floor)});
                    };
                    enum class final_result { success, terminate, fail };
                    auto check_floor = [&](auto&& itv) -> final_result {
                        using enum cyclic_interval_type_t;
                        using itv_type = std::remove_cvref_t<decltype(itv)>;
                        static_assert(itv_type::interval_type() != empty &&
                                      itv_type::interval_type() != entire);
                        constexpr auto infinity = projective_rational{unity{}, zero{}};

                        if constexpr (itv_type::interval_type() == single_point) {
                            // If infinity is the only element in the range, then there is no
                            // further continued fraction coefficients.
                            if (itv.lower_bound() == infinity) {
                                return final_result::terminate;
                            }
                            else {
                                // If a finite rational number is the only element in the range,
                                // then just compute the continued fraction expansion of that
                                // number.
                                update_and_callback(util::div_floor(itv.lower_bound().numerator,
                                                                    itv.lower_bound().denominator));
                                return final_result::success;
                            }
                        }
                        else {
                            auto lower_bound = itv.lower_bound();
                            bool lower_bound_inclusive =
                                (itv_type::left_endpoint_type() == endpoint_type_t::closed);
                            bool upper_bound_inclusive =
                                (itv_type::right_endpoint_type() == endpoint_type_t::closed);
                            if (is_first_) {
                                // Cannot do anything if the range estimate contains the infinity.
                                if (lower_bound == infinity || itv.upper_bound() == infinity ||
                                    cyclic_order(lower_bound, infinity, itv.upper_bound())) {
                                    return final_result::fail;
                                }
                            }
                            else {
                                // Except for the first coefficient, the value must be in (1,infty].
                                // Still, we cannot do anything if infinity is strictly in the
                                // interval or the upper bound is equal to infinity.
                                if (itv.upper_bound() == infinity ||
                                    cyclic_order(lower_bound, infinity, itv.upper_bound())) {
                                    return final_result::fail;
                                }
                                // If lower bound is infinity and is inclusive, then the value might
                                // be infinity.
                                if (lower_bound == infinity) {
                                    if (lower_bound_inclusive &&
                                        (itv.upper_bound() == projective_rational<unity, unity>{} ||
                                         cyclic_order(lower_bound, itv.upper_bound(),
                                                      projective_rational<unity, unity>{}))) {
                                        return final_result::terminate;
                                    }
                                    return final_result::fail;
                                }

                                util::constexpr_assert(itv.upper_bound() !=
                                                       projective_rational<unity, unity>{});
                                if (cyclic_order(lower_bound, projective_rational<unity, unity>{},
                                                 itv.upper_bound())) {
                                    lower_bound.numerator = 1;
                                    lower_bound.denominator = 1;
                                }
                            }

                            // Compare the floor.
                            auto floor_lower = util::div_floor(std::move(lower_bound.numerator),
                                                               std::move(lower_bound.denominator));
                            auto floor_upper =
                                upper_bound_inclusive
                                    ? util::div_floor(itv.upper_bound().numerator,
                                                      itv.upper_bound().denominator)
                                    : (util::div_ceil(itv.upper_bound().numerator,
                                                      itv.upper_bound().denominator) -
                                       1);

                            if (floor_lower == floor_upper) {
                                update_and_callback(
                                    static_cast<decltype(floor_lower)&&>(floor_lower));
                                return final_result::success;
                            }
                            return final_result::fail;
                        }
                    };

                    while (true) {
                        // Step 2. Find the range of the linear fractional transform.
                        auto const result = cf_.current_interval().visit([&](auto&& itv)
                                                                             -> final_result {
                            using enum cyclic_interval_type_t;
                            using itv_type = std::remove_cvref_t<decltype(itv)>;
                            static_assert(itv_type::interval_type() != empty);

                            if constexpr (itv_type::interval_type() == entire) {
                                // Range is the entire RP1, so we cannot proceed.
                                util::constexpr_assert(determinant_sign_ != 0);
                                return final_result::fail;
                            }
                            else if constexpr (itv_type::interval_type() == single_point) {
                                return check_floor(
                                    cyclic_interval<projective_rational<int_type, int_type>,
                                                    single_point>{coeff_(itv.lower_bound())});
                            }
                            else {
                                if (determinant_sign_ == 0) {
                                    return check_floor(
                                        cyclic_interval<projective_rational<int_type, int_type>,
                                                        single_point>{coeff_(itv.lower_bound())});
                                }

                                auto first_end = coeff_(itv.lower_bound());
                                auto second_end = coeff_(itv.upper_bound());

                                if constexpr (itv_type::interval_type() == open ||
                                              itv_type::interval_type() == closed) {
                                    if (determinant_sign_ < 0) {
                                        using std::swap;
                                        swap(first_end, second_end);
                                    }
                                    return check_floor(
                                        cyclic_interval<projective_rational<int_type, int_type>,
                                                        itv_type::interval_type()>{first_end,
                                                                                   second_end});
                                }
                                else if constexpr (itv_type::interval_type() ==
                                                   left_open_right_closed) {
                                    if (determinant_sign_ < 0) {
                                        return check_floor(
                                            cyclic_interval<projective_rational<int_type, int_type>,
                                                            left_closed_right_open>{second_end,
                                                                                    first_end});
                                    }
                                    else {
                                        return check_floor(
                                            cyclic_interval<projective_rational<int_type, int_type>,
                                                            left_open_right_closed>{first_end,
                                                                                    second_end});
                                    }
                                }
                                else {
                                    if (determinant_sign_ < 0) {
                                        return check_floor(
                                            cyclic_interval<projective_rational<int_type, int_type>,
                                                            left_open_right_closed>{second_end,
                                                                                    first_end});
                                    }
                                    else {
                                        return check_floor(
                                            cyclic_interval<projective_rational<int_type, int_type>,
                                                            left_closed_right_open>{first_end,
                                                                                    second_end});
                                    }
                                }
                            }
                        });

                        switch (result) {
                        case final_result::success:
                        case final_result::terminate:
                            return;
                        case final_result::fail:;
                        }
                        cf_.update();
                    }
                }
            };

            template <class ContinuedFractionImplX, class ContinuedFractionImplY,
                      class Unity = unity>
                requires std::is_same_v<typename ContinuedFractionImplX::convergent_type,
                                        typename ContinuedFractionImplY::convergent_type>
            class binary_gosper {
            public:
                using int_type = decltype(ContinuedFractionImplX::convergent_type::numerator);
                using uint_type = decltype(ContinuedFractionImplX::convergent_type::denominator);
                using partial_fraction_type = projective_rational<Unity, int_type>;
                using convergent_type = typename ContinuedFractionImplX::convergent_type;
                using interval_type = variable_shape_cyclic_interval<
                    convergent_type, cyclic_interval_type_t::single_point,
                    cyclic_interval_type_t::left_open_right_closed,
                    cyclic_interval_type_t::left_closed_right_open, cyclic_interval_type_t::entire>;
                using first_internal_continued_fraction_impl_type = ContinuedFractionImplX;
                using second_internal_continued_fraction_impl_type = ContinuedFractionImplY;

            private:
                generator<ContinuedFractionImplX, interval_tracker> xcf_;
                generator<ContinuedFractionImplY, interval_tracker> ycf_;
                bilinear_fractional_transform<int_type> coeff_;
                bool is_first_ = true;

                // Precondition: itv is contained in the domain of transform.
                template <class NumNum, class DenNum, class NumDen, class DenDen,
                          class InputIntervalType>
                static constexpr variable_shape_cyclic_interval<
                    projective_rational<int_type, int_type>, cyclic_interval_type_t::single_point,
                    cyclic_interval_type_t::closed, cyclic_interval_type_t::entire>
                map_cyclic_interval(
                    linear_fractional_transform<NumNum, DenNum, NumDen, DenDen> const& transform,
                    InputIntervalType const& itv, int determinant_sign) {
                    using return_type =
                        variable_shape_cyclic_interval<projective_rational<int_type, int_type>,
                                                       cyclic_interval_type_t::single_point,
                                                       cyclic_interval_type_t::closed,
                                                       cyclic_interval_type_t::entire>;

                    return itv.visit([&transform, determinant_sign](auto&& itv) -> return_type {
                        using enum cyclic_interval_type_t;
                        using value_type = projective_rational<int_type, int_type>;
                        using itv_type = std::remove_cvref_t<decltype(itv)>;
                        static_assert(itv_type::interval_type() != empty);

                        if constexpr (itv_type::interval_type() == entire) {
                            util::constexpr_assert(determinant_sign != 0);
                            return cyclic_interval<value_type, entire>{};
                        }
                        else if constexpr (itv_type::interval_type() == single_point) {
                            return cyclic_interval<value_type, single_point>{
                                transform(itv.lower_bound())};
                        }
                        else {
                            if (determinant_sign != 0) {
                                auto first_end = transform(itv.lower_bound());
                                auto second_end = transform(itv.upper_bound());
                                if (determinant_sign < 0) {
                                    using std::swap;
                                    swap(first_end, second_end);
                                }

                                // Always return closed interval for the sake of simplicity.
                                return cyclic_interval<value_type, closed>{first_end, second_end};
                            }
                            else {
                                return cyclic_interval<value_type, single_point>{
                                    detail::gosper_util<int_type>::range_of_rank1_transform(
                                        transform)};
                            }
                        }
                    });
                }

                // Check a predicate on each component. Assumes itv1, itv2 are both closed.
                template <class InputIntervalType1, class InputIntervalType2, class Pred>
                static constexpr bool check_intersection(InputIntervalType1 const& itv1,
                                                         InputIntervalType2 const& itv2, Pred&& f) {
                    return itv1.visit([&itv2, &f](auto&& itv1) {
                        using enum cyclic_interval_type_t;
                        using value_type = projective_rational<int_type, int_type>;
                        using empty_interval = cyclic_interval<value_type, empty>;
                        using single_point_interval = cyclic_interval<value_type, single_point>;
                        using closed_interval = cyclic_interval<value_type, closed>;

                        constexpr auto itv1_type = itv1.interval_type();

                        if constexpr (itv1_type == empty) {
                            return f(empty_interval{});
                        }
                        else if constexpr (itv1_type == entire) {
                            return itv2.visit(f);
                        }
                        else if constexpr (itv1_type == single_point) {
                            if (detail::gosper_util<int_type>::contains_in_closure(
                                    itv2, itv1.lower_bound())) {
                                return f(single_point_interval{itv1.lower_bound()});
                            }
                            else {
                                return f(empty_interval{});
                            }
                        }
                        else {
                            return itv2.visit([&itv1, &f](auto&& itv2) {
                                constexpr auto itv2_type = itv2.interval_type();

                                if constexpr (itv2_type == empty) {
                                    return f(empty_interval{});
                                }
                                else if constexpr (itv2_type == entire) {
                                    return f(itv1);
                                }
                                else if constexpr (itv2_type == single_point) {
                                    if (detail::gosper_util<int_type>::contains_in_closure(
                                            itv1, itv2.lower_bound())) {
                                        return f(single_point_interval{itv2.lower_bound()});
                                    }
                                    else {
                                        return f(empty_interval{});
                                    }
                                }
                                else {
                                    // 16 different cases.
                                    auto&& a = itv1.lower_bound();
                                    auto&& b = itv1.upper_bound();
                                    auto&& c = itv2.lower_bound();
                                    auto&& d = itv2.upper_bound();
                                    if (a == c) {
                                        if (b == d) {
                                            return f(closed_interval{a, b});
                                        }
                                        else if (cyclic_order(a, b, d)) {
                                            return f(closed_interval{a, b});
                                        }
                                        else {
                                            // [a,d,b]
                                            return f(closed_interval{a, d});
                                        }
                                    }
                                    else if (a == d) {
                                        if (b == c) {
                                            return f(single_point_interval{a}) &&
                                                   f(single_point_interval{b});
                                        }
                                        else if (cyclic_order(a, b, c)) {
                                            return f(single_point_interval{a});
                                        }
                                        else {
                                            // [a,c,b]
                                            return f(single_point_interval{a}) &&
                                                   f(closed_interval{c, b});
                                        }
                                    }
                                    else if (b == c) {
                                        if (cyclic_order(a, b, d)) {
                                            return f(single_point_interval{b});
                                        }
                                        else {
                                            // [a,d,b]
                                            return f(closed_interval{a, d}) &&
                                                   f(single_point_interval{b});
                                        }
                                    }
                                    else if (b == d) {
                                        if (cyclic_order(a, b, c)) {
                                            return f(closed_interval{a, b});
                                        }
                                        else {
                                            // [a,c,b]
                                            return f(closed_interval{c, b});
                                        }
                                    }
                                    else if (cyclic_order(a, b, c)) {
                                        if (cyclic_order(c, d, a)) {
                                            // [a,b,c,d]
                                            return f(cyclic_interval<value_type, empty>{});
                                        }
                                        else if (cyclic_order(a, b, d)) {
                                            // [a,b,d,c]
                                            return f(closed_interval{a, b});
                                        }
                                        else {
                                            // [a,d,b,c]
                                            return f(closed_interval{a, d});
                                        }
                                    }
                                    else {
                                        if (cyclic_order(b, d, a)) {
                                            // [a,c,b,d]
                                            return f(closed_interval{c, b});
                                        }
                                        else if (cyclic_order(a, c, d)) {
                                            // [a,c,d,b]
                                            return f(closed_interval{c, d});
                                        }
                                        else {
                                            // [a,d,c,b]
                                            return f(closed_interval{a, d}) &&
                                                   f(closed_interval{c, b});
                                        }
                                    }
                                }
                            });
                        }
                    });
                }

            public:
                static constexpr auto initial_partial_fraction() {
                    return partial_fraction_type{Unity{}, int_type{0}};
                }
                static constexpr auto initial_interval() {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }

                constexpr binary_gosper(first_internal_continued_fraction_impl_type xcf_impl,
                                        second_internal_continued_fraction_impl_type ycf_impl,
                                        bilinear_fractional_transform<int_type> coeff)
                    : xcf_{std::move(xcf_impl)}, ycf_{std::move(ycf_impl)},
                      coeff_{static_cast<bilinear_fractional_transform<int_type>&&>(coeff)} {
                    // Step 1. Get away from singularities.
                    int det_sign_num_bilinear_form = 0;
                    int rank_num_bilinear_form = 0;
                    int det_sign_den_bilinear_form = 0;
                    int rank_den_bilinear_form = 0;

                    auto calculate_det_sign_rank = [](auto&& a, auto&& b, auto&& c, auto&& d,
                                                      auto& det_sign, auto& rank) {
                        auto const result = a * d <=> b * c;

                        if (result < 0) {
                            det_sign = -1;
                            rank = 2;
                        }
                        if (result > 0) {
                            det_sign = 1;
                            rank = 2;
                        }
                        else {
                            det_sign = 0;
                            rank = (util::is_zero(a) && util::is_zero(b) && util::is_zero(c) &&
                                            util::is_zero(d)
                                        ? 0
                                        : 1);
                        }
                    };

                    calculate_det_sign_rank(coeff_.xnum_ynum_to_num, coeff_.xnum_yden_to_num,
                                            coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num,
                                            det_sign_num_bilinear_form, rank_num_bilinear_form);

                    calculate_det_sign_rank(coeff_.xnum_ynum_to_den, coeff_.xnum_yden_to_den,
                                            coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den,
                                            det_sign_den_bilinear_form, rank_den_bilinear_form);

                    // Degenerate case A = B = 0 is disallowed.
                    util::constexpr_assert(rank_num_bilinear_form != 0 ||
                                           rank_den_bilinear_form != 0);

                    // L = A^T RB - B^T RA = 2sym(ag-ce  ah-cf  \\ bg-de  bh-df).
                    struct sing_det_form_t {
                        int_type a, b, d;
                    } sing_det_form = {((coeff_.xnum_ynum_to_num * coeff_.xden_ynum_to_den -
                                         coeff_.xden_ynum_to_num * coeff_.xnum_ynum_to_den)
                                        << 1),
                                       coeff_.xnum_ynum_to_num * coeff_.xden_yden_to_den +
                                           coeff_.xnum_yden_to_num * coeff_.xden_ynum_to_den -
                                           coeff_.xden_ynum_to_num * coeff_.xnum_yden_to_den -
                                           coeff_.xden_yden_to_num * coeff_.xnum_ynum_to_den,
                                       ((coeff_.xnum_yden_to_num * coeff_.xden_yden_to_den -
                                         coeff_.xden_yden_to_num * coeff_.xnum_yden_to_den)
                                        << 1)};

                    bool globally_well_defined =
                        ((sing_det_form.a * sing_det_form.d) > (sing_det_form.b * sing_det_form.b));

                    if (globally_well_defined) {
                        return;
                    }

                    // Check if there exists [x] in the closure of itv such that Ax and Bx are
                    // parallel.
                    auto has_parallel_output = [globally_well_defined,
                                                &sing_det_form](auto&& itv) -> bool {
                        using enum cyclic_interval_type_t;
                        using itv_type = std::remove_cvref_t<decltype(itv)>;

                        if constexpr (itv_type::interval_type() == empty) {
                            return false;
                        }
                        else if constexpr (itv_type::interval_type() == entire) {
                            return !globally_well_defined;
                        }
                        else if constexpr (itv_type::interval_type() == single_point) {
                            return util::is_zero(sing_det_form.a * itv.lower_bound().numerator *
                                                     itv.lower_bound().numerator +
                                                 ((sing_det_form.b * itv.lower_bound().numerator *
                                                   itv.lower_bound().denominator)
                                                  << 1) +
                                                 sing_det_form.d * itv.lower_bound().denominator *
                                                     itv.lower_bound().denominator);
                        }
                        else {
                            auto Luu = sing_det_form.a * itv.lower_bound().numerator *
                                           itv.lower_bound().numerator +
                                       ((sing_det_form.b * itv.lower_bound().numerator *
                                         itv.lower_bound().denominator)
                                        << 1) +
                                       sing_det_form.d * itv.lower_bound().denominator *
                                           itv.lower_bound().denominator;
                            auto Lvv = sing_det_form.a * itv.upper_bound().numerator *
                                           itv.upper_bound().numerator +
                                       ((sing_det_form.b * itv.upper_bound().numerator *
                                         itv.upper_bound().denominator)
                                        << 1) +
                                       sing_det_form.d * itv.upper_bound().denominator *
                                           itv.upper_bound().denominator;

                            int sign_Luu = util::is_zero(Luu)          ? 0
                                           : util::is_nonnegative(Luu) ? 1
                                                                       : -1;
                            int sign_Lvv = util::is_zero(Lvv)          ? 0
                                           : util::is_nonnegative(Lvv) ? 1
                                                                       : -1;
                            if (sign_Luu * sign_Lvv <= 0) {
                                return true;
                            }

                            auto Luv = sing_det_form.a * itv.lower_bound().numerator *
                                           itv.upper_bound().numerator +
                                       sing_det_form.b * itv.lower_bound().numerator *
                                           itv.upper_bound().denominator +
                                       sing_det_form.b * itv.lower_bound().denominator *
                                           itv.upper_bound().numerator +
                                       sing_det_form.d * itv.lower_bound().denominator *
                                           itv.upper_bound().denominator;
                            int sign_Luv = util::is_zero(Luv)          ? 0
                                           : util::is_nonnegative(Luv) ? 1
                                                                       : -1;
                            return sign_Luu * sign_Luv < 0;
                        }
                    };

                    auto has_singularity = [&] {
                        // Compute the set S = J \cap (T_RA^-1[I] \cup (dom(T_A))^c)
                        // \cap (T_RB^-1[I] \cup (dom(T_B))^c).

                        // A is invertible.
                        if (rank_num_bilinear_form == 2) {
                            // Can ignore the third term from the intersection.
                            // RA^T = (b d \\ -a -c)
                            return check_intersection(
                                ycf_.current_interval(),
                                map_cyclic_interval(
                                    linear_fractional_transform{
                                        coeff_.xnum_yden_to_num, coeff_.xden_yden_to_num,
                                        -coeff_.xnum_ynum_to_num, -coeff_.xden_ynum_to_num},
                                    xcf_.current_interval(), det_sign_num_bilinear_form),
                                has_parallel_output);
                        }
                        // B is invertible.
                        else if (rank_den_bilinear_form == 2) {
                            // Can ignore the second term from the intersection.
                            // RB^T = (f h \\ -e -g)
                            return check_intersection(
                                ycf_.current_interval(),
                                map_cyclic_interval(
                                    linear_fractional_transform{
                                        coeff_.xnum_yden_to_den, coeff_.xden_yden_to_den,
                                        -coeff_.xnum_ynum_to_den, -coeff_.xden_ynum_to_den},
                                    xcf_.current_interval(), det_sign_den_bilinear_form),
                                has_parallel_output);
                        }
                        // A is zero.
                        else if (rank_num_bilinear_form == 0) {
                            util::constexpr_assert(rank_den_bilinear_form == 1);
                            // <u, Bv> = 0 if and only if v in ker B or u is in RB[R^2].
                            return detail::gosper_util<int_type>::contains_in_closure(
                                       xcf_.current_interval(),
                                       detail::gosper_util<int_type>::range_of_rank1_transform(
                                           linear_fractional_transform{
                                               coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den,
                                               -coeff_.xnum_ynum_to_den,
                                               -coeff_.xnum_yden_to_den})) ||
                                   detail::gosper_util<int_type>::contains_in_closure(
                                       ycf_.current_interval(),
                                       detail::gosper_util<int_type>::kernel_of_rank1_transform(
                                           linear_fractional_transform{
                                               coeff_.xnum_ynum_to_den, coeff_.xnum_yden_to_den,
                                               coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den}));
                        }
                        // B is zero.
                        else if (rank_den_bilinear_form == 0) {
                            util::constexpr_assert(rank_num_bilinear_form == 1);
                            // <u, Av> = 0 if and only if v in ker A or u is in RA[R^2].
                            return detail::gosper_util<int_type>::contains_in_closure(
                                       xcf_.current_interval(),
                                       detail::gosper_util<int_type>::range_of_rank1_transform(
                                           linear_fractional_transform{
                                               coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num,
                                               -coeff_.xnum_ynum_to_num,
                                               -coeff_.xnum_yden_to_num})) ||
                                   detail::gosper_util<int_type>::contains_in_closure(
                                       ycf_.current_interval(),
                                       detail::gosper_util<int_type>::kernel_of_rank1_transform(
                                           linear_fractional_transform{
                                               coeff_.xnum_ynum_to_num, coeff_.xnum_yden_to_num,
                                               coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num}));
                        }
                        // Both are of rank 1.
                        else {
                            util::constexpr_assert(rank_num_bilinear_form == 1 &&
                                                   rank_den_bilinear_form == 1);

                            // RA = (c d \\ -a -b)
                            if (detail::gosper_util<int_type>::contains_in_closure(
                                    xcf_.current_interval(),
                                    detail::gosper_util<int_type>::range_of_rank1_transform(
                                        linear_fractional_transform{
                                            coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num,
                                            -coeff_.xnum_ynum_to_num, -coeff_.xnum_yden_to_num}))) {
                                // The second interval is the entire RP1.
                                // RB = (g h \\ -e -f)
                                if (detail::gosper_util<int_type>::contains_in_closure(
                                        xcf_.current_interval(),
                                        detail::gosper_util<int_type>::range_of_rank1_transform(
                                            linear_fractional_transform{
                                                coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den,
                                                -coeff_.xnum_ynum_to_den,
                                                -coeff_.xnum_yden_to_den}))) {
                                    // The third interval is also the entire RP1.
                                    return ycf_.current_interval().visit(has_parallel_output);
                                }
                                else {
                                    // The third interval is a single point.
                                    return check_intersection(
                                        ycf_.current_interval(),
                                        cyclic_interval<projective_rational<int_type, int_type>,
                                                        cyclic_interval_type_t::single_point>{
                                            detail::gosper_util<int_type>::
                                                kernel_of_rank1_transform(
                                                    linear_fractional_transform{
                                                        coeff_.xnum_ynum_to_den,
                                                        coeff_.xnum_yden_to_den,
                                                        coeff_.xden_ynum_to_den,
                                                        coeff_.xden_yden_to_den})},
                                        has_parallel_output);
                                }
                            }
                            else {
                                // The second interval is a single point.
                                // RB = (g h \\ -e -f)
                                if (detail::gosper_util<int_type>::contains_in_closure(
                                        xcf_.current_interval(),
                                        detail::gosper_util<int_type>::range_of_rank1_transform(
                                            linear_fractional_transform{
                                                coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den,
                                                -coeff_.xnum_ynum_to_den,
                                                -coeff_.xnum_yden_to_den}))) {
                                    // The third interval is the entire RP1.
                                    return check_intersection(
                                        ycf_.current_interval(),
                                        cyclic_interval<projective_rational<int_type, int_type>,
                                                        cyclic_interval_type_t::single_point>{
                                            detail::gosper_util<int_type>::
                                                kernel_of_rank1_transform(
                                                    linear_fractional_transform{
                                                        coeff_.xnum_ynum_to_num,
                                                        coeff_.xnum_yden_to_num,
                                                        coeff_.xden_ynum_to_num,
                                                        coeff_.xden_yden_to_num})},
                                        has_parallel_output);
                                }
                                else {
                                    // The third interval is also a single point.
                                    auto first_pt =
                                        detail::gosper_util<int_type>::kernel_of_rank1_transform(
                                            linear_fractional_transform{
                                                coeff_.xnum_ynum_to_num, coeff_.xnum_yden_to_num,
                                                coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num});
                                    auto second_pt =
                                        detail::gosper_util<int_type>::kernel_of_rank1_transform(
                                            linear_fractional_transform{
                                                coeff_.xnum_ynum_to_den, coeff_.xnum_yden_to_den,
                                                coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den});

                                    if (first_pt != second_pt) {
                                        return false;
                                    }
                                    else {
                                        return check_intersection(
                                            ycf_.current_interval(),
                                            cyclic_interval<projective_rational<int_type, int_type>,
                                                            cyclic_interval_type_t::single_point>{
                                                first_pt},
                                            has_parallel_output);
                                    }
                                }
                            }
                        }
                    };

                    while (has_singularity()) {
                        xcf_.update();
                        ycf_.update();
                    }
                }

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    // Step 3. Compare the floor of two endpoints and if they are equal, return.
                    auto update_and_callback = [&](auto&& common_floor) {
                        is_first_ = false;
                        coeff_.translate(-common_floor);
                        coeff_.reflect();
                        callback(partial_fraction_type{Unity{}, std::move(common_floor)});
                    };
                    enum class final_result { success, terminate, fail };
                    auto check_floor = [&](auto&& itv) -> final_result {
                        using enum cyclic_interval_type_t;
                        using itv_type = std::remove_cvref_t<decltype(itv)>;
                        static_assert(itv_type::interval_type() != empty &&
                                      itv_type::interval_type() != entire);
                        constexpr auto infinity = projective_rational{unity{}, zero{}};

                        if constexpr (itv_type::interval_type() == single_point) {
                            // If infinity is the only element in the range, then there is no
                            // further continued fraction coefficients.
                            if (itv.lower_bound() == infinity) {
                                return final_result::terminate;
                            }
                            else {
                                // If a finite rational number is the only element in the range,
                                // then just compute the continued fraction expansion of that
                                // number.
                                update_and_callback(util::div_floor(itv.lower_bound().numerator,
                                                                    itv.lower_bound().denominator));
                                return final_result::success;
                            }
                        }
                        else {
                            auto lower_bound = itv.lower_bound();
                            if (is_first_) {
                                // Cannot do anything if the range estimate contains the infinity.
                                if (lower_bound == infinity || itv.upper_bound() == infinity ||
                                    cyclic_order(lower_bound, infinity, itv.upper_bound())) {
                                    return final_result::fail;
                                }
                            }
                            else {
                                // Except for the first coefficient, the value must be in (1,infty].
                                // Still, we cannot do anything if infinity is strictly in the
                                // interval or the upper bound is equal to infinity.
                                if (itv.upper_bound() == infinity ||
                                    cyclic_order(lower_bound, infinity, itv.upper_bound())) {
                                    return final_result::fail;
                                }
                                // If lower bound is infinity, then the value might be infinity.
                                if (lower_bound == infinity) {
                                    if (itv.upper_bound() == projective_rational<unity, unity>{} ||
                                        cyclic_order(lower_bound, itv.upper_bound(),
                                                     projective_rational<unity, unity>{})) {
                                        return final_result::terminate;
                                    }
                                    return final_result::fail;
                                }

                                util::constexpr_assert(itv.upper_bound() !=
                                                       projective_rational<unity, unity>{});
                                if (cyclic_order(lower_bound, projective_rational<unity, unity>{},
                                                 itv.upper_bound())) {
                                    lower_bound.numerator = 1;
                                    lower_bound.denominator = 1;
                                }
                            }

                            // Compare the floor.
                            auto floor_lower =
                                util::div_floor(lower_bound.numerator, lower_bound.denominator);
                            auto floor_upper = util::div_floor(itv.upper_bound().numerator,
                                                               itv.upper_bound().denominator);
                            if (floor_lower == floor_upper) {
                                update_and_callback(std::move(floor_lower));
                                return final_result::success;
                            }
                            return final_result::fail;
                        }
                    };

                    // Step 2. Find the range of the linear fractional transform.
                    auto compute_range = [this, &check_floor](auto&& x_itv,
                                                              auto&& y_itv) -> final_result {
                        static_assert(x_itv.interval_type() != cyclic_interval_type_t::empty &&
                                      y_itv.interval_type() != cyclic_interval_type_t::empty);

                        if constexpr (x_itv.interval_type() == cyclic_interval_type_t::entire ||
                                      y_itv.interval_type() == cyclic_interval_type_t::entire) {
                            // If one of x,y ranges from the entire RP1, then the range of
                            // f(x,y) must be entire RP1.
                            return final_result::fail;
                        }
                        else {
                            using enum cyclic_interval_type_t;
                            using value_type = projective_rational<int_type, int_type>;

                            value_type const corners[4] = {
                                coeff_(x_itv.lower_bound(), y_itv.lower_bound()), // left-bottom
                                coeff_(x_itv.upper_bound(), y_itv.lower_bound()), // right-bottom
                                coeff_(x_itv.upper_bound(), y_itv.upper_bound()), // right-top
                                coeff_(x_itv.lower_bound(), y_itv.upper_bound())  // left-top
                            };
                            int const edge_directions[4] = {
                                // bottom
                                (x_itv.interval_type() == single_point ? 0 : 1) *
                                    linear_fractional_transform{
                                        coeff_.xnum_ynum_to_num * y_itv.lower_bound().numerator +
                                            coeff_.xnum_yden_to_num *
                                                y_itv.lower_bound().denominator,
                                        coeff_.xden_ynum_to_num * y_itv.lower_bound().numerator +
                                            coeff_.xden_yden_to_num *
                                                y_itv.lower_bound().denominator,
                                        coeff_.xnum_ynum_to_den * y_itv.lower_bound().numerator +
                                            coeff_.xnum_yden_to_den *
                                                y_itv.lower_bound().denominator,
                                        coeff_.xden_ynum_to_den * y_itv.lower_bound().numerator +
                                            coeff_.xden_yden_to_den *
                                                y_itv.lower_bound().denominator}
                                        .determinant_sign(),
                                // right
                                (y_itv.interval_type() == single_point ? 0 : 1) *
                                    linear_fractional_transform{
                                        coeff_.xnum_ynum_to_num * x_itv.upper_bound().numerator +
                                            coeff_.xden_ynum_to_num *
                                                x_itv.upper_bound().denominator,
                                        coeff_.xnum_yden_to_num * x_itv.upper_bound().numerator +
                                            coeff_.xden_yden_to_num *
                                                x_itv.upper_bound().denominator,
                                        coeff_.xnum_ynum_to_den * x_itv.upper_bound().numerator +
                                            coeff_.xden_ynum_to_den *
                                                x_itv.upper_bound().denominator,
                                        coeff_.xnum_yden_to_den * x_itv.upper_bound().numerator +
                                            coeff_.xden_yden_to_den *
                                                x_itv.upper_bound().denominator}
                                        .determinant_sign(),
                                // top
                                (x_itv.interval_type() == single_point ? 0 : -1) *
                                    linear_fractional_transform{
                                        coeff_.xnum_ynum_to_num * y_itv.upper_bound().numerator +
                                            coeff_.xnum_yden_to_num *
                                                y_itv.upper_bound().denominator,
                                        coeff_.xden_ynum_to_num * y_itv.upper_bound().numerator +
                                            coeff_.xden_yden_to_num *
                                                y_itv.upper_bound().denominator,
                                        coeff_.xnum_ynum_to_den * y_itv.upper_bound().numerator +
                                            coeff_.xnum_yden_to_den *
                                                y_itv.upper_bound().denominator,
                                        coeff_.xden_ynum_to_den * y_itv.upper_bound().numerator +
                                            coeff_.xden_yden_to_den *
                                                y_itv.upper_bound().denominator}
                                        .determinant_sign(),
                                // left
                                (y_itv.interval_type() == single_point ? 0 : -1) *
                                    linear_fractional_transform{
                                        coeff_.xnum_ynum_to_num * x_itv.lower_bound().numerator +
                                            coeff_.xden_ynum_to_num *
                                                x_itv.lower_bound().denominator,
                                        coeff_.xnum_yden_to_num * x_itv.lower_bound().numerator +
                                            coeff_.xden_yden_to_num *
                                                x_itv.lower_bound().denominator,
                                        coeff_.xnum_ynum_to_den * x_itv.lower_bound().numerator +
                                            coeff_.xden_ynum_to_den *
                                                x_itv.lower_bound().denominator,
                                        coeff_.xnum_yden_to_den * x_itv.lower_bound().numerator +
                                            coeff_.xden_yden_to_den *
                                                x_itv.lower_bound().denominator}
                                        .determinant_sign()};

                            auto passes_through = [&corners,
                                                   &edge_directions](unsigned int edge_idx,
                                                                     unsigned int corner_idx) {
                                auto const& corner = corners[corner_idx % 4];
                                auto const& first_end = corners[edge_idx % 4];
                                auto const& second_end = corners[(edge_idx + 1) % 4];
                                auto const& edge_direction = edge_directions[edge_idx % 4];

                                if (edge_direction == 0) {
                                    return false;
                                }
                                if (edge_direction > 0) {
                                    return corner == first_end || corner == second_end ||
                                           cyclic_order(first_end, corner, second_end);
                                }
                                else {
                                    return corner == first_end || corner == second_end ||
                                           cyclic_order(second_end, corner, first_end);
                                }
                            };

                            // Cases when f is constant on the region.
                            if (edge_directions[0] == 0 && edge_directions[1] == 0 &&
                                edge_directions[2] == 0 && edge_directions[3] == 0) {
                                return check_floor(
                                    cyclic_interval<value_type const&, single_point>{corners[0]});
                            }

                            int const plus_count = (edge_directions[0] >= 0 ? 1 : 0) +
                                                   (edge_directions[1] >= 0 ? 1 : 0) +
                                                   (edge_directions[2] >= 0 ? 1 : 0) +
                                                   (edge_directions[3] >= 0 ? 1 : 0);

                            // (+,+,+,+) or (-,-,-,-)
                            if (plus_count == 4 || plus_count == 0) {
                                // The range must be RP1.
                                return final_result::fail;
                            }
                            // (+,+,+,-)
                            else if (plus_count == 3) {
                                // Locate the unique negative at the end.
                                unsigned starting_pos = edge_directions[0] < 0   ? 1
                                                        : edge_directions[1] < 0 ? 2
                                                        : edge_directions[2] < 0 ? 3
                                                                                 : 0;

                                // Either [corners[0], corners[3]] or RP1.
                                if (passes_through(starting_pos + 1, starting_pos) ||
                                    passes_through(starting_pos + 2, starting_pos)) {
                                    // The range must be RP1.
                                    return final_result::fail;
                                }
                                else {
                                    return check_floor(cyclic_interval<value_type const&, closed>{
                                        corners[starting_pos], corners[(starting_pos + 3) % 4]});
                                }
                            }
                            // (-,-,-,+)
                            else if (plus_count == 1) {
                                // Locate the unique positive at the end.
                                unsigned int starting_pos = edge_directions[0] >= 0   ? 1
                                                            : edge_directions[1] >= 0 ? 2
                                                            : edge_directions[2] >= 0 ? 3
                                                                                      : 0;

                                // Either [corners[3], corners[0]] or RP1.
                                if (passes_through(starting_pos + 1, starting_pos) ||
                                    passes_through(starting_pos + 2, starting_pos)) {
                                    // The range must be RP1.
                                    return final_result::fail;
                                }
                                else {
                                    return check_floor(cyclic_interval<value_type const&, closed>{
                                        corners[(starting_pos + 3) % 4], corners[starting_pos]});
                                }
                            }
                            // (+,+,-,-)
                            else if (edge_directions[0] != edge_directions[2]) {
                                // Locate two successive positive at the beginning.
                                unsigned int starting_pos =
                                    (edge_directions[0] >= 0 && edge_directions[1] >= 0)   ? 0
                                    : (edge_directions[1] >= 0 && edge_directions[2] >= 0) ? 1
                                    : (edge_directions[2] >= 0 && edge_directions[3] >= 0) ? 2
                                                                                           : 3;

                                // Either [corners[0], corners[2]] or RP1.
                                if (passes_through(starting_pos + 1, starting_pos) ||
                                    passes_through(starting_pos + 2, starting_pos)) {
                                    // The range must be RP1.
                                    return final_result::fail;
                                }
                                else {
                                    return check_floor(cyclic_interval<value_type const&, closed>{
                                        corners[starting_pos], corners[(starting_pos + 2) % 4]});
                                }
                            }
                            // (+,-,+,-)
                            else {
                                // Locate any positive at the beginning.
                                unsigned int starting_pos = edge_directions[0] >= 0 ? 0 : 1;

                                // 6 different cases.
                                if (passes_through(starting_pos + 1, starting_pos)) {
                                    if (passes_through(starting_pos + 2, starting_pos + 1)) {
                                        return check_floor(
                                            cyclic_interval<value_type const&, closed>{
                                                corners[(starting_pos + 2) % 4],
                                                corners[(starting_pos + 3) % 4]});
                                    }
                                    else if (passes_through(starting_pos + 2, starting_pos)) {
                                        return check_floor(
                                            cyclic_interval<value_type const&, closed>{
                                                corners[(starting_pos + 2) % 4],
                                                corners[(starting_pos + 1) % 4]});
                                    }
                                    else {
                                        // The range must be RP1.
                                        return final_result::fail;
                                    }
                                }
                                else if (passes_through(starting_pos + 2, starting_pos)) {
                                    // The range must be RP1.
                                    return final_result::fail;
                                }
                                else if (passes_through(starting_pos + 2, starting_pos + 1)) {
                                    return check_floor(cyclic_interval<value_type const&, closed>{
                                        corners[starting_pos], corners[(starting_pos + 3) % 4]});
                                }
                                else {
                                    return check_floor(cyclic_interval<value_type const&, closed>{
                                        corners[starting_pos], corners[(starting_pos + 1) % 4]});
                                }
                            }
                        }
                    };

                    while (true) {
                        auto const result =
                            xcf_.current_interval().visit([&](auto&& x_itv) -> final_result {
                                return ycf_.current_interval().visit(
                                    [&](auto&& y_itv) -> final_result {
                                        return compute_range(x_itv, y_itv);
                                    });
                            });

                        switch (result) {
                        case final_result::success:
                        case final_result::terminate:
                            return;
                        case final_result::fail:;
                        }

                        xcf_.update();
                        ycf_.update();
                    }
                }
            };
        }
    }
}

#endif
