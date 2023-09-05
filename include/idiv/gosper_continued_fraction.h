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

#ifndef JKJ_HEADER_GOSPER_CONTINUED_FRACTION
#define JKJ_HEADER_GOSPER_CONTINUED_FRACTION

#include "continued_fraction.h"
#include "frac.h"
#include <cstdlib>
#include <type_traits>

namespace jkj {
    namespace cntfrc {
        template <class ContinuedFractionImpl, class Unity = unity,
                  template <class> class... Mixins>
        class unary_gosper;

        template <class ContinuedFractionImpl, class Unity, template <class> class... Mixins>
        struct continued_fraction_traits<unary_gosper<ContinuedFractionImpl, Unity, Mixins...>> {
            using int_type = decltype(ContinuedFractionImpl::convergent_type::numerator);
            using uint_type = decltype(ContinuedFractionImpl::convergent_type::denominator);
            using partial_fraction_type = frac<Unity, int_type>;
            using convergent_type = typename ContinuedFractionImpl::convergent_type;
            using interval_type = variable_shape_cyclic_interval<
                projective_rational<int_type, int_type>, cyclic_interval_type_t::single_point,
                cyclic_interval_type_t::closed, cyclic_interval_type_t::entire>;
        };

        template <class ContinuedFractionImpl, class Unity, template <class> class... Mixins>
        class unary_gosper
            : public continued_fraction_base<unary_gosper<ContinuedFractionImpl, Unity, Mixins...>,
                                             Mixins...> {
            using crtp_base =
                continued_fraction_base<unary_gosper<ContinuedFractionImpl, Unity, Mixins...>,
                                        Mixins...>;
            friend crtp_base;

        public:
            using int_type = typename crtp_base::traits_type::int_type;
            using uint_type = typename crtp_base::traits_type::uint_type;
            using partial_fraction_type = typename crtp_base::traits_type::partial_fraction_type;
            using convergent_type = typename crtp_base::traits_type::convergent_type;
            using interval_type = typename crtp_base::traits_type::interval_type;
            using internal_continued_fraction_impl_type = ContinuedFractionImpl;

        private:
            ContinuedFractionImpl cf_;
            linear_fractional_transform<int_type> coeff_;
            int determinant_sign_ = 0;

            template <class Functor>
            constexpr bool with_next_partial_fraction(Functor&& f) {
                while (true) {
                    // Step 1. Get away from singularities.
                    if (determinant_sign_ == 0) {
                        // The degenerate case a = b = c = d = 0 is disallowed.
                        util::constexpr_assert(!is_zero(coeff_.num_to_num) ||
                                               !is_zero(coeff_.den_to_num));

                        // The unique singularity is at [-b:a]=[-d:c].
                        // If the singularity is contained in the current interval, then refine the
                        // interval.
                        if (cf_.current_interval().contains(projective_rational{
                                invert_sign(coeff_.den_to_num), coeff_.num_to_num})) {
                            // We do not allow the input to be at the singularity.
                            util::constexpr_assert(cf_.current_interval().interval_type() !=
                                                   cyclic_interval_type_t::single_point);

                            cf_.update();
                            continue;
                        }
                    }

                    // Step 2. Find the range of the linear fractional transform.
                    interval_type range = [&]() -> interval_type {
                        // Single point [a:c].
                        if (determinant_sign_ == 0) {
                            return cyclic_interval<projective_rational<int_type, int_type>,
                                                   cyclic_interval_type_t::single_point>{
                                projective_rational<int_type, int_type>{coeff_.num_to_num,
                                                                        coeff_.num_to_den}};
                        }

                        return cf_.current_interval().visit([&](auto&& itv) -> interval_type {
                            using enum cyclic_interval_type_t;
                            static_assert(!interval_type::is_allowed_interval_type(empty));
                            if constexpr (itv.interval_type() == single_point) {
                                return cyclic_interval<projective_rational<int_type, int_type>,
                                                       single_point>{coeff_(itv.lower_bound())};
                            }
                            else if constexpr (itv.interval_type() == entire) {
                                return cyclic_interval<projective_rational<int_type, int_type>,
                                                       entire>{};
                            }
                            else {
                                auto first_end = coeff_(itv.lower_bound());
                                auto second_end = coeff_(itv.upper_bound());
                                if (determinant_sign_ < 0) {
                                    using std::swap;
                                    swap(first_end, second_end);
                                }
                                return cyclic_interval<projective_rational<int_type, int_type>,
                                                       closed>{std::move(first_end),
                                                               std::move(second_end)};
                            }
                        });
                    }();

                    // Step 3. See if the range is a bounded interval inside the real line.
                    if (range.contains(projective_rational{unity{}, zero{}})) {
                        // If infinity is the only point in the range, then terminate.
                        if (range.interval_type() == cyclic_interval_type_t::single_point) {
                            return false;
                        }

                        cf_.update();
                        continue;
                    }

                    // Step 4. Check if the floors of the endpoints agree.
                    int_type lower_floor{};
                    if (range.visit([&](auto&& itv) {
                            using enum cyclic_interval_type_t;
                            if constexpr (itv.interval_type() == empty ||
                                          itv.interval_type() == entire) {
                                util::constexpr_assert(false);
                                return false;
                            }
                            else {
                                lower_floor = div_floor(itv.lower_bound().numerator,
                                                        itv.lower_bound().denominator);
                                if constexpr (itv.interval_type() == single_point) {
                                    return true;
                                }
                                else {
                                    auto upper_floor = div_floor(itv.upper_bound().numerator,
                                                                 itv.upper_bound().denominator);
                                    return lower_floor == upper_floor;
                                }
                            }
                        })) {
                        // Found a new coefficient.
                        coeff_.translate(-lower_floor);
                        coeff_.reflect();
                        determinant_sign_ *= -1;
                        f(partial_fraction_type{Unity{}, static_cast<int_type&&>(lower_floor)});
                        return true;
                    }

                    cf_.update();
                }
            }

        public:
            struct default_mixin_initializer {
                static constexpr partial_fraction_type initial_partial_fraction() {
                    return {Unity{}, int_type{0}};
                }
                static constexpr interval_type initial_interval() {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }
            };

            template <class MixinInitializer = default_mixin_initializer>
            constexpr unary_gosper(ContinuedFractionImpl cf,
                                   linear_fractional_transform<int_type> coeff,
                                   MixinInitializer&& mixin_initializer = {})
                : crtp_base{mixin_initializer}, cf_{static_cast<ContinuedFractionImpl>(cf)},
                  coeff_{static_cast<linear_fractional_transform<int_type>&&>(coeff)} {
                determinant_sign_ =
                    util::strong_order_to_int(coeff_.num_to_num * coeff_.den_to_den <=>
                                              coeff_.den_to_num * coeff_.num_to_den);
            }
        };

        template <class XContinuedFractionImpl, class YContinuedFractionImpl, class Unity = unity,
                  template <class> class... Mixins>
            requires std::is_same_v<typename XContinuedFractionImpl::convergent_type,
                                    typename YContinuedFractionImpl::convergent_type>
        class binary_gosper;

        template <class XContinuedFractionImpl, class YContinuedFractionImpl, class Unity,
                  template <class> class... Mixins>
        struct continued_fraction_traits<
            binary_gosper<XContinuedFractionImpl, YContinuedFractionImpl, Unity, Mixins...>> {
            using int_type = decltype(XContinuedFractionImpl::convergent_type::numerator);
            using uint_type = decltype(XContinuedFractionImpl::convergent_type::denominator);
            using partial_fraction_type = frac<Unity, int_type>;
            using convergent_type = typename XContinuedFractionImpl::convergent_type;
            using interval_type = variable_shape_cyclic_interval<
                convergent_type, cyclic_interval_type_t::single_point,
                cyclic_interval_type_t::left_open_right_closed,
                cyclic_interval_type_t::left_closed_right_open, cyclic_interval_type_t::entire>;
        };

        template <class XContinuedFractionImpl, class YContinuedFractionImpl, class Unity,
                  template <class> class... Mixins>
            requires std::is_same_v<typename XContinuedFractionImpl::convergent_type,
                                    typename YContinuedFractionImpl::convergent_type>
        class binary_gosper
            : public continued_fraction_base<
                  binary_gosper<XContinuedFractionImpl, YContinuedFractionImpl, Unity, Mixins...>,
                  Mixins...> {
            using crtp_base = continued_fraction_base<
                binary_gosper<XContinuedFractionImpl, YContinuedFractionImpl, Unity, Mixins...>,
                Mixins...>;
            friend crtp_base;

        public:
            using int_type = typename crtp_base::traits_type::int_type;
            using uint_type = typename crtp_base::traits_type::uint_type;
            using partial_fraction_type = typename crtp_base::traits_type::partial_fraction_type;
            using convergent_type = typename crtp_base::traits_type::convergent_type;
            using interval_type = typename crtp_base::traits_type::interval_type;
            using first_internal_continued_fraction_impl_type = XContinuedFractionImpl;
            using second_internal_continued_fraction_impl_type = YContinuedFractionImpl;

        private:
            XContinuedFractionImpl xcf_;
            YContinuedFractionImpl ycf_;
            bilinear_fractional_transform<int_type> coeff_;
            int det_sign_num_bilinear_form_ = 0;
            int nullity_num_bilinear_form_ = 0;
            int det_sign_den_bilinear_form_ = 0;
            int nullity_den_bilinear_form_ = 0;
            struct sing_det_form_t {
                int_type a, b, d;
            } sing_det_form_{};
            bool globally_well_defined_ = false;

            constexpr void update_coeff_info() {
                auto calculate_det_sign_nullity = [](auto&& a, auto&& b, auto&& c, auto&& d,
                                                     auto& det_sign, auto& nullity) {
                    auto const result = a * d <=> b * c;

                    if (result < 0) {
                        det_sign = -1;
                        nullity = 0;
                    }
                    if (result > 0) {
                        det_sign = 1;
                        nullity = 0;
                    }
                    else {
                        det_sign = 0;
                        nullity = (is_zero(a) && is_zero(b) && is_zero(c) && is_zero(d) ? 2 : 1);
                    }
                };

                calculate_det_sign_nullity(coeff_.xnum_ynum_to_num, coeff_.xnum_yden_to_num,
                                           coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num,
                                           det_sign_num_bilinear_form_, nullity_num_bilinear_form_);

                calculate_det_sign_nullity(coeff_.xnum_ynum_to_den, coeff_.xnum_yden_to_den,
                                           coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den,
                                           det_sign_den_bilinear_form_, nullity_den_bilinear_form_);

                // L = A^T RB - B^T RA = 2sym(ag-ce  ah-cf  \\ bg-de  bh-df).
                sing_det_form_.a = ((coeff_.xnum_ynum_to_num * coeff_.xden_ynum_to_den -
                                     coeff_.xden_ynum_to_num * coeff_.xnum_ynum_to_den)
                                    << 1);
                sing_det_form_.d = ((coeff_.xnum_yden_to_num * coeff_.xden_yden_to_den -
                                     coeff_.xden_yden_to_num * coeff_.xnum_yden_to_den)
                                    << 1);
                sing_det_form_.b = coeff_.xnum_ynum_to_num * coeff_.xden_yden_to_den +
                                   coeff_.xnum_yden_to_num * coeff_.xden_ynum_to_den -
                                   coeff_.xden_ynum_to_num * coeff_.xnum_yden_to_den -
                                   coeff_.xden_yden_to_num * coeff_.xnum_ynum_to_den;

                globally_well_defined_ =
                    ((sing_det_form_.a * sing_det_form_.d) > (sing_det_form_.b * sing_det_form_.b));

                // Degenerate case A = B = 0 is disallowed.
                util::constexpr_assert(nullity_num_bilinear_form_ != 2 ||
                                       nullity_den_bilinear_form_ != 2);
            }

            // Find the kernel of a rank-1 linear fractional transform.
            template <class NumNum, class DenNum, class NumDen, class DenDen>
            static constexpr projective_rational<int_type, int_type> kernel_of_rank1_transform(
                linear_fractional_transform<NumNum, DenNum, NumDen, DenDen> const& transform) {
                return is_zero(transform.num_to_num) && is_zero(transform.den_to_num)
                           ? projective_rational{-int_type{transform.den_to_den},
                                                 int_type{transform.num_to_den}}
                           : projective_rational{-int_type{transform.den_to_num},
                                                 int_type{transform.num_to_num}};
            }
            // Find the range of a rank-1 linear fractional transform.
            template <class NumNum, class DenNum, class NumDen, class DenDen>
            static constexpr projective_rational<int_type, int_type> range_of_rank1_transform(
                linear_fractional_transform<NumNum, DenNum, NumDen, DenDen> const& transform) {
                return is_zero(transform.num_to_num) && is_zero(transform.num_to_den)
                           ? projective_rational{int_type{transform.den_to_num},
                                                 int_type{transform.den_to_den}}
                           : projective_rational{int_type{transform.num_to_num},
                                                 int_type{transform.num_to_den}};
            }

            // Precondition: itv is contained in the domain of transform.
            template <class NumNum, class DenNum, class NumDen, class DenDen,
                      class InputIntervalType>
            static constexpr variable_shape_cyclic_interval<
                projective_rational<int_type, int_type>, cyclic_interval_type_t::single_point,
                cyclic_interval_type_t::closed, cyclic_interval_type_t::entire>
            map_cyclic_interval(
                linear_fractional_transform<NumNum, DenNum, NumDen, DenDen> const& transform,
                InputIntervalType const& itv, int determinant_sign) {
                using return_type = variable_shape_cyclic_interval<
                    projective_rational<int_type, int_type>, cyclic_interval_type_t::single_point,
                    cyclic_interval_type_t::closed, cyclic_interval_type_t::entire>;

                return itv.visit([&transform, determinant_sign](auto&& itv) -> return_type {
                    using enum cyclic_interval_type_t;
                    using value_type = projective_rational<int_type, int_type>;
                    constexpr auto itv_type = itv.interval_type();
                    static_assert(itv_type != empty);

                    if constexpr (itv_type == entire) {
                        util::constexpr_assert(determinant_sign != 0);
                        return cyclic_interval<value_type, entire>{};
                    }
                    else if constexpr (itv_type == single_point) {
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
                                range_of_rank1_transform(transform)};
                        }
                    }
                });
            }
            template <class NumNum, class DenNum, class NumDen, class DenDen,
                      class InputIntervalType>
            static constexpr variable_shape_cyclic_interval<
                projective_rational<int_type, int_type>, cyclic_interval_type_t::single_point,
                cyclic_interval_type_t::closed, cyclic_interval_type_t::entire>
            map_cyclic_interval(
                linear_fractional_transform<NumNum, DenNum, NumDen, DenDen> const& transform,
                InputIntervalType const& itv) {
                return map_cyclic_interval(
                    transform, itv,
                    util::strong_order_to_int(transform.num_to_num * transform.den_to_den <=>
                                              transform.den_to_num * transform.num_to_den));
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
                        if (itv2.contains(itv1.lower_bound())) {
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
                                if (itv1.contains(itv2.lower_bound())) {
                                    return f(single_point_interval{itv2.lower_bound()});
                                }
                                else {
                                    return f(empty_interval{});
                                }
                            }
                            else {
                                // 14 different cases.
                                auto&& a = itv1.lower_bound();
                                auto&& b = itv1.upper_bound();
                                auto&& c = itv2.lower_bound();
                                auto&& d = itv2.upper_bound();
                                if (a == c) {
                                    if (cyclic_order(a, b, d)) {
                                        return f(closed_interval{a, b});
                                    }
                                    else {
                                        // [a,d,b]
                                        return f(closed_interval{a, d});
                                    }
                                }
                                else if (a == d) {
                                    if (cyclic_order(a, b, c)) {
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
                                        return f(closed_interval{a, d}) && f(closed_interval{c, b});
                                    }
                                }
                            }
                        });
                    }
                });
            }

            // Check if there exists [x] in itv such that Ax and Bx are parallel.
            // Assumes itv is closed.
            template <class Value, cyclic_interval_type_t it>
            constexpr bool has_parallel_output(cyclic_interval<Value, it> const& itv) const {
                using enum cyclic_interval_type_t;
                constexpr auto itv_type = itv.interval_type();

                if constexpr (itv_type == empty) {
                    return false;
                }
                else if constexpr (itv_type == entire) {
                    return !globally_well_defined_;
                }
                else if constexpr (itv_type == single_point) {
                    return is_zero(sing_det_form_.a * itv.lower_bound().numerator *
                                       itv.lower_bound().numerator +
                                   ((sing_det_form_.b * itv.lower_bound().numerator *
                                     itv.lower_bound().denominator)
                                    << 1) +
                                   sing_det_form_.d * itv.lower_bound().denominator *
                                       itv.lower_bound().denominator);
                }
                else {
                    auto Luu = sing_det_form_.a * itv.lower_bound().numerator *
                                   itv.lower_bound().numerator +
                               ((sing_det_form_.b * itv.lower_bound().numerator *
                                 itv.lower_bound().denominator)
                                << 1) +
                               sing_det_form_.d * itv.lower_bound().denominator *
                                   itv.lower_bound().denominator;
                    auto Lvv = sing_det_form_.a * itv.upper_bound().numerator *
                                   itv.upper_bound().numerator +
                               ((sing_det_form_.b * itv.upper_bound().numerator *
                                 itv.upper_bound().denominator)
                                << 1) +
                               sing_det_form_.d * itv.upper_bound().denominator *
                                   itv.upper_bound().denominator;

                    int sign_Luu = is_zero(Luu) ? 0 : is_nonnegative(Luu) ? 1 : -1;
                    int sign_Lvv = is_zero(Lvv) ? 0 : is_nonnegative(Lvv) ? 1 : -1;
                    if (sign_Luu * sign_Lvv <= 0) {
                        return true;
                    }

                    auto Luv = sing_det_form_.a * itv.lower_bound().numerator *
                                   itv.upper_bound().numerator +
                               sing_det_form_.b * itv.lower_bound().numerator *
                                   itv.upper_bound().denominator +
                               sing_det_form_.b * itv.lower_bound().denominator *
                                   itv.upper_bound().numerator +
                               sing_det_form_.d * itv.lower_bound().denominator *
                                   itv.upper_bound().denominator;
                    int sign_Luv = is_zero(Luv) ? 0 : is_nonnegative(Luv) ? 1 : -1;
                    return sign_Luu * sign_Luv < 0;
                }
            }

            constexpr bool has_singularity() const {
                if (globally_well_defined_) {
                    return false;
                }

                auto check_interval = [&](auto&& itv) -> bool { return has_parallel_output(itv); };

                // Compute the set S = J \cap (T_RA^-1[I] \cup (dom(T_A))^c)
                // \cap (T_RB^-1[I] \cup (dom(T_B))^c).
                // We already ruled out the case A = 0 or B = 0 at the start of
                // with_next_partial_fraction().

                // A is invertible.
                if (nullity_num_bilinear_form_ == 0) {
                    // Can ignore the third term from the intersection.
                    // RA^T = (b d \\ -a -c)
                    return check_intersection(
                        ycf_.current_interval(),
                        map_cyclic_interval(linear_fractional_transform{coeff_.xnum_yden_to_num,
                                                                        coeff_.xden_yden_to_num,
                                                                        -coeff_.xnum_ynum_to_num,
                                                                        -coeff_.xden_ynum_to_num},
                                            xcf_.current_interval(), det_sign_num_bilinear_form_),
                        check_interval);
                }
                // B is invertible.
                else if (nullity_den_bilinear_form_ == 0) {
                    // Can ignore the second term from the intersection.
                    // RB^T = (f h \\ -e -g)
                    return check_intersection(
                        ycf_.current_interval(),
                        map_cyclic_interval(linear_fractional_transform{coeff_.xnum_yden_to_den,
                                                                        coeff_.xden_yden_to_den,
                                                                        -coeff_.xnum_ynum_to_den,
                                                                        -coeff_.xden_ynum_to_den},
                                            xcf_.current_interval(), det_sign_den_bilinear_form_),
                        check_interval);
                }
                // Both are not invertible.
                else {
                    util::constexpr_assert(nullity_num_bilinear_form_ == 1 &&
                                           nullity_den_bilinear_form_ == 1);

                    // RA = (c d \\ -a -b)
                    if (xcf_.current_interval().contains(
                            range_of_rank1_transform(linear_fractional_transform{
                                coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num,
                                -coeff_.xnum_ynum_to_num, -coeff_.xnum_yden_to_num}))) {
                        // The second interval is the entire RP1.
                        // RB = (g h \\ -e -f)
                        if (xcf_.current_interval().contains(
                                range_of_rank1_transform(linear_fractional_transform{
                                    coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den,
                                    -coeff_.xnum_ynum_to_den, -coeff_.xnum_yden_to_den}))) {
                            // The third interval is also the entire RP1.
                            return ycf_.current_interval().visit(check_interval);
                        }
                        else {
                            // The third interval is a single point.
                            return check_intersection(
                                ycf_.current_interval(),
                                cyclic_interval<projective_rational<int_type, int_type>,
                                                cyclic_interval_type_t::single_point>{
                                    kernel_of_rank1_transform(linear_fractional_transform{
                                        coeff_.xnum_ynum_to_den, coeff_.xnum_yden_to_den,
                                        coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den})},
                                check_interval);
                        }
                    }
                    else {
                        // The second interval is a single point.
                        // RB = (g h \\ -e -f)
                        if (xcf_.current_interval().contains(
                                range_of_rank1_transform(linear_fractional_transform{
                                    coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den,
                                    -coeff_.xnum_ynum_to_den, -coeff_.xnum_yden_to_den}))) {
                            // The third interval is the entire RP1.
                            return check_intersection(
                                ycf_.current_interval(),
                                cyclic_interval<projective_rational<int_type, int_type>,
                                                cyclic_interval_type_t::single_point>{
                                    kernel_of_rank1_transform(linear_fractional_transform{
                                        coeff_.xnum_ynum_to_num, coeff_.xnum_yden_to_num,
                                        coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num})},
                                check_interval);
                        }
                        else {
                            // The third interval is also a single point.
                            auto first_pt = kernel_of_rank1_transform(linear_fractional_transform{
                                coeff_.xnum_ynum_to_num, coeff_.xnum_yden_to_num,
                                coeff_.xden_ynum_to_num, coeff_.xden_yden_to_num});
                            auto second_pt = kernel_of_rank1_transform(linear_fractional_transform{
                                coeff_.xnum_ynum_to_den, coeff_.xnum_yden_to_den,
                                coeff_.xden_ynum_to_den, coeff_.xden_yden_to_den});

                            if (first_pt != second_pt) {
                                return false;
                            }
                            else {
                                return check_intersection(
                                    ycf_.current_interval(),
                                    cyclic_interval<projective_rational<int_type, int_type>,
                                                    cyclic_interval_type_t::single_point>{first_pt},
                                    check_interval);
                            }
                        }
                    }
                }
            }

            // Assumes that itv1 and itv2 are not disjoint, and are both either single points,
            // closed intervals, or entire RP1.
            template <class InputIntervalType1, class InputIntervalType2, class Functor>
            static constexpr decltype(auto) check_union(InputIntervalType1 const& itv1,
                                                        InputIntervalType2 const& itv2,
                                                        Functor&& f) {
                return itv1.visit([&itv2, &f](auto&& itv1) {
                    using enum cyclic_interval_type_t;
                    using value_type = projective_rational<int_type, int_type>;
                    using closed_interval = cyclic_interval<value_type, closed>;
                    using entire_interval = cyclic_interval<value_type, entire>;

                    constexpr auto itv1_type = itv1.interval_type();

                    if constexpr (itv1_type == empty) {
                        return itv2.visit(f);
                    }
                    else if constexpr (itv1_type == entire) {
                        return f(entire_interval{});
                    }
                    else if constexpr (itv1_type == single_point) {

                        util::constexpr_assert(itv2.contains(itv1.lower_bound()));
                        return itv2.visit(f);
                    }
                    else {
                        return itv2.visit([&itv1, &f](auto&& itv2) {
                            constexpr auto itv2_type = itv2.interval_type();

                            if constexpr (itv2_type == empty) {
                                return f(itv1);
                            }
                            else if constexpr (itv2_type == entire) {
                                return f(entire_interval{});
                            }
                            else if constexpr (itv2_type == single_point) {
                                util::constexpr_assert(itv1.contains(itv2.lower_bound()));
                                return f(itv1);
                            }
                            else {
                                // 14 different cases.
                                auto&& a = itv1.lower_bound();
                                auto&& b = itv1.upper_bound();
                                auto&& c = itv2.lower_bound();
                                auto&& d = itv2.upper_bound();
                                if (a == c) {
                                    if (cyclic_order(a, b, d)) {
                                        return f(closed_interval{a, d});
                                    }
                                    else {
                                        // [a,d,b]
                                        return f(closed_interval{a, b});
                                    }
                                }
                                else if (a == d) {
                                    if (cyclic_order(a, b, c)) {
                                        return f(closed_interval{c, b});
                                    }
                                    else {
                                        // [a,c,b]
                                        return f(entire_interval{});
                                    }
                                }
                                else if (b == c) {
                                    if (cyclic_order(a, b, d)) {
                                        return f(closed_interval{a, d});
                                    }
                                    else {
                                        // [a,d,b]
                                        return f(entire_interval{});
                                    }
                                }
                                else if (b == d) {
                                    if (cyclic_order(a, b, c)) {
                                        return f(closed_interval{c, b});
                                    }
                                    else {
                                        // [a,c,b]
                                        return f(closed_interval{a, b});
                                    }
                                }
                                else if (cyclic_order(a, b, c)) {
                                    util::constexpr_assert(cyclic_order(a, d, c));
                                    if (cyclic_order(a, b, d)) {
                                        // [a,b,d,c]
                                        return f(closed_interval{c, d});
                                    }
                                    else {
                                        // [a,d,b,c]
                                        return f(closed_interval{c, b});
                                    }
                                }
                                else {
                                    if (cyclic_order(b, d, a)) {
                                        // [a,c,b,d]
                                        return f(closed_interval{a, d});
                                    }
                                    else if (cyclic_order(a, c, d)) {
                                        // [a,c,d,b]
                                        return f(closed_interval{a, b});
                                    }
                                    else {
                                        // [a,d,c,b]
                                        return f(entire_interval{});
                                    }
                                }
                            }
                        });
                    }
                });
            }

            template <class Functor>
            constexpr bool with_next_partial_fraction(Functor&& f) {
                // Stops when the denominator coefficients are all zero.
                // (x,y) being not located at a singularity is a precondition, so we do not check if
                // the numerator is also zero.
                if (nullity_den_bilinear_form_ == 2) {
                    return false;
                }

                // If the numerator coefficients are all zero, then we output 0.
                if (nullity_num_bilinear_form_ == 2) {
                    coeff_.reflect();
                    update_coeff_info();
                    f(partial_fraction_type{Unity{}, 0u});
                    return true;
                }

                while (true) {
                    // Step 1. Get away from singularities.
                    if (has_singularity()) {
                        xcf_.update();
                        ycf_.update();
                        continue;
                    }

                    // Step 3. Compare the floor of two endpoints and if they are equal, return.
                    auto update_and_callback = [&](auto&& common_floor) {
                        coeff_.translate(-common_floor);
                        coeff_.reflect();
                        update_coeff_info();
                        f(partial_fraction_type{
                            Unity{}, static_cast<decltype(common_floor)&&>(common_floor)});
                    };
                    enum class final_result { success, retry, fail };
                    auto check_floor = [&](auto&& itv) -> final_result {
                        using enum cyclic_interval_type_t;
                        constexpr auto itv_type = itv.interval_type();
                        static_assert(itv_type != empty);

                        if constexpr (itv_type == entire) {
                            // Cannot do anything if the range estimate is the entire RP1.
                            return final_result::retry;
                        }
                        else if constexpr (itv_type == single_point) {
                            // If infinity is the only element in the range, then there is no
                            // further continued fraction coefficients.
                            if (itv.lower_bound() == projective_rational{unity{}, zero{}}) {
                                return final_result::fail;
                            }
                            else {
                                // If a finite rational number is the only element in the range,
                                // then just compute the continued fraction expansion of that
                                // number.
                                update_and_callback(div_floor(itv.lower_bound().numerator,
                                                              itv.lower_bound().denominator));
                                return final_result::success;
                            }
                        }
                        else {
                            // Cannot do anything if the range estimate contains the infinity.
                            if (itv.lower_bound() == projective_rational{unity{}, zero{}} ||
                                itv.upper_bound() == projective_rational{unity{}, zero{}} ||
                                cyclic_order(itv.lower_bound(),
                                             projective_rational{unity{}, zero{}},
                                             itv.upper_bound())) {
                                return final_result::retry;
                            }

                            // Otherwise, compare the floor.
                            auto floor_lower = div_floor(itv.lower_bound().numerator,
                                                         itv.lower_bound().denominator);
                            auto floor_upper = div_floor(itv.upper_bound().numerator,
                                                         itv.upper_bound().denominator);
                            if (floor_lower == floor_upper) {
                                update_and_callback(
                                    static_cast<decltype(floor_lower)&&>(floor_lower));
                                return final_result::success;
                            }
                            return final_result::retry;
                        }
                    };

                    // Step 2. Find the range of the linear fractional transform.
                    auto compute_range = [this, &check_floor](auto&& x_itv,
                                                              auto&& y_itv) -> final_result {
                        static_assert(x_itv.interval_type() != cyclic_interval_type_t::empty &&
                                      y_itv.interval_type() != cyclic_interval_type_t::empty);

                        if constexpr (x_itv.interval_type() == cyclic_interval_type_t::entire ||
                                      y_itv.interval_type() == cyclic_interval_type_t::entire) {
                            // If one of x,y ranges from the entire RP1, just estimate the range of
                            // f(x,y) as the entire RP1, although that is not a tight estimate.
                            return final_result::retry;
                        }
                        else {
                            auto itv1 = map_cyclic_interval(
                                linear_fractional_transform{
                                    coeff_.xnum_ynum_to_num * y_itv.lower_bound().numerator +
                                        coeff_.xnum_yden_to_num * y_itv.lower_bound().denominator,
                                    coeff_.xden_ynum_to_num * y_itv.lower_bound().numerator +
                                        coeff_.xden_yden_to_num * y_itv.lower_bound().denominator,
                                    coeff_.xnum_ynum_to_den * y_itv.lower_bound().numerator +
                                        coeff_.xnum_yden_to_den * y_itv.lower_bound().denominator,
                                    coeff_.xden_ynum_to_den * y_itv.lower_bound().numerator +
                                        coeff_.xden_yden_to_den * y_itv.lower_bound().denominator},
                                x_itv);
                            auto itv2 = map_cyclic_interval(
                                linear_fractional_transform{
                                    coeff_.xnum_ynum_to_num * x_itv.upper_bound().numerator +
                                        coeff_.xden_ynum_to_num * x_itv.upper_bound().denominator,
                                    coeff_.xnum_yden_to_num * x_itv.upper_bound().numerator +
                                        coeff_.xden_yden_to_num * x_itv.upper_bound().denominator,
                                    coeff_.xnum_ynum_to_den * x_itv.upper_bound().numerator +
                                        coeff_.xden_ynum_to_den * x_itv.upper_bound().denominator,
                                    coeff_.xnum_yden_to_den * x_itv.upper_bound().numerator +
                                        coeff_.xden_yden_to_den * x_itv.upper_bound().denominator},
                                y_itv);
                            auto itv3 = map_cyclic_interval(
                                linear_fractional_transform{
                                    coeff_.xnum_ynum_to_num * y_itv.upper_bound().numerator +
                                        coeff_.xnum_yden_to_num * y_itv.upper_bound().denominator,
                                    coeff_.xden_ynum_to_num * y_itv.upper_bound().numerator +
                                        coeff_.xden_yden_to_num * y_itv.upper_bound().denominator,
                                    coeff_.xnum_ynum_to_den * y_itv.upper_bound().numerator +
                                        coeff_.xnum_yden_to_den * y_itv.upper_bound().denominator,
                                    coeff_.xden_ynum_to_den * y_itv.upper_bound().numerator +
                                        coeff_.xden_yden_to_den * y_itv.upper_bound().denominator},
                                x_itv);
                            auto itv4 = map_cyclic_interval(
                                linear_fractional_transform{
                                    coeff_.xnum_ynum_to_num * x_itv.lower_bound().numerator +
                                        coeff_.xden_ynum_to_num * x_itv.lower_bound().denominator,
                                    coeff_.xnum_yden_to_num * x_itv.lower_bound().numerator +
                                        coeff_.xden_yden_to_num * x_itv.lower_bound().denominator,
                                    coeff_.xnum_ynum_to_den * x_itv.lower_bound().numerator +
                                        coeff_.xden_ynum_to_den * x_itv.lower_bound().denominator,
                                    coeff_.xnum_yden_to_den * x_itv.lower_bound().numerator +
                                        coeff_.xden_yden_to_den * x_itv.lower_bound().denominator},
                                y_itv);

                            return check_union(
                                itv1, itv2, [&itv3, &itv4, &check_floor](auto&& joined_itv12) {
                                    return check_union(
                                        joined_itv12, itv3,
                                        [&itv4, &check_floor](auto&& joined_itv123) {
                                            return check_union(
                                                joined_itv123, itv4,
                                                [&check_floor](auto&& joined_itv1234) {
                                                    return joined_itv1234.visit(check_floor);
                                                });
                                        });
                                });
                        }
                    };

                    switch (xcf_.current_interval().visit([&](auto&& x_itv) -> final_result {
                        return ycf_.current_interval().visit([&](auto&& y_itv) -> final_result {
                            return compute_range(x_itv, y_itv);
                        });
                    })) {
                    case final_result::success:
                        return true;

                    case final_result::fail:
                        return false;

                    case final_result::retry:;
                    }

                    xcf_.update();
                    ycf_.update();
                }
            }

        public:
            struct default_mixin_initializer {
                static constexpr partial_fraction_type initial_partial_fraction() {
                    return {Unity{}, int_type{0}};
                }
                static constexpr interval_type initial_interval() {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }
            };

            template <class MixinInitializer = default_mixin_initializer>
            constexpr binary_gosper(XContinuedFractionImpl xcf, YContinuedFractionImpl ycf,
                                    bilinear_fractional_transform<int_type> coeff,
                                    MixinInitializer&& mixin_initializer = {})
                : crtp_base{mixin_initializer}, xcf_{static_cast<XContinuedFractionImpl>(xcf)},
                  ycf_{static_cast<XContinuedFractionImpl>(ycf)},
                  coeff_{static_cast<bilinear_fractional_transform<int_type>&&>(coeff)} {
                update_coeff_info();
            }
        };
    }
}

#endif
