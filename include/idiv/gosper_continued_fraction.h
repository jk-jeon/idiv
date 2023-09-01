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
            using interval_type =
                variable_shape_cyclic_interval<projective_rational<int_type, int_type>>;
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
                determinant_sign_ = [&] {
                    auto const result = coeff_.num_to_num * coeff_.den_to_den <=>
                                        coeff_.den_to_num * coeff_.num_to_den;

                    if (result < 0) {
                        return -1;
                    }
                    if (result > 0) {
                        return 1;
                    }
                    else {
                        return 0;
                    }
                }();
            }

            template <class Functor>
            constexpr bool with_next_partial_fraction(Functor&& f) {
                while (true) {
                    // Step 1. Get away from singularities.
                    if (determinant_sign_ == 0) {
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
                                {coeff_.num_to_num, coeff_.num_to_den}};
                        }

                        return cf_.current_interval().visit([&](auto&& itv) -> interval_type {
                            using enum cyclic_interval_type_t;
                            if constexpr (itv.interval_type() == empty) {
                                util::constexpr_assert(false);
                                return cyclic_interval<projective_rational<int_type, int_type>,
                                                       empty>{};
                            }
                            else if constexpr (itv.interval_type() == single_point) {
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
        };
    }

#if 0
    template <class Int>
    struct binary_gosper_coeff {
        struct coeff {
            Int const_coeff;
            Int x_coeff;
            Int y_coeff;
            Int xy_coeff;

            template <class PartialFraction>
            constexpr void feed_x(PartialFraction const& partial_fraction) {
                // s/t, (a, b, c, d) |-> (b, sa + tb, d, sc + td)
                const_coeff *= partial_fraction.numerator;
                const_coeff += partial_fraction.denominator * x_coeff;
                y_coeff *= partial_fraction.numerator;
                y_coeff += partial_fraction.denominator * xy_coeff;
                swap(const_coeff, x_coeff);
                swap(y_coeff, xy_coeff);
            }

            template <class PartialFraction>
            constexpr void feed_y(PartialFraction const& partial_fraction) {
                // s/t, (a, b, c, d) |-> (c, d, sa + tc, sb + td)
                const_coeff *= partial_fraction.numerator;
                const_coeff += partial_fraction.denominator * y_coeff;
                x_coeff *= partial_fraction.numerator;
                x_coeff += partial_fraction.denominator * xy_coeff;
                swap(const_coeff, y_coeff);
                swap(x_coeff, xy_coeff);
            }
        };

        coeff numerator;
        coeff denominator;
    };

    template <class ContinuedFractionImplX, class ContinuedFractionImplY, class Unity = unity>
        requires(std::is_same_v<typename ContinuedFractionImplX::partial_fraction_type,
                                typename ContinuedFractionImplY::partial_fraction_type> &&
                 std::is_same_v<typename ContinuedFractionImplX::convergent_type,
                                typename ContinuedFractionImplY::convergent_type>)
    class binary_gosper_continued_fraction {
    public:
        using int_type = decltype(ContinuedFractionImplX::convergent_type::numerator);
        using uint_type = decltype(ContinuedFractionImplX::convergent_type::denominator);
        using partial_fraction_type = frac<Unity, int_type>;
        using convergent_type = typename ContinuedFractionImplX::convergent_type;

    private:
        ContinuedFractionImplX x_cf_;
        ContinuedFractionImplY y_cf_;
        binary_gosper_coeff<int_type> coeff_;
        bool is_x_terminated_ = false;
        bool is_y_terminated_ = false;

        constexpr void progress_x() {
            auto result = x_cf_.next_partial_fraction();
            coeff_.numerator.feed_x(result.partial_fraction);
            coeff_.denominator.feed_x(result.partial_fraction);

            // If the current coefficient is the last one, we have to use a different update
            // formula from now on.
            if (result.is_last) {
                is_x_terminated_ = true;
                using util::swap;
                swap(coeff_.numerator.const_coeff, coeff_.numerator.x_coeff);
                swap(coeff_.numerator.y_coeff, coeff_.numerator.xy_coeff);
                swap(coeff_.denominator.const_coeff, coeff_.denominator.x_coeff);
                swap(coeff_.denominator.y_coeff, coeff_.denominator.xy_coeff);
            }
        }
        constexpr void progress_y() {
            auto result = y_cf_.next_partial_fraction();
            coeff_.numerator.feed_y(result.partial_fraction);
            coeff_.denominator.feed_y(result.partial_fraction);

            // If the current coefficient is the last one, we have to use a different update
            // formula from now on.
            if (result.is_last) {
                is_y_terminated_ = true;
                using util::swap;
                swap(coeff_.numerator.const_coeff, coeff_.numerator.y_coeff);
                swap(coeff_.numerator.x_coeff, coeff_.numerator.xy_coeff);
                swap(coeff_.denominator.const_coeff, coeff_.denominator.y_coeff);
                swap(coeff_.denominator.x_coeff, coeff_.denominator.xy_coeff);
            }
        }

    public:
        constexpr binary_gosper_continued_fraction(ContinuedFractionImplX x_cf,
                                                   ContinuedFractionImplY y_cf,
                                                   binary_gosper_coeff<int_type> coeff)
            : x_cf_{static_cast<ContinuedFractionImplX>(x_cf)},
              y_cf_{static_cast<ContinuedFractionImplY>(y_cf)},
              coeff_{static_cast<binary_gosper_coeff<int_type>&&>(coeff)} {}

        constexpr next_partial_fraction_return<partial_fraction_type> next_partial_fraction() {
            while (true) {
                if (is_x_terminated_ && is_y_terminated_) {
                    // Proceed as in the case of usual rational continued fractions.
                    util::constexpr_assert<util::error_msgs::divide_by_zero>(
                        !is_zero(coeff_.denominator.const_coeff));
                    util::constexpr_assert<util::error_msgs::divide_by_zero>(
                        is_strictly_positive(coeff_.denominator.const_coeff));

                    using std::div;
                    auto div_result =
                        div(coeff_.numerator.const_coeff, abs(coeff_.denominator.const_coeff));
                    coeff_.numerator.const_coeff =
                        static_cast<int_type&&>(coeff_.denominator.const_coeff);
                    coeff_.denominator.const_coeff =
                        int_type(static_cast<uint_type&&>(div_result.rem));

                    // Terminate if all coefficients in the denominator has become zero.
                    return {partial_fraction_type{{}, static_cast<int_type&&>(div_result.quot)},
                            is_zero(coeff_.denominator.const_coeff)};
                }
                else if (is_x_terminated_) {
                    if (!is_zero(coeff_.denominator.const_coeff) &&
                        !is_zero(coeff_.denominator.y_coeff)) {
                        auto const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto y_output =
                            div_floor(coeff_.numerator.y_coeff, coeff_.denominator.y_coeff);

                        if (const_output == y_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.y_coeff -= const_output * coeff_.denominator.y_coeff;

                            using util::swap;
                            swap(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);

                            // Terminate if all coefficients in the denominator has become zero.
                            return {
                                partial_fraction_type{{}, static_cast<int_type&&>(const_output)},
                                is_zero(coeff_.denominator.const_coeff) &&
                                    is_zero(coeff_.denominator.y_coeff)};
                        }
                    }

                    // If two endpoints do not agree, then refine the region.
                    progress_y();
                }
                else if (is_y_terminated_) {
                    if (!is_zero(coeff_.denominator.const_coeff) &&
                        !is_zero(coeff_.denominator.x_coeff)) {
                        auto const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto x_output =
                            div_floor(coeff_.numerator.x_coeff, coeff_.denominator.x_coeff);

                        if (const_output == x_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.x_coeff -= const_output * coeff_.denominator.x_coeff;

                            using util::swap;
                            swap(coeff_.numerator, coeff_.denominator);

                            // Terminate if all coefficients in the denominator has become zero.
                            return {
                                partial_fraction_type{{}, static_cast<int_type&&>(const_output)},
                                is_zero(coeff_.denominator.const_coeff) &&
                                    is_zero(coeff_.denominator.x_coeff)};
                        }
                    }

                    // If two endpoints do not agree, then refine the region.
                    progress_x();
                }
                else {
                    if (!is_zero(coeff_.denominator.const_coeff) &&
                        !is_zero(coeff_.denominator.x_coeff) &&
                        !is_zero(coeff_.denominator.y_coeff) &&
                        !is_zero(coeff_.denominator.xy_coeff)) {
                        auto const_output =
                            div_floor(coeff_.numerator.const_coeff, coeff_.denominator.const_coeff);
                        auto x_output =
                            div_floor(coeff_.numerator.x_coeff, coeff_.denominator.x_coeff);
                        auto y_output =
                            div_floor(coeff_.numerator.y_coeff, coeff_.denominator.y_coeff);
                        auto xy_output =
                            div_floor(coeff_.numerator.xy_coeff, coeff_.denominator.xy_coeff);

                        if (const_output == x_output && x_output == y_output &&
                            y_output == xy_output) {
                            // Found the new coefficient.
                            coeff_.numerator.const_coeff -=
                                const_output * coeff_.denominator.const_coeff;
                            coeff_.numerator.x_coeff -= const_output * coeff_.denominator.x_coeff;
                            coeff_.numerator.y_coeff -= const_output * coeff_.denominator.y_coeff;
                            coeff_.numerator.xy_coeff -= const_output * coeff_.denominator.xy_coeff;

                            using util::swap;
                            swap(coeff_.numerator, coeff_.denominator);

                            // Terminate if all coefficients in the denominator has become zero.
                            return {
                                partial_fraction_type{{}, static_cast<int_type&&>(const_output)},
                                is_zero(coeff_.denominator.const_coeff) &&
                                    is_zero(coeff_.denominator.x_coeff) &&
                                    is_zero(coeff_.denominator.y_coeff) &&
                                    is_zero(coeff_.denominator.xy_coeff)};
                        }
                    }

                    // If four endpoints do not agree, then refine the region.
                    // Compare |b/f - a/e| and |c/g - a/e|, and choose x if the former is larger.
                    // I.e., we choose x if |bge - afg| >= |cfe - afg|.

                    // If e & f are both zero, progressing x will never make them nonzero, so
                    // progress y in this case.
                    if (is_zero(coeff_.denominator.const_coeff) &&
                        is_zero(coeff_.denominator.x_coeff)) {
                        progress_y();
                    }
                    // If e & g are both zero, progressing y will never make them nonzero, so
                    // progress x in this case.
                    else if (is_zero(coeff_.denominator.const_coeff) &&
                             is_zero(coeff_.denominator.y_coeff)) {
                        progress_x();
                    }
                    else {
                        // afg
                        auto const base_point = coeff_.numerator.const_coeff *
                                                coeff_.denominator.x_coeff *
                                                coeff_.denominator.y_coeff;
                        // |bge - afg|
                        auto const const_to_x_length =
                            abs(coeff_.numerator.x_coeff * coeff_.denominator.y_coeff *
                                    coeff_.denominator.const_coeff -
                                base_point);
                        // |bge - afg|
                        auto const const_to_y_length =
                            abs(coeff_.numerator.y_coeff * coeff_.denominator.x_coeff *
                                    coeff_.denominator.const_coeff -
                                base_point);

                        if (const_to_x_length >= const_to_y_length) {
                            progress_x();
                        }
                        else {
                            progress_y();
                        }
                    }
                }
            }
        }
    };
#endif
}

#endif
