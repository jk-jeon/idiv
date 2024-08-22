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

#ifndef JKJ_HEADER_IDIV_CONTINUED_FRACTION_ENGINE_GOSPER
#define JKJ_HEADER_IDIV_CONTINUED_FRACTION_ENGINE_GOSPER

#include "../extended_linear_fractional_mapping.h"
#include "../generator.h"

namespace jkj {
    namespace cntfrc {
        namespace engine {
            namespace detail {
                // If the interval itv gets mapped into a single integer by the floor function, call
                // update_and_callback with that integer and then return true. Otherwise, call
                // refine_interval and the return false.
                template <class VariableShapeInterval, class UpdateAndCallback,
                          class RefineInterval>
                constexpr bool gosper_core(VariableShapeInterval&& itv,
                                           UpdateAndCallback&& update_and_callback,
                                           RefineInterval&& refine_interval) {
                    return itv.visit([&update_and_callback, &refine_interval](auto const& itv_) {
                        using itv_type = std::remove_cvref_t<decltype(itv_)>;
                        constexpr auto infinity = projective_rational<unity, zero>{};

                        using enum cyclic_interval_type_t;
                        if constexpr (itv_type::interval_type() == single_point) {
                            // Stop if we have reached the infinity.
                            if (itv_.lower_bound() != infinity) {
                                update_and_callback(util::div_floor(
                                    itv_.lower_bound().numerator, itv_.lower_bound().denominator));
                            }
                            return true;
                        }
                        else if constexpr (itv_type::interval_type() == entire ||
                                           itv_type::interval_type() == single_complement) {
                            // We cannot do anything further, so refine the interval.
                            refine_interval();
                            return false;
                        }
                        else {
                            // Cannot do anything if the closure of the interval estimate contains
                            // the infinity.
                            if (itv_.lower_bound() == infinity || itv_.upper_bound() == infinity ||
                                cyclic_order(itv_.lower_bound(), infinity, itv_.upper_bound())) {
                                refine_interval();
                                return false;
                            }

                            auto floor_lower_bound = util::div_floor(
                                itv_.lower_bound().numerator, itv_.lower_bound().denominator);
                            auto compare_result = util::strong_order_to_int(
                                (floor_lower_bound + 1) * itv_.upper_bound().denominator <=>
                                itv_.upper_bound().numerator);
                            if (util::is_strictly_negative(itv_.upper_bound().denominator)) {
                                compare_result *= -1;
                            }

                            if (compare_result > 0 ||
                                (compare_result >= 0 &&
                                 itv_type::right_boundary_type() == boundary_type_t::open)) {
                                // floor(lower_bound) + 1 is not in the interval, so there is only
                                // one possible value of the floor.
                                update_and_callback(std::move(floor_lower_bound));
                                return true;
                            }
                            else {
                                // Otherwise, we cannot do anything further, so refine the interval.
                                refine_interval();
                                return false;
                            }
                        }
                    });
                }
            }

            template <class IntervalEstimateProvider>
            class unary_gosper {
            public:
                using decay_type = unary_gosper<
                    typename std::remove_cvref_t<IntervalEstimateProvider>::decay_type>;
                using convergent_type =
                    typename std::remove_cvref_t<decltype(std::declval<IntervalEstimateProvider>()
                                                              .current_interval())>::value_type;
                using int_type = decltype(convergent_type::numerator);
                using partial_fraction_type = projective_rational<unity, int_type>;
                using interval_type =
                    decltype(std::declval<extended_linear_fractional_mapping<int_type>>()
                                 .map_cyclic_interval(
                                     std::declval<IntervalEstimateProvider>().current_interval()));
                using internal_interval_estimate_provider = IntervalEstimateProvider;

            private:
                IntervalEstimateProvider itv_provider_;
                extended_linear_fractional_mapping<int_type> coeff_;
                interval_type complete_fraction_itv_;

                template <class>
                friend class unary_gosper;

            public:
                constexpr unary_gosper(IntervalEstimateProvider itv_provider,
                                       linear_fractional_mapping<int_type> coeff)
                    : itv_provider_{static_cast<IntervalEstimateProvider&&>(itv_provider)},
                      coeff_{std::move(coeff)}, complete_fraction_itv_{coeff_.map_cyclic_interval(
                                                    itv_provider_.current_interval())} {}

                template <class OtherIntervalEstimateProvider>
                    requires(std::is_same_v<unary_gosper, decay_type> &&
                             !std::is_same_v<IntervalEstimateProvider,
                                             OtherIntervalEstimateProvider> &&
                             std::is_constructible_v<IntervalEstimateProvider,
                                                     OtherIntervalEstimateProvider>)
                explicit constexpr unary_gosper(
                    unary_gosper<OtherIntervalEstimateProvider> const& other)
                    : itv_provider_{other.itv_provider_}, coeff_{other.coeff_},
                      complete_fraction_itv_{other.complete_fraction_itv_} {}

                // Copy/move constructor/assignements are allowed only for value-like generators.
                constexpr unary_gosper(unary_gosper const&)
                    requires(std::is_same_v<unary_gosper, decay_type>)
                = default;
                constexpr unary_gosper(unary_gosper&&)
                    requires(std::is_same_v<unary_gosper, decay_type>)
                = default;
                constexpr unary_gosper& operator=(unary_gosper const&)
                    requires(std::is_same_v<unary_gosper, decay_type>)
                = default;
                constexpr unary_gosper& operator=(unary_gosper&&)
                    requires(std::is_same_v<unary_gosper, decay_type>)
                = default;

                constexpr unary_gosper(unary_gosper const&)
                    requires(!std::is_same_v<unary_gosper, decay_type>)
                = delete;
                constexpr unary_gosper(unary_gosper&&)
                    requires(!std::is_same_v<unary_gosper, decay_type>)
                = delete;
                constexpr unary_gosper& operator=(unary_gosper const&)
                    requires(!std::is_same_v<unary_gosper, decay_type>)
                = delete;
                constexpr unary_gosper& operator=(unary_gosper&&)
                    requires(!std::is_same_v<unary_gosper, decay_type>)
                = delete;

                // Make a deep copy.
                constexpr decay_type copy() const { return decay_type{*this}; }

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    auto update_and_callback = [&](auto&& common_floor) {
                        coeff_.translate(-common_floor);
                        coeff_.reflect();
                        complete_fraction_itv_ =
                            extended_linear_fractional_mapping<int_type>{zero{}, unity{}, unity{},
                                                                         -common_floor}
                                .map_cyclic_interval(complete_fraction_itv_);
                        callback.on_next_partial_fraction(partial_fraction_type{
                            unity{}, static_cast<decltype(common_floor)&&>(common_floor)});
                    };
                    auto refine_interval = [this]() {
                        itv_provider_.refine_interval();
                        complete_fraction_itv_ =
                            coeff_.map_cyclic_interval(itv_provider_.current_interval());
                    };
                    while (!detail::gosper_core(complete_fraction_itv_, update_and_callback,
                                                refine_interval)) {
                    }
                }

                constexpr interval_type initial_interval() const {
                    return coeff_.map_cyclic_interval(itv_provider_.current_interval());
                }
            };

            template <class IntervalEstimateProvider>
            unary_gosper(IntervalEstimateProvider,
                         linear_fractional_mapping<
                             decltype(std::remove_cvref_t<
                                      decltype(std::declval<IntervalEstimateProvider>()
                                                   .current_interval())>::value_type::numerator)>)
                -> unary_gosper<IntervalEstimateProvider>;


            template <class IntervalEstimateProviderX, class IntervalEstimateProviderY>
            class binary_gosper {
            public:
                static_assert(std::is_same_v<typename std::remove_cvref_t<
                                                 decltype(std::declval<IntervalEstimateProviderX>()
                                                              .current_interval())>::value_type,
                                             typename std::remove_cvref_t<
                                                 decltype(std::declval<IntervalEstimateProviderY>()
                                                              .current_interval())>::value_type>,
                              "the internal interval estimate providers for binary_gosper must "
                              "have the same value types");

                using decay_type = binary_gosper<
                    typename std::remove_cvref_t<IntervalEstimateProviderX>::decay_type,
                    typename std::remove_cvref_t<IntervalEstimateProviderY>::decay_type>;
                using convergent_type =
                    typename std::remove_cvref_t<decltype(std::declval<IntervalEstimateProviderX>()
                                                              .current_interval())>::value_type;
                using int_type = decltype(convergent_type::numerator);
                using partial_fraction_type = projective_rational<unity, int_type>;
                using interval_type =
                    decltype(std::declval<extended_bilinear_fractional_mapping<int_type>>()
                                 .template map_cyclic_rectangle<true>(
                                     std::declval<IntervalEstimateProviderX>().current_interval(),
                                     std::declval<IntervalEstimateProviderY>().current_interval()));
                using first_internal_interval_estimate_provider = IntervalEstimateProviderX;
                using second_internal_interval_estimate_provider = IntervalEstimateProviderY;

            private:
                IntervalEstimateProviderX itv_provider_x_;
                IntervalEstimateProviderY itv_provider_y_;
                extended_bilinear_fractional_mapping<int_type> coeff_;
                interval_type complete_fraction_itv_;

                template <class, class>
                friend class binary_gosper;

                constexpr interval_type compute_initial_interval_for_complete_fraction() {
                    // Get away from the indeterminacy locus.
                    // This loop will never finish if (x,y) is located inside the indeterminacy
                    // locus.
                    while (coeff_.touches_indeterminacy_locus(itv_provider_x_.current_interval(),
                                                              itv_provider_y_.current_interval())) {
                        itv_provider_x_.refine_interval();
                        itv_provider_y_.refine_interval();
                    }
                    return coeff_.template map_cyclic_rectangle<true>(
                        itv_provider_x_.current_interval(), itv_provider_y_.current_interval());
                }

            public:
                constexpr binary_gosper(IntervalEstimateProviderX itv_provider_x,
                                        IntervalEstimateProviderY itv_provider_y,
                                        bilinear_fractional_mapping<int_type> coeff)
                    : itv_provider_x_{static_cast<IntervalEstimateProviderX&&>(itv_provider_x)},
                      itv_provider_y_{static_cast<IntervalEstimateProviderY&&>(itv_provider_y)},
                      coeff_{std::move(coeff)},
                      complete_fraction_itv_{compute_initial_interval_for_complete_fraction()} {}

                template <class OtherIntervalEstimateProviderX,
                          class OtherIntervalEstimateProviderY>
                    requires(std::is_same_v<binary_gosper, decay_type> &&
                             !(std::is_same_v<IntervalEstimateProviderX,
                                              OtherIntervalEstimateProviderX> &&
                               !std::is_same_v<IntervalEstimateProviderY,
                                               OtherIntervalEstimateProviderY>) &&
                             std::is_constructible_v<IntervalEstimateProviderX,
                                                     OtherIntervalEstimateProviderX> &&
                             std::is_constructible_v<IntervalEstimateProviderY,
                                                     OtherIntervalEstimateProviderY>)
                explicit constexpr binary_gosper(
                    binary_gosper<OtherIntervalEstimateProviderX,
                                  OtherIntervalEstimateProviderY> const& other)
                    : itv_provider_x_{other.itv_provider_x_},
                      itv_provider_y_{other.itv_provider_y_}, coeff_{other.coeff_},
                      complete_fraction_itv_{other.complete_fraction_itv_} {}

                // Make a deep copy.
                constexpr decay_type copy() const { return decay_type{*this}; }

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    auto update_and_callback = [&](auto&& common_floor) {
                        coeff_.translate(-common_floor);
                        coeff_.reflect();
                        complete_fraction_itv_ =
                            extended_linear_fractional_mapping<int_type>{zero{}, unity{}, unity{},
                                                                         -common_floor}
                                .map_cyclic_interval(complete_fraction_itv_);
                        callback.on_next_partial_fraction(partial_fraction_type{
                            unity{}, static_cast<decltype(common_floor)&&>(common_floor)});
                    };
                    auto refine_interval = [this]() {
                        itv_provider_x_.refine_interval();
                        itv_provider_y_.refine_interval();
                        complete_fraction_itv_ = coeff_.template map_cyclic_rectangle<true>(
                            itv_provider_x_.current_interval(), itv_provider_y_.current_interval());
                    };
                    while (!detail::gosper_core(complete_fraction_itv_, update_and_callback,
                                                refine_interval)) {
                    }
                }

                constexpr interval_type initial_interval() const {
                    return coeff_.template map_cyclic_rectangle<true>(
                        itv_provider_x_.current_interval(), itv_provider_y_.current_interval());
                }
            };

            template <class IntervalEstimateProviderX, class IntervalEstimateProviderY>
            binary_gosper(IntervalEstimateProviderX, IntervalEstimateProviderY,
                          bilinear_fractional_mapping<
                              decltype(std::remove_cvref_t<
                                       decltype(std::declval<IntervalEstimateProviderX>()
                                                    .current_interval())>::value_type::numerator)>)
                -> binary_gosper<IntervalEstimateProviderX, IntervalEstimateProviderY>;
        }

        // Helper functions.
        template <class... Mixins, class ContinuedFractionEngine>
        constexpr auto make_unary_gosper_generator(
            ContinuedFractionEngine engine,
            linear_fractional_mapping<decltype(ContinuedFractionEngine::convergent_type::numerator)>
                coeff) {
            return make_generator<Mixins...>(
                engine::unary_gosper{make_generator<interval_estimate_provider>(
                                         static_cast<ContinuedFractionEngine&&>(engine)),
                                     std::move(coeff)});
        }
        template <class... Mixins, class ContinuedFractionEngineX, class ContinuedFractionEngineY>
        constexpr auto make_binary_gosper_generator(
            ContinuedFractionEngineX engine_x, ContinuedFractionEngineY engine_y,
            bilinear_fractional_mapping<
                decltype(ContinuedFractionEngineX::convergent_type::numerator)>
                coeff) {
            return make_generator<Mixins...>(
                engine::binary_gosper{make_generator<interval_estimate_provider>(
                                          static_cast<ContinuedFractionEngineX&&>(engine_x)),
                                      make_generator<interval_estimate_provider>(
                                          static_cast<ContinuedFractionEngineY&&>(engine_y)),
                                      std::move(coeff)});
        }
    }
}

#endif
