// Copyright 2022-2024 Junekey Jeon
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

#ifndef JKJ_HEADER_IDIV_CONTINUED_FRACTION_ENGINE_RATIONAL
#define JKJ_HEADER_IDIV_CONTINUED_FRACTION_ENGINE_RATIONAL

#include "../generator.h"

namespace jkj {
    namespace cntfrc {
        namespace engine {
            template <class Int, class UInt>
            class rational {
            public:
                using partial_fraction_type = projective_rational<unity, Int>;
                using convergent_type = projective_rational<Int, UInt>;
                using interval_type =
                    variable_shape_cyclic_interval<convergent_type,
                                                   cyclic_interval_type_t::single_point>;

            private:
                projective_rational<Int, UInt> fraction_;

            public:
                explicit constexpr rational(projective_rational<Int, UInt> r)
                    : fraction_{std::move(r)} {}

                explicit constexpr rational(Int numerator, UInt denominator)
                    : fraction_{std::move(numerator), std::move(denominator)} {}

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    if (!util::is_zero(fraction_.denominator)) {
                        auto div_result = util::div(fraction_.numerator, fraction_.denominator);
                        fraction_.numerator = static_cast<Int>(std::move(fraction_.denominator));
                        fraction_.denominator = std::move(div_result.rem);

                        callback.on_next_partial_fraction(
                            partial_fraction_type{unity{}, std::move(div_result.quot)});
                    }
                }

                constexpr auto initial_interval() const {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::single_point>{
                        fraction_};
                }
            };

            template <class Int, class UInt>
            rational(projective_rational<Int, UInt>) -> rational<Int, UInt>;

            template <class Int, class UInt>
            rational(Int, UInt) -> rational<Int, UInt>;
        }

        template <class Int, class UInt>
        struct mixin_traits<interval_estimate_provider, engine::rational<Int, UInt>> {
            using interval_type = typename engine::rational<Int, UInt>::interval_type;

            using tracking_data =
                interval_estimate_provider::default_tracking_data<engine::rational<Int, UInt>>;

            class facade : util::noncopyable<>, util::nonmovable<> {
            public:
                template <class Engine>
                explicit constexpr facade(Engine&) noexcept {}

                template <class Advancer>
                constexpr void refine_interval(Advancer&&) {}
            };
        };

        template <class... Mixins, class Int, class UInt>
        constexpr auto make_rational_generator(Int&& numerator, UInt&& denominator) {
            return make_generator<Mixins...>(
                engine::rational{static_cast<Int&&>(numerator), static_cast<UInt&&>(denominator)});
        }

        template <class... Mixins, class Int, class UInt>
        constexpr auto make_rational_generator(projective_rational<Int, UInt> r) {
            return make_generator<Mixins...>(engine::rational{std::move(r)});
        }

        template <class Numerator, class Denominator>
        constexpr auto reduce_fraction(Numerator&& num, Denominator&& den) {
            using signed_denominator_type =
                decltype(util::to_signed(static_cast<Denominator&&>(den)));
            using denominator_type = decltype(util::abs(static_cast<Denominator&&>(den)));
            using signed_numerator_type = decltype(util::to_signed(static_cast<Numerator&&>(num)));
            using numerator_type = std::conditional_t<
                std::is_same_v<signed_numerator_type, std::remove_cvref_t<Numerator>> ||
                    std::is_same_v<signed_denominator_type, std::remove_cvref_t<Denominator>>,
                signed_numerator_type, std::remove_cvref_t<Numerator>>;

            auto cf = [&num, &den] {
                if constexpr (std::is_same_v<denominator_type, std::remove_cvref_t<Denominator>>) {
                    return make_generator<convergent_tracker>(
                        engine::rational<numerator_type, denominator_type>{
                            numerator_type{static_cast<Numerator&&>(num)},
                            denominator_type{static_cast<Denominator&&>(den)}});
                }
                else {
                    if (util::is_strictly_negative(den)) {
                        return make_generator<convergent_tracker>(
                            engine::rational<numerator_type, denominator_type>{
                                -static_cast<numerator_type>(static_cast<Numerator&&>(num)),
                                util::abs(static_cast<Denominator&&>(den))});
                    }
                    else {
                        return make_generator<convergent_tracker>(
                            engine::rational<numerator_type, denominator_type>{
                                static_cast<numerator_type>(static_cast<Numerator&&>(num)),
                                util::abs(static_cast<Denominator&&>(den))});
                    }
                }
            }();

            while (!cf.terminated()) {
                cf.proceed_to_next_partial_fraction();
            }
            return cf.current_convergent();
        }

        template <class ProjectiveRational>
        constexpr auto reduce_fraction(ProjectiveRational&& x) {
            static_assert(
                tmp::is_specialization_of<std::remove_cvref_t<ProjectiveRational>,
                                          projective_rational>(),
                "single parameter overload of cntfrc::reduce_fraction can be called only for "
                "arguments of type projective_rational");

            return reduce_fraction(static_cast<ProjectiveRational&&>(x).numerator,
                                   static_cast<ProjectiveRational&&>(x).denominator);
        }
    }
}

#endif
