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

#ifndef JKJ_HEADER_RATIONAL_CONTINUED_FRACTION
#define JKJ_HEADER_RATIONAL_CONTINUED_FRACTION

#include "continued_fraction.h"

namespace jkj {
    namespace cntfrc {
        namespace impl {
            template <class Int, class UInt, class Unity = unity>
            class rational {
            public:
                using partial_fraction_type = projective_rational<Unity, Int>;
                using convergent_type = projective_rational<Int, UInt>;
                using interval_type = variable_shape_cyclic_interval<
                    convergent_type, cyclic_interval_type_t::single_point,
                    cyclic_interval_type_t::left_open_right_closed,
                    cyclic_interval_type_t::left_closed_right_open, cyclic_interval_type_t::entire>;

            private:
                projective_rational<Int, UInt> fraction_;

            public:
                static constexpr auto initial_partial_fraction() {
                    return partial_fraction_type{Unity{unity{}}, Int{0}};
                }
                static constexpr auto initial_interval() noexcept {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }

                explicit constexpr rational(projective_rational<Int, UInt> r)
                    : fraction_{std::move(r)} {}

                explicit constexpr rational(Int numerator, UInt denominator)
                    : fraction_{std::move(numerator), std::move(denominator)} {}

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    if (!util::is_zero(fraction_.denominator)) {
                        auto div_result = util::div(fraction_.numerator, fraction_.denominator);
                        fraction_.numerator = Int{std::move(fraction_.denominator)};
                        fraction_.denominator = std::move(div_result.rem);

                        callback(partial_fraction_type{Unity{unity{}}, std::move(div_result.quot)});
                    }
                }
            };

            template <class Int, class UInt>
            rational(projective_rational<Int, UInt>) -> rational<Int, UInt, unity>;

            template <class Int, class UInt>
            rational(Int, UInt) -> rational<Int, UInt, unity>;
        }
    }
}

#endif
