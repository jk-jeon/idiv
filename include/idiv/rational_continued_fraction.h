// Copyright 2022 Junekey Jeon
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
#include "frac.h"
#include <cstdlib>

namespace jkj {
    namespace cntfrc {
        namespace impl {
            template <class Int, class UInt, class Unity = unity>
            class rational {
            public:
                using partial_fraction_type = frac<Unity, Int>;
                using convergent_type = projective_rational<Int, UInt>;
                using interval_type = variable_shape_cyclic_interval<
                    convergent_type, cyclic_interval_type_t::single_point,
                    cyclic_interval_type_t::left_open_right_closed,
                    cyclic_interval_type_t::left_closed_right_open, cyclic_interval_type_t::entire>;

            private:
                projective_rational<Int, UInt> fraction_;

            public:
                static constexpr partial_fraction_type initial_partial_fraction() {
                    return {Unity{}, Int{0}};
                }
                static constexpr interval_type initial_interval() noexcept {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }

                explicit constexpr rational(projective_rational<Int, UInt> r)
                    : fraction_{static_cast<Int&&>(r.numerator),
                                static_cast<UInt&&>(r.denominator)} {}

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    if (!is_zero(fraction_.denominator)) {
                        using std::div;
                        auto div_result = div(fraction_.numerator, fraction_.denominator);
                        fraction_.numerator = Int{static_cast<UInt&&>(fraction_.denominator)};
                        fraction_.denominator = static_cast<UInt&&>(div_result.rem);

                        callback(
                            partial_fraction_type{Unity{}, static_cast<Int&&>(div_result.quot)});
                    }
                }
            };
        }
    }
}

#endif
