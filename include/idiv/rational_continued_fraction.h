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
#include "interval.h"
#include <cstdlib>

namespace jkj {
    namespace cntfrc {
        template <class Int, class UInt, class Unity = unity, template <class> class... Mixins>
        class rational_continued_fraction;

        template <class Int, class UInt, class Unity, template <class> class... Mixins>
        struct continued_fraction_traits<rational_continued_fraction<Int, UInt, Unity, Mixins...>> {
            using partial_fraction_type = frac<Unity, Int>;
            using convergent_type = projective_rational<Int, UInt>;
            using interval_type = variable_shape_cyclic_interval<
                convergent_type, cyclic_interval_type_t::single_point,
                cyclic_interval_type_t::left_open_right_closed,
                cyclic_interval_type_t::left_closed_right_open, cyclic_interval_type_t::entire>;
        };

        template <class Int, class UInt, class Unity, template <class> class... Mixins>
        class rational_continued_fraction
            : public continued_fraction_base<
                  rational_continued_fraction<Int, UInt, Unity, Mixins...>, Mixins...> {
            using crtp_base =
                continued_fraction_base<rational_continued_fraction<Int, UInt, Unity, Mixins...>,
                                        Mixins...>;
            friend crtp_base;

        public:
            using partial_fraction_type = typename crtp_base::traits_type::partial_fraction_type;
            using convergent_type = typename crtp_base::traits_type::convergent_type;
            using interval_type = typename crtp_base::traits_type::interval_type;

        private:
            projective_rational<Int, UInt> fraction_;

            template <class Functor>
            constexpr void with_next_partial_fraction(Functor&& f) {
                if (!is_zero(fraction_.denominator)) {
                    using std::div;
                    auto div_result = div(fraction_.numerator, fraction_.denominator);
                    fraction_.numerator = Int{static_cast<UInt&&>(fraction_.denominator)};
                    fraction_.denominator = static_cast<UInt&&>(div_result.rem);

                    f(partial_fraction_type{Unity{}, static_cast<Int&&>(div_result.quot)});
                }
            }

        public:
            struct default_mixin_initializer {
                static constexpr partial_fraction_type initial_partial_fraction() {
                    return {Unity{}, Int{0}};
                }
                static constexpr interval_type initial_interval() noexcept {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }
            };

            template <class MixinInitializer = default_mixin_initializer>
            explicit constexpr rational_continued_fraction(
                projective_rational<Int, UInt> r, MixinInitializer&& mixin_initializer = {})
                : crtp_base{mixin_initializer},
                  fraction_{static_cast<Int&&>(r.numerator), static_cast<UInt&&>(r.denominator)} {}
        };
    }
}

#endif
