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

#ifndef JKJ_HEADER_PRIME_FACTORIZATION
#define JKJ_HEADER_PRIME_FACTORIZATION

#include "frac.h"
#include "util.h"
#include <cstdlib>
#include <vector>

namespace jkj {
    template <class UInt>
    struct prime_factor {
        UInt factor;
        int exponent;

        bool operator==(prime_factor const&) const = default;
    };

    template <class UInt>
    constexpr std::vector<prime_factor<UInt>> prime_factorization(frac<UInt, UInt> n) {
        // I know this is slow, but who cares.
        std::vector<prime_factor<UInt>> prime_factors;

        // Factor out 2 first.
        {
            auto exponent = int(factor_out_power_of_2(n.numerator));
            exponent -= int(factor_out_power_of_2(n.denominator));
            if (exponent != 0) {
                prime_factors.push_back({UInt(2u), exponent});
            }
        }

        UInt divisor{3u};
        while (true) {
            using std::div;
            using util::max;

            auto const& larger_one = max(n.numerator, n.denominator);
            if (divisor * divisor > larger_one) {
                break;
            }

            int exponent = 0;
            while (true) {
                auto div_result = div(n.numerator, divisor);
                if (is_zero(div_result.rem)) {
                    ++exponent;
                    n.numerator = std::move(div_result.quot);
                    continue;
                }
                break;
            }
            while (true) {
                auto div_result = div(n.denominator, divisor);
                if (is_zero(div_result.rem)) {
                    --exponent;
                    n.denominator = std::move(div_result.quot);
                    continue;
                }
                break;
            }

            if (exponent != 0) {
                prime_factors.push_back({divisor, exponent});
            }

            divisor += 2;
        }

        if (n.numerator != 1) {
            prime_factors.push_back({std::move(n.numerator), 1});
        }
        else if (n.denominator != 1) {
            prime_factors.push_back({std::move(n.denominator), -1});
        }

        return prime_factors;
    }
}

#endif