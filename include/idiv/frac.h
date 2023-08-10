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

#ifndef JKJ_HEADER_FRAC
#define JKJ_HEADER_FRAC

#include <compare>

namespace jkj {
    template <class Num, class Den>
    struct frac {
        Num numerator;
        Den denominator;
    };
    template <class Num, class Den>
    frac(Num, Den) -> frac<Num, Den>;


    template <class Num1, class Den1, class Num2, class Den2>
    constexpr bool operator==(frac<Num1, Den1> const& x, frac<Num2, Den2> const& y) {
        return x.numerator * y.denominator == y.numerator * x.denominator;
    }

    template <class Num1, class Den1, class Num2, class Den2>
    constexpr std::strong_ordering operator<=>(frac<Num1, Den1> const& x, frac<Num2, Den2> const& y) {
        return x.numerator * y.denominator <=> y.numerator * x.denominator;
    }

    // Performs no reduction.
    template <class Num1, class Den1, class Num2, class Den2>
    constexpr auto operator+(frac<Num1, Den1> const& x, frac<Num2, Den2> const& y) {
        return frac{x.numerator * y.denominator + y.numerator * x.denominator,
                        x.denominator * y.denominator};
    }

    // Performs no reduction.
    template <class Num1, class Den1, class Num2, class Den2>
    constexpr auto operator-(frac<Num1, Den1> const& x, frac<Num2, Den2> const& y) {
        return frac{x.numerator * y.denominator - y.numerator * x.denominator,
                        x.denominator * y.denominator};
    }

    // Performs no reduction.
    template <class Num1, class Den1, class Num2, class Den2>
    constexpr auto operator*(frac<Num1, Den1> const& x, frac<Num2, Den2> const& y) {
        return frac{x.numerator * y.numerator, x.denominator * y.denominator};
    }

    // Performs no reduction.
    template <class Num1, class Den1, class Num2, class Den2>
    constexpr auto operator/(frac<Num1, Den1> const& x, frac<Num2, Den2> const& y) {
        return frac{x.numerator * y.denominator, x.denominator * y.numerator};
    }
}

#endif