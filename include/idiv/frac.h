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
#include <type_traits>

namespace jkj {
    // Num: supposed to be jkj::bigint::uint_var/uint_const_t/int_var/int_const_t.
    // Den: supposed to be jkj::bigint::uint_var/uint_const_t.
    template <class Num, class Den>
    struct frac {
        Num numerator;
        Den denominator;
    };
    template <class Num, class Den>
    frac(Num&&, Den&&) -> frac<std::remove_cvref_t<Num>, std::remove_cvref_t<Den>>;

    template <class Num1, class Den1, class Num2, class Den2>
    constexpr bool operator==(frac<Num1, Den1> const& x, frac<Num2, Den2> const& y) {
        return x.numerator * y.denominator == y.numerator * x.denominator;
    }

    // Relying on the assumption that the denominator is positive.
    template <class Num1, class Den1, class Num2, class Den2>
    constexpr std::strong_ordering operator<=>(frac<Num1, Den1> const& x,
                                               frac<Num2, Den2> const& y) {
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
    constexpr auto& operator*=(frac<Num1, Den1>& x, frac<Num2, Den2> const& y) {
        x.numerator *= y.numerator;
        x.denominator *= y.denominator;
        return x;
    }

    // Performs no reduction.
    template <class Num1, class Den1, class Num2, class Den2>
    constexpr auto operator/(frac<Num1, Den1> const& x, frac<Num2, Den2> const& y) {
        auto num = x.numerator * y.denominator;
        auto den = x.denominator * y.numerator;
        if (is_strictly_negative(den)) {
            return frac{invert_sign(static_cast<decltype(num)&&>(num)),
                        abs(static_cast<decltype(den)&&>(den))};
        }
        else {
            return frac{static_cast<decltype(num)&&>(num), abs(static_cast<decltype(den)&&>(den))};
        }
    }
}

#endif
