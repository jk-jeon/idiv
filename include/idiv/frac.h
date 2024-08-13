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

#ifndef JKJ_HEADER_IDIV_FRAC
#define JKJ_HEADER_IDIV_FRAC

#include "util.h"

namespace jkj {
    template <class Num, class Den>
    struct frac {
        Num numerator;
        Den denominator;

        frac() = default;

        template <class OtherNum = Num, class OtherDen = Den>
            requires requires(OtherNum num, OtherDen den) {
                Num{num};
                Den{den};
            }
        explicit constexpr frac(OtherNum&& num, OtherDen&& den)
            : numerator{static_cast<OtherNum&&>(num)}, denominator{static_cast<OtherDen&&>(den)} {}

        template <class OtherNum = Num>
            requires requires(OtherNum num) {
                Num{num};
                Den{1u};
            }
        explicit constexpr frac(OtherNum&& num)
            : numerator{static_cast<OtherNum&&>(num)}, denominator{1u} {}

        template <class OtherNum, class OtherDen>
            requires(requires(OtherNum num, OtherDen den) {
                        Num{num};
                        Den{den};
                    } && !(std::is_same_v<Num, OtherNum> && std::is_same_v<Den, OtherDen>))
        explicit constexpr frac(frac<OtherNum, OtherDen> const& other)
            : numerator{other.numerator}, denominator{other.denominator} {}

        template <class OtherNum, class OtherDen>
            requires(requires(OtherNum num, OtherDen den) {
                        Num{num};
                        Den{den};
                    } && !(std::is_same_v<Num, OtherNum> && std::is_same_v<Den, OtherDen>))
        explicit constexpr frac(frac<OtherNum, OtherDen>&& other)
            : numerator{static_cast<OtherNum&&>(other.numerator)},
              denominator{static_cast<OtherDen&&>(other.denominator)} {}
    };
    template <class Num, class Den>
    frac(Num&&, Den&&) -> frac<std::remove_cvref_t<Num>, std::remove_cvref_t<Den>>;
    template <class Num>
    frac(Num&&) -> frac<std::remove_cvref_t<Num>,
                        std::remove_cvref_t<decltype(util::abs(std::declval<Num>()))>>;

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

    template <class Num, class Den>
    constexpr auto operator-(frac<Num, Den> const& x) {
        return frac{-x.numerator, x.denominator};
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
        if (util::is_strictly_negative(den)) {
            return frac{util::invert_sign(std::move(num)), util::abs(std::move(den))};
        }
        else {
            return frac{std::move(num), util::abs(std::move(den))};
        }
    }

    // Performs no reduction.
    template <class Num, class Den>
    constexpr frac<std::remove_cvref_t<Num>, std::remove_cvref_t<Den>>
    make_frac_from_signed(Num&& numerator, Den&& denominator) {
        using return_type = frac<std::remove_cvref_t<Num>, std::remove_cvref_t<Den>>;
        if (util::is_strictly_negative(denominator)) {
            return return_type{util::invert_sign(static_cast<Num&&>(numerator)),
                               util::abs(static_cast<Den&&>(denominator))};
        }
        else {
            return return_type{static_cast<Num&&>(numerator),
                               util::abs(static_cast<Den&&>(denominator))};
        }
    }
}

#endif
