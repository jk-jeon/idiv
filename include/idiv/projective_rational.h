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

#ifndef JKJ_HEADER_PROJECTIVE_RATIONAL
#define JKJ_HEADER_PROJECTIVE_RATIONAL

#include "frac.h"

namespace jkj {
    namespace cntfrc {
        template <class Num, class Den>
        struct projective_rational {
            Num numerator;
            Den denominator;

            projective_rational() = default;

            template <class OtherNum = Num, class OtherDen = Den>
                requires requires(OtherNum num, OtherDen den) {
                    Num{num};
                    Den{den};
                }
            explicit constexpr projective_rational(OtherNum&& num, OtherDen&& den)
                : numerator{static_cast<OtherNum&&>(num)},
                  denominator{static_cast<OtherDen&&>(den)} {}

            template <class OtherNum, class OtherDen>
                requires(requires(OtherNum num, OtherDen den) {
                            Num{num};
                            Den{den};
                        } && !(std::is_same_v<Num, OtherNum> && std::is_same_v<Den, OtherDen>))
            explicit constexpr projective_rational(
                projective_rational<OtherNum, OtherDen> const& other)
                : numerator{other.numerator}, denominator{other.denominator} {}

            template <class OtherNum, class OtherDen>
                requires(requires(OtherNum num, OtherDen den) {
                            Num{num};
                            Den{den};
                        } && !(std::is_same_v<Num, OtherNum> && std::is_same_v<Den, OtherDen>))
            explicit constexpr projective_rational(projective_rational<OtherNum, OtherDen>&& other)
                : numerator{static_cast<OtherNum&&>(other.numerator)},
                  denominator{static_cast<OtherDen&&>(other.denominator)} {}
        };
        template <class Num, class Den>
        projective_rational(Num&&, Den&&)
            -> projective_rational<std::remove_cvref_t<Num>, std::remove_cvref_t<Den>>;

        template <class Num1, class Den1, class Num2, class Den2>
        constexpr bool operator==(projective_rational<Num1, Den1> const& x,
                                  projective_rational<Num2, Den2> const& y) {
            return x.numerator * y.denominator == y.numerator * x.denominator;
        }

        template <class Num1, class Den1, class Num2, class Den2, class Num3, class Den3>
        constexpr bool cyclic_order(projective_rational<Num1, Den1> const& x,
                                    projective_rational<Num2, Den2> const& y,
                                    projective_rational<Num3, Den3> const& z) {
            int const sign1 = util::strong_order_to_int(x.numerator * y.denominator <=>
                                                        x.denominator * y.numerator);
            int const sign2 = util::strong_order_to_int(y.numerator * z.denominator <=>
                                                        y.denominator * z.numerator);
            int const sign3 = util::strong_order_to_int(z.numerator * x.denominator <=>
                                                        z.denominator * x.numerator);

            return sign1 * sign2 * sign3 > 0;
        }

        template <class Num, class Den>
        constexpr frac<Num, Den> project_to_rational(projective_rational<Num, Den> const& x) {
            return frac{x.numerator, x.denominator};
        }
        template <class Num, class Den>
        constexpr frac<Num, Den> project_to_rational(projective_rational<Num, Den>&& x) {
            return frac{std::move(x).numerator, std::move(x).denominator};
        }

        // May use this type to replace constant 0.
        struct zero {
            zero() = default;

            constexpr zero operator*(zero) const noexcept { return {}; }
            template <class T>
            constexpr zero operator*(T&&) const noexcept {
                return {};
            }
            template <class T>
            constexpr zero& operator*=(T&&) noexcept {
                return *this;
            }
            template <class T>
            friend constexpr zero operator*(T&&, zero) noexcept {
                return {};
            }
            template <class T>
            constexpr zero operator/(T&&) const noexcept {
                return {};
            }
            template <class T>
            constexpr zero& operator/=(T&&) noexcept {
                return *this;
            }
            constexpr zero operator+(zero) const noexcept { return {}; }
            template <class T>
            constexpr auto&& operator+(T&& x) const noexcept {
                return static_cast<T&&>(x);
            }
            template <class T>
            friend constexpr auto&& operator+(T&& x, zero) noexcept {
                return static_cast<T&&>(x);
            }
            template <class T>
            friend constexpr auto& operator+=(T& x, zero) noexcept {
                return x;
            }
            template <class T>
            friend constexpr auto&& operator-(T&& x, zero) noexcept {
                return static_cast<T&&>(x);
            }
            template <class T>
            friend constexpr auto& operator-=(T& x, zero) noexcept {
                return x;
            }
            template <class T>
            friend constexpr bool operator==(zero, T&& x) noexcept {
                return util::is_zero(x);
            }
            template <class T>
            friend constexpr auto operator<=>(zero, T&& x) noexcept {
                return util::is_strictly_negative(x) ? std::strong_ordering::greater
                       : util::is_zero(x)            ? std::strong_ordering::equal
                                                     : std::strong_ordering::less;
            }
        };

        // May use this type to replace constant 1.
        struct unity {
            unity() = default;

            constexpr unity operator*(unity) const noexcept { return {}; }
            template <class T>
            constexpr auto&& operator*(T&& x) const noexcept {
                return static_cast<T&&>(x);
            }
            template <class T>
            friend constexpr auto&& operator*(T&& x, unity) noexcept {
                return static_cast<T&&>(x);
            }
            template <class T>
            friend constexpr T& operator*=(T& x, unity) noexcept {
                return x;
            }
            template <class T>
            friend constexpr auto&& operator/(T&& x, unity) noexcept {
                return static_cast<T&&>(x);
            }
            template <class T>
            friend constexpr T& operator/=(T& x, unity) noexcept {
                return x;
            }
        };

        template <class NumNum, class DenNum = NumNum, class NumDen = NumNum, class DenDen = NumNum>
        struct linear_fractional_transform {
            // (ax+b)/(cx+d)
            NumNum num_to_num; // a
            DenNum den_to_num; // b
            NumDen num_to_den; // c
            DenDen den_to_den; // d

            template <class Num, class Den>
            constexpr auto operator()(projective_rational<Num, Den> const& x) const {
                return projective_rational{num_to_num * x.numerator + den_to_num * x.denominator,
                                           num_to_den * x.numerator + den_to_den * x.denominator};
            }

            constexpr int determinant_sign() const {
                return util::strong_order_to_int(num_to_num * den_to_den <=>
                                                 den_to_num * num_to_den);
            }

            // Multiply the matrix (t s \\ 0 t) from left for a number s/t.
            template <class Num, class Den = unity>
            constexpr void translate(Num&& numerator, Den&& denominator = {}) {
                num_to_num *= denominator;
                num_to_num += numerator * num_to_den;
                den_to_num *= denominator;
                den_to_num += numerator * den_to_den;
                num_to_den *= denominator;
                den_to_den *= denominator;
            }

            // Multiply the matrix (0 1 \\ 1 0) from left.
            constexpr void reflect() {
                using std::swap;
                swap(num_to_num, num_to_den);
                swap(den_to_num, den_to_den);
            }
        };
        // This class is supposed to be used both as a temporary and as a stored lvalue, so we
        // do not strip off the reference.
        template <class NumNum, class DenNum, class NumDen, class DenDen>
        linear_fractional_transform(NumNum&&, DenNum&&, NumDen&&, DenDen&&)
            -> linear_fractional_transform<NumNum, DenNum, NumDen, DenDen>;

        template <class Num, class Den = unity>
        constexpr auto linear_fractional_translation(Num&& numerator, Den&& denominator = {}) {
            return linear_fractional_transform{denominator, static_cast<Num&&>(numerator), zero{},
                                               denominator};
        }

        template <class XNumYNumNum, class XNumYDenNum = XNumYNumNum,
                  class XDenYNumNum = XNumYNumNum, class XDenYDenNum = XNumYNumNum,
                  class XNumYNumDen = XNumYNumNum, class XNumYDenDen = XNumYNumNum,
                  class XDenYNumDen = XNumYNumNum, class XDenYDenDen = XNumYNumNum>
        struct bilinear_fractional_transform {
            XNumYNumNum xnum_ynum_to_num;
            XNumYDenNum xnum_yden_to_num;
            XDenYNumNum xden_ynum_to_num;
            XDenYDenNum xden_yden_to_num;
            XNumYNumDen xnum_ynum_to_den;
            XNumYDenDen xnum_yden_to_den;
            XDenYNumDen xden_ynum_to_den;
            XDenYDenDen xden_yden_to_den;

            template <class XNum, class XDen, class YNum, class YDen>
            constexpr auto operator()(projective_rational<XNum, XDen> const& x,
                                      projective_rational<YNum, YDen> const& y) const {
                return projective_rational{xnum_ynum_to_num * x.numerator * y.numerator +
                                               xnum_yden_to_num * x.numerator * y.denominator +
                                               xden_ynum_to_num * x.denominator * y.numerator +
                                               xden_yden_to_num * x.denominator * y.denominator,
                                           xnum_ynum_to_den * x.numerator * y.numerator +
                                               xnum_yden_to_den * x.numerator * y.denominator +
                                               xden_ynum_to_den * x.denominator * y.numerator +
                                               xden_yden_to_den * x.denominator * y.denominator};
            }

            // Multiply the matrix (t s \\ 0 t) from left for a number s/t.
            template <class Num, class Den = unity>
            constexpr void translate(Num&& numerator, Den&& denominator = {}) {
                xnum_ynum_to_num *= denominator;
                xnum_ynum_to_num += numerator * xnum_ynum_to_den;
                xnum_yden_to_num *= denominator;
                xnum_yden_to_num += numerator * xnum_yden_to_den;
                xden_ynum_to_num *= denominator;
                xden_ynum_to_num += numerator * xden_ynum_to_den;
                xden_yden_to_num *= denominator;
                xden_yden_to_num += numerator * xden_yden_to_den;
                xnum_ynum_to_den *= denominator;
                xnum_yden_to_den *= denominator;
                xden_ynum_to_den *= denominator;
                xden_yden_to_den *= denominator;
            }

            // Multiply the matrix (0 1 \\ 1 0) from left.
            constexpr void reflect() {
                using std::swap;
                swap(xnum_ynum_to_num, xnum_ynum_to_den);
                swap(xnum_yden_to_num, xnum_yden_to_den);
                swap(xden_ynum_to_num, xden_ynum_to_den);
                swap(xden_yden_to_num, xden_yden_to_den);
            }
        };
        // This class is supposed to be used both as a temporary and as a stored lvalue, so we
        // do not strip off the reference.
        template <class XNumYNumNum, class XNumYDenNum, class XDenYNumNum, class XDenYDenNum,
                  class XNumYNumDen, class XNumYDenDen, class XDenYNumDen, class XDenYDenDen>
        bilinear_fractional_transform(XNumYNumNum&&, XNumYDenNum&&, XDenYNumNum&&, XDenYDenNum&&,
                                      XNumYNumDen&&, XNumYDenDen&&, XDenYNumDen&&, XDenYDenDen&&)
            -> bilinear_fractional_transform<XNumYNumNum, XNumYDenNum, XDenYNumNum, XDenYDenNum,
                                             XNumYNumDen, XNumYDenDen, XDenYNumDen, XDenYDenDen>;
    }
}

#endif
