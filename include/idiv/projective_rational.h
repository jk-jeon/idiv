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

#include <compare>
#include <concepts>
#include <type_traits>
#include <utility>

namespace jkj {
    namespace cntfrc {
        template <class Num, class Den>
        struct projective_rational {
            Num numerator;
            Den denominator;
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
            return bool((x.numerator * y.denominator > x.denominator * y.numerator) ^
                        (y.numerator * z.denominator > y.denominator * z.numerator) ^
                        (z.numerator * x.denominator > z.denominator * x.numerator));
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
                return is_zero(x);
            }
            template <class T>
            friend constexpr auto operator<=>(zero, T&& x) noexcept {
                return is_strictly_negative(x) ? std::strong_ordering::greater
                       : is_zero(x)            ? std::strong_ordering::equal
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
            constexpr auto operator()(projective_rational<Num, Den> const& x) {
                return projective_rational{num_to_num * x.numerator + den_to_num * x.denominator,
                                           num_to_den * x.numerator + den_to_den * x.denominator};
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
        // This class is supposed to be used both as a temporary and as a stored lvalue, so we do
        // not strip off the reference.
        template <class NumNum, class DenNum, class NumDen, class DenDen>
        linear_fractional_transform(NumNum&&, DenNum&&, NumDen&&, DenDen&&)
            -> linear_fractional_transform<NumNum, DenNum, NumDen, DenDen>;

        template <class Num, class Den = unity>
        constexpr auto linear_fractional_translation(Num&& numerator, Den&& denominator = {}) {
            return linear_fractional_transform{denominator, static_cast<Num&&>(numerator), zero{},
                                               denominator};
        }
    }
}

#endif
