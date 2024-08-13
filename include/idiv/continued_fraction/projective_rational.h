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
        template <class Num, class Den>
        constexpr projective_rational<Num, Den> projectify(frac<Num, Den> const& x) {
            return frac{x.numerator, x.denominator};
        }
        template <class Num, class Den>
        constexpr projective_rational<Num, Den> projectify(frac<Num, Den>&& x) {
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
            constexpr zero& operator*=(zero) noexcept { return *this; }
            template <class T>
            constexpr zero& operator*=(T&&) noexcept {
                return *this;
            }
            template <class T>
            friend constexpr zero operator*(T&&, zero) noexcept {
                return {};
            }
            template <class T>
                requires(
                    requires { T(0); } || requires { T(0u); })
            friend constexpr T& operator*=(T& x, zero) noexcept {
                return x = static_cast<T>(zero{});
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
            friend constexpr bool is_zero(zero) noexcept { return true; }

            template <class T>
                requires requires { T{0}; }
            explicit operator T() const noexcept {
                return T{0};
            }
            template <class T>
                requires(
                    requires { T{0u}; } && !requires { T{0}; })
            explicit operator T() const noexcept {
                return T{0u};
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
            friend constexpr bool is_zero(unity) noexcept { return false; }

            template <class T>
                requires requires { T{1}; }
            explicit operator T() const noexcept {
                return T{1};
            }
            template <class T>
                requires(
                    requires { T{1u}; } && !requires { T{1}; })
            explicit operator T() const noexcept {
                return T{1u};
            }
        };

        template <class NumNum, class DenNum = NumNum, class NumDen = NumNum, class DenDen = NumNum>
        class linear_fractional_mapping {
        protected:
            // (ax+b)/(cx+d)
            NumNum num_to_num_; // a
            DenNum den_to_num_; // b
            NumDen num_to_den_; // c
            DenDen den_to_den_; // d

        public:
            template <class NumNum_, class DenNum_, class NumDen_, class DenDen_>
            constexpr linear_fractional_mapping(NumNum_&& a, DenNum_&& b, NumDen_&& c,
                                                  DenDen_&& d)
                : num_to_num_{static_cast<NumNum_&&>(a)}, den_to_num_{static_cast<DenNum_&&>(b)},
                  num_to_den_{static_cast<NumDen_&&>(c)}, den_to_den_{static_cast<DenDen_&&>(d)} {
                // The degenerate case a = b = c = d = 0 is disallowed.
                util::constexpr_assert(
                    !util::is_zero(num_to_num()) || !util::is_zero(den_to_num()) ||
                    !util::is_zero(num_to_den()) || !util::is_zero(den_to_den()));
            }

            constexpr NumNum const& num_to_num() const& noexcept { return num_to_num_; }
            constexpr DenNum const& den_to_num() const& noexcept { return den_to_num_; }
            constexpr NumDen const& num_to_den() const& noexcept { return num_to_den_; }
            constexpr DenDen const& den_to_den() const& noexcept { return den_to_den_; }

            constexpr NumNum&& num_to_num() && noexcept {
                return static_cast<NumNum&&>(num_to_num_);
            }
            constexpr DenNum&& den_to_num() && noexcept {
                return static_cast<DenNum&&>(den_to_num_);
            }
            constexpr NumDen&& num_to_den() && noexcept {
                return static_cast<NumDen&&>(num_to_den_);
            }
            constexpr DenDen&& den_to_den() && noexcept {
                return static_cast<DenDen&&>(den_to_den_);
            }

            template <class Num, class Den>
            constexpr auto operator()(projective_rational<Num, Den> const& x) const {
                return projective_rational{
                    num_to_num() * x.numerator + den_to_num() * x.denominator,
                    num_to_den() * x.numerator + den_to_den() * x.denominator};
            }

            constexpr int determinant_sign() const {
                return util::strong_order_to_int(num_to_num_ * den_to_den_ <=>
                                                 den_to_num_ * num_to_den_);
            }
        };
        // This class is supposed to be used both as a temporary and as a stored lvalue, so we
        // do not strip off the lvalue reference.
        template <class NumNum, class DenNum, class NumDen, class DenDen>
        linear_fractional_mapping(NumNum&&, DenNum&&, NumDen&&, DenDen&&)
            -> linear_fractional_mapping<NumNum, DenNum, NumDen, DenDen>;

        template <class Num, class Den = unity>
        constexpr auto linear_fractional_translation(Num&& numerator,
                                                     Den&& denominator = Den{unity{}}) {
            return linear_fractional_mapping{denominator, static_cast<Num&&>(numerator), zero{},
                                               denominator};
        }

        template <class XNumYNumNum, class XNumYDenNum = XNumYNumNum,
                  class XDenYNumNum = XNumYNumNum, class XDenYDenNum = XNumYNumNum,
                  class XNumYNumDen = XNumYNumNum, class XNumYDenDen = XNumYNumNum,
                  class XDenYNumDen = XNumYNumNum, class XDenYDenDen = XNumYNumNum>
        class bilinear_fractional_mapping {
        protected:
            // (axy+bx+cy+d)/(exy+fx+gy+h)
            XNumYNumNum xnum_ynum_to_num_; // a
            XNumYDenNum xnum_yden_to_num_; // b
            XDenYNumNum xden_ynum_to_num_; // c
            XDenYDenNum xden_yden_to_num_; // d
            XNumYNumDen xnum_ynum_to_den_; // e
            XNumYDenDen xnum_yden_to_den_; // f
            XDenYNumDen xden_ynum_to_den_; // g
            XDenYDenDen xden_yden_to_den_; // h

        public:
            template <class XNumYNumNum_, class XNumYDenNum_, class XDenYNumNum_,
                      class XDenYDenNum_, class XNumYNumDen_, class XNumYDenDen_,
                      class XDenYNumDen_, class XDenYDenDen_>
            constexpr bilinear_fractional_mapping(XNumYNumNum_&& a, XNumYDenNum_&& b,
                                                    XDenYNumNum_&& c, XDenYDenNum_&& d,
                                                    XNumYNumDen_&& e, XNumYDenDen_&& f,
                                                    XDenYNumDen_&& g, XDenYDenDen_&& h)
                : xnum_ynum_to_num_{static_cast<XNumYNumNum_&&>(a)},
                  xnum_yden_to_num_{static_cast<XNumYDenNum_&&>(b)},
                  xden_ynum_to_num_{static_cast<XDenYNumNum_&&>(c)},
                  xden_yden_to_num_{static_cast<XDenYDenNum_&&>(d)},
                  xnum_ynum_to_den_{static_cast<XNumYNumDen_&&>(e)},
                  xnum_yden_to_den_{static_cast<XNumYDenDen_&&>(f)},
                  xden_ynum_to_den_{static_cast<XDenYNumDen_&&>(g)},
                  xden_yden_to_den_{static_cast<XDenYDenDen_&&>(h)} {
                // The degenerate case a = b = c = d = e = f = g = h = 0 is disallowed.
                util::constexpr_assert(
                    !util::is_zero(num_to_num()) || !util::is_zero(den_to_num()) ||
                    !util::is_zero(num_to_den()) || !util::is_zero(den_to_den()));
            }

            constexpr XNumYNumNum const& xnum_ynum_to_num() const& noexcept {
                return xnum_ynum_to_num_;
            }
            constexpr XNumYDenNum const& xnum_yden_to_num() const& noexcept {
                return xnum_yden_to_num_;
            }
            constexpr XDenYNumNum const& xden_ynum_to_num() const& noexcept {
                return xden_ynum_to_num_;
            }
            constexpr XDenYDenNum const& xden_yden_to_num() const& noexcept {
                return xden_yden_to_num_;
            }
            constexpr XNumYNumDen const& xnum_ynum_to_den() const& noexcept {
                return xnum_ynum_to_den_;
            }
            constexpr XNumYDenDen const& xnum_yden_to_den() const& noexcept {
                return xnum_yden_to_den_;
            }
            constexpr XDenYNumDen const& xden_ynum_to_den() const& noexcept {
                return xden_ynum_to_den_;
            }
            constexpr XDenYDenDen const& xden_yden_to_den() const& noexcept {
                return xden_yden_to_den_;
            }

            constexpr XNumYNumNum&& xnum_ynum_to_num() && noexcept {
                return static_cast<XNumYNumNum&&>(xnum_ynum_to_num_);
            }
            constexpr XNumYDenNum&& xnum_yden_to_num() && noexcept {
                return static_cast<XNumYDenNum&&>(xnum_yden_to_num_);
            }
            constexpr XDenYNumNum&& xden_ynum_to_num() && noexcept {
                return static_cast<XDenYNumNum&&>(xden_ynum_to_num_);
            }
            constexpr XDenYDenNum&& xden_yden_to_num() && noexcept {
                return static_cast<XDenYDenNum&&>(xden_yden_to_num_);
            }
            constexpr XNumYNumDen&& xnum_ynum_to_den() && noexcept {
                return static_cast<XNumYNumDen&&>(xnum_ynum_to_den_);
            }
            constexpr XNumYDenDen&& xnum_yden_to_den() && noexcept {
                return static_cast<XNumYDenDen&&>(xnum_yden_to_den_);
            }
            constexpr XDenYNumDen&& xden_ynum_to_den() && noexcept {
                return static_cast<XDenYNumDen&&>(xden_ynum_to_den_);
            }
            constexpr XDenYDenDen&& xden_yden_to_den() && noexcept {
                return static_cast<XDenYDenDen&&>(xden_yden_to_den_);
            }

            template <class XNum, class XDen, class YNum, class YDen>
            constexpr auto operator()(projective_rational<XNum, XDen> const& x,
                                      projective_rational<YNum, YDen> const& y) const {
                return projective_rational{xnum_ynum_to_num() * x.numerator * y.numerator +
                                               xnum_yden_to_num() * x.numerator * y.denominator +
                                               xden_ynum_to_num() * x.denominator * y.numerator +
                                               xden_yden_to_num() * x.denominator * y.denominator,
                                           xnum_ynum_to_den() * x.numerator * y.numerator +
                                               xnum_yden_to_den() * x.numerator * y.denominator +
                                               xden_ynum_to_den() * x.denominator * y.numerator +
                                               xden_yden_to_den() * x.denominator * y.denominator};
            }

            template <class XNum, class XDen>
            constexpr auto specialize_x(projective_rational<XNum, XDen> const& x) const {
                return linear_fractional_mapping{
                    xnum_ynum_to_num() * x.numerator + xden_ynum_to_num() * x.denominator,
                    xnum_yden_to_num() * x.numerator + xden_yden_to_num() * x.denominator,
                    xnum_ynum_to_den() * x.numerator + xden_ynum_to_den() * x.denominator,
                    xnum_yden_to_den() * x.numerator + xden_yden_to_den() * x.denominator};
            }

            template <class YNum, class YDen>
            constexpr auto specialize_y(projective_rational<YNum, YDen> const& y) const {
                return linear_fractional_mapping{
                    xnum_ynum_to_num() * y.numerator + xnum_yden_to_num() * y.denominator,
                    xden_ynum_to_num() * y.numerator + xden_yden_to_num() * y.denominator,
                    xnum_ynum_to_den() * y.numerator + xnum_yden_to_den() * y.denominator,
                    xden_ynum_to_den() * y.numerator + xden_yden_to_den() * y.denominator};
            }
        };
        // This class is supposed to be used both as a temporary and as a stored lvalue, so we
        // do not strip off the lvalue reference.
        template <class XNumYNumNum, class XNumYDenNum, class XDenYNumNum, class XDenYDenNum,
                  class XNumYNumDen, class XNumYDenDen, class XDenYNumDen, class XDenYDenDen>
        bilinear_fractional_mapping(XNumYNumNum&&, XNumYDenNum&&, XDenYNumNum&&, XDenYDenNum&&,
                                      XNumYNumDen&&, XNumYDenDen&&, XDenYNumDen&&, XDenYDenDen&&)
            -> bilinear_fractional_mapping<XNumYNumNum, XNumYDenNum, XDenYNumNum, XDenYDenNum,
                                             XNumYNumDen, XNumYDenDen, XDenYNumDen, XDenYDenDen>;
    }
}

#endif
