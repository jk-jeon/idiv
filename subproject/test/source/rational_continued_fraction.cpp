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

#include <idiv/bigint.h>
#include <idiv/continued_fraction/generator.h>
#include <idiv/continued_fraction/engine/rational.h>
#include <boost/ut.hpp>

void rational_continued_fraction_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Continued fractions for rational numbers]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        {
            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                cntfrc::engine::rational{convergent_t{156, 179u}});

            auto itr = cf.begin();
            expect(itr->current_convergent() == convergent_t{0, 1u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{1, 1u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{6, 7u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{7, 8u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{27, 31u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{34, 39u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{61, 70u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{156, 179u});
            expect(++itr == cf.end());
        }
        {
            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                cntfrc::engine::rational{convergent_t{-2767, 1982u}});

            auto itr = cf.begin();
            expect(itr->current_convergent() == convergent_t{-2, 1u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-1, 1u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-3, 2u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-4, 3u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-7, 5u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-67, 48u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-74, 53u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-141, 101u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-497, 356u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-1135, 813u});
            expect(++itr != cf.end());
            expect(itr->current_convergent() == convergent_t{-2767, 1982u});
            expect(++itr == cf.end());
        }
    };

    "[reduce_fraction]"_test = [] {
        using jkj::cntfrc::projective_rational;

        // Primitive unsigned/unsigned.
        {
            auto result = cntfrc::reduce_fraction(UINT32_C(22041843), UINT32_C(225055437));
            expect(std::is_same_v<decltype(result),
                                  projective_rational<std::uint_least32_t, std::uint_least32_t>>);
            expect(result.numerator == UINT32_C(16813));
            expect(result.denominator == UINT32_C(171667));
        }

        // Primitive signed/unsigned.
        {
            auto result = cntfrc::reduce_fraction(-INT32_C(22041843), UINT32_C(225055437));
            expect(std::is_same_v<decltype(result),
                                  projective_rational<std::int_least32_t, std::uint_least32_t>>);
            expect(result.numerator == -INT32_C(16813));
            expect(result.denominator == UINT32_C(171667));
        }

        // Primitive unsigned/signed.
        {
            auto result = cntfrc::reduce_fraction(UINT32_C(22041843), -INT32_C(225055437));
            expect(std::is_same_v<decltype(result),
                                  projective_rational<std::int_least32_t, std::uint_least32_t>>);
            expect(result.numerator == -INT32_C(16813));
            expect(result.denominator == UINT32_C(171667));
        }

        // Primitive signed/signed.
        {
            auto result = cntfrc::reduce_fraction(-INT32_C(22041843), -INT32_C(225055437));
            expect(std::is_same_v<decltype(result),
                                  projective_rational<std::int_least32_t, std::uint_least32_t>>);
            expect(result.numerator == INT32_C(16813));
            expect(result.denominator == UINT32_C(171667));
        }

        // jkj::bigint unsigned/unsigned.
        {
            auto num = jkj::bigint::uint_var{UINT32_C(22041843)};
            auto den = jkj::bigint::uint_var{UINT32_C(225055437)};
            auto result = cntfrc::reduce_fraction(num, den);
            expect(
                std::is_same_v<decltype(result),
                               projective_rational<jkj::bigint::uint_var, jkj::bigint::uint_var>>);
            expect(result.numerator == UINT32_C(16813));
            expect(result.denominator == UINT32_C(171667));
        }

        // jkj::bigint signed/unsigned.
        {
            auto num = jkj::bigint::int_var{jkj::bigint::sign_t::negative, UINT32_C(22041843)};
            auto den = jkj::bigint::uint_var{UINT32_C(225055437)};
            auto result = cntfrc::reduce_fraction(num, den);
            expect(
                std::is_same_v<decltype(result),
                               projective_rational<jkj::bigint::int_var, jkj::bigint::uint_var>>);
            expect(result.numerator == -INT32_C(16813));
            expect(result.denominator == UINT32_C(171667));
        }

        // jkj::bigint unsigned/signed.
        {
            auto num = jkj::bigint::uint_var{UINT32_C(22041843)};
            auto den = jkj::bigint::int_var{jkj::bigint::sign_t::negative, UINT32_C(225055437)};
            auto result = cntfrc::reduce_fraction(num, den);
            expect(
                std::is_same_v<decltype(result),
                               projective_rational<jkj::bigint::int_var, jkj::bigint::uint_var>>);
            expect(result.numerator == -INT32_C(16813));
            expect(result.denominator == UINT32_C(171667));
        }

        // jkj::bigint signed/signed.
        {
            auto num = jkj::bigint::int_var{jkj::bigint::sign_t::negative, UINT32_C(22041843)};
            auto den = jkj::bigint::int_var{jkj::bigint::sign_t::negative, UINT32_C(225055437)};
            auto result = cntfrc::reduce_fraction(num, den);
            expect(
                std::is_same_v<decltype(result),
                               projective_rational<jkj::bigint::int_var, jkj::bigint::uint_var>>);
            expect(result.numerator == INT32_C(16813));
            expect(result.denominator == UINT32_C(171667));
        }
    };
}
