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

#include <idiv/wuint.h>
#include <boost/ut.hpp>

void wuint_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Wide unsigned integers]"_test = [] {
        should("add_carry64") = [] {
            unsigned int carry = 0;
            expect(wuint::add_carry64(UINT64_C(0xffff'ffff'ffff'fff0), 0x10, carry) == 0);
            expect(wuint::add_carry64(UINT64_C(0xffff'ffff'ffff'fff0), 0x11, carry) == 2);
            expect(wuint::add_carry64(0, 10, carry) == 11);
            expect(wuint::add_carry64(0, 10, carry) == 10);
        };

        should("sub_borrow64") = [] {
            unsigned int borrow = 0;
            expect(wuint::sub_borrow64(0, 1, borrow) == UINT64_C(0xffff'ffff'ffff'ffff));
            expect(wuint::sub_borrow64(0, 1, borrow) == UINT64_C(0xffff'ffff'ffff'fffe));
            expect(wuint::sub_borrow64(10, 1, borrow) == 8);
            expect(wuint::sub_borrow64(10, 1, borrow) == 9);
        };

        should("uint128::operator+=") = [] {
            wuint::uint128 x{0xffff'ffff, UINT64_C(0xffff'ffff'ffff'ffff)};
            x += UINT64_C(0xffff'ffff'ffff'ffff);
            expect(x.high() == UINT64_C(0x1'0000'0000) &&
                   x.low() == UINT64_C(0xffff'ffff'ffff'fffe));
        };

        should("umul128 and umul128_upper64") = [] {
            constexpr auto x =
                wuint::umul128(UINT64_C(0x12345678'12345678), UINT64_C(0x87654321'87654321));
            expect(x.high() == 0x9a0cd05'83fa2782 && x.low() == 0xeb11e7f5'70b88d78);
            expect(wuint::umul128_upper64(UINT64_C(0x12345678'12345678),
                                          UINT64_C(0x87654321'87654321)) == 0x9a0cd05'83fa2782);
        };
    };
}