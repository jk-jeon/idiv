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

#include <idiv/util.h>
#include <boost/ut.hpp>

void integer_util_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Integer utilities]"_test = [] {
        auto to_signed_from_unsigned = util::to_signed(3u);
        auto to_signed_from_signed = util::to_signed(3);
        auto to_negative_from_unsigned = util::to_negative(2147483647u);
        auto abs_from_unsigned = util::abs(4294967295u);
        auto abs_from_signed = util::abs(static_cast<int>(-2147483648ll));
        auto invert_sign = util::invert_sign(2147483647);
        expect(std::is_same_v<decltype(to_signed_from_unsigned), int> &&
               to_signed_from_unsigned == 3);
        expect(std::is_same_v<decltype(to_signed_from_signed), int> && to_signed_from_signed == 3);
        expect(std::is_same_v<decltype(to_negative_from_unsigned), int> &&
               to_negative_from_unsigned == -2147483647);
        expect(std::is_same_v<decltype(abs_from_unsigned), unsigned int> &&
               abs_from_unsigned == 4294967295u);
        expect(std::is_same_v<decltype(abs_from_signed), unsigned int> &&
               abs_from_signed == 2147483648u);
        expect(std::is_same_v<decltype(invert_sign), int> && invert_sign == -2147483647);

        auto div_uint_uint = util::div(18u, 7u);
        auto div_int_uint_positive = util::div(37, 9u);
        auto div_int_uint_positive_divisible = util::div(36, 9u);
        auto div_int_uint_negative = util::div(-37, 9u);
        auto div_int_uint_negative_divisible = util::div(-36, 9u);
        expect(std::is_same_v<decltype(div_uint_uint.quot), unsigned int> &&
               div_uint_uint.quot == 2u &&
               std::is_same_v<decltype(div_uint_uint.rem), unsigned int> &&
               div_uint_uint.rem == 4u);
        expect(std::is_same_v<decltype(div_int_uint_positive.quot), int> &&
               div_int_uint_positive.quot == 4 &&
               std::is_same_v<decltype(div_int_uint_positive.rem), unsigned int> &&
               div_int_uint_positive.rem == 1u);
        expect(std::is_same_v<decltype(div_int_uint_positive_divisible.quot), int> &&
               div_int_uint_positive_divisible.quot == 4 &&
               std::is_same_v<decltype(div_int_uint_positive_divisible.rem), unsigned int> &&
               div_int_uint_positive_divisible.rem == 0u);
        expect(std::is_same_v<decltype(div_int_uint_negative.quot), int> &&
               div_int_uint_negative.quot == -5 &&
               std::is_same_v<decltype(div_int_uint_negative.rem), unsigned int> &&
               div_int_uint_negative.rem == 8u);
        expect(std::is_same_v<decltype(div_int_uint_negative_divisible.quot), int> &&
               div_int_uint_negative_divisible.quot == -4 &&
               std::is_same_v<decltype(div_int_uint_negative_divisible.rem), unsigned int> &&
               div_int_uint_negative_divisible.rem == 0u);

        auto div_floor_uint_uint = util::div_floor(71u, 9u);
        auto div_floor_uint_uint_divisible = util::div_floor(72u, 9u);
        auto div_floor_int_uint_positive = util::div_floor(71, 9u);
        auto div_floor_int_uint_positive_divisible = util::div_floor(72, 9u);
        auto div_floor_int_uint_negative = util::div_floor(-71, 9u);
        auto div_floor_int_uint_negative_divisible = util::div_floor(-72, 9u);
        auto div_floor_uint_int_positive = util::div_floor(71u, 9);
        auto div_floor_uint_int_positive_divisible = util::div_floor(72u, 9);
        auto div_floor_uint_int_negative = util::div_floor(71u, -9);
        auto div_floor_uint_int_negative_divisible = util::div_floor(72u, -9);
        auto div_floor_int_int_positive = util::div_floor(71, 9);
        auto div_floor_int_int_positive_divisible = util::div_floor(72, 9);
        auto div_floor_int_int_negative = util::div_floor(71, -9);
        auto div_floor_int_int_negative_divisible = util::div_floor(72, -9);
        expect(std::is_same_v<decltype(div_floor_uint_uint), unsigned int> &&
               div_floor_uint_uint == 7u);
        expect(std::is_same_v<decltype(div_floor_uint_uint_divisible), unsigned int> &&
               div_floor_uint_uint_divisible == 8u);
        expect(std::is_same_v<decltype(div_floor_int_uint_positive), int> &&
               div_floor_int_uint_positive == 7);
        expect(std::is_same_v<decltype(div_floor_int_uint_positive_divisible), int> &&
               div_floor_int_uint_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_floor_int_uint_negative), int> &&
               div_floor_int_uint_negative == -8);
        expect(std::is_same_v<decltype(div_floor_int_uint_negative_divisible), int> &&
               div_floor_int_uint_negative_divisible == -8);
        expect(std::is_same_v<decltype(div_floor_uint_int_positive), int> &&
               div_floor_uint_int_positive == 7);
        expect(std::is_same_v<decltype(div_floor_uint_int_positive_divisible), int> &&
               div_floor_uint_int_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_floor_uint_int_negative), int> &&
               div_floor_uint_int_negative == -8);
        expect(std::is_same_v<decltype(div_floor_uint_int_negative_divisible), int> &&
               div_floor_uint_int_negative_divisible == -8);
        expect(std::is_same_v<decltype(div_floor_int_int_positive), int> &&
               div_floor_int_int_positive == 7);
        expect(std::is_same_v<decltype(div_floor_int_int_positive_divisible), int> &&
               div_floor_int_int_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_floor_int_int_negative), int> &&
               div_floor_int_int_negative == -8);
        expect(std::is_same_v<decltype(div_floor_int_int_negative_divisible), int> &&
               div_floor_int_int_negative_divisible == -8);

        auto div_ceil_uint_uint = util::div_ceil(71u, 9u);
        auto div_ceil_uint_uint_divisible = util::div_ceil(72u, 9u);
        auto div_ceil_int_uint_positive = util::div_ceil(71, 9u);
        auto div_ceil_int_uint_positive_divisible = util::div_ceil(72, 9u);
        auto div_ceil_int_uint_negative = util::div_ceil(-71, 9u);
        auto div_ceil_int_uint_negative_divisible = util::div_ceil(-72, 9u);
        auto div_ceil_uint_int_positive = util::div_ceil(71u, 9);
        auto div_ceil_uint_int_positive_divisible = util::div_ceil(72u, 9);
        auto div_ceil_uint_int_negative = util::div_ceil(71u, -9);
        auto div_ceil_uint_int_negative_divisible = util::div_ceil(72u, -9);
        auto div_ceil_int_int_positive = util::div_ceil(71, 9);
        auto div_ceil_int_int_positive_divisible = util::div_ceil(72, 9);
        auto div_ceil_int_int_negative = util::div_ceil(71, -9);
        auto div_ceil_int_int_negative_divisible = util::div_ceil(72, -9);
        expect(std::is_same_v<decltype(div_ceil_uint_uint), unsigned int> &&
               div_ceil_uint_uint == 8u);
        expect(std::is_same_v<decltype(div_ceil_uint_uint_divisible), unsigned int> &&
               div_ceil_uint_uint_divisible == 8u);
        expect(std::is_same_v<decltype(div_ceil_int_uint_positive), int> &&
               div_ceil_int_uint_positive == 8);
        expect(std::is_same_v<decltype(div_ceil_int_uint_positive_divisible), int> &&
               div_ceil_int_uint_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_ceil_int_uint_negative), int> &&
               div_ceil_int_uint_negative == -7);
        expect(std::is_same_v<decltype(div_ceil_int_uint_negative_divisible), int> &&
               div_ceil_int_uint_negative_divisible == -8);
        expect(std::is_same_v<decltype(div_ceil_uint_int_positive), int> &&
               div_ceil_uint_int_positive == 8);
        expect(std::is_same_v<decltype(div_ceil_uint_int_positive_divisible), int> &&
               div_ceil_uint_int_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_ceil_uint_int_negative), int> &&
               div_ceil_uint_int_negative == -7);
        expect(std::is_same_v<decltype(div_ceil_uint_int_negative_divisible), int> &&
               div_ceil_uint_int_negative_divisible == -8);
        expect(std::is_same_v<decltype(div_ceil_int_int_positive), int> &&
               div_ceil_int_int_positive == 8);
        expect(std::is_same_v<decltype(div_ceil_int_int_positive_divisible), int> &&
               div_ceil_int_int_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_ceil_int_int_negative), int> &&
               div_ceil_int_int_negative == -7);
        expect(std::is_same_v<decltype(div_ceil_int_int_negative_divisible), int> &&
               div_ceil_int_int_negative_divisible == -8);

        // max(floor(log2(x / y)), 0).
        expect(util::trunc_floor_log2_div(7u, 1u) == 2);
        expect(util::trunc_floor_log2_div(7u, 2u) == 1);
        expect(util::trunc_floor_log2_div(7u, 3u) == 1);
        expect(util::trunc_floor_log2_div(7u, 4u) == 0);
        expect(util::trunc_floor_log2_div(7u, 5u) == 0);
        expect(util::trunc_floor_log2_div(7u, 6u) == 0);
        expect(util::trunc_floor_log2_div(7u, 7u) == 0);
        expect(util::trunc_floor_log2_div(7u, 8u) == 0);
        expect(util::trunc_floor_log2_div(7u, 9u) == 0);
        expect(util::trunc_floor_log2_div(7u, 10u) == 0);
        expect(util::trunc_floor_log2_div(16u, 1u) == 4);
        expect(util::trunc_floor_log2_div(16u, 2u) == 3);
        expect(util::trunc_floor_log2_div(16u, 3u) == 2);
        expect(util::trunc_floor_log2_div(16u, 4u) == 2);
        expect(util::trunc_floor_log2_div(16u, 5u) == 1);
        expect(util::trunc_floor_log2_div(16u, 6u) == 1);
        expect(util::trunc_floor_log2_div(16u, 7u) == 1);
        expect(util::trunc_floor_log2_div(16u, 8u) == 1);
        expect(util::trunc_floor_log2_div(16u, 9u) == 0);
        expect(util::trunc_floor_log2_div(16u, 10u) == 0);

        // max(ceil(log2(x / y)), 0)
        expect(util::trunc_ceil_log2_div(7u, 1u) == 3);
        expect(util::trunc_ceil_log2_div(7u, 2u) == 2);
        expect(util::trunc_ceil_log2_div(7u, 3u) == 2);
        expect(util::trunc_ceil_log2_div(7u, 4u) == 1);
        expect(util::trunc_ceil_log2_div(7u, 5u) == 1);
        expect(util::trunc_ceil_log2_div(7u, 6u) == 1);
        expect(util::trunc_ceil_log2_div(7u, 7u) == 0);
        expect(util::trunc_ceil_log2_div(7u, 8u) == 0);
        expect(util::trunc_ceil_log2_div(7u, 9u) == 0);
        expect(util::trunc_ceil_log2_div(7u, 10u) == 0);
        expect(util::trunc_ceil_log2_div(16u, 1u) == 4);
        expect(util::trunc_ceil_log2_div(16u, 2u) == 3);
        expect(util::trunc_ceil_log2_div(16u, 3u) == 3);
        expect(util::trunc_ceil_log2_div(16u, 4u) == 2);
        expect(util::trunc_ceil_log2_div(16u, 5u) == 2);
        expect(util::trunc_ceil_log2_div(16u, 6u) == 2);
        expect(util::trunc_ceil_log2_div(16u, 7u) == 2);
        expect(util::trunc_ceil_log2_div(16u, 8u) == 1);
        expect(util::trunc_ceil_log2_div(16u, 9u) == 1);
        expect(util::trunc_ceil_log2_div(16u, 10u) == 1);

        expect(util::pow_uint(10, 0u) == 1);
        expect(util::pow_uint(10, 1u) == 10);
        expect(util::pow_uint(10, 2u) == 100);
        expect(util::pow_uint(10, 3u) == 1'000);
        expect(util::pow_uint(10, 4u) == 10'000);
        expect(util::pow_uint(10, 5u) == 100'000);
        expect(util::pow_uint(10, 6u) == 1'000'000);
        expect(util::pow_uint(10, 7u) == 10'000'000);
        expect(util::pow_uint(10, 8u) == 100'000'000);
        expect(util::pow_uint(10, 9u) == 1'000'000'000);
    };
}

