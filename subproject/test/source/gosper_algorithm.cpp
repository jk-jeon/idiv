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
#include <idiv/continued_fraction/engine/gosper.h>
#include <idiv/continued_fraction/engine/rational.h>
#include <boost/ut.hpp>

#include <vector>
#include <variant>

void gosper_algorithm_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Gosper's algorithms]"_test = [] {
        using rational_continued_fraction_t =
            cntfrc::engine::rational<bigint::int_var, bigint::uint_var>;

        should("unary_gosper") = [] {
            // 481/2245 = (-18*156 + 13*179)/(12*156 - 23*179)
            auto cf1 = cntfrc::make_unary_gosper_generator<cntfrc::convergent_tracker>(
                rational_continued_fraction_t{156, 179u}, {-18, 13, 12, -23});
            auto cf2 = cntfrc::make_generator<cntfrc::convergent_tracker>(
                rational_continued_fraction_t{481, 2245u});

            auto itr1 = cf1.begin();
            auto itr2 = cf2.begin();
            for (; itr2 != cf2.end(); ++itr1, ++itr2) {
                expect(itr1 != cf1.end());
                expect(itr1->current_convergent() == itr2->current_convergent());
            }
            expect(itr1 == cf1.end());
        };

        should("binary_gosper") = [] {
            // Take x = 17/89, y = 31/125, and
            // z = (xy + 4x + 2y + 8)/(2x - 3y + 1) = 2655/182.
            auto cf1 = cntfrc::make_binary_gosper_generator<cntfrc::convergent_tracker>(
                rational_continued_fraction_t{17, 89u}, rational_continued_fraction_t{31, 125u},
                {// numerator
                 1, 4, 2, 8,
                 // denominator
                 0, 2, -3, 1});
            auto cf2 = cntfrc::make_generator<cntfrc::convergent_tracker>(
                rational_continued_fraction_t{2655, 182u});

            auto itr1 = cf1.begin();
            auto itr2 = cf2.begin();
            for (; itr2 != cf2.end(); ++itr1, ++itr2) {
                expect(itr1 != cf1.end());
                expect(itr1->current_convergent() == itr2->current_convergent());
            }
            expect(itr1 == cf1.end());
        };
    };
}
