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
#include <idiv/rational_continued_fraction.h>
#include <boost/ut.hpp>

void rational_continued_fraction_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Continued fractions for rational numbers]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        {
            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                cntfrc::impl::rational{convergent_t{156, 179u}});

            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{0, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{1, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{6, 7u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{7, 8u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{27, 31u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{34, 39u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{61, 70u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{156, 179u});
            expect(cf.update() == false);
        }
        {
            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                cntfrc::impl::rational{convergent_t{-2767, 1982u}});

            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-2, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-1, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-3, 2u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-4, 3u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-7, 5u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-67, 48u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-74, 53u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-141, 101u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-497, 356u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-1135, 813u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-2767, 1982u});
            expect(cf.update() == false);
        }
    };
}

