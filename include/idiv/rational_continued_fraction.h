// Copyright 2022 Junekey Jeon
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

#ifndef JKJ_HEADER_RATIONAL_CONTINUED_FRACTION
#define JKJ_HEADER_RATIONAL_CONTINUED_FRACTION

#include "continued_fraction.h"
#include "frac.h"
#include <cassert>
#include <cstdlib>

namespace jkj {
    template <class Int, class UInt, class Unity = unity>
    class rational_continued_fraction {
        Int prev_error_;
        Int curr_error_;

    public:
        using partial_fraction_type = frac<Unity, Int>;
        using convergent_type = frac<Int, UInt>;

        constexpr rational_continued_fraction(frac<Int, UInt> r)
            : prev_error_{static_cast<Int&&>(r.numerator)},
              curr_error_{static_cast<UInt&&>(r.denominator)} {
            util::constexpr_assert<util::error_msgs::divide_by_zero>(!is_zero(curr_error_));
        }

        constexpr next_partial_fraction_return<partial_fraction_type> next_partial_fraction() {
            util::constexpr_assert<util::error_msgs::divide_by_zero>(!is_zero(curr_error_));

            using std::div;
            auto div_result = div(prev_error_, abs(curr_error_));
            prev_error_ = static_cast<Int&&>(curr_error_);
            curr_error_ = Int(static_cast<UInt&&>(div_result.rem));
            return {{Unity{1}, static_cast<Int&&>(div_result.quot)}, is_zero(curr_error_)};
        }
    };
}

#endif
