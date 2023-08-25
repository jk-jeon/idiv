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

#ifndef JKJ_HEADER_RATIONAL_CONTINUED_FRACTIONS
#define JKJ_HEADER_RATIONAL_CONTINUED_FRACTIONS

#include "continued_fractions.h"
#include <cassert>
#include <cstdlib>

namespace jkj {
    template <class Int, class UInt>
    class rational_continued_fractions
        : public continued_fractions<rational_continued_fractions<Int, UInt>, Int, UInt> {
        using crtp_base = continued_fractions<rational_continued_fractions<Int, UInt>, Int, UInt>;
        friend crtp_base;

        Int prev_error_;
        Int curr_error_;

        constexpr Int compute_next_coefficient() {
            util::constexpr_assert<util::error_msgs::divide_by_zero>(!is_zero(curr_error_));

            using std::div;
            auto div_result = div(prev_error_, abs(curr_error_));
            prev_error_ = static_cast<Int&&>(curr_error_);
            curr_error_ = Int(static_cast<UInt&&>(div_result.rem));

            if (is_zero(curr_error_)) {
                crtp_base::set_terminate_flag();
            }

            return static_cast<Int&&>(div_result.quot);
        }

    public:
        constexpr rational_continued_fractions(frac<Int, UInt> r)
            : prev_error_{static_cast<Int&&>(r.numerator)},
              curr_error_{static_cast<UInt&&>(r.denominator)} {
            util::constexpr_assert<util::error_msgs::divide_by_zero>(!is_zero(curr_error_));
        }
    };
}

#endif
