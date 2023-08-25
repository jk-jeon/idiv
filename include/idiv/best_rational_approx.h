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

#ifndef JKJ_HEADER_BEST_RATIONAL_APPROX
#define JKJ_HEADER_BEST_RATIONAL_APPROX

#include "continued_fractions.h"
#include <cassert>
#include <cstdlib>

namespace jkj {
    template <class Int, class UInt>
    struct best_rational_approx_output {
        frac<Int, UInt> below;
        frac<Int, UInt> above;
    };

    // Find the best rational approximations from below and from above of denominators no more than
    // denominator_upper_bound for the given number x. The number x is given in terms of its
    // continued fractions. The continued fractions calculator cf is assumed to be initialized,
    // i.e., it starts from the first convergent when evaluated.
    template <class ContinuedFractionsImpl, class UInt>
    constexpr best_rational_approx_output<typename ContinuedFractionsImpl::int_type,
                                          typename ContinuedFractionsImpl::uint_type>
    find_best_rational_approx(ContinuedFractionsImpl& cf, UInt const& denominator_upper_bound) {
        util::constexpr_assert<util::error_msgs::no_error_msg>(
            is_strictly_positive(denominator_upper_bound));

        using int_type = typename ContinuedFractionsImpl::int_type;
        using uint_type = typename ContinuedFractionsImpl::uint_type;
        best_rational_approx_output<int_type, uint_type> ret_value;

        // First, find the last convergent whose denominator is bounded above by the given upper
        // bound.
        frac<int_type, uint_type> previous_convergent;
        frac<int_type, uint_type> current_convergent;
        do {
            previous_convergent = cf.previous_convergent();
            current_convergent = cf.current_convergent();

            // Obtain the next convergent.
            if (!cf.update()) {
                // If there is no more convergents, we already obtained the perfect approximation.
                ret_value.below = cf.current_convergent();
                ret_value.above = cf.current_convergent();
                return ret_value;
            }
        } while (cf.current_denominator() <= denominator_upper_bound);

        // If the current convergent is of even index,
        // then the current convergent is the best approximation from below,
        // and the best approximation from above is the last semiconvergent.
        // If the current convergent is of odd index, then the other way around.
        // Note that cf.current_index() is one larger than the index of the current convergent,
        // so we need to reverse the parity.

        auto compute_bounds = [&](auto& major, auto& minor) {
            // The current convergent is the best approximation from below.
            major = current_convergent;

            // The best approximation from above is the last semiconvergent.
            using std::div;
            auto semiconvergent_coeff =
                div(denominator_upper_bound - previous_convergent.denominator,
                    current_convergent.denominator)
                    .quot;

            minor.numerator =
                previous_convergent.numerator + semiconvergent_coeff * current_convergent.numerator;
            minor.denominator = previous_convergent.denominator +
                                semiconvergent_coeff * current_convergent.denominator;
        };

        if (cf.current_index() % 2 == 1) {
            compute_bounds(ret_value.below, ret_value.above);
        }
        else {
            compute_bounds(ret_value.above, ret_value.below);
        }

        return ret_value;
    }
}

#endif
