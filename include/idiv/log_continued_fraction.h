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

#ifndef JKJ_HEADER_LOG_CONTINUED_FRACTION
#define JKJ_HEADER_LOG_CONTINUED_FRACTION

#include "rational_continued_fraction.h"
#include "util.h"

namespace jkj {
    // Computes an arbitrary-precision rational approximation of ln(a/b) for positive integers a, b
    // using the continued fraction expansion of 2atanh(z) = ln((1+z)/(1-z)):
    //
    // 2atanh(z) = 2z / (1 - z^2 / (3 - 4z^2 / (5 - 9z^2 / (7 - 16z^2 / ... ) ) ) )
    //
    // By specializing the general recurrence relation for generalized continued fractions, we
    // obtain
    //
    // A_0 = 0, A_1 = 2z, A_n = (2n-1)A_(n-1) - (n-1)^2 z^2 A_(n-2) for n >= 2,
    // B_0 = 1, B_1 = 1,  B_n = (2n-1)B_(n-1) - (n-1)^2 z^2 B_(n-2) for n >= 2,
    //
    // where A_n/B_n is the nth convergent. By substituting z <- p/q and letting P_n = q^n A_n,
    // Q_n = q^n B_n, we obtain the recurrence relation
    //
    // P_0 = 0, P_1 = 2p, P_n = (2n-1)qP_(n-1) - (n-1)^2 p^2 P_(n-2) for n >= 2,
    // Q_0 = 1, Q_1 = q,  Q_n = (2n-1)qQ_(n-1) - (n-1)^2 p^2 Q_(n-2) for n >= 2.
    //
    // Then it can be shown that the denominators Q_n's are always positive and the sequence
    // (P_n/Q_n)_n of convergents converges to 2atanh(z).

    template <class Int, class UInt>
    class natural_log_calculator {
    public:
        using partial_fraction_type = frac<Int, Int>;
        using convergent_type = frac<Int, UInt>;
        using error_bound_type = frac<UInt, UInt>;

    private:
        frac<Int, Int> current_partial_fraction_{1, 0u};
        UInt p_square_ = 0u;
        UInt two_q_ = 0u;
        int current_index_ = -1;

    public:
        natural_log_calculator(frac<UInt, UInt> const& positive_rational) {
            util::constexpr_assert(!is_zero(positive_rational.denominator) &&
                                   !is_zero(positive_rational.numerator));

            // For given input x, we find the reduced form of z s.t. x = (1+z)/(1-z), i.e.,
            // z = (x-1)/(x+1). Note that z always lies in (-1,1) for x in (0,infty).
            frac<Int, UInt> z = [&] {
                convergent_generator cf{rational_continued_fraction{frac<Int, UInt>{
                    to_signed(positive_rational.numerator) - positive_rational.denominator,
                    positive_rational.numerator + positive_rational.denominator}}};

                while (!cf.is_terminated()) {
                    cf.update();
                }
                return cf.current_convergent();
            }();

            current_partial_fraction_.numerator = (z.numerator << 1);
            current_partial_fraction_.denominator = static_cast<Int>(z.denominator);
            p_square_ = abs(static_cast<Int&&>(z.numerator));
            p_square_ *= p_square_;
            two_q_ = static_cast<UInt&&>(z.denominator <<= 1);
        }

        constexpr next_partial_fraction_return<partial_fraction_type> next_partial_fraction() {
            if (current_index_ == -1) {
                ++current_index_;
                return {{1, 0u}, false};
            }
            else if (current_index_ > 0) {
                current_partial_fraction_.numerator =
                    (-current_index_ * current_index_) * p_square_;
                current_partial_fraction_.denominator += two_q_;
            }
            ++current_index_;
            return {current_partial_fraction_, false};
        }
    };
}

#endif
