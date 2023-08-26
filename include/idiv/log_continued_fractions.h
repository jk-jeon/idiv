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

#ifndef JKJ_HEADER_LOG_CONTINUED_FRACTIONS
#define JKJ_HEADER_LOG_CONTINUED_FRACTIONS

#include "rational_continued_fractions.h"
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
    // Then it can be shown that
    //
    // P_(n-1)/Q_(n-1) - P_n/Q_n = A_(n-1)/B_(n-1) - A_n/B_n = 2 * (-1)^n p^(2n-1)/(Q_(n-1)Q_n).
    //
    // In particular, since Q_n >= (2n-1)qQ_(n-1) >= ... >= ((2n-1)(2n-3) ... 1) * q^n, the sequence
    // (P_n/Q_n)_n converges, in the alternating, "zig-zag" manner, uniformly for all p/q provided
    // |p/q|<1. Therefore, for any n, the limit value 2tanh(z) must lie in between P_(n-1)/Q_(n-1)
    // and P_n/Q_n and the error |P_n/Q_n - 2tanh(z)| is bounded by the difference
    // |P_(n-1)/Q_(n-1) - P_n/Q_n|.

    template <class Int, class UInt>
    class natural_log_calculator {
        frac<Int, UInt> current_convergent_{1, 0u};
        frac<Int, UInt> previous_convergent_{0, 1u};
        frac<UInt, UInt> current_error_bound_{1, 0u};
        UInt p_square_ = 0u;
        UInt two_q_ = 0u;
        UInt coeff_curr_ = 0u;
        UInt coeff_prev_ = 0u;

        int current_index_ = -1;

        // Compute the new convergent.
        void update() {
            coeff_curr_ += two_q_;
            coeff_prev_ = (unsigned(current_index_) * unsigned(current_index_)) * p_square_;
            auto new_numerator = coeff_curr_ * current_convergent_.numerator -
                                 coeff_prev_ * previous_convergent_.numerator;

            auto new_denominator = coeff_curr_ * current_convergent_.denominator -
                                   coeff_prev_ * previous_convergent_.denominator;

            // Update the error bound = 2|p|^(2n-1)((n-1)!)^2 / (Q_(n-1)Q_n).
            current_error_bound_.numerator *= (current_index_ * current_index_);
            current_error_bound_.denominator /= previous_convergent_.denominator;
            current_error_bound_.denominator *= new_denominator;

            // Update convergents.
            previous_convergent_ = static_cast<frac<Int, UInt>&&>(current_convergent_);
            current_convergent_.numerator = static_cast<Int&&>(new_numerator);
            current_convergent_.denominator = static_cast<UInt&&>(new_denominator);

            ++current_index_;
        }

    public:
        using int_type = Int;
        using uint_type = UInt;

        natural_log_calculator(frac<UInt, UInt> const& positive_rational) {
            util::constexpr_assert<util::error_msgs::no_error_msg>(
                !is_zero(positive_rational.denominator) && !is_zero(positive_rational.numerator));

            // For given input x, we find the reduced form of z s.t. x = (1+z)/(1-z), i.e.,
            // z = (x-1)/(x+1). Note that z always lies in (-1,1) for x in (0,infty).
            frac<Int, UInt> z = [&] {
                rational_continued_fractions<Int, UInt> cf{frac<Int, UInt>{
                    to_signed(positive_rational.numerator) - positive_rational.denominator,
                    positive_rational.numerator + positive_rational.denominator}};

                while (!cf.is_terminated()) {
                    cf.update();
                }
                return cf.current_convergent();
            }();

            // Set n = 1.
            previous_convergent_ = frac<Int, UInt>{0, 1u};
            current_convergent_ = frac<Int, UInt>{(z.numerator << 1), z.denominator};

            // error_bound = 2|p|^(2n-1)((n-1)!)^2 / (Q_(n-1)Q_n) where z = p/q.
            // When n = 1, this is equal to |2z| = |2p/q|.
            current_error_bound_ = frac<UInt, UInt>{abs(current_convergent_.numerator),
                                                    current_convergent_.denominator};

            coeff_curr_ = z.denominator;
            p_square_ = abs(static_cast<Int&&>(z.numerator));
            p_square_ *= p_square_;
            two_q_ = static_cast<UInt&&>(z.denominator <<= 1);

            current_index_ = 1;
        }

        frac<Int, UInt> const& current_value() const noexcept { return current_convergent_; }
        frac<UInt, UInt> const& current_error_bound() const noexcept {
            return current_error_bound_;
        }

        // Refine the current value so that the maximum possible error is strictly less than the
        // given bound.
        template <class ErrorValue>
        frac<Int, UInt> const& compute_within_error(ErrorValue const& error_bound) {
            while (current_error_bound() >= error_bound) {
                update();
            }
            return current_value();
        }
    };
}

#endif
