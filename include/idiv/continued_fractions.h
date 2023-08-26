// Copyright 2022-2023 Junekey Jeon
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

#ifndef JKJ_HEADER_CONTINUED_FRACTIONS
#define JKJ_HEADER_CONTINUED_FRACTIONS

#include "frac.h"
#include "util.h"

namespace jkj {
    // An interface for generalized continued fraction calculator for real numbers.
    // Partial numerators and denominators can be signed integers, but the denominator of the
    // convergents are assumed to be always nonnegative.
    template <class Impl, class Int, class UInt>
    class continued_fractions {
        // The (-1)st coefficient is assumed to be 0.
        Int current_coefficient_{0};
        frac<Int, UInt> current_convergent_{1, 0u};
        frac<Int, UInt> previous_convergent_{0, 1u};
        int current_index_ = -1;
        bool terminated_ = false;

    protected:
        constexpr void set_terminate_flag() noexcept { terminated_ = true; }

    public:
        using int_type = Int;
        using uint_type = UInt;

        constexpr int current_index() const noexcept { return current_index_; }

        constexpr Int const& current_coefficient() const noexcept { return current_coefficient_; }

        constexpr frac<Int, UInt> const& current_convergent() const noexcept {
            return current_convergent_;
        }

        constexpr Int const& current_numerator() const noexcept {
            return current_convergent().numerator;
        }

        constexpr UInt const& current_denominator() const noexcept {
            return current_convergent().denominator;
        }

        constexpr frac<Int, UInt> const& previous_convergent() const noexcept {
            return previous_convergent_;
        }

        constexpr Int const& previous_numerator() const noexcept {
            return previous_convergent().numerator;
        }

        constexpr UInt const& previous_denominator() const noexcept {
            return previous_convergent().denominator;
        }

        constexpr bool is_terminated() const noexcept { return terminated_; }

        // Do nothing if the procedure is terminated.
        // Returns true if the update is done,
        // and returns false if the procedure is already terminated.
        constexpr bool update() {
            if (!is_terminated()) {
                frac<Int, UInt> new_output;
                current_coefficient_ = static_cast<Impl&>(*this).compute_next_coefficient();

                auto new_numerator =
                    previous_numerator() + current_coefficient_ * current_numerator();
                auto new_denominator =
                    previous_denominator() + current_coefficient_ * current_denominator();

                util::constexpr_assert(is_nonnegative(new_denominator));

                previous_convergent_ = static_cast<frac<Int, UInt>&&>(current_convergent_);
                current_convergent_.numerator = static_cast<Int&&>(new_numerator);
                current_convergent_.denominator = abs(static_cast<Int&&>(new_denominator));

                ++current_index_;

                return true;
            }
            else {
                return false;
            }
        }
    };
}

#endif
