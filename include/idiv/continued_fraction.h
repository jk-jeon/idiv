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

#ifndef JKJ_HEADER_CONTINUED_FRACTION
#define JKJ_HEADER_CONTINUED_FRACTION

#include "util.h"

namespace jkj {
    // An interface for generalized continued fraction calculator for real numbers.
    // Partial numerators and denominators can be signed integers, but the denominator of the
    // convergents are assumed to be always nonnegative.
    // Given continued fraction expansion
    //
    // b0 + a1 / (b1 + a2 / (b2 + a3 / (b3 + ... ) ) ),
    //
    // we call an/bn the nth "partial fraction", and also we call the fraction obtained by
    // truncating the continued fraction at the nth partial fraction as the nth "convergent".

    template <class PartialFractionType>
    struct next_partial_fraction_return {
        PartialFractionType partial_fraction;
        bool is_last;
    };
    template <class PartialFractionType>
    next_partial_fraction_return(PartialFractionType&&, bool)
        -> next_partial_fraction_return<std::remove_cvref_t<PartialFractionType>>;

    struct unity {
        unity() = default;

        template <class T>
        constexpr unity(T&&) noexcept {}

        template <class T>
        constexpr auto&& operator*(T&& x) const noexcept {
            return static_cast<T&&>(x);
        }
        template <class T>
        friend constexpr auto&& operator*(T&& x, unity) noexcept {
            return static_cast<T&&>(x);
        }
        template <class T>
        friend constexpr T& operator*=(T& x, unity) noexcept {
            return x;
        }
        template <class T>
        friend constexpr auto&& operator/(T&& x, unity) noexcept {
            return static_cast<T&&>(x);
        }
        template <class T>
        friend constexpr T& operator/=(T& x, unity) noexcept {
            return x;
        }
    };

    template <class Impl>
    class convergent_generator_base {
    protected:
        template <class ConvergentType, class PartialFractionType>
        constexpr auto compute_next_convergent(PartialFractionType const& partial_fraction) {
            auto new_numerator = partial_fraction.denominator * current_convergent_numerator() +
                                 partial_fraction.numerator * previous_convergent_numerator();
            auto new_denominator = partial_fraction.denominator * current_convergent_denominator() +
                                   partial_fraction.numerator * previous_convergent_denominator();

            util::constexpr_assert(is_nonnegative(new_denominator));
            return ConvergentType{static_cast<decltype(new_numerator)&&>(new_numerator),
                                  abs(static_cast<decltype(new_denominator)&&>(new_denominator))};
        }

    public:
        constexpr auto const& current_convergent() const noexcept {
            return static_cast<Impl const&>(*this).current_convergent();
        }
        constexpr auto const& current_convergent_numerator() const noexcept {
            return current_convergent().numerator;
        }
        constexpr auto const& current_convergent_denominator() const noexcept {
            return current_convergent().denominator;
        }

        constexpr auto const& previous_convergent() const noexcept {
            return static_cast<Impl const&>(*this).previous_convergent();
        }
        constexpr auto const& previous_convergent_numerator() const noexcept {
            return previous_convergent().numerator;
        }
        constexpr auto const& previous_convergent_denominator() const noexcept {
            return previous_convergent().denominator;
        }
    };

    template <class Impl>
    class convergent_generator : public convergent_generator_base<convergent_generator<Impl>> {
    public:
        using partial_fraction_type = typename std::remove_cvref_t<Impl>::partial_fraction_type;
        using convergent_type = typename std::remove_cvref_t<Impl>::convergent_type;

    private:
        using crtp_base = convergent_generator_base<convergent_generator<Impl>>;

        Impl impl_;
        int current_index_ = -1;
        bool is_terminated_ = false;
        convergent_type current_convergent_{1, 0u};
        convergent_type previous_convergent_{0, 1u};

    public:
        constexpr convergent_generator(Impl&& impl) : impl_{static_cast<Impl&&>(impl)} {}

        constexpr int current_index() const noexcept { return current_index_; }
        constexpr bool is_terminated() const noexcept { return is_terminated_; }

        constexpr convergent_type const& current_convergent() const noexcept {
            return current_convergent_;
        }
        constexpr convergent_type const& previous_convergent() const noexcept {
            return previous_convergent_;
        }

        // Returns true if there are further partial fractions.
        constexpr bool update() {
            if (!is_terminated()) {
                auto result = impl_.next_partial_fraction();

                if (result.is_last) {
                    is_terminated_ = true;
                }
                auto new_convergent = crtp_base::template compute_next_convergent<convergent_type>(
                    result.partial_fraction);
                previous_convergent_ = static_cast<convergent_type&&>(current_convergent_);
                current_convergent_ = static_cast<convergent_type&&>(new_convergent);

                ++current_index_;
            }
            return !is_terminated();
        }
    };
}

#endif
