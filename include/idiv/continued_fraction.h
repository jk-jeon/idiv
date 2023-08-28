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

    // May use this type for the partial numerators of regular continued fractions.
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

    template <class Impl, template <class, class> class... Mixin>
    class convergent_generator
        : public convergent_generator_base<convergent_generator<Impl, Mixin...>>,
          public Mixin<Impl, convergent_generator<Impl, Mixin...>>... {
    public:
        using convergent_type = typename std::remove_cvref_t<Impl>::convergent_type;

    private:
        using crtp_base = convergent_generator_base<convergent_generator<Impl, Mixin...>>;

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
                auto next_convergent = crtp_base::template compute_next_convergent<convergent_type>(
                    result.partial_fraction);

                (static_cast<Mixin<Impl, convergent_generator>&>(*this).update(
                     result.partial_fraction, next_convergent),
                 ...);

                previous_convergent_ = static_cast<convergent_type&&>(current_convergent_);
                current_convergent_ = static_cast<convergent_type&&>(next_convergent);

                ++current_index_;
            }
            return !is_terminated();
        }
    };

    // Keep track of the quantity e = |P_(n-1)/Q_(n-1) - P_n/Q_n| = (|a1| ... |an|)/(Q_(n-1) Q_n).
    // This must be an upper bound on the distance between the current convergent P_n/Q_n and the
    // limiting value of the continued fraction, provided that
    // 1. Q_n's are known to be positive, and
    // 2. the limit exists.
    // Otherwise, e may not have anything to do with the error bound.
    template <class Impl, class ConvergentGenerator>
    class error_bound_tracker {
    public:
        using convergent_type = typename std::remove_cvref_t<Impl>::convergent_type;
        using error_bound_type = typename std::remove_cvref_t<Impl>::error_bound_type;

        friend ConvergentGenerator;

    private:
        error_bound_type current_error_bound_{1u, 0u};

        template <class PartialFractionType>
        constexpr void update(PartialFractionType const& partial_fraction,
                              convergent_type const& next_convergent) {
            decltype(auto) current_denominator =
                static_cast<ConvergentGenerator&>(*this).current_convergent_denominator();
            decltype(auto) previous_denominator =
                static_cast<ConvergentGenerator&>(*this).previous_convergent_denominator();

            // For the initial update.
            if (is_zero(current_denominator)) {
                return;
            }

            current_error_bound_.numerator *= abs(partial_fraction.numerator);
            if (is_zero(current_error_bound_.denominator)) {
                current_error_bound_.denominator = current_denominator;
            }
            else {
                current_error_bound_.denominator /= previous_denominator;
            }
            current_error_bound_.denominator *= next_convergent.denominator;
        }

    public:
        error_bound_type const& current_error_bound() const noexcept {
            return current_error_bound_;
        }

        // Refine the current value so that the maximum possible error is strictly less than the
        // given bound.
        template <class ErrorValue>
        convergent_type const& progress_until(ErrorValue const& error_bound) {
            while (current_error_bound() >= error_bound) {
                static_cast<ConvergentGenerator&>(*this).update();
            }
            return static_cast<ConvergentGenerator&>(*this).current_convergent();
        }
    };
}

#endif
