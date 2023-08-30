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
        using impl_type = Impl;

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
                auto next_convergent = crtp_base::template compute_next_convergent<convergent_type>(
                    result.partial_fraction);
                is_terminated_ = result.is_last;

                (static_cast<Mixin<Impl, convergent_generator>&>(*this).update(
                     result.partial_fraction, next_convergent, impl_),
                 ...);

                previous_convergent_ = static_cast<convergent_type&&>(current_convergent_);
                current_convergent_ = static_cast<convergent_type&&>(next_convergent);

                ++current_index_;
            }
            return !is_terminated();
        }
    };

    // Keep track of an interval where the limit value should live inside.
    template <class Impl, class ConvergentGenerator>
    class interval_tracker {
    public:
        using convergent_type = typename std::remove_cvref_t<Impl>::convergent_type;
        using interval_type = typename std::remove_cvref_t<Impl>::interval_type;
        friend ConvergentGenerator;

    private:
        interval_type current_interval_ = std::remove_cvref_t<Impl>::initial_interval();

        template <class PartialFractionType>
        constexpr void update(PartialFractionType const&, convergent_type const& next_convergent,
                              Impl& impl) {
            auto& self = static_cast<ConvergentGenerator const&>(*this);
            current_interval_ = impl.next_interval(self.previous_convergent(),
                                                   self.current_convergent(), next_convergent);
        }

    public:
        interval_type const& current_interval() const noexcept { return current_interval_; }

        auto current_error_bound() const {
            return current_interval_.upper_bound() - current_interval_.lower_bound();
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
