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

#ifndef JKJ_HEADER_REWINDABLE_CONTINUED_FRACTIONS
#define JKJ_HEADER_REWINDABLE_CONTINUED_FRACTIONS

#include "continued_fractions.h"
#include <type_traits>

namespace jkj {
    // Memorizes all previous convergents and coefficients in a container, and reuse them when
    // rewinded.
    template <class Impl, template <class ValueType, class...> class Container>
    class rewindable_continued_fractions {
    public:
        using int_type = typename Impl::int_type;
        using uint_type = typename Impl::uint_type;

        // Perfect forward everything to the implementation.
        template <class... Args>
        explicit constexpr rewindable_continued_fractions(Args&&... args)
            : impl_(static_cast<Args&&>(args)...) {
            initialize_record();
        }
        template <class Arg>
            requires(std::is_convertible_v<Arg &&, Impl> &&
                     !std::is_base_of_v<rewindable_continued_fractions, std::decay_t<Arg>> &&
                     !std::is_same_v<rewindable_continued_fractions, std::decay_t<Arg>>)
        explicit constexpr rewindable_continued_fractions(Arg&& arg)
            : impl_(static_cast<Arg&&>(arg)) {
            initialize_record();
        }

    private:
        int current_index_ = -1;

        struct record_t {
            int_type coefficient;
            frac<int_type, uint_type> convergent;
        };
        Container<record_t> record_;
        Impl impl_;

        constexpr void initialize_record() {
            // current_coefficient_ = 0;
            // current_convergent_ = {1, 0u};
            // previous_convergent_ = {0, 1u};
            record_.push_back(
                record_t{int_type(0), frac<int_type, uint_type>{int_type(0), uint_type(1u)}});
            record_.push_back(
                record_t{int_type(0), frac<int_type, uint_type>{int_type(1), uint_type(0u)}});
        }

    public:
        constexpr int current_index() const noexcept { return current_index_; }

        constexpr int_type const& current_coefficient() const noexcept {
            return record_[current_index_ + 2].coefficient;
        }

        constexpr frac<int_type, uint_type> const& current_convergent() const noexcept {
            return record_[current_index_ + 2].convergent;
        }

        constexpr int_type const& current_numerator() const noexcept {
            return current_convergent().numerator;
        }

        constexpr uint_type const& current_denominator() const noexcept {
            return current_convergent().denominator;
        }

        constexpr frac<int_type, uint_type> const& previous_convergent() const noexcept {
            return record_[current_index_ + 1].convergent;
        }

        constexpr int_type const& previous_numerator() const noexcept {
            return previous_convergent().numerator;
        }

        constexpr uint_type const& previous_denominator() const noexcept {
            return previous_convergent().denominator;
        }

        constexpr bool is_terminated() const noexcept {
            return impl_.is_terminated() && current_index_ == impl_.current_index();
        }

        // Do nothing if the procedure is terminated.
        // Returns true if the update is done,
        // and returns false if the procedure is already terminated.
        constexpr bool update() {
            if (!is_terminated()) {
                // Progress the implementation only if the indices match.
                if (current_index_ == impl_.current_index()) {
                    record_.reserve(record_.size() + 1);
                    impl_.update();
                    // Record the new states.
                    record_.push_back(
                        record_t{impl_.current_coefficient(), impl_.current_convergent()});
                }
                ++current_index_;

                return true;
            }
            else {
                return false;
            }
        }

        // Go back to the initial state/
        constexpr void rewind() noexcept { current_index_ = -1; }
    };
}

#endif