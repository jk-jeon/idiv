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

#ifndef JKJ_HEADER_CACHING_CONTINUED_FRACTION
#define JKJ_HEADER_CACHING_CONTINUED_FRACTION

#include "continued_fraction.h"
#include <type_traits>
#include <vector>

namespace jkj {
    namespace cntfrc {
        // Memorizes all previous partial fractions and convergents in a container, and reuse them
        // when rewinded.
        template <class Impl>
        class caching_continued_fraction {
        public:
            using partial_fraction_type = typename std::remove_cvref_t<Impl>::partial_fraction_type;
            using convergent_type = typename std::remove_cvref_t<Impl>::convergent_type;

        private:
            Impl impl_;
            int current_index_ = -1;
            bool terminated_ = false;

            struct record_t {
                partial_fraction_type partial_fraction;
                convergent_type convergent;
            };
            std::vector<record_t> record_;

            constexpr void initialize_record() {
                // previous_previous_convergent_ = {1, 0u};
                // previous_convergent_ = {0, 1u};
                // current_convergent_ = {1, 0u};
                record_.push_back(record_t{partial_fraction_type{}, convergent_type{1, 0u}});
                record_.push_back(record_t{partial_fraction_type{}, convergent_type{0, 1u}});
                record_.push_back(record_t{partial_fraction_type{}, convergent_type{1, 0u}});
            }

        public:
            constexpr caching_continued_fraction(Impl&& impl) : impl_{static_cast<Impl&&>(impl)} {
                initialize_record();
            }

            constexpr int current_index() const noexcept { return current_index_; }
            constexpr bool terminated() const noexcept {
                return terminated_ && current_index_ + 4 == record_.size();
            }

            constexpr convergent_type const& current_convergent() const noexcept {
                return record_[current_index_ + 3].convergent;
            }
            constexpr auto const& current_convergent_numerator() const noexcept {
                return current_convergent().numerator;
            }
            constexpr auto const& current_convergent_denominator() const noexcept {
                return current_convergent().denominator;
            }

            constexpr convergent_type const& previous_convergent() const noexcept {
                return record_[current_index_ + 2].convergent;
            }
            constexpr auto const& previous_convergent_numerator() const noexcept {
                return previous_convergent().numerator;
            }
            constexpr auto const& previous_convergent_denominator() const noexcept {
                return previous_convergent().denominator;
            }

            constexpr convergent_type const& previous_previous_convergent() const noexcept {
                return record_[current_index_ + 1].convergent;
            }
            constexpr auto const& previous_previous_convergent_numerator() const noexcept {
                return previous_previous_convergent().numerator;
            }
            constexpr auto const& previous_previous_convergent_denominator() const noexcept {
                return previous_previous_convergent().denominator;
            }

            // Returns true if there are further partial fractions.
            constexpr bool update() {
                // Progress the implementation only if the indices match.
                if (current_index_ + 4 == record_.size()) {
                    if (!terminated_) {
                        record_.reserve(record_.size() + 1);
                        terminated_ = !impl_.update();
                        record_.push_back(
                            record_t{impl_.current_partial_fraction(), impl_.current_convergent()});
                        ++current_index_;
                    }
                    return !terminated_;
                }
                else {
                    ++current_index_;
                    return true;
                }
            }

            // Go back to the initial state.
            constexpr void rewind() noexcept { current_index_ = -1; }
        };
    }
}

#endif
