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

#ifndef JKJ_HEADER_CACHING_GENERATOR
#define JKJ_HEADER_CACHING_GENERATOR

#include "continued_fraction.h"
#include <type_traits>
#include <vector>

namespace jkj {
    namespace cntfrc {
        // Memorizes all previous partial fractions, convergents, and intervals in a container, and
        // reuse them when rewinded.
        template <class ContinuedFractionGenerator>
        class caching_generator {
        public:
            using impl_type = typename std::remove_cvref_t<ContinuedFractionGenerator>::impl_type;
            using decay_type = caching_generator<
                typename std::remove_cvref_t<ContinuedFractionGenerator>::decay_type>;
            using partial_fraction_type =
                typename std::remove_cvref_t<ContinuedFractionGenerator>::partial_fraction_type;
            using convergent_type =
                typename std::remove_cvref_t<ContinuedFractionGenerator>::convergent_type;
            using interval_type =
                typename std::remove_cvref_t<ContinuedFractionGenerator>::interval_type;

            template <class OtherGenerator>
            friend class caching_generator;

        private:
            template <class T>
            struct dummy {
                template <class... U>
                constexpr dummy(U&&...) noexcept {}
            };

            struct record_partial_fraction {
                partial_fraction_type partial_fraction;
            };
            static constexpr bool has_partial_fraction_tracker =
                std::remove_cvref_t<ContinuedFractionGenerator>::template is_implementing_mixins<
                    partial_fraction_tracker>();
            using partial_fraction_record_fragment =
                std::conditional_t<has_partial_fraction_tracker, record_partial_fraction,
                                   dummy<record_partial_fraction>>;
            constexpr partial_fraction_record_fragment snapshot_partial_fraction() const {
                if constexpr (has_partial_fraction_tracker) {
                    return {cf_.current_partial_fraction()};
                }
                else {
                    return {};
                }
            }

            struct record_convergent {
                convergent_type convergent;
            };
            static constexpr bool has_convergent_tracker = std::remove_cvref_t<
                ContinuedFractionGenerator>::template is_implementing_mixins<convergent_tracker>();
            using convergent_record_fragment =
                std::conditional_t<has_convergent_tracker, record_convergent,
                                   dummy<record_convergent>>;
            constexpr convergent_record_fragment snapshot_convergent() const {
                if constexpr (has_convergent_tracker) {
                    return {cf_.current_convergent()};
                }
                else {
                    return {};
                }
            }

            struct record_interval {
                interval_type interval;
            };
            static constexpr bool has_interval_tracker = std::remove_cvref_t<
                ContinuedFractionGenerator>::template is_implementing_mixins<interval_tracker>();
            using interval_record_fragment =
                std::conditional_t<has_interval_tracker, record_interval, dummy<record_interval>>;
            constexpr interval_record_fragment snapshot_interval() const {
                if constexpr (has_interval_tracker) {
                    return {cf_.current_interval()};
                }
                else {
                    return {};
                }
            }

            ContinuedFractionGenerator cf_;
            int current_index_ = -1;
            bool terminated_ = false;

            struct record_t : partial_fraction_record_fragment,
                              convergent_record_fragment,
                              interval_record_fragment {};
            std::vector<record_t> record_;

            constexpr record_t snapshot() const {
                return {
                    {snapshot_partial_fraction()}, {snapshot_convergent()}, {snapshot_interval()}};
            }


            constexpr void initialize_record() {
                // previous_previous_convergent_ = {1, 0u};
                // previous_convergent_ = {0, 1u};
                // current_convergent_ = {1, 0u};
                record_.push_back(record_t{{snapshot_partial_fraction()},
                                           {convergent_type{1, 0u}},
                                           {snapshot_interval()}});
                record_.push_back(record_t{{snapshot_partial_fraction()},
                                           {convergent_type{0, 1u}},
                                           {snapshot_interval()}});
                record_.push_back(record_t{{snapshot_partial_fraction()},
                                           {convergent_type{1, 0u}},
                                           {snapshot_interval()}});
            }

            struct decay_copy_key {};
            template <class OtherGenerator>
            constexpr caching_generator(decay_copy_key,
                                        caching_generator<OtherGenerator> const& other)
                : cf_{other.cf_.copy()}, current_index_{other.current_index_},
                  terminated_{other.terminated_}, record_(other.record_) {}

        public:
            // caching_generator implements:
            // - index_tracker,
            // - partial_fraction_tracker if cf_ implements partial_fraction_tracker,
            // - convergent_tracker and previous_previous_convergent_tracker if cf_ implements
            //   convergent_tracker, and
            // - interval_tracker if cf_ implements interval_tracker.
            template <template <class, class> class... QueriedMixins>
            static constexpr bool is_implementing_mixins() noexcept {
                using list = tmp::join<
                    tmp::typelist<detail::mixin_type_wrapper<index_tracker>>,
                    std::conditional_t<
                        has_partial_fraction_tracker,
                        tmp::typelist<detail::mixin_type_wrapper<partial_fraction_tracker>>,
                        tmp::typelist<>>,
                    std::conditional_t<has_convergent_tracker,
                                       tmp::typelist<detail::mixin_type_wrapper<convergent_tracker>,
                                                     detail::mixin_type_wrapper<
                                                         previous_previous_convergent_tracker>>,
                                       tmp::typelist<>>,
                    std::conditional_t<has_interval_tracker,
                                       tmp::typelist<detail::mixin_type_wrapper<interval_tracker>>,
                                       tmp::typelist<>>>;
                return (... && tmp::is_in<detail::mixin_type_wrapper<QueriedMixins>>(list{}));
            }

            constexpr caching_generator(ContinuedFractionGenerator&& impl)
                : cf_{static_cast<ContinuedFractionGenerator&&>(impl)} {
                initialize_record();
            }

            // Make an identical deep copy whose internal generator is a deep copy of the one used
            // by this. The returned copy has independent state from this, but has the same record.
            // In particular, calling rewind() on both the copy and this will make them to produce
            // the same output.
            constexpr decay_type copy() const { return decay_type{*this}; }

            constexpr int current_index() const noexcept { return current_index_; }

            constexpr partial_fraction_type const& current_partial_fraction() const noexcept
                requires has_partial_fraction_tracker
            {
                return record_[std::size_t(current_index_ + 3)].partial_fraction;
            }

            constexpr convergent_type const& current_convergent() const noexcept
                requires has_convergent_tracker
            {
                return record_[std::size_t(current_index_ + 3)].convergent;
            }
            constexpr auto const& current_convergent_numerator() const noexcept
                requires has_convergent_tracker
            {
                return current_convergent().numerator;
            }
            constexpr auto const& current_convergent_denominator() const noexcept
                requires has_convergent_tracker
            {
                return current_convergent().denominator;
            }

            constexpr convergent_type const& previous_convergent() const noexcept
                requires has_convergent_tracker
            {
                return record_[std::size_t(current_index_ + 2)].convergent;
            }
            constexpr auto const& previous_convergent_numerator() const noexcept
                requires has_convergent_tracker
            {
                return previous_convergent().numerator;
            }
            constexpr auto const& previous_convergent_denominator() const noexcept
                requires has_convergent_tracker
            {
                return previous_convergent().denominator;
            }

            constexpr convergent_type const& previous_previous_convergent() const noexcept
                requires has_convergent_tracker
            {
                return record_[std::size_t(current_index_ + 1)].convergent;
            }
            constexpr auto const& previous_previous_convergent_numerator() const noexcept
                requires has_convergent_tracker
            {
                return previous_previous_convergent().numerator;
            }
            constexpr auto const& previous_previous_convergent_denominator() const noexcept
                requires has_convergent_tracker
            {
                return previous_previous_convergent().denominator;
            }

            constexpr auto const& current_interval() const noexcept
                requires has_interval_tracker
            {
                // current_interval() can change after termination.
                if (terminated()) {
                    return cf_.current_interval();
                }
                else {
                    return record_[std::size_t(current_index_ + 3)].interval;
                }
            }

            // Returns true if there are further partial fractions.
            constexpr bool proceed_to_next_partial_fraction() {
                // Progress the implementation only if the indices match.
                if (std::size_t(current_index_ + 4) == record_.size()) {
                    if (!terminated_) {
                        record_.reserve(record_.size() + 1);
                        terminated_ = !cf_.proceed_to_next_partial_fraction();
                        if (!terminated_) {
                            record_.push_back(snapshot());
                            ++current_index_;
                        }
                    }
                    return !terminated_;
                }
                else {
                    ++current_index_;
                    return true;
                }
            }

            constexpr bool terminated() const noexcept {
                return terminated_ && std::size_t(current_index_ + 4) == record_.size();
            }

            // Go back to the initial state.
            constexpr void rewind() noexcept { current_index_ = -1; }
        };

        // A convenient factory function.
        template <class Generator>
        constexpr auto make_caching_generator(Generator&& cf) {
            return caching_generator<std::remove_cvref_t<Generator>>{static_cast<Generator&&>(cf)};
        }
    }
}

#endif
