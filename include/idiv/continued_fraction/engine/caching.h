// Copyright 2023-2024 Junekey Jeon
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

#ifndef JKJ_HEADER_IDIV_CONTINUED_FRACTION_ENGINE_CACHING
#define JKJ_HEADER_IDIV_CONTINUED_FRACTION_ENGINE_CACHING

#include "../generator.h"
#include <vector>

namespace jkj {
    namespace cntfrc {
        namespace engine {
            template <class InternalEngine>
                requires(std::is_object_v<InternalEngine>)
            class caching {
            public:
                using partial_fraction_type = typename InternalEngine::partial_fraction_type;
                using required_mixins = tmp::typelist<index_tracker>;

                using decay_type = caching<decay_type_of<InternalEngine>>;

                static constexpr bool is_engine_rewindable = true;

            private:
                std::vector<partial_fraction_type> records_{};
                InternalEngine impl_;

                class callback_type {
                    std::vector<partial_fraction_type>& records_;

                    friend caching;

                    constexpr callback_type(std::vector<partial_fraction_type>& records)
                        : records_{records} {}

                public:
                    constexpr void
                    on_next_partial_fraction(partial_fraction_type const& next_partial_fraction) {
                        records_.push_back(next_partial_fraction);
                    }
                };

            public:
                explicit constexpr caching(InternalEngine internal_engine)
                    : impl_{static_cast<InternalEngine&&>(internal_engine)} {}

                template <class Callback>
                constexpr void with_zeroth_partial_fraction(Callback&& callback) {
                    auto callback_to_chain = callback_type{records_};
                    auto callback_arg_for_internal = callback.chain_front(callback_to_chain);
                    if constexpr (requires {
                                      impl_.with_zeroth_partial_fraction(callback_arg_for_internal);
                                  }) {
                        impl_.with_zeroth_partial_fraction(callback_arg_for_internal);
                    }
                    else {
                        impl_.with_next_partial_fraction(callback_arg_for_internal);
                    }
                }

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    util::constexpr_assert(callback.current_state().current_index() >= 0);

                    // Draw the new partial fraction only when necessary.
                    if (std::size_t(callback.current_state().current_index() + 1) ==
                        records_.size()) {
                        impl_.with_next_partial_fraction(
                            callback.chain_front(callback_type{records_}));
                    }
                    // Otherwise, use the stored partial fraction.
                    else {
                        callback.on_next_partial_fraction(
                            records_[std::size_t(callback.current_state().current_index() + 1)]);
                    }
                }

                constexpr decltype(auto) initial_interval()
                    requires requires { impl_.initial_interval(); }
                {
                    return impl_.initial_interval();
                }

                constexpr util::span<partial_fraction_type const> records() const noexcept {
                    return {records_.data(), records_.size()};
                }

                constexpr InternalEngine const& internal_engine() const noexcept { return impl_; }

                template <class Mixin>
                constexpr facade_of<Mixin, InternalEngine> create_facade_from_internal_engine() {
                    return facade_of<Mixin, InternalEngine>{impl_};
                }
            };

            template <class InternalEngine>
            caching(InternalEngine) -> caching<InternalEngine>;
        }

        // Specialize mixin_traits for mixins that engine::caching understands. Any other
        // mixin not listed here cannot be used along with engine::caching unless they
        // provide appropriate specializations of mixin_traits for engine::caching.
        template <class InternalEngine>
        struct mixin_traits<index_tracker, engine::caching<InternalEngine>> {
            class tracking_data {
                using partial_fraction_type = typename InternalEngine::partial_fraction_type;

                int current_index_ = -1;

            public:
                constexpr tracking_data(engine::caching<InternalEngine> const&) noexcept {}

                template <class Generator>
                constexpr void on_next_partial_fraction(Generator const&,
                                                        partial_fraction_type const&) noexcept {
                    ++current_index_;
                }

                template <class Generator>
                constexpr int current_index(Generator const&) const noexcept {
                    return current_index_;
                }

                template <class Generator>
                constexpr void rewind(Generator const&) noexcept {
                    util::constexpr_assert(current_index_ != -1);
                    current_index_ = 0;
                }
            };
        };

        template <class InternalEngine>
        struct mixin_traits<partial_fraction_tracker, engine::caching<InternalEngine>> {
            using required_mixins = tmp::typelist<index_tracker>;
            using after_than = tmp::typelist<index_tracker>;

            class tracking_data {
                using partial_fraction_type = typename InternalEngine::partial_fraction_type;

            public:
                constexpr tracking_data(engine::caching<InternalEngine> const&) noexcept {}

                template <class Generator>
                constexpr void on_next_partial_fraction(Generator const&,
                                                        partial_fraction_type const&) noexcept {}

                template <class Generator>
                constexpr partial_fraction_type const&
                current_partial_fraction(Generator const& gen) const noexcept {
                    auto const& records = gen.engine().records();
                    util::constexpr_assert(gen.current_index() >= 0 &&
                                           std::size_t(gen.current_index()) < records.size());
                    return records[std::size_t(gen.current_index())];
                }
            };
        };

        template <class InternalEngine>
        struct mixin_traits<convergent_tracker, engine::caching<InternalEngine>> {
            using convergent_type =
                typename mixin_traits<convergent_tracker, InternalEngine>::convergent_type;
            using required_mixins = tmp::typelist<index_tracker>;
            using after_than = tmp::typelist<index_tracker>;

            class tracking_data {
                using partial_fraction_type = typename InternalEngine::partial_fraction_type;

                // The 0th element corresponds to the (-2)nd convergent.
                std::vector<convergent_type> records_{convergent_type{zero{}, unity{}},
                                                      convergent_type{unity{}, zero{}}};

            public:
                constexpr tracking_data(engine::caching<InternalEngine> const&) noexcept {}

                template <class Generator>
                constexpr void on_next_partial_fraction(
                    Generator const& gen,
                    [[maybe_unused]] partial_fraction_type const& next_partial_fraction) {
                    // Don't do anything if the input is not new.
                    if (std::size_t(gen.current_index() + 2) < records_.size()) {
                        return;
                    }
                    // Otherwise, compute the new convergent as usual.
                    else {
                        util::constexpr_assert(gen.current_index() >= 0);
                        auto const& current = records_[std::size_t(gen.current_index() + 1)];
                        auto const& previous = records_[std::size_t(gen.current_index())];
                        auto next_numerator =
                            next_partial_fraction.denominator * current.numerator +
                            next_partial_fraction.numerator * previous.numerator;
                        auto next_denominator =
                            next_partial_fraction.denominator * current.denominator +
                            next_partial_fraction.numerator * previous.denominator;

                        records_.push_back(convergent_type{std::move(next_numerator),
                                                           util::abs(std::move(next_denominator))});
                    }
                }

                template <class Generator>
                constexpr convergent_type const&
                current_convergent(Generator const& gen) const noexcept {
                    util::constexpr_assert(gen.current_index() >= -1 &&
                                           std::size_t(gen.current_index() + 2) < records_.size());
                    return records_[std::size_t(gen.current_index() + 2)];
                }
                template <class Generator>
                constexpr convergent_type const&
                previous_convergent(Generator const& gen) const noexcept {
                    util::constexpr_assert(gen.current_index() >= -1 &&
                                           std::size_t(gen.current_index() + 2) < records_.size());
                    return records_[std::size_t(gen.current_index() + 1)];
                }

                // Open this for previous_previous_convergent_tracker.
                constexpr std::vector<convergent_type> const& records() const noexcept {
                    return records_;
                }
            };
        };

        template <class InternalEngine>
        struct mixin_traits<previous_previous_convergent_tracker, engine::caching<InternalEngine>> {
            using required_mixins = tmp::typelist<index_tracker, convergent_tracker>;
            using after_than = tmp::typelist<index_tracker, convergent_tracker>;

            class tracking_data {
                using partial_fraction_type = typename InternalEngine::partial_fraction_type;
                using convergent_type = typename InternalEngine::convergent_type;

            public:
                constexpr tracking_data(engine::caching<InternalEngine> const&) noexcept {}

                template <class Generator>
                constexpr void on_next_partial_fraction(Generator const&,
                                                        partial_fraction_type const&) {}

                template <class Generator>
                constexpr convergent_type const&
                previous_previous_convergent(Generator const& gen) const noexcept {
                    auto const& records =
                        tracking_data_accessor<convergent_tracker>::tracking_data<Generator>(*this)
                            .records();
                    util::constexpr_assert(gen.current_index() >= -1 &&
                                           std::size_t(gen.current_index()) < records.size());

                    if (gen.current_index() == -1) {
                        // (-3)rd convergent is infinity by convention, and it is equal to the
                        // (-1)st convergent, so just return the (-1)st convergent.
                        return records[1];
                    }
                    else {
                        return records[gen.current_index()];
                    }
                }
            };
        };

        template <class InternalEngine>
        struct mixin_traits<interval_estimate_provider, engine::caching<InternalEngine>> {
            using interval_type =
                typename mixin_traits<interval_estimate_provider, InternalEngine>::interval_type;

            class tracking_data {
                using partial_fraction_type = typename InternalEngine::partial_fraction_type;

                tracking_data_of<interval_estimate_provider, InternalEngine>
                    internal_tracking_data_;

            public:
                constexpr tracking_data(engine::caching<InternalEngine> const& engine)
                    : internal_tracking_data_{engine.internal_engine()} {}

                // These two functions (as well as the absent final_update()) do not need to do
                // anything because the tracking data attached to the internal generator already
                // handles all the necessary updates.
                template <class Generator>
                constexpr void on_next_partial_fraction(
                    Generator const& gen,
                    partial_fraction_type const& next_partial_fraction) noexcept {
                    if constexpr (requires {
                                      internal_tracking_data_.on_next_partial_fraction(
                                          gen, next_partial_fraction);
                                  }) {
                        internal_tracking_data_.on_next_partial_fraction(gen,
                                                                         next_partial_fraction);
                    }
                }
                template <class Generator>
                constexpr void
                on_next_interval(Generator const& gen,
                                 interval_type const& next_partial_fraction) noexcept {
                    if constexpr (requires {
                                      internal_tracking_data_.on_next_interval(
                                          gen, next_partial_fraction);
                                  }) {

                        internal_tracking_data_.on_next_interval(gen, next_partial_fraction);
                    }
                }

                template <class Generator>
                constexpr decltype(auto) current_interval(Generator const& gen) const noexcept {
                    return internal_tracking_data_.current_interval(gen);
                }
            };

            class facade {
                using partial_fraction_type = typename InternalEngine::partial_fraction_type;

                facade_of<interval_estimate_provider, InternalEngine> internal_facade_;

            public:
                using required_mixins =
                    decltype(find_sorted_required_mixin_list<InternalEngine,
                                                             interval_estimate_provider>());

                explicit constexpr facade(engine::caching<InternalEngine>& engine) noexcept
                    : internal_facade_{engine.template create_facade_from_internal_engine<
                          interval_estimate_provider>()} {}

                template <class Advancer>
                constexpr void refine_interval(Advancer&& advancer) {
                    internal_facade_.refine_interval(advancer);
                }
            };

            using required_mixins = typename facade::required_mixins;
        };

        template <class... Mixins, class Engine>
        constexpr auto make_caching_generator(Engine&& engine) {
            return generator<engine::caching<std::remove_cvref_t<Engine>>, Mixins...>{
                engine::caching<std::remove_cvref_t<Engine>>{static_cast<Engine&&>(engine)}};
        }
    }
}

#endif
