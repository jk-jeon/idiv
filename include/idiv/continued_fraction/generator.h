// Copyright 2022-2024 Junekey Jeon
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

#ifndef JKJ_HEADER_IDIV_CONTINUED_FRACTION_GENERATOR
#define JKJ_HEADER_IDIV_CONTINUED_FRACTION_GENERATOR

#include "../interval.h"
#include "../tmp.h"
#include "projective_rational.h"
#include <iterator>

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

    // The main interface type is the generator class template. It takes a type parameter Engine and
    // a list of type parameters Mixins. A mixin is a small individual feature that the user wants
    // to "mix-into" the generator instance. There are several predefined mixin classes, but users
    // can add their owns as well. The Engine parameter specifies "the engine class", the class with
    // actual continued fraction implementation, which will be instantiated as a subobject of
    // the generator class. The generator class only provides const access to this subobject.
    //
    // Upon construction, the generator class interacts with the engine class and computes the 0th
    // partial fraction. By keep calling proceed_to_next_partial_fraction() member function of the
    // generator class, it moves to the next partial fraction one-by-one, until the expansion ends.
    // If proceed_to_next_partial_fraction() fails to produce a new partial fraction (as the
    // expansion reaches the end), it returns false. Otherwise, it returns true. Or, the user can
    // call terminated() member function of the generator to check if the expansion is known to be
    // ended, i.e., whether or not the last call to proceed_to_next_partial_fraction() produced a
    // new partial fraction. Calling proceed_to_next_partial_fraction() when a previous call to it
    // returned false has no side effect, and it just simply returns false again. By convention,
    // infinity is considered to be the unique number with empty continued fraction expansion. In
    // this case, terminated() always returns true and calling proceed_to_next_partial_fraction()
    // has no effect.
    //
    // Alternatively, there are begin() and end() member functions which makes the generator class
    // an input range. Note that incrementing any iterator obtained by begin() permanently modifies
    // the generator. It is not a forward range, and generally the iteration can be done only once.
    //
    // The way newly computed partial fractions are reported back to the caller is through callback.
    // The caller can provide possibly multiple callback arguments to
    // proceed_to_next_partial_fraction(). When the generator succeeds in computing a new partial
    // fraction, it will call on_next_partial_fraction() member function of the callback parameters,
    // one-by-one, from the first callback parameter to the last callback parameter, in order,
    // passing the newly computed partial fraction as the argument.
    //
    // Callback is also the way of communication between the generator and the engine. When the user
    // calls proceed_to_next_partial_fraction() from the generator, the generator then calls
    // with_next_partial_fraction() member function of its engine, with a callback argument. The
    // engine is supposed to call on_next_partial_fraction() member function of this callback
    // parameter, which then will eventually call the on_next_partial_fraction() member functions of
    // the callback paramters the user supplied to proceed_to_next_partial_fraction(). If the engine
    // does not call on_next_partial_fraction() member function of the callback parameter passed
    // into with_next_partial_fraction() until it returns, the generator considers it as the end of
    // the continued fraction expansion, thus it sets the termination flag. Note that this callback
    // parameter of with_next_partial_fraction() member function must be a template paramter (so
    // with_next_partial_fraction() must be a template function) because the engine does not know
    // about the generator.
    //
    // By default, the generator does not do any state-keeping other than the termination flag and
    // what the engine internally does. However, the user can add such a feature by specifying
    // mixins. For instance, partial_fraction_tracker remembers the lastly computed partial fraction
    // (with the convention that the (-1)st partial fraction is 0). Once partial_fraction_tracker is
    // specified, now the generator has a member function current_partial_fraction() which retrieves
    // this information. Similarly, convergent_tracker remembers two lastly computed convergents
    // (with the convention that the (-1)st convergent is infinity and the (-2)nd convergent is 0).
    // Specifying convergent_tracker results in the generator to have current_convergent() and
    // previous_convergent() member functions (together with several others small helper functions),
    // which retreives this information.
    //
    // A mixin is consisting of 5 parts: proxy mixin, callback mixin, generator mixin, tracking data
    // and engine facade.
    //
    // The proxy is a lightweight reference object connected to the generator that exposes some
    // selected state information of the generator. The proxy is constructed in two ways: by
    // dereferencing the iterator, or by calling current_state() member function of the callback
    // parameter passed into with_next_partial_fraction(). By default, terminated() is the only
    // accessible member function of the proxy, but each mixin can add whatever additional member
    // functions into the proxy by specifying them in its proxy mixin. The proxy class of the
    // generator then inherts from all the proxy mixin classes of the specified mixins.
    //
    // Similarly, the callback object passed into with_next_partial_fraction() member function of
    // the engine class inherits from all the callback mixins of the specified mixins. In this way
    // mixins can add different kinds of communication channels between the generator and the
    // engine. Finally, the generator class itself inherits from the generator mixin classes of the
    // specified mixins. For some mixins (index_tracker, partial_fraction_tracker,
    // convergent_tracker, previous_previous_convergent_tracker), the proxy mixin and the generator
    // mixin are the same classes, which means the same information is provided either through the
    // proxy or directly from the generator. But this is not mandatory and the mixin author can
    // choose what to do.
    //
    // The tracking data of a mixin is the class that holds the actual state information that the
    // mixin cares about. For instance, the tracking data for the convergent_tracker mixin stores
    // the current and the previous convergent. The generator contains all of the tracking data
    // associated to the specified mixins as its subobjects. Whenever on_next_partial_fraction()
    // member function of the callback passed to with_next_partial_fraction() is called,
    // on_next_partial_fraction() member function of each of the tracking data subobjects is called.
    // On the other hand, if on_next_partial_fraction() is never called until
    // with_next_partial_fraction() returns, then final_update() member function of each of the
    // tracking data subobjects is called.
    //
    // Finally, the engine facade for a mixin is an interface between the engine and the mixin. The
    // key point here is that obtaining the facade is the only way for a mixin to get a mutable
    // access to the engine.
    //
    // Some mixins may properly function only when some other mixins coexist inside the same
    // generator. Such a dependency, if needed, is supposed to be specified by specializing a class
    // template mixin_traits. All the mixins required by the specified mixins are automatically and
    // transitively added to the list of mixins when the generator is instantiated, even if not
    // explicltly requested by the user. Some mixins may have ordering constraints; that is, to
    // correctly function they may require certain constraints on the order of the calls to
    // on_next_partial_fraction()/final_update() member functions of the tracking data. Such a
    // constraint is supposed to be also specified by specializing the template mixin_traits.
    // For instance, mixin_traits for previous_previous_convergent_tracker specifies that
    // convergent_tracker is required whenever previous_previous_convergent_tracker is used, and
    // on_next_partial_fraction() of the tracking data associated to
    // previous_previous_convergent_tracker must be called after that of convergent_tracker
    // returned.
    //
    // Also, some engine may require some mixins to exist in the same generator. Such a dependency
    // can be also specified if needed by defining a member type alias required_mixins in the
    // implementation class. Any mixins specified there are automatically (and transitively) added
    // to the list of mixins when the generator is instantiated, even if not explicitly requested by
    // the user. The implementation class can also specify an additional ordering constraint on
    // mixins by defining a member type alias mixin_ordering_constraints. When the generator is
    // instantiated, all of these ordering constraints are taken into account and it automatically
    // sets up the proper order of calling on_next_partial_fraction() members of the tracking data.
    // Note that the ordering constraint should not form a cycle. In such a case, the compilation
    // will fail with a diagnostic.

    namespace cntfrc {
        namespace detail {
            template <class T>
            struct decay_type_of_impl {
                using type = std::remove_cvref_t<T>;
            };

            template <class T>
                requires requires { typename T::decay_type; }
            struct decay_type_of_impl<T> {
                using type = typename T::decay_type;
            };
        }

        template <class T>
        using decay_type_of = typename detail::decay_type_of_impl<T>::type;

        // Obtain proxy/callback/generator mixin class template specializations from a mixin class
        // Mixin. When Mixin does not provide proxy/callback/generator mixins, then the
        // specialization of the following template is returned instead.
        template <class Mixin>
        struct empty_mixin {};

        namespace detail {
            template <class Mixin, class Engine, class Proxy>
            struct proxy_mixin_of_impl;
            template <class Mixin, class Engine, class Proxy>
                requires requires { typename Mixin::template proxy_mixin<Engine, Proxy>; }
            struct proxy_mixin_of_impl<Mixin, Engine, Proxy> {
                using type = typename Mixin::template proxy_mixin<Engine, Proxy>;
            };
            template <class Mixin, class Engine, class Proxy>
                requires(!requires { typename Mixin::template proxy_mixin<Engine, Proxy>; })
            struct proxy_mixin_of_impl<Mixin, Engine, Proxy> {
                using type = empty_mixin<Mixin>;
            };

            template <class Mixin, class Engine, class Callback>
            struct callback_mixin_of_impl;
            template <class Mixin, class Engine, class Callback>
                requires requires { typename Mixin::template callback_mixin<Engine, Callback>; }
            struct callback_mixin_of_impl<Mixin, Engine, Callback> {
                using type = typename Mixin::template callback_mixin<Engine, Callback>;
            };
            template <class Mixin, class Engine, class Callback>
                requires(!requires { typename Mixin::template callback_mixin<Engine, Callback>; })
            struct callback_mixin_of_impl<Mixin, Engine, Callback> {
                using type = empty_mixin<Mixin>;
            };

            template <class Mixin, class Engine, class Generator>
            struct generator_mixin_of_impl;
            template <class Mixin, class Engine, class Generator>
                requires requires { typename Mixin::template generator_mixin<Engine, Generator>; }
            struct generator_mixin_of_impl<Mixin, Engine, Generator> {
                using type = typename Mixin::template generator_mixin<Engine, Generator>;
            };
            template <class Mixin, class Engine, class Generator>
                requires(!requires { typename Mixin::template generator_mixin<Engine, Generator>; })
            struct generator_mixin_of_impl<Mixin, Engine, Generator> {
                using type = empty_mixin<Mixin>;
            };
        }
        template <class Mixin, class Engine, class Proxy>
        using proxy_mixin_of = typename detail::proxy_mixin_of_impl<Mixin, Engine, Proxy>::type;
        template <class Mixin, class Engine, class Callback>
        using callback_mixin_of =
            typename detail::callback_mixin_of_impl<Mixin, Engine, Callback>::type;
        template <class Mixin, class Engine, class Generator>
        using generator_mixin_of =
            typename detail::generator_mixin_of_impl<Mixin, Engine, Generator>::type;

        // Inspect if std::remove_cvref_t<T> contains specified mixins.
        // T can be either proxy, callback, or generator.
        template <class T, class... QueriedMixins>
        constexpr bool has_mixins() noexcept {
            return tmp::is_contained_in(tmp::typelist<QueriedMixins...>{},
                                        typename std::remove_cvref_t<T>::included_mixins{});
        }
        template <class T, class... QueriedMixins>
        constexpr bool has_mixins(tmp::typelist<QueriedMixins...>) noexcept {
            return has_mixins<T, QueriedMixins...>();
        }

        enum class callback_type_tag_t { normal, advancer };

        template <class... Callbacks>
        struct callback_ref_chain;

        template <>
        struct callback_ref_chain<> {
            template <class PartialFraction>
            constexpr void on_next_partial_fraction(PartialFraction const&) noexcept {}

            constexpr void final_update() noexcept {}

            template <class Callback>
            constexpr callback_ref_chain<std::remove_reference_t<Callback>>
            chain_front(Callback&& callback) const noexcept {
                return {callback, *this};
            }
        };

        template <class FirstCallback, class... RemainingCallbacks>
        struct callback_ref_chain<FirstCallback, RemainingCallbacks...>
            : callback_ref_chain<RemainingCallbacks...> {
        private:
            FirstCallback& first_;

        public:
            constexpr callback_ref_chain(
                FirstCallback& first, callback_ref_chain<RemainingCallbacks...> remaining) noexcept
                : callback_ref_chain<RemainingCallbacks...>{remaining}, first_{first} {}

            constexpr callback_ref_chain(FirstCallback& first,
                                         RemainingCallbacks&... remaining) noexcept
                : callback_ref_chain<RemainingCallbacks...>{remaining...}, first_{first} {}

            constexpr callback_ref_chain(callback_ref_chain const&) noexcept = default;

            template <class PartialFraction>
            constexpr void on_next_partial_fraction(PartialFraction const& next_partial_fraction) {
                if constexpr (requires {
                                  first_.on_next_partial_fraction(next_partial_fraction);
                              }) {
                    first_.on_next_partial_fraction(next_partial_fraction);
                }
                callback_ref_chain<RemainingCallbacks...>::on_next_partial_fraction(
                    next_partial_fraction);
            }

            constexpr void final_update() {
                if constexpr (requires { first_.final_update(); }) {
                    first_.final_update();
                }
                callback_ref_chain<RemainingCallbacks...>::final_update();
            }

            template <class Callback>
            constexpr callback_ref_chain<std::remove_reference_t<Callback>, FirstCallback,
                                         RemainingCallbacks...>
            chain_front(Callback&& callback) const noexcept {
                return {callback, *this};
            }
        };

        using empty_callback = callback_ref_chain<>;

        // Use this class if the mixin's tracking data does not have any data to track.
        template <class Mixin, class Engine>
        class empty_tracking_data {
            using partial_fraction_type = typename Engine::partial_fraction_type;

        public:
            constexpr empty_tracking_data(Engine const&) noexcept {}
        };

        struct empty_facade {};

        template <class Mixin, class Engine>
        struct mixin_traits {
            using tracking_data = typename Mixin::template default_tracking_data<Engine>;
        };

        namespace detail {
            template <class Mixin, class Engine>
            struct tracking_data_of_impl;
            template <class Mixin, class Engine>
                requires requires { typename mixin_traits<Mixin, Engine>::tracking_data; }
            struct tracking_data_of_impl<Mixin, Engine> {
                using type = typename mixin_traits<Mixin, Engine>::tracking_data;
            };
            template <class Mixin, class Engine>
                requires(
                    !requires { typename mixin_traits<Mixin, Engine>::tracking_data; } &&
                    requires { typename Mixin::template default_tracking_data<Engine>; })
            struct tracking_data_of_impl<Mixin, Engine> {
                using type = typename Mixin::template default_tracking_data<Engine>;
            };
            template <class Mixin, class Engine>
                requires(
                    !requires { typename mixin_traits<Mixin, Engine>::tracking_data; } &&
                    !requires { typename Mixin::template default_tracking_data<Engine>; })
            struct tracking_data_of_impl<Mixin, Engine> {
                using type = empty_tracking_data<Mixin, Engine>;
            };

            template <class Mixin, class Engine>
            struct facade_of_impl;
            template <class Mixin, class Engine>
                requires requires { typename mixin_traits<Mixin, Engine>::facade; }
            struct facade_of_impl<Mixin, Engine> {
                using type = typename mixin_traits<Mixin, Engine>::facade;
            };
            template <class Mixin, class Engine>
                requires(
                    !requires { typename mixin_traits<Mixin, Engine>::facade; } &&
                    requires { typename Mixin::template default_facade<Engine>; })
            struct facade_of_impl<Mixin, Engine> {
                using type = typename Mixin::template default_facade<Engine>;
            };
            template <class Mixin, class Engine>
                requires(
                    !requires { typename mixin_traits<Mixin, Engine>::facade; } &&
                    !requires { typename Mixin::template default_facade<Engine>; })
            struct facade_of_impl<Mixin, Engine> {
                using type = empty_facade;
            };
        }
        template <class Mixin, class Engine>
        using tracking_data_of = typename detail::tracking_data_of_impl<Mixin, Engine>::type;
        template <class Mixin, class Engine>
        using facade_of = typename detail::facade_of_impl<Mixin, Engine>::type;

        // Obtain the tracking data instance associated to Mixin from one of the mixins included in
        // Derived. Derived must be one of proxy/callback/generator, and Base must be one of
        // proxy_mixin/callback_mixin/generator_mixin included in Derived. This helper function can
        // be used to access the tracking data associated to the mixin trying to access the data, or
        // the one associated to other mixins included in Derived.
        template <class Mixin>
        struct tracking_data_accessor {
            template <class Derived, class Base>
            static constexpr decltype(auto) tracking_data(Base&& base) noexcept {
                static_assert(has_mixins<Derived, Mixin>(),
                              "the requested mixin is not present in Derived");

                static_assert(
                    tmp::is_in<std::remove_cvref_t<Base>>(
                        tmp::map<typename Derived::included_mixins, proxy_mixin_of,
                                 typename Derived::engine_type, Derived>{}) ||
                        tmp::is_in<std::remove_cvref_t<Base>>(
                            tmp::map<typename Derived::included_mixins, callback_mixin_of,
                                     typename Derived::engine_type, Derived>{}) ||
                        tmp::is_in<std::remove_cvref_t<Base>>(
                            tmp::map<typename Derived::included_mixins, generator_mixin_of,
                                     typename Derived::engine_type, Derived>{}),
                    "the type trying to access the tracking data is not the "
                    "proxy_mixin/callback_mixin/generator_mixin associated to one of the "
                    "mixins included in Derived");


                return static_cast<tmp::forward_cvref<Base&&, Derived>>(base)
                    .template tracking_data<Mixin>();
            }
        };

        struct generator_accessor {
            template <class Derived, class Base>
            static constexpr decltype(auto) generator(Base&& base) noexcept {
                static_assert(
                    tmp::is_in<std::remove_cvref_t<Base>>(
                        tmp::map<typename Derived::included_mixins, proxy_mixin_of,
                                 typename Derived::engine_type, Derived>{}) ||
                        tmp::is_in<std::remove_cvref_t<Base>>(
                            tmp::map<typename Derived::included_mixins, callback_mixin_of,
                                     typename Derived::engine_type, Derived>{}) ||
                        tmp::is_in<std::remove_cvref_t<Base>>(
                            tmp::map<typename Derived::included_mixins, generator_mixin_of,
                                     typename Derived::engine_type, Derived>{}),
                    "the type trying to obtain the generator object is not the "
                    "proxy_mixin/callback_mixin/generator_mixin associated to one of the "
                    "mixins included in Derived");

                return static_cast<tmp::forward_cvref<Base&&, Derived>>(base).generator();
            }

            template <class Mixin, class Generator, class GeneratorMixin>
            static constexpr decltype(auto) engine_facade(GeneratorMixin&& base) noexcept {
                static_assert(has_mixins<Generator, Mixin>(),
                              "the requested mixin is not present in Generator");

                static_assert(
                    std::is_same_v<
                        std::remove_cvref_t<GeneratorMixin>,
                        generator_mixin_of<Mixin, typename Generator::engine_type, Generator>>,
                    "the type trying to access the engine facade is not the "
                    "generator_mixin associated to Mixin");

                return static_cast<tmp::forward_cvref<GeneratorMixin&&, Generator>>(base)
                    .template engine_facade<Mixin>();
            }

            template <class Generator, class GeneratorMixin, class... Callbacks>
            static constexpr auto create_advancer(GeneratorMixin&& base,
                                                  Callbacks&&... callbacks) noexcept {
                static_assert(
                    tmp::is_in<std::remove_cvref_t<GeneratorMixin>>(
                        tmp::map<typename Generator::included_mixins, generator_mixin_of,
                                 typename Generator::engine_type, Generator>{}),
                    "the type trying to access the engine facade is not the "
                    "generator_mixin associated to one of the mixins included in Generator");

                return static_cast<Generator&>(base).create_advancer(
                    static_cast<Callbacks&&>(callbacks)...);
            }

            template <class Generator, class GeneratorMixin>
            static constexpr void reset_termination_flag(GeneratorMixin&& base) noexcept {
                static_assert(
                    tmp::is_in<std::remove_cvref_t<GeneratorMixin>>(
                        tmp::map<typename Generator::included_mixins, generator_mixin_of,
                                 typename Generator::engine_type, Generator>{}),
                    "the type trying to access the engine facade is not the "
                    "generator_mixin associated to one of the mixins included in Generator");

                return static_cast<Generator&>(base).reset_termination_flag();
            }
        };

        namespace detail {
            // A direct edge from Before to After representing a mixin ordering constraint.
            template <class Before, class After>
            struct mixin_ordering_constraint_edge {};

            // First -> Second -> ... -> Last
            template <class ConstraintList, class First>
            constexpr ConstraintList
            linear_chain_mixin_ordering_constraint_impl(ConstraintList,
                                                        std::type_identity<First>) noexcept {
                return {};
            }
            template <class ConstraintList, class First, class Second, class... Remaining>
            constexpr auto
            linear_chain_mixin_ordering_constraint_impl(ConstraintList, std::type_identity<First>,
                                                        std::type_identity<Second>,
                                                        std::type_identity<Remaining>...) noexcept {
                return linear_chain_mixin_ordering_constraint_impl(
                    tmp::push_back<ConstraintList, mixin_ordering_constraint_edge<First, Second>>{},
                    std::type_identity<Second>{}, std::type_identity<Remaining>{}...);
            }
        }

        namespace mixin_ordering_constraint {
            // Before::on_next_partial_fraction() should be called before
            // After::on_next_partial_fraction().
            template <class Before, class After>
            using before_after =
                tmp::typelist<detail::mixin_ordering_constraint_edge<Before, After>>;

            // First::on_next_partial_fraction() and then Second::on_next_partial_fraction() and
            // then Remaining::on_next_partial_fraction()...
            template <class First, class Second, class... Remaining>
            using linear_chain = decltype(detail::linear_chain_mixin_ordering_constraint_impl(
                tmp::typelist<>{}, std::type_identity<First>{}, std::type_identity<Second>{},
                std::type_identity<Remaining>{}...));

            // E.g. mixin_ordering_constraints =
            //                       constraint_list<before_after<index_tracker,
            //                                                    interval_estimate_provider>,
            //                                       linear_chain<partial_fraction_tracker,
            //                                                    convergent_tracaker,
            //                                                    interval_estimate_provider>>;
            template <class... Constraints>
            using constraint_list = tmp::join<Constraints...>;
        }

        namespace detail {
            // Topological sort for the order dependency graph of mixins.
            template <std::size_t node_count>
            struct topological_sort_output {
                bool succeed = false;
                util::array<std::size_t, node_count> sorted_indices{};
            };
            struct graph_edge {
                std::size_t from_idx;
                std::size_t to_idx;
            };

            // Ignore any edge from/to a node with index >= node_count.
            // Such an edge can present if an ordering constraint is specified for mixins that are
            // not actually specified/required.
            template <std::size_t node_count, std::size_t edge_count>
            constexpr topological_sort_output<node_count>
            topological_sort(util::array<graph_edge, edge_count> const& edges) noexcept {
                topological_sort_output<node_count> result;

                if constexpr (node_count == 0) {
                    result.succeed = true;
                    return result;
                }
                else if constexpr (edge_count == 0) {
                    for (std::size_t i = 0; i < node_count; ++i) {
                        result.sorted_indices[i] = i;
                    }
                    result.succeed = true;
                    return result;
                }
                else {
                    std::size_t sorted_node_count = 0;

                    util::array<bool, node_count> nonfree_node_marks{};
                    util::array<bool, edge_count> deleted_edge_marks{};
                    for (std::size_t i = 0; i < edge_count; ++i) {
                        if (edges[i].from_idx >= node_count || edges[i].to_idx >= node_count) {
                            deleted_edge_marks[i] = true;
                        }
                    }

                    // Find all nodes with no incoming edge (called free nodes) and put them into a
                    // stack.
                    for (std::size_t i = 0; i < edge_count; ++i) {
                        if (!deleted_edge_marks[i]) {
                            nonfree_node_marks[edges[i].to_idx] = true;
                        }
                    }
                    util::array<std::size_t, node_count> free_node_stack{};
                    std::size_t free_node_count = 0;
                    for (std::size_t i = 0; i < node_count; ++i) {
                        if (nonfree_node_marks[i]) {
                            continue;
                        }
                        free_node_stack[free_node_count++] = i;
                    }

                    // Loop until there is no more free nodes.
                    while (free_node_count > 0) {
                        auto const free_node_idx = free_node_stack[--free_node_count];
                        result.sorted_indices[sorted_node_count++] = free_node_idx;

                        // Loop through all edges starting from the chosen free node.
                        for (std::size_t edge_idx = 0; edge_idx < edge_count; ++edge_idx) {
                            if (deleted_edge_marks[edge_idx]) {
                                continue;
                            }
                            if (edges[edge_idx].from_idx != free_node_idx) {
                                continue;
                            }

                            // Delete the chosen edge, and if this is the last incoming edge into
                            // the target node, then put that node into the free node stack.
                            deleted_edge_marks[edge_idx] = true;
                            bool last_incoming_edge = true;
                            for (std::size_t i = 0; i < edge_count; ++i) {
                                if (deleted_edge_marks[i]) {
                                    continue;
                                }
                                if (edges[i].to_idx == edges[edge_idx].to_idx) {
                                    last_incoming_edge = false;
                                    break;
                                }
                            }

                            if (last_incoming_edge) {
                                free_node_stack[free_node_count++] = edges[edge_idx].to_idx;
                            }
                        } // Loop over outcoming edges.
                    }     // Loop over free nodes.

                    // Sort is successful if there is no remaining edge. Otherwise, there is a
                    // cycle.
                    for (std::size_t i = 0; i < edge_count; ++i) {
                        if (!deleted_edge_marks[i]) {
                            // Fail.
                            return result;
                        }
                    }
                    // Success.
                    result.succeed = true;
                    return result;
                }
            }

            template <class T>
            struct get_required_mixins;
            template <class T>
                requires requires { typename T::required_mixins; }
            struct get_required_mixins<T> {
                using type = typename T::required_mixins;
            };
            template <class T>
                requires(!requires { typename T::required_mixins; })
            struct get_required_mixins<T> {
                using type = tmp::typelist<>;
            };

            template <class T>
            struct get_before_than;
            template <class T>
                requires requires { typename T::before_than; }
            struct get_before_than<T> {
                using type = typename T::before_than;
            };
            template <class T>
                requires(!requires { typename T::before_than; })
            struct get_before_than<T> {
                using type = tmp::typelist<>;
            };

            template <class T>
            struct get_after_than;
            template <class T>
                requires requires { typename T::after_than; }
            struct get_after_than<T> {
                using type = typename T::after_than;
            };
            template <class T>
                requires(!requires { typename T::after_than; })
            struct get_after_than<T> {
                using type = tmp::typelist<>;
            };

            template <class T>
            struct get_mixin_ordering_constraints;
            template <class T>
                requires requires { typename T::mixin_ordering_constraints; }
            struct get_mixin_ordering_constraints<T> {
                using type = typename T::mixin_ordering_constraints;
            };
            template <class T>
                requires(!requires { typename T::mixin_ordering_constraints; })
            struct get_mixin_ordering_constraints<T> {
                using type = tmp::typelist<>;
            };

            // Construct an array of graph_edge from a list of mixin_ordering_constraint_edge's.
            template <class... Mixins, class... Constraints>
            constexpr auto
            convert_local_mixin_ordering_constraints_list(tmp::typelist<Mixins...>,
                                                          tmp::typelist<Constraints...>) noexcept {
                if constexpr (sizeof...(Constraints) == 0) {
                    return util::array<graph_edge, 0>{};
                }
                else {
                    using list = tmp::typelist<Mixins...>;
                    auto constraint_to_graph_edge =
                        []<class Before, class After>(
                            mixin_ordering_constraint_edge<Before, After>) {
                            return graph_edge{tmp::find_first_index<Before>(list{}),
                                              tmp::find_first_index<After>(list{})};
                        };
                    return util::array<graph_edge, sizeof...(Constraints)>{
                        {constraint_to_graph_edge(Constraints{})...}};
                }
            }

            // Construct an array of graph_edge from the ordering constraints of TargetMixin
            // specified in it's traits specialization.
            template <class TargetMixin, class... Mixins, class... Before, class... After>
            constexpr auto convert_single_mixin_traits(tmp::typelist<Mixins...>,
                                                       tmp::typelist<Before...>,
                                                       tmp::typelist<After...>) noexcept {
                constexpr std::size_t total_size = sizeof...(Before) + sizeof...(After);

                if constexpr (total_size == 0) {
                    return util::array<graph_edge, 0>{};
                }
                else {
                    using list = tmp::typelist<Mixins...>;
                    constexpr std::size_t idx_of_target =
                        tmp::find_first_index<TargetMixin>(list{});

                    return util::array<graph_edge, total_size>{
                        {{idx_of_target, tmp::find_first_index<Before>(list{})}...,
                         {tmp::find_first_index<After>(list{}), idx_of_target}...}};
                }
            }

            template <std::size_t... sizes>
            constexpr auto merge_graphs(util::array<graph_edge, sizes> const&... graphs) noexcept {
                auto impl = [](auto& global_index, auto& ret_value, auto&& graph) {
                    if constexpr (std::remove_cvref_t<decltype(graph)>::size() != 0) {
                        for (std::size_t idx = 0; idx < graph.size(); ++idx, ++global_index) {
                            ret_value[global_index] = graph[idx];
                        }
                    }
                };
                std::size_t global_index = 0;
                util::array<graph_edge, (std::size_t(0) + ... + sizes)> ret_value{};
                (impl(global_index, ret_value, graphs), ...);
                return ret_value;
            }

            // Construct an array of graph_edge from the ordering constraints of all mixins in
            // dMixins specified in their mixin_traits specializations.
            template <class Engine, class... Mixins>
            constexpr auto convert_multiple_mixin_traits(tmp::typelist<Mixins...> list) noexcept {
                return merge_graphs(convert_single_mixin_traits<Mixins>(
                    list, typename get_before_than<mixin_traits<Mixins, Engine>>::type{},
                    typename get_after_than<mixin_traits<Mixins, Engine>>::type{})...);
            }

            // Find the transitive closure of required mixins.
            template <class Engine, class... Mixins>
            constexpr auto get_transitive_required_mixin_list(tmp::typelist<Mixins...>) noexcept {
                using original_list = tmp::typelist<Mixins...>;
                using augmented_list = tmp::remove_duplicate<
                    tmp::join<original_list, decltype(typename get_required_mixins<
                                                      mixin_traits<Mixins, Engine>>::type{})...>>;

                if constexpr (std::is_same_v<augmented_list, original_list>) {
                    return augmented_list{};
                }
                else {
                    return get_transitive_required_mixin_list<Engine>(augmented_list{});
                }
            }

            // Main compile-time function for computing the mixin list.
            template <class Engine, class LocalConstraintList, class... RequiredMixins,
                      class... AdditionalMixins>
            constexpr auto
            find_sorted_required_mixin_list_impl(tmp::typelist<RequiredMixins...>,
                                                 tmp::typelist<AdditionalMixins...>) noexcept {
                // Collect all required mixins transitively.
                using mixin_list = decltype(get_transitive_required_mixin_list<Engine>(
                    tmp::remove_duplicate<tmp::typelist<RequiredMixins..., AdditionalMixins...>>{}

                    ));

                // Collect all ordering constraints and convert them into an array of graph_edge.
                constexpr auto mixin_dependency_graph =
                    merge_graphs(convert_local_mixin_ordering_constraints_list(
                                     mixin_list{}, LocalConstraintList{}),
                                 convert_multiple_mixin_traits<Engine>(mixin_list{}));

                // Get a topologically sorted array of mixin indices pointing into
                // mixin_list. Any ordering constraint pointing from/to a mixin not
                // included in mixin_list is ignored.
                constexpr auto sorted_mixin_indices =
                    topological_sort<mixin_list::size>(mixin_dependency_graph);

                static_assert(sorted_mixin_indices.succeed,
                              "mixin's ordering constraint should not form a cycle");

                // Convert the index array into a tmp::typelist and return.
                if constexpr (mixin_list::size == 0) {
                    return tmp::typelist<>{};
                }
                else {
                    return [&sorted_mixin_indices]<std::size_t... I>(std::index_sequence<I...>) {
                        return tmp::typelist<
                            tmp::get_type<sorted_mixin_indices.sorted_indices[I], mixin_list>...>{};
                    }(std::make_index_sequence<mixin_list::size>{});
                }
            }
        }

        // Find the transitive closure of required mixins of the engine and the provided list of
        // mixins, without duplicate, topologically sorted the list according to the ordering
        // constraints given between the mixins, and then return the resulting typelist.
        template <class Engine, class... AdditionalMixins>
        constexpr auto find_sorted_required_mixin_list() noexcept {
            return detail::find_sorted_required_mixin_list_impl<
                Engine, typename detail::get_mixin_ordering_constraints<Engine>::type>(
                typename detail::get_required_mixins<Engine>::type{},
                tmp::typelist<AdditionalMixins...>{});
        }

        namespace detail {
            template <class Engine>
            class generator_base {
                Engine engine_;

            public:
                constexpr generator_base(Engine&& engine)
                    : engine_{static_cast<Engine&&>(engine)} {}

                constexpr Engine& engine() & noexcept { return engine_; }
                constexpr Engine const& engine() const& noexcept { return engine_; }
                constexpr Engine&& engine() && noexcept { return std::move(*this).engine_; }
            };

            template <class Mixin, class Engine>
            class tracking_data_wrapper {
                using tracking_data_type = tracking_data_of<Mixin, Engine>;
                tracking_data_type data_;

                template <class, class>
                friend class tracking_data_wrapper;

            public:
                constexpr tracking_data_wrapper(Engine const& engine) : data_{engine} {}

                template <class OtherEngine>
                constexpr tracking_data_wrapper(
                    tracking_data_wrapper<Mixin, OtherEngine> const& other)
                    : data_{other.data_} {}

                constexpr tracking_data_type& get() & noexcept { return data_; }
                constexpr tracking_data_type const& get() const& noexcept { return data_; }
                constexpr tracking_data_type&& get() && noexcept { return std::move(*this).data_; }
            };

            // The actual implementation of the generator class.
            template <class Engine, class... Mixins>
                requires(std::is_object_v<Engine>)
            class generator_impl
                : private generator_base<Engine>,
                  private tracking_data_wrapper<Mixins, Engine>...,
                  public generator_mixin_of<Mixins, Engine, generator_impl<Engine, Mixins...>>... {
                bool terminated_ = false;

                template <class Mixin>
                using tracking_data_type = tracking_data_of<Mixin, Engine>;

            public:
                using engine_type = Engine;
                using included_mixins = tmp::typelist<Mixins...>;
                using decay_type = generator_impl<decay_type_of<Engine>, Mixins...>;
                using partial_fraction_type = typename Engine::partial_fraction_type;

                template <class... Callbacks>
                explicit constexpr generator_impl(Engine engine_in, Callbacks&&... callbacks)
                    : generator_base<Engine>{static_cast<Engine&&>(engine_in)},
                      tracking_data_wrapper<Mixins, Engine>{engine()}... {
                    proceed_to_zeroth_partial_fraction(callbacks...);
                }

                template <class OtherEngine>
                    requires(std::is_same_v<generator_impl, decay_type> &&
                             !std::is_same_v<Engine, OtherEngine> &&
                             std::is_constructible_v<Engine, OtherEngine>)
                explicit constexpr generator_impl(
                    generator_impl<OtherEngine, Mixins...> const& other)
                    : generator_base<Engine>{other.engine()},
                      tracking_data_wrapper<Mixins, Engine>{
                          static_cast<tracking_data_wrapper<Mixins, OtherEngine> const&>(
                              other)}... {}

                // Copy/move constructor/assignements are allowed only for value-like generators.
                constexpr generator_impl(generator_impl const&)
                    requires(std::is_same_v<generator_impl, decay_type>)
                = default;
                constexpr generator_impl(generator_impl&&)
                    requires(std::is_same_v<generator_impl, decay_type>)
                = default;
                constexpr generator_impl& operator=(generator_impl const&)
                    requires(std::is_same_v<generator_impl, decay_type>)
                = default;
                constexpr generator_impl& operator=(generator_impl&&)
                    requires(std::is_same_v<generator_impl, decay_type>)
                = default;

                constexpr generator_impl(generator_impl const&)
                    requires(!std::is_same_v<generator_impl, decay_type>)
                = delete;
                constexpr generator_impl(generator_impl&&)
                    requires(!std::is_same_v<generator_impl, decay_type>)
                = delete;
                constexpr generator_impl& operator=(generator_impl const&)
                    requires(!std::is_same_v<generator_impl, decay_type>)
                = delete;
                constexpr generator_impl& operator=(generator_impl&&)
                    requires(!std::is_same_v<generator_impl, decay_type>)
                = delete;

                // Make a copy of the generator. The internal continued fraction engine is copied in
                // value even if it has a reference-like part.
                constexpr decay_type copy() const { return decay_type{*this}; }

                constexpr bool terminated() const noexcept { return terminated_; }

                constexpr Engine const& engine() const noexcept {
                    return generator_base<Engine>::engine();
                }

                // Returns true if succeeded obtaining a further partial fraction.
                template <class... Callbacks>
                constexpr bool proceed_to_next_partial_fraction(Callbacks&&... callbacks) {
                    if (!terminated_) {
                        callback_type<callback_type_tag_t::normal,
                                      std::remove_reference_t<Callbacks>...>
                            callback_arg_for_engine{
                                *this, callback_ref_chain<std::remove_reference_t<Callbacks>...>{
                                           callbacks...}};

                        terminated_ = true;
                        modifiable_engine().with_next_partial_fraction(callback_arg_for_engine);

                        if (terminated_) {
                            final_update(callback_ref_chain<Callbacks...>{callbacks...});
                        }
                    }
                    return !terminated_;
                }

                class iterator;
                class sentinel;
                class pointer;

                class proxy_type final : util::noncopyable<proxy_type>,
                                         util::nonmovable<proxy_type>,
                                         public proxy_mixin_of<Mixins, Engine, proxy_type>... {
                    generator_impl const& gen_;

                    template <class Mixin>
                    friend struct cntfrc::tracking_data_accessor;

                    friend generator_accessor;
                    friend generator_impl;
                    friend iterator;
                    friend pointer;

                    explicit constexpr proxy_type(generator_impl const& gen) noexcept : gen_{gen} {}

                    template <class Mixin>
                    constexpr tracking_data_type<Mixin> const& tracking_data() const noexcept {
                        return gen_.tracking_data<Mixin>();
                    }

                    constexpr generator_impl const& generator() const noexcept { return gen_; }

                public:
                    using engine_type = Engine;
                    using included_mixins = tmp::typelist<Mixins...>;

                    constexpr bool terminated() const noexcept { return gen_.terminated(); }
                };

                class sentinel final {
                    friend generator_impl;

                    constexpr sentinel() noexcept = default;
                };

                class iterator final {
                    generator_impl* gen_ptr_; // non-null.

                    friend generator_impl;
                    friend sentinel;

                    explicit constexpr iterator(generator_impl& gen) : gen_ptr_{&gen} {}

                public:
                    // std::weakly_incrementable.
                    using difference_type = std::ptrdiff_t;
                    constexpr iterator& operator++() {
                        gen_ptr_->proceed_to_next_partial_fraction();
                        return *this;
                    }
                    constexpr iterator operator++(int) {
                        iterator ret_value = *this;
                        ++*this;
                        return ret_value;
                    }

                    // std::input_or_output_iterator.
                    constexpr proxy_type operator*() const noexcept {
                        return proxy_type{*gen_ptr_};
                    }

                    // std::indirectly_readable.
                    using value_type = proxy_type;

                    // std::input_iterator.
                    using iterator_concept = std::input_iterator_tag;

                    constexpr pointer operator->() const noexcept { return pointer{gen_ptr_}; }

                    friend constexpr bool operator==(iterator const& it, sentinel) noexcept {
                        return it.gen_ptr_->terminated();
                    }
                };

                class pointer final : util::noncopyable<>, util::nonmovable<> {
                    proxy_type proxy_;

                    friend iterator;

                    explicit constexpr pointer(generator_impl* gen_ptr) : proxy_{*gen_ptr} {}

                public:
                    constexpr proxy_type const* operator->() const noexcept { return &proxy_; }
                };

                constexpr iterator begin() noexcept { return iterator{*this}; }
                constexpr sentinel end() const noexcept { return sentinel{}; }

            private:
                template <class Mixin>
                friend struct cntfrc::tracking_data_accessor;

                friend struct cntfrc::generator_accessor;

                template <class Mixin>
                constexpr tracking_data_type<Mixin>& tracking_data() noexcept {
                    return static_cast<tracking_data_wrapper<Mixin, Engine>&>(*this).get();
                }
                template <class Mixin>
                constexpr tracking_data_type<Mixin> const& tracking_data() const noexcept {
                    return static_cast<tracking_data_wrapper<Mixin, Engine> const&>(*this).get();
                }

                constexpr engine_type& modifiable_engine() noexcept {
                    return generator_base<Engine>::engine();
                }

                template <class Mixin>
                constexpr decltype(auto) engine_facade() & {
                    return facade_of<Mixin, Engine>{generator_base<Engine>::engine()};
                }
                template <class Mixin>
                constexpr decltype(auto) engine_facade() const& {
                    return facade_of<Mixin, Engine>{generator_base<Engine>::engine()};
                }
                template <class Mixin>
                constexpr decltype(auto) engine_facade() && {
                    return facade_of<Mixin, Engine>{
                        static_cast<generator_base<Engine>&&>(*this).engine()};
                }

                constexpr generator_impl const& generator() const noexcept { return *this; }
                constexpr proxy_type current_state() const noexcept { return proxy_type{*this}; }
                constexpr void reset_termination_flag() noexcept { terminated_ = false; }

                template <callback_type_tag_t callback_type_tag, class... ChainedCallbacks>
                class callback_type
                    : util::noncopyable<callback_type<callback_type_tag, ChainedCallbacks...>>,
                      util::nonmovable<callback_type<callback_type_tag, ChainedCallbacks...>>,
                      public callback_mixin_of<
                          Mixins, Engine,
                          callback_type<callback_type_tag, ChainedCallbacks...>>... {
                private:
                    generator_impl& gen_;
                    callback_ref_chain<ChainedCallbacks...> chained_callbacks_;

                    template <class Mixin>
                    friend struct cntfrc::tracking_data_accessor;

                    friend struct cntfrc::generator_accessor;

                    friend generator_impl;

                    template <callback_type_tag_t, class...>
                    friend class callback_type;

                    constexpr callback_type(
                        generator_impl& gen,
                        callback_ref_chain<ChainedCallbacks...> chained_callbacks) noexcept
                        : gen_{gen}, chained_callbacks_{chained_callbacks} {}

                    // Provides non-const access.
                    template <class Mixin>
                    constexpr tracking_data_type<Mixin>& tracking_data() const noexcept {
                        return static_cast<tracking_data_wrapper<Mixin, Engine>&>(gen_).get();
                    }

                    constexpr generator_impl const& generator() const noexcept { return gen_; }

                public:
                    using engine_type = Engine;
                    using included_mixins = tmp::typelist<Mixins...>;

                    constexpr void
                    on_next_partial_fraction(partial_fraction_type const& next_partial_fraction) {
                        gen_.terminated_ = false;
                        gen_.on_next_partial_fraction(chained_callbacks_, next_partial_fraction);
                    }

                    constexpr proxy_type current_state() const noexcept {
                        return gen_.current_state();
                    }

                    template <class Callback>
                    constexpr callback_type<callback_type_tag, std::remove_reference_t<Callback>,
                                            ChainedCallbacks...>
                    chain_front(Callback&& callback) noexcept {
                        return {gen_, chained_callbacks_.chain_front(callback)};
                    }

                    constexpr void final_update()
                        requires(callback_type_tag == callback_type_tag_t::advancer)
                    {
                        gen_.terminated_ = true;
                        gen_.final_update(chained_callbacks_);
                    }

                    constexpr bool proceed_to_next_partial_fraction()
                        requires(callback_type_tag == callback_type_tag_t::advancer)
                    {
                        return gen_.proceed_to_next_partial_fraction(chained_callbacks_);
                    }
                };

                template <class... Callbacks>
                constexpr callback_type<callback_type_tag_t::advancer,
                                        std::remove_reference_t<Callbacks>...>
                create_advancer(Callbacks&&... callbacks) {
                    return {*this, callback_ref_chain<std::remove_reference_t<Callbacks>...>(
                                       callbacks...)};
                }

                // Update the generator with a new partial fraction obtained.
                template <class Callback>
                constexpr void
                on_next_partial_fraction(Callback&& callback,
                                         partial_fraction_type const& next_partial_fraction) {
                    auto invoke_on_next_partial_fraction =
                        [this, &next_partial_fraction](auto&& mixin_tracking_data) {
                            if constexpr (requires {
                                              mixin_tracking_data.on_next_partial_fraction(
                                                  generator(), next_partial_fraction);
                                          }) {
                                mixin_tracking_data.on_next_partial_fraction(generator(),
                                                                             next_partial_fraction);
                            }
                        };
                    (invoke_on_next_partial_fraction(tracking_data<Mixins>()), ...);

                    callback.on_next_partial_fraction(next_partial_fraction);
                }

                // Finalize the generator when the continued fraction expansion reached the end.
                template <class Callback>
                constexpr void final_update(Callback&& callback) {
                    auto invoke_final_update = [this](auto&& mixin_tracking_data) {
                        if constexpr (requires { mixin_tracking_data.final_update(generator()); }) {
                            mixin_tracking_data.final_update(generator());
                        }
                    };
                    (invoke_final_update(tracking_data<Mixins>()), ...);

                    if constexpr (requires { callback.final_update(); }) {
                        callback.final_update();
                    }
                }

                // Proceed to the 0th partial fraction.
                // If the engine provides a special routine for the 0th partial fraction, then
                // call it. Otherwise, call with_next_partial_fraction().
                template <class... Callbacks>
                constexpr void proceed_to_zeroth_partial_fraction(Callbacks&&... callbacks) {
                    callback_type<callback_type_tag_t::normal,
                                  std::remove_reference_t<Callbacks>...>
                        callback_arg_for_engine{
                            *this, callback_ref_chain<std::remove_reference_t<Callbacks>...>{
                                       callbacks...}};

                    terminated_ = true;
                    if constexpr (requires {
                                      modifiable_engine().with_zeroth_partial_fraction(
                                          callback_arg_for_engine);
                                  }) {
                        modifiable_engine().with_zeroth_partial_fraction(callback_arg_for_engine);
                    }
                    else {
                        modifiable_engine().with_next_partial_fraction(callback_arg_for_engine);
                    }
                    if (terminated_) {
                        final_update(callback_ref_chain<Callbacks...>{callbacks...});
                    }
                }
            };

            template <class Engine, class MixinList>
            struct generator_type;

            template <class Engine, class... Mixins>
            struct generator_type<Engine, tmp::typelist<Mixins...>> {
                using type = generator_impl<Engine, Mixins...>;
            };
        }

        // The main interface type.
        // Automatically add all transitively required mixins into the list of mixins, and sort
        // the list according to the imposed ordering constraints. Then generator derives from
        // the resulting list of mixins.
        // Engine must be a non-reference type, but it may be a reference-like type or may have a
        // reference-like subobject. This is because a separate instance of Engine is never
        // considered a "complete object", rather only when packaged into a generator it becomes
        // such.
        template <class Engine, class... Mixins>
        using generator = typename detail::generator_type<
            Engine, decltype(find_sorted_required_mixin_list<Engine, Mixins...>())>::type;

        // A convenient factory function.
        template <class... Mixins, class Engine>
        constexpr auto make_generator(Engine&& engine) {
            return generator<std::remove_cvref_t<Engine>, Mixins...>{static_cast<Engine&&>(engine)};
        }

        // Mixin: stores the current index of the continued fraction expansion. The index is
        // retrievable from the proxy as well as the generator.
        struct index_tracker {
        private:
            template <class Engine, class Derived>
            struct mixin_impl {
                constexpr int current_index() const {
                    return tracking_data_accessor<index_tracker>::template tracking_data<Derived>(
                               *this)
                        .current_index(generator_accessor::generator<Derived>(*this));
                }
            };

            class default_tracking_data_impl {
                int current_index_ = -1;

            public:
                template <class Engine>
                constexpr default_tracking_data_impl(Engine const&) noexcept {}

                template <class Generator, class PartialFractionType>
                constexpr void on_next_partial_fraction(Generator const&,
                                                        PartialFractionType const&) noexcept {
                    ++current_index_;
                }

                template <class Generator>
                constexpr int current_index(Generator const&) const noexcept {
                    return current_index_;
                }
            };

        public:
            template <class Engine, class Proxy>
            using proxy_mixin = mixin_impl<Engine, Proxy>;

            template <class Engine, class Generator>
            using generator_mixin = mixin_impl<Engine, Generator>;

            template <class Engine>
            using default_tracking_data = default_tracking_data_impl;
        };

        // Mixin: stores the current partial fraction of the continued fraction expansion. The
        // stored partial fraction can be retrieved only from the proxy.
        struct partial_fraction_tracker {
        private:
            template <class Engine, class Derived>
            struct mixin_impl {
                using partial_fraction_type = typename Engine::partial_fraction_type;

                constexpr partial_fraction_type current_partial_fraction() const {
                    return tracking_data_accessor<partial_fraction_tracker>::template tracking_data<
                               Derived>(*this)
                        .current_partial_fraction(generator_accessor::generator<Derived>(*this));
                }
            };

            template <class PartialFractionType>
            class default_tracking_data_impl {
                PartialFractionType current_partial_fraction_{unity{}, zero{}};

            public:
                template <class Engine>
                constexpr default_tracking_data_impl(Engine const&) {}

                template <class Generator>
                constexpr void
                on_next_partial_fraction(Generator const&,
                                         PartialFractionType const& next_partial_fraction) {
                    current_partial_fraction_ = next_partial_fraction;
                }

                template <class Generator>
                constexpr PartialFractionType const&
                current_partial_fraction(Generator const&) const noexcept {
                    return current_partial_fraction_;
                }
            };

        public:
            template <class Engine, class Proxy>
            using proxy_mixin = mixin_impl<Engine, Proxy>;

            template <class Engine, class Generator>
            using generator_mixin = mixin_impl<Engine, Generator>;

            template <class Engine>
            using default_tracking_data =
                default_tracking_data_impl<typename Engine::partial_fraction_type>;
        };

        // Mixin: stores the previous and the current convergents of the continued fraction
        // expansion. The stored convergents are retrievable from the proxy as well as the
        // generator.
        struct convergent_tracker {
        private:
            template <class Engine, class Derived>
            struct mixin_impl {
                static_assert(
                    requires {
                        typename mixin_traits<convergent_tracker, Engine>::convergent_type;
                    }, "the continued fraction engine should provide the type for storing "
                       "convergents through mixin_traits");
                using convergent_type =
                    typename mixin_traits<convergent_tracker, Engine>::convergent_type;

                constexpr convergent_type const& current_convergent() const noexcept {
                    return tracking_data_accessor<convergent_tracker>::template tracking_data<
                               Derived>(*this)
                        .current_convergent(generator_accessor::generator<Derived>(*this));
                }
                constexpr auto const& current_convergent_numerator() const noexcept {
                    return current_convergent().numerator;
                }
                constexpr auto const& current_convergent_denominator() const noexcept {
                    return current_convergent().denominator;
                }

                constexpr convergent_type const& previous_convergent() const noexcept {
                    return tracking_data_accessor<convergent_tracker>::template tracking_data<
                               Derived>(*this)
                        .previous_convergent(generator_accessor::generator<Derived>(*this));
                }
                constexpr auto const& previous_convergent_numerator() const noexcept {
                    return previous_convergent().numerator;
                }
                constexpr auto const& previous_convergent_denominator() const noexcept {
                    return previous_convergent().denominator;
                }
            };

            template <class ConvergentType>
            class default_tracking_data_impl {
                ConvergentType current_convergent_{unity{}, zero{}};
                ConvergentType previous_convergent_{zero{}, unity{}};

            public:
                template <class Engine>
                constexpr default_tracking_data_impl(Engine const&) {}

                template <class Generator, class PartialFractionType>
                constexpr void
                on_next_partial_fraction(Generator const&,
                                         PartialFractionType const& next_partial_fraction) {
                    auto next_numerator =
                        next_partial_fraction.denominator * current_convergent_.numerator +
                        next_partial_fraction.numerator * previous_convergent_.numerator;
                    auto next_denominator =
                        next_partial_fraction.denominator * current_convergent_.denominator +
                        next_partial_fraction.numerator * previous_convergent_.denominator;

                    previous_convergent_ = std::move(current_convergent_);
                    current_convergent_ = ConvergentType{std::move(next_numerator),
                                                         util::abs(std::move(next_denominator))};
                }

                template <class Generator>
                constexpr ConvergentType const&
                current_convergent(Generator const&) const noexcept {
                    return current_convergent_;
                }
                template <class Generator>
                constexpr ConvergentType const&
                previous_convergent(Generator const&) const noexcept {
                    return previous_convergent_;
                }
            };

        public:
            template <class Engine, class Proxy>
            using proxy_mixin = mixin_impl<Engine, Proxy>;

            template <class Engine, class Generator>
            using generator_mixin = mixin_impl<Engine, Generator>;

            template <class Engine>
            using default_tracking_data = default_tracking_data_impl<
                typename mixin_traits<convergent_tracker, Engine>::convergent_type>;
        };
        template <class Engine>
        struct mixin_traits<convergent_tracker, Engine> {
            using convergent_type = typename Engine::convergent_type;
        };

        // Mixin: stores the previous previous convergent of the continued fraction expansion.
        // Depends on convergent_tracker's presence, so including this mixin results in tracking
        // the last three convergents. The stored convergent is retrievable only from the proxy.
        struct previous_previous_convergent_tracker {
        private:
            template <class Engine, class Derived>
            struct mixin_impl {
            private:
                using convergent_type =
                    typename mixin_traits<convergent_tracker, Engine>::convergent_type;

            public:
                constexpr convergent_type const& previous_previous_convergent() const noexcept {
                    return tracking_data_accessor<previous_previous_convergent_tracker>::
                        template tracking_data<Derived>(*this)
                            .previous_previous_convergent(
                                generator_accessor::generator<Derived>(*this));
                }
                constexpr auto const& previous_previous_convergent_numerator() const noexcept {
                    return previous_previous_convergent().numerator;
                }
                constexpr auto const& previous_previous_convergent_denominator() const noexcept {
                    return previous_previous_convergent().denominator;
                }
            };

            template <class ConvergentType>
            class default_tracking_data_impl {
                ConvergentType previous_previous_convergent_{zero{}, zero{}};

            public:
                template <class Engine>
                constexpr default_tracking_data_impl(Engine const&) {}

                template <class Generator, class PartialFractionType>
                constexpr void on_next_partial_fraction(Generator const& gen,
                                                        PartialFractionType const&) {
                    previous_previous_convergent_ = gen.previous_convergent();
                }

                template <class Generator>
                constexpr ConvergentType const&
                previous_previous_convergent(Generator const&) const noexcept {
                    return previous_previous_convergent_;
                }
            };

        public:
            template <class Engine, class Proxy>
            using proxy_mixin = mixin_impl<Engine, Proxy>;

            template <class Engine, class Generator>
            using generator_mixin = mixin_impl<Engine, Generator>;

            template <class Engine>
            using default_tracking_data = default_tracking_data_impl<
                typename mixin_traits<convergent_tracker, Engine>::convergent_type>;
        };
        template <class Engine>
        struct mixin_traits<previous_previous_convergent_tracker, Engine> {
            using required_mixins = tmp::typelist<convergent_tracker>;
            using before_than = tmp::typelist<convergent_tracker>;
        };

        // Mixin: provides interface to treat the generator as an interval estimate provider.
        struct interval_estimate_provider {
        private:
            template <class IntervalType>
            class default_tracking_data_impl {
                IntervalType current_interval_;

            public:
                template <class Engine>
                constexpr default_tracking_data_impl(Engine const& engine)
                    : current_interval_{engine.initial_interval()} {}

                template <class Generator>
                constexpr void on_next_interval(Generator const&, IntervalType next_interval) {
                    current_interval_ = static_cast<IntervalType&&>(next_interval);
                }

                template <class Generator>
                constexpr IntervalType const& current_interval(Generator const&) const noexcept {
                    return current_interval_;
                }
            };

            template <class PartialFractionType, class IntervalType>
            class default_facade_impl {
            public:
                using required_mixins = tmp::typelist<index_tracker, convergent_tracker>;

                static_assert(
                    std::is_same_v<decltype(PartialFractionType::numerator), unity>,
                    "inspecting the partial fraction type concludes that the continued "
                    "fraction is not regular; the default facade for "
                    "interval_estimate_provider only works for regular continued fractions");

                template <class Engine>
                constexpr default_facade_impl(Engine&) noexcept {}

                template <class Advancer>
                constexpr void refine_interval(Advancer&& advancer) {
                    if constexpr (IntervalType::allowed_interval_types() ==
                                  util::array<cyclic_interval_type_t, 1>{
                                      cyclic_interval_type_t::single_point}) {
                        // Nothing to do if the only possible interval type is the single point.
                    }
                    else {
                        auto&& state = advancer.current_state();
                        if (state.terminated()) {
                            // Nothing to do if there is no further continued fraction expansion.
                        }
                        else if (advancer.proceed_to_next_partial_fraction()) {
                            // Use convergents. Note that this is valid only for regular continued
                            // fractions.
                            if (state.current_index() >= 1) {
                                if (state.current_index() % 2 == 0) {
                                    advancer.on_next_interval(
                                        cyclic_interval<
                                            typename IntervalType::value_type,
                                            cyclic_interval_type_t::left_closed_right_open>{
                                            state.current_convergent(),
                                            state.previous_convergent()});
                                }
                                else {
                                    advancer.on_next_interval(
                                        cyclic_interval<
                                            typename IntervalType::value_type,
                                            cyclic_interval_type_t::left_open_right_closed>{
                                            state.previous_convergent(),
                                            state.current_convergent()});
                                }
                            }
                        }
                        else {
                            if constexpr (IntervalType::is_allowed_interval_type(
                                              cyclic_interval_type_t::single_point)) {
                                advancer.on_next_interval(
                                    cyclic_interval<typename IntervalType::value_type,
                                                    cyclic_interval_type_t::single_point>{
                                        state.current_convergent()});
                            }
                        }
                    }
                }
            };

        public:
            template <class Engine, class Callback>
            struct callback_mixin {
                static_assert(
                    requires {
                        typename mixin_traits<interval_estimate_provider, Engine>::interval_type;
                    }, "the continued fraction engine should provide the type for storing "
                       "interval estimates");

                using interval_type =
                    typename mixin_traits<interval_estimate_provider, Engine>::interval_type;

                constexpr void on_next_interval(interval_type next_interval) {
                    tracking_data_accessor<interval_estimate_provider>::template tracking_data<
                        Callback>(*this)
                        .on_next_interval(generator_accessor::generator<Callback>(*this),
                                          static_cast<interval_type&&>(next_interval));
                }
            };

            template <class Engine, class Generator>
            class generator_mixin {
            public:
                using interval_type =
                    typename mixin_traits<interval_estimate_provider, Engine>::interval_type;

                template <class... Callbacks>
                constexpr void refine_interval(Callbacks&&... callbacks) {
                    generator_accessor::engine_facade<interval_estimate_provider, Generator>(*this)
                        .refine_interval(
                            generator_accessor::create_advancer<Generator>(*this, callbacks...));
                }

                constexpr interval_type const& current_interval() const noexcept {
                    return tracking_data_accessor<
                               interval_estimate_provider>::template tracking_data<Generator>(*this)
                        .current_interval(generator_accessor::generator<Generator>(*this));
                }

                // Refine the current value so that the maximum possible error is no more than
                // the given bound.
                template <class ErrorValue>
                constexpr void refine_interval_until(ErrorValue const& error_bound) {
                    // Translate the lower bound by error_bound and see if it belongs to the
                    // interval.
                    while (current_interval().visit([&](auto&& itv) {
                        using itv_type = std::remove_cvref_t<decltype(itv)>;
                        static_assert(itv_type::interval_type() != cyclic_interval_type_t::empty);
                        if constexpr (itv_type::interval_type() ==
                                      cyclic_interval_type_t::single_point) {
                            return false;
                        }
                        else if constexpr (itv_type::interval_type() ==
                                           cyclic_interval_type_t::entire) {
                            return true;
                        }
                        else {
                            auto translated = linear_fractional_translation(
                                error_bound.numerator, error_bound.denominator)(itv.lower_bound());

                            if (!itv.contains(translated)) {
                                return false;
                            }

                            if (itv.interval_type() == cyclic_interval_type_t::closed &&
                                translated == itv.upper_bound()) {
                                return false;
                            }

                            return true;
                        }
                    })) {
                        refine_interval();
                    }
                }
            };

            template <class Engine>
            using default_tracking_data = default_tracking_data_impl<
                typename mixin_traits<interval_estimate_provider, Engine>::interval_type>;

            template <class Engine>
            using default_facade = default_facade_impl<
                typename Engine::partial_fraction_type,
                typename mixin_traits<interval_estimate_provider, Engine>::interval_type>;
        };
        template <class Engine>
        struct mixin_traits<interval_estimate_provider, Engine> {
            using interval_type = typename Engine::interval_type;
            using facade = interval_estimate_provider::default_facade<Engine>;
            using required_mixins = typename facade::required_mixins;
        };

        // Mixin: provides interface for rewinding the generator back to its initial partial
        // fraction.
        struct rewinder {
            template <class Engine, class Generator>
            class generator_mixin {
            public:
                constexpr void rewind() {
                    // Since there is no way to call rewind() before returning from the
                    // generator's constructor, current_index() == -1 means the number is equal
                    // to infinity and there is no partial fraction. In such a case, rewind()
                    // does nothing.
                    if (generator_accessor::generator<Generator>(*this).current_index() != -1) {
                        // Call rewind() for the engine.
                        generator_accessor::engine_facade<rewinder, Generator>(*this).rewind(
                            generator_accessor::generator<Generator>(*this));

                        // Call rewind() for the tracking data associated to each of the
                        // included mixins.
                        call_rewind_for_tracking_data(typename Generator::included_mixins{});

                        generator_accessor::reset_termination_flag<Generator>(*this);
                    }
                }

            private:
                constexpr void call_rewind_for_tracking_data(tmp::typelist<>) noexcept {}

                template <class FirstMixin, class... RemainingMixins>
                constexpr void
                call_rewind_for_tracking_data(tmp::typelist<FirstMixin, RemainingMixins...>) {
                    auto&& tracking_data =
                        tracking_data_accessor<FirstMixin>::template tracking_data<Generator>(
                            *this);
                    auto&& generator = generator_accessor::generator<Generator>(*this);
                    if constexpr (requires { tracking_data.rewind(generator); }) {
                        tracking_data.rewind(generator);
                    }
                    call_rewind_for_tracking_data(tmp::typelist<RemainingMixins...>{});
                }
            };

            template <class Engine>
            class default_facade {
            public:
                static_assert(mixin_traits<rewinder, Engine>::is_engine_rewindable,
                              "rewinder mixin requires the engine to provide the facade for it");

                constexpr default_facade(Engine&) noexcept {}

                template <class Generator>
                constexpr void rewind(Generator const&) noexcept {}
            };
        };
        template <class Engine>
        struct mixin_traits<rewinder, Engine> {
            static constexpr bool is_engine_rewindable = requires {
                Engine::is_engine_rewindable;
                requires(Engine::is_engine_rewindable == true);
            };
            using required_mixins = tmp::typelist<index_tracker>;
        };
    }
}

#endif
