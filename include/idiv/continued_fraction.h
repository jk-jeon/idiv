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

#include "interval.h"
#include "projective_rational.h"
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

    // The main interface type is the generator class template. It takes a type parameter Impl and a
    // list of template template parameters Mixins. A mixin is a small individual feature that the
    // user wants to "mix-into" the generator instance. The generator gets the specified feature by
    // deriving from the corresponding mixin type. There are several predefined mixin templates, but
    // users can add their own as well. The Impl parameter specifies the class with an actual
    // continued fraction implementation which will be instantiated as a data member of the
    // generator class.
    //
    // When update() member function of the generator class is called, the
    // generator class calls the with_next_partial_fraction() member function of the implementation.
    // This function is supposed to take one callback parameter whose function call operator will be
    // called with the newly computed partial fraction is computed. The callback parameter also has
    // a member function get_generator() which returns a const reference to the containing generator
    // object. Since the implementation class does not know about the containing generator class,
    // the type of the callback should be a template parameter.
    //
    // When the callback is called inside with_next_partial_fraction(), it then calls the update()
    // member function of each of the mixins the generator is deriving from, with two arguments, one
    // for the passed new partial fraction, and another for the reference to the instance of the
    // implementation class. After that, the update() function returns true.
    //
    // If the callback is not called and with_next_partial_fraction() returns, then the generator
    // considers that there is no more partial fractions remaining and the continued fraction
    // expansion is done. It then sets the termination flag so that now its terminated() member
    // function returns true. Also, it calls the final_update() member function of each of the
    // mixins with the reference to the implementation object. After that, the update() function
    // returns false. Once termination flag is set, now calling update() does not do anything.
    //
    // Some mixins may properly function only when some other mixins coexist inside the same
    // generator. Such a dependency, if needed, is supposed to be specified by specializing a class
    // template mixin_traits. All the mixins required by the specified mixins are automatically and
    // transitively added to the list of mixins the generator is deriving from, even if not
    // explicltly requested by the user. Some mixins may have ordering constraints; that is, to
    // correctly function they may require certain constraints on the order of the calls to their
    // update() member functions. Such a constraint is supposed to be also specified by specializing
    // the template mixin_traits.
    //
    // Each specialization of mixin_traits may have three member types: required_mixins, before_than
    // and after_than. required_mixins specifies additional mixins required for the proper
    // functioning of the mixin, and the other two specifies that the update() for the mixin should
    // be called before/after the update() of all mixins specified there. Any list of mixins should
    // be given by specializing the struct mixin_list. If a mixin appearing in
    // before_than/after_than is not in the list of mixins the generator is actually deriving
    // from, then it is silently ignored. Any of these three types required_mixins, before_than and
    // after_than can be missing in the specialization of mixin_traits, which has the same effect as
    // specifying an empty mixin_list.
    //
    // Also, some implementations may require some mixins to exist in the same generator. Such a
    // dependency can be also specified if needed by defining a member type alias required_mixins in
    // the implementation class. Any mixins specified there are automatically (and transitively)
    // added to the list of mixins the generator is deriving from, even if not explicitly requested
    // by the user. The implementation class can also specify an additional ordering constraint on
    // mixins by defining a member type alias mixin_ordering_constraints. This alias should be
    // a specialization of the type alias template mixin_ordering_constraint::constraint_list.
    // There are several possible types of parameters for this template alias:
    //   - mixin_ordering_constraint::before_after<Before, After>: the update() of the mixin Before
    //   should be called before the update() of the mixin After.
    //   - mixin_ordering_constraint::linear_chain<Mixin1, Mixin2, ... , MixinN>: the update() of
    //   any mixin appearing first should be called before the update() of any other mixin appearing
    //   after it.
    //
    // Note that the ordering constraint should not form a cycle. In such a case, the compilation
    // will fail with a diagnostic.

    namespace cntfrc {
        template <template <class, class> class... Mixins>
        struct mixin_list {};

        template <template <class, class> class Mixin>
        struct mixin_traits {};

        namespace detail {
            // Template template parameters work not so nicely with metaprogramming, so we wrap
            // mixins into a unique type.
            template <template <class, class> class Mixin>
            struct mixin_type_wrapper {};

            // The alias template type in the above mixin_type_wrapper is a *different* template
            // from Mixin, thus mixin_traits<mixin_type_wrapper<Mixin>::template type> is a
            // *different* type from mixin_traits<Mixin>. To correctly point to the latter from
            // mixin_type_wrapper<Mixin>, we use a helper alias.
            template <template <class, class> class Mixin>
            constexpr mixin_traits<Mixin>
            traits_from_wrapped_mixin_helper(mixin_type_wrapper<Mixin>) noexcept {
                return {};
            }
            template <class WrappedMixin>
            using traits_from_wrapped_mixin =
                decltype(traits_from_wrapped_mixin_helper(WrappedMixin{}));

            // A direct edge from Before to After representing a mixin ordering constraint.
            template <template <class, class> class Before, template <class, class> class After>
            struct mixin_ordering_constraint_edge {};

            // First -> Second -> ... -> Last
            template <class ConstraintList, template <class, class> class First>
            constexpr ConstraintList
            linear_chain_mixin_ordering_constraint_impl(ConstraintList,
                                                        mixin_type_wrapper<First>) noexcept {
                return {};
            }
            template <class ConstraintList, template <class, class> class First,
                      template <class, class> class Second,
                      template <class, class> class... Remaining>
            constexpr auto
            linear_chain_mixin_ordering_constraint_impl(ConstraintList, mixin_type_wrapper<First>,
                                                        mixin_type_wrapper<Second>,
                                                        mixin_type_wrapper<Remaining>...) noexcept {
                return linear_chain_mixin_ordering_constraint_impl(
                    tmp::push_back<ConstraintList, mixin_ordering_constraint_edge<First, Second>>{},
                    mixin_type_wrapper<Second>{}, mixin_type_wrapper<Remaining>{}...);
            }
        }

        namespace mixin_ordering_constraint {
            // Before::update() should be called before After::update().
            template <template <class, class> class Before, template <class, class> class After>
            using before_after =
                tmp::typelist<detail::mixin_ordering_constraint_edge<Before, After>>;

            // First::update() and then Second::update() and then Remaining::update()...
            template <template <class, class> class First, template <class, class> class Second,
                      template <class, class> class... Remaining>
            using linear_chain = decltype(detail::linear_chain_mixin_ordering_constraint_impl(
                tmp::typelist<>{}, detail::mixin_type_wrapper<First>{},
                detail::mixin_type_wrapper<Second>{}, detail::mixin_type_wrapper<Remaining>{}...));

            // E.g. mixin_ordering_constraints =
            //                       constraint_list<before_after<index_tracker, interval_tracker>,
            //                                       linear_chain<partial_fraction_tracker,
            //                                                    convergent_tracaker,
            //                                                    interval_tracker>>;
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
                using type = mixin_list<>;
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
                using type = mixin_list<>;
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
                using type = mixin_list<>;
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
            template <class... WrappedMixins, class... Constraints>
            constexpr auto
            convert_local_mixin_ordering_constraints_list(tmp::typelist<WrappedMixins...>,
                                                          tmp::typelist<Constraints...>) noexcept {
                if constexpr (sizeof...(Constraints) == 0) {
                    return util::array<graph_edge, 0>{};
                }
                else {
                    using list = tmp::typelist<WrappedMixins...>;
                    auto constraint_to_graph_edge =
                        []<template <class, class> class Before,
                           template <class, class> class After>(
                            mixin_ordering_constraint_edge<Before, After>) {
                            return graph_edge{
                                tmp::find_first_index<mixin_type_wrapper<Before>>(list{}),
                                tmp::find_first_index<mixin_type_wrapper<After>>(list{})};
                        };
                    return util::array<graph_edge, sizeof...(Constraints)>{
                        {constraint_to_graph_edge(Constraints{})...}};
                }
            }

            template <class WrappedTargetMixin, class... WrappedMixins,
                      template <class, class> class... Before,
                      template <class, class> class... After>
            constexpr auto convert_single_mixin_traits(tmp::typelist<WrappedMixins...>,
                                                       mixin_list<Before...>,
                                                       mixin_list<After...>) noexcept {
                constexpr std::size_t total_size = sizeof...(Before) + sizeof...(After);

                if constexpr (total_size == 0) {
                    return util::array<graph_edge, 0>{};
                }
                else {
                    using list = tmp::typelist<WrappedMixins...>;
                    constexpr std::size_t idx_of_target =
                        tmp::find_first_index<WrappedTargetMixin>(list{});

                    return util::array<graph_edge, total_size>{
                        {{idx_of_target,
                          tmp::find_first_index<mixin_type_wrapper<Before>>(list{})}...,
                         {tmp::find_first_index<mixin_type_wrapper<After>>(list{}),
                          idx_of_target}...}};
                }
            }

            template <std::size_t... sizes>
            constexpr auto merge_graphs(util::array<graph_edge, sizes> const&... graphs) noexcept {
                std::size_t global_index = 0;
                util::array<graph_edge, (std::size_t(0) + ... + sizes)> ret_value{};
                // I have no idea why the compiler gets mad if I remove this branch...
                if constexpr (ret_value.size() != 0) {
                    auto impl = [&](auto&& graph) {
                        if constexpr (graph.size() != 0) {
                            for (std::size_t idx = 0; idx < graph.size(); ++idx, ++global_index) {
                                ret_value[global_index] = graph[idx];
                            }
                        }
                    };
                    (impl(graphs), ...);
                }
                return ret_value;
            }

            template <class... WrappedMixins>
            constexpr auto
            convert_multiple_mixin_traits(tmp::typelist<WrappedMixins...> list) noexcept {
                return merge_graphs(convert_single_mixin_traits<WrappedMixins>(
                    list,
                    typename get_before_than<traits_from_wrapped_mixin<WrappedMixins>>::type{},
                    typename get_after_than<traits_from_wrapped_mixin<WrappedMixins>>::type{})...);
            }

            // Convert mixin_list<> into tmp::typelist<mixin_type_wrapper<>> for metaprogramming.
            template <template <class, class> class... Mixins>
            constexpr tmp::typelist<mixin_type_wrapper<Mixins>...>
            wrap_mixin_list(mixin_list<Mixins...>) noexcept {
                return {};
            }

            // Convert tmp::typelist<mixin_type_wrapper<>> back into mixin_list<>.
            template <template <class, class> class... Mixins>
            constexpr mixin_list<Mixins...>
            unwrap_mixin_list(tmp::typelist<mixin_type_wrapper<Mixins>...>) noexcept {
                return {};
            }

            // Find the transitive closure of required mixins.
            template <class... WrappedMixins>
            constexpr auto
            get_transitive_required_mixin_list(tmp::typelist<WrappedMixins...>) noexcept {
                using original_list = tmp::typelist<WrappedMixins...>;
                using augmented_list = tmp::remove_duplicate<tmp::join<
                    original_list, decltype(wrap_mixin_list(
                                       typename get_required_mixins<
                                           traits_from_wrapped_mixin<WrappedMixins>>::type{}))...>>;

                if constexpr (std::is_same_v<augmented_list, original_list>) {
                    return augmented_list{};
                }
                else {
                    return get_transitive_required_mixin_list(augmented_list{});
                }
            }

            // Main compile-time function for computing the mixin list.
            template <class LocalConstraintList, template <class, class> class... RequiredMixins,
                      template <class, class> class... AdditionalMixins>
            constexpr auto find_sorted_mixin_list_impl(mixin_list<RequiredMixins...>,
                                                       mixin_list<AdditionalMixins...>) noexcept {
                // Collect all required mixins transitively.
                using wrapped_mixin_list = decltype(get_transitive_required_mixin_list(
                    tmp::remove_duplicate<tmp::typelist<mixin_type_wrapper<RequiredMixins>...,
                                                        mixin_type_wrapper<AdditionalMixins>...>>{}

                    ));

                // Collect all ordering constraints and convert them into an array of graph_edge.
                constexpr auto mixin_dependency_graph =
                    merge_graphs(convert_local_mixin_ordering_constraints_list(
                                     wrapped_mixin_list{}, LocalConstraintList{}),
                                 convert_multiple_mixin_traits(wrapped_mixin_list{}));

                // Get a topological sorted array of mixin indices pointing into wrapped_mixin_list.
                // Any ordering constraint pointing from/to a mixin not included in
                // wrapped_mixin_list is ignored.
                constexpr auto sorted_mixin_indices =
                    topological_sort<wrapped_mixin_list::size>(mixin_dependency_graph);

                static_assert(sorted_mixin_indices.succeed,
                              "mixin's ordering constraint should not form a cylce");

                // Convert the index array into a tmp::typelist and return.
                return [&sorted_mixin_indices]<std::size_t... I>(std::index_sequence<I...>) {
                    return tmp::typelist<tmp::get_type<sorted_mixin_indices.sorted_indices[I],
                                                       wrapped_mixin_list>...>{};
                }(std::make_index_sequence<wrapped_mixin_list::size>{});
            }

            template <class Impl, template <class, class> class... AdditionalMixins>
            constexpr auto find_sorted_mixin_list() noexcept {
                return find_sorted_mixin_list_impl<
                    typename get_mixin_ordering_constraints<Impl>::type>(
                    typename get_required_mixins<Impl>::type{}, mixin_list<AdditionalMixins...>{});
            }

            template <class Impl, template <class, class> class... Mixins>
            class generator_impl : public Mixins<Impl, generator_impl<Impl, Mixins...>>... {
            public:
                using impl_type = Impl;
                using partial_fraction_type = typename Impl::partial_fraction_type;
                using convergent_type = typename Impl::convergent_type;
                using interval_type = typename Impl::interval_type;

            private:
                Impl impl_;
                bool terminated_ = false;

                struct callback_type {
                private:
                    generator_impl& gen_;

                    friend generator_impl;
                    explicit constexpr callback_type(generator_impl& gen) noexcept : gen_{gen} {}

                public:
                    constexpr generator_impl const& get_generator() const noexcept { return gen_; }

                    constexpr void operator()(partial_fraction_type const& next_partial_fraction) {
                        (static_cast<Mixins<Impl, generator_impl>&>(gen_).update(
                             next_partial_fraction, gen_.impl_),
                         ...);
                        gen_.terminated_ = false;
                    }
                };

            public:
                explicit constexpr generator_impl(Impl impl)
                    : Mixins<Impl, generator_impl>{impl}..., impl_{static_cast<Impl&&>(impl)} {}

                // Returns true if succeeded obtaining a further partial fraction.
                constexpr bool update() {
                    if (!terminated_) {
                        terminated_ = true;
                        impl_.with_next_partial_fraction(callback_type{*this});

                        if (terminated_) {
                            auto invoke_final_update = [this](auto&& mixin) {
                                if constexpr (requires { mixin.final_update(impl_); }) {
                                    mixin.final_update(impl_);
                                }
                            };
                            (invoke_final_update(static_cast<Mixins<Impl, generator_impl>&>(*this)),
                             ...);
                        }
                    }
                    return !terminated_;
                }

                constexpr bool terminated() const noexcept { return terminated_; }
            };

            template <class Impl, class MixinList>
            struct get_generator_type;

            template <class Impl, template <class, class> class... Mixins>
            struct get_generator_type<Impl, mixin_list<Mixins...>> {
                using type = generator_impl<Impl, Mixins...>;
            };
        }

        // The main interface type.
        // Automatically add all transitively required mixins into the list of mixins, and sort the
        // list according to the imposed ordering constraints. Then generator derives from the
        // resulting list of mixins.
        template <class Impl, template <class, class> class... Mixins>
        using generator = typename detail::get_generator_type<
            Impl, decltype(detail::unwrap_mixin_list(
                      detail::find_sorted_mixin_list<Impl, Mixins...>()))>::type;

        // A convenient factory function.
        template <template <class, class> class... Mixins, class Impl>
        constexpr auto make_generator(Impl&& impl) {
            return generator<std::remove_cvref_t<Impl>, Mixins...>{static_cast<Impl&&>(impl)};
        }

        // Mixin: stores the current index of the continued fraction expansion.
        template <class Impl, class Generator>
        class index_tracker {
            using partial_fraction_type = typename Impl::partial_fraction_type;
            int current_index_ = -1;

            friend Generator;

            explicit constexpr index_tracker(Impl const&) noexcept {}

            constexpr void update(partial_fraction_type const&, Impl const&) noexcept {
                ++current_index_;
            }

        public:
            constexpr int current_index() const noexcept { return current_index_; }
        };

        // Mixin: stores the current partial fraction of the continued fraction expansion.
        template <class Impl, class Generator>
        class partial_fraction_tracker {
            using partial_fraction_type = typename Impl::partial_fraction_type;

            partial_fraction_type current_partial_fraction_;

            friend Generator;

            explicit constexpr partial_fraction_tracker(Impl const& impl)
                : current_partial_fraction_{impl.initial_partial_fraction()} {}

            constexpr void update(partial_fraction_type const& next_partial_fraction, Impl const&) {
                current_partial_fraction_ = next_partial_fraction;
            }

        public:
            constexpr partial_fraction_type const& current_partial_fraction() const noexcept {
                return current_partial_fraction_;
            }
        };

        // Mixin: stores the previous and the current convergents of the continued fraction
        // expansion.
        template <class Impl, class Generator>
        class convergent_tracker {
            using partial_fraction_type = typename Impl::partial_fraction_type;
            using convergent_type = typename Impl::convergent_type;

            convergent_type current_convergent_{1, 0u};
            convergent_type previous_convergent_{0, 1u};

            friend Generator;

            explicit constexpr convergent_tracker(Impl const&) noexcept {}

            constexpr void update(partial_fraction_type const& next_partial_fraction, Impl const&) {
                auto next_numerator =
                    next_partial_fraction.denominator * current_convergent_numerator() +
                    next_partial_fraction.numerator * previous_convergent_numerator();
                auto next_denominator =
                    next_partial_fraction.denominator * current_convergent_denominator() +
                    next_partial_fraction.numerator * previous_convergent_denominator();

                previous_convergent_ = static_cast<convergent_type&&>(current_convergent_);
                current_convergent_ =
                    convergent_type{std::move(next_numerator), abs(std::move(next_denominator))};
            }

        public:
            constexpr convergent_type const& current_convergent() const noexcept {
                return current_convergent_;
            }
            constexpr auto const& current_convergent_numerator() const noexcept {
                return current_convergent().numerator;
            }
            constexpr auto const& current_convergent_denominator() const noexcept {
                return current_convergent().denominator;
            }

            constexpr convergent_type const& previous_convergent() const noexcept {
                return previous_convergent_;
            }
            constexpr auto const& previous_convergent_numerator() const noexcept {
                return previous_convergent().numerator;
            }
            constexpr auto const& previous_convergent_denominator() const noexcept {
                return previous_convergent().denominator;
            }
        };

        // Mixin: stores the last three convergents (including the current one) of the continued
        // fraction expansion.
        template <class Impl, class Generator>
        class convergent_triple_tracker {
            using partial_fraction_type = typename Impl::partial_fraction_type;
            using convergent_type = typename Impl::convergent_type;

            convergent_type previous_previous_convergent_{1, 0u};

            friend Generator;

            explicit constexpr convergent_triple_tracker(Impl const&) noexcept {}

            constexpr void update(partial_fraction_type const&, Impl const&) {
                previous_previous_convergent_ =
                    static_cast<Generator const&>(*this).previous_convergent();
            }

        public:
            constexpr convergent_type const& previous_previous_convergent() const noexcept {
                return previous_previous_convergent_;
            }
            constexpr auto const& previous_previous_convergent_numerator() const noexcept {
                return previous_previous_convergent().numerator;
            }
            constexpr auto const& previous_previous_convergent_denominator() const noexcept {
                return previous_previous_convergent().denominator;
            }
        };
        template <>
        struct mixin_traits<convergent_triple_tracker> {
            using required_mixins = mixin_list<convergent_tracker>;
            using before_than = mixin_list<convergent_tracker>;
        };

        // Mixin: stores the current estimate of an interval where the true value should be in.
        template <class Impl, class Generator>
        class interval_tracker {
            using partial_fraction_type = typename Impl::partial_fraction_type;
            using convergent_type = typename Impl::convergent_type;
            using interval_type = typename Impl::interval_type;

            interval_type current_interval_;

            friend Generator;

            explicit constexpr interval_tracker(Impl const& impl)
                : current_interval_{impl.initial_interval()} {}

            constexpr void update(partial_fraction_type const&, Impl& impl) {
                auto const& self = static_cast<Generator const&>(*this);
                if constexpr (requires { impl.next_interval(self); }) {
                    current_interval_ = impl.next_interval(self);
                }
                else {
                    // Use convergents. Note that this is valid only for regular continued
                    // fractions.
                    if (self.current_index() >= 1) {
                        if (self.current_index() % 2 == 0) {
                            current_interval_ =
                                cyclic_interval<typename interval_type::value_type,
                                                cyclic_interval_type_t::left_closed_right_open>{
                                    self.current_convergent(), self.previous_convergent()};
                        }
                        else {
                            current_interval_ =
                                cyclic_interval<typename interval_type::value_type,
                                                cyclic_interval_type_t::left_open_right_closed>{
                                    self.previous_convergent(), self.current_convergent()};
                        }
                    }
                }
            }
            constexpr void final_update(Impl const&) {
                current_interval_ = cyclic_interval<typename interval_type::value_type,
                                                    cyclic_interval_type_t::single_point>{
                    static_cast<Generator const&>(*this).current_convergent()};
            }

        public:
            interval_type const& current_interval() const noexcept { return current_interval_; }

            // Refine the current value so that the maximum possible error is no more than the
            // given bound.
            template <class ErrorValue>
            convergent_type const& progress_until(ErrorValue const& error_bound) {
                // Translate the lower bound by error_bound and see if it belongs to the interval.
                while (current_interval().visit([&](auto&& itv) {
                    static_assert(itv.interval_type() != cyclic_interval_type_t::empty);
                    if constexpr (itv.interval_type() == cyclic_interval_type_t::single_point) {
                        return false;
                    }
                    else if constexpr (itv.interval_type() == cyclic_interval_type_t::entire) {
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
                    static_cast<Generator&>(*this).update();
                }
                return static_cast<Generator const&>(*this).current_convergent();
            }
        };
        template <>
        struct mixin_traits<interval_tracker> {
            using required_mixins = mixin_list<index_tracker, convergent_tracker>;
            using after_than = mixin_list<index_tracker, convergent_tracker>;
        };
    }
}

#endif
