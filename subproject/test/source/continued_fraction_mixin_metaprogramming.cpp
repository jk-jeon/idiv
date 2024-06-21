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

#include <idiv/continued_fraction.h>
#include <boost/ut.hpp>

void continued_fraction_mixin_metaprogramming_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Metaprogramming stuffs for continued fraction mixins]"_test = [] {
        should("topological_sort") = [] {
            {
                //       --- 2 -
                //      /       \
                //      |  - 3 - 6 -
                //      | /         \
                // 0 -- 1 ------ 7 - 9
                //      |\        \ /
                //      |  - 4 - 8 -
                //      \       /   \
                //       --- 5 -     10
                util::array<cntfrc::detail::graph_edge, 14> edges{{{0, 1},
                                                                   {1, 2},
                                                                   {1, 3},
                                                                   {1, 7},
                                                                   {1, 4},
                                                                   {1, 5},
                                                                   {2, 6},
                                                                   {3, 6},
                                                                   {4, 8},
                                                                   {5, 8},
                                                                   {6, 9},
                                                                   {7, 9},
                                                                   {8, 9},
                                                                   {7, 10}}};
                auto result = cntfrc::detail::topological_sort<11>(edges);
                expect(result.succeed == true);
                expect(result.sorted_indices[0] == 0);
                expect(result.sorted_indices[1] == 1);
                expect(result.sorted_indices[2] == 5);
                expect(result.sorted_indices[3] == 4);
                expect(result.sorted_indices[4] == 8);
                expect(result.sorted_indices[5] == 7);
                expect(result.sorted_indices[6] == 10);
                expect(result.sorted_indices[7] == 3);
                expect(result.sorted_indices[8] == 2);
                expect(result.sorted_indices[9] == 6);
                expect(result.sorted_indices[10] == 9);
            }
            {
                // Cyclic case.
                util::array<cntfrc::detail::graph_edge, 14> edges{{{0, 1}, {1, 2}, {2, 3}, {3, 0}}};
                auto result = cntfrc::detail::topological_sort<4>(edges);
                expect(result.succeed == false);
            }
        };

        should("get_transitive_required_mixin_list") = [] {
            auto list = cntfrc::detail::get_transitive_required_mixin_list(
                tmp::typelist<cntfrc::detail::mixin_type_wrapper<cntfrc::interval_tracker>,
                              cntfrc::detail::mixin_type_wrapper<cntfrc::convergent_tracker>>{});

            expect(std::is_same_v<
                   decltype(list),
                   tmp::typelist<cntfrc::detail::mixin_type_wrapper<cntfrc::interval_tracker>,
                                 cntfrc::detail::mixin_type_wrapper<cntfrc::convergent_tracker>,
                                 cntfrc::detail::mixin_type_wrapper<cntfrc::index_tracker>>>);
        };

        should("find_sorted_mixin_list") = [] {
            struct dummy_type {
                using required_mixins =
                    cntfrc::mixin_list<cntfrc::index_tracker, cntfrc::partial_fraction_tracker,
                                       cntfrc::convergent_tracker>;
                using mixin_ordering_constraints =
                    cntfrc::mixin_ordering_constraint::constraint_list<
                        cntfrc::mixin_ordering_constraint::before_after<cntfrc::index_tracker,
                                                                        cntfrc::interval_tracker>,
                        cntfrc::mixin_ordering_constraint::before_after<cntfrc::convergent_tracker,
                                                                        cntfrc::interval_tracker>>;
            };
            auto sorted_wrapped_mixin_list =
                cntfrc::detail::find_sorted_mixin_list<dummy_type, cntfrc::interval_tracker>();
            expect(std::is_same_v<
                   decltype(sorted_wrapped_mixin_list),
                   tmp::typelist<
                       cntfrc::detail::mixin_type_wrapper<cntfrc::convergent_tracker>,
                       cntfrc::detail::mixin_type_wrapper<cntfrc::partial_fraction_tracker>,
                       cntfrc::detail::mixin_type_wrapper<cntfrc::index_tracker>,
                       cntfrc::detail::mixin_type_wrapper<cntfrc::interval_tracker>>>);
        };
    };
}