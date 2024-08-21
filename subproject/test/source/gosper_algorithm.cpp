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

#include <idiv/bigint.h>
#include <idiv/continued_fraction/engine/gosper.h>
#include <idiv/continued_fraction/engine/rational.h>
#include <boost/ut.hpp>

#include <vector>
#include <variant>

void gosper_algorithm_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Gosper's algorithms]"_test = [] {
        should("extended_linear_fractional_mapping") = [] {
            using enum cyclic_interval_type_t;
            using value_type = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
            using mapping_type = cntfrc::extended_linear_fractional_mapping<bigint::int_var>;

            // Return type of map_cyclic_interval.
            {
                expect(std::is_same_v<decltype(std::declval<mapping_type>().map_cyclic_interval(
                                          std::declval<cyclic_interval<value_type, empty>>())),
                                      cyclic_interval<value_type, empty>>);
                expect(
                    std::is_same_v<decltype(std::declval<mapping_type>().map_cyclic_interval(
                                       std::declval<cyclic_interval<value_type, single_point>>())),
                                   cyclic_interval<value_type, single_point>>);
                expect(
                    std::is_same_v<decltype(std::declval<mapping_type>().map_cyclic_interval(
                                       std::declval<cyclic_interval<value_type, open>>())),
                                   variable_shape_cyclic_interval<value_type, single_point, open>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_interval(
                           std::declval<cyclic_interval<value_type, left_open_right_closed>>())),
                       variable_shape_cyclic_interval<value_type, single_point,
                                                      left_open_right_closed,
                                                      left_closed_right_open>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_interval(
                           std::declval<cyclic_interval<value_type, left_closed_right_open>>())),
                       variable_shape_cyclic_interval<value_type, single_point,
                                                      left_open_right_closed,
                                                      left_closed_right_open>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_interval(
                           std::declval<cyclic_interval<value_type, closed>>())),
                       variable_shape_cyclic_interval<value_type, single_point, closed>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_interval(
                           std::declval<cyclic_interval<value_type, single_complement>>())),
                       variable_shape_cyclic_interval<value_type, single_point,
                                                      single_complement>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_interval(
                           std::declval<cyclic_interval<value_type, entire>>())),
                       variable_shape_cyclic_interval<value_type, single_point, entire>>);

                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_interval(
                           std::declval<variable_shape_cyclic_interval<
                               value_type, single_point, closed, single_complement>>())),
                       variable_shape_cyclic_interval<value_type, single_point, closed,
                                                      single_complement>>);
                expect(std::is_same_v<decltype(std::declval<mapping_type>().map_cyclic_interval(
                                          std::declval<variable_shape_cyclic_interval<
                                              value_type, empty, open, entire>>())),
                                      variable_shape_cyclic_interval<value_type, empty,
                                                                     single_point, open, entire>>);
                expect(
                    std::is_same_v<decltype(std::declval<mapping_type>().map_cyclic_interval(
                                       std::declval<variable_shape_cyclic_interval<
                                           value_type, left_open_right_closed, closed>>())),
                                   variable_shape_cyclic_interval<value_type, single_point,
                                                                  left_open_right_closed,
                                                                  left_closed_right_open, closed>>);
            }

            // Positive determinant case.
            {
                cntfrc::extended_linear_fractional_mapping<bigint::int_var> f{1, 2, -1, 3};
                expect(f.determinant_sign() == 1);
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{7, 4}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{3, 1});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{7, 4u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{3, 1});


                // Return value of map_cyclic_interval.
                {
                    auto const zero = value_type{0, 1u};
                    auto const infinity = value_type{1, 0u};

                    expect(f.map_cyclic_interval(cyclic_interval<value_type, empty>{}) ==
                           cyclic_interval<value_type, empty>{});
                    expect(f.map_cyclic_interval(cyclic_interval<value_type, single_point>{zero}) ==
                           cyclic_interval<value_type, single_point>{value_type{2, 3u}});
                    expect(
                        f.map_cyclic_interval(cyclic_interval<value_type, open>{zero, infinity}) ==
                        cyclic_interval<value_type, open>{value_type{2, 3u}, value_type{-1, 1u}});
                    expect(
                        f.map_cyclic_interval(cyclic_interval<value_type, left_open_right_closed>{
                            zero, infinity}) == cyclic_interval<value_type, left_open_right_closed>{
                                                    value_type{2, 3u}, value_type{-1, 1u}});
                    expect(
                        f.map_cyclic_interval(cyclic_interval<value_type, left_closed_right_open>{
                            zero, infinity}) == cyclic_interval<value_type, left_closed_right_open>{
                                                    value_type{2, 3u}, value_type{-1, 1u}});
                    expect(f.map_cyclic_interval(cyclic_interval<value_type, closed>{
                               zero, infinity}) == cyclic_interval<value_type, closed>{
                                                       value_type{2, 3u}, value_type{-1, 1u}});

                    expect(f.map_cyclic_interval(
                               cyclic_interval<value_type, single_complement>{zero}) ==
                           cyclic_interval<value_type, single_complement>{value_type{2, 3u}});
                    expect(f.map_cyclic_interval(cyclic_interval<value_type, entire>{}) ==
                           cyclic_interval<value_type, entire>{});
                }
            }

            // Negative determinant case.
            {
                cntfrc::extended_linear_fractional_mapping<bigint::int_var> f{1, 2, 3, 0};
                expect(f.determinant_sign() == -1);
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{7, 4}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{5, 7});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{7, 4u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{5, 7});

                // Return value of map_cyclic_interval.
                {
                    auto const zero = value_type{0, 1u};
                    auto const infinity = value_type{1, 0u};

                    expect(f.map_cyclic_interval(cyclic_interval<value_type, empty>{}) ==
                           cyclic_interval<value_type, empty>{});
                    expect(f.map_cyclic_interval(cyclic_interval<value_type, single_point>{zero}) ==
                           cyclic_interval<value_type, single_point>{infinity});
                    expect(
                        f.map_cyclic_interval(cyclic_interval<value_type, open>{zero, infinity}) ==
                        cyclic_interval<value_type, open>{value_type{1, 3u}, infinity});
                    expect(
                        f.map_cyclic_interval(cyclic_interval<value_type, left_open_right_closed>{
                            zero, infinity}) == cyclic_interval<value_type, left_closed_right_open>{
                                                    value_type{1, 3u}, infinity});
                    expect(
                        f.map_cyclic_interval(cyclic_interval<value_type, left_closed_right_open>{
                            zero, infinity}) == cyclic_interval<value_type, left_open_right_closed>{
                                                    value_type{1, 3u}, infinity});
                    expect(f.map_cyclic_interval(
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, closed>{value_type{1, 3u}, infinity});

                    expect(f.map_cyclic_interval(cyclic_interval<value_type, single_complement>{
                               zero}) == cyclic_interval<value_type, single_complement>{infinity});
                    expect(f.map_cyclic_interval(cyclic_interval<value_type, entire>{}) ==
                           cyclic_interval<value_type, entire>{});
                }
            }

            // Singular case.
            {
                cntfrc::extended_linear_fractional_mapping<bigint::int_var> f{1, 2, 2, 4};
                expect(f.determinant_sign() == 0);
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-2, 1}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{1, 2});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{-2, 1u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{1, 2});

                // Return value of map_cyclic_interval.
                {
                    auto const singularity = value_type{-2, 1u};
                    auto const zero = value_type{0, 1u};
                    auto const constant_value = value_type{1, 2u};

                    expect(f.map_cyclic_interval(cyclic_interval<value_type, empty>{}) ==
                           cyclic_interval<value_type, empty>{});
                    expect(f.map_cyclic_interval(
                               cyclic_interval<value_type, single_point>{singularity}) ==
                           cyclic_interval<value_type, single_point>{constant_value});
                    expect(f.map_cyclic_interval(
                               cyclic_interval<value_type, open>{singularity, zero}) ==
                           cyclic_interval<value_type, single_point>{constant_value});
                    expect(
                        f.map_cyclic_interval(cyclic_interval<value_type, left_open_right_closed>{
                            singularity, zero}) ==
                        cyclic_interval<value_type, single_point>{constant_value});
                    expect(
                        f.map_cyclic_interval(cyclic_interval<value_type, left_closed_right_open>{
                            singularity, zero}) ==
                        cyclic_interval<value_type, single_point>{constant_value});
                    expect(f.map_cyclic_interval(
                               cyclic_interval<value_type, closed>{singularity, zero}) ==
                           cyclic_interval<value_type, single_point>{constant_value});

                    expect(f.map_cyclic_interval(
                               cyclic_interval<value_type, single_complement>{singularity}) ==
                           cyclic_interval<value_type, single_point>{constant_value});
                    expect(f.map_cyclic_interval(cyclic_interval<value_type, entire>{}) ==
                           cyclic_interval<value_type, single_point>{constant_value});
                }
            }
        };

        using rational_continued_fraction_t =
            cntfrc::engine::rational<bigint::int_var, bigint::uint_var>;

        should("unary_gosper") = [] {
            // 481/2245 = (-18*156 + 13*179)/(12*156 - 23*179)
            auto cf1 = cntfrc::make_unary_gosper_generator<cntfrc::convergent_tracker>(
                rational_continued_fraction_t{156, 179u}, {-18, 13, 12, -23});
            auto cf2 = cntfrc::make_generator<cntfrc::convergent_tracker>(
                rational_continued_fraction_t{481, 2245u});

            auto itr1 = cf1.begin();
            auto itr2 = cf2.begin();
            for (; itr2 != cf2.end(); ++itr1, ++itr2) {
                expect(itr1 != cf1.end());
                expect(itr1->current_convergent() == itr2->current_convergent());
            }
            expect(itr1 == cf1.end());
        };

        should("binary_gosper") = [] {
            // Take x = 17/89, y = 31/125, and
            // z = (xy + 4x + 2y + 8)/(2x - 3y + 1) = 2655/182.
            auto cf1 = cntfrc::make_binary_gosper_generator<cntfrc::convergent_tracker>(
                rational_continued_fraction_t{17, 89u}, rational_continued_fraction_t{31, 125u},
                {// numerator
                 1, 4, 2, 8,
                 // denominator
                 0, 2, -3, 1});
            auto cf2 = cntfrc::make_generator<cntfrc::convergent_tracker>(
                rational_continued_fraction_t{2655, 182u});

            auto itr1 = cf1.begin();
            auto itr2 = cf2.begin();
            for (; itr2 != cf2.end(); ++itr1, ++itr2) {
                expect(itr1 != cf1.end());
                expect(itr1->current_convergent() == itr2->current_convergent());
            }
            expect(itr1 == cf1.end());
        };
    };
}
