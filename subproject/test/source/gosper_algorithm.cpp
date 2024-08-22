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
                mapping_type f{1, 2, -1, 3};
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
                mapping_type f{1, 2, 3, 0};
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
                mapping_type f{1, 2, 2, 4};
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

        should("extended_bilinear_fractional_mapping") = [] {
            using enum cyclic_interval_type_t;
            using value_type = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
            using mapping_type = cntfrc::extended_bilinear_fractional_mapping<bigint::int_var>;

            // Return type of map_cyclic_rectangle<false>.
            {
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                           std::declval<
                               variable_shape_cyclic_interval<value_type, empty, single_point>>(),
                           std::declval<
                               variable_shape_cyclic_interval<value_type, single_point>>())),
                       variable_shape_cyclic_interval<value_type, empty, single_point>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                           std::declval<variable_shape_cyclic_interval<value_type, open, entire>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, single_point,
                                                                       single_complement>>())),
                       variable_shape_cyclic_interval<value_type, single_point, open,
                                                      single_complement, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                           std::declval<
                               variable_shape_cyclic_interval<value_type, single_point, entire>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, single_point,
                                                                       single_complement>>())),
                       variable_shape_cyclic_interval<value_type, empty, single_point,
                                                      single_complement, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                           std::declval<
                               variable_shape_cyclic_interval<value_type, single_point, closed>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, single_point,
                                                                       closed, entire>>())),
                       variable_shape_cyclic_interval<value_type, empty, single_point, closed,
                                                      entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                           std::declval<variable_shape_cyclic_interval<value_type, open, closed>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, entire>>())),
                       variable_shape_cyclic_interval<value_type, single_point, open, closed,
                                                      entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                           std::declval<variable_shape_cyclic_interval<value_type, open>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, open>>())),
                       variable_shape_cyclic_interval<value_type, single_point, open,
                                                      single_complement, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<variable_shape_cyclic_interval<value_type, open>>(),
                           std::declval<variable_shape_cyclic_interval<value_type,
                                                                       left_open_right_closed>>())),
                       variable_shape_cyclic_interval<
                           value_type, single_point, open, left_open_right_closed,
                           left_closed_right_open, single_complement, entire>>);
                expect(
                    std::is_same_v<decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                                       std::declval<variable_shape_cyclic_interval<
                                           value_type, left_open_right_closed>>(),
                                       std::declval<variable_shape_cyclic_interval<
                                           value_type, left_open_right_closed>>())),
                                   variable_shape_cyclic_interval<
                                       value_type, single_point, open, left_open_right_closed,
                                       left_closed_right_open, closed, single_complement, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                           std::declval<variable_shape_cyclic_interval<
                               value_type, open, left_open_right_closed, closed, entire>>(),
                           std::declval<
                               variable_shape_cyclic_interval<value_type, single_point>>())),
                       variable_shape_cyclic_interval<value_type, single_point, open,
                                                      left_open_right_closed,
                                                      left_closed_right_open, closed, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle(
                           std::declval<variable_shape_cyclic_interval<
                               value_type, open, left_open_right_closed, closed, entire>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, open>>())),
                       variable_shape_cyclic_interval<
                           value_type, single_point, open, left_open_right_closed,
                           left_closed_right_open, closed, single_complement, entire>>);
            }
            // Return type of map_cyclic_rectangle<true>.
            {
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<
                               variable_shape_cyclic_interval<value_type, empty, single_point>>(),
                           std::declval<
                               variable_shape_cyclic_interval<value_type, single_point>>())),
                       variable_shape_cyclic_interval<value_type, empty, single_point>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<variable_shape_cyclic_interval<value_type, open, entire>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, single_point,
                                                                       single_complement>>())),
                       variable_shape_cyclic_interval<value_type, single_point, open,
                                                      single_complement, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<
                               variable_shape_cyclic_interval<value_type, single_point, entire>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, single_point,
                                                                       single_complement>>())),
                       variable_shape_cyclic_interval<value_type, single_point, single_complement,
                                                      entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<
                               variable_shape_cyclic_interval<value_type, single_point, closed>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, single_point,
                                                                       closed, entire>>())),
                       variable_shape_cyclic_interval<value_type, single_point, closed, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<variable_shape_cyclic_interval<value_type, open, closed>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, entire>>())),
                       variable_shape_cyclic_interval<value_type, single_point, open, closed,
                                                      entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<variable_shape_cyclic_interval<value_type, open>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, open>>())),
                       variable_shape_cyclic_interval<value_type, single_point, open,
                                                      single_complement, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<variable_shape_cyclic_interval<value_type, open>>(),
                           std::declval<variable_shape_cyclic_interval<value_type,
                                                                       left_open_right_closed>>())),
                       variable_shape_cyclic_interval<
                           value_type, single_point, open, left_open_right_closed,
                           left_closed_right_open, single_complement, entire>>);
                expect(
                    std::is_same_v<decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                                       std::declval<variable_shape_cyclic_interval<
                                           value_type, left_open_right_closed>>(),
                                       std::declval<variable_shape_cyclic_interval<
                                           value_type, left_open_right_closed>>())),
                                   variable_shape_cyclic_interval<
                                       value_type, single_point, open, left_open_right_closed,
                                       left_closed_right_open, single_complement, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<variable_shape_cyclic_interval<
                               value_type, open, left_open_right_closed, closed, entire>>(),
                           std::declval<
                               variable_shape_cyclic_interval<value_type, single_point>>())),
                       variable_shape_cyclic_interval<value_type, single_point, open,
                                                      left_open_right_closed,
                                                      left_closed_right_open, closed, entire>>);
                expect(std::is_same_v<
                       decltype(std::declval<mapping_type>().map_cyclic_rectangle<true>(
                           std::declval<variable_shape_cyclic_interval<
                               value_type, open, left_open_right_closed, closed, entire>>(),
                           std::declval<variable_shape_cyclic_interval<value_type, open>>())),
                       variable_shape_cyclic_interval<
                           value_type, single_point, open, left_open_right_closed,
                           left_closed_right_open, closed, single_complement, entire>>);
            }

            // Constant case.
            {
                mapping_type f{// numerator.
                               1, 2, 0, 3,
                               // denominator.
                               2, 4, 0, 6};
                expect(f.kind() == mapping_type::kind_t::constant);
                expect(f.number_of_points_in_indeterminacy_locus() == 0);

                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-15, 11},
                         cntfrc::projective_rational<bigint::int_var, bigint::int_var>{1, 5}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{1, 2});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{-15, 11u},
                         cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{1, 5u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{1, 2});

                auto const singularity_x = value_type{-15, 11u};
                auto const singularity_y = value_type{1, 5u};
                auto const zero = value_type{0, 1u};
                auto const constant_value = value_type{1, 2u};

                // Return value of map_cyclic_rectangle.
                {
                    auto test = [&f](auto const& itv_x, auto const& itv_y,
                                     auto const& expected_answer) {
                        auto answer1 = f.map_cyclic_rectangle<false>(itv_x, itv_y);
                        auto answer2 = f.map_cyclic_rectangle<true>(itv_x, itv_y);
                        return answer1 == expected_answer && answer1 == answer2;
                    };

                    expect(test(cyclic_interval<value_type, single_point>{singularity_x},
                                cyclic_interval<value_type, single_point>{singularity_y},
                                cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(cyclic_interval<value_type, open>{zero, singularity_x},
                                cyclic_interval<value_type, open>{zero, singularity_y},
                                cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(cyclic_interval<value_type, open>{zero, singularity_x},
                                cyclic_interval<value_type, closed>{zero, singularity_y},
                                cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(cyclic_interval<value_type, closed>{zero, singularity_x},
                                cyclic_interval<value_type, closed>{zero, singularity_y},
                                cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(
                        cyclic_interval<value_type, left_open_right_closed>{zero, singularity_x},
                        cyclic_interval<value_type, left_closed_right_open>{zero, singularity_y},
                        cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(
                        cyclic_interval<value_type, left_open_right_closed>{zero, singularity_x},
                        cyclic_interval<value_type, open>{zero, singularity_y},
                        cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(
                        cyclic_interval<value_type, closed>{zero, singularity_x},
                        cyclic_interval<value_type, left_closed_right_open>{zero, singularity_y},
                        cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, open>{zero, singularity_y},
                                cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(cyclic_interval<value_type, open>{zero, singularity_x},
                                cyclic_interval<value_type, single_complement>{zero},
                                cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(cyclic_interval<value_type, closed>{zero, singularity_x},
                                cyclic_interval<value_type, single_complement>{zero},
                                cyclic_interval<value_type, single_point>{constant_value}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, single_point>{constant_value}));
                }
            }

            // Function-of-only-x case.
            {
                // ((3x - 1)(2y + 7))/((x + 4)(2y + 7)) = (6xy + 21x - 2y - 7)/(2xy + 7x + 8y + 28).
                mapping_type f{// numerator.
                               6, 21, -2, -7,
                               // denominator.
                               2, 7, 8, 28};
                expect(f.kind() == mapping_type::kind_t::nonconstant_unary_function_of_x);
                expect(f.number_of_points_in_indeterminacy_locus() == 0);

                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{1, 3},
                         cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-7, 2}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{0, 1});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{1, 3u},
                         cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{-7, 2u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{0, 1});

                auto const singularity_y = value_type{-7, 2u};
                auto const zero = value_type{0, 1u};
                auto const infinity = value_type{1, 0u};

                // Return value of map_cyclic_rectangle.
                {
                    auto test = [&f](auto const& itv_x, auto const& itv_y,
                                     auto const& expected_answer) {
                        auto answer1 = f.map_cyclic_rectangle<false>(itv_x, itv_y);
                        auto answer2 = f.map_cyclic_rectangle<true>(itv_x, itv_y);
                        return answer1 == expected_answer && answer1 == answer2;
                    };

                    expect(test(cyclic_interval<value_type, single_point>{zero},
                                cyclic_interval<value_type, single_point>{singularity_y},
                                cyclic_interval<value_type, single_point>{value_type{-1, 4u}}));
                    expect(test(
                        cyclic_interval<value_type, open>{zero, infinity},
                        cyclic_interval<value_type, open>{zero, singularity_y},
                        cyclic_interval<value_type, open>{value_type{-1, 4u}, value_type{3, 1u}}));
                    expect(test(
                        cyclic_interval<value_type, open>{zero, infinity},
                        cyclic_interval<value_type, closed>{zero, singularity_y},
                        cyclic_interval<value_type, open>{value_type{-1, 4u}, value_type{3, 1u}}));
                    expect(test(cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, closed>{zero, singularity_y},
                                cyclic_interval<value_type, closed>{value_type{-1, 4u},
                                                                    value_type{3, 1u}}));
                    expect(test(
                        cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                        cyclic_interval<value_type, left_closed_right_open>{zero, singularity_y},
                        cyclic_interval<value_type, left_open_right_closed>{value_type{-1, 4u},
                                                                            value_type{3, 1u}}));
                    expect(test(cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                                cyclic_interval<value_type, open>{zero, singularity_y},
                                cyclic_interval<value_type, left_open_right_closed>{
                                    value_type{-1, 4u}, value_type{3, 1u}}));
                    expect(test(
                        cyclic_interval<value_type, closed>{zero, infinity},
                        cyclic_interval<value_type, left_closed_right_open>{zero, singularity_y},
                        cyclic_interval<value_type, closed>{value_type{-1, 4u},
                                                            value_type{3, 1u}}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, open>{zero, singularity_y},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(
                        cyclic_interval<value_type, open>{zero, infinity},
                        cyclic_interval<value_type, single_complement>{zero},
                        cyclic_interval<value_type, open>{value_type{-1, 4u}, value_type{3, 1u}}));
                    expect(test(cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, single_complement>{zero},
                                cyclic_interval<value_type, closed>{value_type{-1, 4u},
                                                                    value_type{3, 1u}}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{}));
                }
            }

            // Function-of-only-y case.
            {
                // ((3x - 1)(2y + 7))/((3x - 1)(y - 3)) = (6xy + 21x - 2y - 7)/(3xy - 9x - y + 3).
                mapping_type f{// numerator.
                               6, 21, -2, -7,
                               // denominator.
                               3, -9, -1, 3};
                expect(f.kind() == mapping_type::kind_t::nonconstant_unary_function_of_y);
                expect(f.number_of_points_in_indeterminacy_locus() == 0);

                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{1, 3},
                         cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-7, 2}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{0, 1});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{1, 3u},
                         cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{-7, 2u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{0, 1});

                auto const singularity_x = value_type{1, 3u};
                auto const zero = value_type{0, 1u};
                auto const infinity = value_type{1, 0u};

                // Return value of map_cyclic_rectangle.
                {
                    auto test = [&f](auto const& itv_x, auto const& itv_y,
                                     auto const& expected_answer) {
                        auto answer1 = f.map_cyclic_rectangle<false>(itv_x, itv_y);
                        auto answer2 = f.map_cyclic_rectangle<true>(itv_x, itv_y);
                        return answer1 == expected_answer && answer1 == answer2;
                    };

                    expect(test(cyclic_interval<value_type, single_point>{singularity_x},
                                cyclic_interval<value_type, single_point>{zero},
                                cyclic_interval<value_type, single_point>{value_type{-7, 3u}}));
                    expect(test(
                        cyclic_interval<value_type, open>{zero, singularity_x},
                        cyclic_interval<value_type, open>{zero, infinity},
                        cyclic_interval<value_type, open>{value_type{2, 1u}, value_type{-7, 3u}}));
                    expect(test(cyclic_interval<value_type, open>{zero, singularity_x},
                                cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, closed>{value_type{2, 1u},
                                                                    value_type{-7, 3u}}));
                    expect(test(cyclic_interval<value_type, closed>{zero, singularity_x},
                                cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, closed>{value_type{2, 1u},
                                                                    value_type{-7, 3u}}));
                    expect(test(
                        cyclic_interval<value_type, left_open_right_closed>{zero, singularity_x},
                        cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                        cyclic_interval<value_type, left_open_right_closed>{value_type{2, 1u},
                                                                            value_type{-7, 3u}}));
                    expect(test(
                        cyclic_interval<value_type, left_open_right_closed>{zero, singularity_x},
                        cyclic_interval<value_type, open>{zero, infinity},
                        cyclic_interval<value_type, open>{value_type{2, 1u}, value_type{-7, 3u}}));
                    expect(test(cyclic_interval<value_type, closed>{zero, singularity_x},
                                cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                                cyclic_interval<value_type, left_open_right_closed>{
                                    value_type{2, 1u}, value_type{-7, 3u}}));
                    expect(test(
                        cyclic_interval<value_type, entire>{},
                        cyclic_interval<value_type, open>{zero, infinity},
                        cyclic_interval<value_type, open>{value_type{2, 1u}, value_type{-7, 3u}}));
                    expect(
                        test(cyclic_interval<value_type, open>{zero, singularity_x},
                             cyclic_interval<value_type, single_complement>{zero},
                             cyclic_interval<value_type, single_complement>{value_type{-7, 3u}}));
                    expect(
                        test(cyclic_interval<value_type, closed>{zero, singularity_x},
                             cyclic_interval<value_type, single_complement>{zero},
                             cyclic_interval<value_type, single_complement>{value_type{-7, 3u}}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{}));
                }
            }

            // Generic case with no point in the indeterminacy locus.
            {
                // (x - y)/(xy + 1).
                mapping_type f{// numerator.
                               0, 1, -1, 0,
                               // denominator.
                               1, 0, 0, 1};
                expect(f.kind() == mapping_type::kind_t::generic);
                expect(f.number_of_points_in_indeterminacy_locus() == 0);

                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{2, 1},
                         cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-7, 1}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-9, 13});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{2, 1u},
                         cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{-7, 1u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-9, 13u});

                auto const zero = value_type{0, 1u};
                auto const one = value_type{1, 1u};
                auto const minus_one = value_type{-1, 1u};
                auto const infinity = value_type{1, 0u};

                // Return value of map_cyclic_rectangle.
                {
                    auto test = [&f](auto const& itv_x, auto const& itv_y,
                                     auto const& expected_answer) {
                        auto answer1 = f.map_cyclic_rectangle<false>(itv_x, itv_y);
                        auto answer2 = f.map_cyclic_rectangle<true>(itv_x, itv_y);
                        return answer1 == expected_answer && answer1 == answer2;
                    };

                    expect(test(cyclic_interval<value_type, single_point>{zero},
                                cyclic_interval<value_type, single_point>{zero},
                                cyclic_interval<value_type, single_point>{zero}));
                    expect(test(cyclic_interval<value_type, open>{zero, infinity},
                                cyclic_interval<value_type, open>{zero, infinity},
                                cyclic_interval<value_type, single_complement>{infinity}));
                    expect(test(cyclic_interval<value_type, open>{zero, infinity},
                                cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, single_complement>{infinity}));
                    expect(test(cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                                cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                                cyclic_interval<value_type, open>{zero, infinity},
                                cyclic_interval<value_type, single_complement>{infinity}));
                    expect(test(cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, open>{zero, infinity},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, open>{zero, infinity},
                                cyclic_interval<value_type, single_complement>{zero},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, closed>{zero, infinity},
                                cyclic_interval<value_type, single_complement>{zero},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{}));

                    expect(test(cyclic_interval<value_type, single_point>{one},
                                cyclic_interval<value_type, single_point>{one},
                                cyclic_interval<value_type, single_point>{zero}));
                    expect(test(cyclic_interval<value_type, open>{one, infinity},
                                cyclic_interval<value_type, open>{one, infinity},
                                cyclic_interval<value_type, open>{minus_one, one}));
                    expect(test(cyclic_interval<value_type, open>{one, infinity},
                                cyclic_interval<value_type, closed>{one, infinity},
                                cyclic_interval<value_type, open>{minus_one, one}));
                    expect(test(cyclic_interval<value_type, closed>{one, infinity},
                                cyclic_interval<value_type, closed>{one, infinity},
                                cyclic_interval<value_type, closed>{minus_one, one}));
                    expect(
                        test(cyclic_interval<value_type, left_open_right_closed>{one, infinity},
                             cyclic_interval<value_type, left_closed_right_open>{one, infinity},
                             cyclic_interval<value_type, left_open_right_closed>{minus_one, one}));
                    expect(test(cyclic_interval<value_type, left_open_right_closed>{one, infinity},
                                cyclic_interval<value_type, open>{one, infinity},
                                cyclic_interval<value_type, open>{minus_one, one}));
                    expect(
                        test(cyclic_interval<value_type, closed>{one, infinity},
                             cyclic_interval<value_type, left_closed_right_open>{one, infinity},
                             cyclic_interval<value_type, left_open_right_closed>{minus_one, one}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, open>{one, infinity},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, open>{one, infinity},
                                cyclic_interval<value_type, single_complement>{one},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, closed>{one, infinity},
                                cyclic_interval<value_type, single_complement>{one},
                                cyclic_interval<value_type, entire>{}));
                    expect(test(cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{},
                                cyclic_interval<value_type, entire>{}));
                }
            }

            // Generic case with one point in the indeterminacy locus.
            {
                // xy/(x + y).
                mapping_type f{// numerator.
                               1, 0, 0, 0,
                               // denominator.
                               0, 1, 1, 0};
                expect(f.kind() == mapping_type::kind_t::generic);
                expect(f.number_of_points_in_indeterminacy_locus() == 1);

                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{2, 1},
                         cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-7, 1}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{14, 5});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{2, 1u},
                         cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{-7, 1u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{14, 5u});

                auto const zero = value_type{0, 1u};
                auto const one = value_type{1, 1u};
                auto const minus_one = value_type{-1, 1u};
                auto const infinity = value_type{1, 0u};

                expect(f.is_indeterminacy_locus(zero, zero) == true);
                expect(f.is_indeterminacy_locus(zero, infinity) == false);
                expect(f.is_indeterminacy_locus(infinity, zero) == false);
                expect(f.is_indeterminacy_locus(infinity, infinity) == false);

                // Indeterminacy locus intersection check.
                {
                    // Left-bottom corner.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    // Right-bottom corner.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    // Right-top corner.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, open>{infinity, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);

                    // Left-top corner.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, open>{infinity, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);

                    // Bottom edge.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{zero, one}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, open>{zero, one}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, closed>{zero, one}) == true);

                    // Right edge.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, single_point>{zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, single_point>{zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{minus_one,
                                                                                   one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{minus_one,
                                                                                   one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, single_complement>{infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, entire>{}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, single_point>{zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{minus_one,
                                                                                   one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{minus_one,
                                                                                   one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, single_complement>{infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, entire>{}) == true);

                    // Top edge.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);

                    // Left edge.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, single_point>{zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, single_point>{zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{minus_one,
                                                                                   one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{minus_one,
                                                                                   one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, single_complement>{infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, entire>{}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, single_point>{zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{minus_one,
                                                                                   one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{minus_one,
                                                                                   one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, single_complement>{infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, entire>{}) == true);

                    // Interior.
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, one}) == true);

                    // Exterior.
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{one, minus_one},
                               cyclic_interval<value_type, open>{one, minus_one}) == false);

                    // Others.
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{zero},
                               cyclic_interval<value_type, entire>{}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, single_complement>{zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{zero},
                               cyclic_interval<value_type, single_complement>{zero}) == false);
                }

                // Return value of map_cyclic_rectangle.
                {
                    // Left-bottom corner.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{zero, infinity},
                                               cyclic_interval<value_type, open>{zero, infinity}) ==
                        cyclic_interval<value_type, open>{zero, infinity});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{zero, infinity},
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity}) ==
                        cyclic_interval<value_type, open>{zero, infinity});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{zero, infinity},
                            cyclic_interval<value_type, left_closed_right_open>{zero, infinity}) ==
                        cyclic_interval<value_type, left_closed_right_open>{zero, infinity});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, left_closed_right_open>{zero, infinity});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) ==
                           cyclic_interval<value_type, open>{zero, infinity});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity}) ==
                        cyclic_interval<value_type, left_open_right_closed>{zero, infinity});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                            cyclic_interval<value_type, left_closed_right_open>{zero, infinity}) ==
                        cyclic_interval<value_type, left_closed_right_open>{zero, infinity});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, closed>{zero, infinity});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) ==
                           cyclic_interval<value_type, left_closed_right_open>{zero, infinity});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity}) ==
                        cyclic_interval<value_type, left_closed_right_open>{zero, infinity});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                            cyclic_interval<value_type, left_closed_right_open>{zero, infinity}) ==
                        cyclic_interval<value_type, left_closed_right_open>{zero, infinity});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, left_closed_right_open>{zero, infinity});

                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{zero, infinity},
                                               cyclic_interval<value_type, open>{zero, infinity}) ==
                        cyclic_interval<value_type, left_closed_right_open>{zero, infinity});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, closed>{zero, infinity},
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity}) ==
                        cyclic_interval<value_type, closed>{zero, infinity});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, closed>{zero, infinity},
                            cyclic_interval<value_type, left_closed_right_open>{zero, infinity}) ==
                        cyclic_interval<value_type, left_closed_right_open>{zero, infinity});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, closed>{zero, infinity});

                    // Right-bottom corner.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{infinity, zero},
                                               cyclic_interval<value_type, open>{zero, infinity}) ==
                        cyclic_interval<value_type, single_complement>{zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{infinity, zero},
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity}) ==
                        cyclic_interval<value_type, single_complement>{zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, entire>{});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) ==
                           cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, entire>{});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) ==
                           cyclic_interval<value_type, single_complement>{zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity}) ==
                        cyclic_interval<value_type, single_complement>{zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, entire>{});

                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{infinity, zero},
                                               cyclic_interval<value_type, open>{zero, infinity}) ==
                        cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) ==
                           cyclic_interval<value_type, entire>{});

                    // Right-top corner.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{infinity, zero},
                                               cyclic_interval<value_type, open>{infinity, zero}) ==
                        cyclic_interval<value_type, open>{infinity, zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{infinity, zero},
                            cyclic_interval<value_type, left_open_right_closed>{infinity, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{infinity, zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{infinity, zero},
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero}) ==
                        cyclic_interval<value_type, open>{infinity, zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) ==
                           cyclic_interval<value_type, left_open_right_closed>{infinity, zero});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, open>{infinity, zero}) ==
                           cyclic_interval<value_type, left_open_right_closed>{infinity, zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                            cyclic_interval<value_type, left_open_right_closed>{infinity, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{infinity, zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{infinity, zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) ==
                           cyclic_interval<value_type, left_open_right_closed>{infinity, zero});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, open>{infinity, zero}) ==
                           cyclic_interval<value_type, open>{infinity, zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                            cyclic_interval<value_type, left_open_right_closed>{infinity, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{infinity, zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero}) ==
                        cyclic_interval<value_type, left_closed_right_open>{infinity, zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) ==
                           cyclic_interval<value_type, closed>{infinity, zero});

                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{infinity, zero},
                                               cyclic_interval<value_type, open>{infinity, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{infinity, zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, closed>{infinity, zero},
                            cyclic_interval<value_type, left_open_right_closed>{infinity, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{infinity, zero});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, closed>{infinity, zero},
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero}) ==
                        cyclic_interval<value_type, closed>{infinity, zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) ==
                           cyclic_interval<value_type, closed>{infinity, zero});

                    // Left-top corner.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{zero, infinity},
                                               cyclic_interval<value_type, open>{infinity, zero}) ==
                        cyclic_interval<value_type, single_complement>{zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   infinity, zero}) == cyclic_interval<value_type, entire>{});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{zero, infinity},
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero}) ==
                        cyclic_interval<value_type, single_complement>{zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) ==
                           cyclic_interval<value_type, entire>{});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, open>{infinity, zero}) ==
                           cyclic_interval<value_type, single_complement>{zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   infinity, zero}) == cyclic_interval<value_type, entire>{});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                            cyclic_interval<value_type, left_closed_right_open>{infinity, zero}) ==
                        cyclic_interval<value_type, single_complement>{zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) ==
                           cyclic_interval<value_type, entire>{});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, open>{infinity, zero}) ==
                           cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   infinity, zero}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   infinity, zero}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) ==
                           cyclic_interval<value_type, entire>{});

                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{zero, infinity},
                                               cyclic_interval<value_type, open>{infinity, zero}) ==
                        cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   infinity, zero}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   infinity, zero}) == cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) ==
                           cyclic_interval<value_type, entire>{});

                    // Bottom edge.
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, single_point>{zero},
                                                  cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, single_point>{zero});
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, one},
                                                  cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, entire>{});

                    // Right edge.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, zero},
                                               cyclic_interval<value_type, single_point>{zero}) ==
                        cyclic_interval<value_type, single_point>{zero});
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, zero},
                                               cyclic_interval<value_type, open>{minus_one, one}) ==
                        cyclic_interval<value_type, entire>{});

                    // Top edge.
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, single_point>{zero});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, entire>{});

                    // Left edge.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{zero, one},
                                               cyclic_interval<value_type, single_point>{zero}) ==
                        cyclic_interval<value_type, single_point>{zero});
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{zero, one},
                                               cyclic_interval<value_type, open>{minus_one, one}) ==
                        cyclic_interval<value_type, entire>{});

                    // Interior.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, one},
                                               cyclic_interval<value_type, open>{minus_one, one}) ==
                        cyclic_interval<value_type, entire>{});

                    // Exterior.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{one, minus_one},
                                               cyclic_interval<value_type, open>{one, minus_one}) ==
                        cyclic_interval<value_type, open>{value_type{1, 2u}, value_type{-1, 2u}});

                    // Others.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, single_complement>{zero},
                                               cyclic_interval<value_type, entire>{}) ==
                        cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, single_complement>{zero}) ==
                           cyclic_interval<value_type, entire>{});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, single_complement>{zero},
                               cyclic_interval<value_type, single_complement>{zero}) ==
                           cyclic_interval<value_type, single_complement>{zero});
                }
            }

            // Generic case with two points in the indeterminacy locus.
            {
                // (x - y)/(x + y).
                mapping_type f{// numerator.
                               0, 1, -1, 0,
                               // denominator.
                               0, 1, 1, 0};
                expect(f.kind() == mapping_type::kind_t::generic);
                expect(f.number_of_points_in_indeterminacy_locus() == 2);

                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::int_var>{2, 1},
                         cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-7, 1}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-9, 5});
                expect(f(cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{2, 1u},
                         cntfrc::projective_rational<bigint::int_var, bigint::uint_var>{-7, 1u}) ==
                       cntfrc::projective_rational<bigint::int_var, bigint::int_var>{-9, 5u});

                auto const zero = value_type{0, 1u};
                auto const one = value_type{1, 1u};
                auto const two = value_type{2, 1u};
                auto const minus_one = value_type{-1, 1u};
                auto const infinity = value_type{1, 0u};

                expect(f.is_indeterminacy_locus(zero, zero) == true);
                expect(f.is_indeterminacy_locus(zero, infinity) == false);
                expect(f.is_indeterminacy_locus(infinity, zero) == false);
                expect(f.is_indeterminacy_locus(infinity, infinity) == true);

                // Indeterminacy locus intersection check.
                {
                    // Left-bottom corner.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, open>{zero, one}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, closed>{zero, one}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, closed>{zero, one}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);

                    // Right-bottom corner.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, open>{zero, one}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{zero, one}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{zero, one}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{zero, one}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{zero, one}) == true);

                    // Right-top corner.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);

                    // Left-top corner.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{infinity,
                                                                                   zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{infinity,
                                                                                   zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);

                    // Left-bottom and right-top corners.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, infinity},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    // Right-bottom and left-top corners.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, open>{zero, infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{
                                   zero, infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{zero, infinity}) == true);

                    // Bottom edge.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{zero, one}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, one},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, open>{zero, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, open>{zero, one}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, one},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, closed>{zero, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, closed>{zero, one}) == true);

                    // Right edge.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, single_point>{zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, single_point>{zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{minus_one,
                                                                                   one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{minus_one,
                                                                                   one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, single_complement>{infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, entire>{}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, single_point>{zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{minus_one,
                                                                                   one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{minus_one,
                                                                                   one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, single_complement>{infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, entire>{}) == true);

                    // Top edge.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, open>{minus_one, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, closed>{minus_one, zero}) == true);

                    // Left edge.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, single_point>{zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, single_point>{zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{minus_one,
                                                                                   one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{minus_one,
                                                                                   one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, single_complement>{infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, entire>{}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, single_point>{zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{minus_one,
                                                                                   one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{minus_one,
                                                                                   one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, single_complement>{infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, entire>{}) == true);

                    // Left and right edges.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, open>{two, one}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, open>{two, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{two, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{two, one}) ==
                           false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, closed>{two, one}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, single_complement>{infinity}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{infinity, zero},
                               cyclic_interval<value_type, entire>{}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, open>{two, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_open_right_closed>{two, one}) ==
                           true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, left_closed_right_open>{two, one}) ==
                           true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, closed>{two, one}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, single_complement>{infinity}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{infinity, zero},
                               cyclic_interval<value_type, entire>{}) == true);

                    // Bottom and top edges.
                    expect(f.touches_indeterminacy_locus(
                               cyclic_interval<value_type, open>{two, one},
                               cyclic_interval<value_type, open>{infinity, zero}) == true);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{two, one},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{two, one},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{two, one},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{two, one},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, open>{infinity, zero}) == false);

                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{two, one},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_open_right_closed>{two, one},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, left_closed_right_open>{two, one},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, closed>{two, one},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, single_complement>{infinity},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, entire>{},
                               cyclic_interval<value_type, closed>{infinity, zero}) == true);

                    // Both in interior.
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{two, one},
                               cyclic_interval<value_type, open>{two, one}) == true);

                    // One in interior, another in exterior.
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, one}) == true);

                    // Both in exterior.
                    expect(f.intersects_indeterminacy_locus(
                               cyclic_interval<value_type, open>{one, two},
                               cyclic_interval<value_type, open>{one, two}) == false);
                }

                // Return value of map_cyclic_rectangle.
                {
                    // Left-bottom corner.
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, open>{zero, one},
                                                  cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           cyclic_interval<value_type, open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           cyclic_interval<value_type, left_open_right_closed>{minus_one, one});
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, open>{zero, one},
                                                  cyclic_interval<value_type, closed>{zero, one}) ==
                           cyclic_interval<value_type, left_open_right_closed>{minus_one, one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           cyclic_interval<value_type, open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           cyclic_interval<value_type, left_open_right_closed>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, closed>{zero, one}) ==
                           cyclic_interval<value_type, left_open_right_closed>{minus_one, one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, left_closed_right_open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           cyclic_interval<value_type, left_closed_right_open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           cyclic_interval<value_type, closed>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, closed>{zero, one}) ==
                           cyclic_interval<value_type, closed>{minus_one, one});

                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{zero, one},
                                                  cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, left_closed_right_open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           cyclic_interval<value_type, left_closed_right_open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           cyclic_interval<value_type, closed>{minus_one, one});
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{zero, one},
                                                  cyclic_interval<value_type, closed>{zero, one}) ==
                           cyclic_interval<value_type, closed>{minus_one, one});

                    // Right-bottom corner.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, zero},
                                               cyclic_interval<value_type, open>{zero, one}) ==
                        cyclic_interval<value_type, open>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           cyclic_interval<value_type, open>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           cyclic_interval<value_type, left_closed_right_open>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, zero},
                                               cyclic_interval<value_type, closed>{zero, one}) ==
                        cyclic_interval<value_type, left_closed_right_open>{one, minus_one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, left_open_right_closed>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           cyclic_interval<value_type, left_open_right_closed>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           cyclic_interval<value_type, closed>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{zero, one}) ==
                           cyclic_interval<value_type, closed>{one, minus_one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, open>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           cyclic_interval<value_type, open>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           cyclic_interval<value_type, left_closed_right_open>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{zero, one}) ==
                           cyclic_interval<value_type, left_closed_right_open>{one, minus_one});

                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{minus_one, zero},
                                               cyclic_interval<value_type, open>{zero, one}) ==
                        cyclic_interval<value_type, left_open_right_closed>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_open_right_closed>{zero, one}) ==
                           cyclic_interval<value_type, left_open_right_closed>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, left_closed_right_open>{zero, one}) ==
                           cyclic_interval<value_type, closed>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{minus_one, zero},
                                               cyclic_interval<value_type, closed>{zero, one}) ==
                        cyclic_interval<value_type, closed>{one, minus_one});

                    // Right-top corner.
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, open>{minus_one, one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{minus_one, zero},
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{minus_one, one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{minus_one, zero},
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero}) ==
                        cyclic_interval<value_type, open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, left_open_right_closed>{minus_one, one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, left_closed_right_open>{minus_one, one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero}) ==
                        cyclic_interval<value_type, closed>{minus_one, one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero}) ==
                        cyclic_interval<value_type, left_closed_right_open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, closed>{minus_one, one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, open>{minus_one, one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{minus_one, one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero}) ==
                        cyclic_interval<value_type, open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, left_open_right_closed>{minus_one, one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, left_closed_right_open>{minus_one, one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, closed>{minus_one, zero},
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero}) ==
                        cyclic_interval<value_type, closed>{minus_one, one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, closed>{minus_one, zero},
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero}) ==
                        cyclic_interval<value_type, left_closed_right_open>{minus_one, one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{minus_one, zero},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, closed>{minus_one, one});

                    // Left-top corner.
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, open>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{zero, one},
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero}) ==
                        cyclic_interval<value_type, left_closed_right_open>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, open>{zero, one},
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero}) ==
                        cyclic_interval<value_type, open>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, left_closed_right_open>{one, minus_one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, open>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{zero, one},
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero}) ==
                        cyclic_interval<value_type, left_closed_right_open>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_open_right_closed>{zero, one},
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero}) ==
                        cyclic_interval<value_type, open>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_open_right_closed>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, left_closed_right_open>{one, minus_one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, left_open_right_closed>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{zero, one},
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero}) ==
                        cyclic_interval<value_type, closed>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, left_closed_right_open>{zero, one},
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, left_closed_right_open>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, closed>{one, minus_one});

                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, left_open_right_closed>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, closed>{zero, one},
                            cyclic_interval<value_type, left_open_right_closed>{minus_one, zero}) ==
                        cyclic_interval<value_type, closed>{one, minus_one});
                    expect(
                        f.map_cyclic_rectangle(
                            cyclic_interval<value_type, closed>{zero, one},
                            cyclic_interval<value_type, left_closed_right_open>{minus_one, zero}) ==
                        cyclic_interval<value_type, left_open_right_closed>{one, minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, closed>{zero, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, closed>{one, minus_one});

                    // Bottom edge.
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, single_point>{zero},
                                                  cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, single_point>{minus_one});
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, one},
                                                  cyclic_interval<value_type, open>{zero, one}) ==
                           cyclic_interval<value_type, single_complement>{one});
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, one},
                                                  cyclic_interval<value_type, closed>{zero, one}) ==
                           cyclic_interval<value_type, entire>{});

                    // Right edge.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, zero},
                                               cyclic_interval<value_type, single_point>{zero}) ==
                        cyclic_interval<value_type, single_point>{one});
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, zero},
                                               cyclic_interval<value_type, open>{minus_one, one}) ==
                        cyclic_interval<value_type, single_complement>{minus_one});
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{minus_one, zero},
                                               cyclic_interval<value_type, open>{minus_one, one}) ==
                        cyclic_interval<value_type, entire>{});

                    // Top edge.
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, single_point>{zero},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, single_point>{minus_one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, open>{minus_one, zero}) ==
                           cyclic_interval<value_type, single_complement>{one});
                    expect(f.map_cyclic_rectangle(
                               cyclic_interval<value_type, open>{minus_one, one},
                               cyclic_interval<value_type, closed>{minus_one, zero}) ==
                           cyclic_interval<value_type, entire>{});

                    // Left edge.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{zero, one},
                                               cyclic_interval<value_type, single_point>{zero}) ==
                        cyclic_interval<value_type, single_point>{one});
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{zero, one},
                                               cyclic_interval<value_type, open>{minus_one, one}) ==
                        cyclic_interval<value_type, single_complement>{minus_one});
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, closed>{zero, one},
                                               cyclic_interval<value_type, open>{minus_one, one}) ==
                        cyclic_interval<value_type, entire>{});

                    // Both in interior.
                    expect(f.map_cyclic_rectangle(cyclic_interval<value_type, open>{two, one},
                                                  cyclic_interval<value_type, open>{two, one}) ==
                           cyclic_interval<value_type, entire>{});

                    // One in interior, another in exterior.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{minus_one, one},
                                               cyclic_interval<value_type, open>{minus_one, one}) ==
                        cyclic_interval<value_type, entire>{});

                    // Both in exterior.
                    expect(
                        f.map_cyclic_rectangle(cyclic_interval<value_type, open>{one, two},
                                               cyclic_interval<value_type, open>{one, two}) ==
                        cyclic_interval<value_type, open>{value_type{-1, 3u}, value_type{1, 3u}});
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
