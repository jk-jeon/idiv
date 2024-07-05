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

#include <idiv/best_rational_approx.h>
#include <idiv/bigint.h>
#include <idiv/rational_continued_fraction.h>
#include <boost/ut.hpp>

void best_rational_approx_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Best rational approximation]"_test = [] {
        using rational_t = frac<bigint::int_var, bigint::uint_var>;

        should("best_rational_approx") = [] {
            auto perform_test = [](bigint::int_var const& numerator,
                                   bigint::uint_var const& denominator, std::size_t nmax) {
                auto cf = cntfrc::make_generator<cntfrc::index_tracker,
                                                 cntfrc::previous_previous_convergent_tracker>(
                    cntfrc::impl::rational{numerator, denominator});

                auto result = idiv::find_best_rational_approx(cf, nmax);

                auto from_below = rational_t{util::div_floor(numerator, denominator), 1u};
                auto from_above = rational_t{util::div_ceil(numerator, denominator), 1u};
                for (std::size_t i = 1; i <= nmax; ++i) {
                    auto low = rational_t{util::div_floor(i * numerator, denominator), unsigned(i)};
                    auto high = rational_t{util::div_ceil(i * numerator, denominator), unsigned(i)};

                    if (from_below < low) {
                        from_below = low;
                    }
                    if (high < from_above) {
                        from_above = high;
                    }
                }

                expect(result.below == from_below);
                expect(result.above == from_above);
            };
            // Effectively rational case.
            perform_test(137, 1290u, 1500);
            // Effectively irrational case.
            perform_test(6614777, 12961230u, 1500);
        };

        should("find_floor_quotient_range") = [] {
            auto perform_test = [](bigint::int_var const& numerator,
                                   bigint::uint_var const& denominator, std::size_t nmax) {
                auto cf = cntfrc::make_generator<cntfrc::index_tracker,
                                                 cntfrc::previous_previous_convergent_tracker>(
                    cntfrc::impl::rational{numerator, denominator});

                auto result = idiv::find_floor_quotient_range(cf, nmax);

                auto from_below = rational_t{util::div_floor(numerator, denominator), 1u};
                auto from_above = rational_t{util::div_floor(numerator, denominator) + 1, 1u};
                for (std::size_t i = 1; i <= nmax; ++i) {
                    auto low = rational_t{util::div_floor(i * numerator, denominator), unsigned(i)};
                    auto high =
                        rational_t{util::div_floor(i * numerator, denominator) + 1, unsigned(i)};

                    if (from_below < low) {
                        from_below = low;
                    }
                    if (high < from_above) {
                        from_above = high;
                    }
                }

                expect(result.lower_bound() == from_below);
                expect(result.upper_bound() == from_above);
            };
            // Effectively rational case.
            perform_test(137, 1290u, 1500);
            // Effectively irrational case.
            perform_test(6614777, 12961230u, 1500);
        };

        should("find_extremizers_of_fractional_part") = [] {
            auto perform_test = [](bigint::int_var const& numerator,
                                   bigint::uint_var const& denominator, std::size_t nmax) {
                auto cf = cntfrc::make_generator<cntfrc::index_tracker,
                                                 cntfrc::previous_previous_convergent_tracker>(
                    cntfrc::impl::rational{numerator, denominator});

                auto result = idiv::find_extremizers_of_fractional_part(cf, nmax);
                expect(result.smallest_minimizer >= 1 && result.smallest_minimizer <= nmax);
                expect(result.largest_maximizer >= 1 && result.largest_maximizer <= nmax);

                auto smallest_remainder = denominator - 1u;
                auto smallest_minimizer = nmax;
                auto largest_remainder = bigint::uint_var{0u};
                auto largest_maximizer = std::size_t(1);
                for (std::size_t i = 1; i <= nmax; ++i) {
                    auto remainder = i * numerator;
                    remainder.long_division(denominator);

                    if (smallest_remainder > remainder) {
                        smallest_remainder = util::abs(remainder);
                        smallest_minimizer = i;
                    }
                    if (largest_remainder <= remainder) {
                        largest_remainder = util::abs(remainder);
                        largest_maximizer = i;
                    }
                }

                expect(result.smallest_minimizer == smallest_minimizer);
                expect(result.largest_maximizer == largest_maximizer);
            };
            // Effectively rational case.
            perform_test(137, 1290u, 1500);
            // Effectively irrational case.
            perform_test(6614777, 12961230u, 1500);
        };
    };
}
