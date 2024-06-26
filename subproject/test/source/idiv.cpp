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

#include <idiv/idiv.h>
#include <boost/ut.hpp>
#include <format>
#include <iostream>

namespace jkj {
    namespace bigint {
        static std::ostream& operator<<(std::ostream& out, jkj::bigint::uint_var const& n) {
            auto digits = n.to_decimal();
            if (digits.empty()) {
                out << "0";
            }
            else {
                auto itr = digits.cbegin();
                out << std::format("{}", *itr);

                for (++itr; itr != digits.cend(); ++itr) {
                    out << std::format("{:019d}", *itr);
                }
            }

            return out;
        }

        static std::ostream& operator<<(std::ostream& out, jkj::bigint::int_var const& n) {
            if (n.is_strictly_negative()) {
                out << "-";
            }
            out << n.abs();
            return out;
        }
    }
}

void idiv_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[idiv]"_test = [] {
        using projective_rational_t =
            cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using nrange_t = interval<bigint::int_var, interval_type_t::bounded_closed>;

        should("find_optimal_multiply_shift") = [] {
            auto perform_test = [](bigint::int_var const& numerator,
                                   bigint::uint_var const& denominator, std::size_t nmax) {
                auto cf = cntfrc::make_generator<cntfrc::index_tracker,
                                                 cntfrc::previous_previous_convergent_tracker>(
                    cntfrc::impl::rational{numerator, denominator});

                auto result = idiv::find_optimal_multiply_shift(cf, nmax);

                std::size_t k = 0;
                bigint::uint_var m = 0u;
                while (true) {
                    bool success = true;
                    for (std::size_t n = 1; n <= nmax; ++n) {
                        if (util::div_floor(n * numerator, denominator) != ((n * m) >> k)) {
                            success = false;
                            break;
                        }
                    }
                    if (!success) {
                        if (++m == bigint::uint_var::power_of_2(k)) {
                            m = 0u;
                            ++k;
                        }
                    }
                    else {
                        break;
                    }
                }
                expect(k == result.shift_amount && m == result.multiplier);
            };
            // Effectively rational case.
            perform_test(17, 129u, 150);
            // Effectively irrational case.
            perform_test(6614777, 12961230u, 150);
        };

        should("find_simultaneous_multiply_add_shift") = [] {
            auto perform_test = [](projective_rational_t const& x, projective_rational_t const& y,
                                   nrange_t const& nrange) {
                auto xcf = cntfrc::make_generator(cntfrc::impl::rational{x});
                auto ycf = cntfrc::make_generator(cntfrc::impl::rational{y});

                auto result = idiv::find_simultaneous_multiply_add_shift(xcf, ycf, nrange);

                auto const xdyd = x.denominator * y.denominator;
                auto const xnyd = x.numerator * y.denominator;
                auto const ynxd = y.numerator * x.denominator;
                for (auto n = nrange.lower_bound(); n <= nrange.upper_bound(); ++n) {
                    auto true_answer = util::div_floor(xnyd * n + ynxd, xdyd);
                    auto answer = (n * result.multiplier + result.adder) >> result.shift_amount;

                    expect(answer == true_answer) << "n = " << n;
                }
            };
            perform_test(projective_rational_t{17, 129u}, projective_rational_t{39, 176u},
                         nrange_t{-150, 150});
            perform_test(projective_rational_t{1'936'274, 6'432'163u},
                         projective_rational_t{-4'206'456, 33'668'149u}, nrange_t{-1000, 1000});
            // Case when the shift amount should be increased.
            perform_test(projective_rational_t{432, 611u}, projective_rational_t{46, 47u},
                         nrange_t{0, 273});
        };

        should("find_extrema_of_fractional_part (two unknowns)") = [] {
            auto perform_test = [](projective_rational_t const& x, projective_rational_t const& y,
                                   nrange_t const& nrange) {
                auto xcf = cntfrc::make_generator(cntfrc::impl::rational{x});
                auto ycf = cntfrc::make_generator(cntfrc::impl::rational{y});

                auto result = idiv::find_extrema_of_fractional_part(xcf, ycf, nrange);
                expect(result.smallest_minimizer >= nrange.lower_bound() &&
                       result.smallest_minimizer <= nrange.upper_bound());
                expect(result.largest_maximizer >= nrange.lower_bound() &&
                       result.largest_maximizer <= nrange.upper_bound());

                auto const xdyd = x.denominator * y.denominator;
                auto const xnyd = x.numerator * y.denominator;
                auto const ynxd = y.numerator * x.denominator;
                auto smallest_remainder = xdyd - 1u;
                auto smallest_minimizer = nrange.upper_bound();
                auto largest_remainder = bigint::uint_var{0u};
                auto largest_maximizer = nrange.lower_bound();
                for (bigint::int_var n = nrange.lower_bound(); n <= nrange.upper_bound(); ++n) {
                    auto remainder = n * xnyd + ynxd;
                    remainder.long_division(xdyd);

                    if (smallest_remainder > remainder) {
                        smallest_remainder = util::abs(remainder);
                        smallest_minimizer = n;
                    }
                    if (largest_remainder <= remainder) {
                        largest_remainder = util::abs(remainder);
                        largest_maximizer = n;
                    }
                }

                expect(result.smallest_minimizer == smallest_minimizer);
                expect(result.largest_maximizer == largest_maximizer);
            };
            perform_test(projective_rational_t{17, 129u}, projective_rational_t{39, 176u},
                         nrange_t{-150, 150});
            perform_test(projective_rational_t{1'936'274, 6'432'163u},
                         projective_rational_t{-4'206'456, 33'668'149u}, nrange_t{-1000, 1000});
        };

        should("find_maximizer/minimizer_of_floor_subtract_quotient_positive_range") = [] {
            auto perform_test = [](projective_rational_t const& x, projective_rational_t const& y,
                                   projective_rational_t const& zeta, nrange_t const& nrange) {
                auto xcf = cntfrc::make_generator(cntfrc::impl::rational{x});
                auto ycf = cntfrc::make_generator(cntfrc::impl::rational{y});
                auto zetacf = cntfrc::make_generator(cntfrc::impl::rational{zeta});

                auto const maximizer_computed =
                    idiv::find_maximizer_of_floor_subtract_quotient_positive_range(xcf, ycf, zetacf,
                                                                                   nrange);
                expect(maximizer_computed >= nrange.lower_bound() &&
                       maximizer_computed <= nrange.upper_bound());

                auto const minimizer_computed =
                    idiv::find_minimizer_of_floor_subtract_quotient_positive_range(xcf, ycf, zetacf,
                                                                                   nrange);
                expect(minimizer_computed >= nrange.lower_bound() &&
                       minimizer_computed <= nrange.upper_bound());

                auto const xdyd = x.denominator * y.denominator;
                auto const xnyd = x.numerator * y.denominator;
                auto const ynxd = y.numerator * x.denominator;

                auto maximizer = maximizer_computed;
                auto max_value = make_frac_from_signed(
                    zeta.denominator * util::div_floor(maximizer * xnyd + ynxd, xdyd) -
                        zeta.numerator,
                    maximizer * zeta.denominator);

                auto minimizer = minimizer_computed;
                auto min_value = make_frac_from_signed(
                    zeta.denominator * util::div_floor(minimizer * xnyd + ynxd, xdyd) -
                        zeta.numerator,
                    minimizer * zeta.denominator);

                for (bigint::int_var n = nrange.lower_bound(); n <= nrange.upper_bound(); ++n) {
                    auto const value = make_frac_from_signed(
                        zeta.denominator * util::div_floor(n * xnyd + ynxd, xdyd) - zeta.numerator,
                        n * zeta.denominator);

                    if (value > max_value) {
                        maximizer = n;
                        max_value = value;
                    }

                    if (value < min_value) {
                        minimizer = n;
                        min_value = value;
                    }
                }

                expect(maximizer_computed == maximizer);
                expect(minimizer_computed == minimizer);
            };
            perform_test(projective_rational_t{19'282, 23'129u}, projective_rational_t{98, 519u},
                         projective_rational_t{661, 8'331u}, nrange_t{5'123, 12'150});
        };
    };
}