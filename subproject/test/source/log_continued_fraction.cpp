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
#include <idiv/gosper_continued_fraction.h>
#include <idiv/log_continued_fraction.h>
#include <boost/ut.hpp>

void log_continued_fraction_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Continued fractions for natural and general log]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using unsigned_frac_t = frac<bigint::uint_var, bigint::uint_var>;

        should("natural_log_calculator") = [] {
            unsigned_frac_t error_bound{
                1u, bigint::uint_var(bigint::decimal_uint_const_v<1'000'000, 0, 0, 0, 0, 0>)};

            {
                // Compute ln(2) up to 100 digits.
                auto nlc =
                    cntfrc::make_generator<cntfrc::index_tracker, cntfrc::partial_fraction_tracker,
                                           cntfrc::convergent_tracker, cntfrc::interval_tracker>(
                        cntfrc::impl::natural_log_calculator<bigint::int_var, bigint::uint_var>{
                            unsigned_frac_t{2u, 1u}});

                auto const approx_ln2 = nlc.progress_until(error_bound);
                auto const digits = div_floor(
                    approx_ln2.numerator * bigint::decimal_uint_const_v<100'000, 0, 0, 0, 0, 0>,
                    approx_ln2.denominator);

                expect(
                    digits ==
                    bigint::decimal_uint_const_v<
                        69'314, UINT64_C(7'180'559'945'309'417'232),
                        UINT64_C(1'214'581'765'680'755'001), UINT64_C(3'436'025'525'412'068'000),
                        UINT64_C(9'493'393'621'969'694'715), UINT64_C(6'058'633'269'964'186'875)>);
            }
        };

        should("natural_log") = [] {
            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                cntfrc::impl::natural_log<bigint::int_var, bigint::uint_var>{
                    unsigned_frac_t{3u, 1u}});

            // First 15 convergents of ln(3).
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{1, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{11, 10u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{78, 71u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{713, 649u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{1'504, 1'369u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{3'721, 3'387u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{5'225, 4'756u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{19'396, 17'655u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{24'621, 22'411u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{807'268, 734'807u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{1'639'157, 1'492'025u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{28'672'937, 26'099'232u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{30'312'094, 27'591'257u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{483'354'347, 439'968'087u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{513'666'441, 467'559'344u});
        };

        should("general_log") = [] {
            // Rational cases.
            {
                // log(175616/91125) / log(3136/2025)
                // = log(2^9*7^4 / 3^6*5^3) / log(2^6*7^2 / 3^4*5^2)
                // = 3/2.
                auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                    cntfrc::impl::general_log<bigint::int_var, bigint::uint_var>{
                        unsigned_frac_t{3136u, 2025u}, unsigned_frac_t{175'616u, 91'125u}});

                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{1, 1u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{3, 2u});
                expect(cf.update() == false);
            }
            {
                // c = 151238619319231523311 / 51098098609801923098592931,
                // 151238619319231523311 = 8 * 2^64 + 3664666729555110383,
                // 51098098609801923098592931 = 2770033 * 2^64 + 8783072032707069603.
                // a = c^17, b=c^9, so log_a(b) = 9/17.
                auto c = unsigned_frac_t{
                    bigint::uint_var{UINT64_C(8), UINT64_C(3664666729555110383)},
                    bigint::uint_var{UINT64_C(2770033), UINT64_C(8783072032707069603)}};
                auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                    cntfrc::impl::general_log<bigint::int_var, bigint::uint_var>{
                        util::pow_uint(c, 17u), util::pow_uint(c, 9u)});

                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{0, 1u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{1, 1u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{1, 2u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{9, 17u});
                expect(cf.update() == false);
            }
            // Irrational case.
            {
                auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                    cntfrc::impl::general_log<bigint::int_var, bigint::uint_var>{
                        unsigned_frac_t{2u, 1u}, unsigned_frac_t{4u, 3u}});

                // First 20 convergents of log2(4/3).
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{0, 1u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{1, 2u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{2, 5u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{5, 12u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{17, 41u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{22, 53u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{127, 306u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{276, 665u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{6'475, 15'601u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{13'226, 31'867u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{32'927, 79'335u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{46'153, 111'202u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{79'080, 190'537u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{4'395'553, 10'590'737u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{4'474'633, 10'781'274u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{22'294'085, 53'715'833u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{71'356'888, 171'928'773u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{93'650'973, 225'644'606u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{165'007'861, 397'573'379u});
                expect(cf.update() == true);
                expect(cf.current_convergent() ==
                       convergent_t{INT64_C(2'568'768'888), UINT64_C(6'189'245'291)});
            }
        };

        should("additional_unary_gosper") = [] {
            using continued_fraction_t = cntfrc::impl::unary_gosper<
                cntfrc::impl::natural_log_calculator<bigint::int_var, bigint::uint_var>>;

            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                continued_fraction_t{continued_fraction_t::internal_continued_fraction_impl_type{
                                         unsigned_frac_t{4u, 3u}},
                                     {0, 7, 2, 1}});

            // First 15 convergents of 7/(2ln(4/3) + 1).
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{4, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{9, 2u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{31, 7u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{40, 9u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{471, 106u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{1'924, 433u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{10'091, 2'271u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{12'015, 2'704u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{34'121, 7'679u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{592'072, 133'247u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{4'770'697, 1'073'655u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{5'362'769, 1'206'902u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{20'859'004, 4'694'361u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{276'529'821, 62'233'595u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{297'388'825, 66'927'956u});
        };

        should("additional_binary_gosper") = [] {
            using log_calculator =
                cntfrc::impl::natural_log_calculator<bigint::int_var, bigint::uint_var>;
            using continued_fraction_t =
                cntfrc::impl::binary_gosper<log_calculator, log_calculator>;

            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                continued_fraction_t{log_calculator{unsigned_frac_t{176u, 39u}},
                                     log_calculator{unsigned_frac_t{95u, 771u}},
                                     {// numerator
                                      0, 0, -4, 1,
                                      // denominator
                                      7, 3, -1, 0}});

            // First 20 convergents of
            // (-4ln(95/771) + 1)/(7ln(176/39)ln(95/771) + 3ln(176/39) - ln(95/771)).
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-1, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-1, 2u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-2, 3u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-3, 5u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-20, 33u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-163, 269u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-346, 571u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-509, 840u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-855, 1411u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-2'219, 3'662u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-3'074, 5'073u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-75'995, 125'414u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-79'069, 130'487u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-155'064, 255'901u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-1'474'645, 2'433'596u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-7'528'289, 12'423'881u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-9'002'934, 14'857'477u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-16'531'223, 27'281'358u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-25'534'157, 42'138'835u});
        };
    };
}