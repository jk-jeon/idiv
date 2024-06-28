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
#include <boost/ut.hpp>

void bigint_test() {
    using namespace boost::ut;
    using namespace jkj;

    "[Arbitrary-precision integers (bigint)]"_test = [] {
        should("static_block_holder utils") = [] {
            bigint::static_block_holder<7> a{1, 2, 3, 4, 5, 6, 7};
            expect(std::is_same_v<decltype(bigint::detail::slice<4>(a)),
                                  bigint::static_block_holder<4>> &&
                   bigint::detail::slice<4>(a)[3] == 4);

            expect(bigint::detail::count_blocks_excluding_leading_zeros(
                       bigint::static_block_holder<16>{1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0,
                                                       0}) == 9);

            auto b = bigint::detail::remove_leading_zero_blocks<bigint::static_block_holder<16>{
                1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 0, 0, 0}>();
            expect(std::is_same_v<decltype(b), bigint::static_block_holder<9>> && b[5] == 6);

            expect(bigint::detail::reverse_and_validate(b)[2] == 7);
        };

        should("uint_var comparisons") = [] {
            expect(bigint::uint_var{1, 2, 3, 4} == bigint::uint_var{0, 0, 0, 0, 1, 2, 3, 4});
            expect(bigint::uint_var{1, 2, 3, 4} != bigint::uint_var{1, 2, 3, 3});
            expect(bigint::uint_var{1, 2, 3, 4} > bigint::uint_var{1, 2, 3, 3});
            expect(bigint::uint_var{1, 2, 3, 4} > bigint::uint_var{1, 2, 3});
            expect(bigint::uint_var{1, 3, 3} > bigint::uint_var{1, 2, 6});
            expect(bigint::uint_var{1, 2, 3} < bigint::uint_var{1, 2, 3, 3});
            expect(bigint::uint_var{1, 2, 3} < bigint::uint_var{1, 2, 5});
            expect(bigint::uint_var{1, 5, 3} < bigint::uint_var{1, 7, 1});
            expect(bigint::uint_var{1, 1, 1} <= bigint::uint_var{1, 1, 1});
            expect(bigint::uint_var{1, 1, 1} <= bigint::uint_var{3, 1, 1});
            expect(bigint::uint_var{1, 2, 3} <= bigint::uint_var{1, 3, 1});
            expect(bigint::uint_var{11, 61, 31} >= bigint::uint_var{11, 61, 30});
            expect(bigint::uint_var{11, 91, 31} >= bigint::uint_var{11, 51, 51});
            expect(bigint::uint_var{116, 76, 31} >= bigint::uint_var{116, 76, 31});

            expect(bigint::uint_var{0, 0, 0, 12345678} == 12345678u);
            expect(bigint::uint_var{0, 0, 0, 123456789} != 12345678u);
            expect(bigint::uint_var{1, 12345678} != 12345678u);
            expect(bigint::uint_var{} == 0u);
            expect(bigint::uint_var{} != 1u);
            expect(bigint::uint_var{1, 2} != 1u);
        };

        should("uint_var member functions") = [] {
            expect(bigint::uint_var{0, 0, 0, 0, 1, 2, 3, 4}.number_of_blocks() == 4);
            expect(bigint::uint_var{0, 0, 0, 0, 1, 2, 3, 4}[3] ==
                   1); // Recall the LSB is located first.
            expect(bigint::uint_var{}.is_zero());
            expect(!bigint::uint_var{1, 2, 3}.is_zero());
            {
                auto x = bigint::uint_var{UINT64_C(0xffff'ffff'ffff'ffff),
                                          UINT64_C(0xffff'ffff'ffff'ffff),
                                          UINT64_C(0xffff'ffff'ffff'ffff)};
                x += bigint::uint_var{UINT64_C(0xffff'ffff'ffff'ffff),
                                      UINT64_C(0xffff'ffff'ffff'ffff)};
                expect(x == bigint::uint_const_v<1, 0, UINT64_C(0xffff'ffff'ffff'ffff),
                                                 UINT64_C(0xffff'ffff'ffff'fffe)>);
            }
            {
                auto x = bigint::uint_var{UINT64_C(0xffff'ffff'ffff'ffff),
                                          UINT64_C(0xffff'ffff'ffff'ffff)};
                x += bigint::uint_const_v<UINT64_C(0xffff'ffff'ffff'ffff),
                                          UINT64_C(0xffff'ffff'ffff'ffff),
                                          UINT64_C(0xffff'ffff'ffff'ffff)>;
                expect(x == bigint::uint_const_v<1, 0, UINT64_C(0xffff'ffff'ffff'ffff),
                                                 UINT64_C(0xffff'ffff'ffff'fffe)>);
            }
            {
                auto x = bigint::uint_var{
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                    UINT64_C(0xffff'ffff'ffff'f000)};
                x += 0x1000;
                expect(x == bigint::uint_const_v<1, 0, 0, 0, 0, 0>);
            }
            {
                auto x = bigint::uint_var{};
                ++x;
                expect(x == 1u);
                expect(x++ == 1u);
                expect(x == 2u);
            }
            {
                auto x = bigint::uint_var{
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff)};
                ++x;
                expect(x == bigint::uint_const_v<1, 0, 0, 0, 0, 0, 0, 0, 0>);
                expect(x++ == bigint::uint_const_v<1, 0, 0, 0, 0, 0, 0, 0, 0>);
                expect(x == bigint::uint_const_v<1, 0, 0, 0, 0, 0, 0, 0, 1>);
            }
            {
                auto x = bigint::uint_var{1, 0, 0, 0, 0, 0, 0, 0, 0};
                x -= UINT64_C(0xffff'ffff'ffff'ffff);
                expect(x == bigint::uint_const_v<
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), 1>);
                x -= UINT64_C(0xffff'ffff'ffff'ffff);
                expect(x == bigint::uint_const_v<
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'fffe), 2>);
                x -= bigint::uint_const_v<UINT64_C(0xffff'ffff'ffff'ffff),
                                          UINT64_C(0xffff'ffff'ffff'ffff)>;
                expect(x == bigint::uint_const_v<
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'fffe),
                                UINT64_C(0xffff'ffff'ffff'fffe), 3>);
                x -= bigint::uint_const_v<
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'fffe), 0, 3>;
                expect(x == bigint::uint_const_v<UINT64_C(0xffff'ffff'ffff'fffe), 0>);
                x -= bigint::uint_const_v<UINT64_C(0xffff'ffff'ffff'fffd), 1>;
                expect(x == UINT64_C(0xffff'ffff'ffff'ffff));
                x -= UINT64_C(0xffff'ffff'ffff'ffff);
                expect(x == 0u);
            }
            {
                auto x = bigint::uint_var{1, 0, 0, 0, 0, 0, 0, 0, 0};
                --x;
                expect(x == bigint::uint_const_v<
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff)>);
                expect(x-- ==
                       bigint::uint_const_v<
                           UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                           UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                           UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                           UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff)>);
                expect(x == bigint::uint_const_v<
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                                UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'fffe)>);
            }
            {
                auto x = bigint::uint_var{
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff),
                    UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff)};
                x *= UINT64_C(0xffff'ffff'ffff'ffff);
                expect(x ==
                       bigint::uint_const_v<
                           UINT64_C(0xffff'ffff'ffff'fffe), UINT64_C(0xffff'ffff'ffff'ffff),
                           UINT64_C(0xffff'ffff'ffff'ffff), UINT64_C(0xffff'ffff'ffff'ffff), 1>);
                auto y = bigint::uint_var{1, 0};
                y *= UINT64_C(0xffff'ffff'ffff'ffff);
                expect(y == bigint::uint_const_v<UINT64_C(0xffff'ffff'ffff'ffff), 0>);
            }
            {
                auto x =
                    bigint::uint_var{UINT64_C(0x12345678'12345678), UINT64_C(0x12345678'12345678),
                                     UINT64_C(0x87654321'87654321)};
                auto q = x.long_division(bigint::uint_const_v<UINT64_C(0x87654321'87654321),
                                                              UINT64_C(0x12345678'12345678)>);
                expect(q == UINT64_C(0x226b'9021'46b3'cd4f));
                expect(x == bigint::uint_const_v<UINT64_C(0x3aad'98ad'581b'da5e),
                                                 UINT64_C(0xd786'bdca'ba18'7c19)>);
            }
            // Not a member but a hidden friend.
            expect(bit_width(bigint::uint_var{UINT64_C(0x0000'8000'0000'0000),
                                              UINT64_C(12345678900123456789)}) == 112);

            expect(bigint::uint_var::power_of_2(175) ==
                   bigint::uint_const_v<UINT64_C(0x0000'8000'0000'0000), 0, 0>);

            {
                auto x = bigint::uint_var{UINT64_C(0x0000'1234'dead'beef),
                                          UINT64_C(0xdead'beef'dead'beef),
                                          UINT64_C(0xdead'beef'dead'beef)};
                x <<= 4;
                expect(x == bigint::uint_const_v<UINT64_C(0x0001'234d'eadb'eefd),
                                                 UINT64_C(0xeadb'eefd'eadb'eefd),
                                                 UINT64_C(0xeadb'eefd'eadb'eef0)>);
                x <<= 12;
                expect(x == bigint::uint_const_v<UINT64_C(0x1234'dead'beef'dead),
                                                 UINT64_C(0xbeef'dead'beef'dead),
                                                 UINT64_C(0xbeef'dead'beef'0000)>);
                x <<= 128;
                expect(x == bigint::uint_const_v<UINT64_C(0x1234'dead'beef'dead),
                                                 UINT64_C(0xbeef'dead'beef'dead),
                                                 UINT64_C(0xbeef'dead'beef'0000), 0, 0>);
                x <<= 16;
                expect(x == bigint::uint_const_v<0x1234, UINT64_C(0xdead'beef'dead'beef),
                                                 UINT64_C(0xdead'beef'dead'beef),
                                                 UINT64_C(0xdead'beef'0000'0000), 0, 0>);
                x >>= 144;
                expect(x == bigint::uint_const_v<UINT64_C(0x1234'dead'beef'dead),
                                                 UINT64_C(0xbeef'dead'beef'dead),
                                                 UINT64_C(0xbeef'dead'beef'0000)>);
                x >>= 192;
                expect(x == 0u);
            }

            {
                auto x = bigint::uint_var{
                    UINT64_C(0x0000'0012'dead'beef), UINT64_C(0xdead'beef'dead'beef),
                    UINT64_C(0xdead'beef'dead'beef), UINT64_C(0xdead'beef'0000'0000),
                    UINT64_C(0x0000'0000'0000'0000), UINT64_C(0x0000'0000'0000'0000),
                    UINT64_C(0x0000'0000'0000'0000), UINT64_C(0x0000'0000'0000'0000)};
                expect(x.factor_out_power_of_2() == 288);
                expect(x == bigint::uint_const_v<0x12, UINT64_C(0xdead'beef'dead'beef),
                                                 UINT64_C(0xdead'beef'dead'beef),
                                                 UINT64_C(0xdead'beef'dead'beef)>);
            }

            {
                auto x =
                    bigint::uint_var{UINT64_C(1234567890123456789), UINT64_C(1234567890123456789),
                                     UINT64_C(1234567890123456789), UINT64_C(1234567890123456789),
                                     UINT64_C(1234567890123456789), UINT64_C(1234567890123456789),
                                     UINT64_C(1234567890123456789), UINT64_C(1234567890123456789),
                                     UINT64_C(1234567890123456789), UINT64_C(1234567890123456789),
                                     UINT64_C(1234567890123456789), UINT64_C(1234567890123456789),
                                     UINT64_C(1234567890123456789), UINT64_C(1234567890123456789),
                                     UINT64_C(1234567890123456789)};
                auto decimal = std::vector<bigint::block_type>{652,
                                                               UINT64_C(2'154'646'679'190'934'311),
                                                               UINT64_C(6'475'070'707'159'753'481),
                                                               UINT64_C(9'644'072'700'111'562'292),
                                                               UINT64_C(6'078'271'397'402'072'764),
                                                               UINT64_C(7'177'842'004'041'080'699),
                                                               UINT64_C(1'294'679'995'335'186'520),
                                                               UINT64_C(3'227'147'820'766'060'402),
                                                               UINT64_C(7'809'184'506'609'574'266),
                                                               UINT64_C(6'140'449'254'819'569'255),
                                                               UINT64_C(3'017'536'223'147'787'056),
                                                               UINT64_C(8'218'105'215'718'557'812),
                                                               UINT64_C(2'032'926'231'820'575),
                                                               UINT64_C(6'519'632'540'758'027'403),
                                                               UINT64_C(1'665'352'115'267'003'556),
                                                               UINT64_C(6'326'492'966'752'452'885)};
                expect(x.to_decimal() == decimal);
                expect(bigint::uint_var::from_decimal(decimal) == x);
                expect(
                    x ==
                    bigint::decimal_uint_const_v<
                        652, UINT64_C(2'154'646'679'190'934'311),
                        UINT64_C(6'475'070'707'159'753'481), UINT64_C(9'644'072'700'111'562'292),
                        UINT64_C(6'078'271'397'402'072'764), UINT64_C(7'177'842'004'041'080'699),
                        UINT64_C(1'294'679'995'335'186'520), UINT64_C(3'227'147'820'766'060'402),
                        UINT64_C(7'809'184'506'609'574'266), UINT64_C(6'140'449'254'819'569'255),
                        UINT64_C(3'017'536'223'147'787'056), UINT64_C(8'218'105'215'718'557'812),
                        UINT64_C(2'032'926'231'820'575), UINT64_C(6'519'632'540'758'027'403),
                        UINT64_C(1'665'352'115'267'003'556), UINT64_C(6'326'492'966'752'452'885)>);
            }
        };

        should("unsigned binary operations") = [] {
            {
                auto x = bigint::uint_var{
                    UINT64_C(0x7777'7777'7777'7777), UINT64_C(0x7777'7777'7777'7777),
                    UINT64_C(0x7777'7777'7777'7777), UINT64_C(0x7777'7777'7777'7777)};
                auto y = bigint::uint_var{
                    UINT64_C(0x4444'4444'4444'4444), UINT64_C(0x4444'4444'4444'4444),
                    UINT64_C(0x4444'4444'4444'4444), UINT64_C(0x4444'4444'4444'4444)};
                auto sum = bigint::uint_var{
                    UINT64_C(0xbbbb'bbbb'bbbb'bbbb), UINT64_C(0xbbbb'bbbb'bbbb'bbbb),
                    UINT64_C(0xbbbb'bbbb'bbbb'bbbb), UINT64_C(0xbbbb'bbbb'bbbb'bbbb)};
                auto diff = bigint::uint_var{
                    UINT64_C(0x3333'3333'3333'3333), UINT64_C(0x3333'3333'3333'3333),
                    UINT64_C(0x3333'3333'3333'3333), UINT64_C(0x3333'3333'3333'3333)};
                auto sum1 = bigint::uint_var{
                    UINT64_C(0x7777'7777'7777'7777), UINT64_C(0x7777'7777'7777'7777),
                    UINT64_C(0x7777'7777'7777'7777), UINT64_C(0xbbbb'bbbb'bbbb'bbbb)};
                auto diff1 = bigint::uint_var{
                    UINT64_C(0x7777'7777'7777'7777), UINT64_C(0x7777'7777'7777'7777),
                    UINT64_C(0x7777'7777'7777'7777), UINT64_C(0x3333'3333'3333'3333)};
                auto prod = bigint::uint_var{
                    UINT64_C(0x1fdb'9753'0eca'8641), UINT64_C(0xfdb9'7530'eca8'641f),
                    UINT64_C(0xdb97'530e'ca86'41fd), UINT64_C(0xb975'30ec'a864'1fdb),
                    UINT64_C(0x579b'e024'68ac'f135), UINT64_C(0x79be'0246'8acf'1357),
                    UINT64_C(0x9be0'2468'acf1'3579), UINT64_C(0xbe02'468a'cf13'579c)};
                auto prod1 = bigint::uint_var{
                    UINT64_C(0x1fdb'9753'0eca'8641), UINT64_C(0xdddd'dddd'dddd'dddd),
                    UINT64_C(0xdddd'dddd'dddd'dddd), UINT64_C(0xdddd'dddd'dddd'dddd),
                    UINT64_C(0xbe02'468a'cf13'579c)};
                auto quot1 = bigint::uint_var{1, UINT64_C(0xc000'0000'0000'0001),
                                              UINT64_C(0xc000'0000'0000'0001),
                                              UINT64_C(0xc000'0000'0000'0001)};

                expect(x + y == sum);
                expect(bigint::uint_var(x) + y == sum);
                expect(x + bigint::uint_var(y) == sum);
                expect(bigint::uint_var(x) + bigint::uint_var(y) == sum);
                expect(x + UINT64_C(0x4444'4444'4444'4444) == sum1);
                expect(bigint::uint_var(x) + UINT64_C(0x4444'4444'4444'4444) == sum1);
                expect(UINT64_C(0x4444'4444'4444'4444) + x == sum1);
                expect(UINT64_C(0x4444'4444'4444'4444) + bigint::uint_var(x) == sum1);

                expect(x - y == diff);
                expect(bigint::uint_var(x) - y == diff);
                expect(x - bigint::uint_var(y) == diff);
                expect(bigint::uint_var(x) - bigint::uint_var(y) == diff);
                expect(x - UINT64_C(0x4444'4444'4444'4444) == diff1);
                expect(bigint::uint_var(x) - UINT64_C(0x4444'4444'4444'4444) == diff1);

                expect(x * y == prod);
                expect(bigint::uint_var(x) * y == prod);
                expect(x * bigint::uint_var(y) == prod);
                expect(bigint::uint_var(x) * bigint::uint_var(y) == prod);
                expect(x * UINT64_C(0x4444'4444'4444'4444) == prod1);
                expect(bigint::uint_var(x) * UINT64_C(0x4444'4444'4444'4444) == prod1);
                expect(UINT64_C(0x4444'4444'4444'4444) * x == prod1);
                expect(UINT64_C(0x4444'4444'4444'4444) * bigint::uint_var(x) == prod1);

                expect(x / y == 1u);
                expect(bigint::uint_var(x) / y == 1u);
                expect(x / bigint::uint_var(y) == 1u);
                expect(bigint::uint_var(x) / bigint::uint_var(y) == 1u);
                expect(x / UINT64_C(0x4444'4444'4444'4444) == quot1);
                expect(bigint::uint_var(x) / UINT64_C(0x4444'4444'4444'4444) == quot1);

                expect(x % y == diff);
                expect(bigint::uint_var(x) % y == diff);
                expect(x % bigint::uint_var(y) == diff);
                expect(bigint::uint_var(x) % bigint::uint_var(y) == diff);
                expect(x % UINT64_C(0x4444'4444'4444'4444) == UINT64_C(0x3333'3333'3333'3333));
                expect(bigint::uint_var(x) % UINT64_C(0x4444'4444'4444'4444) ==
                       UINT64_C(0x3333'3333'3333'3333));
            }
            {
                auto x = bigint::uint_var{UINT64_C(0x0010'dead'c0de'beef), 0, 0, 0, 0, 0};
                auto y = bigint::uint_var{UINT64_C(0x0000'0060'0000'0000),
                                          UINT64_C(0xdead'c0de'dead'c0de),
                                          UINT64_C(0xdead'c0de'dead'c0de)};
                auto z = bigint::uint_var{UINT64_C(0x0000'0043'7ab7'037a),
                                          UINT64_C(0xfbbc'0000'0000'0000), 0};

                expect(bigint::trunc_floor_log2_div(x, y) == 205);
                expect(bigint::trunc_ceil_log2_div(x, y) == 206);
                expect(bigint::trunc_floor_log2_div(x, z) == 206);
                expect(bigint::trunc_ceil_log2_div(x, z) == 206);
            }
        };
    };
}
