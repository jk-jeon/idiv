#include <idiv/bigint.h>
#include <idiv/gosper_continued_fraction.h>
#include <idiv/log_continued_fraction.h>
#include <idiv/rational_continued_fraction.h>
#include <boost/ut.hpp>

int main() {
    using namespace boost::ut;
    using namespace jkj;

    "[wuint]"_test = [] {
        should("add_carry64") = [] {
            unsigned int carry = 0;
            expect(wuint::add_carry64(UINT64_C(0xffff'ffff'ffff'fff0), 0x10, carry) == 0);
            expect(wuint::add_carry64(UINT64_C(0xffff'ffff'ffff'fff0), 0x11, carry) == 2);
            expect(wuint::add_carry64(0, 10, carry) == 11);
            expect(wuint::add_carry64(0, 10, carry) == 10);
        };

        should("sub_borrow64") = [] {
            unsigned int borrow = 0;
            expect(wuint::sub_borrow64(0, 1, borrow) == UINT64_C(0xffff'ffff'ffff'ffff));
            expect(wuint::sub_borrow64(0, 1, borrow) == UINT64_C(0xffff'ffff'ffff'fffe));
            expect(wuint::sub_borrow64(10, 1, borrow) == 8);
            expect(wuint::sub_borrow64(10, 1, borrow) == 9);
        };

        should("uint128::operator+=") = [] {
            wuint::uint128 x{0xffff'ffff, UINT64_C(0xffff'ffff'ffff'ffff)};
            x += UINT64_C(0xffff'ffff'ffff'ffff);
            expect(x.high() == UINT64_C(0x1'0000'0000) &&
                   x.low() == UINT64_C(0xffff'ffff'ffff'fffe));
        };

        should("umul128 and umul128_upper64") = [] {
            constexpr auto x =
                wuint::umul128(UINT64_C(0x12345678'12345678), UINT64_C(0x87654321'87654321));
            expect(x.high() == 0x9a0cd05'83fa2782 && x.low() == 0xeb11e7f5'70b88d78);
            expect(wuint::umul128_upper64(UINT64_C(0x12345678'12345678),
                                          UINT64_C(0x87654321'87654321)) == 0x9a0cd05'83fa2782);
        };
    };

    "[bigint]"_test = [] {
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

    "[prime_factorization]"_test = [] {
        using uint_type = bigint::uint_var;
        auto prime_factors = prime_factorization(
            frac<uint_type, uint_type>{UINT64_C(185'712'372'396), UINT64_C(198'759'128'733)});

        expect(prime_factors[0] == prime_factor<uint_type>{2u, 2});
        expect(prime_factors[1] == prime_factor<uint_type>{11u, -1});
        expect(prime_factors[2] == prime_factor<uint_type>{37u, 1});
        expect(prime_factors[3] == prime_factor<uint_type>{79u, 1});
        expect(prime_factors[4] == prime_factor<uint_type>{853u, 1});
        expect(prime_factors[5] == prime_factor<uint_type>{2'069u, 1});
        expect(prime_factors[6] == prime_factor<uint_type>{5'861u, -1});
        expect(prime_factors[7] == prime_factor<uint_type>{342'547u, -1});
        expect(prime_factors.size() == 8);
    };

    "[rational_continued_fraction]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        {
            auto cf = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
                cntfrc::rational_continued_fraction{convergent_t{156, 179u}});

            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{0, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{1, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{6, 7u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{7, 8u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{27, 31u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{34, 39u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{61, 70u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{156, 179u});
            expect(cf.update() == false);
        }
        {
            auto cf = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
                cntfrc::rational_continued_fraction{convergent_t{-2767, 1982u}});

            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-2, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-1, 1u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-3, 2u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-4, 3u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-7, 5u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-67, 48u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-74, 53u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-141, 101u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-497, 356u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-1135, 813u});
            expect(cf.update() == true);
            expect(cf.current_convergent() == convergent_t{-2767, 1982u});
            expect(cf.update() == false);
        }
    };

    "[unary_gosper]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using unary_gosper_t = cntfrc::unary_gosper<
            cntfrc::continued_fraction<
                cntfrc::rational_continued_fraction<bigint::int_var, bigint::uint_var>,
                cntfrc::index_tracker, cntfrc::convergent_tracker, cntfrc::interval_tracker>,
            cntfrc::unity>;

        // 481/2245 = (-18*156 + 13*179)/(12*156 - 23*179)
        auto cf1 =
            cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(unary_gosper_t{
                unary_gosper_t::internal_continued_fraction_impl_type{convergent_t{156, 179u}},
                {-18, 13, 12, -23}});
        auto cf2 = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
            cntfrc::rational_continued_fraction<bigint::int_var, bigint::uint_var>{
                convergent_t{481, 2245u}});

        while (!cf2.terminated()) {
            expect(cf1.update() == cf2.update());
            expect(cf1.current_convergent() == cf2.current_convergent());
        }
        expect(cf1.update() == cf2.update());
    };

    "[binary_gosper]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using rational_continued_fraction_t =
            cntfrc::rational_continued_fraction<bigint::int_var, bigint::uint_var>;
        using binary_gosper_t = cntfrc::binary_gosper<
            cntfrc::continued_fraction<rational_continued_fraction_t, cntfrc::index_tracker,
                                       cntfrc::convergent_tracker, cntfrc::interval_tracker>,
            cntfrc::continued_fraction<rational_continued_fraction_t, cntfrc::index_tracker,
                                       cntfrc::convergent_tracker, cntfrc::interval_tracker>>;

        // Take x = 17/89, y = 31/125, and
        // z = (xy + 4x + 2y + 8)/(2x - 3y + 1) = 2655/182.
        auto cf1 = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
            binary_gosper_t{rational_continued_fraction_t{convergent_t{17, 89u}},
                            rational_continued_fraction_t{convergent_t{31, 125u}},
                            {// numerator
                             1, 4, 2, 8,
                             // denominator
                             0, 2, -3, 1}});
        auto cf2 = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
            rational_continued_fraction_t{convergent_t{2655, 182u}});

        while (!cf2.terminated()) {
            expect(cf1.update() == cf2.update());
            expect(cf1.current_convergent() == cf2.current_convergent());
        }
        expect(cf1.update() == cf2.update());
    };

    "[log_continued_fraction]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using unsigned_frac_t = frac<bigint::uint_var, bigint::uint_var>;

        should("natural_log_calculator") = [] {
            unsigned_frac_t error_bound{
                1u, bigint::uint_var(bigint::decimal_uint_const_v<1'000'000, 0, 0, 0, 0, 0>)};

            {
                // Compute ln(2) up to 100 digits.
                auto nlc = cntfrc::make_continued_fraction_generator<
                    cntfrc::index_tracker, cntfrc::partial_fraction_tracker,
                    cntfrc::convergent_tracker, cntfrc::interval_tracker>(
                    cntfrc::natural_log_calculator<bigint::int_var, bigint::uint_var>{
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
            auto cf = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
                cntfrc::natural_log<bigint::int_var, bigint::uint_var>{unsigned_frac_t{3u, 1u}});

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
            // Rational case.
            {
                // log(175616/91125) / log(3136/2025)
                // = log(2^9*7^4 / 3^6*5^3) / log(2^6*7^2 / 3^4*5^2)
                // = 3/2.
                auto cf = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
                    cntfrc::general_log<bigint::int_var, bigint::uint_var>{
                        unsigned_frac_t{3136u, 2025u}, unsigned_frac_t{175'616u, 91'125u}});

                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{1, 1u});
                expect(cf.update() == true);
                expect(cf.current_convergent() == convergent_t{3, 2u});
                expect(cf.update() == false);
            }
            // Irrational case.
            {
                auto cf = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
                    cntfrc::general_log<bigint::int_var, bigint::uint_var>{
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
            using continued_fraction_t = cntfrc::unary_gosper<cntfrc::continued_fraction<
                cntfrc::natural_log_calculator<bigint::int_var, bigint::uint_var>,
                cntfrc::index_tracker, cntfrc::partial_fraction_tracker, cntfrc::convergent_tracker,
                cntfrc::interval_tracker>>;

            auto cf = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
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
                cntfrc::natural_log_calculator<bigint::int_var, bigint::uint_var>;
            using continued_fraction_t = cntfrc::binary_gosper<
                cntfrc::continued_fraction<log_calculator, cntfrc::index_tracker,
                                           cntfrc::partial_fraction_tracker,
                                           cntfrc::convergent_tracker, cntfrc::interval_tracker>,
                cntfrc::continued_fraction<log_calculator, cntfrc::index_tracker,
                                           cntfrc::partial_fraction_tracker,
                                           cntfrc::convergent_tracker, cntfrc::interval_tracker>>;

            auto cf = cntfrc::make_continued_fraction_generator<cntfrc::convergent_tracker>(
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
