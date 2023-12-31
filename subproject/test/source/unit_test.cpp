#include <idiv/idiv.h>
#include <idiv/log_continued_fraction.h>
#include <idiv/best_rational_approx.h>
#include <boost/ut.hpp>
#include <format>
#include <iostream>

namespace jkj {
    namespace bigint {
        std::ostream& operator<<(std::ostream& out, jkj::bigint::uint_var const& n) {
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

        std::ostream& operator<<(std::ostream& out, jkj::bigint::int_var const& n) {
            if (n.is_strictly_negative()) {
                out << "-";
            }
            out << n.abs();
            return out;
        }
    }
}

int main() {
    using namespace boost::ut;
    using namespace jkj;

    "[integer_util]"_test = [] {
        auto to_signed_from_unsigned = util::to_signed(3u);
        auto to_signed_from_signed = util::to_signed(3);
        auto to_negative_from_unsigned = util::to_negative(2147483647u);
        auto abs_from_unsigned = util::abs(4294967295u);
        auto abs_from_signed = util::abs(static_cast<int>(-2147483648ll));
        auto invert_sign = util::invert_sign(2147483647);
        expect(std::is_same_v<decltype(to_signed_from_unsigned), int> &&
               to_signed_from_unsigned == 3);
        expect(std::is_same_v<decltype(to_signed_from_signed), int> && to_signed_from_signed == 3);
        expect(std::is_same_v<decltype(to_negative_from_unsigned), int> &&
               to_negative_from_unsigned == -2147483647);
        expect(std::is_same_v<decltype(abs_from_unsigned), unsigned int> &&
               abs_from_unsigned == 4294967295u);
        expect(std::is_same_v<decltype(abs_from_signed), unsigned int> &&
               abs_from_signed == 2147483648u);
        expect(std::is_same_v<decltype(invert_sign), int> && invert_sign == -2147483647);

        auto div_uint_uint = util::div(18u, 7u);
        auto div_int_uint_positive = util::div(37, 9u);
        auto div_int_uint_positive_divisible = util::div(36, 9u);
        auto div_int_uint_negative = util::div(-37, 9u);
        auto div_int_uint_negative_divisible = util::div(-36, 9u);
        expect(std::is_same_v<decltype(div_uint_uint.quot), unsigned int> &&
               div_uint_uint.quot == 2u &&
               std::is_same_v<decltype(div_uint_uint.rem), unsigned int> &&
               div_uint_uint.rem == 4u);
        expect(std::is_same_v<decltype(div_int_uint_positive.quot), int> &&
               div_int_uint_positive.quot == 4 &&
               std::is_same_v<decltype(div_int_uint_positive.rem), unsigned int> &&
               div_int_uint_positive.rem == 1u);
        expect(std::is_same_v<decltype(div_int_uint_positive_divisible.quot), int> &&
               div_int_uint_positive_divisible.quot == 4 &&
               std::is_same_v<decltype(div_int_uint_positive_divisible.rem), unsigned int> &&
               div_int_uint_positive_divisible.rem == 0u);
        expect(std::is_same_v<decltype(div_int_uint_negative.quot), int> &&
               div_int_uint_negative.quot == -5 &&
               std::is_same_v<decltype(div_int_uint_negative.rem), unsigned int> &&
               div_int_uint_negative.rem == 8u);
        expect(std::is_same_v<decltype(div_int_uint_negative_divisible.quot), int> &&
               div_int_uint_negative_divisible.quot == -4 &&
               std::is_same_v<decltype(div_int_uint_negative_divisible.rem), unsigned int> &&
               div_int_uint_negative_divisible.rem == 0u);

        auto div_floor_uint_uint = util::div_floor(71u, 9u);
        auto div_floor_uint_uint_divisible = util::div_floor(72u, 9u);
        auto div_floor_int_uint_positive = util::div_floor(71, 9u);
        auto div_floor_int_uint_positive_divisible = util::div_floor(72, 9u);
        auto div_floor_int_uint_negative = util::div_floor(-71, 9u);
        auto div_floor_int_uint_negative_divisible = util::div_floor(-72, 9u);
        auto div_floor_uint_int_positive = util::div_floor(71u, 9);
        auto div_floor_uint_int_positive_divisible = util::div_floor(72u, 9);
        auto div_floor_uint_int_negative = util::div_floor(71u, -9);
        auto div_floor_uint_int_negative_divisible = util::div_floor(72u, -9);
        auto div_floor_int_int_positive = util::div_floor(71, 9);
        auto div_floor_int_int_positive_divisible = util::div_floor(72, 9);
        auto div_floor_int_int_negative = util::div_floor(71, -9);
        auto div_floor_int_int_negative_divisible = util::div_floor(72, -9);
        expect(std::is_same_v<decltype(div_floor_uint_uint), unsigned int> &&
               div_floor_uint_uint == 7u);
        expect(std::is_same_v<decltype(div_floor_uint_uint_divisible), unsigned int> &&
               div_floor_uint_uint_divisible == 8u);
        expect(std::is_same_v<decltype(div_floor_int_uint_positive), int> &&
               div_floor_int_uint_positive == 7);
        expect(std::is_same_v<decltype(div_floor_int_uint_positive_divisible), int> &&
               div_floor_int_uint_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_floor_int_uint_negative), int> &&
               div_floor_int_uint_negative == -8);
        expect(std::is_same_v<decltype(div_floor_int_uint_negative_divisible), int> &&
               div_floor_int_uint_negative_divisible == -8);
        expect(std::is_same_v<decltype(div_floor_uint_int_positive), int> &&
               div_floor_uint_int_positive == 7);
        expect(std::is_same_v<decltype(div_floor_uint_int_positive_divisible), int> &&
               div_floor_uint_int_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_floor_uint_int_negative), int> &&
               div_floor_uint_int_negative == -8);
        expect(std::is_same_v<decltype(div_floor_uint_int_negative_divisible), int> &&
               div_floor_uint_int_negative_divisible == -8);
        expect(std::is_same_v<decltype(div_floor_int_int_positive), int> &&
               div_floor_int_int_positive == 7);
        expect(std::is_same_v<decltype(div_floor_int_int_positive_divisible), int> &&
               div_floor_int_int_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_floor_int_int_negative), int> &&
               div_floor_int_int_negative == -8);
        expect(std::is_same_v<decltype(div_floor_int_int_negative_divisible), int> &&
               div_floor_int_int_negative_divisible == -8);

        auto div_ceil_uint_uint = util::div_ceil(71u, 9u);
        auto div_ceil_uint_uint_divisible = util::div_ceil(72u, 9u);
        auto div_ceil_int_uint_positive = util::div_ceil(71, 9u);
        auto div_ceil_int_uint_positive_divisible = util::div_ceil(72, 9u);
        auto div_ceil_int_uint_negative = util::div_ceil(-71, 9u);
        auto div_ceil_int_uint_negative_divisible = util::div_ceil(-72, 9u);
        auto div_ceil_uint_int_positive = util::div_ceil(71u, 9);
        auto div_ceil_uint_int_positive_divisible = util::div_ceil(72u, 9);
        auto div_ceil_uint_int_negative = util::div_ceil(71u, -9);
        auto div_ceil_uint_int_negative_divisible = util::div_ceil(72u, -9);
        auto div_ceil_int_int_positive = util::div_ceil(71, 9);
        auto div_ceil_int_int_positive_divisible = util::div_ceil(72, 9);
        auto div_ceil_int_int_negative = util::div_ceil(71, -9);
        auto div_ceil_int_int_negative_divisible = util::div_ceil(72, -9);
        expect(std::is_same_v<decltype(div_ceil_uint_uint), unsigned int> &&
               div_ceil_uint_uint == 8u);
        expect(std::is_same_v<decltype(div_ceil_uint_uint_divisible), unsigned int> &&
               div_ceil_uint_uint_divisible == 8u);
        expect(std::is_same_v<decltype(div_ceil_int_uint_positive), int> &&
               div_ceil_int_uint_positive == 8);
        expect(std::is_same_v<decltype(div_ceil_int_uint_positive_divisible), int> &&
               div_ceil_int_uint_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_ceil_int_uint_negative), int> &&
               div_ceil_int_uint_negative == -7);
        expect(std::is_same_v<decltype(div_ceil_int_uint_negative_divisible), int> &&
               div_ceil_int_uint_negative_divisible == -8);
        expect(std::is_same_v<decltype(div_ceil_uint_int_positive), int> &&
               div_ceil_uint_int_positive == 8);
        expect(std::is_same_v<decltype(div_ceil_uint_int_positive_divisible), int> &&
               div_ceil_uint_int_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_ceil_uint_int_negative), int> &&
               div_ceil_uint_int_negative == -7);
        expect(std::is_same_v<decltype(div_ceil_uint_int_negative_divisible), int> &&
               div_ceil_uint_int_negative_divisible == -8);
        expect(std::is_same_v<decltype(div_ceil_int_int_positive), int> &&
               div_ceil_int_int_positive == 8);
        expect(std::is_same_v<decltype(div_ceil_int_int_positive_divisible), int> &&
               div_ceil_int_int_positive_divisible == 8);
        expect(std::is_same_v<decltype(div_ceil_int_int_negative), int> &&
               div_ceil_int_int_negative == -7);
        expect(std::is_same_v<decltype(div_ceil_int_int_negative_divisible), int> &&
               div_ceil_int_int_negative_divisible == -8);
    };

    "[tmp]"_test = [] {
        using typelist = tmp::typelist<int, double, float, char, float, double, unsigned int, long>;
        should("find_first_index") = [] {
            expect(tmp::find_first_index<int>(typelist{}) == 0);
            expect(tmp::find_first_index<double>(typelist{}) == 1);
            expect(tmp::find_first_index<float>(typelist{}) == 2);
            expect(tmp::find_first_index<char>(typelist{}) == 3);
            expect(tmp::find_first_index<unsigned int>(typelist{}) == 6);
            expect(tmp::find_first_index<long>(typelist{}) == 7);
            expect(tmp::find_first_index<short>(typelist{}) == 8);
        };

        should("get_type") = [] {
            expect(std::is_same_v<tmp::get_type<0, typelist>, int>);
            expect(std::is_same_v<tmp::get_type<1, typelist>, double>);
            expect(std::is_same_v<tmp::get_type<2, typelist>, float>);
            expect(std::is_same_v<tmp::get_type<3, typelist>, char>);
            expect(std::is_same_v<tmp::get_type<4, typelist>, float>);
            expect(std::is_same_v<tmp::get_type<5, typelist>, double>);
            expect(std::is_same_v<tmp::get_type<6, typelist>, unsigned int>);
            expect(std::is_same_v<tmp::get_type<7, typelist>, long>);
        };

        should("back_sublist") = [] {
            expect(std::is_same_v<tmp::back_sublist<0, typelist>, tmp::typelist<>>);
            expect(std::is_same_v<tmp::back_sublist<1, typelist>, tmp::typelist<long>>);
            expect(
                std::is_same_v<tmp::back_sublist<2, typelist>, tmp::typelist<unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<3, typelist>,
                                  tmp::typelist<double, unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<4, typelist>,
                                  tmp::typelist<float, double, unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<5, typelist>,
                                  tmp::typelist<char, float, double, unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<6, typelist>,
                                  tmp::typelist<float, char, float, double, unsigned int, long>>);
            expect(std::is_same_v<
                   tmp::back_sublist<7, typelist>,
                   tmp::typelist<double, float, char, float, double, unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<8, typelist>, typelist>);
        };

        should("remove_duplicate") = [] {
            expect(std::is_same_v<tmp::remove_duplicate<typelist>,
                                  tmp::typelist<int, double, float, char, unsigned int, long>>);
        };

        should("push_back") = [] {
            expect(std::is_same_v<tmp::push_back<typelist, unsigned short>,
                                  tmp::typelist<int, double, float, char, float, double,
                                                unsigned int, long, unsigned short>>);
        };

        should("join") = [] {
            expect(std::is_same_v<
                   tmp::join<typelist, tmp::typelist<unsigned short, unsigned short>, typelist>,
                   tmp::typelist<int, double, float, char, float, double, unsigned int, long,
                                 unsigned short, unsigned short, int, double, float, char, float,
                                 double, unsigned int, long>>);
        };

        should("filter") = [] {
            auto predicate = [](auto arg) {
                return std::is_integral_v<typename decltype(arg)::type>;
            };

            expect(std::is_same_v<tmp::filter<typelist, decltype(predicate)>,
                                  tmp::typelist<int, char, unsigned int, long>>);
        };
    };

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

    "[topological_sort]"_test = [] {
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

    "[get_transitive_required_mixin_list]"_test = [] {
        auto list = cntfrc::detail::get_transitive_required_mixin_list(
            tmp::typelist<cntfrc::detail::mixin_type_wrapper<cntfrc::interval_tracker>,
                          cntfrc::detail::mixin_type_wrapper<cntfrc::convergent_tracker>>{});

        expect(std::is_same_v<
               decltype(list),
               tmp::typelist<cntfrc::detail::mixin_type_wrapper<cntfrc::interval_tracker>,
                             cntfrc::detail::mixin_type_wrapper<cntfrc::convergent_tracker>,
                             cntfrc::detail::mixin_type_wrapper<cntfrc::index_tracker>>>);
    };

    "[find_sorted_mixin_list]"_test = [] {
        struct dummy_type {
            using required_mixins =
                cntfrc::mixin_list<cntfrc::index_tracker, cntfrc::partial_fraction_tracker,
                                   cntfrc::convergent_tracker>;
            using mixin_ordering_constraints = cntfrc::mixin_ordering_constraint::constraint_list<
                cntfrc::mixin_ordering_constraint::before_after<cntfrc::index_tracker,
                                                                cntfrc::interval_tracker>,
                cntfrc::mixin_ordering_constraint::before_after<cntfrc::convergent_tracker,
                                                                cntfrc::interval_tracker>>;
        };
        auto sorted_wrapped_mixin_list =
            cntfrc::detail::find_sorted_mixin_list<dummy_type, cntfrc::interval_tracker>();
        expect(std::is_same_v<
               decltype(sorted_wrapped_mixin_list),
               tmp::typelist<cntfrc::detail::mixin_type_wrapper<cntfrc::convergent_tracker>,
                             cntfrc::detail::mixin_type_wrapper<cntfrc::partial_fraction_tracker>,
                             cntfrc::detail::mixin_type_wrapper<cntfrc::index_tracker>,
                             cntfrc::detail::mixin_type_wrapper<cntfrc::interval_tracker>>>);
    };

    "[rational_continued_fraction]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        {
            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                cntfrc::impl::rational{convergent_t{156, 179u}});

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
            auto cf = cntfrc::make_generator<cntfrc::convergent_tracker>(
                cntfrc::impl::rational{convergent_t{-2767, 1982u}});

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
        using unary_gosper_t =
            cntfrc::impl::unary_gosper<cntfrc::impl::rational<bigint::int_var, bigint::uint_var>,
                                       cntfrc::unity>;

        // 481/2245 = (-18*156 + 13*179)/(12*156 - 23*179)
        auto cf1 = cntfrc::make_generator<cntfrc::convergent_tracker>(unary_gosper_t{
            unary_gosper_t::internal_continued_fraction_impl_type{convergent_t{156, 179u}},
            {-18, 13, 12, -23}});
        auto cf2 = cntfrc::make_generator<cntfrc::convergent_tracker>(
            cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{convergent_t{481, 2245u}});

        while (!cf2.terminated()) {
            expect(cf1.update() == cf2.update());
            expect(cf1.current_convergent() == cf2.current_convergent());
        }
        expect(cf1.update() == cf2.update());
    };

    "[binary_gosper]"_test = [] {
        using convergent_t = cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using rational_continued_fraction_t =
            cntfrc::impl::rational<bigint::int_var, bigint::uint_var>;
        using binary_gosper_t = cntfrc::impl::binary_gosper<rational_continued_fraction_t,
                                                            rational_continued_fraction_t>;

        // Take x = 17/89, y = 31/125, and
        // z = (xy + 4x + 2y + 8)/(2x - 3y + 1) = 2655/182.
        auto cf1 = cntfrc::make_generator<cntfrc::convergent_tracker>(
            binary_gosper_t{rational_continued_fraction_t{convergent_t{17, 89u}},
                            rational_continued_fraction_t{convergent_t{31, 125u}},
                            {// numerator
                             1, 4, 2, 8,
                             // denominator
                             0, 2, -3, 1}});
        auto cf2 = cntfrc::make_generator<cntfrc::convergent_tracker>(
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
            // Rational case.
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

    "[best_rational_approx]"_test = [] {
        using projective_rational_t =
            cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using rational_t = frac<bigint::int_var, bigint::uint_var>;
        auto perform_test = [](bigint::int_var const& numerator,
                               bigint::uint_var const& denominator, std::size_t nmax) {
            auto cf = cntfrc::make_generator<cntfrc::index_tracker,
                                             cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{
                    projective_rational_t{numerator, denominator}});

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

    "[find_floor_quotient_range]"_test = [] {
        using projective_rational_t =
            cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using rational_t = frac<bigint::int_var, bigint::uint_var>;
        auto perform_test = [](bigint::int_var const& numerator,
                               bigint::uint_var const& denominator, std::size_t nmax) {
            auto cf = cntfrc::make_generator<cntfrc::index_tracker,
                                             cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{
                    projective_rational_t{numerator, denominator}});

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

    "[find_extrema_of_fractional_part]"_test = [] {
        using projective_rational_t =
            cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        auto perform_test = [](bigint::int_var const& numerator,
                               bigint::uint_var const& denominator, std::size_t nmax) {
            auto cf = cntfrc::make_generator<cntfrc::index_tracker,
                                             cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{
                    projective_rational_t{numerator, denominator}});

            auto result = idiv::find_extrema_of_fractional_part(cf, nmax);
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

    "[find_optimal_multiply_shift]"_test = [] {
        using projective_rational_t =
            cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        auto perform_test = [](bigint::int_var const& numerator,
                               bigint::uint_var const& denominator, std::size_t nmax) {
            auto cf = cntfrc::make_generator<cntfrc::index_tracker,
                                             cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{
                    projective_rational_t{numerator, denominator}});

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

    "[find_suboptimal_multiply_add_shift]"_test = [] {
        using projective_rational_t =
            cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using nrange_t = interval<bigint::int_var, interval_type_t::bounded_closed>;
        auto perform_test = [](projective_rational_t const& x, projective_rational_t const& y,
                               nrange_t const& nrange) {
            auto xcf = cntfrc::make_generator<cntfrc::index_tracker,
                                              cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{x});
            auto ycf = cntfrc::make_generator<cntfrc::index_tracker,
                                              cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{y});

            auto result = idiv::find_suboptimal_multiply_add_shift(xcf, ycf, nrange);

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

    "[find_extrema_of_fractional_part (two unknowns)]"_test = [] {
        using projective_rational_t =
            cntfrc::projective_rational<bigint::int_var, bigint::uint_var>;
        using nrange_t = interval<bigint::int_var, interval_type_t::bounded_closed>;
        auto perform_test = [](projective_rational_t const& x, projective_rational_t const& y,
                               nrange_t const& nrange) {
            auto xcf = cntfrc::make_generator<cntfrc::index_tracker,
                                              cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{x});
            auto ycf = cntfrc::make_generator<cntfrc::index_tracker,
                                              cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational<bigint::int_var, bigint::uint_var>{y});

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
}
