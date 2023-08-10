#include "idiv/idiv.h"
#include <iostream>

int main() {
    std::uint64_t numerator = 7;
    std::uint64_t denominator = 18;
    std::uint64_t nmax = 0xffff'ffff;

    auto info1 = jkj::idiv::convert_to_multiply_shift_effectively_rational({numerator, denominator},
                                                                           nmax);

    auto info2 = jkj::idiv::convert_to_multiply_add_shift_effectively_rational(
        {numerator, denominator}, nmax,
        jkj::big_uint::var{0xffff'ffff'ffff'ffff, 0xffff'ffff'ffff'ffff});

    std::cout << "     Number = " << numerator << " / " << denominator
              << "\n\n Multiplier = " << info1.multiplier[0]
              << "\n      Shift = " << info1.shift_amount
              << "\n\n Multiplier = " << info2->multiplier[0]
              << "\n      Adder = " << (info2->adder.is_zero() ? 0 : info2->adder[0])
              << "\n      Shift = " << info2->shift_amount << "\n";
}