#include "idiv/idiv.h"
#include <format>
#include <iostream>

std::ostream& operator<<(std::ostream& out, jkj::big_uint::var const& n) {
    auto n_copy = n;
    auto const divisor = jkj::big_uint::var{UINT64_C(1000'0000'0000'0000'0000)};
    std::vector<std::uint64_t> nineteen_digits;
    while (!n_copy.is_zero()) {
        auto q = n_copy.long_division(divisor);
        if (n_copy.is_zero()) {
            nineteen_digits.push_back(0);
        }
        else {
            nineteen_digits.push_back(n_copy[0]);
        }
        n_copy = std::move(q);
    }

    if (nineteen_digits.empty()) {
        out << "0";
    }
    else {
        auto itr = nineteen_digits.crbegin();
        out << std::format("{}", *itr);

        for (++itr; itr != nineteen_digits.crend(); ++itr) {
            out << std::format("{:19}", *itr);
        }
    }

    return out;
}

int main() {
    jkj::big_uint::var numerator{7};
    jkj::big_uint::var denominator{18};
    jkj::big_uint::var nmax = 0xffff'ffff;

    auto info1 =
        jkj::idiv::convert_to_multiply_shift_effectively_rational({numerator, denominator}, nmax);

    auto info2 = jkj::idiv::convert_to_multiply_add_shift_effectively_rational(
        {numerator, denominator}, nmax,
        jkj::big_uint::var{0xffff'ffff'ffff'ffff, 0xffff'ffff'ffff'ffff});

    std::cout << "     Number = " << numerator << " / " << denominator
              << "\n\n Multiplier = " << info1.multiplier
              << "\n      Shift = " << info1.shift_amount
              << "\n\n Multiplier = " << info2->multiplier << "\n      Adder = " << info2->adder
              << "\n      Shift = " << info2->shift_amount << "\n";
}