#include "idiv/idiv.h"
#include <format>
#include <iostream>

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

int main() {
    using frac = jkj::frac<jkj::bigint::int_var, jkj::bigint::uint_var>;
    frac x{7, 18};
    std::cout << "      Number = " << x.numerator << " / " << x.denominator;

    using continued_fractions_calc_type = jkj::caching_continued_fractions<
        jkj::rational_continued_fractions<jkj::bigint::int_var, jkj::bigint::uint_var>,
        std::vector>;
    continued_fractions_calc_type cf{x};

    jkj::bigint::uint_var nmax = 0xffff'ffff;
    std::cout << "\n       n_max = " << nmax;

    auto info1 = jkj::idiv::convert_to_multiply_shift(cf, nmax);

    std::cout << "\n\n[Multiply-and-shift method]\n"
              << "  Multiplier = " << info1.multiplier << "\n       Shift = " << info1.shift_amount;

    
    jkj::bigint::uint_var max_allowed = jkj::bigint::uint_var{0xffff'ffff'ffff'ffff};
    auto info2 = jkj::idiv::convert_to_multiply_add_shift_effectively_rational(
        {x.numerator.abs(), x.denominator}, nmax, max_allowed);

    std::cout << "\n\n[Multiply-add-and-shift method (with max_allowed = " << max_allowed << ")]\n";
    if (info2.succeeded) {
        std::cout << "  Multiplier = " << info2.multiplier << "\n       Adder = " << info2.adder
                  << "\n       Shift = " << info2.shift_amount << "\n";
    }
    else {
        std::cout << "Failed to find any solution.\n";
    }
    
}
