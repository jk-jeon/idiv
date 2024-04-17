#include "idiv/idiv.h"
#include "idiv/rational_continued_fraction.h"
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

int main() {
    using convergent_t =
        jkj::cntfrc::projective_rational<jkj::bigint::int_var, jkj::bigint::uint_var>;
    convergent_t x{7, 18u};
    std::cout << "      Number = " << x.numerator << " / " << x.denominator;

    auto cf = jkj::cntfrc::make_generator<jkj::cntfrc::index_tracker,
                                          jkj::cntfrc::previous_previous_convergent_tracker>(
        jkj::cntfrc::impl::rational{x});
    jkj::bigint::uint_var nmax = 0xffff'ffff;
    std::cout << "\n       n_max = " << nmax;

    auto info1 = jkj::idiv::find_optimal_multiply_shift(cf, nmax);

    std::cout << "\n\n[Multiply-and-shift method]\n"
              << "  Multiplier = " << info1.multiplier << "\n       Shift = " << info1.shift_amount
              << "\n\n";

#if 0				
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
#endif // 0
}
