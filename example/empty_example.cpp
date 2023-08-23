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
    jkj::bigint::decimal_uint_const_t<3, 7, 11> x;
    auto y = to_negative(x) % jkj::bigint::uint_const_v<6, 11>;
    auto q = to_negative(x) / jkj::bigint::uint_const_v<6, 11>;
    jkj::bigint::uint_var yy{y};
    std::cout << yy << "\n" << jkj::bigint::int_var{q};
}