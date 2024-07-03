#include <idiv/idiv.h>
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

struct multiply_add_shift_info {
    bool succeeded = false;
    jkj::bigint::uint_var multiplier = {};
    jkj::bigint::uint_var adder = {};
    unsigned int shift_amount = 0;
};

constexpr inline multiply_add_shift_info convert_to_multiply_add_shift_effectively_rational(
    jkj::cntfrc::projective_rational<jkj::bigint::int_var, jkj::bigint::uint_var> const& x,
    jkj::bigint::uint_var const& nmax, jkj::bigint::uint_var const& max_allowed_value) {
    jkj::util::constexpr_assert<jkj::util::error_msgs::divide_by_zero>(!x.denominator.is_zero());
    jkj::util::constexpr_assert(x.denominator <= nmax);

    auto internal_continued_fractions_calculator =
        jkj::cntfrc::make_generator<jkj::cntfrc::partial_fraction_tracker,
                                    jkj::cntfrc::convergent_tracker, jkj::cntfrc::interval_tracker>(
            jkj::cntfrc::impl::rational<jkj::bigint::int_var, jkj::bigint::uint_var>{x});

    multiply_add_shift_info ret_value;
    auto continued_fractions_calculator =
        jkj::cntfrc::caching_generator<decltype(internal_continued_fractions_calculator)>{
            std::move(internal_continued_fractions_calculator)};

    jkj::bigint::uint_var n_L0, n_U0;
    if (x.denominator != 1u) {
        // Find the largest multiple of x.denominator <= nmax.
        n_L0 = (nmax / x.denominator) * x.denominator;

        // Compute the modular inverse of -x.numerator.
        auto const mod_inv =
            jkj::idiv::find_best_rational_approx(continued_fractions_calculator, x.denominator - 1u)
                .above.denominator;

        // v = floor((nmax - b) / q) * q + b.
        n_U0 = ((nmax - mod_inv) / x.denominator) * x.denominator;
        n_U0 += mod_inv;
    }
    else {
        n_L0 = nmax;
        n_U0 = nmax;
    }

    using ufrac = jkj::frac<jkj::bigint::uint_var, jkj::bigint::uint_var>;

    ufrac zeta_max{0u, 1u}, zeta_Lmax{0u, 1u}, zeta_Umax{0u, 1u}, zeta_min;
    jkj::bigint::uint_var n_L1 = 0, n_U1 = 0;
    jkj::bigint::uint_var floor_n_L0_x, floor_n_U0_x_p1;

    while (zeta_max.numerator != zeta_max.denominator) {
        multiply_add_shift_info candidate;

        // Update zeta_Lmax if necessary.
        if (zeta_max == zeta_Lmax) {
            n_L0 += n_L1;
            floor_n_L0_x = (n_L0 * x.numerator.abs()) / x.denominator;

            if (n_L0 == nmax) {
                zeta_Lmax = ufrac{1u, 1u};
            }
            else {
                auto const new_nmax = nmax - n_L0;
                continued_fractions_calculator.rewind();
                auto const best_approx =
                    jkj::idiv::find_best_rational_approx(continued_fractions_calculator, new_nmax)
                        .below;
                auto const largest_multiplier = new_nmax / best_approx.denominator;
                n_L1 = largest_multiplier * best_approx.denominator;

                zeta_Lmax = ufrac{(n_L1 * floor_n_L0_x) -
                                      n_L0 * (largest_multiplier * best_approx.numerator.abs()),
                                  n_L1};

                // Truncate to 1 if necessary.
                if (zeta_Lmax.numerator > zeta_Lmax.denominator) {
                    zeta_Lmax = ufrac{1u, 1u};
                }
            }
        }

        // Update zeta_Umax if necessary.
        if (zeta_max == zeta_Umax) {
            n_U0 -= n_U1;
            floor_n_U0_x_p1 = ((n_U0 * x.numerator.abs()) / x.denominator) + 1u;

            if (n_U0 == 1) {
                zeta_Umax = ufrac{1u, 1u};
            }
            else {
                auto const new_nmax = n_U0 - 1u;
                continued_fractions_calculator.rewind();
                auto const best_approx =
                    jkj::idiv::find_best_rational_approx(continued_fractions_calculator, new_nmax)
                        .below;
                auto const largest_multiplier = new_nmax / best_approx.denominator;
                n_U1 = largest_multiplier * best_approx.denominator;

                zeta_Umax = ufrac{(n_U1 * floor_n_U0_x_p1) -
                                      n_U0 * (largest_multiplier * best_approx.numerator.abs()),
                                  n_U1};

                // Truncate to 1 if necessary.
                if (zeta_Umax.numerator > zeta_Umax.denominator) {
                    zeta_Umax = ufrac{1u, 1u};
                }
            }
        }

        zeta_min = static_cast<ufrac&&>(zeta_max);
        zeta_max = zeta_Lmax < zeta_Umax ? zeta_Lmax : zeta_Umax;

        auto const left_end = ufrac{floor_n_L0_x * zeta_max.denominator - zeta_max.numerator,
                                    n_L0 * zeta_max.denominator};
        auto const right_end = ufrac{floor_n_U0_x_p1 * zeta_min.denominator - zeta_min.numerator,
                                     n_U0 * zeta_min.denominator};

        {
            auto delta = ufrac{right_end.numerator * left_end.denominator,
                               left_end.denominator * right_end.denominator};
            auto numerator_diff = right_end.denominator * left_end.numerator;

            // If the interval is empty, move to the next subinterval for zeta.
            if (delta.numerator <= numerator_diff) {
                continue;
            }

            delta.numerator -= numerator_diff;
            jkj::util::constexpr_assert(delta.denominator >= delta.numerator);
            candidate.shift_amount = trunc_floor_log2_div(delta.denominator, delta.numerator);
        }

        candidate.multiplier =
            ((left_end.numerator << candidate.shift_amount) / left_end.denominator) + 1u;

        // If t goes out of the interval, then increase k0.
        if (candidate.multiplier * right_end.denominator >=
            (right_end.numerator << candidate.shift_amount)) {
            ++candidate.shift_amount;
            candidate.multiplier =
                ((left_end.numerator << candidate.shift_amount) / left_end.denominator) + 1u;
        }
        else {
            candidate.shift_amount -= candidate.multiplier.factor_out_power_of_2();
        }

        // Truncate zeta0 from 0 to avoid underflow.
        auto zeta0 = ufrac{(floor_n_L0_x << candidate.shift_amount),
                           jkj::bigint::uint_var::power_of_2(candidate.shift_amount)};
        {
            auto numerator_diff = n_L0 * candidate.multiplier;
            if (zeta0.numerator > numerator_diff) {
                zeta0.numerator -= numerator_diff;
            }
            else {
                zeta0.numerator = 0;
            }
        }

        while (candidate.multiplier * nmax + zeta0.numerator <= max_allowed_value) {
            // Loop over all numerators in the interval.
            while (true) {
                // These branches do not touch candidate.multiplier and
                // candidate.shift_amount, but may modify candidate.adder.
                if (zeta0 >= zeta_min) {
                    candidate.adder = zeta0.numerator;

                    // Check admissibility of xi with respect to zeta.
                    if (candidate.multiplier * n_U0 + candidate.adder <
                        (floor_n_U0_x_p1 << candidate.shift_amount)) {
                        // Found.
                        candidate.succeeded = true;
                        break;
                    }
                }
                else {
                    auto const delta_zeta = zeta_min - zeta0;
                    candidate.adder =
                        zeta0.numerator + div_ceil((delta_zeta.numerator << candidate.shift_amount),
                                                   delta_zeta.denominator);

                    // Check zeta < zeta_max.
                    if (candidate.adder * zeta_max.denominator <
                        (zeta_max.numerator << candidate.shift_amount)) {
                        // Check the max_allowed_value constraint.
                        if (candidate.multiplier * nmax + candidate.adder <= max_allowed_value) {
                            // Check admissibility of xi with respect to zeta.
                            if (candidate.multiplier * n_U0 + candidate.adder <
                                (floor_n_U0_x_p1 << candidate.shift_amount)) {
                                // Found.
                                candidate.succeeded = true;
                                break;
                            }
                        }
                    }
                }

                // Try the next numerator.
                ++candidate.multiplier;
                if (candidate.multiplier * right_end.denominator >=
                    (right_end.numerator << candidate.shift_amount)) {
                    break;
                }

                if (zeta0.numerator > n_L0) {
                    zeta0.numerator -= n_L0;
                }
                else {
                    zeta0.numerator = 0;
                }

                if (candidate.multiplier * nmax + zeta0.numerator > max_allowed_value) {
                    break;
                }
            }

            if (candidate.succeeded) {
                // If this is the first success, record it and move to the next
                // subinterval.
                if (!ret_value.succeeded) {
                    ret_value = candidate;
                }
                // Otherwise, compare it with the previous best one and replace if
                // appropriate.
                else {
                    if (ret_value.shift_amount > candidate.shift_amount) {
                        ret_value = candidate;
                    }
                    else if (ret_value.shift_amount == candidate.shift_amount) {
                        if (ret_value.multiplier > candidate.multiplier) {
                            ret_value = candidate;
                        }
                        else if (ret_value.multiplier == candidate.multiplier) {
                            if (ret_value.adder > candidate.adder) {
                                ret_value = candidate;
                            }
                        }
                    }
                }

                break;
            }

            // Increase k0 and recompute t, zeta0.
            ++candidate.shift_amount;
            candidate.multiplier =
                ((left_end.numerator << candidate.shift_amount) / left_end.denominator) + 1u;

            zeta0.numerator = (floor_n_L0_x << candidate.shift_amount);
            zeta0.denominator <<= 1;
            {
                auto numerator_diff = n_L0 * candidate.multiplier;
                if (zeta0.numerator > numerator_diff) {
                    zeta0.numerator -= numerator_diff;
                }
                else {
                    zeta0.numerator = 0;
                }
            }
        }
    }

    return ret_value;
}

int main() {
    using convergent_t =
        jkj::cntfrc::projective_rational<jkj::bigint::int_var, jkj::bigint::uint_var>;
    convergent_t x{3, 7u};
    std::cout << "      Number = " << x.numerator << " / " << x.denominator;

    auto cf = jkj::cntfrc::make_generator<jkj::cntfrc::index_tracker,
                                          jkj::cntfrc::previous_previous_convergent_tracker>(
        jkj::cntfrc::impl::rational{x});
    jkj::bigint::uint_var nmax = 0xffff'ffff;
    std::cout << "\n       n_max = " << nmax;

    using positive_rational_t = jkj::frac<jkj::bigint::uint_var, jkj::bigint::uint_var>;

    auto info1 = jkj::idiv::find_optimal_multiply_shift(cf, nmax);

    std::cout << "\n\n[Multiply-and-shift method]\n"
              << "  Multiplier = " << info1.multiplier << "\n       Shift = " << info1.shift_amount
              << "\n\n";

    jkj::bigint::uint_var max_allowed = jkj::bigint::uint_var{0xffff'ffff'ffff'ffff};
    auto info2 = convert_to_multiply_add_shift_effectively_rational(x, nmax, max_allowed);

    std::cout << "\n\n[Multiply-add-and-shift method (with max_allowed = " << max_allowed << ")]\n";
    if (info2.succeeded) {
        std::cout << "  Multiplier = " << info2.multiplier << "\n       Adder = " << info2.adder
                  << "\n       Shift = " << info2.shift_amount << "\n";
    }
    else {
        std::cout << "Failed to find any solution.\n";
    }
}
