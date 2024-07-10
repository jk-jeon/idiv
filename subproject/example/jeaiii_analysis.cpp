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

// Given a nonzero real number x and a range [nmin:nmax] of positive integers, find the interval I
// such that n = floor(floor(n xi)/x) holds for all n in [nmin:nmax] if and only if xi in I. The
// interval might be empty.
template <class ContinuedFractionGenerator>
jkj::variable_shape_interval<jkj::frac<jkj::bigint::int_var, jkj::bigint::uint_var>,
                             jkj::interval_type_t::empty,
                             jkj::interval_type_t::bounded_left_open_right_closed>
jeaiii_analysis_floor(
    ContinuedFractionGenerator&& xcf,
    jkj::interval<jkj::bigint::int_var, jkj::interval_type_t::bounded_closed> const& nrange) {
    // Note that n = floor(floor(n xi)/x) holds if and only if
    // nx <= floor(n xi) < (n+1)x, if and only if
    // ceil(nx)/n <= xi < ceil((n+1)x)/n, if and only if
    // floor(n(-x) + (-x))/n < -xi <= floor(n(-x))/n.
    // Hence, we find the maximum of the LHS and the minimum of the RHS.

    namespace cntfrc = jkj::cntfrc;

    // To compute ceil(nx) and ceil((n+1)x).
    auto approx_x = jkj::idiv::find_best_rational_approx(
                        cntfrc::make_generator<cntfrc::index_tracker,
                                               cntfrc::previous_previous_convergent_tracker>(
                            xcf.copy_internal_implementation()),
                        nrange.upper_bound() + 1u)
                        .above;

    auto minus_xcf = cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                            cntfrc::interval_tracker>(
        cntfrc::impl::unary_gosper{static_cast<ContinuedFractionGenerator&&>(xcf), {-1, 0, 0, 1}});
    auto zero_cf = cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                          cntfrc::interval_tracker>(
        cntfrc::impl::rational<jkj::bigint::int_var, jkj::bigint::uint_var>{0, 1u});

    auto maximizer = jkj::idiv::find_maximizer_of_floor_subtract_quotient_positive_range(
        minus_xcf.copy(), minus_xcf.copy(), zero_cf.copy(), nrange);
    auto minimizer = jkj::idiv::find_minimizer_of_floor_subtract_quotient_positive_range(
        minus_xcf, zero_cf.copy(), zero_cf, nrange);

    // ceil(nx)/n.
    auto lower_bound =
        jkj::idiv::find_best_rational_approx(
            cntfrc::make_generator<cntfrc::index_tracker,
                                   cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational{
                    jkj::util::div_ceil(minimizer * approx_x.numerator, approx_x.denominator),
                    jkj::util::abs(minimizer)}),
            approx_x.denominator)
            .above;

    // ceil((n+1)x)/n.
    auto upper_bound =
        jkj::idiv::find_best_rational_approx(
            cntfrc::make_generator<cntfrc::index_tracker,
                                   cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational{
                    jkj::util::div_ceil((maximizer + 1) * approx_x.numerator, approx_x.denominator),
                    jkj::util::abs(maximizer)}),
            approx_x.denominator)
            .above;

    if (lower_bound < upper_bound) {
        return jkj::interval<jkj::frac<jkj::bigint::int_var, jkj::bigint::uint_var>,
                             jkj::interval_type_t::bounded_left_open_right_closed>{
            std::move(lower_bound), std::move(upper_bound)};
    }
    else {
        return jkj::interval<jkj::frac<jkj::bigint::int_var, jkj::bigint::uint_var>,
                             jkj::interval_type_t::empty>{};
    }
}

// Given a nonzero real number x and a range [nmin:nmax] of positive integers, find the interval I
// such that n = floor((floor(n xi) + 1)/x) holds for all n in [nmin:nmax] if and only if xi in I.
// The interval might be empty.
template <class ContinuedFractionGenerator>
jkj::variable_shape_interval<jkj::frac<jkj::bigint::int_var, jkj::bigint::uint_var>,
                             jkj::interval_type_t::empty,
                             jkj::interval_type_t::bounded_left_open_right_closed>
jeaiii_analysis_floor_plus_one(
    ContinuedFractionGenerator&& xcf,
    jkj::interval<jkj::bigint::int_var, jkj::interval_type_t::bounded_closed> const& nrange) {
    // Note that n = floor((floor(n xi) + 1)/x) holds if and only if
    // nx <= floor(n xi) + 1 < (n+1)x, if and only if
    // (ceil(nx) - 1)/n <= xi < (ceil((n+1)x) - 1)/n, if and only if
    // (floor(n(-x) + (-x)) + 1)/n < -xi <= (floor(n(-x)) + 1)/n.
    // Hence, we find the maximum of the LHS and the minimum of the RHS.

    namespace cntfrc = jkj::cntfrc;

    // To compute ceil(nx) and ceil((n+1)x).
    auto approx_x = jkj::idiv::find_best_rational_approx(
                        cntfrc::make_generator<cntfrc::index_tracker,
                                               cntfrc::previous_previous_convergent_tracker>(
                            xcf.copy_internal_implementation()),
                        nrange.upper_bound() + 1u)
                        .above;

    auto minus_xcf = cntfrc::make_generator<cntfrc::previous_previous_convergent_tracker,
                                            cntfrc::interval_tracker>(
        cntfrc::impl::unary_gosper{static_cast<ContinuedFractionGenerator&&>(xcf), {-1, 0, 0, 1}});
    auto zero_cf = cntfrc::make_generator<cntfrc::interval_tracker>(
        cntfrc::impl::rational<jkj::bigint::int_var, jkj::bigint::uint_var>{0, 1u});

    auto minus_one_cf =
        cntfrc::make_generator<cntfrc::index_tracker, cntfrc::previous_previous_convergent_tracker>(
            cntfrc::impl::rational<jkj::bigint::int_var, jkj::bigint::uint_var>{-1, 1u});

    auto maximizer = jkj::idiv::find_maximizer_of_floor_subtract_quotient_positive_range(
        minus_xcf.copy(), minus_xcf.copy(), minus_one_cf.copy(), nrange);
    auto minimizer = jkj::idiv::find_minimizer_of_floor_subtract_quotient_positive_range(
        minus_xcf, zero_cf, minus_one_cf, nrange);

    // (ceil(nx) - 1)/n.
    auto lower_bound =
        jkj::idiv::find_best_rational_approx(
            cntfrc::make_generator<cntfrc::index_tracker,
                                   cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational{
                    jkj::util::div_ceil(minimizer * approx_x.numerator, approx_x.denominator) - 1,
                    jkj::util::abs(minimizer)}),
            approx_x.denominator)
            .above;

    // (ceil((n+1)x) - 1)/n.
    auto upper_bound =
        jkj::idiv::find_best_rational_approx(
            cntfrc::make_generator<cntfrc::index_tracker,
                                   cntfrc::previous_previous_convergent_tracker>(
                cntfrc::impl::rational{jkj::util::div_ceil((maximizer + 1) * approx_x.numerator,
                                                           approx_x.denominator) -
                                           1,
                                       jkj::util::abs(maximizer)}),
            approx_x.denominator)
            .above;

    if (lower_bound < upper_bound) {
        return jkj::interval<jkj::frac<jkj::bigint::int_var, jkj::bigint::uint_var>,
                             jkj::interval_type_t::bounded_left_open_right_closed>{
            std::move(lower_bound), std::move(upper_bound)};
    }
    else {
        return jkj::interval<jkj::frac<jkj::bigint::int_var, jkj::bigint::uint_var>,
                             jkj::interval_type_t::empty>{};
    }
}

int main() {
    while (true) {
        std::cout << "Input x = 2^D / 10^k:\n";
        std::size_t D;
        while (true) {
            std::cout << "  D: ";
            std::cin >> D;
            if (std::cin) {
                break;
            }
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<int>::max(), '\n');
        }
        std::size_t k;
        while (true) {
            std::cout << "  k: ";
            std::cin >> k;
            if (std::cin) {
                break;
            }
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<int>::max(), '\n');
        }

        std::cout << "Input the range of n:\n";
        std::uint64_t nmin, nmax;
        while (true) {
            while (true) {
                std::cout << "  nmin: ";
                std::cin >> nmin;
                if (std::cin && nmin > 0) {
                    break;
                }
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<int>::max(), '\n');
            }
            while (true) {
                std::cout << "  nmax: ";
                std::cin >> nmax;
                if (std::cin && nmax > 0) {
                    break;
                }
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<int>::max(), '\n');
            }
            if (nmin <= nmax) {
                break;
            }
        }
        std::cout << "\n";

        namespace cntfrc = jkj::cntfrc;

        {
            std::cout << "[Analysis for the floor case]\n";
            auto const admissible_range = jeaiii_analysis_floor(
                cntfrc::make_generator<cntfrc::interval_tracker>(cntfrc::impl::rational{
                    jkj::util::to_signed(jkj::bigint::uint_var::power_of_2(D)),
                    jkj::util::pow_uint(jkj::bigint::uint_var{10}, k)}),
                jkj::interval<jkj::bigint::int_var, jkj::interval_type_t::bounded_closed>{
                    jkj::bigint::int_var{jkj::bigint::sign_t::positive, nmin},
                    jkj::bigint::int_var{jkj::bigint::sign_t::positive, nmax}});

            std::cout << "The interval of admissible xi is ";
            admissible_range.visit([&](auto&& itv) {
                using itv_type = std::remove_cvref_t<decltype(itv)>;
                if constexpr (itv_type::interval_type() == jkj::interval_type_t::empty) {
                    std::cout << "empty.\n\n\n";
                }
                else {
                    static_assert(itv_type::interval_type() ==
                                  jkj::interval_type_t::bounded_left_open_right_closed);

                    std::cout << "(" << itv.lower_bound().numerator << " / "
                              << itv.lower_bound().denominator << ", "
                              << itv.upper_bound().numerator << " / "
                              << itv.upper_bound().denominator << "].\n";

                    auto multiply_shift_result = jkj::idiv::find_optimal_multiply_shift(itv);
                    std::cout
                        << "The smallest admissible magic number and the associated shift amount: ("
                        << multiply_shift_result.multiplier << ", "
                        << multiply_shift_result.shift_amount << ").\n\n\n";
                }
            });
        }

        {
            std::cout << "[Analysis for the floor-plus-1 case]\n";
            auto const admissible_range = jeaiii_analysis_floor_plus_one(
                cntfrc::make_generator<cntfrc::interval_tracker>(cntfrc::impl::rational{
                    jkj::util::to_signed(jkj::bigint::uint_var::power_of_2(D)),
                    jkj::util::pow_uint(jkj::bigint::uint_var{10}, k)}),
                jkj::interval<jkj::bigint::int_var, jkj::interval_type_t::bounded_closed>{
                    jkj::bigint::int_var{jkj::bigint::sign_t::positive, nmin},
                    jkj::bigint::int_var{jkj::bigint::sign_t::positive, nmax}});

            std::cout << "The interval of admissible xi is ";
            admissible_range.visit([&](auto&& itv) {
                using itv_type = std::remove_cvref_t<decltype(itv)>;
                if constexpr (itv_type::interval_type() == jkj::interval_type_t::empty) {
                    std::cout << "empty.\n\n\n";
                }
                else {
                    static_assert(itv_type::interval_type() ==
                                  jkj::interval_type_t::bounded_left_open_right_closed);

                    std::cout << "(" << itv.lower_bound().numerator << " / "
                              << itv.lower_bound().denominator << ", "
                              << itv.upper_bound().numerator << " / "
                              << itv.upper_bound().denominator << "].\n";

                    auto multiply_shift_result = jkj::idiv::find_optimal_multiply_shift(itv);
                    std::cout
                        << "The smallest admissible magic number and the associated shift amount: ("
                        << multiply_shift_result.multiplier << ", "
                        << multiply_shift_result.shift_amount << ").\n\n\n";
                }
            });
        }
    }
}
