// Copyright 2023 Junekey Jeon
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

#ifndef JKJ_HEADER_LOG_CONTINUED_FRACTION
#define JKJ_HEADER_LOG_CONTINUED_FRACTION

#include "rational_continued_fraction.h"
#include "gosper_continued_fraction.h"

namespace jkj {
    namespace cntfrc {
        namespace impl {
            // Computes an arbitrary-precision rational approximation of ln(a/b) for positive
            // integers
            // a, b using the continued fraction expansion of 2atanh(z) = ln((1+z)/(1-z)):
            //
            // 2atanh(z) = 2z / (1 - z^2 / (3 - 4z^2 / (5 - 9z^2 / (7 - 16z^2 / ... ) ) ) )
            //
            // By specializing the general recurrence relation for generalized continued fractions,
            // we obtain
            //
            // A_0 = 0, A_1 = 2z, A_n = (2n-1)A_(n-1) - (n-1)^2 z^2 A_(n-2) for n >= 2,
            // B_0 = 1, B_1 = 1,  B_n = (2n-1)B_(n-1) - (n-1)^2 z^2 B_(n-2) for n >= 2,
            //
            // where A_n/B_n is the nth convergent. By substituting z <- p/q and letting P_n = q^n
            // A_n, Q_n = q^n B_n, we obtain the recurrence relation
            //
            // P_0 = 0, P_1 = 2p, P_n = (2n-1)qP_(n-1) - (n-1)^2 p^2 P_(n-2) for n >= 2,
            // Q_0 = 1, Q_1 = q,  Q_n = (2n-1)qQ_(n-1) - (n-1)^2 p^2 Q_(n-2) for n >= 2.
            //
            // Then it can be shown that the denominators Q_n's are always positive and the sequence
            // (P_n/Q_n)_n of convergents converges to 2atanh(z).

            template <class Int, class UInt>
            class natural_log_calculator {
            public:
                using partial_fraction_type = projective_rational<Int, Int>;
                using convergent_type = projective_rational<Int, UInt>;
                using interval_type = variable_shape_cyclic_interval<
                    convergent_type, cyclic_interval_type_t::single_point,
                    cyclic_interval_type_t::open, cyclic_interval_type_t::entire>;

                using required_mixins =
                    mixin_list<index_tracker, partial_fraction_tracker, convergent_tracker>;
                using mixin_ordering_constraints = mixin_ordering_constraint::constraint_list<
                    mixin_ordering_constraint::before_after<index_tracker, interval_tracker>,
                    mixin_ordering_constraint::before_after<convergent_tracker, interval_tracker>>;

            private:
                // Used as a temporary storage for the quantity |p|/q for initial indices.
                frac<UInt, UInt> current_error_bound_{1u, 0u};
                UInt p_square_ = 0u;
                UInt two_q_ = 0u;
                UInt previous_denominator = 0u;
                bool is_z_negative_ = false;

            public:
                static constexpr auto initial_partial_fraction() {
                    return partial_fraction_type{Int{1}, Int{0}};
                }
                static constexpr auto initial_interval() {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }

                explicit constexpr natural_log_calculator(
                    frac<UInt, UInt> const& positive_rational) {
                    util::constexpr_assert(!util::is_zero(positive_rational.denominator) &&
                                           !util::is_zero(positive_rational.numerator));

                    // For given input x, we find the reduced form of z s.t. x = (1+z)/(1-z), i.e.,
                    // z = (x-1)/(x+1). Note that z always lies in (-1,1) for x in (0,infty).
                    auto z = [&] {
                        auto cf = make_generator<convergent_tracker>(
                            rational{projective_rational<Int, UInt>{
                                util::to_signed(positive_rational.numerator) -
                                    positive_rational.denominator,
                                positive_rational.numerator + positive_rational.denominator}});

                        while (!cf.terminated()) {
                            cf.update();
                        }
                        return cf.current_convergent();
                    }();

                    is_z_negative_ = util::is_strictly_negative(z.numerator);
                    current_error_bound_.numerator = util::abs(z.numerator);
                    current_error_bound_.denominator = z.denominator;
                    two_q_ = (z.denominator << 1);
                    p_square_ = util::abs(static_cast<Int&&>(z.numerator));
                    p_square_ *= p_square_;
                }

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    auto const& gen = callback.get_generator();
                    if (gen.current_index() == -1) {
                        callback(partial_fraction_type{1, 0u});
                    }
                    else if (gen.current_index() == 0) {
                        current_error_bound_.numerator <<= 1;
                        callback(partial_fraction_type{
                            is_z_negative_ ? util::to_negative(static_cast<UInt&&>(
                                                 current_error_bound_.numerator))
                                           : util::to_signed(static_cast<UInt&&>(
                                                 current_error_bound_.numerator)),
                            Int{static_cast<UInt&&>(current_error_bound_.denominator)}});
                    }
                    else {
                        auto result = gen.current_partial_fraction();
                        result.numerator = (-gen.current_index() * gen.current_index()) * p_square_;
                        result.denominator += two_q_;
                        callback(static_cast<decltype(result)&&>(result));
                    }
                }

                template <class ContinuedFractionGenerator>
                constexpr interval_type next_interval(ContinuedFractionGenerator const& gen) {
                    util::constexpr_assert(gen.current_index() >= 0);
                    // When p = 0.
                    if (util::is_zero(p_square_)) {
                        return cyclic_interval<convergent_type,
                                               cyclic_interval_type_t::single_point>{
                            convergent_type{0, 1u}};
                    }

                    if (gen.current_index() == 0) {
                        // |p|/q -> 2|p|q/(q^2 - p^2)
                        auto error_bound = current_error_bound_;
                        error_bound.numerator *= error_bound.denominator;
                        error_bound.numerator <<= 1;
                        error_bound.denominator *= error_bound.denominator;
                        error_bound.denominator -= p_square_;

                        if (is_z_negative_) {
                            return cyclic_interval<convergent_type, cyclic_interval_type_t::open>{
                                convergent_type{
                                    util::to_negative(static_cast<UInt&&>(error_bound.numerator)),
                                    static_cast<UInt&&>(error_bound.denominator)},
                                convergent_type{0, 1u}};
                        }
                        else {
                            return cyclic_interval<convergent_type, cyclic_interval_type_t::open>{
                                convergent_type{0, 1u},
                                convergent_type{
                                    util::to_signed(static_cast<UInt&&>(error_bound.numerator)),
                                    static_cast<UInt&&>(error_bound.denominator)}};
                        }
                    }
                    else if (gen.current_index() == 1) {
                        // 2|p|/q -> 2|p|^3/(2q(q^2 - p^2))
                        previous_denominator = gen.current_convergent_denominator();
                        current_error_bound_.numerator =
                            abs(gen.current_convergent_numerator()) * p_square_;
                        current_error_bound_.denominator = gen.current_convergent_denominator();
                        current_error_bound_.denominator *= current_error_bound_.denominator;
                        current_error_bound_.denominator -= p_square_;
                        current_error_bound_.denominator *= two_q_;
                    }
                    else if (gen.current_index() >= 2) {
                        current_error_bound_.numerator *=
                            (unsigned(gen.current_index() - 1) * unsigned(gen.current_index()));
                        current_error_bound_.numerator *= p_square_;
                        current_error_bound_.denominator /=
                            (unsigned(gen.current_index()) * previous_denominator);
                        current_error_bound_.denominator *= (unsigned(gen.current_index() + 1) *
                                                             gen.current_convergent_denominator());
                        previous_denominator = gen.previous_convergent_denominator();
                    }

                    if (is_z_negative_) {
                        auto lower_bound = linear_fractional_translation(
                            util::to_negative(current_error_bound_.numerator),
                            current_error_bound_.denominator)(gen.current_convergent());
                        return cyclic_interval<convergent_type, cyclic_interval_type_t::open>{
                            convergent_type{
                                util::is_strictly_negative(lower_bound.denominator)
                                    ? util::invert_sign(std::move(lower_bound.numerator))
                                    : std::move(lower_bound.numerator),
                                std::move(lower_bound.denominator)},
                            gen.current_convergent()};
                    }
                    else {
                        return cyclic_interval<convergent_type, cyclic_interval_type_t::open>{
                            gen.current_convergent(),
                            linear_fractional_translation(current_error_bound_.numerator,
                                                          current_error_bound_.denominator)(
                                gen.current_convergent())};
                    }
                }
            };

            template <class Int, class UInt, class Unity = unity>
            class natural_log {
                using internal_impl_type = natural_log_calculator<Int, UInt>;
                using impl_type = unary_gosper<internal_impl_type, Unity>;

                impl_type impl_;

            public:
                using partial_fraction_type = typename impl_type::partial_fraction_type;
                using convergent_type = typename impl_type::convergent_type;
                using interval_type = typename impl_type::interval_type;

                static constexpr partial_fraction_type initial_partial_fraction() {
                    return {Unity{}, Int{0}};
                }
                static constexpr interval_type initial_interval() {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }

                explicit constexpr natural_log(frac<UInt, UInt> const& positive_rational)
                    : impl_{internal_impl_type{positive_rational}, {1, 0, 0, 1}} {}

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    impl_.with_next_partial_fraction(callback);
                }
            };

            template <class Int, class UInt, class Unity = unity>
            class general_log {
                using rational_impl_type = rational<Int, UInt, Unity>;
                using irrational_internal_impl_type = natural_log_calculator<Int, UInt>;
                using irrational_impl_type = binary_gosper<irrational_internal_impl_type,
                                                           irrational_internal_impl_type, Unity>;

            public:
                using partial_fraction_type = typename rational_impl_type::partial_fraction_type;
                using convergent_type = typename rational_impl_type::convergent_type;
                using interval_type = typename rational_impl_type::interval_type;

                static_assert(std::is_same_v<partial_fraction_type,
                                             typename irrational_impl_type::partial_fraction_type>);
                static_assert(std::is_same_v<convergent_type,
                                             typename irrational_impl_type::convergent_type>);
                static_assert(
                    std::is_same_v<interval_type, typename irrational_impl_type::interval_type>);

            private:
                rational_impl_type rational_impl_;
                irrational_impl_type irrational_impl_;
                bool is_rational_;

                struct check_rational_return {
                    bool is_rational = false;
                    projective_rational<Int, UInt> result{1, 0u};
                };

                static constexpr check_rational_return check_rational(frac<UInt, UInt> const& base,
                                                                      frac<UInt, UInt> const& x) {
                    // Both base and x must lie in the interval (0, infinity).
                    util::constexpr_assert(!util::is_zero(base.numerator) &&
                                           !util::is_zero(base.denominator));
                    util::constexpr_assert(!util::is_zero(x.numerator) &&
                                           !util::is_zero(x.denominator));

                    // Log with base 1 naturally corresponds to infinity.
                    if (base.numerator == base.denominator) {
                        return {true, projective_rational<Int, UInt>{1, 0u}};
                    }

                    // This algorithm is a refinement of a suggestion by Seok-Hyeong Lee.
                    // Suppose we want to see if log_a b = m/n holds, which means a^m = b^n.
                    // If a = p_1^e_1 ... p_d^e_d is the prime factorization of a, then
                    // p_1^(me_1) ... p_d^(me_d) must be the prime factorization of a^m = b^n.
                    // This means that p_1, ... p_d should be the complete list of prime factors
                    // appearing in the prime factorization of b^n, and since b^n is an nth power,
                    // every exponent me_i must be a multiple of n. However, since gcd(m,n) = 1, it
                    // follows that n divides every e_i. In other words, a = c^n must hold for some
                    // rational number c, which then implies b = c^m.
                    // Now, write a = p/q, b = r/s and c = t/u, then we must have
                    // p/q = t^n/u^n and r/s = t^m/u^m. Since both of p/q and t^n/u^n are already of
                    // their reduced fraction form, we must have p = t^n and q = u^n.
                    // Similarly, we must have r = t^m and s = u^m if m >= 0, and
                    // r = u^-m and s = t^-m if m < 0.

                    // So first, we decide which one is the case, that is, to see the sign of
                    // log_a b. Since log_a b >= 0 if and only if
                    // either a > 1 and b >= 1 or a < 1 and b <= 1,
                    // we simply inspect these inequalities. When log_a b < 0, replace b by 1/b so
                    // that we always have m >= 0, p = t^n, r = t^m, q = u^n, and s = u^m.
                    bool const log_is_nonnegative =
                        (base.numerator > base.denominator && x.numerator >= x.denominator) ||
                        (base.numerator < base.denominator && x.numerator <= x.denominator);
                    auto p = base.numerator;
                    auto q = base.denominator;
                    auto r = log_is_nonnegative ? x.numerator : x.denominator;
                    auto s = log_is_nonnegative ? x.denominator : x.numerator;

                    // Next, we determine which one between m and n is larger. Since m/n >= 1
                    // if and only if either a > 1 and b >= a or a < 1 and b <= a,
                    // we simply inspect these inequalities.
                    // When m < n, swap a and b so that we always have m >= n.
                    bool const numerator_is_greater_than_or_equal_to = [&] {
                        auto ps = p * s;
                        auto qr = q * r;
                        return (p > q && ps <= qr) || (p < q && ps >= qr);
                    }();
                    if (!numerator_is_greater_than_or_equal_to) {
                        using util::swap;
                        swap(p, r);
                        swap(q, s);
                    }

                    // We now solve p = t^n, r = t^m.
                    projective_rational<UInt, UInt> convergent{1u, 0u}, previous_convergent{0u, 1u};
                    while (true) {
                        // Since m >= n, p must divide r, so we factor out the maximum power of p
                        // from r and write r = p^d * p', r' = p, n' = m - dn, and m' = n so that
                        // p' = t^n' and r' = t^m' hold.
                        // Then we will substitute p <- p', r <- r', m <- m', n <- n', and then
                        // iterate this procedure until we reach to (m',n') = (1,0).
                        // Note that reconstructing (m, n) from (m',n') can be done by multiplying
                        // the coefficient matrix [d 1;1 0], so we cummulatively multiply [d 1;1 0]
                        // to right at each iteration step.
                        // The coefficient d should be continued fraction coefficient of m/n, and
                        // the coefficient matrix we are keeping track of is nothing but the matrix
                        // consisting of two consecutive convergents.
                        UInt d = 0u;
                        while (true) {
                            auto div_result = util::div(r, p);
                            if (!util::is_zero(div_result.rem)) {
                                using util::swap;
                                swap(r, p);

                                swap(convergent, previous_convergent);
                                convergent.numerator += d * previous_convergent.numerator;
                                convergent.denominator += d * previous_convergent.denominator;

                                break;
                            }

                            ++d;
                            r = static_cast<decltype(div_result.quot)&&>(div_result.quot);
                        }

                        // If the invariant "p divides r" fails to hold, then we conclude that
                        // log_a b must be irrational.
                        if (util::is_zero(d)) {
                            return {};
                        }

                        // If it is impossible to write s = q^d * q', s' = q, then we conclude that
                        // log_a b must be irrational.
                        {
                            auto q_to_d = util::pow_uint(q, d);
                            auto div_result = util::div(s, q_to_d);
                            if (!util::is_zero(div_result.rem)) {
                                return {};
                            }
                            s = static_cast<UInt&&>(q);
                            q = static_cast<decltype(div_result.quot)&&>(div_result.quot);
                        }

                        // Since m and n are coprime, we must have m' = 1 if n' = 0, which is the
                        // termination condition.
                        if (p == 1u) {
                            // At this point, q must be equal to 1 as well. Otherwise, we conclude
                            // that log_a b must be irrational.
                            if (q != 1u) {
                                return {};
                            }
                            break;
                        }
                    }

                    // From (m',n') = (1,0), we obtain (m,n).
                    if (!numerator_is_greater_than_or_equal_to) {
                        using util::swap;
                        swap(convergent.numerator, convergent.denominator);
                    }
                    return {true,
                            projective_rational<Int, UInt>{
                                log_is_nonnegative
                                    ? util::to_signed(static_cast<UInt&&>(convergent.numerator))
                                    : util::to_negative(static_cast<UInt&&>(convergent.numerator)),
                                static_cast<UInt&&>(convergent.denominator)}};
                }

                explicit constexpr general_log(frac<UInt, UInt> const& base,
                                               frac<UInt, UInt> const& x,
                                               check_rational_return rational_check)
                    : rational_impl_{rational_check.result},
                      irrational_impl_{irrational_internal_impl_type{x},
                                       irrational_internal_impl_type{base},
                                       // (0xy + 1x + 0y + 0) / (0xy + 0x + 1y + 0)
                                       {0, 1, 0, 0, 0, 0, 1, 0}},
                      is_rational_{rational_check.is_rational} {}

            public:
                static constexpr partial_fraction_type initial_partial_fraction() {
                    return {Unity{}, Int{0}};
                }
                static constexpr interval_type initial_interval() {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
                }

                explicit constexpr general_log(frac<UInt, UInt> const& base,
                                               frac<UInt, UInt> const& x)
                    : general_log{base, x, check_rational(base, x)} {}

                template <class Callback>
                constexpr void with_next_partial_fraction(Callback&& callback) {
                    if (is_rational_) {
                        rational_impl_.with_next_partial_fraction(callback);
                    }
                    else {
                        irrational_impl_.with_next_partial_fraction(callback);
                    }
                }
            };
        }
    }
}

#endif
