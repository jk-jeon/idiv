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
#include "prime_factorization.h"

namespace jkj {
    namespace cntfrc {
        // Computes an arbitrary-precision rational approximation of ln(a/b) for positive integers
        // a, b using the continued fraction expansion of 2atanh(z) = ln((1+z)/(1-z)):
        //
        // 2atanh(z) = 2z / (1 - z^2 / (3 - 4z^2 / (5 - 9z^2 / (7 - 16z^2 / ... ) ) ) )
        //
        // By specializing the general recurrence relation for generalized continued fractions, we
        // obtain
        //
        // A_0 = 0, A_1 = 2z, A_n = (2n-1)A_(n-1) - (n-1)^2 z^2 A_(n-2) for n >= 2,
        // B_0 = 1, B_1 = 1,  B_n = (2n-1)B_(n-1) - (n-1)^2 z^2 B_(n-2) for n >= 2,
        //
        // where A_n/B_n is the nth convergent. By substituting z <- p/q and letting P_n = q^n A_n,
        // Q_n = q^n B_n, we obtain the recurrence relation
        //
        // P_0 = 0, P_1 = 2p, P_n = (2n-1)qP_(n-1) - (n-1)^2 p^2 P_(n-2) for n >= 2,
        // Q_0 = 1, Q_1 = q,  Q_n = (2n-1)qQ_(n-1) - (n-1)^2 p^2 Q_(n-2) for n >= 2.
        //
        // Then it can be shown that the denominators Q_n's are always positive and the sequence
        // (P_n/Q_n)_n of convergents converges to 2atanh(z).

        template <class Int, class UInt>
        class natural_log_calculator {
        public:
            using partial_fraction_type = frac<Int, Int>;
            using convergent_type = projective_rational<Int, UInt>;
            using interval_type = variable_shape_cyclic_interval<
                convergent_type, cyclic_interval_type_t::single_point, cyclic_interval_type_t::open,
                cyclic_interval_type_t::entire>;

        private:
            // Used as a temporary storage for the quantity |p|/q for initial indices.
            frac<UInt, UInt> current_error_bound_{1u, 0u};
            UInt p_square_ = 0u;
            UInt two_q_ = 0u;
            UInt previous_denominator = 0u;
            bool is_z_negative_ = false;

        public:
            static constexpr partial_fraction_type initial_partial_fraction() {
                return {Int{1}, Int{0}};
            }
            static constexpr interval_type initial_interval() {
                return cyclic_interval<convergent_type, cyclic_interval_type_t::entire>{};
            }

            explicit constexpr natural_log_calculator(frac<UInt, UInt> const& positive_rational) {
                util::constexpr_assert(!is_zero(positive_rational.denominator) &&
                                       !is_zero(positive_rational.numerator));

                // For given input x, we find the reduced form of z s.t. x = (1+z)/(1-z), i.e.,
                // z = (x-1)/(x+1). Note that z always lies in (-1,1) for x in (0,infty).
                auto z = [&] {
                    auto cf = make_continued_fraction_generator<convergent_tracker>(
                        rational_continued_fraction{projective_rational<Int, UInt>{
                            to_signed(positive_rational.numerator) - positive_rational.denominator,
                            positive_rational.numerator + positive_rational.denominator}});

                    while (!cf.terminated()) {
                        cf.update();
                    }
                    return cf.current_convergent();
                }();

                is_z_negative_ = is_strictly_negative(z.numerator);
                current_error_bound_.numerator = abs(z.numerator);
                current_error_bound_.denominator = z.denominator;
                two_q_ = (z.denominator << 1);
                p_square_ = abs(static_cast<Int&&>(z.numerator));
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
                        is_z_negative_
                            ? to_negative(static_cast<UInt&&>(current_error_bound_.numerator))
                            : to_signed(static_cast<UInt&&>(current_error_bound_.numerator)),
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
                if (is_zero(p_square_)) {
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::single_point>{
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
                            convergent_type{to_negative(static_cast<UInt&&>(error_bound.numerator)),
                                            static_cast<UInt&&>(error_bound.denominator)},
                            convergent_type{0, 1u}};
                    }
                    else {
                        return cyclic_interval<convergent_type, cyclic_interval_type_t::open>{
                            convergent_type{0, 1u},
                            convergent_type{to_signed(static_cast<UInt&&>(error_bound.numerator)),
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
                    current_error_bound_.denominator *=
                        (unsigned(gen.current_index() + 1) * gen.current_convergent_denominator());
                    previous_denominator = gen.previous_convergent_denominator();
                }

                if (is_z_negative_) {
                    auto lower_bound = linear_fractional_translation(
                        to_negative(current_error_bound_.numerator),
                        current_error_bound_.denominator)(gen.current_convergent());
                    return cyclic_interval<convergent_type, cyclic_interval_type_t::open>{
                        convergent_type{is_strictly_negative(lower_bound.denominator)
                                            ? invert_sign(std::move(lower_bound.numerator))
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
            using impl_type = unary_gosper<
                continued_fraction<internal_impl_type, index_tracker, partial_fraction_tracker,
                                   convergent_tracker, interval_tracker>,
                Unity>;

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
            using rational_impl_type = rational_continued_fraction<Int, UInt, Unity>;
            using irrational_internal_impl_type = natural_log_calculator<Int, UInt>;
            using irrational_impl_type = binary_gosper<
                continued_fraction<irrational_internal_impl_type, index_tracker,
                                   partial_fraction_tracker, convergent_tracker, interval_tracker>,
                continued_fraction<irrational_internal_impl_type, index_tracker,
                                   partial_fraction_tracker, convergent_tracker, interval_tracker>,
                Unity>;

        public:
            using partial_fraction_type = typename rational_impl_type::partial_fraction_type;
            using convergent_type = typename rational_impl_type::convergent_type;
            using interval_type = typename rational_impl_type::interval_type;

            static_assert(std::is_same_v<partial_fraction_type,
                                         typename irrational_impl_type::partial_fraction_type>);
            static_assert(
                std::is_same_v<convergent_type, typename irrational_impl_type::convergent_type>);
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
                util::constexpr_assert(!is_zero(base.numerator) && !is_zero(base.denominator));
                util::constexpr_assert(!is_zero(x.numerator) && !is_zero(x.denominator));

                // Log with base 1 does not make sense.
                util::constexpr_assert(base.numerator != base.denominator);

                auto prime_factors_base = prime_factorization(base);
                auto prime_factors_x = prime_factorization(x);

                if (prime_factors_x.size() != prime_factors_base.size()) {
                    return {};
                }
                util::constexpr_assert(prime_factors_x.size() != 0);

                auto ratio = projective_rational<int, int>{prime_factors_x[0].exponent,
                                                           prime_factors_base[0].exponent};

                for (std::size_t idx = 1; idx < prime_factors_x.size(); ++idx) {
                    if (prime_factors_x[idx].factor != prime_factors_base[idx].factor) {
                        return {};
                    }

                    if (ratio != projective_rational<int, int>{prime_factors_x[idx].exponent,
                                                               prime_factors_base[idx].exponent}) {
                        return {};
                    }
                }

                return {true, projective_rational<Int, UInt>{
                                  ratio.denominator > 0 ? ratio.numerator : -ratio.numerator,
                                  ratio.denominator > 0 ? unsigned(ratio.denominator)
                                                        : unsigned(-ratio.denominator)}};
            }

            explicit constexpr general_log(frac<UInt, UInt> const& base, frac<UInt, UInt> const& x,
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

            explicit constexpr general_log(frac<UInt, UInt> const& base, frac<UInt, UInt> const& x)
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

#endif
