// Copyright 2022-2024 Junekey Jeon
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

#ifndef JKJ_HEADER_IDIV_BEST_RATIONAL_APPROX
#define JKJ_HEADER_IDIV_BEST_RATIONAL_APPROX

#include "continued_fraction/generator.h"

namespace jkj {
    namespace idiv {
        namespace detail {
            // Common framework for two functions below.
            template <class SemiconvergentType, class ReturnType, class ContinuedFractionGenerator,
                      class UInt, class AfterTerminate, class ConvergentWrapper>
            constexpr ReturnType find_best_rational_approx_impl(ContinuedFractionGenerator&& cf,
                                                                UInt const& nmax,
                                                                AfterTerminate&& after_terminate,
                                                                ConvergentWrapper&& wrapper) {
                auto get_last_semiconvergent = [&cf](auto const& denominator_bound) {
                    auto inverted_semiconvergent_coeff =
                        util::div_ceil(cf.current_convergent_denominator() - denominator_bound,
                                       cf.previous_convergent_denominator());

                    return SemiconvergentType{cf.current_convergent_numerator() -
                                                  inverted_semiconvergent_coeff *
                                                      cf.previous_convergent_numerator(),
                                              util::abs(cf.current_convergent_denominator() -
                                                        inverted_semiconvergent_coeff *
                                                            cf.previous_convergent_denominator())};
                };

                // First, find the last convergent and the last semiconvergent whose denominator is
                // bounded above by nmax.
                for (auto const& state : cf) {
                    if (state.current_convergent_denominator() > nmax) {
                        // If there the last convergent is still not a perfect approximation, then
                        // we return the last semiconvergent and the convergent as the bounds. Which
                        // one is the lower bound and which one is the upper bound is determined by
                        // the parity of state.current_index(). Note that state.current_index() is
                        // the index of the first convergent with the denominator strictly larger
                        // than nmax, so if this index is even, then the semiconvergent is the lower
                        // bound and the convergent is the upepr bound, and if the index is odd,
                        // then the other way around.
                        if (state.current_index() % 2 == 0) {
                            return ReturnType{get_last_semiconvergent(nmax),
                                              wrapper(state.previous_convergent())};
                        }
                        else {
                            return ReturnType{wrapper(state.previous_convergent()),
                                              get_last_semiconvergent(nmax)};
                        }
                    }
                }

                return after_terminate(get_last_semiconvergent);
            }
        }

        template <class ProjectiveRational>
        struct best_rational_approx_output {
            ProjectiveRational below;
            ProjectiveRational above;
        };

        // For a given real number x and a nonnegative integer nmax, find the best rational
        // approximation of it from below and from above whose denominators are no more than nmax.
        // The number x is specified in terms of a continued fraction generator giving its continued
        // fraction expansion. The generator needs to have index_tracker and convergent_tracker
        // within it, and it also needs to be at its initial stage, i.e., any iterator obtained from
        // it has never been advanced. After the function returns, the generator is terminated if x
        // is rational and its denominator is at most nmax.
        template <class ContinuedFractionGenerator, class UInt>
        constexpr auto find_best_rational_approx(ContinuedFractionGenerator&& cf,
                                                 UInt const& nmax) {
            static_assert(
                cntfrc::has_mixins<ContinuedFractionGenerator, cntfrc::index_tracker,
                                   cntfrc::convergent_tracker>(),
                "the passed continued fraction generator must implement index_tracker and "
                "convergent_tracker");

            using convergent_type =
                typename std::remove_cvref_t<ContinuedFractionGenerator>::convergent_type;
            using return_type = best_rational_approx_output<convergent_type>;

            return detail::find_best_rational_approx_impl<convergent_type, return_type>(
                cf, nmax,
                [&](auto&&) -> return_type {
                    return {cf.current_convergent(), cf.current_convergent()};
                },
                [&](auto&& convergent) -> decltype(auto) {
                    return static_cast<decltype(convergent)&&>(convergent);
                });
        }

        // For a given real number x and a positive integer nmax, find the interval
        // [max_n floor(nx)/n, min_n (floor(nx)+1)/n)), where n ranges from {1, ... , nmax}. The
        // number x is specified in terms of a continued fraction generator giving its continued
        // fraction expansion. The generator needs to have index_tracker and convergent_tracker
        // within it, and it also needs to be at its initial stage, i.e., any iterator obtained from
        // it has never been advanced. After the function returns, the generator is terminated if x
        // is rational and its denominator is at most nmax.
        template <class ContinuedFractionGenerator, class UInt>
        constexpr auto find_floor_quotient_range(ContinuedFractionGenerator&& cf,
                                                 UInt const& nmax) {
            static_assert(
                cntfrc::has_mixins<ContinuedFractionGenerator, cntfrc::index_tracker,
                                   cntfrc::convergent_tracker>(),
                "the passed continued fraction generator must implement index_tracker and "
                "convergent_tracker");

            util::constexpr_assert(util::is_strictly_positive(nmax));

            using convergent_type =
                typename std::remove_cvref_t<ContinuedFractionGenerator>::convergent_type;
            using rational_type = decltype(project_to_rational(std::declval<convergent_type>()));
            using return_type =
                interval<rational_type, interval_type_t::bounded_left_closed_right_open>;

            return detail::find_best_rational_approx_impl<rational_type, return_type>(
                cf, nmax,
                [&](auto&& get_last_semiconvergent) -> return_type {
                    // If we reach to the perfect approximation, then we have to find the largest
                    // positive integer v <= nmax such that vp == -1 (mod q), where x = p/q. Then
                    // the lower bound is p/q, while the upper bound is ((vp+1)/q) / v.
                    // To compute v, we find the modular inverse b of -p. This can be done by
                    // observing that the best rational approximation from above whose denominator
                    // is strictly less than q must be precisely ((bp+1)/q) / b. Then
                    //        v = b + floor((nmax - b)/q)q and
                    // (vp+1)/q = (bp+1)/q + floor((nmax - b)/q)p.

                    // If we ended at an even convergent, the last convergent is the best rational
                    // approximation from above. Otherwise, the last semiconvergent is the best
                    // rational approximation from above.

                    auto upper_bound =
                        cf.current_index() % 2 == 0
                            ? project_to_rational(cf.previous_convergent())
                            : get_last_semiconvergent(cf.current_convergent_denominator() - 1);

                    // At this point, upper_bound is ((bp+1)/q) / b, so we adjust the numerator and
                    // the denominator by floor((nmax - b)/q).
                    auto max_quotient =
                        (nmax - upper_bound.denominator) / cf.current_convergent().denominator;
                    upper_bound.numerator += max_quotient * cf.current_convergent().numerator;
                    upper_bound.denominator += max_quotient * cf.current_convergent().denominator;

                    return return_type{project_to_rational(cf.current_convergent()),
                                       std::move(upper_bound)};
                },
                [&](auto&& convergent) -> rational_type {
                    return project_to_rational(static_cast<decltype(convergent)&&>(convergent));
                });
        }

        template <class Int>
        struct extremizers_of_fractional_part {
            Int smallest_minimizer;
            Int largest_maximizer;
        };

        // For a given real number x and a positive integer nmax, find
        // min argmin_n (nx - floor(nx)) and max argmax_n (nx - floor(nx)). The number x is
        // specified in terms of a continued fraction generator giving its continued fraction
        // expansion. The generator needs to have index_tracker and convergent_tracker within it,
        // and it also needs to be at its initial stage, i.e., any iterator obtained from it has
        // never been advanced. After the function returns, the generator is terminated if x is
        // rational and its denominator is at most nmax.
        template <class ContinuedFractionGenerator, class UInt>
        constexpr auto find_extremizers_of_fractional_part(ContinuedFractionGenerator&& cf,
                                                           UInt const& nmax) {
            static_assert(
                cntfrc::has_mixins<ContinuedFractionGenerator, cntfrc::index_tracker,
                                   cntfrc::convergent_tracker>(),
                "the passed continued fraction generator must implement index_tracker and "
                "convergent_tracker");

            util::constexpr_assert(util::is_strictly_positive(nmax));

            using convergent_type =
                typename std::remove_cvref_t<ContinuedFractionGenerator>::convergent_type;
            using return_type =
                extremizers_of_fractional_part<decltype(convergent_type::denominator)>;

            auto floor_quotient_range =
                find_floor_quotient_range(std::forward<ContinuedFractionGenerator>(cf), nmax);
            return return_type{floor_quotient_range.lower_bound().denominator,
                               floor_quotient_range.upper_bound().denominator};
        }
    }
}

#endif
