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


#ifndef JKJ_HEADER_IDIV_BIG_UINT
#define JKJ_HEADER_IDIV_BIG_UINT

#include "wuint.h"
#include <bit>
#include <compare>
#include <limits>
#include <vector>

namespace jkj {
    namespace big_uint {
        using block_type = std::uint64_t;
        static constexpr std::size_t number_of_bits_in_block = 64;
        static constexpr block_type largest_pow10_in_block = UINT64_C(1000'0000'0000'0000'0000);
        static_assert(std::numeric_limits<block_type>::max() >= largest_pow10_in_block);

        // Implementations of arithemetic operations.
        namespace detail {
            // All of these functions require the parameters to have nonzero sizes.

            // Precondition: no leading zero blocks.
            template <class X, class Y>
            constexpr std::strong_ordering comparison_impl(X const& x, Y const& y) noexcept {
                if (x.size() < y.size()) {
                    return std::strong_ordering::less;
                }
                else if (x.size() > y.size()) {
                    return std::strong_ordering::greater;
                }

                for (std::size_t idx_p1 = x.size(); idx_p1 > 0; --idx_p1) {
                    if (x[idx_p1 - 1] < y[idx_p1 - 1]) {
                        return std::strong_ordering::less;
                    }
                    else if (x[idx_p1 - 1] > y[idx_p1 - 1]) {
                        return std::strong_ordering::greater;
                    }
                }
                return std::strong_ordering::equal;
            }

            // Precondition: x.size() >= y.size().
            // Modifies x in-place. Returns true if a leading block with the carry bit should be
            // appended.
            template <class X, class Y>
            constexpr decltype(auto) add_impl(X&& x, Y const& y) {
                util::constexpr_assert<util::error_msgs::no_error_msg>(x.size() >= y.size());

                unsigned int carry = 0;
                for (std::size_t idx = 0; idx < y.size(); ++idx) {
                    auto with_carry = x[idx] + carry;
                    unsigned int first_carry = (with_carry < x[idx]) ? 1 : 0;

                    auto const y_block = y[idx];
                    x[idx] = with_carry + y_block;
                    carry = first_carry | ((x[idx] < y_block) ? 1 : 0);
                }

                if (carry != 0) {
                    for (std::size_t idx = y.size(); idx < x.size(); ++idx) {
                        ++x[idx];
                        if (x[idx] != 0) {
                            return false;
                        }
                    }
                    // Carry.
                    return true;
                }
                return false;
            }

            // Precondition: x.size() >= y.size() and x >= y.
            // Modifies x in-place. x may contain leading zeros after the function has executed.
            template <class X, class Y>
            constexpr void subtract_impl(X&& x, Y const& y) {
                util::constexpr_assert<util::error_msgs::underflow>(x.size() >= y.size());

                unsigned int borrow = 0;
                for (std::size_t idx = 0; idx < y.size(); ++idx) {
                    auto with_borrow = x[idx] - borrow;
                    unsigned int first_borrow = (with_borrow > x[idx]) ? 1 : 0;

                    auto y_block = y[idx];
                    x[idx] = with_borrow - y_block;
                    borrow = first_borrow | ((x[idx] > with_borrow) ? 1 : 0);
                }

                if (borrow != 0) {
                    util::constexpr_assert<util::error_msgs::underflow>(x.size() > y.size());
                    for (std::size_t idx = y.size(); idx < x.size(); ++idx) {
                        --x[idx];
                        if (x[idx] != std::numeric_limits<block_type>::max()) {
                            return;
                        }
                    }
                    util::constexpr_assert<util::error_msgs::underflow>(
                        x.back() != std::numeric_limits<block_type>::max());
                }
            }

            // Precondition: result is zero-initialized, and has the size x.size() + y.size().
            // result should not alias with x, y.
            template <class X, class Y, class Result>
            constexpr void multiply_impl(X const& x, Y const& y, Result& result) {
                for (std::size_t y_idx = 0; y_idx < y.size(); ++y_idx) {
                    // Compute y[y_idx] * x and accumulate it into the result
                    for (std::size_t x_idx = 0; x_idx < x.size(); ++x_idx) {
                        auto mul = wuint::umul128(x[x_idx], y[y_idx]);

                        // Add the first half
                        result[x_idx + y_idx] += mul.low();
                        unsigned int carry = result[x_idx + y_idx] < mul.low() ? 1 : 0;

                        // Add the second half
                        auto with_carry = mul.high() + carry;
                        carry = with_carry < mul.high() ? 1 : 0;
                        result[x_idx + y_idx + 1] += with_carry;

                        // If there is carry,
                        if (result[x_idx + y_idx + 1] < with_carry) {
                            // Propagate.
                            util::constexpr_assert<util::error_msgs::no_error_msg>(
                                x_idx + y_idx + 2 < result.size());
                            for (auto idx = x_idx + y_idx + 2; idx < result.size(); ++idx) {
                                ++result[idx];
                                if (result[idx] != 0) {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // - x: dividend. Should not have any leading zero blocks. Modified in-place and becomes
            // the remainder after the function is executed.
            // - y: divisor. Should not have any leading zero blocks.
            // - q: quotient. Should be zero-initialized and has the size at least x.size() -
            // y.size() + 1.
            // Returns the number of blocks in the remainder.
            template <class X, class Y, class Quotient>
            constexpr std::size_t long_division_impl(X&& x, Y const& y, Quotient& q) {
                util::constexpr_assert<util::error_msgs::no_error_msg>(q.size() + y.size() >=
                                                                       x.size() + 1);

                auto x_size = x.size();
                auto x_leading_one_pos_in_leading_block = std::size_t(std::bit_width(x.back()));
                auto x_leading_one_pos =
                    x_size * number_of_bits_in_block + x_leading_one_pos_in_leading_block;

                auto const y_leading_one_pos_in_leading_block =
                    std::size_t(std::bit_width(y.back()));
                auto const y_leading_one_pos =
                    y.size() * number_of_bits_in_block + y_leading_one_pos_in_leading_block;

                // Return if x is smaller than y.
                while (x_leading_one_pos >= y_leading_one_pos) {
                    // Compare x and y after aligning their leading 1's.
                    // If x is smaller, then decrease the shift amount by 1.
                    std::size_t block_shift, bit_shift;
                    {
                        auto const total_shift = x_leading_one_pos - y_leading_one_pos;
                        block_shift = total_shift / number_of_bits_in_block;
                        bit_shift = total_shift % number_of_bits_in_block;
                    }

                    if (bit_shift != 0) {
                        auto const should_shift_less = [&] {
                            // The first block of y.
                            if (y.size() + block_shift < x_size) {
                                auto const y_block =
                                    (y.back() >> (number_of_bits_in_block - bit_shift));
                                auto const x_block = x[y.size() + block_shift];
                                if (x_block < y_block) {
                                    return true;
                                }
                                else if (x_block > y_block) {
                                    return false;
                                }
                            }

                            // Middle blocks of y.
                            for (std::size_t idx = y.size() - 1; idx > 0; --idx) {
                                auto const y_block =
                                    (y[idx] << bit_shift) |
                                    (y[idx - 1] >> (number_of_bits_in_block - bit_shift));
                                auto const x_block = x[idx + block_shift];

                                if (x_block < y_block) {
                                    return true;
                                }
                                else if (x_block > y_block) {
                                    return false;
                                }
                            }

                            // The last block of y.
                            return x[block_shift] < (y[0] << bit_shift);
                        }();

                        bit_shift = should_shift_less ? bit_shift - 1 : bit_shift;
                    }
                    else {
                        auto const should_shift_less = [&] {
                            for (std::size_t idx_p1 = y.size(); idx_p1 > 0; --idx_p1) {
                                auto const x_block = x[idx_p1 - 1 + block_shift];
                                auto const y_block = y[idx_p1 - 1];

                                if (x_block < y_block) {
                                    return true;
                                }
                                else if (x_block > y_block) {
                                    return false;
                                }
                            }
                            return false;
                        }();

                        if (should_shift_less) {
                            // Return if x is smaller than y.
                            if (block_shift == 0) {
                                break;
                            }

                            bit_shift = number_of_bits_in_block - 1;
                            --block_shift;
                        }
                    }

                    // Subtract the aligned y from x.
                    if (bit_shift != 0) {
                        // The last block of y.
                        auto y_block = (y[0] << bit_shift);
                        unsigned int borrow = x[block_shift] < y_block ? 1 : 0;
                        x[block_shift] -= y_block;

                        // Middle blocks of y.
                        for (std::size_t idx = 1; idx < y.size(); ++idx) {
                            auto with_borrow = x[idx + block_shift] - borrow;
                            unsigned int first_borrow =
                                (with_borrow > x[idx + block_shift]) ? 1 : 0;

                            y_block = (y[idx] << bit_shift) |
                                      (y[idx - 1] >> (number_of_bits_in_block - bit_shift));
                            x[idx + block_shift] = with_borrow - y_block;
                            borrow = first_borrow | ((x[idx + block_shift] > with_borrow) ? 1 : 0);
                        }

                        // The first block of y.
                        if (y.size() + block_shift < x_size) {
                            auto with_borrow = x[y.size() + block_shift] - borrow;
                            unsigned int first_borrow =
                                (with_borrow > x[y.size() + block_shift]) ? 1 : 0;

                            auto const y_block =
                                (y.back() >> (number_of_bits_in_block - bit_shift));
                            util::constexpr_assert<util::error_msgs::no_error_msg>(with_borrow >=
                                                                                   y_block);
                            x[y.size() + block_shift] = with_borrow - y_block;
                            borrow =
                                first_borrow | ((x[y.size() + block_shift] > with_borrow) ? 1 : 0);
                        }
                        util::constexpr_assert<util::error_msgs::no_error_msg>(borrow == 0);
                    }
                    else {
                        unsigned int borrow = 0;
                        for (std::size_t idx = 0; idx < y.size(); ++idx) {
                            auto with_borrow = x[idx + block_shift] - borrow;
                            unsigned int first_borrow =
                                (with_borrow > x[idx + block_shift]) ? 1 : 0;

                            x[idx + block_shift] = with_borrow - y[idx];
                            borrow = first_borrow | ((x[idx + block_shift] > with_borrow) ? 1 : 0);
                        }

                        if (borrow != 0) {
                            util::constexpr_assert<util::error_msgs::no_error_msg>(
                                y.size() + block_shift < x_size);
                            util::constexpr_assert<util::error_msgs::no_error_msg>(
                                x[y.size() + block_shift] != 0);
                            --x[y.size() + block_shift];
                        }
                    }

                    // Update quotient.
                    q[block_shift] |= (block_type(1) << bit_shift);

                    // Recalculate the leading 1 position of x.
                    while (x[x_size - 1] == 0) {
                        if (--x_size == 0) {
                            // x is divisible by y.
                            return x_size;
                        }
                    }
                    x_leading_one_pos_in_leading_block = std::size_t(std::bit_width(x[x_size - 1]));
                    x_leading_one_pos =
                        x_size * number_of_bits_in_block + x_leading_one_pos_in_leading_block;
                }

                return x_size;
            }
        }

        // Implememtation details for big_uint::constant: big unsigned integer constants whose
        // values are encoded in their types.

        // The first block is the least significant block and the last block is the most
        // significant block.
        template <std::size_t N>
        using static_block_holder = util::array<block_type, N>;

        namespace detail {
            // Reverse the order of the array.
            template <std::size_t N>
            constexpr static_block_holder<N> reverse(static_block_holder<N> const& arr) noexcept {
                static_block_holder<N> result{};
                for (std::size_t idx = 0; idx < N; ++idx) {
                    result[idx] = arr[N - idx - 1];
                }
                return result;
            }

            // Take the first N blocks and discard the remaining.
            template <std::size_t N, class Container>
            constexpr static_block_holder<N> slice(Container const& arr) noexcept {
                util::constexpr_assert<util::error_msgs::index_out_of_range>(N <= arr.size());
                if constexpr (N != 0) {
                    static_block_holder<N> sliced{};
                    for (std::size_t idx = 0; idx < N; ++idx) {
                        sliced[idx] = arr[idx];
                    }
                    return sliced;
                }
                else {
                    return static_block_holder<0>{};
                }
            }

            // Count the number of blocks excluding leading zero blocks.
            template <std::size_t N>
            constexpr std::size_t
            count_blocks_excluding_leading_zeros(static_block_holder<N> arr) noexcept {
                if constexpr (N == 0) {
                    return 0;
                }
                else {
                    std::size_t ret_value = N;
                    for (std::size_t idx_p1 = N; idx_p1 > 0; --idx_p1) {
                        if (arr[idx_p1 - 1] == 0) {
                            --ret_value;
                        }
                        else {
                            break;
                        }
                    }
                    return ret_value;
                }
            }

            namespace adl_guard {
                template <std::size_t N, static_block_holder<N> arr>
                struct constant_impl_impl {
                    using block_holder_type = static_block_holder<N>;
                    static constexpr std::size_t number_of_blocks() noexcept {
                        return block_holder_type::size();
                    }
                    static constexpr block_holder_type blocks = arr;
                };

                template <static_block_holder<1> arr>
                struct constant_impl_impl<1, arr> {
                    using block_holder_type = static_block_holder<1>;
                    static constexpr std::size_t number_of_blocks() noexcept { return 1; }
                    static constexpr block_holder_type blocks = arr;

                    constexpr operator block_type() const noexcept { return arr[0]; }
                };

                template <static_block_holder<0> arr>
                struct constant_impl_impl<0, arr> {
                    using block_holder_type = static_block_holder<0>;
                    static constexpr std::size_t number_of_blocks() noexcept { return 0; }
                    static constexpr block_holder_type blocks = arr;

                    constexpr operator block_type() const noexcept { return 0; }
                };
            }

            template <auto arr>
            using constant_impl = adl_guard::constant_impl_impl<arr.size(), arr>;

            // Remove leading zero blocks from arr.
            template <auto arr>
            constexpr auto remove_leading_zero_blocks() noexcept {
                return slice<count_blocks_excluding_leading_zeros(arr)>(arr);
            }

            namespace adl_guard {
                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr std::strong_ordering operator<=>(constant_impl<x>,
                                                           constant_impl<y>) noexcept {
                    return comparison_impl(x, y);
                }
                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr bool operator==(constant_impl<x>, constant_impl<y>) noexcept {
                    return comparison_impl(x, y) == 0;
                }

                template <std::size_t M, static_block_holder<M> x>
                constexpr std::strong_ordering operator<=>(constant_impl<x>,
                                                           block_type y) noexcept {
                    return comparison_impl(x, static_block_holder<1>{y});
                }
                template <std::size_t M, static_block_holder<M> x>
                constexpr std::strong_ordering operator==(constant_impl<x>, block_type y) noexcept {
                    return comparison_impl(x, static_block_holder<1>{y}) == 0;
                }

                template <std::size_t N, static_block_holder<N> y>
                constexpr std::strong_ordering operator<=>(block_type x,
                                                           constant_impl<y>) noexcept {
                    return comparison_impl(static_block_holder<1>{x}, y);
                }
                template <std::size_t N, static_block_holder<N> y>
                constexpr std::strong_ordering operator==(block_type x, constant_impl<y>) noexcept {
                    return comparison_impl(static_block_holder<1>{x}, y) == 0;
                }

                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr auto operator+(constant_impl<x>, constant_impl<y>) noexcept {
                    if constexpr (x.size() >= y.size()) {
                        if constexpr (y.size() == 0) {
                            return constant_impl<x>{};
                        }
                        else {
                            constexpr auto result = [] {
                                // First, compute the addition without the last carry.
                                auto result = util::apply(
                                    [](auto... elems) {
                                        return static_block_holder<x.size() + 1>{elems..., 0};
                                    },
                                    x);
                                // Add the last carry if needed.
                                if (add_impl(result, y)) {
                                    result.back() = 1;
                                }
                                return result;
                            }();

                            // The result we get in the above may contain a leading zero, so remove
                            // it if needed.
                            if constexpr (result.back() == 0) {
                                return constant_impl<slice<x.size()>(result)>{};
                            }
                            else {
                                return constant_impl<result>{};
                            }
                        }
                    }
                    else {
                        return constant_impl<y>{} + constant_impl<x>{};
                    }
                }

                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr auto operator-(constant_impl<x>, constant_impl<y>) noexcept {
                    constexpr auto result = [] {
                        auto result = x;
                        // This line will generate an error if underflow occurs. The compiler will
                        // generate an error message about the call to the function
                        // 'assert_failed<error_msgs::underflow>'.
                        subtract_impl(result, y);
                        return result;
                    }();
                    return constant_impl<remove_leading_zero_blocks<result>()>{};
                }

                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr auto operator*(constant_impl<x>, constant_impl<y>) noexcept {
                    if constexpr (M == 0 || N == 0) {
                        return constant_impl<static_block_holder<0>{}>{};
                    }
                    else {
                        constexpr auto result = [] {
                            static_block_holder<M + N> result{};
                            multiply_impl(x, y, result);
                            return result;
                        }();

                        // The result we get in the above may contain a leading zero, so remove it
                        // if needed.
                        if constexpr (result.back() == 0) {
                            return constant_impl<slice<result.size() - 1>(result)>{};
                        }
                        else {
                            return constant_impl<result>{};
                        }
                    }
                }

                template <class Quotient, class Remainder>
                struct constant_div_t;

                template <std::size_t M, std::size_t N, static_block_holder<M> quotient,
                          static_block_holder<N> remainder>
                struct constant_div_t<constant_impl<quotient>, constant_impl<remainder>> {
                    constant_impl<quotient> quot;
                    constant_impl<remainder> rem;
                };

                template <std::size_t M, std::size_t N>
                struct div_intermediate_result {
                    static_block_holder<M> quot;
                    static_block_holder<N> rem;
                };

                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr auto div(constant_impl<x>, constant_impl<y>) noexcept {
                    // The divisor shall not be zero.
                    static_assert(N != 0);

                    // When the divisor has strictly more blocks than the dividend.
                    if constexpr (N > M) {
                        // Quotient is zero, remainder is the dividend.
                        return constant_div_t<constant_impl<static_block_holder<0>{}>,
                                              constant_impl<x>>{};
                    }
                    else {
                        constexpr div_intermediate_result<M - N + 1, M> result = [] {
                            static_block_holder<M - N + 1> quot{};
                            static_block_holder<M> dividend = x;
                            long_division_impl(dividend, y, quot);
                            return div_intermediate_result<M - N + 1, M>{quot, dividend};
                        }();

                        // Remove leading zero blocks.
                        constexpr auto quot = remove_leading_zero_blocks<result.quot>();
                        constexpr auto rem = remove_leading_zero_blocks<result.rem>();

                        return constant_div_t<constant_impl<quot>, constant_impl<rem>>{};
                    }
                }

                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr auto operator/(constant_impl<x> xx, constant_impl<y> yy) {
                    return div(xx, yy).quot;
                }

                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr auto operator%(constant_impl<x> xx, constant_impl<y> yy) {
                    return div(xx, yy).rem;
                }
            }

            // Tools for conversion into/from decimals.
            namespace adl_guard {
                template <std::size_t N, static_block_holder<N> arr, bool enable>
                struct decimal_constant_impl_impl;

                template <std::size_t N, static_block_holder<N> arr>
                struct decimal_constant_impl_impl<N, arr, true> {
                    using block_holder_type = static_block_holder<N>;
                    static constexpr std::size_t number_of_blocks = block_holder_type::size();
                    static constexpr block_holder_type blocks = arr;
                };

                template <static_block_holder<1> arr>
                struct decimal_constant_impl_impl<1, arr, true> {
                    using block_holder_type = static_block_holder<1>;
                    static constexpr std::size_t number_of_blocks = block_holder_type::size();
                    static constexpr block_holder_type blocks = arr;

                    constexpr operator block_type() const noexcept { return arr[0]; }
                };

                template <static_block_holder<0> arr>
                struct decimal_constant_impl_impl<0, arr, true> {
                    using block_holder_type = static_block_holder<0>;
                    static constexpr std::size_t number_of_blocks = block_holder_type::size();
                    static constexpr block_holder_type blocks = arr;

                    constexpr operator block_type() const noexcept { return 0; }
                };
            }

            template <std::size_t N>
            constexpr bool is_valid_decimal(static_block_holder<N> const& arr) noexcept {
                for (std::size_t i = 0; i < N; ++i) {
                    if (arr[i] >= largest_pow10_in_block) {
                        return false;
                    }
                }
                return true;
            }

            constexpr bool is_valid_decimal(static_block_holder<0> const& arr) noexcept {
                return true;
            }

            template <auto arr>
            using decimal_constant_impl =
                adl_guard::decimal_constant_impl_impl<arr.size(), arr, is_valid_decimal(arr)>;

            namespace adl_guard {
                template <std::size_t N, static_block_holder<N> x>
                constexpr auto to_decimal(constant_impl<x>) noexcept {
                    if constexpr (N == 0) {
                        return decimal_constant_impl<static_block_holder<0>{}>{};
                    }
                    else {
                        // Compute an upper bound on the number of 19-digits blocks that will be
                        // generated. The exact number of blocks is floor(64N/(19 log2(10))) + 1,
                        // and since 72/71 > 64/(19 log2(10)) holds, the following is an upper
                        // bound.
                        constexpr std::size_t number_of_19_digits_blocks = ((72 * N) / 71) + 1;

                        // In the below, the last element of result stores the actual number of
                        // blocks.
                        constexpr static_block_holder<number_of_19_digits_blocks + 1> result = [] {
                            constexpr static_block_holder<1> divisor{largest_pow10_in_block};

                            static_block_holder<number_of_19_digits_blocks + 1> result{};
                            std::size_t idx = 0;
                            std::size_t remaining_binary_blocks = N;
                            static_block_holder<N> dividend = x;
                            static_block_holder<N> quotient{};

                            while (remaining_binary_blocks != 0) {
                                long_division_impl(
                                    util::span<block_type>{dividend.data, remaining_binary_blocks},
                                    divisor, quotient);
                                result[idx++] = dividend.front();
                                dividend = quotient;
                                remaining_binary_blocks =
                                    count_blocks_excluding_leading_zeros(dividend);
                                quotient = static_block_holder<N>{};
                            }

                            // Record the number of blocks and return.
                            result[number_of_19_digits_blocks] = idx;
                            return result;
                        }();

                        // Return the correct slice of the result obtained above.
                        return decimal_constant_impl<slice<result[number_of_19_digits_blocks]>(
                            result)>{};
                    }
                }

                template <std::size_t N, static_block_holder<N> x>
                constexpr auto from_decimal(decimal_constant_impl<x>) noexcept {
                    if constexpr (N == 0) {
                        return constant_impl<static_block_holder<0>{}>{};
                    }
                    else {
                        // The number of blocks required is at most N.
                        constexpr static_block_holder<N + 1> result = [] {
                            constexpr static_block_holder<1> multiplier{largest_pow10_in_block};

                            static_block_holder<N + 1> result{};
                            result[0] = x[N - 1];
                            for (std::size_t idx_p1 = 1; idx_p1 > 0; --idx_p1) {
                                auto current_copy = slice<N>(result);
                                result = static_block_holder<N + 1>{};
                                multiply_impl(current_copy, multiplier, result);
                                auto const has_carry =
                                    add_impl(result, static_block_holder<1>{x[idx_p1 - 1]});
                                util::constexpr_assert<util::error_msgs::no_error_msg>(!has_carry);
                            }

                            util::constexpr_assert<util::error_msgs::no_error_msg>(result[N] == 0);
                            return result;
                        }();

                        return constant_impl<remove_leading_zero_blocks<result>()>{};
                    }
                }

                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr std::strong_ordering operator<=>(decimal_constant_impl<x>,
                                                           decimal_constant_impl<y>) noexcept {
                    return comparison_impl(x, y);
                }
                template <std::size_t M, std::size_t N, static_block_holder<M> x,
                          static_block_holder<N> y>
                constexpr bool operator==(decimal_constant_impl<x>,
                                          decimal_constant_impl<y>) noexcept {
                    return comparison_impl(x, y) == 0;
                }

                template <std::size_t M, static_block_holder<M> x>
                constexpr std::strong_ordering operator<=>(decimal_constant_impl<x>,
                                                           block_type y) noexcept {
                    return comparison_impl(x, static_block_holder<1>{y});
                }
                template <std::size_t M, static_block_holder<M> x>
                constexpr std::strong_ordering operator==(decimal_constant_impl<x>,
                                                          block_type y) noexcept {
                    return comparison_impl(x, static_block_holder<1>{y}) == 0;
                }

                template <std::size_t N, static_block_holder<N> y>
                constexpr std::strong_ordering operator<=>(block_type x,
                                                           decimal_constant_impl<y>) noexcept {
                    return comparison_impl(static_block_holder<1>{x}, y);
                }
                template <std::size_t N, static_block_holder<N> y>
                constexpr std::strong_ordering operator==(block_type x,
                                                          decimal_constant_impl<y>) noexcept {
                    return comparison_impl(static_block_holder<1>{x}, y) == 0;
                }
            }
        }

        // The first parameter is the most significant block.
        // E.g., constant<3, 7, 11> represents 3 * 2^128 + 7 * 2^64 + 11.
        // All leading zero blocks will be automatically removed, so for e.g., constant<0, 0, 1> is
        // the same type as constant<1>.
        template <block_type... blocks>
        using constant = detail::constant_impl<detail::remove_leading_zero_blocks<detail::reverse(
            static_block_holder<sizeof...(blocks)>{blocks...})>()>;

        // A sequence of decimal numbers consisting of 19-digits.
        // E.g., decimal_constant<1234, 1234> represents 1234'000'0000'0000'0000'1234.
        // The order convention is same as above.
        template <block_type... blocks>
        using decimal_constant =
            detail::decimal_constant_impl<detail::remove_leading_zero_blocks<detail::reverse(
                static_block_holder<sizeof...(blocks)>{blocks...})>()>;

        // Dynamically-sized big unsigned integers.
        class var {
            // Least significant element first.
            std::vector<block_type> blocks;

            constexpr void remove_leading_zero_blocks() {
                auto itr = blocks.end();
                for (; itr != blocks.begin(); --itr) {
                    if (*(itr - 1) != 0) {
                        break;
                    }
                }
                blocks.erase(itr, blocks.end());
            }

        public:
            // blocks is empty if and only if it represents 0.
            var() = default;
            constexpr var(block_type n) {
                if (n != 0) {
                    blocks.push_back(n);
                }
            }

            // Remove leading zeros and reverse the order.
            template <class Iter, class = typename std::iterator_traits<Iter>::iterator_category>
            constexpr var(Iter first, Iter last) {
                while (*first == 0) {
                    if (++first == last) {
                        return;
                    }
                }

                blocks.assign(first, last);
                for (std::size_t idx = 0; idx < blocks.size() / 2; ++idx) {
                    auto const temp = blocks[idx];
                    blocks[idx] = blocks[blocks.size() - idx - 1];
                    blocks[blocks.size() - idx - 1] = temp;
                }
            }

            // Remove leading zeros and reverse the order.
            explicit constexpr var(std::initializer_list<block_type> list)
                : var(list.begin(), list.end()) {}

            constexpr std::size_t number_of_blocks() const noexcept { return blocks.size(); }
            constexpr block_type operator[](std::size_t idx) const { return blocks[idx]; }

            constexpr bool is_zero() const noexcept { return blocks.empty(); }
            constexpr bool is_even() const noexcept {
                if (blocks.empty()) {
                    return true;
                }
                else {
                    return blocks[0] % 2 == 0;
                }
            }

            constexpr bool operator==(var const& n) const noexcept { return blocks == n.blocks; }
            constexpr std::strong_ordering operator<=>(var const& n) const noexcept {
                return detail::comparison_impl(blocks, n.blocks);
            }
            constexpr bool operator==(block_type n) const noexcept {
                return (blocks.size() == 0 && n == 0) || (blocks.size() == 1 && blocks[0] == n);
            }
            constexpr std::strong_ordering operator<=>(block_type n) const noexcept {
                if (blocks.size() > 1) {
                    return std::strong_ordering::greater;
                }
                else if (blocks.size() == 1) {
                    return blocks[0] <=> n;
                }
                else {
                    return 0 <=> n;
                }
            }
            friend constexpr std::strong_ordering operator<=>(block_type x, var const& y) noexcept {
                if (y.blocks.size() > 1) {
                    return std::strong_ordering::less;
                }
                else if (y.blocks.size() == 1) {
                    return x <=> y.blocks[0];
                }
                else {
                    return x <=> 0;
                }
            }

            constexpr var& operator+=(var const& n) {
                if (blocks.size() < n.blocks.size()) {
                    blocks.resize(n.blocks.size(), 0);
                }
                if (detail::add_impl(blocks, n.blocks)) {
                    blocks.push_back(1);
                }
                return *this;
            }
            constexpr var& operator+=(block_type n) {
                if (is_zero()) {
                    blocks.push_back(n);
                    return *this;
                }

                blocks[0] += n;

                // If carry happens,
                if (blocks[0] < n) {
                    // Propagate carry.
                    for (std::size_t idx = 1; idx < blocks.size(); ++idx) {
                        ++blocks[idx];
                        if (blocks[idx] != 0) {
                            return *this;
                        }
                    }
                    blocks.push_back(1);
                }
                return *this;
            }
            template <class T>
            constexpr var operator+(T const& n) const {
                auto r = *this;
                return r += n;
            }
            friend constexpr var operator+(block_type n, var const& m) { return m + n; }
            constexpr var& operator++() {
                *this += 1;
                return *this;
            }
            constexpr var operator++(int) {
                auto temp = *this;
                *this += 1;
                return temp;
            }

            // Precondition: n should be strictly smaller than or equal to the current number.
            constexpr var& operator-=(var const& n) {
                if (n.is_zero()) {
                    return *this;
                }
                detail::subtract_impl(blocks, n.blocks);
                remove_leading_zero_blocks();
                return *this;
            }
            constexpr var operator-(var const& n) const {
                auto r = *this;
                return r -= n;
            }

            // Precondition: *this should be nonzero.
            constexpr var& operator--() {
                util::constexpr_assert<util::error_msgs::underflow>(!is_zero());

                for (std::size_t idx = 0; idx < blocks.size(); ++idx) {
                    --blocks[idx];
                    if (blocks[idx] != std::numeric_limits<block_type>::max()) {
                        break;
                    }
                }
                remove_leading_zero_blocks();
                return *this;
            }

            constexpr var& operator*=(block_type n) {
                if (n == 0) {
                    blocks.clear();
                    return *this;
                }

                block_type carry = 0;
                for (std::size_t idx = 0; idx < blocks.size(); ++idx) {
                    auto mul = wuint::umul128(blocks[idx], n);
                    blocks[idx] = mul.low() + carry;
                    carry = mul.high() + (blocks[idx] < mul.low() ? 1 : 0);
                }
                if (carry != 0) {
                    blocks.push_back(carry);
                }

                return *this;
            }
            constexpr var operator*(block_type n) {
                auto r = *this;
                return r *= n;
            }
            friend constexpr var operator*(var const& x, var const& y) {
                var result;
                if (x.is_zero() || y.is_zero()) {
                    return result;
                }
                result.blocks.resize(x.blocks.size() + y.blocks.size(), 0);

                detail::multiply_impl(x.blocks, y.blocks, result.blocks);
                if (result.blocks.back() == 0) {
                    result.blocks.pop_back();
                }

                return result;
            }
            constexpr var& operator*=(var const& y) {
                auto result = *this * y;
                *this = result;
                return *this;
            }

            // Perform long division
            // *this becomes the remainder, returns the quotient
            // Precondition: n != 0
            constexpr var long_division(var const& n) {
                // The divisor shall not be zero.
                util::constexpr_assert<util::error_msgs::divide_by_zero>(!n.is_zero());

                var quotient;

                // When the divisor has strictly more blocks than the dividend.
                if (n.blocks.size() > blocks.size()) {
                    // Quotient is zero, remainder is the dividend.
                    return quotient;
                }

                quotient.blocks.resize(blocks.size() - n.blocks.size() + 1, 0);
                blocks.resize(detail::long_division_impl(blocks, n.blocks, quotient.blocks));
                quotient.remove_leading_zero_blocks();

                return quotient;
            }

            constexpr var operator/(var const& n) const {
                auto temp = *this;
                return temp.long_division(n);
            }
            constexpr var operator/(block_type n) const {
                auto temp = *this;
                return temp.long_division(n);
            }
            constexpr var operator%(var const& n) const {
                auto temp = *this;
                temp.long_division(n);
                return temp;
            }
            constexpr var operator%(block_type n) const {
                auto temp = *this;
                temp.long_division(n);
                return temp;
            }

            friend constexpr std::size_t bit_width(var const& n) {
                if (n.is_zero()) {
                    return 0;
                }
                return (n.blocks.size() - 1) * number_of_bits_in_block +
                       std::size_t(std::bit_width(n.blocks.back()));
            }

            static constexpr var power_of_2(std::size_t exp) {
                var ret_value;
                ret_value.blocks.resize((exp / number_of_bits_in_block) + 1, 0);
                ret_value.blocks.back() = (block_type(1) << (exp % number_of_bits_in_block));
                return ret_value;
            }

            constexpr bool is_power_of_2() const noexcept {
                if (is_zero()) {
                    return false;
                }
                if (std::bit_floor(blocks.back()) != blocks.back()) {
                    return false;
                }
                for (std::size_t idx = 0; idx < blocks.size() - 1; ++idx) {
                    if (blocks[idx] != 0) {
                        return false;
                    }
                }
                return true;
            }

            constexpr var& operator<<=(std::size_t k) {
                if (k == 0 || is_zero()) {
                    return *this;
                }

                auto const new_leading_one_pos = bit_width(*this) + k;
                auto const prev_number_of_blocks = number_of_blocks();
                blocks.resize((new_leading_one_pos + (number_of_bits_in_block - 1)) /
                              number_of_bits_in_block);

                auto const block_shift = k / number_of_bits_in_block;
                auto const bit_shift = k % number_of_bits_in_block;

                if (bit_shift != 0) {
                    // The first block.
                    if (prev_number_of_blocks + block_shift < blocks.size()) {
                        blocks[prev_number_of_blocks + block_shift] =
                            (blocks[prev_number_of_blocks - 1] >>
                             (number_of_bits_in_block - bit_shift));
                    }
                    // Middle blocks.
                    for (std::size_t idx = prev_number_of_blocks - 1; idx > 0; --idx) {
                        auto const block =
                            ((blocks[idx] << bit_shift) |
                             (blocks[idx - 1] >> (number_of_bits_in_block - bit_shift)));
                        blocks[idx + block_shift] = block;
                    }
                    // The last block.
                    blocks[block_shift] = (blocks[0] << bit_shift);
                }
                else {
                    for (std::size_t idx_p1 = prev_number_of_blocks; idx_p1 > 0; --idx_p1) {
                        blocks[idx_p1 - 1 + block_shift] = blocks[idx_p1 - 1];
                    }
                }

                for (std::size_t idx = 0; idx < block_shift; ++idx) {
                    blocks[idx] = 0;
                }

                return *this;
            }
            constexpr var operator<<(std::size_t k) const {
                auto r = *this;
                return r <<= k;
            }

            constexpr var& operator>>=(std::size_t k) {
                if (k == 0) {
                    return *this;
                }

                auto const prev_leading_one_pos = bit_width(*this);
                if (k >= prev_leading_one_pos) {
                    blocks.clear();
                    return *this;
                }
                auto const new_leading_one_pos = prev_leading_one_pos - k;
                auto const new_number_of_blocks =
                    (new_leading_one_pos + (number_of_bits_in_block - 1)) / number_of_bits_in_block;

                auto const block_shift = k / number_of_bits_in_block;
                auto const bit_shift = k % number_of_bits_in_block;

                if (bit_shift != 0) {
                    // Middle blocks.
                    for (std::size_t idx = 0; idx + block_shift + 1 < number_of_blocks(); ++idx) {
                        blocks[idx] = ((blocks[idx + block_shift] >> bit_shift) |
                                       (blocks[idx + block_shift + 1]
                                        << (number_of_bits_in_block - bit_shift)));
                    }
                    // The last block.
                    blocks[number_of_blocks() - 1 - block_shift] =
                        (blocks[number_of_blocks() - 1] >> bit_shift);
                }
                else {
                    for (std::size_t idx = 0; idx < number_of_blocks() - block_shift; ++idx) {
                        blocks[idx] = blocks[idx + block_shift];
                    }
                }

                blocks.resize(new_number_of_blocks);
                return *this;
            }
            constexpr var operator>>(std::size_t k) const {
                auto r = *this;
                return r >>= k;
            }

            // Find the largest power of 2 dividing *this, divide *this by that power, and return
            // the exponent.
            constexpr std::size_t factor_out_power_of_2() {
                util::constexpr_assert<util::error_msgs::no_error_msg>(!is_zero());

                std::size_t trailing_zero_blocks = 0;
                for (; trailing_zero_blocks < blocks.size(); ++trailing_zero_blocks) {
                    if (blocks[trailing_zero_blocks] != 0) {
                        break;
                    }
                }

                // Shift all blocks.
                for (std::size_t idx = trailing_zero_blocks; idx < blocks.size(); ++idx) {
                    blocks[idx - trailing_zero_blocks] = blocks[idx];
                }
                blocks.resize(blocks.size() - trailing_zero_blocks);

                auto const trailing_zero_bits = std::size_t(std::countr_zero(blocks.front()));
                if (trailing_zero_bits != 0) {
                    for (std::size_t idx = 0; idx < blocks.size() - 1; ++idx) {
                        blocks[idx] =
                            ((blocks[idx] >> trailing_zero_bits) |
                             (blocks[idx + 1] << (number_of_bits_in_block - trailing_zero_bits)));
                    }
                    blocks.back() >>= trailing_zero_bits;

                    if (blocks.back() == 0) {
                        blocks.pop_back();
                    }
                }

                return trailing_zero_blocks * number_of_bits_in_block + trailing_zero_bits;
            }

            // Computes max(floor(log2(x / y)), 0).
            // Precondition: x, y are not zero.
            friend constexpr std::size_t trunc_floor_log2_div(var const& x, var const& y) {
                util::constexpr_assert<util::error_msgs::divide_by_zero>(!x.is_zero() &&
                                                                         !y.is_zero());

                auto const x_leading_one_pos = bit_width(x);
                auto const y_leading_one_pos = bit_width(y);

                if (y_leading_one_pos >= x_leading_one_pos) {
                    return 0;
                }

                auto const total_shift = x_leading_one_pos - y_leading_one_pos;
                auto const block_shift = total_shift / number_of_bits_in_block;
                auto const bit_shift = total_shift % number_of_bits_in_block;

                if (bit_shift != 0) {
                    // The first block of y.
                    if (y.number_of_blocks() + block_shift < x.number_of_blocks()) {
                        if (x[y.number_of_blocks() + block_shift] <
                            (y[y.number_of_blocks() - 1] >>
                             (number_of_bits_in_block - bit_shift))) {
                            return total_shift - 1;
                        }
                    }
                    // Middle blocks of y.
                    for (std::size_t idx = y.number_of_blocks() - 1; idx > 0; --idx) {
                        auto const y_block =
                            ((y[idx] << bit_shift) |
                             (y[idx - 1] >> (number_of_bits_in_block - bit_shift)));
                        if (x[idx + block_shift] < y_block) {
                            return total_shift - 1;
                        }
                    }
                    // The last block of y.
                    if (x[block_shift] < (y[0] << bit_shift)) {
                        return total_shift - 1;
                    }
                }
                else {
                    for (std::size_t idx_p1 = y.number_of_blocks(); idx_p1 > 0; --idx_p1) {
                        if (x[idx_p1 - 1 + block_shift] < y[idx_p1 - 1]) {
                            return total_shift - 1;
                        }
                    }
                }
                return total_shift;
            }
        };

        struct var_div_t {
            var quot;
            var rem;
        };
        inline constexpr var_div_t div(var const& x, var const& y) {
            var_div_t ret;
            ret.rem = x;
            ret.quot = ret.rem.long_division(y);
            return ret;
        }
        inline constexpr var_div_t div(var&& x, var const& y) {
            var_div_t ret;
            ret.rem = static_cast<var&&>(x);
            ret.quot = ret.rem.long_division(y);
            return ret;
        }

        inline constexpr var div_ceil(var const& x, var const& y) {
            auto dividend = x;
            auto quotient = dividend.long_division(y);
            if (!dividend.is_zero()) {
                ++quotient;
            }
            return quotient;
        }
        inline constexpr var div_ceil(var&& x, var const& y) {
            auto quotient = x.long_division(y);
            if (!x.is_zero()) {
                ++quotient;
            }
            return quotient;
        }
    }
}

#endif