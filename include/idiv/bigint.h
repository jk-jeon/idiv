// Copyright 2023-2024 Junekey Jeon
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
#include <concepts>
#include <limits>
#include <vector>

namespace jkj {
    namespace bigint {
        // Big integers will be represented as a list of 64-bit blocks.
        // std::uint_least64_t is used to represent each block.
        // Even if std::uint_least64_t is not exactly of 64-bit, we implictly assume that many of
        // functions in this header is implicitly assuming that each block does not contain anything
        // larger than 64-bits.
        using block_type = std::uint_least64_t;
        using signed_block_type = std::int_least64_t;
        static constexpr std::size_t number_of_bits_in_block = 64;
        static constexpr block_type largest_pow10_in_block = UINT64_C(1000'0000'0000'0000'0000);
        static_assert(std::numeric_limits<block_type>::max() >= largest_pow10_in_block);

        // Represents sign.
        using sign_t = util::sign_t;

        constexpr sign_t invert_sign(sign_t s) noexcept {
            return s == sign_t::positive ? sign_t::negative : sign_t::positive;
        }

        // Concepts used for operator overloading.
        // If std::uint_fast64_t has more than 64 digits, then block_type itself ironically doesn't
        // satisfy this concept. However, this saves a lot of headaches so let us just do it like
        // this.
        template <class T>
        concept convertible_to_block_type =
            std::unsigned_integral<T> &&
            (std::numeric_limits<T>::digits <= number_of_bits_in_block);
        template <class T>
        concept convertible_to_signed_block_type =
            std::signed_integral<T> && (std::numeric_limits<T>::digits <= number_of_bits_in_block);

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

            template <class X, class Y>
            constexpr bool equal_impl(X const& x, Y const& y) noexcept {
                if (x.size() != y.size()) {
                    return false;
                }

                for (std::size_t idx = 0; idx < x.size(); ++idx) {
                    if (x[idx] != y[idx]) {
                        return false;
                    }
                }
                return true;
            }

            // Precondition: x.size() >= y.size(), either x and y point at the same position or have
            // no overlap.
            // Modifies x in-place. Returns true if a leading block with the carry bit
            // should be appended.
            template <class X, class Y>
            constexpr decltype(auto) add_impl(X&& x, Y const& y) noexcept {
                util::constexpr_assert(x.size() >= y.size());

                unsigned int carry = 0;
                for (std::size_t idx = 0; idx < y.size(); ++idx) {
                    x[idx] = wuint::add_carry64(x[idx], y[idx], carry);
                }

                if (carry != 0) {
                    std::size_t idx = y.size();
                    while (true) {
                        // Carry.
                        if (idx == x.size()) {
                            return true;
                        }

                        if (x[idx] == wuint::uint64_mask) {
                            x[idx] = 0;
                            ++idx;
                        }
                        else {
                            ++x[idx];
                            break;
                        }
                    }
                }
                return false;
            }

            // Precondition: x.size() >= y.size(), x >= y, either x and y point at the same position
            // or have no overlap.
            // Modifies x in-place. x may contain leading zeros after the
            // function has executed.
            template <class X, class Y>
            constexpr void subtract_impl(X&& x, Y const& y) noexcept {
                util::constexpr_assert<util::error_msgs::underflow>(x.size() >= y.size());

                unsigned int borrow = 0;
                for (std::size_t idx = 0; idx < y.size(); ++idx) {
                    x[idx] = wuint::sub_borrow64(x[idx], y[idx], borrow);
                }

                if (borrow != 0) {
                    util::constexpr_assert<util::error_msgs::underflow>(x.size() > y.size());
                    std::size_t idx = y.size();
                    while (x[idx] == 0) {
                        x[idx] = wuint::uint64_mask;
                        ++idx;
                        util::constexpr_assert<util::error_msgs::underflow>(idx != x.size());
                    }
                    --x[idx];
                }
            }

            // Precondition: result is zero-initialized, and has the size x.size() + y.size().
            // result should not alias with x, y.
            template <class X, class Y, class Result>
            constexpr void multiply_impl(X const& x, Y const& y, Result& result) noexcept {
                for (std::size_t y_idx = 0; y_idx < y.size(); ++y_idx) {
                    // Compute y[y_idx] * x and accumulate it into the result
                    for (std::size_t x_idx = 0; x_idx < x.size(); ++x_idx) {
                        auto mul = wuint::umul128(x[x_idx], y[y_idx]);

                        // Add the first half
                        unsigned int carry = 0;
                        result[x_idx + y_idx] =
                            wuint::add_carry64(result[x_idx + y_idx], mul.low(), carry);

                        // Add the second half
                        result[x_idx + y_idx + 1] =
                            wuint::add_carry64(result[x_idx + y_idx + 1], mul.high(), carry);

                        // If there is carry,
                        if (carry != 0) {
                            // Propagate.
                            auto idx = x_idx + y_idx + 2;
                            while (result[idx] == wuint::uint64_mask) {
                                result[idx] = 0;
                                ++idx;
                                util::constexpr_assert(idx < result.size());
                            }
                            ++result[idx];
                        }
                    }
                }
            }

            // - x: dividend. Should not have any leading zero blocks. Modified in-place and becomes
            // the remainder after the function is executed.
            // - y: divisor. Should not have any leading zero blocks.
            // - q: quotient. Should be zero-initialized and has the size at least x.size() -
            // y.size() + 1.
            // All of x, y, q should not alias with each other.
            // Returns the number of blocks in the remainder.
            template <class X, class Y, class Quotient>
            constexpr std::size_t long_division_impl(X&& x, Y const& y, Quotient& q) noexcept {
                util::constexpr_assert(q.size() + y.size() >= x.size() + 1);

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
                        unsigned int borrow;
                        {
                            auto const y_block = (y[0] << bit_shift);
                            borrow = x[block_shift] < y_block ? 1 : 0;
                            x[block_shift] -= y_block;
                        }

                        // Middle blocks of y.
                        for (std::size_t idx = 1; idx < y.size(); ++idx) {
                            auto with_borrow = x[idx + block_shift] - borrow;
                            unsigned int first_borrow =
                                (with_borrow > x[idx + block_shift]) ? 1 : 0;

                            auto const y_block =
                                (y[idx] << bit_shift) |
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
                            util::constexpr_assert(with_borrow >= y_block);
                            x[y.size() + block_shift] = with_borrow - y_block;
                            borrow =
                                first_borrow | ((x[y.size() + block_shift] > with_borrow) ? 1 : 0);
                        }
                        util::constexpr_assert(borrow == 0);
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
                            util::constexpr_assert(y.size() + block_shift < x_size);
                            util::constexpr_assert(x[y.size() + block_shift] != 0);
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

        // Static array of blocks.
        // The first block is the least significant block and the last block is the most
        // significant block.
        template <std::size_t N>
        using static_block_holder = util::array<block_type, N>;

        // Forward declarations.
        class uint_view;
        class int_view;
        namespace detail {
            template <auto arr>
            struct uint_const_impl;

            template <sign_t sign_, auto arr>
            struct int_const_impl;
        }
        class uint_var;
        class int_var;

        // std::span-like view type over uint_const_t / uint_var.
        class uint_view {
            util::span<block_type const> blocks_;

            explicit constexpr uint_view(util::span<block_type const> blocks) noexcept
                : blocks_{blocks} {}

        public:
            friend class int_view;
            template <auto arr>
            friend struct detail::uint_const_impl;
            template <sign_t sign_, auto arr>
            friend struct detail::int_const_impl;
            friend class uint_var;
            friend class int_var;

            constexpr auto blocks() const noexcept { return blocks_; }
            constexpr std::size_t number_of_blocks() const noexcept { return blocks_.size(); }
            constexpr block_type operator[](std::size_t idx) const { return blocks_[idx]; }

            constexpr bool is_zero() const noexcept { return blocks_.empty(); }
            constexpr bool is_even() const noexcept {
                if (blocks_.empty()) {
                    return true;
                }
                else {
                    return blocks_[0] % 2 == 0;
                }
            }

            static constexpr uint_view make_view_from_single_block(block_type const& n) noexcept {
                util::constexpr_assert<util::error_msgs::overflow>(n <= wuint::uint64_mask);
                return uint_view{{&n, 1}};
            }
        };

        // Signed version of uint_view.
        // Internally uses sign-magnitude representation, but does not allow negative zero.
        class int_view {
            util::span<block_type const> blocks_;
            sign_t sign_;

            explicit constexpr int_view(sign_t sign, util::span<block_type const> blocks) noexcept
                : blocks_{blocks}, sign_{sign} {}

        public:
            friend class uint_view;
            template <auto arr>
            friend struct detail::uint_const_impl;
            template <sign_t sign_, auto arr>
            friend struct detail::int_const_impl;
            friend class uint_var;
            friend class int_var;

            friend constexpr int_view to_signed(uint_view n) noexcept;
            friend constexpr int_view to_negative(uint_view n) noexcept;

            constexpr auto blocks() const noexcept { return blocks_; }
            constexpr std::size_t number_of_blocks() const noexcept { return blocks_.size(); }
            constexpr block_type operator[](std::size_t idx) const { return blocks_[idx]; }

            constexpr bool is_zero() const noexcept { return blocks_.empty(); }
            constexpr bool is_even() const noexcept { return uint_view{blocks_}.is_even(); }

            constexpr sign_t sign() const noexcept { return sign_; }
            constexpr bool is_strictly_positive() const noexcept {
                return !is_zero() && sign_ == sign_t::positive;
            }
            constexpr bool is_strictly_negative() const noexcept {
                return sign_ == sign_t::negative;
            }
            constexpr bool is_nonnegative() const noexcept { return sign_ == sign_t::positive; }
            constexpr bool is_nonpositive() const noexcept {
                return is_zero() || sign_ == sign_t::negative;
            }

            // Absolute value.
            constexpr uint_view abs() const noexcept { return uint_view{blocks_}; }
        };

        // Implememtation details for bigint::uint_const: big unsigned integer constants whose
        // values are encoded in their types.
        namespace detail {
            // Take the first N blocks and discard the remaining.
            template <std::size_t N, class Range>
            constexpr static_block_holder<N> slice(Range const& arr) noexcept {
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

            // Remove leading zero blocks from arr.
            template <auto arr>
            constexpr auto remove_leading_zero_blocks() noexcept {
                return slice<count_blocks_excluding_leading_zeros(arr)>(arr);
            }

            // Reverse the order of the array and check all of elements are within 64-bits.
            template <std::size_t N>
            constexpr static_block_holder<N>
            reverse_and_validate(static_block_holder<N> const& arr) noexcept {
                if constexpr (N == 0) {
                    return {};
                }
                else {
                    static_block_holder<N> result = arr;
                    for (std::size_t idx = 0; idx < N / 2; ++idx) {
                        auto const temp = result[idx];
                        result[idx] = result[N - idx - 1];
                        result[N - idx - 1] = temp;

                        util::constexpr_assert(result[idx] <= wuint::uint64_mask &&
                                               result[N - idx - 1] <= wuint::uint64_mask);
                    }
                    return result;
                }
            }

            template <auto arr>
            struct uint_const_impl;

            // Does not allow leading zeros.
            template <std::size_t N, static_block_holder<N> arr>
            struct uint_const_impl<arr> {
                static_assert(count_blocks_excluding_leading_zeros<N>(arr) == N);
                using block_holder_type = static_block_holder<N>;
                static constexpr block_holder_type const& blocks() noexcept { return blocks_; }
                static constexpr std::size_t number_of_blocks() noexcept { return N; }

                template <class = void>
                    requires(N <= 1)
                constexpr explicit operator block_type() const noexcept {
                    if constexpr (N == 0) {
                        return 0;
                    }
                    else {
                        return arr[0];
                    }
                }

                static constexpr bool is_zero() noexcept { return N == 0; }
                static constexpr bool is_even() noexcept {
                    if constexpr (N == 0) {
                        return true;
                    }
                    else {
                        return blocks()[0] % 2 == 0;
                    }
                }

                template <class = void>
                    requires(N != 0)
                constexpr operator uint_view() const noexcept {
                    return uint_view{{blocks().data(), blocks().size()}};
                }

            private:
                static constexpr block_holder_type blocks_ = arr;
            };

            template <std::size_t N, static_block_holder<N> arr>
                requires(N != 0)
            constexpr uint_view to_view(uint_const_impl<arr> x) noexcept {
                return uint_view(x);
            }
        }

        // The first parameter is the most significant block.
        // E.g., uint_const_t<3, 7, 11> represents 3 * 2^128 + 7 * 2^64 + 11.
        // All leading zero blocks will be automatically removed, so for e.g., uint_const_t<0, 0, 1>
        // is the same type as uint_const_t<1>.
        template <block_type... blocks>
        using uint_const_t =
            detail::uint_const_impl<detail::remove_leading_zero_blocks<detail::reverse_and_validate(
                static_block_holder<sizeof...(blocks)>{blocks...})>()>;

        namespace detail {
            template <class T>
            struct is_uint_const_impl : std::false_type {};

            template <std::size_t N, static_block_holder<N> arr>
            struct is_uint_const_impl<uint_const_impl<arr>> : std::true_type {};
        }

        template <class T>
        concept uint_const = detail::is_uint_const_impl<T>::value;


        // Implememtation details for bigint::int_const: big signed integer constants whose
        // values are encoded in their types.
        namespace detail {
            template <sign_t sign, auto arr>
            struct int_const_impl;

            // Internally uses sign-magnitude representation, but does not allow negative zero.
            template <sign_t sign_, std::size_t N, static_block_holder<N> arr>
            struct int_const_impl<sign_, arr> {
                static_assert(N > 0 || sign_ == sign_t::positive);
                static_assert(count_blocks_excluding_leading_zeros<N>(arr) == N);

                using block_holder_type = static_block_holder<N>;
                static constexpr block_holder_type const& blocks() noexcept { return blocks_; }
                static constexpr std::size_t number_of_blocks() noexcept { return N; }

            private:
                static constexpr bool is_overflow_safe() noexcept {
                    // Recall that signed-to-unsigned conversion is defined, and the result should
                    // be as if bit-pattern is reinterpreted from a 2's complement signed integer
                    // to the unsigned integer of the same bit-width.
                    if constexpr (N == 0) {
                        return true;
                    }
                    else if constexpr (N > 1) {
                        return false;
                    }
                    else if constexpr (sign_ == sign_t::positive) {
                        return arr[0] <= static_cast<block_type>(
                                             std::numeric_limits<signed_block_type>::max());
                    }
                    else {
                        return (block_type(0) - arr[0]) >=
                               static_cast<block_type>(
                                   std::numeric_limits<signed_block_type>::min());
                    }
                }

            public:
                template <class = void>
                    requires(is_overflow_safe())
                constexpr explicit operator signed_block_type() const noexcept {
                    if constexpr (N == 0) {
                        return 0;
                    }
                    else {
                        return arr[0];
                    }
                }

                static constexpr bool is_zero() noexcept { return N == 0; }
                static constexpr bool is_even() noexcept {
                    if constexpr (N == 0) {
                        return true;
                    }
                    else {
                        return blocks()[0] % 2 == 0;
                    }
                }

                static constexpr sign_t sign() noexcept { return sign_; }
                static constexpr bool is_strictly_positive() noexcept {
                    return !is_zero() && sign_ == sign_t::positive;
                }
                static constexpr bool is_strictly_negative() noexcept {
                    return sign_ == sign_t::negative;
                }
                static constexpr bool is_nonnegative() noexcept {
                    return sign_ == sign_t::positive;
                }
                static constexpr bool is_nonpositive() noexcept {
                    return is_zero() || sign_ == sign_t::negative;
                }

                // Absolute value.
                static constexpr uint_const_impl<arr> abs() noexcept { return {}; }

                template <class = void>
                    requires(N != 0)
                constexpr operator int_view() const noexcept {
                    return int_view{sign_, {blocks().data(), blocks().size()}};
                }

            private:
                static constexpr block_holder_type blocks_ = arr;
            };

            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
                requires(N != 0)
            constexpr int_view to_view(int_const_impl<sign, arr> x) noexcept {
                return int_view(x);
            }
        }

        using detail::to_view;

        // The first parameter is the sign and the second parameter is the most significant block.
        // E.g., int_const_t<sign_t::negative, 3, 7, 11> represents -(3 * 2^128 + 7 * 2^64 + 11).
        // All leading zero blocks will be automatically removed, so for e.g.,
        // int_const_t<sign_t::negative, 0, 0, 1> is the same type as int_const_t<sign_t::negative,
        // 1>, and the positive zero and the negative zero yield the same type.
        template <sign_t sign = sign_t::positive, block_type... blocks>
        using int_const_t =
            detail::int_const_impl<sizeof...(blocks) == 0 ? sign_t::positive : sign,
                                   detail::remove_leading_zero_blocks<detail::reverse_and_validate(
                                       static_block_holder<sizeof...(blocks)>{blocks...})>()>;

        namespace detail {
            template <class T>
            struct is_int_const_impl : std::false_type {};

            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            struct is_int_const_impl<int_const_impl<sign, arr>> : std::true_type {};
        }

        template <class T>
        concept int_const = detail::is_int_const_impl<T>::value;

        // std::bit_width counterpart for uint_view.
        constexpr std::size_t bit_width(uint_view n) noexcept {
            if (n.is_zero()) {
                return 0;
            }
            return (n.blocks().size() - 1) * number_of_bits_in_block +
                   std::size_t(std::bit_width(n.blocks().back()));
        }

        // Comparison operators.

        constexpr bool operator==(uint_view x, uint_view y) noexcept {
            return detail::equal_impl(x.blocks(), y.blocks());
        }
        constexpr bool operator==(uint_view x, int_view y) noexcept {
            if (y.is_nonnegative()) {
                return detail::equal_impl(x.blocks(), y.blocks());
            }
            return false;
        }
        constexpr bool operator==(uint_view x, convertible_to_block_type auto y) noexcept {
            if (x.is_zero()) {
                return y == 0;
            }
            else if (x.number_of_blocks() == 1) {
                return x[0] == y;
            }
            else {
                return false;
            }
        }
        constexpr bool operator==(uint_view x, convertible_to_signed_block_type auto y) noexcept {
            if (y >= 0) {
                return x == static_cast<block_type>(y);
            }
            return false;
        }
        constexpr bool operator==(uint_view x, uint_const_t<>) noexcept { return x.is_zero(); }
        constexpr bool operator==(uint_view x, int_const_t<>) noexcept { return x.is_zero(); }
        constexpr bool operator==(int_view x, int_view y) noexcept {
            if (x.is_nonnegative()) {
                if (y.is_nonnegative()) {
                    return x.abs() == y.abs();
                }
                return false;
            }
            return y.is_strictly_negative() && x.abs() == y.abs();
        }
        constexpr bool operator==(int_view x, convertible_to_block_type auto y) noexcept {
            if (x.is_nonnegative()) {
                return x.abs() == y;
            }
            return false;
        }
        constexpr bool operator==(int_view x, convertible_to_signed_block_type auto y) noexcept {
            // Signed-to-unsigned conversion is well-defined.
            if (x.is_nonnegative()) {
                if (y >= 0) {
                    return x.abs() == static_cast<block_type>(y);
                }
                return false;
            }
            return y < 0 && (block_type(0) - x[0]) == static_cast<block_type>(y);
        }
        constexpr bool operator==(int_view x, uint_const_t<>) noexcept { return x.is_zero(); }
        constexpr bool operator==(int_view x, int_const_t<>) noexcept { return x.is_zero(); }
        constexpr bool operator==(uint_const_t<>, convertible_to_block_type auto y) noexcept {
            return y == 0;
        }
        constexpr bool operator==(uint_const_t<>,
                                  convertible_to_signed_block_type auto y) noexcept {
            return y == 0;
        }
        constexpr bool operator==(uint_const_t<>, uint_const_t<>) noexcept { return true; }
        constexpr bool operator==(uint_const_t<>, int_const_t<>) noexcept { return true; }
        constexpr bool operator==(int_const_t<>, convertible_to_block_type auto y) noexcept {
            return y == 0;
        }
        constexpr bool operator==(int_const_t<>, convertible_to_signed_block_type auto y) noexcept {
            return y == 0;
        }


        constexpr std::strong_ordering operator<=>(uint_view x, uint_view y) noexcept {
            return detail::comparison_impl(x.blocks(), y.blocks());
        }
        constexpr std::strong_ordering operator<=>(uint_view x, int_view y) noexcept {
            if (y.is_nonnegative()) {
                return detail::comparison_impl(x.blocks(), y.blocks());
            }
            return std::strong_ordering::greater;
        }
        constexpr std::strong_ordering operator<=>(uint_view x,
                                                   convertible_to_block_type auto y) noexcept {
            if (x.number_of_blocks() == 0) {
                return y == 0 ? std::strong_ordering::equal : std::strong_ordering::less;
            }
            else if (x.number_of_blocks() > 1) {
                return std::strong_ordering::greater;
            }
            else {
                return x[0] <=> y;
            }
        }
        constexpr std::strong_ordering
        operator<=>(uint_view x, convertible_to_signed_block_type auto y) noexcept {
            if (y >= 0) {
                return x <=> static_cast<block_type>(y);
            }
            return std::strong_ordering::greater;
        }
        constexpr std::strong_ordering operator<=>(uint_view x, uint_const_t<>) noexcept {
            return x.number_of_blocks() == 0 ? std::strong_ordering::equal
                                             : std::strong_ordering::greater;
        }
        constexpr std::strong_ordering operator<=>(uint_view x, int_const_t<>) noexcept {
            return x <=> uint_const_t<>{};
        }
        constexpr std::strong_ordering operator<=>(int_view x, int_view y) noexcept {
            if (x.is_nonnegative()) {
                if (y.is_nonnegative()) {
                    return x.abs() <=> y.abs();
                }
                return std::strong_ordering::greater;
            }
            if (y.is_nonnegative()) {
                return std::strong_ordering::less;
            }
            return y.abs() <=> x.abs();
        }
        constexpr std::strong_ordering operator<=>(int_view x,
                                                   convertible_to_block_type auto y) noexcept {
            if (x.is_nonnegative()) {
                return x.abs() <=> y;
            }
            return std::strong_ordering::less;
        }
        constexpr std::strong_ordering
        operator<=>(int_view x, convertible_to_signed_block_type auto y) noexcept {
            if (x.is_nonnegative()) {
                if (y >= 0) {
                    return x.abs() <=> static_cast<block_type>(y);
                }
                return std::strong_ordering::greater;
            }
            if (y >= 0) {
                return std::strong_ordering::less;
            }
            return ((block_type(0) - static_cast<block_type>(y)) & wuint::uint64_mask) <=> x.abs();
        }
        constexpr std::strong_ordering operator<=>(int_view x, uint_const_t<>) noexcept {
            if (x.is_strictly_positive()) {
                return std::strong_ordering::greater;
            }
            else if (x.is_strictly_negative()) {
                return std::strong_ordering::less;
            }
            return std::strong_ordering::equal;
        }
        constexpr std::strong_ordering operator<=>(int_view x, int_const_t<>) noexcept {
            return x <=> uint_const_t<>{};
        }
        constexpr std::strong_ordering operator<=>(uint_const_t<>,
                                                   convertible_to_block_type auto y) noexcept {
            return 0 <=> y;
        }
        constexpr std::strong_ordering
        operator<=>(uint_const_t<>, convertible_to_signed_block_type auto y) noexcept {
            return 0 <=> y;
        }
        constexpr std::strong_ordering operator<=>(uint_const_t<>, uint_const_t<>) noexcept {
            return std::strong_ordering::equal;
        }
        constexpr std::strong_ordering operator<=>(uint_const_t<>, int_const_t<>) noexcept {
            return std::strong_ordering::equal;
        }
        constexpr std::strong_ordering operator<=>(int_const_t<>,
                                                   convertible_to_block_type auto y) noexcept {
            return 0 <=> y;
        }
        constexpr std::strong_ordering
        operator<=>(int_const_t<>, convertible_to_signed_block_type auto y) noexcept {
            return 0 <=> y;
        }


        // Dynamically-sized big unsigned integers.
        // Internally stores a list of blocks of bits, where the least significant block comes
        // first. The stored list having no leading zero is a class invariant. The entire class is
        // mostly strongly exception-safe:
        // - every const member function is noexcept and does not actually change anything,
        // - in any non-const member function, the only source of exceptions is memory allocation
        // faillure, and if it throws, the previous state of the class is retained.
        // - Some binary operations involving rvalue-reference parameters can modify some of those
        // parameters, with or without the presene of exceptions. But they always leave them in a
        // well-defined state.
        // - Any operations should work just fine even if operands alias to each other.
        //
        // Note that, however, violation of preconditions may terminate (through assert) or lead to
        // undefined behavior (if assert does not terminate).
        class uint_var {
            std::vector<block_type> blocks_;

            constexpr void remove_leading_zero_blocks() noexcept {
                auto itr = blocks_.end();
                for (; itr != blocks_.begin(); --itr) {
                    if (*(itr - 1) != 0) {
                        break;
                    }
                }
                blocks_.erase(itr, blocks_.end());
            }

        public:
            friend class int_var;

            // blocks_ is empty if and only if it represents 0.
            uint_var() = default;
            constexpr uint_var(block_type n) {
                util::constexpr_assert<util::error_msgs::overflow>(n <= wuint::uint64_mask);
                if (n != 0) {
                    blocks_.resize(1, n);
                }
            }

            // Remove leading zeros and reverse the order.
            template <class Iter, class = typename std::iterator_traits<Iter>::iterator_category>
            explicit constexpr uint_var(Iter first, Iter last) {
                while (*first == 0) {
                    if (++first == last) {
                        return;
                    }
                }

                blocks_.assign(first, last);
                for (std::size_t idx = 0; idx < blocks_.size() / 2; ++idx) {
                    auto const temp = blocks_[idx];
                    blocks_[idx] = blocks_[blocks_.size() - idx - 1];
                    blocks_[blocks_.size() - idx - 1] = temp;
                    util::constexpr_assert<util::error_msgs::overflow>(
                        blocks_[idx] <= wuint::uint64_mask &&
                        blocks_[blocks_.size() - idx - 1] <= wuint::uint64_mask);
                }
            }

            // Remove leading zeros and reverse the order.
            explicit constexpr uint_var(std::initializer_list<block_type> list)
                : uint_var(list.begin(), list.end()) {}

            // Copy the blocks.
            constexpr uint_var(uint_view n) {
                blocks_.assign(n.blocks().cbegin(), n.blocks().cend());
            }
            constexpr uint_var(uint_const_t<>) {}

            constexpr util::span<block_type const> blocks() const noexcept {
                return {blocks_.data(), blocks_.size()};
            }
            constexpr std::size_t number_of_blocks() const noexcept { return blocks_.size(); }
            constexpr block_type operator[](std::size_t idx) const noexcept { return blocks_[idx]; }

            // Obtain view.
            constexpr operator uint_view() const noexcept { return uint_view{blocks()}; }

            constexpr bool is_zero() const noexcept { return blocks_.empty(); }
            constexpr bool is_even() const noexcept { return uint_view(*this).is_even(); }

            constexpr uint_var& operator+=(uint_view n) {
                if (n.is_zero()) {
                    return *this;
                }
                else if (blocks_.size() < n.blocks().size()) {
                    blocks_.reserve(n.blocks().size() + 1);
                    blocks_.resize(n.blocks().size(), 0);
                }
                else {
                    blocks_.reserve(blocks_.size() + 1);
                }

                if (detail::add_impl(blocks_, n.blocks())) {
                    blocks_.push_back(1);
                }
                return *this;
            }
            constexpr uint_var& operator+=(block_type n) {
                util::constexpr_assert<util::error_msgs::overflow>(n <= wuint::uint64_mask);
                if (is_zero()) {
                    if (n != 0) {
                        blocks_.resize(1, n);
                    }
                    return *this;
                }
                blocks_.reserve(blocks_.size() + 1);

                blocks_[0] += n;

                // If carry happens,
                if (blocks_[0] < n) {
                    // Propagate carry.
                    for (std::size_t idx = 1; idx < blocks_.size(); ++idx) {
                        ++blocks_[idx];
                        if (blocks_[idx] != 0) {
                            return *this;
                        }
                    }
                    blocks_.push_back(1);
                }
                return *this;
            }
            constexpr uint_var& operator+=(uint_const_t<>) { return *this; }

            constexpr uint_var& operator++() {
                if (is_zero()) {
                    blocks_.resize(1, 1);
                    return *this;
                }
                blocks_.reserve(blocks_.size() + 1);

                std::size_t idx = 0;
                while (blocks_[idx] == wuint::uint64_mask) {
                    blocks_[idx] = 0;
                    ++idx;

                    if (idx == number_of_blocks()) {
                        blocks_.push_back(1);
                        return *this;
                    }
                }
                ++blocks_[idx];
                return *this;
            }
            constexpr uint_var operator++(int) & {
                auto temp = *this;
                ++*this;
                return temp;
            }
            constexpr uint_var operator++(int) && {
                auto temp = std::move(*this);
                ++*this;
                return temp;
            }

            // Precondition: n should be strictly smaller than or equal to the current number.
            constexpr uint_var& operator-=(uint_view n) noexcept {
                if (n.is_zero()) {
                    return *this;
                }
                detail::subtract_impl(blocks_, n.blocks());
                remove_leading_zero_blocks();
                return *this;
            }
            constexpr uint_var& operator-=(block_type n) noexcept {
                util::constexpr_assert<util::error_msgs::overflow>(n <= wuint::uint64_mask);
                jkj::util::constexpr_assert<jkj::util::error_msgs::underflow>(!is_zero());

                if (blocks_[0] >= n) {
                    blocks_[0] -= n;
                    if (blocks_.size() == 1 && blocks_[0] == 0) {
                        blocks_.clear();
                    }
                }
                else {
                    jkj::util::constexpr_assert<jkj::util::error_msgs::underflow>(
                        number_of_blocks() > 1);

                    blocks_[0] -= n;
                    std::size_t idx = 1;
                    while (blocks_[idx] == 0) {
                        blocks_[idx] = std::numeric_limits<block_type>::max();
                        ++idx;
                    }
                    jkj::util::constexpr_assert<jkj::util::error_msgs::underflow>(
                        number_of_blocks() > idx);
                    --blocks_[idx];
                    if (idx == number_of_blocks() - 1 && blocks_[idx] == 0) {
                        blocks_.pop_back();
                    }
                }

                return *this;
            }
            constexpr uint_var& operator-=(uint_const_t<>) noexcept { return *this; }

            // Precondition: *this should be nonzero.
            constexpr uint_var& operator--() noexcept {
                util::constexpr_assert<util::error_msgs::underflow>(!is_zero());

                std::size_t idx = 0;
                while (blocks_[idx] == 0) {
                    blocks_[idx] = wuint::uint64_mask;
                    ++idx;
                }
                --blocks_[idx];
                if (idx == number_of_blocks() - 1 && blocks_[idx] == 0) {
                    blocks_.pop_back();
                }

                return *this;
            }
            constexpr uint_var operator--(int) & {
                auto temp = *this;
                --*this;
                return temp;
            }
            constexpr uint_var operator--(int) && {
                auto temp = std::move(*this);
                --*this;
                return temp;
            }

            constexpr uint_var& operator*=(uint_const_t<>) noexcept {
                blocks_.clear();
                return *this;
            }
            constexpr uint_var& operator*=(uint_const_t<1u>) noexcept { return *this; }
            constexpr uint_var& operator*=(block_type n) {
                util::constexpr_assert<util::error_msgs::overflow>(n <= wuint::uint64_mask);
                if (n == 0) {
                    blocks_.clear();
                    return *this;
                }
                blocks_.reserve(blocks_.size() + 1);

                block_type high_bits = 0;
                for (std::size_t idx = 0; idx < blocks_.size(); ++idx) {
                    auto mul = wuint::umul128(blocks_[idx], n);
                    unsigned int carry = 0;
                    blocks_[idx] = wuint::add_carry64(mul.low(), high_bits, carry);
                    high_bits = wuint::add_carry64(mul.high(), 0, carry);
                }
                if (high_bits != 0) {
                    blocks_.push_back(high_bits);
                }

                return *this;
            }

            friend constexpr uint_var operator*(uint_view x, uint_view y);

            // Perform long division.
            // *this becomes the remainder, returns the quotient.
            // Precondition: n != 0
            constexpr uint_var long_division(uint_view n) {
                // The divisor shall not be zero.
                util::constexpr_assert<util::error_msgs::divide_by_zero>(!n.is_zero());

                uint_var quotient;

                // Avoid aliasing issue. There are only handful of places where uint_view can get
                // constructed, so this should be enough.
                if (blocks_.data() == n.blocks().data()) {
                    // *this and n are the same, so the quotient is 1 and the remainder is 0.
                    quotient.blocks_.push_back(1);
                    blocks_.clear();
                    return quotient;
                }
                // When the divisor has strictly more blocks than the dividend.
                else if (n.blocks().size() > blocks_.size()) {
                    // Quotient is zero, remainder is the dividend.
                    return quotient;
                }

                quotient.blocks_.resize(blocks_.size() - n.blocks().size() + 1, 0);
                blocks_.resize(detail::long_division_impl(blocks_, n.blocks(), quotient.blocks_));
                quotient.remove_leading_zero_blocks();

                return quotient;
            }
            constexpr uint_var long_division(uint_const_t<1u>) {
                auto quotient = std::move(*this);
                util::constexpr_assert(is_zero());
                return quotient;
            }

            constexpr uint_var& operator%=(uint_view y) {
                long_division(y);
                return *this;
            }
            constexpr uint_var& operator%=(block_type y) {
                long_division(uint_view::make_view_from_single_block(y));
                return *this;
            }
            constexpr uint_var& operator%=(uint_const_t<1u>) noexcept {
                blocks_.clear();
                return *this;
            }
            constexpr uint_var& operator%=(uint_const_t<>) const = delete;

            static constexpr uint_var power_of_2(std::size_t exp) {
                uint_var ret_value;
                ret_value.blocks_.resize((exp / number_of_bits_in_block) + 1, 0);
                ret_value.blocks_.back() = (block_type(1) << (exp % number_of_bits_in_block));
                return ret_value;
            }

            constexpr bool is_power_of_2() const noexcept {
                if (is_zero()) {
                    return false;
                }
                if (std::bit_floor(blocks_.back()) != blocks_.back()) {
                    return false;
                }
                for (std::size_t idx = 0; idx < blocks_.size() - 1; ++idx) {
                    if (blocks_[idx] != 0) {
                        return false;
                    }
                }
                return true;
            }

            constexpr uint_var& operator<<=(std::size_t k) {
                if (k == 0 || is_zero()) {
                    return *this;
                }

                auto const new_leading_one_pos = bit_width(*this) + k;
                auto const prev_number_of_blocks = number_of_blocks();
                blocks_.resize((new_leading_one_pos + (number_of_bits_in_block - 1)) /
                               number_of_bits_in_block);

                auto const block_shift = k / number_of_bits_in_block;
                auto const bit_shift = k % number_of_bits_in_block;

                if (bit_shift != 0) {
                    // The first block.
                    if (prev_number_of_blocks + block_shift < blocks_.size()) {
                        blocks_[prev_number_of_blocks + block_shift] =
                            (blocks_[prev_number_of_blocks - 1] >>
                             (number_of_bits_in_block - bit_shift));
                    }
                    // Middle blocks.
                    for (std::size_t idx = prev_number_of_blocks - 1; idx > 0; --idx) {
                        auto const block =
                            ((blocks_[idx] << bit_shift) |
                             (blocks_[idx - 1] >> (number_of_bits_in_block - bit_shift))) &
                            wuint::uint64_mask;
                        blocks_[idx + block_shift] = block;
                    }
                    // The last block.
                    blocks_[block_shift] = (blocks_[0] << bit_shift) & wuint::uint64_mask;
                }
                else {
                    for (std::size_t idx_p1 = prev_number_of_blocks; idx_p1 > 0; --idx_p1) {
                        blocks_[idx_p1 - 1 + block_shift] = blocks_[idx_p1 - 1];
                    }
                }

                for (std::size_t idx = 0; idx < block_shift; ++idx) {
                    blocks_[idx] = 0;
                }

                return *this;
            }
            constexpr uint_var operator<<(std::size_t k) const& {
                auto r = *this;
                r <<= k;
                return r;
            }
            constexpr uint_var operator<<(std::size_t k) && {
                auto r = std::move(*this);
                r <<= k;
                return r;
            }

            constexpr uint_var& operator>>=(std::size_t k) noexcept {
                if (k == 0) {
                    return *this;
                }

                auto const prev_leading_one_pos = bit_width(*this);
                if (k >= prev_leading_one_pos) {
                    blocks_.clear();
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
                        blocks_[idx] = ((blocks_[idx + block_shift] >> bit_shift) |
                                        (blocks_[idx + block_shift + 1]
                                         << (number_of_bits_in_block - bit_shift))) &
                                       wuint::uint64_mask;
                    }
                    // The last block.
                    blocks_[number_of_blocks() - 1 - block_shift] =
                        (blocks_[number_of_blocks() - 1] >> bit_shift);
                }
                else {
                    for (std::size_t idx = 0; idx < number_of_blocks() - block_shift; ++idx) {
                        blocks_[idx] = blocks_[idx + block_shift];
                    }
                }

                blocks_.resize(new_number_of_blocks);
                return *this;
            }
            constexpr uint_var operator>>(std::size_t k) const& {
                auto r = *this;
                r >>= k;
                return r;
            }
            constexpr uint_var operator>>(std::size_t k) && {
                auto r = std::move(*this);
                r >>= k;
                return r;
            }

            // Find the largest power of 2 dividing *this, divide *this by that power, and return
            // the exponent.
            constexpr std::size_t factor_out_power_of_2() noexcept {
                util::constexpr_assert(!is_zero());

                std::size_t trailing_zero_blocks = 0;
                for (; trailing_zero_blocks < blocks_.size(); ++trailing_zero_blocks) {
                    if (blocks_[trailing_zero_blocks] != 0) {
                        break;
                    }
                }

                // Shift all blocks.
                for (std::size_t idx = trailing_zero_blocks; idx < blocks_.size(); ++idx) {
                    blocks_[idx - trailing_zero_blocks] = blocks_[idx];
                }
                blocks_.resize(blocks_.size() - trailing_zero_blocks);

                auto const trailing_zero_bits = std::size_t(std::countr_zero(blocks_.front()));
                if (trailing_zero_bits != 0) {
                    for (std::size_t idx = 0; idx < blocks_.size() - 1; ++idx) {
                        blocks_[idx] =
                            ((blocks_[idx] >> trailing_zero_bits) |
                             (blocks_[idx + 1] << (number_of_bits_in_block - trailing_zero_bits))) &
                            wuint::uint64_mask;
                    }
                    blocks_.back() >>= trailing_zero_bits;

                    if (blocks_.back() == 0) {
                        blocks_.pop_back();
                    }
                }

                return trailing_zero_blocks * number_of_bits_in_block + trailing_zero_bits;
            }

            // Return a list consisting of blocks of decimal digits of the stored number.
            // The list starts from the most significant digits, and each block after the first
            // block is consisting of 19 digits.
            constexpr std::vector<block_type> to_decimal() const {
                std::vector<std::uint64_t> digits;
                auto n_copy = *this;
                uint_var divisor{largest_pow10_in_block};
                while (!n_copy.is_zero()) {
                    auto q = n_copy.long_division(divisor);
                    if (n_copy.is_zero()) {
                        digits.push_back(0);
                    }
                    else {
                        digits.push_back(n_copy[0]);
                    }
                    n_copy = std::move(q);
                }

                // Reverse the list to place the most significant block at the first position.
                for (std::size_t idx = 0; idx < digits.size() / 2; ++idx) {
                    auto const temp = digits[idx];
                    digits[idx] = digits[digits.size() - idx - 1];
                    digits[digits.size() - idx - 1] = temp;
                }

                return digits;
            }

            // Converts a list of blocks of decimal digits into uint_var.
            // The list starts from the most significant digits, and each block is consisting of 19
            // digits. It may be better to use std::ranges concepts, but pulling the header <ranges>
            // would be too much of hassle.
            template <class Range>
            static constexpr uint_var from_decimal(Range const& digits) {
                uint_var result;
                uint_var multiplier = 1;
                for (std::size_t idx_p1 = digits.size(); idx_p1 > 0; --idx_p1) {
                    util::constexpr_assert(digits[idx_p1 - 1] < largest_pow10_in_block);
                    result += digits[idx_p1 - 1] * multiplier;
                    multiplier *= largest_pow10_in_block;
                }

                return result;
            }
        };

        constexpr uint_view to_view(uint_var const& x) noexcept { return uint_view(x); }

        constexpr uint_var operator+(uint_view x, uint_view y) {
            auto r = uint_var(x);
            r += y;
            return r;
        }
        constexpr uint_var operator+(uint_var&& x, uint_view y) { return std::move(x += y); }
        constexpr uint_var operator+(uint_view x, uint_var&& y) { return std::move(y += x); }
        constexpr uint_var operator+(uint_var&& x, uint_var&& y) { return std::move(x += y); }
        constexpr uint_var operator+(uint_view x, convertible_to_block_type auto y) {
            auto r = uint_var(x);
            r += y;
            return r;
        }
        constexpr uint_var operator+(uint_var&& x, convertible_to_block_type auto y) {
            return std::move(x += y);
        }
        constexpr uint_var operator+(uint_view x, uint_const_t<>) { return uint_var(x); }
        constexpr uint_var operator+(uint_var&& x, uint_const_t<>) { return std::move(x); }
        constexpr uint_var operator+(convertible_to_block_type auto x, uint_view y) {
            auto r = uint_var(y);
            r += x;
            return r;
        }
        constexpr uint_var operator+(convertible_to_block_type auto x, uint_var&& y) {
            return std::move(y += x);
        }
        constexpr uint_var operator+(uint_const_t<>, uint_view y) { return uint_var(y); }
        constexpr uint_var operator+(uint_const_t<>, uint_var&& y) { return std::move(y); }

        constexpr uint_var operator-(uint_view x, uint_view y) {
            auto r = uint_var(x);
            r -= y;
            return r;
        }
        constexpr uint_var operator-(uint_var&& x, uint_view y) { return std::move(x -= y); }
        constexpr uint_var operator-(uint_view x, convertible_to_block_type auto y) {
            auto r = uint_var(x);
            r -= y;
            return r;
        }
        constexpr uint_var operator-(uint_var&& x, convertible_to_block_type auto y) {
            return std::move(x -= y);
        }
        constexpr uint_var operator-(uint_view x, uint_const_t<>) { return uint_var(x); }
        constexpr uint_var operator-(uint_var&& x, uint_const_t<>) { return std::move(x); }
        constexpr uint_var operator-(convertible_to_block_type auto x, uint_view y) {
            auto r = uint_var(x);
            r -= y;
            return r;
        }

        constexpr uint_var operator*(uint_view x, uint_view y) {
            uint_var r;
            if (x.is_zero() || y.is_zero()) {
                return r;
            }
            r.blocks_.resize(x.blocks().size() + y.blocks().size(), 0);

            detail::multiply_impl(x.blocks(), y.blocks(), r.blocks_);
            if (r.blocks_.back() == 0) {
                r.blocks_.pop_back();
            }

            return r;
        }
        constexpr uint_var operator*(uint_view x, convertible_to_block_type auto y) {
            auto r = uint_var(x);
            r *= y;
            return r;
        }
        constexpr uint_var operator*(uint_var&& x, convertible_to_block_type auto y) {
            return std::move(x *= y);
        }
        constexpr uint_var operator*(uint_view x, uint_const_t<1u>) { return uint_var(x); }
        constexpr uint_var operator*(uint_view, uint_const_t<>) { return {}; }
        constexpr uint_var operator*(uint_var&& x, uint_const_t<1u>) { return std::move(x); }
        constexpr uint_var operator*(convertible_to_block_type auto x, uint_view y) {
            auto r = uint_var(y);
            r *= x;
            return r;
        }
        constexpr uint_var operator*(convertible_to_block_type auto x, uint_var&& y) {
            return std::move(y *= x);
        }
        constexpr uint_var operator*(uint_const_t<1u>, uint_view y) { return uint_var(y); }
        constexpr uint_var operator*(uint_const_t<1u>, uint_var&& y) { return std::move(y); }
        constexpr uint_var operator*(uint_const_t<>, uint_view) { return {}; }

        constexpr uint_var& operator*=(uint_var& x, uint_view y) {
            auto r = x * y;
            x = std::move(r);
            return x;
        }

        constexpr uint_var operator/(uint_view x, uint_view y) {
            auto r = uint_var(x);
            return r.long_division(y);
        }
        constexpr uint_var operator/(uint_var&& x, uint_view y) { return x.long_division(y); }
        constexpr uint_var operator/(uint_view x, convertible_to_block_type auto y) {
            auto r = uint_var(x);
            return r.long_division(uint_view::make_view_from_single_block(y));
        }
        constexpr uint_var operator/(uint_var&& x, convertible_to_block_type auto y) {
            return x.long_division(uint_view::make_view_from_single_block(y));
        }
        constexpr uint_var operator/(uint_view x, uint_const_t<1u>) { return uint_var(x); }
        constexpr uint_var operator/(uint_var&& x, uint_const_t<1u>) { return std::move(x); }
        constexpr uint_var operator/(uint_view, uint_const_t<>) = delete;

        constexpr uint_var& operator/=(uint_var& x, uint_view y) {
            auto r = x / y;
            x = std::move(r);
            return x;
        }
        constexpr uint_var& operator/=(uint_var& x, convertible_to_block_type auto y) {
            auto r = x / y;
            x = std::move(r);
            return x;
        }
        constexpr uint_var& operator/=(uint_var& x, uint_const_t<1u>) { return x; }
        constexpr uint_var& operator/=(uint_var& x, uint_const_t<>) = delete;

        constexpr uint_var operator%(uint_view x, uint_view y) {
            auto r = uint_var(x);
            r.long_division(y);
            return r;
        }
        constexpr uint_var operator%(uint_var&& x, uint_view y) {
            x.long_division(y);
            return std::move(x);
        }
        constexpr uint_var operator%(uint_view x, convertible_to_block_type auto y) {
            auto r = uint_var(x);
            r.long_division(uint_view::make_view_from_single_block(y));
            return r;
        }
        constexpr uint_var operator%(uint_var&& x, convertible_to_block_type auto y) {
            x.long_division(uint_view::make_view_from_single_block(y));
            return std::move(x);
        }
        constexpr uint_var operator%(uint_view, uint_const_t<1u>) { return {}; }
        constexpr uint_var operator%(uint_view, uint_const_t<>) = delete;

        using uint_var_div_t = util::div_t<uint_var>;
        constexpr uint_var_div_t div(uint_view x, uint_view y) {
            uint_var_div_t ret;
            ret.rem = uint_var(x);
            ret.quot = ret.rem.long_division(y);
            return ret;
        }
        constexpr uint_var_div_t div(uint_view x, uint_const_t<1u>) { return {uint_var(x), {}}; }
        constexpr uint_var_div_t div(uint_view x, uint_const_t<>) = delete;
        constexpr uint_var_div_t div(uint_var&& x, uint_view y) {
            uint_var_div_t ret;
            ret.quot = x.long_division(y);
            ret.rem = std::move(x);
            return ret;
        }
        constexpr uint_var_div_t div(uint_var&& x, uint_const_t<1u>) { return {std::move(x), {}}; }

        constexpr uint_var div_floor(uint_view x, uint_view y) { return x / y; }
        constexpr uint_var div_floor(uint_view x, uint_const_t<1u>) { return uint_var(x); }
        constexpr uint_var div_floor(uint_view x, uint_const_t<>) = delete;
        constexpr uint_var div_floor(uint_var&& x, uint_view y) { return x.long_division(y); }
        constexpr uint_var div_floor(uint_var&& x, uint_const_t<1u>) { return std::move(x); }
        constexpr uint_var div_floor(uint_var&& x, uint_const_t<>) = delete;

        constexpr uint_var div_ceil(uint_view x, uint_view y) {
            auto dividend = uint_var(x);
            auto quotient = dividend.long_division(y);
            if (!dividend.is_zero()) {
                ++quotient;
            }
            return quotient;
        }
        constexpr uint_var div_ceil(uint_view x, uint_const_t<1u>) { return uint_var(x); }
        constexpr uint_var div_ceil(uint_view x, uint_const_t<>) = delete;
        constexpr uint_var div_ceil(uint_var&& x, uint_view y) {
            auto quotient = x.long_division(y);
            if (!x.is_zero()) {
                ++quotient;
            }
            return quotient;
        }
        constexpr uint_var div_ceil(uint_var&& x, uint_const_t<1u>) { return std::move(x); }

        // Computes max(floor(log2(x / y)), 0).
        // Precondition: x, y are not zero.
        constexpr std::size_t trunc_floor_log2_div(uint_view x, uint_view y) noexcept {
            util::constexpr_assert<util::error_msgs::divide_by_zero>(!x.is_zero() && !y.is_zero());

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
                    auto const x_block = x[y.number_of_blocks() + block_shift];
                    auto const y_block =
                        (y[y.number_of_blocks() - 1] >> (number_of_bits_in_block - bit_shift));

                    if (x_block < y_block) {
                        return total_shift - 1;
                    }
                    else if (x_block > y_block) {
                        return total_shift;
                    }
                }
                // Middle blocks of y.
                for (std::size_t idx = y.number_of_blocks() - 1; idx > 0; --idx) {
                    auto const x_block = x[idx + block_shift];
                    auto const y_block = ((y[idx] << bit_shift) |
                                          (y[idx - 1] >> (number_of_bits_in_block - bit_shift))) &
                                         wuint::uint64_mask;
                    if (x_block < y_block) {
                        return total_shift - 1;
                    }
                    else if (x_block < y_block) {
                        return total_shift;
                    }
                }
                // The last block of y.
                {
                    auto const x_block = x[block_shift];
                    auto const y_block = ((y[0] << bit_shift) & wuint::uint64_mask);

                    if (x_block < y_block) {
                        return total_shift - 1;
                    }
                }
            }
            else {
                for (std::size_t idx_p1 = y.number_of_blocks(); idx_p1 > 0; --idx_p1) {
                    auto const x_block = x[idx_p1 - 1 + block_shift];
                    auto const y_block = y[idx_p1 - 1];
                    if (x_block < y_block) {
                        return total_shift - 1;
                    }
                    else if (x_block > y_block) {
                        return total_shift;
                    }
                }
            }
            return total_shift;
        }

        // Computes max(ceil(log2(x / y)), 0).
        // Precondition: x, y are not zero.
        constexpr std::size_t trunc_ceil_log2_div(uint_view x, uint_view y) noexcept {
            util::constexpr_assert<util::error_msgs::divide_by_zero>(!x.is_zero() && !y.is_zero());

            auto const x_leading_one_pos = bit_width(x);
            auto const y_leading_one_pos = bit_width(y);

            if (y_leading_one_pos > x_leading_one_pos) {
                return 0;
            }

            auto const total_shift = x_leading_one_pos - y_leading_one_pos;
            auto const block_shift = total_shift / number_of_bits_in_block;
            auto const bit_shift = total_shift % number_of_bits_in_block;

            if (bit_shift != 0) {
                // The first block of y.
                if (y.number_of_blocks() + block_shift < x.number_of_blocks()) {
                    auto const x_block = x[y.number_of_blocks() + block_shift];
                    auto const y_block =
                        (y[y.number_of_blocks() - 1] >> (number_of_bits_in_block - bit_shift));

                    if (x_block > y_block) {
                        return total_shift + 1;
                    }
                    else if (x_block < y_block) {
                        return total_shift;
                    }
                }
                // Middle blocks of y.
                for (std::size_t idx = y.number_of_blocks() - 1; idx > 0; --idx) {
                    auto const x_block = x[idx + block_shift];
                    auto const y_block = ((y[idx] << bit_shift) |
                                          (y[idx - 1] >> (number_of_bits_in_block - bit_shift))) &
                                         wuint::uint64_mask;
                    if (x_block > y_block) {
                        return total_shift + 1;
                    }
                    else if (x_block < y_block) {
                        return total_shift;
                    }
                }
                // The last block of y.
                {
                    auto const x_block = x[block_shift];
                    auto const y_block = ((y[0] << bit_shift) & wuint::uint64_mask);

                    if (x_block > y_block) {
                        return total_shift + 1;
                    }
                }
            }
            else {
                for (std::size_t idx_p1 = y.number_of_blocks(); idx_p1 > 0; --idx_p1) {
                    auto const x_block = x[idx_p1 - 1 + block_shift];
                    auto const y_block = y[idx_p1 - 1];
                    if (x_block > y_block) {
                        return total_shift + 1;
                    }
                    else if (x_block < y_block) {
                        return total_shift;
                    }
                }
            }
            return total_shift;
        }

        // Find the largest power of 2 dividing n, divide n by that power, and return
        // the exponent.
        constexpr std::size_t factor_out_power_of_2(uint_var& n) noexcept {
            return n.factor_out_power_of_2();
        }

        // Operations on uint_const.
        namespace detail {
            // Used for converting uint_var into uint_const.
            template <std::size_t N>
            struct block_holder_size_pair {
                static_block_holder<N> arr{};
                std::size_t size = 0;

                constexpr block_holder_size_pair(uint_view n) noexcept {
                    for (std::size_t idx = 0; idx < n.number_of_blocks(); ++idx) {
                        arr[idx] = n[idx];
                    }
                    size = n.number_of_blocks();
                }
            };

            template <auto pair>
            using block_holder_size_pair_to_uint_const =
                uint_const_impl<slice<pair.size>(pair.arr)>;

            template <std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator+(uint_const_impl<x> xx, uint_const_impl<y> yy) noexcept {
                if constexpr (M == 0) {
                    return yy;
                }
                else {
                    constexpr auto result = block_holder_size_pair<util::max(M, N) + 1>{
                        to_view(decltype(xx){}) + decltype(yy){}};
                    return block_holder_size_pair_to_uint_const<result>{};
                }
            }

            template <std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
                requires(M >= N)
            constexpr auto operator-(uint_const_impl<x> xx, uint_const_impl<y> yy) noexcept {
                if constexpr (M == 0) {
                    return uint_const_t<>{};
                }
                else {
                    constexpr auto result = block_holder_size_pair<util::max(M, N)>{
                        to_view(decltype(xx){}) - decltype(yy){}};
                    return block_holder_size_pair_to_uint_const<result>{};
                }
            }

            template <std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator*(uint_const_impl<x> xx, uint_const_impl<y> yy) noexcept {
                if constexpr (M == 0) {
                    return uint_const_t<>{};
                }
                else {
                    constexpr auto result =
                        block_holder_size_pair<M + N>{to_view(decltype(xx){}) * decltype(yy){}};
                    return block_holder_size_pair_to_uint_const<result>{};
                }
            }

            template <class Quotient, class Remainder>
            struct uint_const_div_t;

            template <std::size_t M, std::size_t N, static_block_holder<M> quotient,
                      static_block_holder<N> remainder>
            struct uint_const_div_t<uint_const_impl<quotient>, uint_const_impl<remainder>> {
                uint_const_impl<quotient> quot;
                uint_const_impl<remainder> rem;
            };

            template <std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
                requires(N != 0)
            constexpr auto div(uint_const_impl<x> xx, uint_const_impl<y> yy) noexcept {
                if constexpr (M == 0) {
                    return uint_const_div_t<uint_const_t<>, uint_const_t<>>{};
                }
                else {
                    constexpr auto result = [] {
                        auto result_var = div(to_view(decltype(xx){}), to_view(decltype(yy){}));
                        auto quot =
                            block_holder_size_pair<M - util::min(M, N) + 1>(result_var.quot);
                        auto rem = block_holder_size_pair<util::min(M, N)>(result_var.rem);

                        struct intermediate_result {
                            decltype(quot) quot_;
                            decltype(rem) rem_;
                        };
                        return intermediate_result{quot, rem};
                    }();
                    return uint_const_div_t<block_holder_size_pair_to_uint_const<result.quot_>,
                                            block_holder_size_pair_to_uint_const<result.rem_>>{};
                }
            }

            template <std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
                requires(N != 0)
            constexpr auto operator/(uint_const_impl<x> xx, uint_const_impl<y> yy) noexcept {
                return div(xx, yy).quot;
            }

            template <std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
                requires(N != 0)
            constexpr auto operator%(uint_const_impl<x> xx, uint_const_impl<y> yy) noexcept {
                return div(xx, yy).rem;
            }

            template <std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator<<(uint_const_impl<x> xx, uint_const_impl<y>) noexcept {
                static_assert(N <= 1, "jkj::idiv: shift amount too large");
                if constexpr (N == 0) {
                    return xx;
                }
                else {
                    constexpr auto result =
                        block_holder_size_pair<M + (y[0] + number_of_bits_in_block - 1) /
                                                       number_of_bits_in_block>{
                            to_view(decltype(xx){}) << y[0]};
                    return block_holder_size_pair_to_uint_const<result>{};
                }
            }

            template <std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator>>(uint_const_impl<x> xx, uint_const_impl<y>) noexcept {
                static_assert(N <= 1, "jkj::idiv: shift amount too large");
                if constexpr (N == 0) {
                    return xx;
                }
                else {
                    constexpr auto result =
                        block_holder_size_pair<M>{to_view(decltype(xx){}) >> y[0]};
                    return block_holder_size_pair_to_uint_const<result>{};
                }
            }

            // Return a list consisting of blocks of decimal digits of the stored number.
            // The list starts from the most significant digits, and each block after the first
            // block is consisting of 19 digits.
            template <std::size_t M, static_block_holder<M> x>
            constexpr auto to_decimal(uint_const_impl<x>) noexcept {
                // Compute an upper bound on the number of 19-digits blocks that will be
                // generated. An obvious upper bound is floor(64M/(19 log2(10))) + 1, and
                // since 72/71 > 64/(19 log2(10)) holds, the following is an upper bound.
                constexpr std::size_t number_of_19_digits_blocks = ((72 * M) / 71) + 1;
                constexpr auto result = block_holder_size_pair<number_of_19_digits_blocks>(
                    uint_var{uint_const_impl<x>{}}.to_decimal());

                return slice<result.size>(result.arr);
            }

            // Converts a list of blocks of decimal digits into uint_var.
            // The list starts from the most significant digits, and each block is consisting of 19
            // digits.
            template <auto arr>
            constexpr auto from_decimal() noexcept {
                // Compute an upper bound on the number of blocks that will be
                // generated. An obvious upper bound is floor(19N * log2(10) / 64)) + 1, and
                // since 143/145 > 19 log2(10)/64 holds, the following is an upper bound.
                constexpr std::size_t number_of_blocks = ((143 * arr.size()) / 145) + 1;
                constexpr auto result =
                    block_holder_size_pair<number_of_blocks>(uint_var::from_decimal(arr));

                return block_holder_size_pair_to_uint_const<result>();
            }
        }

        using detail::to_decimal;
        using detail::from_decimal;

        template <block_type... blocks>
        inline constexpr uint_const_t<blocks...> uint_const_v{};

        // Convert a sequence of decimal digit blocks, each consisting of 19-digits, into a
        // uint_const. E.g., decimal_uint_const_t<1234, 1234> represents
        // 1234'000'0000'0000'0000'1234. Any leading zeros are ignored.
        template <block_type... decimal_digit_blocks>
        using decimal_uint_const_t =
            decltype(detail::from_decimal<static_block_holder<sizeof...(decimal_digit_blocks)>{
                         decimal_digit_blocks...}>());

        template <block_type... decimal_digit_blocks>
        inline constexpr decimal_uint_const_t<decimal_digit_blocks...> decimal_uint_const_v{};


        // Signed version of uint_var.
        // Internally uses sign-magnitude representation, but does not allow negative zero.
        // The entire class is mostly strongly exception-safe:
        // - every const member function is noexcept and does not actually change anything,
        // - in any non-const member function, the only source of exceptions is memory allocation
        // faillure, and if it throws, the previous state of the class is retained.
        // - Some binary operations involving rvalue-reference parameters can modify some of those
        // parameters, with or without the presene of exceptions. But they always leave them in a
        // well-defined state.
        // - Any operations should work just fine even if operands alias to each other.
        //
        // Note that, however, violation of preconditions may terminate (through assert) or lead to
        // undefined behavior (if assert does not terminate).
        class int_var {
            // Absolute value.
            uint_var abs_;
            // Sign.
            sign_t sign_ = sign_t::positive;

        public:
            friend class uint_var;

            int_var() = default;
            explicit constexpr int_var(sign_t sign, uint_view abs) : abs_(abs), sign_(sign) {
                if (abs_.is_zero()) {
                    sign_ = sign_t::positive;
                }
            }
            explicit constexpr int_var(sign_t sign, uint_var&& abs)
                : abs_(std::move(abs)), sign_(sign) {
                if (abs_.is_zero()) {
                    sign_ = sign_t::positive;
                }
            }

            explicit constexpr int_var(uint_view abs) : int_var(sign_t::positive, abs) {}
            explicit constexpr int_var(uint_var&& abs)
                : int_var(sign_t::positive, std::move(abs)) {}
            constexpr int_var(int_view n) : int_var(n.sign(), n.abs()) {}

            constexpr int_var(signed_block_type n) {
                // Signed-to-unsigned conversion is well-defined.
                if (n > 0) {
                    auto n_abs = static_cast<block_type>(n);
                    util::constexpr_assert<util::error_msgs::overflow>(n_abs <= wuint::uint64_mask);
                    abs_.blocks_.push_back(n_abs);
                }
                else if (n < 0) {
                    auto n_abs = block_type(0) - static_cast<block_type>(n);
                    util::constexpr_assert<util::error_msgs::overflow>(n_abs <= wuint::uint64_mask);
                    sign_ = sign_t::negative;
                    abs_.blocks_.push_back(n_abs);
                }
            }

            constexpr int_var(uint_const_t<>) {}
            constexpr int_var(int_const_t<>) {}

            // Obtain view.
            constexpr operator int_view() const noexcept { return int_view{sign_, abs_.blocks()}; }

            constexpr util::span<block_type const> blocks() const noexcept { return abs_.blocks(); }
            constexpr std::size_t number_of_blocks() const noexcept {
                return abs_.number_of_blocks();
            }
            constexpr block_type operator[](std::size_t idx) const { return abs_[idx]; }

            constexpr bool is_zero() const noexcept { return abs_.is_zero(); }
            constexpr bool is_even() const noexcept { return abs_.is_even(); }

            constexpr sign_t sign() const noexcept { return sign_; }
            constexpr bool is_strictly_positive() const noexcept {
                return !is_zero() && sign_ == sign_t::positive;
            }
            constexpr bool is_strictly_negative() const noexcept {
                return sign_ == sign_t::negative;
            }
            constexpr bool is_nonnegative() const noexcept { return sign_ == sign_t::positive; }
            constexpr bool is_nonpositive() const noexcept {
                return is_zero() || sign_ == sign_t::negative;
            }

            // Absolute value.
            constexpr uint_var const& abs() const& noexcept { return abs_; }
            constexpr uint_var&& abs() && noexcept { return std::move(abs_); }

            constexpr int_var& invert_sign() noexcept {
                if (!is_zero()) {
                    sign_ = ::jkj::bigint::invert_sign(sign_);
                }
                return *this;
            }

            constexpr int_var operator-() const& {
                return int_var(::jkj::bigint::invert_sign(sign_), abs_);
            }

            constexpr int_var operator-() && {
                return int_var(::jkj::bigint::invert_sign(sign_), std::move(abs_));
            }

            constexpr int_var& operator+=(uint_view n) {
                if (is_nonnegative()) {
                    abs_ += n;
                }
                else {
                    if (abs_ > n) {
                        abs_ -= n;
                    }
                    else {
                        abs_ = n - abs_;
                        sign_ = sign_t::positive;
                    }
                }
                return *this;
            }
            constexpr int_var& operator+=(int_view n) {
                if (sign_ == n.sign()) {
                    abs_ += n.abs();
                }
                else {
                    auto comp = abs_ <=> n.abs();

                    if (comp < 0) {
                        abs_ = n.abs() - abs_;
                        sign_ = ::jkj::bigint::invert_sign(sign_);
                    }
                    else if (comp > 0) {
                        abs_ -= n.abs();
                    }
                    else {
                        abs_.blocks_.clear();
                        sign_ = sign_t::positive;
                    }
                }
                return *this;
            }
            constexpr int_var& operator+=(convertible_to_block_type auto n) {
                if (is_nonnegative()) {
                    abs_ += n;
                }
                else {
                    if (abs_ > n) {
                        abs_ -= n;
                    }
                    else {
                        // abs_ should be consisting of exactly one block.
                        abs_.blocks_[0] = n - abs_.blocks_[0];
                        if (abs_.blocks_[0] == 0) {
                            abs_.blocks_.clear();
                        }
                        sign_ = sign_t::positive;
                    }
                }
                return *this;
            }
            constexpr int_var& operator+=(convertible_to_signed_block_type auto n) {
                // Signed-to-unsigned conversion is well-defined.
                auto const n_abs = n >= 0 ? static_cast<block_type>(n)
                                          : block_type(0) - static_cast<block_type>(n);

                if (is_nonnegative() == (n >= 0)) {
                    abs_ += n_abs;
                }
                else {
                    if (abs_ >= n_abs) {
                        abs_ -= n_abs;
                        if (abs_.is_zero()) {
                            sign_ = sign_t::positive;
                        }
                    }
                    else {
                        // abs_ should be either zero or consisting of a single block.
                        if (is_zero()) {
                            abs_.blocks_.resize(1, n_abs);
                        }
                        else {
                            abs_.blocks_[0] = n_abs - abs_.blocks_[0];
                        }
                        sign_ = ::jkj::bigint::invert_sign(sign_);
                    }
                }
                return *this;
            }
            constexpr int_var& operator+=(uint_const_t<>) { return *this; }
            constexpr int_var& operator+=(int_const_t<>) { return *this; }

            constexpr int_var& operator++() {
                if (is_nonnegative()) {
                    ++abs_;
                }
                else {
                    --abs_;
                    if (abs_.is_zero()) {
                        sign_ = sign_t::positive;
                    }
                }
                return *this;
            }
            constexpr int_var operator++(int) & {
                auto temp = *this;
                ++*this;
                return temp;
            }
            constexpr int_var operator++(int) && {
                auto temp = std::move(*this);
                ++*this;
                return temp;
            }

            constexpr int_var& operator-=(uint_view n) {
                if (is_strictly_negative()) {
                    abs_ += n;
                }
                else {
                    if (abs_ >= n) {
                        abs_ -= n;
                        if (abs_.is_zero()) {
                            sign_ = sign_t::positive;
                        }
                    }
                    else {
                        abs_ = n - abs_;
                        sign_ = sign_t::negative;
                    }
                }
                return *this;
            }
            constexpr int_var& operator-=(int_view n) {
                if (n.is_zero()) {
                    return *this;
                }
                return *this += int_view{::jkj::bigint::invert_sign(n.sign()), n.blocks_};
            }
            constexpr int_var& operator-=(convertible_to_block_type auto n) {
                if (is_strictly_negative()) {
                    abs_ += n;
                }
                else {
                    if (abs_ >= n) {
                        abs_ -= n;
                        if (abs_.is_zero()) {
                            sign_ = sign_t::positive;
                        }
                    }
                    else {
                        if (is_zero()) {
                            abs_.blocks_.resize(1, n);
                        }
                        else {
                            abs_.blocks_[0] = n - abs_.blocks_[0];
                        }
                        sign_ = sign_t::negative;
                    }
                }
                return *this;
            }
            constexpr int_var& operator-=(convertible_to_signed_block_type auto n) {
                // Signed-to-unsigned conversion is well-defined.
                auto const n_abs = n >= 0 ? static_cast<block_type>(n)
                                          : block_type(0) - static_cast<block_type>(n);

                if (is_strictly_negative() == (n >= 0)) {
                    abs_ += n_abs;
                }
                else {
                    if (abs_ >= n_abs) {
                        abs_ -= n_abs;
                        if (abs_.is_zero()) {
                            sign_ = sign_t::positive;
                        }
                    }
                    else {
                        // abs_ should be either zero or consisting of a single block.
                        if (is_zero()) {
                            abs_.blocks_.resize(1, n_abs);
                        }
                        else {
                            abs_.blocks_[0] = n_abs - abs_.blocks_[0];
                        }
                        sign_ = ::jkj::bigint::invert_sign(sign_);
                    }
                }
                return *this;
            }
            constexpr int_var& operator-=(uint_const_t<>) { return *this; }
            constexpr int_var& operator-=(int_const_t<>) { return *this; }

            constexpr int_var& operator--() {
                if (is_strictly_negative()) {
                    ++abs_;
                }
                else if (is_zero()) {
                    abs_.blocks_.resize(1, 1);
                    sign_ = sign_t::negative;
                }
                else {
                    --abs_;
                    if (abs_.is_zero()) {
                        sign_ = sign_t::positive;
                    }
                }
                return *this;
            }
            constexpr int_var operator--(int) & {
                auto temp = *this;
                --*this;
                return temp;
            }
            constexpr int_var operator--(int) && {
                auto temp = std::move(*this);
                --*this;
                return temp;
            }

            constexpr int_var& operator*=(convertible_to_block_type auto n) {
                abs_ *= n;
                if (abs_.is_zero()) {
                    sign_ = sign_t::positive;
                }
                return *this;
            }
            constexpr int_var& operator*=(convertible_to_signed_block_type auto n) {
                // Signed-to-unsigned conversion is well-defined.
                if (n >= 0) {
                    abs_ *= static_cast<block_type>(n);
                }
                else {
                    abs_ *= (block_type(0) - static_cast<block_type>(n));
                    sign_ = ::jkj::bigint::invert_sign(sign_);
                }

                if (abs_.is_zero()) {
                    sign_ = sign_t::positive;
                }
                return *this;
            }
            constexpr int_var& operator*=(uint_const_t<1u>) { return *this; }
            constexpr int_var& operator*=(int_const_t<sign_t::positive, 1u>) { return *this; }
            constexpr int_var& operator*=(int_const_t<sign_t::negative, 1u>) {
                invert_sign();
                return *this;
            }
            constexpr int_var& operator*=(uint_const_t<>) {
                abs_.blocks_.clear();
                sign_ = sign_t::positive;
                return *this;
            }
            constexpr int_var& operator*=(int_const_t<>) { return (*this) *= uint_const_t<>{}; }

            // Perform long division.
            // *this becomes the remainder, returns the quotient.
            // The divisor is assumed to be positive.
            // Precondition: n != 0
            constexpr int_var long_division(uint_view n) {
                util::constexpr_assert<util::error_msgs::divide_by_zero>(!n.is_zero());

                // Do the usual unsigned division if the dividend is nonnegative.
                if (is_nonnegative()) {
                    auto quotient = abs_.long_division(n);
                    sign_ = sign_t::positive;
                    return int_var{sign_t::positive, std::move(quotient)};
                }
                // Otherwise, use the identities
                // floor(-p/q) = -ceil(p/q) = -floor((p + q - 1)/q) and
                // -p - floor(-p/q)q = -p + floor((p + q - 1)/q)q
                // = (q - 1) - ((p + q - 1) - floor((p + q - 1)/q)q).
                else {
                    auto q_m1 = uint_var(n);
                    --q_m1;

                    auto dividend = abs_;
                    dividend += q_m1;

                    auto quotient = dividend.long_division(n);
                    abs_ = std::move(q_m1);
                    abs_ -= dividend;
                    sign_ = sign_t::positive;

                    return int_var{sign_t::negative, std::move(quotient)};
                }
            }
            constexpr int_var long_division(uint_const_t<1u>) {
                auto quotient = std::move(*this);
                util::constexpr_assert(abs_.blocks_.empty());
                sign_ = sign_t::positive;
                return quotient;
            }
            constexpr int_var long_division(uint_const_t<>) = delete;

            constexpr int_var& operator<<=(std::size_t k) {
                abs_ <<= k;
                return *this;
            }

            constexpr int_var operator<<(std::size_t k) const& {
                auto r = *this;
                r <<= k;
                return r;
            }
            constexpr int_var operator<<(std::size_t k) && {
                auto r = std::move(*this);
                r <<= k;
                return r;
            }

            // Performs arithemtic shift: x |-> floor(x/2^k).
            constexpr int_var& operator>>=(std::size_t k) {
                if (is_nonnegative()) {
                    abs_ >>= k;
                    return *this;
                }
                else {
                    // floor(-x/2^k) = -ceil(x/2^k), so check if the lower k-bits are all zero.
                    if (is_zero()) {
                        return *this;
                    }

                    bool divisible = [&] {
                        std::size_t bits_to_check = k;
                        for (std::size_t idx = 0; idx < number_of_blocks(); ++idx) {
                            if (bits_to_check >= number_of_bits_in_block) {
                                if (blocks()[idx] != 0) {
                                    return false;
                                }
                                bits_to_check -= number_of_bits_in_block;
                            }
                            else {
                                if (static_cast<std::size_t>(std::countr_zero(blocks()[idx])) <
                                    bits_to_check) {
                                    return false;
                                }
                                return true;
                            }
                        }
                        return bits_to_check == 0;
                    }();

                    abs_ >>= k;
                    if (!divisible) {
                        ++abs_;
                    }
                    return *this;
                }
            }

            constexpr int_var operator>>(std::size_t k) const& {
                auto r = *this;
                r >>= k;
                return r;
            }
            constexpr int_var operator>>(std::size_t k) && {
                auto r = std::move(*this);
                r >>= k;
                return r;
            }

            // Find the largest power of 2 dividing *this, divide *this by that power, and return
            // the exponent.
            constexpr std::size_t factor_out_power_of_2() noexcept {
                return abs_.factor_out_power_of_2();
            }
        };

        constexpr int_view to_view(int_var const& x) noexcept { return int_view(x); }

        // Conversion between unsigned entities and signed entities.
        constexpr int_view to_signed(uint_view n) noexcept {
            return int_view{sign_t::positive, n.blocks()};
        };
        constexpr int_view to_negative(uint_view n) noexcept {
            return int_view{n.is_zero() ? sign_t::positive : sign_t::negative, n.blocks()};
        };
        constexpr uint_view abs(uint_view n) noexcept { return n; }
        constexpr uint_view abs(int_view n) noexcept { return n.abs(); }

        constexpr int_var to_signed(uint_var const& n) { return int_var(n); }
        constexpr int_var to_negative(uint_var const& n) { return int_var(sign_t::negative, n); }
        constexpr uint_var abs(uint_var const& n) { return n; }
        constexpr uint_var abs(int_var const& n) { return n.abs(); }

        constexpr int_var to_signed(uint_var&& n) noexcept { return int_var(std::move(n)); }
        constexpr int_var to_negative(uint_var&& n) noexcept {
            return int_var{sign_t::negative, std::move(n)};
        }
        constexpr uint_var abs(uint_var&& n) noexcept { return std::move(n); }
        constexpr uint_var abs(int_var&& n) noexcept { return std::move(n).abs(); }

        namespace detail {
            template <std::size_t N, static_block_holder<N> arr>
            constexpr int_const_impl<sign_t::positive, arr>
            to_signed(uint_const_impl<arr>) noexcept {
                return {};
            }
            template <std::size_t N, static_block_holder<N> arr>
            constexpr int_const_impl<(arr.size() == 0 ? sign_t::positive : sign_t::negative), arr>
            to_negative(uint_const_impl<arr>) noexcept {
                return {};
            }
            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            constexpr uint_const_impl<arr> abs(int_const_impl<sign, arr>) noexcept {
                return {};
            }
        }
        using detail::to_signed;
        using detail::to_negative;
        using detail::abs;

        // Some inspection functions.
        constexpr bool is_zero(uint_view n) noexcept { return n.is_zero(); }
        constexpr bool is_even(uint_view n) noexcept { return n.is_even(); }
        constexpr bool is_strictly_positive(uint_view n) noexcept { return !n.is_zero(); }
        constexpr bool is_strictly_negative(uint_view) noexcept { return false; }
        constexpr bool is_nonnegative(uint_view) noexcept { return true; }
        constexpr bool is_nonpositive(uint_view n) noexcept { return n.is_zero(); }
        constexpr sign_t sign(uint_view) noexcept { return sign_t::positive; }

        constexpr bool is_zero(int_view n) noexcept { return n.is_zero(); }
        constexpr bool is_even(int_view n) noexcept { return n.is_even(); }
        constexpr bool is_strictly_positive(int_view n) noexcept {
            return n.is_strictly_positive();
        }
        constexpr bool is_strictly_negative(int_view n) noexcept {
            return n.is_strictly_negative();
        }
        constexpr bool is_nonnegative(int_view n) noexcept { return n.is_nonnegative(); }
        constexpr bool is_nonpositive(int_view n) noexcept { return n.is_nonpositive(); }
        constexpr sign_t sign(int_view n) noexcept { return n.sign(); }

        namespace detail {
            template <std::size_t N, static_block_holder<N> arr>
            constexpr bool is_zero(uint_const_impl<arr> n) noexcept {
                return n.is_zero();
            }
            template <std::size_t N, static_block_holder<N> arr>
            constexpr bool is_even(uint_const_impl<arr> n) noexcept {
                return n.is_even();
            }
            template <std::size_t N, static_block_holder<N> arr>
            constexpr bool is_strictly_positive(uint_const_impl<arr> n) noexcept {
                return !n.is_zero();
            }
            template <std::size_t N, static_block_holder<N> arr>
            constexpr bool is_strictly_negative(uint_const_impl<arr>) noexcept {
                return false;
            }
            template <std::size_t N, static_block_holder<N> arr>
            constexpr bool is_nonnegative(uint_const_impl<arr>) noexcept {
                return true;
            }
            template <std::size_t N, static_block_holder<N> arr>
            constexpr bool is_nonpositive(uint_const_impl<arr> n) noexcept {
                return n.is_zero();
            }
            template <std::size_t N, static_block_holder<N> arr>
            constexpr sign_t sign(uint_const_impl<arr>) noexcept {
                return sign_t::positive;
            }

            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            constexpr bool is_zero(int_const_impl<sign, arr> n) noexcept {
                return n.is_zero();
            }
            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            constexpr bool is_even(int_const_impl<sign, arr> n) noexcept {
                return n.is_even();
            }
            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            constexpr bool is_strictly_positive(int_const_impl<sign, arr> n) noexcept {
                return n.is_strictly_positive();
            }
            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            constexpr bool is_strictly_negative(int_const_impl<sign, arr> n) noexcept {
                return n.is_strictly_negative();
            }
            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            constexpr bool is_nonnegative(int_const_impl<sign, arr> n) noexcept {
                return n.is_nonnegative();
            }
            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            constexpr bool is_nonpositive(int_const_impl<sign, arr> n) noexcept {
                return n.is_nonpositive();
            }
            template <sign_t sign_, std::size_t N, static_block_holder<N> arr>
            constexpr sign_t sign(int_const_impl<sign_, arr> n) noexcept {
                return n.sign();
            }
        }
        using detail::is_zero;
        using detail::is_even;
        using detail::is_strictly_positive;
        using detail::is_strictly_negative;
        using detail::is_nonnegative;
        using detail::is_nonpositive;
        using detail::sign;

        // Sign inversion functions.
        constexpr int_var invert_sign(int_view n) {
            return int_var(::jkj::bigint::invert_sign(n.sign()), uint_var(n.abs()));
        }
        constexpr int_var invert_sign(int_var&& n) {
            return int_var(::jkj::bigint::invert_sign(n.sign()), std::move(n).abs());
        }
        namespace detail {
            template <sign_t sign, std::size_t N, static_block_holder<N> arr>
            constexpr int_const_impl<::jkj::bigint::invert_sign(sign), arr>
            invert_sign(int_const_impl<sign, arr>) noexcept {
                return {};
            }
            template <sign_t sign, static_block_holder<0> arr>
            constexpr int_const_impl<sign_t::positive, arr>
            invert_sign(int_const_impl<sign, arr>) noexcept {
                return {};
            }
        }
        using detail::invert_sign;


        constexpr int_var operator+(uint_view x, int_view y) {
            auto r = int_var(y);
            r += x;
            return r;
        }
        constexpr int_var operator+(uint_view x, int_var&& y) { return std::move(y += x); }
        constexpr int_var operator+(uint_view x, convertible_to_signed_block_type auto y) {
            auto r = int_var(y);
            r += x;
            return r;
        }
        constexpr int_var operator+(uint_view x, int_const_t<>) { return int_var(uint_var(x)); }

        constexpr int_var operator+(uint_var&& x, int_view y) {
            auto r = to_signed(std::move(x));
            r += y;
            return r;
        }
        constexpr int_var operator+(uint_var&& x, int_var&& y) { return std::move(y += x); }
        constexpr int_var operator+(uint_var&& x, convertible_to_signed_block_type auto y) {
            auto r = to_signed(std::move(x));
            r += y;
            return r;
        }
        constexpr int_var operator+(uint_var&& x, int_const_t<>) { return to_signed(std::move(x)); }

        constexpr int_var operator+(int_view x, uint_view y) {
            auto r = int_var(x);
            r += y;
            return r;
        }
        constexpr int_var operator+(int_view x, uint_var&& y) {
            auto r = to_signed(std::move(y));
            r += x;
            return r;
        }
        constexpr int_var operator+(int_view x, int_view y) {
            auto r = int_var(x);
            r += y;
            return r;
        }
        constexpr int_var operator+(int_view x, int_var&& y) { return std::move(y += x); }
        constexpr int_var operator+(int_view x, convertible_to_block_type auto y) {
            auto r = int_var(x);
            r += y;
            return r;
        }
        constexpr int_var operator+(int_view x, convertible_to_signed_block_type auto y) {
            auto r = int_var(x);
            r += y;
            return r;
        }
        constexpr int_var operator+(int_view x, uint_const_t<>) { return int_var(x); }
        constexpr int_var operator+(int_view x, int_const_t<>) { return int_var(x); }

        constexpr int_var operator+(int_var&& x, uint_view y) { return std::move(x += y); }
        constexpr int_var operator+(int_var&& x, uint_var&& y) { return std::move(x += y); }
        constexpr int_var operator+(int_var&& x, int_view y) { return std::move(x += y); }
        constexpr int_var operator+(int_var&& x, int_var&& y) { return std::move(x += y); }
        constexpr int_var operator+(int_var&& x, convertible_to_block_type auto y) {
            return std::move(x += y);
        }
        constexpr int_var operator+(int_var&& x, convertible_to_signed_block_type auto y) {
            return std::move(x += y);
        }
        constexpr int_var operator+(int_var&& x, uint_const_t<>) { return std::move(x); }
        constexpr int_var operator+(int_var&& x, int_const_t<>) { return std::move(x); }

        constexpr int_var operator+(convertible_to_block_type auto x, int_view y) { return y + x; }
        constexpr int_var operator+(convertible_to_block_type auto x, int_var&& y) {
            return std::move(y += x);
        }

        constexpr int_var operator+(convertible_to_signed_block_type auto x, uint_view y) {
            return y + x;
        }
        constexpr int_var operator+(convertible_to_signed_block_type auto x, uint_var&& y) {
            return std::move(y) + x;
        }
        constexpr int_var operator+(convertible_to_signed_block_type auto x, int_view y) {
            return y + x;
        }
        constexpr int_var operator+(convertible_to_signed_block_type auto x, int_var&& y) {
            return std::move(y) + x;
        }

        constexpr int_var operator+(uint_const_t<>, int_view y) { return int_var(y); }
        constexpr int_var operator+(uint_const_t<>, int_var&& y) { return std::move(y); }

        constexpr int_var operator+(int_const_t<>, uint_view y) { return int_var(uint_var(y)); }
        constexpr int_var operator+(int_const_t<>, uint_var&& y) { return to_signed(std::move(y)); }
        constexpr int_var operator+(int_const_t<>, int_view y) { return int_var(y); }
        constexpr int_var operator+(int_const_t<>, int_var&& y) { return std::move(y); }

        constexpr int_var operator-(uint_view x, int_view y) {
            auto r = int_var(uint_var(x));
            r -= y;
            return r;
        }
        constexpr int_var operator-(uint_view x, int_var&& y) {
            y -= x;
            y.invert_sign();
            return std::move(y);
        }
        constexpr int_var operator-(uint_view x, convertible_to_signed_block_type auto y) {
            auto r = int_var(uint_var(x));
            r -= y;
            return r;
        }
        constexpr int_var operator-(uint_view x, int_const_t<>) { return int_var(uint_var(x)); }

        constexpr int_var operator-(uint_var&& x, int_view y) {
            auto r = to_signed(std::move(x));
            r -= y;
            return r;
        }
        constexpr int_var operator-(uint_var&& x, int_var&& y) {
            y -= x;
            y.invert_sign();
            return std::move(y);
        }
        constexpr int_var operator-(uint_var&& x, convertible_to_signed_block_type auto y) {
            auto r = to_signed(std::move(x));
            r -= y;
            return r;
        }
        constexpr int_var operator-(uint_var&& x, int_const_t<>) { return to_signed(std::move(x)); }

        constexpr int_var operator-(int_view x, uint_view y) {
            auto r = int_var(x);
            r -= y;
            return r;
        }
        constexpr int_var operator-(int_view x, uint_var&& y) {
            auto r = to_signed(std::move(y));
            r -= x;
            r.invert_sign();
            return r;
        }
        constexpr int_var operator-(int_view x, int_view y) {
            auto r = int_var(x);
            r -= y;
            return r;
        }
        constexpr int_var operator-(int_view x, int_var&& y) {
            y -= x;
            y.invert_sign();
            return std::move(y);
        }
        constexpr int_var operator-(int_view x, convertible_to_block_type auto y) {
            auto r = int_var(x);
            r -= y;
            return r;
        }
        constexpr int_var operator-(int_view x, convertible_to_signed_block_type auto y) {
            auto r = int_var(x);
            r -= y;
            return r;
        }
        constexpr int_var operator-(int_view x, uint_const_t<>) { return int_var(x); }
        constexpr int_var operator-(int_view x, int_const_t<>) { return int_var(x); }

        constexpr int_var operator-(int_var&& x, uint_view y) { return std::move(x -= y); }
        constexpr int_var operator-(int_var&& x, uint_var&& y) { return std::move(x -= y); }
        constexpr int_var operator-(int_var&& x, int_view y) { return std::move(x -= y); }
        constexpr int_var operator-(int_var&& x, int_var&& y) { return std::move(x -= y); }
        constexpr int_var operator-(int_var&& x, convertible_to_block_type auto y) {
            return std::move(x -= y);
        }
        constexpr int_var operator-(int_var&& x, convertible_to_signed_block_type auto y) {
            return std::move(x -= y);
        }
        constexpr int_var operator-(int_var&& x, uint_const_t<>) { return std::move(x); }
        constexpr int_var operator-(int_var&& x, int_const_t<>) { return std::move(x); }

        constexpr int_var operator-(convertible_to_block_type auto x, int_view y) {
            auto r = int_var(y);
            r -= x;
            r.invert_sign();
            return r;
        }
        constexpr int_var operator-(convertible_to_block_type auto x, int_var&& y) {
            y -= x;
            y.invert_sign();
            return std::move(y);
        }

        constexpr int_var operator-(convertible_to_signed_block_type auto x, uint_view y) {
            auto r = int_var(x);
            r -= y;
            return r;
        }
        constexpr int_var operator-(convertible_to_signed_block_type auto x, uint_var&& y) {
            auto r = to_signed(std::move(y));
            r -= x;
            r.invert_sign();
            return r;
        }
        constexpr int_var operator-(convertible_to_signed_block_type auto x, int_view y) {
            auto r = int_var(x);
            r -= y;
            return r;
        }
        constexpr int_var operator-(convertible_to_signed_block_type auto x, int_var&& y) {
            y -= x;
            y.invert_sign();
            return std::move(y);
        }

        constexpr int_var operator-(uint_const_t<>, int_view y) {
            auto r = int_var(y);
            r.invert_sign();
            return r;
        }
        constexpr int_var operator-(uint_const_t<>, int_var&& y) {
            y.invert_sign();
            return std::move(y);
        }

        constexpr int_var operator-(int_const_t<>, uint_view y) {
            return int_var(sign_t::negative, uint_var(y));
        }
        constexpr int_var operator-(int_const_t<>, uint_var&& y) {
            return to_negative(std::move(y));
        }
        constexpr int_var operator-(int_const_t<>, int_view y) {
            auto r = int_var(y);
            r.invert_sign();
            return r;
        }
        constexpr int_var operator-(int_const_t<>, int_var&& y) {
            y.invert_sign();
            return std::move(y);
        }

        constexpr int_var operator*(uint_view x, int_view y) {
            return int_var(y.sign(), x * y.abs());
        }
        constexpr int_var operator*(uint_view x, convertible_to_signed_block_type auto y) {
            return int_var(y >= 0 ? sign_t::positive : sign_t::negative,
                           x * (y >= 0 ? static_cast<block_type>(y)
                                       : (wuint::uint64_mask - static_cast<block_type>(y) + 1)));
        }
        constexpr int_var operator*(uint_view x, int_const_t<sign_t::positive, 1u>) {
            return int_var(sign_t::positive, uint_var(x));
        }
        constexpr int_var operator*(uint_view x, int_const_t<sign_t::negative, 1u>) {
            return int_var(sign_t::negative, uint_var(x));
        }
        constexpr int_var operator*(uint_view, int_const_t<>) { return {}; }

        constexpr int_var operator*(uint_var&& x, convertible_to_signed_block_type auto y) {
            auto r = to_signed(std::move(x));
            r *= y;
            return r;
        }
        constexpr int_var operator*(uint_var&& x, int_const_t<sign_t::positive, 1u>) {
            return int_var(sign_t::positive, std::move(x));
        }
        constexpr int_var operator*(uint_var&& x, int_const_t<sign_t::negative, 1u>) {
            return int_var(sign_t::negative, std::move(x));
        }
        constexpr int_var operator*(uint_var&&, int_const_t<>) { return {}; }

        constexpr int_var operator*(int_view x, uint_view y) { return y * x; }
        constexpr int_var operator*(int_view x, int_view y) {
            return int_var(x.sign() == y.sign() ? sign_t::positive : sign_t::negative,
                           x.abs() * y.abs());
        }
        constexpr int_var operator*(int_view x, convertible_to_block_type auto y) {
            auto r = int_var(x);
            r *= y;
            return r;
        }
        constexpr int_var operator*(int_view x, convertible_to_signed_block_type auto y) {
            auto r = int_var(x);
            r *= y;
            return r;
        }
        constexpr int_var operator*(int_view x, uint_const_t<1u>) { return int_var(x); }
        constexpr int_var operator*(int_view x, int_const_t<sign_t::positive, 1u>) {
            return int_var(x);
        }
        constexpr int_var operator*(int_view x, int_const_t<sign_t::negative, 1u>) {
            return int_var(x).invert_sign();
        }
        constexpr int_var operator*(int_view, uint_const_t<>) { return {}; }
        constexpr int_var operator*(int_view, int_const_t<>) { return {}; }

        constexpr int_var operator*(int_var&& x, convertible_to_block_type auto y) {
            return std::move(x *= y);
        }
        constexpr int_var operator*(int_var&& x, convertible_to_signed_block_type auto y) {
            return std::move(x *= y);
        }
        constexpr int_var operator*(int_var&& x, uint_const_t<1u>) { return std::move(x); }
        constexpr int_var operator*(int_var&& x, int_const_t<sign_t::positive, 1u>) {
            return std::move(x);
        }
        constexpr int_var operator*(int_var&& x, int_const_t<sign_t::negative, 1u>) {
            return std::move(x.invert_sign());
        }

        constexpr int_var operator*(convertible_to_block_type auto x, int_view y) { return y * x; }
        constexpr int_var operator*(convertible_to_block_type auto x, int_var&& y) {
            return std::move(y) * x;
        }

        constexpr int_var operator*(convertible_to_signed_block_type auto x, uint_view y) {
            return y * x;
        }
        constexpr int_var operator*(convertible_to_signed_block_type auto x, uint_var&& y) {
            return std::move(y) * x;
        }
        constexpr int_var operator*(convertible_to_signed_block_type auto x, int_view y) {
            return y * x;
        }
        constexpr int_var operator*(convertible_to_signed_block_type auto x, int_var&& y) {
            return std::move(y) * x;
        }

        constexpr int_var operator*(uint_const_t<1u>, int_view y) { return int_var(y); }
        constexpr int_var operator*(uint_const_t<1u>, int_var&& y) { return std::move(y); }

        constexpr int_var operator*(uint_const_t<>, int_view) { return {}; }

        constexpr int_var operator*(int_const_t<sign_t::positive, 1u>, uint_view y) {
            return int_var(sign_t::positive, uint_var(y));
        }
        constexpr int_var operator*(int_const_t<sign_t::positive, 1u>, uint_var&& y) {
            return int_var(sign_t::positive, std::move(y));
        }
        constexpr int_var operator*(int_const_t<sign_t::positive, 1u>, int_view y) {
            return int_var(y);
        }
        constexpr int_var operator*(int_const_t<sign_t::positive, 1u>, int_var&& y) {
            return std::move(y);
        }

        constexpr int_var operator*(int_const_t<sign_t::negative, 1u>, uint_view y) {
            return int_var(sign_t::negative, uint_var(y));
        }
        constexpr int_var operator*(int_const_t<sign_t::negative, 1u>, uint_var&& y) {
            return int_var(sign_t::negative, std::move(y));
        }
        constexpr int_var operator*(int_const_t<sign_t::negative, 1u>, int_view y) {
            return int_var(y).invert_sign();
        }
        constexpr int_var operator*(int_const_t<sign_t::negative, 1u>, int_var&& y) {
            return std::move(y.invert_sign());
        }

        constexpr int_var operator*(int_const_t<>, uint_view) { return {}; }
        constexpr int_var operator*(int_const_t<>, int_view) { return {}; }

        constexpr int_var& operator*=(int_var& x, uint_view y) {
            auto r = x * y;
            x = std::move(r);
            return x;
        }
        constexpr int_var& operator*=(int_var& x, int_view y) {
            auto r = x * y;
            x = std::move(r);
            return x;
        }

        // We do not define division by signed integers, since it has ambiguous semantics.

        constexpr int_var operator/(int_view x, uint_view y) {
            auto r = int_var(x);
            return r.long_division(y);
        }
        constexpr int_var operator/(int_var&& x, uint_view y) { return x.long_division(y); }
        constexpr int_var operator/(int_view x, convertible_to_block_type auto y) {
            auto r = int_var(x);
            return r.long_division(uint_view::make_view_from_single_block(y));
        }
        constexpr int_var operator/(int_var&& x, convertible_to_block_type auto y) {
            return x.long_division(uint_view::make_view_from_single_block(y));
        }
        constexpr int_var operator/(int_view x, uint_const_t<1u>) { return int_var(x); }
        constexpr int_var operator/(int_var&& x, uint_const_t<1u>) { return std::move(x); }
        constexpr int_var operator/(int_view x, uint_const_t<>) = delete;

        constexpr int_var& operator/=(int_var& x, uint_view y) {
            auto r = x / y;
            x = std::move(r);
            return x;
        }
        constexpr int_var& operator/=(int_var& x, convertible_to_block_type auto y) {
            auto r = x / y;
            x = std::move(r);
            return x;
        }
        constexpr int_var& operator/=(int_var& x, uint_const_t<1u>) { return x; }
        constexpr int_var& operator/=(int_var& x, uint_const_t<>) = delete;

        constexpr int_var operator%(int_view x, uint_view y) {
            auto r = int_var(x);
            r.long_division(y);
            return r;
        }
        constexpr int_var operator%(int_var&& x, uint_view y) {
            x.long_division(y);
            return std::move(x);
        }
        constexpr int_var operator%(int_view x, convertible_to_block_type auto y) {
            auto r = int_var(x);
            r.long_division(y);
            return r;
        }
        constexpr int_var operator%(int_var&& x, convertible_to_block_type auto y) {
            x.long_division(uint_view::make_view_from_single_block(y));
            return std::move(x);
        }
        constexpr int_var operator%(int_view, uint_const_t<1u>) { return {}; }
        constexpr int_var operator%(int_view, uint_const_t<>) = delete;

        using int_var_div_t = util::div_t<int_var, uint_var>;
        constexpr int_var_div_t div(int_view x, uint_view y) {
            int_var_div_t ret;
            auto temp = int_var(x);
            ret.quot = temp.long_division(y);
            ret.rem = std::move(temp).abs();
            return ret;
        }
        constexpr int_var_div_t div(int_var&& x, uint_view y) {
            int_var_div_t ret;
            ret.quot = x.long_division(y);
            ret.rem = std::move(x).abs();
            return ret;
        }
        constexpr int_var_div_t div(int_view x, uint_const_t<1u>) { return {int_var(x), {}}; }
        constexpr int_var_div_t div(int_var&& x, uint_const_t<1u>) { return {std::move(x), {}}; }
        constexpr int_var_div_t div(int_view, uint_const_t<>) = delete;

        constexpr int_var div_floor(uint_view x, int_view y) {
            if (y.is_nonnegative()) {
                return int_var(sign_t::positive, div_floor(x, y.abs()));
            }
            else {
                return int_var(sign_t::negative, div_ceil(x, y.abs()));
            }
        }
        constexpr int_var div_floor(uint_var&& x, int_view y) {
            if (y.is_nonnegative()) {
                return int_var(sign_t::positive, div_floor(std::move(x), y.abs()));
            }
            else {
                return int_var(sign_t::negative, div_ceil(std::move(x), y.abs()));
            }
        }
        constexpr int_var div_floor(uint_view x, int_const_t<sign_t::positive, 1u>) {
            return int_var{x};
        }
        constexpr int_var div_floor(uint_view x, int_const_t<sign_t::negative, 1u>) {
            return int_var{sign_t::negative, x};
        }
        constexpr int_var div_floor(uint_var&& x, int_const_t<sign_t::positive, 1u>) {
            return int_var{std::move(x)};
        }
        constexpr int_var div_floor(uint_var&& x, int_const_t<sign_t::negative, 1u>) {
            return int_var{sign_t::negative, std::move(x)};
        }
        constexpr int_var div_floor(uint_view x, int_const_t<>) = delete;

        constexpr int_var div_floor(int_view x, uint_view y) { return x / y; }
        constexpr int_var div_floor(int_var&& x, uint_view y) { return x.long_division(y); }
        constexpr int_var div_floor(int_view x, uint_const_t<1u>) { return int_var(x); }
        constexpr int_var div_floor(int_var&& x, uint_const_t<1u>) { return std::move(x); }
        constexpr int_var div_floor(int_view x, uint_const_t<>) = delete;

        constexpr int_var div_ceil(uint_view x, int_view y) {
            if (y.is_nonnegative()) {
                return int_var(sign_t::positive, div_ceil(x, y.abs()));
            }
            else {
                return int_var(sign_t::negative, div_floor(x, y.abs()));
            }
        }
        constexpr int_var div_ceil(uint_var&& x, int_view y) {
            if (y.is_nonnegative()) {
                return int_var(sign_t::positive, div_ceil(std::move(x), y.abs()));
            }
            else {
                return int_var(sign_t::negative, div_floor(std::move(x), y.abs()));
            }
        }
        constexpr int_var div_ceil(uint_view x, int_const_t<sign_t::positive, 1u>) {
            return int_var{x};
        }
        constexpr int_var div_ceil(uint_view x, int_const_t<sign_t::negative, 1u>) {
            return int_var{sign_t::negative, x};
        }
        constexpr int_var div_ceil(uint_var&& x, int_const_t<sign_t::positive, 1u>) {
            return int_var{std::move(x)};
        }
        constexpr int_var div_ceil(uint_var&& x, int_const_t<sign_t::negative, 1u>) {
            return int_var{sign_t::negative, std::move(x)};
        }
        constexpr int_var div_ceil(uint_view x, int_const_t<>) = delete;

        constexpr int_var div_ceil(int_view x, uint_view y) {
            if (x.is_nonnegative()) {
                return int_var{sign_t::positive, div_ceil(x.abs(), y)};
            }
            else {
                return int_var{sign_t::negative, div_floor(x.abs(), y)};
            }
        }
        constexpr int_var div_ceil(int_var&& x, uint_view y) {
            if (x.is_nonnegative()) {
                return int_var{sign_t::positive, div_ceil(std::move(x).abs(), y)};
            }
            else {
                return int_var{sign_t::negative, div_floor(std::move(x).abs(), y)};
            }
        }
        constexpr int_var div_ceil(int_view x, uint_const_t<1u>) { return int_var(x); }
        constexpr int_var div_ceil(int_var&& x, uint_const_t<1u>) { return std::move(x); }
        constexpr int_var div_ceil(int_view x, uint_const_t<>) = delete;

        constexpr int_var div_floor(int_view x, int_view y) {
            if (x.sign() == y.sign()) {
                return int_var{sign_t::positive, div_floor(x.abs(), y.abs())};
            }
            else {
                return int_var{sign_t::negative, div_ceil(x.abs(), y.abs())};
            }
        }
        constexpr int_var div_floor(int_var&& x, int_view y) {
            if (x.sign() == y.sign()) {
                return int_var{sign_t::positive, div_floor(std::move(x).abs(), y.abs())};
            }
            else {
                return int_var{sign_t::negative, div_ceil(std::move(x).abs(), y.abs())};
            }
        }
        constexpr int_var div_floor(int_view x, int_const_t<sign_t::positive, 1u>) {
            return int_var{x};
        }
        constexpr int_var div_floor(int_var&& x, int_const_t<sign_t::positive, 1u>) {
            return std::move(x);
        }
        constexpr int_var div_floor(int_view x, int_const_t<sign_t::negative, 1u>) {
            return int_var{x}.invert_sign();
        }
        constexpr int_var div_floor(int_var&& x, int_const_t<sign_t::negative, 1u>) {
            return std::move(x.invert_sign());
        }
        constexpr int_var div_floor(int_view x, int_const_t<>) = delete;

        constexpr int_var div_ceil(int_view x, int_view y) {
            if (x.sign() == y.sign()) {
                return int_var{sign_t::positive, div_ceil(x.abs(), y.abs())};
            }
            else {
                return int_var{sign_t::negative, div_floor(x.abs(), y.abs())};
            }
        }
        constexpr int_var div_ceil(int_var&& x, int_view y) {
            if (x.sign() == y.sign()) {
                return int_var{sign_t::positive, div_ceil(std::move(x).abs(), y.abs())};
            }
            else {
                return int_var{sign_t::negative, div_floor(std::move(x).abs(), y.abs())};
            }
        }
        constexpr int_var div_ceil(int_view x, int_const_t<sign_t::positive, 1u>) {
            return int_var{x};
        }
        constexpr int_var div_ceil(int_var&& x, int_const_t<sign_t::positive, 1u>) {
            return std::move(x);
        }
        constexpr int_var div_ceil(int_view x, int_const_t<sign_t::negative, 1u>) {
            return int_var{x}.invert_sign();
        }
        constexpr int_var div_ceil(int_var&& x, int_const_t<sign_t::negative, 1u>) {
            return std::move(x.invert_sign());
        }
        constexpr int_var div_ceil(int_view x, int_const_t<>) = delete;

        // Find the largest power of 2 dividing n, divide n by that power, and return
        // the exponent.
        constexpr std::size_t factor_out_power_of_2(int_var& n) noexcept {
            return n.factor_out_power_of_2();
        }

        // Operations on int_const.
        namespace detail {
            // Used for converting int_var into int_const.
            // Used for converting uint_var into uint_const.
            template <std::size_t N>
            struct block_holder_size_sign_triple : public block_holder_size_pair<N> {
                sign_t sign;

                constexpr block_holder_size_sign_triple(int_view n) noexcept
                    : block_holder_size_pair<N>(n.abs()), sign(n.sign()) {}
            };

            template <auto triple>
            using block_holder_size_sign_triple_to_int_const =
                int_const_impl<triple.sign, slice<triple.size>(triple.arr)>;

            template <sign_t s, std::size_t M, static_block_holder<M> x>
            constexpr auto operator-(int_const_impl<s, x>) noexcept {
                return int_const_impl<
                    (x.size() == 0 ? sign_t::positive : ::jkj::bigint::invert_sign(s)), x>{};
            }

            template <sign_t s, sign_t t, std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator+(int_const_impl<s, x> xx, int_const_impl<t, y> yy) noexcept {
                if constexpr (M == 0) {
                    return yy;
                }
                else {
                    constexpr auto result = block_holder_size_sign_triple<util::max(M, N) + 1>{
                        to_view(decltype(xx){}) + decltype(yy){}};
                    return block_holder_size_sign_triple_to_int_const<result>{};
                }
            }

            template <sign_t s, sign_t t, std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator-(int_const_impl<s, x> xx, int_const_impl<t, y> yy) noexcept {
                if constexpr (M == 0) {
                    return int_const_t<>{};
                }
                else {
                    constexpr auto result = block_holder_size_sign_triple<util::max(M, N) + 1>{
                        to_view(decltype(xx){}) - decltype(yy){}};
                    return block_holder_size_sign_triple_to_int_const<result>{};
                }
            }

            template <sign_t s, sign_t t, std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator*(int_const_impl<s, x> xx, int_const_impl<t, y> yy) noexcept {
                if constexpr (M == 0) {
                    return int_const_t<>{};
                }
                else {
                    constexpr auto result = block_holder_size_sign_triple<M + N>{
                        to_view(decltype(xx){}) * decltype(yy){}};
                    return block_holder_size_sign_triple_to_int_const<result>{};
                }
            }

            template <class Quotient, class Remainder>
            struct int_const_div_t;

            template <sign_t s, std::size_t M, std::size_t N, static_block_holder<M> quotient,
                      static_block_holder<N> remainder>
            struct int_const_div_t<int_const_impl<s, quotient>, uint_const_impl<remainder>> {
                int_const_impl<s, quotient> quot;
                uint_const_impl<remainder> rem;
            };

            template <sign_t s, std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
                requires(N != 0)
            constexpr auto div(int_const_impl<s, x> xx, uint_const_impl<y> yy) noexcept {
                if constexpr (M == 0) {
                    return int_const_div_t<int_const_t<>, uint_const_t<>>{};
                }
                else {
                    constexpr auto result = [] {
                        auto result_var = div(to_view(decltype(xx){}), to_view(decltype(yy){}));
                        auto quot =
                            block_holder_size_sign_triple<util::max(M, N) - util::min(M, N) + 2>(
                                result_var.quot);
                        auto rem = block_holder_size_pair<N>(result_var.rem);

                        struct intermediate_result {
                            decltype(quot) quot_;
                            decltype(rem) rem_;
                        };
                        return intermediate_result{quot, rem};
                    }();
                    return int_const_div_t<block_holder_size_sign_triple_to_int_const<result.quot_>,
                                           block_holder_size_pair_to_uint_const<result.rem_>>{};
                }
            }

            template <sign_t s, std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
                requires(N != 0)
            constexpr auto operator/(int_const_impl<s, x> xx, uint_const_impl<y> yy) noexcept {
                return div(xx, yy).quot;
            }

            template <sign_t s, std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
                requires(N != 0)
            constexpr auto operator%(int_const_impl<s, x> xx, uint_const_impl<y> yy) noexcept {
                return div(xx, yy).rem;
            }

            template <sign_t s, std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator<<(int_const_impl<s, x> xx, uint_const_impl<y>) noexcept {
                static_assert(N <= 1, "jkj::idiv: shift amount too large");
                if constexpr (N == 0) {
                    return xx;
                }
                else {
                    constexpr auto result =
                        block_holder_size_sign_triple<M + (y[0] + number_of_bits_in_block - 1) /
                                                              number_of_bits_in_block>{
                            to_view(decltype(xx){}) << y[0]};
                    return block_holder_size_sign_triple_to_int_const<result>{};
                }
            }

            template <sign_t s, std::size_t M, std::size_t N, static_block_holder<M> x,
                      static_block_holder<N> y>
            constexpr auto operator>>(int_const_impl<s, x> xx, uint_const_impl<y>) noexcept {
                static_assert(N <= 1, "jkj::idiv: shift amount too large");
                if constexpr (N == 0) {
                    return xx;
                }
                else {
                    constexpr auto result =
                        block_holder_size_sign_triple<M>{to_view(decltype(xx){}) >> y[0]};
                    return block_holder_size_sign_triple_to_int_const<result>{};
                }
            }
        }

        // Specializations for swap.
        constexpr void swap(uint_var& x, uint_var& y) noexcept { util::swap(x, y); }
        constexpr void swap(int_var& x, int_var& y) noexcept { util::swap(x, y); }
    }
}

#endif
