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

#ifndef JKJ_HEADER_IDIV_WUINT
#define JKJ_HEADER_IDIV_WUINT

#include "util.h"
#include <cstdint>
#include <limits>

namespace jkj {
    namespace wuint {
        inline constexpr std::uint_least64_t uint64_mask = UINT64_C(0xffff'ffff'ffff'ffff);
        inline constexpr std::uint_least32_t uint32_mask = UINT32_C(0xffff'ffff);

        // Add-with-carry implementation assuming 64-bit words.
        // std::uint_least64_t may not be of 64-bit, but this should work anyway.
        constexpr std::uint_least64_t add_carry64(std::uint_least64_t a, std::uint_least64_t b,
                                                  unsigned int& carry) noexcept {
            [[maybe_unused]] auto generic_impl = [&] {
                std::uint_least64_t with_carry = ((a & uint64_mask) + carry) & uint64_mask;
                unsigned int first_carry = (with_carry < (a & uint64_mask)) ? 1 : 0;

                with_carry += (b & uint64_mask);
                with_carry &= uint64_mask;
                carry = first_carry | ((with_carry < (b & uint64_mask)) ? 1 : 0);

                return with_carry;
            };

            JKJ_IF_NOT_CONSTEVAL {
                if constexpr (std::numeric_limits<std::uint_least64_t>::digits == 64) {
#if JKJ_HAS_BUILTIN(__builtin_addcll) || JKJ_HAS_BUILTIN(__builtin_addcl)
    #if JKJ_HAS_BUILTIN(__builtin_addcll)
                    if constexpr (std::is_same_v<std::uint_least64_t, unsigned long long>) {
                        unsigned long long carryout{};
                        auto result = __builtin_addcll(a, b, carry, &carryout);
                        carry = static_cast<unsigned int>(carryout);
                        return result;
                    }
    #endif
    #if JKJ_HAS_BUILTIN(__builtin_addcl)
                    if constexpr (std::is_same_v<std::uint_least64_t, unsigned long>) {
                        unsigned long carryout{};
                        auto result = __builtin_addcl(a, b, carry, &carryout);
                        carry = static_cast<unsigned int>(carryout);
                        return result;
                    }
    #endif
#elif JKJ_HAS_BUILTIN(__builtin_ia32_addcarryx_u64)
                    unsigned long long result{};
                    auto carryout = __builtin_ia32_addcarryx_u64(static_cast<unsigned char>(carry),
                                                                 a, b, &result);
                    carry = static_cast<unsigned int>(carryout);
                    return result;
#elif JKJ_HAS_BUILTIN(__builtin_ia32_addcarry_u64)
                    unsigned long long result{};
                    auto carryout = __builtin_ia32_addcarry_u64(static_cast<unsigned char>(carry),
                                                                a, b, &result);
                    carry = static_cast<unsigned int>(carryout);
                    return result;
#elif defined(_MSC_VER) && defined(_M_X64)
                    auto carryout = _addcarry_u64(static_cast<unsigned char>(carry), a, b, &a);
                    carry = static_cast<unsigned int>(carryout);
                    return a;
#endif
                }
            }
            return generic_impl();
        }

        // Subtract-with-borrow implementation assuming 64-bit words.
        // std::uint_least64_t may not be of 64-bit, but this should work anyway.
        constexpr std::uint_least64_t sub_borrow64(std::uint_least64_t a, std::uint_least64_t b,
                                                   unsigned int& borrow) noexcept {
            [[maybe_unused]] auto generic_impl = [&] {
                std::uint_least64_t with_borrow = ((a & uint64_mask) - borrow) & uint64_mask;
                unsigned int first_borrow = (with_borrow > (a & uint64_mask)) ? 1 : 0;

                a = with_borrow - (b & uint64_mask);
                a &= uint64_mask;
                borrow = first_borrow | ((a > with_borrow) ? 1 : 0);

                return a;
            };

            JKJ_IF_NOT_CONSTEVAL {
                if constexpr (std::numeric_limits<std::uint_least64_t>::digits == 64) {
#if JKJ_HAS_BUILTIN(__builtin_subcll) || JKJ_HAS_BUILTIN(__builtin_subcl)
    #if JKJ_HAS_BUILTIN(__builtin_subcll)
                    if constexpr (std::is_same_v<std::uint_least64_t, unsigned long long>) {
                        unsigned long long borrowout{};
                        auto result = __builtin_subcll(a, b, borrow, &borrowout);
                        borrow = static_cast<unsigned int>(borrowout);
                        return result;
                    }
    #endif
    #if JKJ_HAS_BUILTIN(__builtin_subcl)
                    if constexpr (std::is_same_v<std::uint_least64_t, unsigned long>) {
                        unsigned long borrowout{};
                        auto result = __builtin_subcl(a, b, borrow, &borrowout);
                        borrow = static_cast<unsigned int>(borrowout);
                        return result;
                    }
    #endif
#elif JKJ_HAS_BUILTIN(__builtin_ia32_subborrow_u64)
                    unsigned long long result{};
                    auto borrowout = __builtin_ia32_subborrow_u64(
                        static_cast<unsigned char>(borrow), a, b, &result);
                    borrow = static_cast<unsigned int>(borrowout);
                    return result;
#elif defined(_MSC_VER) && defined(_M_X64)
                    auto borrowout = _subborrow_u64(static_cast<unsigned char>(borrow), a, b, &a);
                    borrow = static_cast<unsigned int>(borrowout);
                    return a;
#endif
                }
            }
            return generic_impl();
        }

        // Compilers might support built-in 128-bit integer types. However, it seems
        // that emulating them with a pair of 64-bit integers actually produces a better
        // code, so we avoid using those built-ins. That said, they are still useful for
        // implementing 64-bit x 64-bit -> 128-bit multiplication.

        // clang-format off
#if defined(__SIZEOF_INT128__)
		// To silence "error: ISO C++ does not support '__int128' for 'type name'
		// [-Wpedantic]"
#if defined(__GNUC__)
			__extension__
#endif
				using builtin_uint128_t = unsigned __int128;
#endif
        // clang-format on

        struct uint128 {
            uint128() = default;

            constexpr uint128(std::uint_least64_t high, std::uint_least64_t low) noexcept
                : high_{high & uint64_mask}, low_{low & uint64_mask} {}

            constexpr std::uint_least64_t high() const noexcept { return high_; }
            constexpr std::uint_least64_t low() const noexcept { return low_; }

            constexpr uint128& operator+=(std::uint_least64_t n) & noexcept {
                unsigned int carry = 0;
                low_ = add_carry64(low_, n, carry);
                high_ = add_carry64(high_, 0, carry);
                return *this;
            }

        private:
            std::uint_least64_t high_;
            std::uint_least64_t low_;
        };

        constexpr inline std::uint_least64_t umul64(std::uint_least32_t x,
                                                    std::uint_least32_t y) noexcept {
#if defined(_MSC_VER) && defined(_M_IX86)
            JKJ_IF_NOT_CONSTEVAL { return __emulu(x, y); }
#else
            return (x * std::uint_least64_t(y)) & uint64_mask;
#endif
        }

        // Get 128-bit result of multiplication of two 64-bit unsigned integers.
        constexpr inline uint128 umul128(std::uint_least64_t x, std::uint_least64_t y) noexcept {
            [[maybe_unused]] auto generic_impl = [&]() -> uint128 {
                auto a = std::uint_least32_t((x >> 32) & uint32_mask);
                auto b = std::uint_least32_t(x & uint32_mask);
                auto c = std::uint_least32_t((y >> 32) & uint32_mask);
                auto d = std::uint_least32_t(y & uint32_mask);

                auto ac = umul64(a, c);
                auto bc = umul64(b, c);
                auto ad = umul64(a, d);
                auto bd = umul64(b, d);

                auto intermediate = (bd >> 32) + (ad & uint32_mask) + (bc & uint32_mask);

                return {ac + (intermediate >> 32) + (ad >> 32) + (bc >> 32),
                        (intermediate << 32) + (bd & uint32_mask)};
            };
            JKJ_IF_CONSTEVAL { return generic_impl(); }
#if defined(__SIZEOF_INT128__)
            auto result = builtin_uint128_t(x & uint64_mask) * builtin_uint128_t(y & uint64_mask);
            return {std::uint_least64_t(result >> 64), std::uint_least64_t(result) & uint64_mask};
#elif defined(_MSC_VER) && defined(_M_X64)
            std::uint_least64_t high;
            auto low = _umul128(x, y, &high);
            return {high, low};
#else
            return generic_impl();
#endif
        }

        constexpr inline std::uint_least64_t umul128_upper64(std::uint_least64_t x,
                                                             std::uint_least64_t y) noexcept {
            [[maybe_unused]] auto generic_impl = [&]() -> std::uint_least64_t {
                auto a = std::uint_least32_t((x >> 32) & uint32_mask);
                auto b = std::uint_least32_t(x & uint32_mask);
                auto c = std::uint_least32_t((y >> 32) & uint32_mask);
                auto d = std::uint_least32_t(y & uint32_mask);

                auto ac = umul64(a, c);
                auto bc = umul64(b, c);
                auto ad = umul64(a, d);
                auto bd = umul64(b, d);

                auto intermediate = (bd >> 32) + (ad & uint32_mask) + (bc & uint32_mask);

                return ac + (intermediate >> 32) + (ad >> 32) + (bc >> 32);
            };
            JKJ_IF_CONSTEVAL { return generic_impl(); }
#if defined(__SIZEOF_INT128__)
            auto result = builtin_uint128_t(x & uint64_mask) * builtin_uint128_t(y & uint64_mask);
            return std::uint_least64_t(result >> 64);
#elif defined(_MSC_VER) && defined(_M_X64)
            return __umulh(x, y);
#else
            return generic_impl();
#endif
        }
    }
}

#endif
