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


#ifndef JKJ_HEADER_IDIV_WUINT
#define JKJ_HEADER_IDIV_WUINT

#include "core.h"
#include <cstdint>

namespace jkj {
    namespace wuint {
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

            std::uint64_t high_;
            std::uint64_t low_;

            constexpr uint128(std::uint64_t high, std::uint64_t low) noexcept
                : high_{high}, low_{low} {}

            constexpr std::uint64_t high() const noexcept { return high_; }
            constexpr std::uint64_t low() const noexcept { return low_; }

            constexpr uint128& operator+=(std::uint64_t n) & noexcept {
                [[maybe_unused]] auto const generic_impl = [&] {
                    auto sum = low_ + n;
                    high_ += (sum < low_ ? 1 : 0);
                    low_ = sum;
                };
                JKJ_IF_CONSTEVAL {
                    generic_impl();
                    return *this;
                }
#if JKJ_HAS_BUILTIN(__builtin_addcll)
                unsigned long long carry{};
                low_ = __builtin_addcll(low_, n, 0, &carry);
                high_ = __builtin_addcll(high_, 0, carry, &carry);
#elif JKJ_HAS_BUILTIN(__builtin_ia32_addcarryx_u64)
                unsigned long long result{};
                auto carry = __builtin_ia32_addcarryx_u64(0, low_, n, &result);
                low_ = result;
                __builtin_ia32_addcarryx_u64(carry, high_, 0, &result);
                high_ = result;
#elif defined(_MSC_VER) && defined(_M_X64)
                auto carry = _addcarry_u64(0, low_, n, &low_);
                _addcarry_u64(carry, high_, 0, &high_);
#else
                generic_impl();
#endif
                return *this;
            }
        };

        constexpr inline std::uint64_t umul64(std::uint32_t x, std::uint32_t y) noexcept {
#if defined(_MSC_VER) && defined(_M_IX86)
            JKJ_IF_NOT_CONSTEVAL { return __emulu(x, y); }
#else
            return x * std::uint64_t(y);
#endif
        }

        // Get 128-bit result of multiplication of two 64-bit unsigned integers.
        constexpr inline uint128 umul128(std::uint64_t x, std::uint64_t y) noexcept {
            [[maybe_unused]] auto generic_impl = [&]() -> uint128 {
                auto a = std::uint32_t(x >> 32);
                auto b = std::uint32_t(x);
                auto c = std::uint32_t(y >> 32);
                auto d = std::uint32_t(y);

                auto ac = umul64(a, c);
                auto bc = umul64(b, c);
                auto ad = umul64(a, d);
                auto bd = umul64(b, d);

                auto intermediate = (bd >> 32) + std::uint32_t(ad) + std::uint32_t(bc);

                return {ac + (intermediate >> 32) + (ad >> 32) + (bc >> 32),
                        (intermediate << 32) + std::uint32_t(bd)};
            };
            JKJ_IF_CONSTEVAL { return generic_impl(); }
#if defined(__SIZEOF_INT128__)
            auto result = builtin_uint128_t(x) * builtin_uint128_t(y);
            return {std::uint64_t(result >> 64), std::uint64_t(result)};
#elif defined(_MSC_VER) && defined(_M_X64)
            uint128 result;
            result.low_ = _umul128(x, y, &result.high_);
            return result;
#else
            return generic_impl();
#endif
        }

        constexpr inline std::uint64_t umul128_upper64(std::uint64_t x, std::uint64_t y) noexcept {
            [[maybe_unused]] auto generic_impl = [&]() -> std::uint64_t {
                auto a = std::uint32_t(x >> 32);
                auto b = std::uint32_t(x);
                auto c = std::uint32_t(y >> 32);
                auto d = std::uint32_t(y);

                auto ac = umul64(a, c);
                auto bc = umul64(b, c);
                auto ad = umul64(a, d);
                auto bd = umul64(b, d);

                auto intermediate = (bd >> 32) + std::uint32_t(ad) + std::uint32_t(bc);

                return ac + (intermediate >> 32) + (ad >> 32) + (bc >> 32);
            };
            JKJ_IF_CONSTEVAL { return generic_impl(); }
#if defined(__SIZEOF_INT128__)
            auto result = builtin_uint128_t(x) * builtin_uint128_t(y);
            return std::uint64_t(result >> 64);
#elif defined(_MSC_VER) && defined(_M_X64)
            return __umulh(x, y);
#else
            return generic_impl();
#endif
        }
    }
}

#endif