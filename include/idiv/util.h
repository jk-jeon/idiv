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


#ifndef JKJ_HEADER_IDIV_CORE
#define JKJ_HEADER_IDIV_CORE

#include <cassert>
#include <cstddef>
#include <type_traits>
#include <utility>

#if defined(__cpp_if_consteval) && __cpp_is_consteval >= 202106L
    #define JKJ_IF_CONSTEVAL if consteval
    #define JKJ_IF_NOT_CONSTEVAL if !consteval
    #define JKJ_CAN_BRANCH_ON_CONSTEVAL 1
#elif defined(__cpp_lib_is_constant_evaluated) && __cpp_lib_is_constant_evaluated >= 201811L
    #define JKJ_IF_CONSTEVAL if (std::is_constant_evaluated())
    #define JKJ_IF_NOT_CONSTEVAL if (!std::is_constant_evaluated())
    #define JKJ_CAN_BRANCH_ON_CONSTEVAL 1
#else
    #define JKJ_IF_CONSTEVAL if constexpr (false)
    #define JKJ_IF_NOT_CONSTEVAL if constexpr (true)
    #define JKJ_CAN_BRANCH_ON_CONSTEVAL 0
#endif

#if defined(__has_builtin)
    #define JKJ_HAS_BUILTIN(x) __has_builtin(x)
#else
    #define JKJ_HAS_BUILTIN(x) false
#endif

#if defined(_MSC_VER)
    #include <intrin.h>
#endif

namespace jkj {
    namespace util {
        // For compile-time error-handling.
        // Note that this function cannot be called in constexpr context.
        template <char const* error_msg>
        void assert_failed() {}

        // A replacement for assert to use in constexpr functions.
        template <char const* error_msg>
        constexpr void constexpr_assert(bool cond) {
            // Use the assert macro for the runtime.
            assert(cond);
            // Call a non-constexpr function to trigger an error for the compile-time.
            if (!cond) {
                assert_failed<error_msg>();
            }
        }

        namespace error_msgs {
            inline constexpr char no_error_msg[] = "no error message";
            inline constexpr char index_out_of_range[] = "index out of range";
            inline constexpr char overflow[] = "overflow";
            inline constexpr char underflow[] = "underflow";
            inline constexpr char divide_by_zero[] = "divide by zero";
        }

        template <class T>
        constexpr T const& min(T const& x, T const& y) {
            return x >= y ? y : x;
        }

        template <class T>
        constexpr T const& max(T const& x, T const& y) {
            return x >= y ? x : y;
        }

        // A minimal implementation of std::array.
        template <class T, std::size_t N>
        struct array {
            using value_type = T;

            constexpr value_type* data() noexcept { return data_; }
            constexpr value_type const* data() const noexcept { return data_; }

            constexpr value_type& front() noexcept { return data_[0]; }
            constexpr value_type const& front() const noexcept { return data_[0]; }

            constexpr value_type& back() noexcept { return data_[N - 1]; }
            constexpr value_type const& back() const noexcept { return data_[N - 1]; }

            constexpr value_type* begin() noexcept { return data_; }
            constexpr value_type const* begin() const noexcept { return data_; }
            constexpr value_type const* cbegin() const noexcept { return data_; }

            constexpr value_type* end() noexcept { return data_ + N; }
            constexpr value_type const* end() const noexcept { return data_ + N; }
            constexpr value_type const* cend() const noexcept { return data_ + N; }

            constexpr value_type& operator[](std::size_t idx) & noexcept {
                constexpr_assert<error_msgs::index_out_of_range>(idx < N);
                return data_[idx];
            }
            constexpr value_type const& operator[](std::size_t idx) const& noexcept {
                constexpr_assert<error_msgs::index_out_of_range>(idx < N);
                return data_[idx];
            }

            static constexpr std::size_t size() noexcept { return N; }

            value_type data_[N];
        };

        // Zero-sized arrays.
        template <class T>
        struct array<T, 0> {
            using block_type = T;
            static constexpr std::size_t size() noexcept { return 0; }
        };

        // A minimal implementation of std::apply for array.
        namespace detail {
            template <class Functor, class Array, std::size_t... indices>
            constexpr decltype(auto) apply_impl(Functor&& f, Array&& arr,
                                                std::index_sequence<indices...>) noexcept {
                return static_cast<Functor&&>(f)(arr[indices]...);
            }
        }
        template <class Functor, class T, std::size_t N>
        constexpr decltype(auto) apply(Functor&& f, array<T, N> const& arr) noexcept {
            return detail::apply_impl(static_cast<Functor&&>(f), arr,
                                      std::make_index_sequence<N>{});
        }
        template <class Functor, class T>
        constexpr decltype(auto) apply(Functor&& f, array<T, 0> const&) noexcept {
            return static_cast<Functor&&>(f)();
        }

        // A minimal implementation of std::span.
        template <class T>
        class span {
        public:
            using value_type = T;

            constexpr span(value_type* ptr, std::size_t size) noexcept : ptr_{ptr}, size_{size} {}

            constexpr value_type* data() const noexcept { return ptr_; }

            constexpr value_type& front() const noexcept {
                constexpr_assert<error_msgs::index_out_of_range>(size_ != 0);
                return ptr_[0];
            }

            constexpr value_type& back() const noexcept {
                constexpr_assert<error_msgs::index_out_of_range>(size_ != 0);
                return ptr_[size_ - 1];
            }

            constexpr value_type const* begin() const noexcept { return ptr_; }
            constexpr value_type const* cbegin() const noexcept { return ptr_; }

            constexpr value_type const* end() const noexcept { return ptr_ + size_; }
            constexpr value_type const* cend() const noexcept { return ptr_ + size_; }

            constexpr value_type& operator[](std::size_t idx) const noexcept {
                constexpr_assert<error_msgs::index_out_of_range>(idx < size_);
                return ptr_[idx];
            }

            constexpr std::size_t size() const noexcept { return size_; }
            constexpr bool empty() const noexcept { return size_ == 0; }

        private:
            value_type* ptr_;
            std::size_t size_;
        };
    }
}

#endif