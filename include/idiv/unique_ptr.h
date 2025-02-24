// Copyright 2024 Junekey Jeon
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

#ifndef JKJ_HEADER_IDIV_UNIQUE_PTR
#define JKJ_HEADER_IDIV_UNIQUE_PTR

#include <cstddef>
#include <type_traits>

namespace jkj {
    namespace util {
        // Minimal implementation of std::unique_ptr since std::unique_ptr is not constexpr in
        // C++20. It does not support array and custom deleter since it is not used.
        template <class T>
            requires(!std::is_array_v<T>)
        class unique_ptr {
            T* ptr_;

            template <class U>
                requires(!std::is_array_v<U>)
            friend class unique_ptr;

        public:
            using pointer = T*;
            using element_type = T;

            constexpr unique_ptr() noexcept : ptr_{nullptr} {}
            constexpr unique_ptr(std::nullptr_t) noexcept : ptr_{nullptr} {}
            explicit constexpr unique_ptr(pointer ptr) noexcept : ptr_{ptr} {}
            constexpr unique_ptr(unique_ptr&& other) noexcept : ptr_{other.release()} {}
            template <class U>
                requires(!std::is_array_v<U> &&
                         std::is_convertible_v<typename unique_ptr<U>::pointer, pointer>)
            constexpr unique_ptr(unique_ptr<U>&& other) noexcept : ptr_{other.release()} {}

            constexpr unique_ptr& operator=(std::nullptr_t) noexcept {
                reset();
                return *this;
            }
            template <class U>
                requires(!std::is_array_v<U> &&
                         std::is_convertible_v<typename unique_ptr<U>::pointer, pointer>)
            constexpr unique_ptr& operator=(unique_ptr<U>&& other) noexcept {
                reset(other.release());
                return *this;
            }
            constexpr unique_ptr& operator=(unique_ptr&& other) noexcept {
                reset(other.release());
                return *this;
            }

            unique_ptr(const unique_ptr&) = delete;
            unique_ptr& operator=(const unique_ptr&) = delete;

            constexpr void swap(unique_ptr& other) noexcept {
                using std::swap;
                swap(ptr_, other.ptr_);
            }

            constexpr ~unique_ptr() noexcept { delete ptr_; }

            constexpr T& operator*() const noexcept { return *ptr_; }
            constexpr pointer operator->() const noexcept { return ptr_; }
            constexpr pointer get() const noexcept { return ptr_; }
            explicit constexpr operator bool() const noexcept { return ptr_ != nullptr; }

            constexpr pointer release() noexcept {
                auto tmp = ptr_;
                ptr_ = nullptr;
                return tmp;
            }

            constexpr void reset(pointer ptr = nullptr) noexcept {
                auto tmp = ptr_;
                ptr_ = ptr;
                delete tmp;
            }
        };

        template <class T, class... Args>
            requires(!std::is_array_v<T>)
        constexpr unique_ptr<T> make_unique(Args&&... args) {
            return unique_ptr<T>{new T(static_cast<Args&&>(args)...)};
        }
    }
}

#endif