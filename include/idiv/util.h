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


#ifndef JKJ_HEADER_IDIV_UTIL
#define JKJ_HEADER_IDIV_UTIL

#include <cassert>
#include <concepts>
#include <cstddef>
#include <limits>
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

        namespace error_msgs {
            inline constexpr char no_error_msg[] = "no error message";
            inline constexpr char index_out_of_range[] = "index out of range";
            inline constexpr char overflow[] = "overflow";
            inline constexpr char underflow[] = "underflow";
            inline constexpr char divide_by_zero[] = "divide by zero";
        }

        // A replacement for assert to use in constexpr functions.
        template <char const* error_msg = error_msgs::no_error_msg>
        constexpr void constexpr_assert(bool cond) {
            // Use the assert macro for the runtime.
            assert(cond);
            // Call a non-constexpr function to trigger an error for the compile-time.
            if (!cond) {
                assert_failed<error_msg>();
            }
        }

        template <class T>
        constexpr T const& min(T const& x, T const& y) {
            return x >= y ? y : x;
        }

        template <class T>
        constexpr T const& max(T const& x, T const& y) {
            return x >= y ? x : y;
        }

        template <class T>
        constexpr void swap(T& x, T& y) noexcept {
            auto temp = static_cast<T&&>(x);
            x = static_cast<T&&>(y);
            y = static_cast<T&&>(temp);
        }

        constexpr int strong_order_to_int(std::strong_ordering r) noexcept {
            if (r > 0) {
                return 1;
            }
            else if (r < 0) {
                return -1;
            }
            else {
                return 0;
            }
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

        // Some utilities for dealing with primitive integer types.
        enum class sign_t : bool { positive = false, negative = true };
        template <class MaybeSigned, class Remainder = MaybeSigned>
        struct div_t {
            MaybeSigned quot;
            Remainder rem;
        };
        namespace detail {
            // Conversion between unsigned entities and signed entities.
            constexpr auto to_signed(std::unsigned_integral auto n) noexcept {
                using unsigned_t = decltype(n);
                using signed_t = std::make_signed_t<unsigned_t>;
                // There should be no UB here.
                util::constexpr_assert(
                    n <= static_cast<unsigned_t>(std::numeric_limits<signed_t>::max()));
                return static_cast<signed_t>(n);
            }
            constexpr auto to_signed(std::signed_integral auto n) noexcept { return n; }
            struct to_signed_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return to_signed(static_cast<decltype(n)&&>(n));
                }
            };

            constexpr auto to_negative(std::unsigned_integral auto n) noexcept {
                using unsigned_t = decltype(n);
                using signed_t = std::make_signed_t<unsigned_t>;
                // There should be no UB here.
                util::constexpr_assert(
                    static_cast<unsigned_t>(-n) >=
                    static_cast<unsigned_t>(std::numeric_limits<signed_t>::min()));
                return static_cast<signed_t>(-n);
            }
            struct to_negative_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return to_negative(static_cast<decltype(n)&&>(n));
                }
            };

            constexpr auto abs(std::unsigned_integral auto n) noexcept { return n; }
            constexpr auto abs(std::signed_integral auto n) noexcept {
                using unsigned_t = std::make_unsigned_t<decltype(n)>;
                return n >= 0 ? static_cast<unsigned_t>(n)
                              : static_cast<unsigned_t>(-static_cast<unsigned_t>(n));
            }
            struct abs_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return abs(static_cast<decltype(n)&&>(n));
                }
            };

            // Some inspection functions.
            constexpr bool is_zero(std::unsigned_integral auto n) noexcept { return n == 0; }
            constexpr bool is_zero(std::signed_integral auto n) noexcept { return n == 0; }
            struct is_zero_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return is_zero(static_cast<decltype(n)&&>(n));
                }
            };

            constexpr bool is_even(std::unsigned_integral auto n) noexcept { return n % 2 == 0; }
            constexpr bool is_even(std::signed_integral auto n) noexcept { return n % 2 == 0; }
            struct is_even_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return is_even(static_cast<decltype(n)&&>(n));
                }
            };

            constexpr bool is_strictly_positive(std::unsigned_integral auto n) noexcept {
                return n > 0;
            }
            constexpr bool is_strictly_positive(std::signed_integral auto n) noexcept {
                return n > 0;
            }
            struct is_strictly_positive_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return is_strictly_positive(static_cast<decltype(n)&&>(n));
                }
            };

            constexpr bool is_strictly_negative(std::unsigned_integral auto) noexcept {
                return false;
            }
            constexpr bool is_strictly_negative(std::signed_integral auto n) noexcept {
                return n < 0;
            }
            struct is_strictly_negative_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return is_strictly_negative(static_cast<decltype(n)&&>(n));
                }
            };

            constexpr bool is_nonnegative(std::unsigned_integral auto) noexcept { return true; }
            constexpr bool is_nonnegative(std::signed_integral auto n) noexcept { return n >= 0; }
            struct is_nonnegative_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return is_nonnegative(static_cast<decltype(n)&&>(n));
                }
            };

            constexpr bool is_nonpositive(std::unsigned_integral auto n) noexcept { return n == 0; }
            constexpr bool is_nonpositive(std::signed_integral auto n) noexcept { return n <= 0; }
            struct is_nonpositive_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return is_nonpositive(static_cast<decltype(n)&&>(n));
                }
            };

            constexpr sign_t sign(std::unsigned_integral auto) noexcept { return sign_t::positive; }
            constexpr sign_t sign(std::signed_integral auto n) noexcept {
                return n >= 0 ? sign_t::positive : sign_t::negative;
            }
            struct sign_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return sign(static_cast<decltype(n)&&>(n));
                }
            };

            // Sign inversion functions.
            constexpr auto invert_sign(std::signed_integral auto n) noexcept {
                using signed_t = decltype(n);
                using unsigned_t = std::make_unsigned_t<signed_t>;
                auto const unsigned_negative_n =
                    static_cast<unsigned_t>(-static_cast<unsigned_t>(n));

                // There should be no UB here.
                util::constexpr_assert(
                    (n >= 0 && unsigned_negative_n >=
                                   static_cast<unsigned_t>(std::numeric_limits<signed_t>::min())) ||
                    (n < 0 && unsigned_negative_n <=
                                  static_cast<unsigned_t>(std::numeric_limits<signed_t>::max())));
                return static_cast<signed_t>(unsigned_negative_n);
            }
            struct invert_sign_impl {
                constexpr decltype(auto) operator()(auto&& n) const {
                    return invert_sign(static_cast<decltype(n)&&>(n));
                }
            };

            // Divisions.
            template <std::unsigned_integral UInt>
            constexpr div_t<UInt> div(UInt x, UInt y) noexcept {
                return {static_cast<UInt>(x / y), static_cast<UInt>(x % y)};
            }
            // Follows the convention 0 <= r < q, where q is the divisor and r is the remainder.
            template <std::signed_integral Int, std::unsigned_integral UInt>
            constexpr auto div(Int x, UInt y) noexcept {
                using common_uint = std::common_type_t<std::make_unsigned_t<Int>, UInt>;

                auto const abs_x = static_cast<common_uint>(abs(x));
                auto div_result = div(abs_x, y);
                if (x < 0 && div_result.rem != 0) {
                    ++div_result.quot;
                    div_result.rem = static_cast<common_uint>(y - div_result.rem);
                }
                return div_t<std::make_signed_t<common_uint>, common_uint>{
                    to_negative(div_result.quot), div_result.rem};
            }
            struct div_impl {
                constexpr decltype(auto) operator()(auto&& x, auto&& y) const {
                    return div(static_cast<decltype(x)&&>(x), static_cast<decltype(y)&&>(y));
                }
            };

            template <std::unsigned_integral UInt>
            constexpr UInt div_floor(UInt x, UInt y) noexcept {
                return static_cast<UInt>(x / y);
            }
            template <std::signed_integral Int, std::unsigned_integral UInt>
            constexpr auto div_floor(Int x, UInt y) noexcept {
                return div(x, y).quot;
            }
            template <std::unsigned_integral UInt, std::signed_integral Int>
            constexpr auto div_floor(UInt x, Int y) noexcept {
                using common_uint = std::common_type_t<UInt, std::make_unsigned_t<Int>>;

                if (y > 0) {
                    return static_cast<common_uint>(x / static_cast<common_uint>(y));
                }
                else {
                    auto const abs_y = static_cast<common_uint>(abs(y));
                    auto div_result = div(x, abs_y);
                    if (div_result.rem != 0) {
                        ++div_result.quot;
                    }
                    return to_negative(div_result.quot);
                }
            }
            template <std::signed_integral Int>
            constexpr Int div_floor(Int x, Int y) noexcept {
                auto const abs_x = abs(x);
                auto const abs_y = abs(y);
                if (sign(x) == sign(y)) {
                    return abs_x / abs_y;
                }
                else {
                    auto div_result = div(abs_x, abs_y);
                    if (div_result.rem != 0) {
                        ++div_result.quot;
                    }
                    return to_negative(div_result.quot);
                }
            }
            struct div_floor_impl {
                constexpr decltype(auto) operator()(auto&& x, auto&& y) const {
                    return div_floor(static_cast<decltype(x)&&>(x), static_cast<decltype(y)&&>(y));
                }
            };

            template <std::unsigned_integral UInt>
            constexpr UInt div_ceil(UInt x, UInt y) noexcept {
                auto div_result = div(x, y);
                if (div_result.rem != 0) {
                    ++div_result.quot;
                }
                return div_result.quot;
            }
            template <std::signed_integral Int, std::unsigned_integral UInt>
            constexpr auto div_ceil(Int x, UInt y) noexcept {
                using common_uint = std::common_type_t<std::make_unsigned_t<Int>, UInt>;
                auto const abs_x = static_cast<common_uint>(abs(x));
                if (x > 0) {
                    return to_signed(div_ceil(abs_x, y));
                }
                else {
                    return to_negative(abs_x / y);
                }
            }
            template <std::unsigned_integral UInt, std::signed_integral Int>
            constexpr auto div_ceil(UInt x, Int y) noexcept {
                using common_uint = std::common_type_t<UInt, std::make_unsigned_t<Int>>;
                auto const abs_y = static_cast<common_uint>(abs(y));
                if (y > 0) {
                    return to_signed(div_ceil(x, abs_y));
                }
                else {
                    return to_negative(x / abs_y);
                }
            }
            template <std::signed_integral Int>
            constexpr Int div_ceil(Int x, Int y) noexcept {
                auto const abs_x = abs(x);
                auto const abs_y = abs(y);
                if (sign(x) == sign(y)) {
                    return to_signed(div_ceil(abs_x, abs_y));
                }
                else {
                    return to_negative(abs_x / abs_y);
                }
            }
            struct div_ceil_impl {
                constexpr decltype(auto) operator()(auto&& x, auto&& y) const {
                    return div_ceil(static_cast<decltype(x)&&>(x), static_cast<decltype(y)&&>(y));
                }
            };
        }
        inline constexpr auto to_signed = detail::to_signed_impl{};
        inline constexpr auto to_negative = detail::to_negative_impl{};
        inline constexpr auto abs = detail::abs_impl{};
        inline constexpr auto is_zero = detail::is_zero_impl{};
        inline constexpr auto is_even = detail::is_even_impl{};
        inline constexpr auto is_strictly_positive = detail::is_strictly_positive_impl{};
        inline constexpr auto is_strictly_negative = detail::is_strictly_negative_impl{};
        inline constexpr auto is_nonnegative = detail::is_nonnegative_impl{};
        inline constexpr auto is_nonpositive = detail::is_nonpositive_impl{};
        inline constexpr auto sign = detail::sign_impl{};
        inline constexpr auto invert_sign = detail::invert_sign_impl{};
        inline constexpr auto div = detail::div_impl{};
        inline constexpr auto div_floor = detail::div_floor_impl{};
        inline constexpr auto div_ceil = detail::div_ceil_impl{};
    }

    // Some metaprogramming utilities.
    namespace tmp {
        template <class... Types>
        struct typelist {
            static constexpr std::size_t size = sizeof...(Types);
        };

        template <class Type, class... Types>
        constexpr std::size_t find_first_index(typelist<Types...>) noexcept {
            bool found = false;
            auto impl = [&](auto arg) -> std::size_t {
                if (!found) {
                    if (std::is_same_v<Type, typename decltype(arg)::type>) {
                        found = true;
                        return 0;
                    }
                    return 1;
                }
                return 0;
            };
            return (std::size_t(0) + ... + impl(std::type_identity<Types>{}));
        }

        namespace detail {
            template <auto...>
            struct value_placeholder {
                template <class T>
                constexpr value_placeholder(T&&) noexcept {}
            };

            template <std::size_t... I, class... Types>
            constexpr auto get_type_helper(std::index_sequence<I...>, typelist<Types...>) noexcept {
                return [](value_placeholder<I>..., auto nth, auto...) {
                    return nth;
                }(std::type_identity<Types>{}...);
            }

            template <class T, std::size_t index>
            struct indexed_type_placeholder {
                using type = T;

                template <class U>
                constexpr indexed_type_placeholder(U&&) noexcept {}
            };

            template <std::size_t... I, class... Types>
            constexpr auto back_sublist_helper(std::index_sequence<I...>,
                                               typelist<Types...>) noexcept {
                return [](value_placeholder<I>..., auto... args) {
                    return typelist<typename decltype(args)::type...>{};
                }(std::type_identity<Types>{}...);
            }
        }

        template <std::size_t N, class Typelist>
        using get_type = typename decltype(detail::get_type_helper(std::make_index_sequence<N>{},
                                                                   Typelist{}))::type;

        template <std::size_t N, class Typelist>
        using back_sublist = decltype(detail::back_sublist_helper(
            std::make_index_sequence<Typelist::size - N>{}, Typelist{}));

        namespace detail {
            template <auto prefix_sum, class... Types>
            constexpr auto prefix_sum_compaction(typelist<Types...>) noexcept {
                if constexpr (sizeof...(Types) == 0) {
                    return typelist<>{};
                }
                else {
                    using list = typelist<Types...>;
                    constexpr std::size_t size = prefix_sum[sizeof...(Types) - 1];
                    constexpr auto index_array = [] {
                        util::array<std::size_t, size> result{};
                        for (std::size_t i = sizeof...(Types); i > 0; --i) {
                            result[prefix_sum[i - 1] - 1] = i - 1;
                        }
                        return result;
                    }();

                    return [&index_array]<std::size_t... I>(std::index_sequence<I...>) {
                        return typelist<get_type<index_array[I], list>...>{};
                    }(std::make_index_sequence<size>{});
                }
            }

            template <class... Types>
            constexpr auto remove_duplicate_impl(typelist<Types...>) noexcept {
                using list = typelist<Types...>;
                constexpr auto prefix_sum = [] {
                    std::size_t count = 0;
                    std::size_t index = 0;
                    auto impl = [&](auto arg) {
                        if (find_first_index<typename decltype(arg)::type>(list{}) == index++) {
                            return ++count;
                        }
                        else {
                            return count;
                        }
                    };
                    return util::array<std::size_t, sizeof...(Types)>{
                        impl(std::type_identity<Types>{})...};
                }();

                return prefix_sum_compaction<prefix_sum>(list{});
            }
        }

        // Guranteed to preserve the order.
        template <class Typelist>
        using remove_duplicate = decltype(detail::remove_duplicate_impl(Typelist{}));

        template <class... Types>
        constexpr bool has_duplicate(typelist<Types...>) noexcept {
            using list = typelist<Types...>;
            std::size_t index = 0;
            auto impl = [&](auto arg) {
                if (find_first_index<typename decltype(arg)::type>(list{}) == index++) {
                    return false;
                }
                else {
                    return true;
                }
            };
            return (impl(std::type_identity<Types>{}) || ...);
        }

        namespace detail {
            template <class Type, class... Types>
            constexpr typelist<Types..., Type> push_back_impl(typelist<Types...>) noexcept {
                return {};
            }

            constexpr typelist<> join_impl() noexcept { return {}; }
            template <class... Types>
            constexpr typelist<Types...> join_impl(typelist<Types...>) noexcept {
                return {};
            }
            template <class... Types1, class... Types2>
            constexpr typelist<Types1..., Types2...> join_impl(typelist<Types1...>,
                                                               typelist<Types2...>) noexcept {
                return {};
            }
            template <class... Types1, class... Types2, class... Typelists>
            constexpr auto join_impl(typelist<Types1...> first, typelist<Types2...> second,
                                     Typelists... remaining) noexcept {
                return join_impl(join_impl(first, second), join_impl(remaining...));
            }
        }

        template <class Typelist, class Type>
        using push_back = decltype(detail::push_back_impl<Type>(Typelist{}));

        template <class... Typelists>
        using join = decltype(detail::join_impl(Typelists{}...));

        namespace detail {
            template <class Predicate, class... Types>
            constexpr auto filter_impl(typelist<Types...>) noexcept {
                using list = typelist<Types...>;
                constexpr auto prefix_sum = [] {
                    std::size_t count = 0;
                    auto impl = [&](auto arg) {
                        if (Predicate{}(arg)) {
                            return ++count;
                        }
                        else {
                            return count;
                        }
                    };
                    return util::array<std::size_t, sizeof...(Types)>{
                        impl(std::type_identity<Types>{})...};
                }();

                return prefix_sum_compaction<prefix_sum>(list{});
            }
        }

        // Guranteed to preserve the order.
        // The predicate is evaluated in the form Predicate{}(std::type_identity<T>{}) on type T.
        template <class Typelist, class Predicate>
        using filter = decltype(detail::filter_impl<Predicate>(Typelist{}));
    }
}

#endif
