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

#ifndef JKJ_HEADER_INTERVAL
#define JKJ_HEADER_INTERVAL

#include "util.h"
#include <concepts>

namespace jkj {
    enum class boundary_type_t : bool { open, closed };

    enum class interval_type_t {
        empty,                          // empty
        bounded_open,                   // (a,b)
        bounded_left_open_right_closed, // (a,b]
        bounded_left_closed_right_open, // [a,b)
        bounded_closed,                 // [a,b]
        bounded_below_open,             // (a,infty)
        bounded_below_closed,           // [a,infty)
        bounded_above_open,             // (-infty,b)
        bounded_above_closed,           // (-infty,b]
        entire                          // (-infty,infty)
    };

    template <std::totally_ordered Value, interval_type_t it>
    struct interval;

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::empty> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept { return interval_type_t::empty; }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        template <class T>
        static constexpr bool contains(T const&) noexcept {
            return false;
        }
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_open> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept { return interval_type_t::bounded_open; }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T, class U>
            requires(std::is_constructible_v<Value, T> && std::is_constructible_v<Value, U>)
        explicit constexpr interval(T&& lower_bound, U&& upper_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)},
              upper_bound_{static_cast<U&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return lower_bound_ < x && x < upper_bound_;
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_left_open_right_closed> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_left_open_right_closed;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T, class U>
            requires(std::is_constructible_v<Value, T> && std::is_constructible_v<Value, U>)
        explicit constexpr interval(T&& lower_bound, U&& upper_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)},
              upper_bound_{static_cast<U&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return lower_bound_ < x && x <= upper_bound_;
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_left_closed_right_open> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_left_closed_right_open;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::closed; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T, class U>
            requires(std::is_constructible_v<Value, T> && std::is_constructible_v<Value, U>)
        explicit constexpr interval(T&& lower_bound, U&& upper_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)},
              upper_bound_{static_cast<U&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return lower_bound_ <= x && x < upper_bound_;
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_closed> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept { return interval_type_t::bounded_closed; }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::closed; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T, class U>
            requires(std::is_constructible_v<Value, T> && std::is_constructible_v<Value, U>)
        explicit constexpr interval(T&& lower_bound, U&& upper_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)},
              upper_bound_{static_cast<U&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return lower_bound_ <= x && x <= upper_bound_;
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_below_open> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_below_open;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        template <class T>
            requires(!std::is_base_of_v<interval, T> && std::is_constructible_v<Value, T>)
        explicit constexpr interval(T&& lower_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return lower_bound_ < x;
        }

    private:
        Value lower_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_below_closed> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_below_closed;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::closed; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        template <class T>
            requires(!std::is_base_of_v<interval, T> && std::is_constructible_v<Value, T>)
        explicit constexpr interval(T&& lower_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return lower_bound_ <= x;
        }

    private:
        Value lower_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_above_open> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_above_open;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T>
            requires(!std::is_base_of_v<interval, T> && std::is_constructible_v<Value, T>)
        explicit constexpr interval(T&& upper_bound)
            : upper_bound_{static_cast<T&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return x < upper_bound_;
        }

    private:
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_above_closed> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_above_closed;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T>
            requires(!std::is_base_of_v<interval, T> && std::is_constructible_v<Value, T>)
        explicit constexpr interval(T&& upper_bound)
            : upper_bound_{static_cast<T&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return x <= upper_bound_;
        }

    private:
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::entire> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept { return interval_type_t::entire; }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        template <class T>
        static constexpr bool contains(T const&) {
            return true;
        }
    };

    enum class cyclic_interval_type_t {
        empty,                  // empty
        single_point,           // {a}
        open,                   // (a,b)
        left_open_right_closed, // (a,b]
        left_closed_right_open, // [a,b)
        closed,                 // [a,b]
        entire                  // entire
    };

    template <class T>
    concept cyclically_ordered = std::equality_comparable<T> && requires(T a, T b, T c) {
        { cyclic_order(a, b, c) } -> std::same_as<bool>;
    };

    template <cyclically_ordered Value, cyclic_interval_type_t it>
    struct cyclic_interval;

    template <cyclically_ordered Value>
    struct cyclic_interval<Value, cyclic_interval_type_t::empty> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <cyclic_interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept { return cyclic_interval_type_t::empty; }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        template <class T>
        static constexpr bool contains(T const&) noexcept {
            return false;
        }
    };

    template <cyclically_ordered Value>
    struct cyclic_interval<Value, cyclic_interval_type_t::single_point> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <cyclic_interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return cyclic_interval_type_t::single_point;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::closed; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return point_; }
        constexpr auto& lower_bound() & noexcept { return point_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(point_); }

        constexpr auto const& upper_bound() const& noexcept { return point_; }
        constexpr auto& upper_bound() & noexcept { return point_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(point_); }

        template <class T>
            requires(!std::is_base_of_v<cyclic_interval, T> && std::is_constructible_v<Value, T>)
        explicit constexpr cyclic_interval(T&& point) : point_{static_cast<T&&>(point)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return x == point_;
        }

    private:
        Value point_;
    };

    template <cyclically_ordered Value>
    struct cyclic_interval<Value, cyclic_interval_type_t::open> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <cyclic_interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept { return cyclic_interval_type_t::open; }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T, class U>
            requires(std::is_constructible_v<Value, T> && std::is_constructible_v<Value, U>)
        explicit constexpr cyclic_interval(T&& lower_bound, U&& upper_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)},
              upper_bound_{static_cast<U&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return cyclic_order(lower_bound_, x, upper_bound_);
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <cyclically_ordered Value>
    struct cyclic_interval<Value, cyclic_interval_type_t::left_open_right_closed> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <cyclic_interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return cyclic_interval_type_t::left_open_right_closed;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::open; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T, class U>
            requires(std::is_constructible_v<Value, T> && std::is_constructible_v<Value, U>)
        explicit constexpr cyclic_interval(T&& lower_bound, U&& upper_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)},
              upper_bound_{static_cast<U&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return x == upper_bound_ || cyclic_order(lower_bound_, x, upper_bound_);
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <cyclically_ordered Value>
    struct cyclic_interval<Value, cyclic_interval_type_t::left_closed_right_open> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <cyclic_interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept {
            return cyclic_interval_type_t::left_closed_right_open;
        }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::closed; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T, class U>
            requires(std::is_constructible_v<Value, T> && std::is_constructible_v<Value, U>)
        explicit constexpr cyclic_interval(T&& lower_bound, U&& upper_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)},
              upper_bound_{static_cast<U&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return x == lower_bound_ || cyclic_order(lower_bound_, x, upper_bound_);
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <cyclically_ordered Value>
    struct cyclic_interval<Value, cyclic_interval_type_t::closed> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <cyclic_interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept { return cyclic_interval_type_t::closed; }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::closed; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        template <class T, class U>
            requires(std::is_constructible_v<Value, T> && std::is_constructible_v<Value, U>)
        explicit constexpr cyclic_interval(T&& lower_bound, U&& upper_bound)
            : lower_bound_{static_cast<T&&>(lower_bound)},
              upper_bound_{static_cast<U&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const {
            return x == lower_bound_ || x == upper_bound_ ||
                   cyclic_order(lower_bound_, x, upper_bound_);
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <cyclically_ordered Value>
    struct cyclic_interval<Value, cyclic_interval_type_t::entire> {
        using value_type = std::remove_cvref_t<Value>;

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            return static_cast<Functor&&>(f)(*this);
        }
        // Returns true if the visitation was successful.
        template <cyclic_interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if constexpr (it == interval_type()) {
                static_cast<Functor&&>(f)(*this);
                return true;
            }
            return false;
        }

        static constexpr auto interval_type() noexcept { return cyclic_interval_type_t::entire; }
        static constexpr auto left_boundary_type() noexcept { return boundary_type_t::closed; }
        static constexpr auto right_boundary_type() noexcept { return boundary_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        template <class T>
        static constexpr bool contains(T const&) noexcept {
            return true;
        }
    };

    // Variable-shape intervals implementation details.
    namespace detail {
        // Sort and remove duplications from a list of enum values.
        template <class Enum, std::size_t N>
        constexpr util::array<Enum, N>
        sort_and_remove_duplication_impl(util::array<Enum, N> arr) noexcept {
            using underlying = std::underlying_type_t<Enum>;
            std::size_t size = N;
            for (std::size_t i = 0; i < size; ++i) {
                // Find the smallest element.
                std::size_t min_idx = i;
                for (std::size_t j = i + 1; j < size; ++j) {
                    if (static_cast<underlying>(arr[min_idx]) > static_cast<underlying>(arr[j])) {
                        min_idx = j;
                    }
                }

                // Swap that smallest element with the pivot.
                using std::swap;
                swap(arr[i], arr[min_idx]);

                // Send all duplicates of that smallest element to the end of the list.
                for (std::size_t j = i + 1; j < size; ++j) {
                    if (static_cast<underlying>(arr[i]) == static_cast<underlying>(arr[j])) {
                        swap(arr[j], arr[size - 1]);
                        --size;
                    }
                }
            }

            return arr;
        }
        template <class Enum, std::size_t N>
        constexpr std::size_t count_without_duplicates(util::array<Enum, N> sorted) noexcept {
            for (std::size_t i = 1; i < N; ++i) {
                if (sorted[i - 1] >= sorted[i]) {
                    return i;
                }
            }
            return N;
        }

        template <class Enum, Enum... vals>
        constexpr auto sort_and_remove_duplication() noexcept {
            constexpr auto sorted =
                sort_and_remove_duplication_impl(util::array<Enum, sizeof...(vals)>{vals...});
            constexpr auto size = count_without_duplicates(sorted);

            util::array<Enum, size> return_value{};
            for (std::size_t i = 0; i < size; ++i) {
                return_value[i] = sorted[i];
            }
            return return_value;
        }

        template <class Value, class Enum, template <class, Enum> class StaticIntervalTemplate,
                  auto allowed_interval_types>
        class variable_shape_interval_impl;

        template <class Value, class Enum, template <class, Enum> class StaticIntervalTemplate,
                  std::size_t N, util::array<Enum, N> allowed_interval_types_arg>
        class variable_shape_interval_impl<Value, Enum, StaticIntervalTemplate,
                                           allowed_interval_types_arg> {
            static constexpr auto allowed_interval_types_ = allowed_interval_types_arg;

            static_assert(allowed_interval_types_.size() != 0,
                          "no allowed interval type is specified.");

            using storage_type = decltype([]<std::size_t... I>(std::index_sequence<I...>) {
                return util::variant_storage<
                    StaticIntervalTemplate<Value, allowed_interval_types_[I]>...>{};
            }(std::make_index_sequence<allowed_interval_types_.size()>{}));

            Enum interval_type_;
            storage_type storage_;

            // Check if all the alternatives provide lower_bound().
            static constexpr bool lower_bound_exists() noexcept {
                return []<std::size_t... I>(std::index_sequence<I...>) {
                    return ((requires(StaticIntervalTemplate<Value, allowed_interval_types_[I]> x) {
                                x.lower_bound();
                            }) && ...);
                }(std::make_index_sequence<N>{});
            }
            // Check if all the alternatives provide upper_bound().
            static constexpr bool upper_bound_exists() noexcept {
                return []<std::size_t... I>(std::index_sequence<I...>) {
                    return ((requires(StaticIntervalTemplate<Value, allowed_interval_types_[I]> x) {
                                x.upper_bound();
                            }) && ...);
                }(std::make_index_sequence<N>{});
            }

        public:
            template <class, class Enum_, template <class, Enum_> class, auto>
            friend class variable_shape_interval_impl;

            static constexpr auto const& allowed_interval_types() noexcept {
                return allowed_interval_types_;
            }
            static constexpr bool is_allowed_interval_type(Enum it) noexcept {
                for (std::size_t idx = 0; idx < allowed_interval_types_.size(); ++idx) {
                    if (it == allowed_interval_types_[idx]) {
                        return true;
                    }
                }
                return false;
            }
            template <std::size_t NOther>
            static constexpr bool are_allowed_interval_types(
                util::array<Enum, NOther> const& allowed_interval_types_arg_other) noexcept {
                for (std::size_t idx = 0; idx < NOther; ++idx) {
                    if (!is_allowed_interval_type(allowed_interval_types_arg_other[idx])) {
                        return false;
                    }
                }
                return true;
            }

            using value_type = std::remove_cvref_t<Value>;

            constexpr variable_shape_interval_impl(variable_shape_interval_impl const& itv)
                : interval_type_{itv.interval_type()} {
                itv.visit([this](auto const& itv_) {
                    storage_.template construct<std::remove_cvref_t<decltype(itv_)>>(itv_);
                });
            }

            constexpr variable_shape_interval_impl(variable_shape_interval_impl&& itv)
                : interval_type_{itv.interval_type()} {
                itv.visit([this](auto&& itv_) {
                    storage_.template construct<std::remove_cvref_t<decltype(itv_)>>(
                        std::move(itv_));
                });
            }

            template <class T, Enum it>
                requires std::is_constructible_v<Value, T>
            constexpr variable_shape_interval_impl(StaticIntervalTemplate<T, it> itv)
                : interval_type_{it} {
                static_assert(is_allowed_interval_type(it),
                              "specified interval type not is allowed");
                storage_.template construct<
                    StaticIntervalTemplate<Value, decltype(itv)::interval_type()>>(std::move(itv));
            }

            template <class T, std::size_t NOther,
                      util::array<Enum, NOther> allowed_interval_types_arg_other>
                requires std::is_constructible_v<Value, T>
            constexpr variable_shape_interval_impl(
                variable_shape_interval_impl<T, Enum, StaticIntervalTemplate,
                                             allowed_interval_types_arg_other> const& itv)
                : interval_type_{itv.interval_type()} {
                static_assert(are_allowed_interval_types(allowed_interval_types_arg_other),
                              "one of the possible interval type is not allowed");

                itv.visit([this](auto const& itv_) {
                    storage_.template construct<StaticIntervalTemplate<
                        Value, std::remove_cvref_t<decltype(itv_)>::interval_type()>>(itv_);
                });
            }

            template <class T, std::size_t NOther,
                      util::array<Enum, NOther> allowed_interval_types_arg_other>
                requires std::is_constructible_v<Value, T>
            constexpr variable_shape_interval_impl(
                variable_shape_interval_impl<T, Enum, StaticIntervalTemplate,
                                             allowed_interval_types_arg_other>&& itv)
                : interval_type_{itv.interval_type()} {
                static_assert(are_allowed_interval_types(allowed_interval_types_arg_other),
                              "one of the possible interval type is not allowed");

                itv.visit([this](auto&& itv_) {
                    storage_.template construct<StaticIntervalTemplate<
                        Value, std::remove_cvref_t<decltype(itv_)>::interval_type()>>(
                        std::move(itv_));
                });
            }

            constexpr ~variable_shape_interval_impl() {
                visit([](auto&& itv) { std::destroy_at(&itv); });
            }

            constexpr variable_shape_interval_impl&
            operator=(variable_shape_interval_impl const& itv) {
                itv.visit([this](auto const& itv_) {
                    storage_.template assign_from<allowed_interval_types_>(interval_type_, itv_);
                });
                interval_type_ = itv.interval_type();
                return *this;
            }

            constexpr variable_shape_interval_impl& operator=(variable_shape_interval_impl&& itv) {
                itv.visit([this](auto&& itv_) {
                    storage_.template assign_from<allowed_interval_types_>(interval_type_,
                                                                           std::move(itv_));
                });
                interval_type_ = itv.interval_type();
                return *this;
            }

            template <class T, Enum it>
                requires std::is_constructible_v<Value, T>
            constexpr variable_shape_interval_impl& operator=(StaticIntervalTemplate<T, it> itv) {
                static_assert(is_allowed_interval_type(it),
                              "specified interval type is not allowed");
                storage_.template assign_from<allowed_interval_types_>(interval_type_,
                                                                       std::move(itv));
                interval_type_ = itv.interval_type();
                return *this;
            }

            template <class T, std::size_t NOther,
                      util::array<Enum, NOther> allowed_interval_types_arg_other>
                requires std::is_constructible_v<Value, T>
            constexpr variable_shape_interval_impl&
            operator=(variable_shape_interval_impl<T, Enum, StaticIntervalTemplate,
                                                   allowed_interval_types_arg_other> const& itv) {
                static_assert(are_allowed_interval_types(allowed_interval_types_arg_other),
                              "one of the possible interval type is not allowed");
                itv.visit([this](auto const& itv_) {
                    storage_.template assign_from<allowed_interval_types_>(interval_type_, itv_);
                });
                interval_type_ = itv.interval_type();
                return *this;
            }

            template <class T, std::size_t NOther,
                      util::array<Enum, NOther> allowed_interval_types_arg_other>
                requires std::is_constructible_v<Value, T>
            constexpr variable_shape_interval_impl&
            operator=(variable_shape_interval_impl<T, Enum, StaticIntervalTemplate,
                                                   allowed_interval_types_arg_other>&& itv) {
                static_assert(are_allowed_interval_types(allowed_interval_types_arg_other),
                              "one of the possible interval type is not allowed");
                itv.visit([this](auto&& itv_) {
                    storage_.template assign_from<allowed_interval_types_>(interval_type_,
                                                                           std::move(itv_));
                });
                interval_type_ = itv.interval_type();
                return *this;
            }

            template <class Functor>
            constexpr decltype(auto) visit(Functor&& f) const {
                return storage_.template visit<allowed_interval_types_>(interval_type_,
                                                                        static_cast<Functor&&>(f));
            }

            // Returns true if the visitation was successful.
            template <Enum it, class Functor>
            constexpr bool visit_if(Functor&& f) const {
                if (interval_type_ == it) {
                    static_cast<Functor&&>(f)(
                        storage_.template get<StaticIntervalTemplate<Value, it>>());
                    return true;
                }
                return false;
            }

            constexpr Enum interval_type() const noexcept { return interval_type_; }
            constexpr auto left_boundary_type() const noexcept {
                return visit([](auto&& itv) { return itv.left_boundary_type(); });
            }
            constexpr auto right_boundary_type() const noexcept {
                return visit([](auto&& itv) { return itv.right_boundary_type(); });
            }
            template <class T>
            constexpr bool contains(T const& x) const {
                return visit([&x](auto&& itv) { return itv.contains(x); });
            }

            // Call success with lower_bound() if the current interval type supports lower_bound(),
            // and call fail otherwise.
            template <class SuccessFunctor, class FailFunctor>
            constexpr decltype(auto) with_lower_bound(SuccessFunctor&& success,
                                                      FailFunctor&& fail) const {
                return visit([&](auto&& itv) -> decltype(auto) {
                    using itv_type = std::remove_cvref_t<decltype(itv)>;
                    if constexpr (requires(
                                      StaticIntervalTemplate<Value, itv_type::interval_type()> x) {
                                      x.lower_bound();
                                  }) {
                        return static_cast<SuccessFunctor&&>(success)(itv.lower_bound());
                    }
                    else {
                        return static_cast<FailFunctor&&>(fail)();
                    }
                });
            }

            // Call success with upper_bound() if the current interval type supports upper_bound(),
            // and call fail otherwise.
            template <class SuccessFunctor, class FailFunctor>
            constexpr decltype(auto) with_upper_bound(SuccessFunctor&& success,
                                                      FailFunctor&& fail) const {
                return visit([&](auto&& itv) -> decltype(auto) {
                    using itv_type = std::remove_cvref_t<decltype(itv)>;
                    if constexpr (requires(
                                      StaticIntervalTemplate<Value, itv_type::interval_type()> x) {
                                      x.upper_bound();
                                  }) {
                        return static_cast<SuccessFunctor&&>(success)(itv.upper_bound());
                    }
                    else {
                        return static_cast<FailFunctor&&>(fail)();
                    }
                });
            }

            constexpr value_type const& lower_bound() const noexcept
                requires(lower_bound_exists())
            {
                return visit([&](auto&& itv) -> decltype(auto) { return itv.lower_bound(); });
            }

            constexpr value_type const& upper_bound() const noexcept
                requires(upper_bound_exists())
            {
                return visit([&](auto&& itv) -> decltype(auto) { return itv.upper_bound(); });
            }
        };
    }

    // Variable shape intervals. Possible interval types can be specified. If no type is specified,
    // any type is allowed.
    template <std::totally_ordered Value, interval_type_t... allowed_interval_types>
    using variable_shape_interval = std::conditional_t<
        sizeof...(allowed_interval_types) == 0,
        detail::variable_shape_interval_impl<
            Value, interval_type_t, interval,
            util::array<interval_type_t, 10>{
                interval_type_t::empty, interval_type_t::bounded_open,
                interval_type_t::bounded_left_open_right_closed,
                interval_type_t::bounded_left_closed_right_open, interval_type_t::bounded_closed,
                interval_type_t::bounded_below_open, interval_type_t::bounded_below_closed,
                interval_type_t::bounded_above_open, interval_type_t::bounded_above_closed,
                interval_type_t::entire}>,
        detail::variable_shape_interval_impl<
            Value, interval_type_t, interval,
            detail::sort_and_remove_duplication<interval_type_t, allowed_interval_types...>()>>;

    // Variable shape cyclic intervals. Possible cyclic interval types can be specified. If no type
    // is specified, any type is allowed.
    template <cyclically_ordered Value, cyclic_interval_type_t... allowed_interval_types>
    using variable_shape_cyclic_interval = std::conditional_t<
        sizeof...(allowed_interval_types) == 0,
        detail::variable_shape_interval_impl<
            Value, cyclic_interval_type_t, cyclic_interval,
            util::array<cyclic_interval_type_t, 7>{
                cyclic_interval_type_t::empty, cyclic_interval_type_t::single_point,
                cyclic_interval_type_t::open, cyclic_interval_type_t::left_closed_right_open,
                cyclic_interval_type_t::left_closed_right_open, cyclic_interval_type_t::closed,
                cyclic_interval_type_t::entire}>,
        detail::variable_shape_interval_impl<
            Value, cyclic_interval_type_t, cyclic_interval,
            detail::sort_and_remove_duplication<cyclic_interval_type_t,
                                                allowed_interval_types...>()>>;
}

#endif
