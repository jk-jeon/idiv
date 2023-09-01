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

#ifndef JKJ_HEADER_INTERVAL
#define JKJ_HEADER_INTERVAL

#include "util.h"
#include <concepts>

namespace jkj {
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
    enum class endpoint_type_t { empty, open, closed };

    template <std::totally_ordered Value, interval_type_t it>
    struct interval;

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::empty> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept { return interval_type_t::empty; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::empty; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::empty; }

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
        static constexpr auto interval_type() noexcept { return interval_type_t::bounded_open; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)},
              upper_bound_{static_cast<Value&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return lower_bound_ < x && x < upper_bound_;
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_left_open_right_closed> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_left_open_right_closed;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)},
              upper_bound_{static_cast<Value&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return lower_bound_ < x && x <= upper_bound_;
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_left_closed_right_open> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_left_closed_right_open;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::closed; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)},
              upper_bound_{static_cast<Value&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return lower_bound_ <= x && x < upper_bound_;
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_closed> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept { return interval_type_t::bounded_closed; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::closed; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)},
              upper_bound_{static_cast<Value&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return lower_bound_ <= x && x <= upper_bound_;
        }

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_below_open> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_below_open;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        constexpr interval(Value lower_bound) : lower_bound_{static_cast<Value&&>(lower_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return lower_bound_ < x;
        }

    private:
        Value lower_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_below_closed> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_below_closed;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::closed; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        constexpr interval(Value lower_bound) : lower_bound_{static_cast<Value&&>(lower_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return lower_bound_ <= x;
        }

    private:
        Value lower_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_above_open> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_above_open;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr interval(Value upper_bound) : upper_bound_{static_cast<Value&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return x < upper_bound_;
        }

    private:
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_above_closed> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_above_closed;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr interval(Value upper_bound) : upper_bound_{static_cast<Value&&>(upper_bound)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return x <= upper_bound_;
        }

    private:
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::entire> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept { return interval_type_t::entire; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

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

    template <std::totally_ordered Value>
        requires(std::is_nothrow_move_assignable_v<Value>)
    class variable_shape_interval {
        interval_type_t interval_type_;
        Value lower_bound_{};
        Value upper_bound_{};

    public:
        using value_type = std::remove_cvref_t<Value>;
        template <interval_type_t it>
        constexpr variable_shape_interval(interval<Value, it> itv) noexcept : interval_type_{it} {
            if constexpr (requires { itv.lower_bound(); }) {
                lower_bound_ = static_cast<interval<Value, it>&&>(itv).lower_bound();
            }
            if constexpr (requires { itv.upper_bound(); }) {
                upper_bound_ = static_cast<interval<Value, it>&&>(itv).upper_bound();
            }
        }

        template <interval_type_t it>
        constexpr variable_shape_interval& operator=(interval<Value, it> itv) noexcept {
            interval_type_ = it;
            if constexpr (requires { itv.lower_bound(); }) {
                lower_bound_ = static_cast<interval<Value, it>&&>(itv).lower_bound();
            }
            if constexpr (requires { itv.upper_bound(); }) {
                upper_bound_ = static_cast<interval<Value, it>&&>(itv).upper_bound();
            }
            return *this;
        }

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            using enum interval_type_t;
            switch (interval_type_) {
            case empty:
                return f(interval<Value const&, empty>{});
            case bounded_open:
                return f(interval<Value const&, bounded_open>{lower_bound_, upper_bound_});
            case bounded_left_open_right_closed:
                return f(interval<Value const&, bounded_left_open_right_closed>{lower_bound_,
                                                                                upper_bound_});
            case bounded_left_closed_right_open:
                return f(interval<Value const&, bounded_left_closed_right_open>{lower_bound_,
                                                                                upper_bound_});
            case bounded_closed:
                return f(interval<Value const&, bounded_closed>{lower_bound_, upper_bound_});
            case bounded_below_open:
                return f(interval<Value const&, bounded_below_open>{lower_bound_});
            case bounded_below_closed:
                return f(interval<Value const&, bounded_below_closed>{lower_bound_});
            case bounded_above_open:
                return f(interval<Value const&, bounded_above_open>{upper_bound_});
            case bounded_above_closed:
                return f(interval<Value const&, bounded_above_closed>{upper_bound_});
            case entire:;
            }
            return f(interval<Value const&, entire>{});
        }

        // Returns true if the visitation was successful.
        template <interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if (interval_type_ == it) {
                using enum interval_type_t;
                if constexpr (it == empty || it == entire) {
                    f(interval<Value const&, it>{});
                }
                else if constexpr (it == bounded_below_open || it == bounded_below_closed) {
                    f(interval<Value const&, it>{lower_bound_});
                }
                else if constexpr (it == bounded_above_open || it == bounded_above_closed) {
                    f(interval<Value const&, it>{upper_bound_});
                }
                else {
                    f(interval<Value const&, it>{lower_bound_, upper_bound_});
                }
                return true;
            }
            return false;
        }

        constexpr interval_type_t interval_type() const noexcept { return interval_type_; }
        constexpr auto left_endpoint_type() noexcept {
            return visit([](auto&& itv) { return itv.left_endpoint_type(); });
        }
        constexpr auto right_endpoint_type() noexcept {
            return visit([](auto&& itv) { return itv.right_endpoint_type(); });
        }
        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return visit([&x](auto&& itv) { return itv.contains(x); });
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
        static constexpr auto interval_type() noexcept { return cyclic_interval_type_t::empty; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::empty; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::empty; }

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
        static constexpr auto interval_type() noexcept {
            return cyclic_interval_type_t::single_point;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::closed; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return point_; }
        constexpr auto& lower_bound() & noexcept { return point_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(point_); }

        constexpr auto const& upper_bound() const& noexcept { return point_; }
        constexpr auto& upper_bound() & noexcept { return point_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(point_); }

        constexpr cyclic_interval(Value point) : point_{static_cast<Value&&>(point)} {}

        template <class T>
        constexpr bool contains(T const& x) const noexcept {
            return x == point_;
        }

    private:
        Value point_;
    };

    template <cyclically_ordered Value>
    struct cyclic_interval<Value, cyclic_interval_type_t::open> {
        using value_type = std::remove_cvref_t<Value>;
        static constexpr auto interval_type() noexcept { return cyclic_interval_type_t::open; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr cyclic_interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)},
              upper_bound_{static_cast<Value&&>(upper_bound)} {}

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
        static constexpr auto interval_type() noexcept {
            return cyclic_interval_type_t::left_open_right_closed;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr cyclic_interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)},
              upper_bound_{static_cast<Value&&>(upper_bound)} {}

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
        static constexpr auto interval_type() noexcept {
            return cyclic_interval_type_t::left_closed_right_open;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::closed; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr cyclic_interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)},
              upper_bound_{static_cast<Value&&>(upper_bound)} {}

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
        static constexpr auto interval_type() noexcept { return cyclic_interval_type_t::closed; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::closed; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr cyclic_interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)},
              upper_bound_{static_cast<Value&&>(upper_bound)} {}

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
        static constexpr auto interval_type() noexcept { return cyclic_interval_type_t::entire; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::empty; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::empty; }

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

    template <cyclically_ordered Value>
        requires(std::is_nothrow_move_assignable_v<Value>)
    class variable_shape_cyclic_interval {
        cyclic_interval_type_t interval_type_;
        Value lower_bound_{};
        Value upper_bound_{};

    public:
        using value_type = std::remove_cvref_t<Value>;
        template <cyclic_interval_type_t it>
        constexpr variable_shape_cyclic_interval(cyclic_interval<Value, it> itv) noexcept
            : interval_type_{it} {
            if constexpr (requires { itv.lower_bound(); }) {
                lower_bound_ = static_cast<cyclic_interval<Value, it>&&>(itv).lower_bound();
            }
            if constexpr (requires { itv.upper_bound(); }) {
                upper_bound_ = static_cast<cyclic_interval<Value, it>&&>(itv).upper_bound();
            }
        }

        template <cyclic_interval_type_t it>
        constexpr variable_shape_cyclic_interval&
        operator=(cyclic_interval<Value, it> itv) noexcept {
            interval_type_ = it;
            if constexpr (requires { itv.lower_bound(); }) {
                lower_bound_ = static_cast<cyclic_interval<Value, it>&&>(itv).lower_bound();
            }
            if constexpr (requires { itv.upper_bound(); }) {
                upper_bound_ = static_cast<cyclic_interval<Value, it>&&>(itv).upper_bound();
            }
            return *this;
        }

        template <class Functor>
        constexpr decltype(auto) visit(Functor&& f) const {
            using enum cyclic_interval_type_t;
            switch (interval_type_) {
            case empty:
                return f(cyclic_interval<Value const&, empty>{});
            case single_point:
                return f(cyclic_interval<Value const&, single_point>{lower_bound_});
            case open:
                return f(cyclic_interval<Value const&, open>{lower_bound_, upper_bound_});
            case left_open_right_closed:
                return f(cyclic_interval<Value const&, left_open_right_closed>{lower_bound_,
                                                                               upper_bound_});
            case left_closed_right_open:
                return f(cyclic_interval<Value const&, left_closed_right_open>{lower_bound_,
                                                                               upper_bound_});
            case closed:
                return f(cyclic_interval<Value const&, closed>{lower_bound_, upper_bound_});
            case entire:;
            }
            return f(cyclic_interval<Value const&, entire>{});
        }

        // Returns true if the visitation was successful.
        template <cyclic_interval_type_t it, class Functor>
        constexpr bool visit_if(Functor&& f) const {
            if (interval_type_ == it) {
                using enum cyclic_interval_type_t;
                if constexpr (it == empty || it == entire) {
                    f(cyclic_interval<Value const&, it>{});
                }
                else if constexpr (it == single_point) {
                    f(cyclic_interval<Value const&, it>{lower_bound_});
                }
                else {
                    f(cyclic_interval<Value const&, it>{lower_bound_, upper_bound_});
                }
                return true;
            }
            return false;
        }

        constexpr cyclic_interval_type_t interval_type() const noexcept { return interval_type_; }
        constexpr auto left_endpoint_type() noexcept {
            return visit([](auto&& itv) { return itv.left_endpoint_type(); });
        }
        constexpr auto right_endpoint_type() noexcept {
            return visit([](auto&& itv) { return itv.right_endpoint_type(); });
        }
        template <class T>
        constexpr bool contains(T const& x) const {
            return visit([&x](auto&& itv) { return itv.contains(x); });
        }
    };
}

#endif
