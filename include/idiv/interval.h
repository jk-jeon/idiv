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
    enum class endpoint_type_t { empty, open, closed, infinity };

    template <std::totally_ordered Value, interval_type_t it>
    struct interval;

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::empty> {
        static constexpr auto interval_type() noexcept { return interval_type_t::empty; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::empty; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::empty; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_open> {
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

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_left_open_right_closed> {
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

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_left_closed_right_open> {
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

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_closed> {
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

    private:
        Value lower_bound_;
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_below_open> {
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_below_open;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::open; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::infinity; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        constexpr interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)} {}

    private:
        Value lower_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_below_closed> {
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_below_closed;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::closed; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::infinity; }

        constexpr auto const& lower_bound() const& noexcept { return lower_bound_; }
        constexpr auto& lower_bound() & noexcept { return lower_bound_; }
        constexpr auto&& lower_bound() && noexcept { return static_cast<Value&&>(lower_bound_); }

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;

        constexpr interval(Value lower_bound, Value upper_bound)
            : lower_bound_{static_cast<Value&&>(lower_bound)} {}

    private:
        Value lower_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_above_open> {
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_above_open;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::infinity; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::open; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr interval(Value lower_bound, Value upper_bound)
            : upper_bound_{static_cast<Value&&>(upper_bound)} {}

    private:
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::bounded_above_closed> {
        static constexpr auto interval_type() noexcept {
            return interval_type_t::bounded_above_closed;
        }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::infinity; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::closed; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept { return upper_bound_; }
        constexpr auto& upper_bound() & noexcept { return upper_bound_; }
        constexpr auto&& upper_bound() && noexcept { return static_cast<Value&&>(upper_bound_); }

        constexpr interval(Value lower_bound, Value upper_bound)
            : upper_bound_{static_cast<Value&&>(upper_bound)} {}

    private:
        Value upper_bound_;
    };

    template <std::totally_ordered Value>
    struct interval<Value, interval_type_t::entire> {
        static constexpr auto interval_type() noexcept { return interval_type_t::entire; }
        static constexpr auto left_endpoint_type() noexcept { return endpoint_type_t::infinity; }
        static constexpr auto right_endpoint_type() noexcept { return endpoint_type_t::infinity; }

        constexpr auto const& lower_bound() const& noexcept = delete;
        constexpr auto& lower_bound() & noexcept = delete;
        constexpr auto&& lower_bound() && noexcept = delete;

        constexpr auto const& upper_bound() const& noexcept = delete;
        constexpr auto& upper_bound() & noexcept = delete;
        constexpr auto&& upper_bound() && noexcept = delete;
    };
}

#endif