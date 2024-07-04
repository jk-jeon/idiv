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

#include <idiv/idiv.h>
#include <idiv/log_continued_fraction.h>
#include <idiv/type_erased_continued_fraction.h>
#include <daw/json/daw_json_link.h>
#include <stdexcept>

/////////////////////////////////////////////////////////////////////////////////
// [Example input JSON file]
//
// {
//   "x" : {
//     "type" : "rational",
//     "params" : {
//       "numerator" : -614,
//       "denominator" : 36899
//     }
//   },
//   "y" : {
//     "type" : "unary_gosper",
//     "params" : {
//       "x" : {
//         "type" : "binary_gosper",
//         "params" : {
//           "x" : {
//             "type" : "general_log",
//             "params" : {
//               "base" : {
//                 "numerator" : 10,
//                 "denominator" : 3
//               },
//               "number" : {
//                 "numerator" : 5,
//                 "denominator" : 2
//               }
//             }
//           },
//           "y" : {
//             "type" : "natural_log",
//             "params" : {
//               "number" : {
//                 "numerator" : 18,
//                 "denominator" : 7
//               }
//             }
//           },
//           "coefficients" : [
//             3, 6, 17, -6,
//             -5, 13, 9, 1
//           ]
//         }
//       },
//       "coefficients" : [
//         -7, 14,
//         0, 3
//       ]
//     }
//   },
//   "ranges" : [
//     { "min" : -100, "max" : 100 },
//     { "min" : -100000, "max" : -700 },
//     { "min" : 300, "max" : 100000 }
//   ]
// }

using partial_fraction_type =
    jkj::cntfrc::projective_rational<jkj::cntfrc::unity, jkj::bigint::int_var>;
using convergent_type =
    jkj::cntfrc::projective_rational<jkj::bigint::int_var, jkj::bigint::uint_var>;
using interval_type =
    jkj::variable_shape_cyclic_interval<convergent_type, jkj::cyclic_interval_type_t::single_point,
                                        jkj::cyclic_interval_type_t::left_open_right_closed,
                                        jkj::cyclic_interval_type_t::left_closed_right_open,
                                        jkj::cyclic_interval_type_t::entire>;

using continued_fraction_generator = jkj::cntfrc::generator<
    jkj::cntfrc::impl::type_erased<partial_fraction_type, convergent_type, interval_type>,
    jkj::cntfrc::previous_previous_convergent_tracker, jkj::cntfrc::interval_tracker>;

template <>
struct daw::json::json_data_contract<
    jkj::cntfrc::impl::rational<jkj::bigint::int_var, jkj::bigint::uint_var>> {
    using type = json_member_list<json_number<"numerator", std::int_least64_t>,
                                  json_number<"denominator", std::uint_least64_t>>;
};

template <>
struct daw::json::json_data_contract<jkj::frac<jkj::bigint::uint_var, jkj::bigint::uint_var>> {
    using type = json_member_list<json_number<"numerator", std::uint_least64_t>,
                                  json_number<"denominator", std::uint_least64_t>>;
};

template <>
struct daw::json::json_data_contract<
    jkj::cntfrc::impl::natural_log<jkj::bigint::int_var, jkj::bigint::uint_var>> {
    using type = json_member_list<
        json_class<"number", ::jkj::frac<::jkj::bigint::uint_var, ::jkj::bigint::uint_var>>>;
};

template <>
struct daw::json::json_data_contract<
    jkj::cntfrc::impl::general_log<jkj::bigint::int_var, jkj::bigint::uint_var>> {
    using type = json_member_list<
        json_class<"base", ::jkj::frac<::jkj::bigint::uint_var, ::jkj::bigint::uint_var>>,
        json_class<"number", ::jkj::frac<::jkj::bigint::uint_var, ::jkj::bigint::uint_var>>>;
};

template <>
struct daw::json::json_data_contract<
    jkj::cntfrc::linear_fractional_transform<jkj::bigint::int_var>> {
    using type = json_tuple_member_list<std::int_least64_t, std::int_least64_t, std::int_least64_t,
                                        std::int_least64_t>;
};

template <>
struct daw::json::json_data_contract<
    jkj::cntfrc::impl::unary_gosper<continued_fraction_generator>> {
    using type = json_member_list<
        json_class<"x", continued_fraction_generator>,
        json_class<"coefficients",
                   ::jkj::cntfrc::linear_fractional_transform<::jkj::bigint::int_var>>>;
};

template <>
struct daw::json::json_data_contract<
    jkj::cntfrc::bilinear_fractional_transform<jkj::bigint::int_var>> {
    using type = json_tuple_member_list<std::int_least64_t, std::int_least64_t, std::int_least64_t,
                                        std::int_least64_t, std::int_least64_t, std::int_least64_t,
                                        std::int_least64_t, std::int_least64_t>;
};

template <>
struct daw::json::json_data_contract<
    jkj::cntfrc::impl::binary_gosper<continued_fraction_generator, continued_fraction_generator>> {
    using type = json_member_list<
        json_class<"x", continued_fraction_generator>,
        json_class<"y", continued_fraction_generator>,
        json_class<"coefficients",
                   ::jkj::cntfrc::bilinear_fractional_transform<::jkj::bigint::int_var>>>;
};

template <>
struct daw::json::json_data_contract<continued_fraction_generator> {
    using type =
        json_member_list<json_string<"type">, daw::json::json_raw<"params", std::string_view>>;

    struct constructor_t {
        continued_fraction_generator operator()(std::string_view type,
                                                std::string_view params_raw_json) const {
            using impl_type = ::jkj::cntfrc::impl::type_erased<partial_fraction_type,
                                                               convergent_type, interval_type>;

            if (type == "rational") {
                return continued_fraction_generator{from_json<
                    ::jkj::cntfrc::impl::rational<::jkj::bigint::int_var, ::jkj::bigint::uint_var>>(
                    params_raw_json)};
            }
            else if (type == "natural_log") {
                return continued_fraction_generator{from_json<::jkj::cntfrc::impl::natural_log<
                    ::jkj::bigint::int_var, ::jkj::bigint::uint_var>>(params_raw_json)};
            }
            else if (type == "general_log") {
                return continued_fraction_generator{from_json<::jkj::cntfrc::impl::general_log<
                    ::jkj::bigint::int_var, ::jkj::bigint::uint_var>>(params_raw_json)};
            }
            else if (type == "unary_gosper") {
                return continued_fraction_generator{
                    from_json<::jkj::cntfrc::impl::unary_gosper<continued_fraction_generator>>(
                        params_raw_json)};
            }
            else if (type == "binary_gosper") {
                return continued_fraction_generator{from_json<::jkj::cntfrc::impl::binary_gosper<
                    continued_fraction_generator, continued_fraction_generator>>(params_raw_json)};
            }
            throw std::invalid_argument{"unknown type of continued fraction"};
        }
    };
};

int main() {}