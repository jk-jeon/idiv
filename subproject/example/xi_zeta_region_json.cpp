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
//
// 
// [Corresponding output]
//
// When zeta = 80765 / 36899:
//   (-317 - zeta) / 19182 <= xi <= (297 - zeta) / -17717
// 
// When zeta in (80765 / 36899, 366865 / 167608):
//   (-1545 - zeta) / 92980 <= xi <= (1525 - zeta) / -91515
// 
// When zeta = 366865 / 167608:
//   (-1545 - zeta) / 92980 <= xi < (1525 - zeta) / -91515
// 
// When zeta in (366865 / 167608, 360025 / 164483):
//   (-1545 - zeta) / 92980 <= xi < (-1264 - zeta) / 76093
// 
// When zeta = 360025 / 164483:
//   (-1545 - zeta) / 92980 < xi < (-1264 - zeta) / 76093
// 
// When zeta in (360025 / 164483, 80766 / 36899):
//   (1192 - zeta) / -71503 < xi < (-1264 - zeta) / 76093

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

template <>
struct daw::json::json_data_contract<
    jkj::interval<jkj::bigint::int_var, jkj::interval_type_t::bounded_closed>> {
    using type = json_member_list<json_number<"min", std::int_least64_t>,
                                  json_number<"max", std::int_least64_t>>;
};

struct input_params {
    continued_fraction_generator x;
    continued_fraction_generator y;
    std::vector<jkj::interval<jkj::bigint::int_var, jkj::interval_type_t::bounded_closed>> ranges;
};

template <>
struct daw::json::json_data_contract<input_params> {
    using type = json_member_list<
        json_class<"x", continued_fraction_generator>,
        json_class<"y", continued_fraction_generator>,
        json_array<"ranges", ::jkj::interval<::jkj::bigint::int_var,
                                             ::jkj::interval_type_t::bounded_closed>>>;
};

#include <fstream>
#include <sstream>

namespace jkj {
    namespace bigint {
        static std::ostream& operator<<(std::ostream& out, jkj::bigint::uint_var const& n) {
            auto digits = n.to_decimal();
            if (digits.empty()) {
                out << "0";
            }
            else {
                auto itr = digits.cbegin();
                out << std::format("{}", *itr);

                for (++itr; itr != digits.cend(); ++itr) {
                    out << std::format("{:019d}", *itr);
                }
            }

            return out;
        }

        static std::ostream& operator<<(std::ostream& out, jkj::bigint::int_var const& n) {
            if (n.is_strictly_negative()) {
                out << "-";
            }
            out << n.abs();
            return out;
        }
    }
}

int main() {
    auto const params = [] {
        std::ifstream input_file{"xi_zeta_region_json_example_input.json"};
        std::stringstream strstream;
        strstream << input_file.rdbuf();
        return daw::json::from_json<input_params>(strstream.str());
    }();

    auto const result = jkj::idiv::find_xi_zeta_region(params.x, params.y, params.ranges);

    for (auto const& elementary_region : result) {
        if (jkj::util::is_zero(elementary_region.xi_left_endpoint_denominator)) {
            elementary_region.zeta_range.visit([](auto const& itv) {
                using itv_type = std::remove_cvref_t<decltype(itv)>;

                if constexpr (itv_type::interval_type() == jkj::interval_type_t::bounded_open) {
                    std::cout << itv.lower_bound().numerator << " / "
                              << itv.lower_bound().denominator << " <= zeta <"
                              << itv.upper_bound().numerator << " / "
                              << itv.upper_bound().denominator << ")\n\n";
                }
            });
            break;
        }

        std::cout << "When zeta ";
        elementary_region.zeta_range.visit([](auto const& itv) {
            using itv_type = std::remove_cvref_t<decltype(itv)>;

            if constexpr (itv_type::interval_type() == jkj::interval_type_t::bounded_open) {
                std::cout << "in (" << itv.lower_bound().numerator << " / "
                          << itv.lower_bound().denominator << ", " << itv.upper_bound().numerator
                          << " / " << itv.upper_bound().denominator << "):\n";
            }
            else if constexpr (itv_type::interval_type() == jkj::interval_type_t::bounded_closed) {
                std::cout << "= " << itv.lower_bound().numerator << " / "
                          << itv.lower_bound().denominator << ":\n";
            }
            else {
                std::cout << "in (-inf, inf):\n";
            }
        });

        std::cout << "  (" << elementary_region.xi_left_endpoint_numerator << " - zeta) / "
                  << elementary_region.xi_left_endpoint_denominator;

        if (elementary_region.xi_left_endpoint_included) {
            std::cout << " <= xi ";
        }
        else {
            std::cout << " < xi ";
        }
        if (elementary_region.xi_right_endpoint_included) {
            std::cout << "<= ";
        }
        else {
            std::cout << "< ";
        }

        std::cout << "(" << elementary_region.xi_right_endpoint_numerator << " - zeta) / "
                  << elementary_region.xi_right_endpoint_denominator << "\n\n";
    }
}