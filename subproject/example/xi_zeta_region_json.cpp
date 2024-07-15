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
//   "x": {
//     "type": "rational",
//     "params": {
//       "numerator": -614,
//       "denominator": 36899
//     }
//   },
//   "y": {
//     "type": "unary_gosper",
//     "params": {
//       "x": {
//         "type": "binary_gosper",
//         "params": {
//           "x": {
//             "type": "general_log",
//             "params": {
//               "base": {
//                 "numerator": 10,
//                 "denominator": 3
//               },
//               "number": {
//                 "numerator": 5,
//                 "denominator": 2
//               }
//             }
//           },
//           "y": {
//             "type": "natural_log",
//             "params": {
//               "number": {
//                 "numerator": 18,
//                 "denominator": 7
//               }
//             }
//           },
//           "coefficients": [
//             3,
//             6,
//             17,
//             -6,
//             -5,
//             13,
//             9,
//             1
//           ]
//         }
//       },
//       "coefficients": [
//         -7,
//         14,
//         0,
//         3
//       ]
//     }
//   },
//   "constraints": [
//     {
//       "affine_coefficients": {
//         "linear_coefficients": [
//           {
//             "numerator": 1,
//             "denominator": 1
//           },
//           {
//             "numerator": 0,
//             "denominator": 1
//           },
//           {
//             "numerator": 0,
//             "denominator": 1
//           },
//           {
//             "numerator": 1,
//             "denominator": 1
//           }
//         ],
//         "constant_coefficient_x": {
//           "numerator": 0,
//           "denominator": 1
//         },
//         "constant_coefficient_y": {
//           "numerator": 0,
//           "denominator": 1
//         }
//       },
//       "range": {
//         "min": -100,
//         "max": 100
//       }
//     },
//     {
//       "affine_coefficients": {
//         "linear_coefficients": [
//           {
//             "numerator": 1,
//             "denominator": 1
//           },
//           {
//             "numerator": 0,
//             "denominator": 1
//           },
//           {
//             "numerator": 0,
//             "denominator": 1
//           },
//           {
//             "numerator": 1,
//             "denominator": 1
//           }
//         ],
//         "constant_coefficient_x": {
//           "numerator": 0,
//           "denominator": 1
//         },
//         "constant_coefficient_y": {
//           "numerator": 0,
//           "denominator": 1
//         }
//       },
//       "range": {
//         "min": -100000,
//         "max": -700
//       }
//     },
//     {
//       "affine_coefficients": {
//         "linear_coefficients": [
//           {
//             "numerator": 1,
//             "denominator": 1
//           },
//           {
//             "numerator": 0,
//             "denominator": 1
//           },
//           {
//             "numerator": 0,
//             "denominator": 1
//           },
//           {
//             "numerator": 1,
//             "denominator": 1
//           }
//         ],
//         "constant_coefficient_x": {
//           "numerator": 0,
//           "denominator": 1
//         },
//         "constant_coefficient_y": {
//           "numerator": 0,
//           "denominator": 1
//         }
//       },
//       "range": {
//         "min": 300,
//         "max": 100000
//       }
//     }
//   ]
// }
//
//
// [Corresponding output]
//
// When -2737 / 164483 < xi < -614 / 36899,
//   (-92980 / 1) xi + (-1545 / 1) <= zeta < (71503 / 1) xi + (1192 / 1)
//
// When xi = -614 / 36899,
//   80765 / 36899 <= zeta < 80766 / 36899
//
// When -614 / 36899 < xi < -2789 / 167608,
//   (91515 / 1) xi + (1525 / 1) <= zeta < (-76093 / 1) xi + (-1264 / 1)
//

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
struct daw::json::json_data_contract<jkj::frac<jkj::bigint::int_var, jkj::bigint::uint_var>> {
    using type = json_member_list<json_number<"numerator", std::int_least64_t>,
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
    jkj::idiv::floor_constraint_spec::affine_transform::linear_transform> {
    using type =
        json_tuple_member_list<::jkj::frac<::jkj::bigint::int_var, ::jkj::bigint::uint_var>,
                               ::jkj::frac<::jkj::bigint::int_var, ::jkj::bigint::uint_var>,
                               ::jkj::frac<::jkj::bigint::int_var, ::jkj::bigint::uint_var>,
                               ::jkj::frac<::jkj::bigint::int_var, ::jkj::bigint::uint_var>>;
};

template <>
struct daw::json::json_data_contract<jkj::idiv::floor_constraint_spec::affine_transform> {
    using type = json_member_list<
        json_class<"linear_coefficients",
                   ::jkj::idiv::floor_constraint_spec::affine_transform::linear_transform>,
        json_class<"constant_coefficient_x",
                   ::jkj::frac<::jkj::bigint::int_var, ::jkj::bigint::uint_var>>,
        json_class<"constant_coefficient_y",
                   ::jkj::frac<::jkj::bigint::int_var, ::jkj::bigint::uint_var>>>;
};

template <>
struct daw::json::json_data_contract<
    jkj::interval<jkj::bigint::int_var, jkj::interval_type_t::bounded_closed>> {
    using type = json_member_list<json_number<"min", std::int_least64_t>,
                                  json_number<"max", std::int_least64_t>>;
};

template <>
struct daw::json::json_data_contract<jkj::idiv::floor_constraint_spec> {
    using type = json_member_list<
        json_class<"affine_coefficients", ::jkj::idiv::floor_constraint_spec::affine_transform>,
        json_class<"range", ::jkj::interval<::jkj::bigint::int_var,
                                            ::jkj::interval_type_t::bounded_closed>>>;
};

struct input_params {
    continued_fraction_generator x;
    continued_fraction_generator y;
    std::vector<jkj::idiv::floor_constraint_spec> constraints;
};

template <>
struct daw::json::json_data_contract<input_params> {
    using type = json_member_list<json_class<"x", continued_fraction_generator>,
                                  json_class<"y", continued_fraction_generator>,
                                  json_array<"constraints", ::jkj::idiv::floor_constraint_spec>>;
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
    auto params = [] {
        std::ifstream input_file{"xi_zeta_region_json_example_input.json"};
        std::stringstream strstream;
        strstream << input_file.rdbuf();
        return daw::json::from_json<input_params>(strstream.str());
    }();

    auto const result = jkj::idiv::find_xi_zeta_region(params.x, params.y, params.constraints);

    std::visit(
        [](auto const& region) {
            static constexpr auto region_type = std::remove_cvref_t<decltype(region)>::region_type;

            namespace xi_zeta_region = jkj::idiv::xi_zeta_region;

            if constexpr (region_type == jkj::idiv::xi_zeta_region::region_type_t::entire_plane) {
                std::cout << "Entire plane.\n\n";
            }
            else if constexpr (region_type == xi_zeta_region::region_type_t::single_point) {
                std::cout << "(xi, zeta) = (" << region.xi.numerator << " / "
                          << region.xi.denominator << ", " << region.zeta.numerator << " / "
                          << region.zeta.denominator << ")\n\n";
            }
            else if constexpr (region_type == xi_zeta_region::region_type_t::line_segment) {
                std::cout << "For t in ";
                if (region.left_boundary_type == xi_zeta_region::boundary_type_t::inclusive) {
                    std::cout << "[";
                }
                else {
                    std::cout << "(";
                }
                std::cout << "0, 1";
                if (region.right_boundary_type == xi_zeta_region::boundary_type_t::inclusive) {
                    std::cout << "]";
                }
                else {
                    std::cout << "),\n    xi = (";
                }
                std::cout << region.base_point_xi.numerator << " / "
                          << region.base_point_xi.denominator << ") + t("
                          << region.direction_vector_xi.numerator << " / "
                          << region.direction_vector_xi.denominator << "),\n zeta = ("
                          << region.base_point_zeta.numerator << " / "
                          << region.base_point_zeta.denominator << ") + t("
                          << region.direction_vector_zeta.numerator << " / "
                          << region.direction_vector_zeta.denominator << ")\n\n";
            }
            else if constexpr (region_type ==
                               xi_zeta_region::region_type_t::infinite_parallelogram) {
                if (jkj::util::is_zero(region.value_gap)) {
                    std::cout << region.xi_coeff << " xi + " << region.zeta_coeff
                              << " zeta = " << region.min_value;
                }
                else {
                    std::cout << region.min_value;
                    if (region.lower_boundary_type == xi_zeta_region::boundary_type_t::inclusive) {
                        std::cout << " <= ";
                    }
                    else {
                        std::cout << " < ";
                    }
                    std::cout << region.xi_coeff << " xi + " << region.zeta_coeff << " zeta";
                    if (region.upper_boundary_type == xi_zeta_region::boundary_type_t::inclusive) {
                        std::cout << " <= ";
                    }
                    else {
                        std::cout << " < ";
                    }
                    std::cout << region.min_value + region.value_gap << "\n\n";
                }
            }
            else {
                region.for_each_vertical_slice([](auto const& slice) {
                    static constexpr auto slice_type =
                        std::remove_cvref_t<decltype(slice)>::slice_type;

                    if constexpr (slice_type ==
                                  xi_zeta_region::bounded_polygon::slice_type_t::single_point) {
                        std::cout << "When xi = " << slice.xi.numerator << " / "
                                  << slice.xi.denominator << ",\n  zeta = " << slice.zeta.numerator
                                  << " / " << slice.zeta.denominator << "\n\n";
                    }
                    else if constexpr (slice_type == xi_zeta_region::bounded_polygon::slice_type_t::
                                                         vertical_line_segment) {
                        std::cout << "When xi = " << slice.xi.numerator << " / "
                                  << slice.xi.denominator << ",\n  "
                                  << slice.zeta_range.lower_bound().numerator << " / "
                                  << slice.zeta_range.lower_bound().denominator;
                        if (slice.zeta_range.left_endpoint_type() == jkj::endpoint_type_t::open) {
                            std::cout << " < zeta ";
                        }
                        else {
                            std::cout << " <= zeta ";
                        }
                        if (slice.zeta_range.right_endpoint_type() == jkj::endpoint_type_t::open) {
                            std::cout << "< ";
                        }
                        else {
                            std::cout << "<= ";
                        }
                        std::cout << slice.zeta_range.upper_bound().numerator << " / "
                                  << slice.zeta_range.upper_bound().denominator << "\n\n";
                    }
                    else {
                        std::cout << "When " << slice.xi_range.lower_bound().numerator << " / "
                                  << slice.xi_range.lower_bound().denominator << " < xi < "
                                  << slice.xi_range.upper_bound().numerator << " / "
                                  << slice.xi_range.upper_bound().denominator << ",\n  ("
                                  << slice.lower_boundary_linear_coeff.numerator << " / "
                                  << slice.lower_boundary_linear_coeff.denominator << ") xi + ("
                                  << slice.lower_boundary_constant_coeff.numerator << " / "
                                  << slice.lower_boundary_constant_coeff.denominator << ")";
                        if (slice.lower_boundary_type ==
                            xi_zeta_region::boundary_type_t::inclusive) {
                            std::cout << " <= ";
                        }
                        else {
                            std::cout << " < ";
                        }
                        std::cout << "zeta";
                        if (slice.upper_boundary_type ==
                            xi_zeta_region::boundary_type_t::inclusive) {
                            std::cout << " <= ";
                        }
                        else {
                            std::cout << " < ";
                        }
                        std::cout << "(" << slice.upper_boundary_linear_coeff.numerator << " / "
                                  << slice.upper_boundary_linear_coeff.denominator << ") xi + ("
                                  << slice.upper_boundary_constant_coeff.numerator << " / "
                                  << slice.upper_boundary_constant_coeff.denominator << ")\n\n";
                    }
                });
            }
        },
        result);
}