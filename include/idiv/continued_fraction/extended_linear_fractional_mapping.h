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

#ifndef JKJ_HEADER_EXTENDED_LINEAR_FRACTIONAL_MAPPING
#define JKJ_HEADER_EXTENDED_LINEAR_FRACTIONAL_MAPPING

#include "engine/rational.h"

namespace jkj {
    namespace cntfrc {
        namespace detail {
            // Find the kernel of a rank-1 linear fractional transform.
            template <class Int, class NumNum, class DenNum, class NumDen, class DenDen>
            constexpr projective_rational<Int, Int> kernel_of_rank1_transform(
                linear_fractional_mapping<NumNum, DenNum, NumDen, DenDen> const& transform) {
                return util::is_zero(transform.num_to_num()) &&
                               util::is_zero(transform.den_to_num())
                           ? projective_rational{Int{-transform.den_to_den()},
                                                 Int{transform.num_to_den()}}
                           : projective_rational{Int{-transform.den_to_num()},
                                                 Int{transform.num_to_num()}};
            }
            template <class Int, class NumNum, class DenNum, class NumDen, class DenDen>
            constexpr projective_rational<Int, Int> kernel_of_rank1_transform(
                linear_fractional_mapping<NumNum, DenNum, NumDen, DenDen>&& transform) {
                return util::is_zero(transform.num_to_num()) &&
                               util::is_zero(transform.den_to_num())
                           ? projective_rational{Int{-std::move(transform).den_to_den()},
                                                 Int{std::move(transform).num_to_den()}}
                           : projective_rational{Int{-std::move(transform).den_to_num()},
                                                 Int{std::move(transform).num_to_num()}};
            }
            // Find the range of a rank-1 linear fractional transform.
            template <class Int, class NumNum, class DenNum, class NumDen, class DenDen>
            constexpr projective_rational<Int, Int> range_of_rank1_transform(
                linear_fractional_mapping<NumNum, DenNum, NumDen, DenDen> const& transform) {
                return util::is_zero(transform.num_to_num()) &&
                               util::is_zero(transform.num_to_den())
                           ? projective_rational{Int{transform.den_to_num()},
                                                 Int{transform.den_to_den()}}
                           : projective_rational{Int{transform.num_to_num()},
                                                 Int{transform.num_to_den()}};
            }
            template <class Int, class NumNum, class DenNum, class NumDen, class DenDen>
            constexpr projective_rational<Int, Int> range_of_rank1_transform(
                linear_fractional_mapping<NumNum, DenNum, NumDen, DenDen>&& transform) {
                return util::is_zero(transform.num_to_num()) &&
                               util::is_zero(transform.num_to_den())
                           ? projective_rational{Int{std::move(transform).den_to_num()},
                                                 Int{std::move(transform).den_to_den()}}
                           : projective_rational{Int{std::move(transform).num_to_num()},
                                                 Int{std::move(transform).num_to_den()}};
            }

            template <class ProjectiveRational>
            constexpr auto sign_normalize(ProjectiveRational&& x)
                -> projective_rational<std::remove_cvref_t<decltype(x.numerator)>,
                                       decltype(util::abs(x.denominator))> {
                using return_type = projective_rational<std::remove_cvref_t<decltype(x.numerator)>,
                                                        decltype(util::abs(x.denominator))>;
                if (util::is_strictly_negative(x.denominator)) {
                    return return_type{-static_cast<ProjectiveRational&&>(x).numerator,
                                       util::abs(static_cast<ProjectiveRational&&>(x).denominator)};
                }
                else {
                    return return_type{static_cast<ProjectiveRational&&>(x).numerator,
                                       util::abs(static_cast<ProjectiveRational&&>(x).denominator)};
                }
            }
        }

        template <class Int>
            requires(std::is_object_v<Int>)
        class extended_linear_fractional_mapping : private linear_fractional_mapping<Int> {
            using base_type = linear_fractional_mapping<Int>;
            using int_type = Int;
            using uint_type = decltype(util::abs(std::declval<Int>()));
            int determinant_sign_;

            template <auto allowed_interval_types>
            static constexpr auto compute_map_cyclic_interval_return_type() noexcept {
                static_assert(allowed_interval_types.size() != 0);

                constexpr auto flags = [] {
                    util::array<bool, variable_shape_cyclic_interval<projective_rational<
                                          int_type, uint_type>>::allowed_interval_types()
                                          .size()>
                        ret_value{};

                    for (std::size_t idx = 0; idx < allowed_interval_types.size(); ++idx) {
                        switch (allowed_interval_types[idx]) {
                            using enum cyclic_interval_type_t;

                        case empty:
                            ret_value[std::size_t(empty)] = true;
                            break;

                        case single_point:
                            ret_value[std::size_t(single_point)] = true;
                            break;

                        case open:
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(open)] = true;
                            break;

                        case left_open_right_closed:
                        case left_closed_right_open:
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(left_open_right_closed)] = true;
                            ret_value[std::size_t(left_closed_right_open)] = true;
                            break;

                        case closed:
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(closed)] = true;
                            break;

                        case single_complement:
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            break;

                        case entire:
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;
                        }
                    }
                    return ret_value;
                }();

                constexpr std::size_t number_of_possible_types = [flags] {
                    std::size_t ret_value = 0;
                    for (std::size_t idx = 0; idx < flags.size(); ++idx) {
                        if (flags[idx]) {
                            ++ret_value;
                        }
                    }
                    return ret_value;
                }();

                constexpr auto possible_types = [flags] {
                    util::array<std::size_t, number_of_possible_types> ret_value{};
                    std::size_t target_idx = 0;
                    for (std::size_t idx = 0; idx < flags.size(); ++idx) {
                        if (flags[idx]) {
                            ret_value[target_idx++] = idx;
                        }
                    }
                    return ret_value;
                }();

                constexpr auto ret_value =
                    []<std::size_t... possible_types_>(std::index_sequence<possible_types_...>) {
                        return std::type_identity<variable_shape_cyclic_interval<
                            projective_rational<int_type, uint_type>,
                            cyclic_interval_type_t(possible_types_)...>>{};
                    }(tmp::unpack_array<possible_types>{});

                return ret_value;
            }

            template <auto allowed_interval_types>
            using map_cyclic_interval_return_type =
                typename decltype(compute_map_cyclic_interval_return_type<
                                  allowed_interval_types>())::type;

            template <class ReturnType, class StaticShapeCyclicInterval>
            constexpr ReturnType
            map_cyclic_interval_impl(StaticShapeCyclicInterval const& itv) const {
                using enum cyclic_interval_type_t;
                using itv_type = StaticShapeCyclicInterval;

                if constexpr (itv_type::interval_type() == empty) {
                    return cyclic_interval<projective_rational<int_type, uint_type>, empty>{};
                }
                else if constexpr (itv_type::interval_type() == entire) {
                    return cyclic_interval<projective_rational<int_type, uint_type>, entire>{};
                }
                else if constexpr (itv_type::interval_type() == single_point) {
                    return cyclic_interval<projective_rational<int_type, uint_type>, single_point>{
                        detail::sign_normalize(base_type::operator()(itv.lower_bound()))};
                }
                else if constexpr (itv_type::interval_type() == single_complement) {
                    return cyclic_interval<projective_rational<int_type, uint_type>,
                                           single_complement>{
                        detail::sign_normalize(base_type::operator()(itv.lower_bound()))};
                }
                else {
                    auto lower_bound =
                        detail::sign_normalize(base_type::operator()(itv.lower_bound()));
                    auto upper_bound =
                        detail::sign_normalize(base_type::operator()(itv.upper_bound()));

                    if (determinant_sign_ < 0) {
                        using std::swap;
                        swap(lower_bound, upper_bound);
                    }

                    if constexpr (itv_type::interval_type() == open) {
                        return cyclic_interval<projective_rational<int_type, uint_type>, open>{
                            std::move(lower_bound), std::move(upper_bound)};
                    }
                    else if constexpr (itv_type::interval_type() == closed) {
                        return cyclic_interval<projective_rational<int_type, uint_type>, closed>{
                            std::move(lower_bound), std::move(upper_bound)};
                    }
                    else {
                        if ((itv_type::interval_type() == left_open_right_closed &&
                             determinant_sign_ > 0) ||
                            (itv_type::interval_type() == left_closed_right_open &&
                             determinant_sign_ < 0)) {
                            return cyclic_interval<projective_rational<int_type, uint_type>,
                                                   left_open_right_closed>{std::move(lower_bound),
                                                                           std::move(upper_bound)};
                        }
                        else {
                            return cyclic_interval<projective_rational<int_type, uint_type>,
                                                   left_closed_right_open>{std::move(lower_bound),
                                                                           std::move(upper_bound)};
                        }
                    }
                }
            }

        public:
            enum class kind_t { generic, constant };

            template <class NumNum, class DenNum, class NumDen, class DenDen>
            constexpr extended_linear_fractional_mapping(NumNum&& a, DenNum&& b, NumDen&& c,
                                                         DenDen&& d)
                : base_type{static_cast<NumNum&&>(a), static_cast<DenNum&&>(b),
                            static_cast<NumDen&&>(c), static_cast<DenDen&&>(d)},
                  determinant_sign_{base_type::determinant_sign()} {
                // Falls back to a constant function if the rank is 1.
                if (determinant_sign_ == 0) {
                    auto reduced_fraction = reduce_fraction(
                        detail::range_of_rank1_transform<Int>(static_cast<base_type&&>(*this)));

                    static_cast<base_type&>(*this) =
                        base_type{std::move(reduced_fraction).numerator, zero{},
                                  std::move(reduced_fraction).denominator, zero{}};
                }
            }

            template <class NumNum, class DenNum, class NumDen, class DenDen>
            constexpr extended_linear_fractional_mapping(
                linear_fractional_mapping<NumNum, DenNum, NumDen, DenDen> const& lft)
                : extended_linear_fractional_mapping(lft.num_to_num(), lft.den_to_num(),
                                                     lft.num_to_den(), lft.den_to_den()) {}

            template <class NumNum, class DenNum, class NumDen, class DenDen>
            constexpr extended_linear_fractional_mapping(
                linear_fractional_mapping<NumNum, DenNum, NumDen, DenDen>&& lft)
                : extended_linear_fractional_mapping(
                      std::move(lft).num_to_num(), std::move(lft).den_to_num(),
                      std::move(lft).num_to_den(), std::move(lft).den_to_den()) {}

            base_type const& coefficients() const noexcept {
                return static_cast<base_type const&>(*this);
            }
            using base_type::num_to_num;
            using base_type::den_to_num;
            using base_type::num_to_den;
            using base_type::den_to_den;

            template <class Num, class Den>
            constexpr projective_rational<Int, Int>
            operator()(projective_rational<Num, Den> const& x) const {
                if (determinant_sign_ != 0) {
                    return projective_rational<Int, Int>{base_type::operator()(x)};
                }
                else {
                    return projective_rational<Int, Int>{num_to_num(), num_to_den()};
                }
            }

            constexpr int determinant_sign() const noexcept { return determinant_sign_; }
            constexpr kind_t kind() const noexcept {
                return determinant_sign_ == 0 ? kind_t::constant : kind_t::generic;
            }

            // Compose the linear fractional transform given by the matrix (t s;0 t), which
            // corresponds to the translation map x |-> x + s/t. This map is undefined if
            // s/t = x = infinity, but when s/t = infinity, we can continuously extend the function
            // into the constant function x |-> infinity.
            // !! NOTE !! This function is NOT exception-safe.
            template <class Num, class Den = unity>
            constexpr void translate(Num&& numerator, Den&& denominator = Den{unity{}}) {
                if (util::is_zero(denominator)) {
                    static_cast<base_type&>(*this) = base_type{unity{}, zero{}, zero{}, zero{}};
                    determinant_sign_ = 0;
                }
                else {
                    base_type::num_to_num_ *= denominator;
                    base_type::num_to_num_ += numerator * num_to_den();
                    base_type::den_to_num_ *= denominator;
                    base_type::den_to_num_ += numerator * den_to_den();
                    base_type::num_to_den_ *= denominator;
                    base_type::den_to_den_ *= denominator;
                }
            }

            // Compose the linear fractional transform given by the matrix (0 1;1 0), which
            // corresponds to the reciprocal map x |-> 1/x.
            // !! NOTE !! This function is exception-safe ONLY IF Int is no-throw swappable.
            constexpr void reflect() noexcept(std::is_nothrow_swappable_v<Int>) {
                using std::swap;
                swap(base_type::num_to_num_, base_type::num_to_den_);
                swap(base_type::den_to_num_, base_type::den_to_den_);
                determinant_sign_ *= -1;
            }

            // Compute the image of a cyclic interval.
            template <class VariableShapeCyclicInterval>
            constexpr map_cyclic_interval_return_type<
                VariableShapeCyclicInterval::allowed_interval_types()>
            map_cyclic_interval(VariableShapeCyclicInterval const& itv) const {
                using return_type = map_cyclic_interval_return_type<
                    VariableShapeCyclicInterval::allowed_interval_types()>;

                if constexpr (return_type::is_allowed_interval_type(
                                  cyclic_interval_type_t::single_point)) {
                    // Rank = 1 means a constant function.
                    if (determinant_sign_ == 0) {
                        if constexpr (return_type::is_allowed_interval_type(
                                          cyclic_interval_type_t::empty)) {
                            if (itv.interval_type() == cyclic_interval_type_t::empty) {
                                return cyclic_interval<projective_rational<int_type, uint_type>,
                                                       cyclic_interval_type_t::empty>{};
                            }
                        }

                        return cyclic_interval<projective_rational<int_type, uint_type>,
                                               cyclic_interval_type_t::single_point>{
                            detail::sign_normalize(projective_rational{base_type::num_to_num(),
                                                                       base_type::num_to_den()})};
                    }
                    else {
                        return itv.visit([this](auto&& itv_) {
                            return this->map_cyclic_interval_impl<return_type>(itv_);
                        });
                    }
                }
                // The only case where single-point image is impossible is when the input interval
                // is empty.
                else {
                    return cyclic_interval<projective_rational<int_type, uint_type>,
                                           cyclic_interval_type_t::empty>{};
                }
            }
        };

        template <class Int>
            requires(std::is_object_v<Int>)
        class extended_bilinear_fractional_mapping : private bilinear_fractional_mapping<Int> {
        public:
            enum class kind_t {
                generic,
                nonconstant_unary_function_of_x,
                nonconstant_unary_function_of_y,
                constant
            };

        private:
            using base_type = bilinear_fractional_mapping<Int>;
            using int_type = Int;
            using uint_type = decltype(util::abs(std::declval<Int>()));

            // L1 = B^T R A - A^T R B.
            // L2 = B R A^T - A R B^T.
            struct determinant_form_t {
                Int a, b, d;
            };
            determinant_form_t det_form1_ = {Int{zero{}}, Int{zero{}}, Int{zero{}}};
            determinant_form_t det_form2_ = {Int{zero{}}, Int{zero{}}, Int{zero{}}};

            int numerator_determinant_sign_;
            int denominator_determinant_sign_;

            kind_t kind_ = kind_t::generic;

            int number_of_points_in_indeterminacy_locus_ = 0;

            template <bool away_from_indeterminacy_locus, auto allowed_interval_types_x,
                      auto allowed_interval_types_y>
            static constexpr auto compute_map_cyclic_rectangle_return_type() noexcept {
                static_assert(allowed_interval_types_x.size() != 0 &&
                              allowed_interval_types_y.size() != 0);

                constexpr auto total_possible_allowed_interval_types =
                    variable_shape_cyclic_interval<
                        projective_rational<int_type, uint_type>>::allowed_interval_types()
                        .size();

                constexpr auto pack = [](cyclic_interval_type_t first,
                                         cyclic_interval_type_t second) {
                    struct index_pair {
                        std::size_t larger;
                        std::size_t smaller;
                    };

                    auto first_idx = std::size_t(first);
                    auto second_idx = std::size_t(second);
                    auto idx_pair = first_idx > second_idx ? index_pair{first_idx, second_idx}
                                                           : index_pair{second_idx, first_idx};

                    return idx_pair.smaller + (idx_pair.larger * (idx_pair.larger + 1)) / 2;
                };

                constexpr auto allowed_combinations = [pack] {
                    struct allowed_combination_flags_t {
                        util::array<bool, (total_possible_allowed_interval_types *
                                           (total_possible_allowed_interval_types + 1)) /
                                              2>
                            flags;
                        std::size_t number_of_set;
                    };

                    constexpr auto allowed_combination_flags = [pack] {
                        allowed_combination_flags_t ret_value{};
                        for (std::size_t idx_x = 0; idx_x < allowed_interval_types_x.size();
                             ++idx_x) {
                            for (std::size_t idx_y = 0; idx_y < allowed_interval_types_y.size();
                                 ++idx_y) {
                                auto combination_idx = pack(allowed_interval_types_x[idx_x],
                                                            allowed_interval_types_y[idx_y]);
                                if (!ret_value.flags[combination_idx]) {
                                    ret_value.flags[combination_idx] = true;
                                    ++ret_value.number_of_set;
                                }
                            }
                        }
                        return ret_value;
                    }();

                    util::array<std::size_t, allowed_combination_flags.number_of_set> ret_value{};
                    std::size_t idx = 0;
                    for (std::size_t combination_idx = 0;
                         combination_idx < allowed_combination_flags.flags.size();
                         ++combination_idx) {
                        if (allowed_combination_flags.flags[combination_idx]) {
                            ret_value[idx++] = combination_idx;
                        }
                    }
                    return ret_value;
                }();

                constexpr auto flags = [pack, allowed_combinations] {
                    util::array<bool, variable_shape_cyclic_interval<projective_rational<
                                          int_type, uint_type>>::allowed_interval_types()
                                          .size()>
                        ret_value{};

                    for (std::size_t idx = 0; idx < allowed_combinations.size(); ++idx) {
                        switch (allowed_combinations[idx]) {
                            using enum cyclic_interval_type_t;
                        case pack(empty, empty):
                        case pack(empty, single_point):
                        case pack(empty, open):
                        case pack(empty, left_open_right_closed):
                        case pack(empty, left_closed_right_open):
                        case pack(empty, closed):
                        case pack(empty, single_complement):
                        case pack(empty, entire):
                            ret_value[std::size_t(empty)] = true;
                            break;

                        case pack(single_point, single_point):
                            if (!away_from_indeterminacy_locus) {
                                ret_value[std::size_t(empty)] = true;
                            }
                            ret_value[std::size_t(single_point)] = true;
                            break;

                        case pack(single_point, open):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(open)] = true;
                            break;

                        case pack(single_point, left_open_right_closed):
                        case pack(single_point, left_closed_right_open):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(left_open_right_closed)] = true;
                            ret_value[std::size_t(left_closed_right_open)] = true;
                            break;

                        case pack(single_point, closed):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(closed)] = true;
                            break;

                        case pack(single_point, single_complement):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            break;

                        case pack(single_point, entire):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(open, open):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(open)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(open, left_open_right_closed):
                        case pack(open, left_closed_right_open):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(open)] = true;
                            ret_value[std::size_t(left_open_right_closed)] = true;
                            ret_value[std::size_t(left_closed_right_open)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(open, closed):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(open)] = true;
                            ret_value[std::size_t(left_open_right_closed)] = true;
                            ret_value[std::size_t(left_closed_right_open)] = true;
                            ret_value[std::size_t(closed)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(open, single_complement):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(open)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(open, entire):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(open)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(left_open_right_closed, left_open_right_closed):
                        case pack(left_open_right_closed, left_closed_right_open):
                        case pack(left_closed_right_open, left_closed_right_open):
                            if (!away_from_indeterminacy_locus) {
                                ret_value[std::size_t(closed)] = true;
                            }
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(open)] = true;
                            ret_value[std::size_t(left_open_right_closed)] = true;
                            ret_value[std::size_t(left_closed_right_open)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(left_open_right_closed, closed):
                        case pack(left_closed_right_open, closed):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(left_open_right_closed)] = true;
                            ret_value[std::size_t(left_closed_right_open)] = true;
                            ret_value[std::size_t(open)] = true;
                            ret_value[std::size_t(closed)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(left_open_right_closed, single_complement):
                        case pack(left_closed_right_open, single_complement):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(left_open_right_closed)] = true;
                            ret_value[std::size_t(left_closed_right_open)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(left_open_right_closed, entire):
                        case pack(left_closed_right_open, entire):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(left_open_right_closed)] = true;
                            ret_value[std::size_t(left_closed_right_open)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(closed, closed):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(closed)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(closed, single_complement):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(closed)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(closed, entire):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(closed)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(single_complement, single_complement):
                        case pack(single_complement, entire):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(single_complement)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        case pack(entire, entire):
                            ret_value[std::size_t(single_point)] = true;
                            ret_value[std::size_t(entire)] = true;
                            break;

                        default:
                            util::constexpr_assert(false);
                        }
                    }
                    return ret_value;
                }();

                constexpr std::size_t number_of_possible_types = [flags] {
                    std::size_t ret_value = 0;
                    for (std::size_t idx = 0; idx < flags.size(); ++idx) {
                        if (flags[idx]) {
                            ++ret_value;
                        }
                    }
                    return ret_value;
                }();

                constexpr auto possible_types = [flags] {
                    util::array<std::size_t, number_of_possible_types> ret_value{};
                    std::size_t target_idx = 0;
                    for (std::size_t idx = 0; idx < flags.size(); ++idx) {
                        if (flags[idx]) {
                            ret_value[target_idx++] = idx;
                        }
                    }
                    return ret_value;
                }();

                constexpr auto ret_value =
                    []<std::size_t... possible_types_>(std::index_sequence<possible_types_...>) {
                        return std::type_identity<variable_shape_cyclic_interval<
                            projective_rational<int_type, uint_type>,
                            cyclic_interval_type_t(possible_types_)...>>{};
                    }(tmp::unpack_array<possible_types>{});

                return ret_value;
            }

            template <bool away_from_indeterminacy_locus, auto allowed_interval_types_x,
                      auto allowed_interval_types_y>
            using map_cyclic_rectangle_return_type =
                typename decltype(compute_map_cyclic_rectangle_return_type<
                                  away_from_indeterminacy_locus, allowed_interval_types_x,
                                  allowed_interval_types_y>())::type;

            template <class VariableShapeCyclicInterval>
            static constexpr auto compute_closure_type() noexcept {
                using enum cyclic_interval_type_t;
                constexpr bool is_empty_allowed =
                    VariableShapeCyclicInterval::is_allowed_interval_type(empty);
                constexpr bool is_single_point_allowed =
                    VariableShapeCyclicInterval::is_allowed_interval_type(single_point);
                constexpr bool is_closed_allowed =
                    VariableShapeCyclicInterval::is_allowed_interval_type(open) ||
                    VariableShapeCyclicInterval::is_allowed_interval_type(left_open_right_closed) ||
                    VariableShapeCyclicInterval::is_allowed_interval_type(left_closed_right_open) ||
                    VariableShapeCyclicInterval::is_allowed_interval_type(closed);
                constexpr bool is_entire_allowed =
                    VariableShapeCyclicInterval::is_allowed_interval_type(single_complement) ||
                    VariableShapeCyclicInterval::is_allowed_interval_type(entire);

                constexpr std::size_t number_of_allowed_interval_types =
                    std::size_t(is_empty_allowed) + std::size_t(is_single_point_allowed) +
                    std::size_t(is_closed_allowed) + std::size_t(is_entire_allowed);
                static_assert(number_of_allowed_interval_types != 0);

                constexpr auto possible_types = [] {
                    util::array<std::size_t, number_of_allowed_interval_types> ret_value{};
                    std::size_t idx = 0;
                    if constexpr (is_empty_allowed) {
                        ret_value[idx++] = std::size_t(empty);
                    }
                    if constexpr (is_single_point_allowed) {
                        ret_value[idx++] = std::size_t(single_point);
                    }
                    if constexpr (is_closed_allowed) {
                        ret_value[idx++] = std::size_t(closed);
                    }
                    if constexpr (is_entire_allowed) {
                        ret_value[idx++] = std::size_t(entire);
                    }
                    return ret_value;
                }();

                return []<std::size_t... possible_types_>(std::index_sequence<possible_types_...>) {
                    return std::type_identity<variable_shape_cyclic_interval<
                        typename VariableShapeCyclicInterval::value_type,
                        cyclic_interval_type_t(possible_types_)...>>{};
                }(tmp::unpack_array<possible_types>{});
            }

            template <class VariableShapeCyclicInterval, class StaticShapeCyclicInterval>
            static constexpr
                typename decltype(compute_closure_type<VariableShapeCyclicInterval>())::type
                closure_impl(StaticShapeCyclicInterval const& itv) {
                using enum cyclic_interval_type_t;
                using itv_type = StaticShapeCyclicInterval;
                using value_type = typename itv_type::value_type;

                if constexpr (itv_type::interval_type() == empty ||
                              itv_type::interval_type() == single_point ||
                              itv_type::interval_type() == closed ||
                              itv_type::interval_type() == entire) {
                    return itv;
                }
                else if constexpr (itv_type::interval_type() == single_complement) {
                    return cyclic_interval<value_type, entire>{};
                }
                else {
                    return cyclic_interval<value_type, closed>{itv.lower_bound(),
                                                               itv.upper_bound()};
                }
            }

            template <class ProjectiveRational>
            constexpr Int evaluate_det_form1(ProjectiveRational&& y) const {
                return det_form1_.a * y.numerator * y.numerator +
                       ((det_form1_.b * y.numerator * y.denominator) << 1) +
                       det_form1_.d * y.denominator * y.denominator;
            }
            template <class ProjectiveRational>
            constexpr Int evaluate_det_form2(ProjectiveRational&& y) const {
                return det_form2_.a * y.numerator * y.numerator +
                       ((det_form2_.b * y.numerator * y.denominator) << 1) +
                       det_form2_.d * y.denominator * y.denominator;
            }

            // Check existence of y in the given interval such that L1(y,y) = 0,
            // assuming det L1 <= 0.
            template <class LowerBound, class UpperBound>
            constexpr bool has_parallel_point(LowerBound&& u, boundary_type_t left_boundary_type,
                                              UpperBound&& v,
                                              boundary_type_t right_boundary_type) const {
                auto const Luu = util::strong_order_to_int(
                    evaluate_det_form1(static_cast<LowerBound&&>(u)) <=> 0);
                auto const Lvv = util::strong_order_to_int(
                    evaluate_det_form1(static_cast<UpperBound&&>(v)) <=> 0);
                auto const Luv =
                    util::strong_order_to_int((det_form1_.a * u.numerator * v.numerator +
                                               det_form1_.b * u.numerator * v.denominator +
                                               det_form1_.b * u.denominator * v.numerator +
                                               det_form1_.d * u.denominator * v.denominator) <=> 0);

                // If y exists in the open interval, then always return true.
                if (Luu * Lvv < 0 || Luu * Luv < 0 || Lvv * Luv < 0) {
                    return true;
                }

                // If y does not exist in the closed interval, then always return false.
                if (Luu * Lvv > 0 && Luu * Luv > 0) {
                    return false;
                }

                // Otherwise, the only possibility is that y is on the boundary.
                if (Luu == 0 && left_boundary_type == boundary_type_t::closed) {
                    return true;
                }
                if (Lvv == 0 && right_boundary_type == boundary_type_t::closed) {
                    return true;
                }
                return false;
            }

            template <class FirstLowerBound, class FirstUpperBound, class SecondLowerBound,
                      class SecondUpperBound>
            constexpr bool intersects_indeterminacy_locus_interval_pair(
                FirstLowerBound&& a, boundary_type_t first_left_boundary_type,   //
                FirstUpperBound&& b, boundary_type_t first_right_boundary_type,  //
                SecondLowerBound&& c, boundary_type_t second_left_boundary_type, //
                SecondUpperBound&& d, boundary_type_t second_right_boundary_type) const {
                // 16 different cases.
                using enum boundary_type_t;
                if (a == c) {
                    if (b == d) {
                        // [a,b] cap [c,d] = [a,b].
                        return has_parallel_point(
                            static_cast<FirstLowerBound&&>(a),
                            first_left_boundary_type && second_left_boundary_type,
                            static_cast<FirstUpperBound&&>(b),
                            first_right_boundary_type && second_right_boundary_type);
                    }
                    // [a,b,d].
                    else if (cyclic_order(a, b, d)) {
                        // [a,b] cap [c,d] = [a,b].
                        return has_parallel_point(
                            static_cast<FirstLowerBound&&>(a),
                            first_left_boundary_type && second_left_boundary_type,
                            static_cast<FirstUpperBound&&>(b), first_right_boundary_type);
                    }
                    // [a,d,b].
                    else {
                        // [a,b] cap [c,d] = [a,d].
                        return has_parallel_point(
                            static_cast<FirstLowerBound&&>(a),
                            first_left_boundary_type && second_left_boundary_type,
                            static_cast<SecondUpperBound&&>(d), second_right_boundary_type);
                    }
                }
                else if (a == d) {
                    if (b == c) {
                        // [a,b] cap [c,d] = {a} cup {b}.
                        bool result = false;
                        if ((first_left_boundary_type && second_right_boundary_type) == closed) {
                            result = (result || util::is_zero(evaluate_det_form1(
                                                    static_cast<FirstLowerBound&&>(a))));
                        }
                        if ((first_right_boundary_type && second_left_boundary_type) == closed) {
                            result = (result || util::is_zero(evaluate_det_form1(
                                                    static_cast<FirstUpperBound&&>(b))));
                        }
                        return result;
                    }
                    // [a,b,c].
                    else if (cyclic_order(a, b, c)) {
                        // [a,b] cap [c,d] = {a}.
                        return ((first_left_boundary_type && second_right_boundary_type) ==
                                closed) &&
                               util::is_zero(evaluate_det_form1(static_cast<FirstLowerBound&&>(a)));
                    }
                    // [a,c,b].
                    else {
                        // [a,b] cap [c,d] = {a} cup [c,b].
                        bool result = false;
                        if ((first_left_boundary_type && second_right_boundary_type) == closed) {
                            result = (result || util::is_zero(evaluate_det_form1(
                                                    static_cast<FirstLowerBound&&>(a))));
                        }
                        result = (result || has_parallel_point(static_cast<SecondLowerBound&&>(c),
                                                               second_left_boundary_type,
                                                               static_cast<FirstUpperBound&&>(b),
                                                               first_right_boundary_type));
                        return result;
                    }
                }
                else if (b == d) {
                    // [a,b,c].
                    if (cyclic_order(a, b, c)) {
                        // [a,b] cap [c,d] = [a,b].
                        return has_parallel_point(
                            static_cast<FirstLowerBound&&>(a), first_left_boundary_type,
                            static_cast<FirstUpperBound&&>(b),
                            first_right_boundary_type && second_right_boundary_type);
                    }
                    // [a,c,b].
                    else {
                        // [a,b] cap [c,d] = [c,b].
                        return has_parallel_point(
                            static_cast<SecondLowerBound&&>(c), second_left_boundary_type,
                            static_cast<FirstUpperBound&&>(b),
                            first_right_boundary_type && second_right_boundary_type);
                    }
                }
                else if (b == c) {
                    // [a,b,d].
                    if (cyclic_order(a, b, d)) {
                        // [a,b] cap [c,d] = {b}.
                        return ((first_right_boundary_type && second_left_boundary_type) ==
                                closed) &&
                               util::is_zero(evaluate_det_form1(static_cast<FirstUpperBound&&>(b)));
                    }
                    // [a,d,b].
                    else {
                        // [a,b] cap [c,d] = {b} cup [a,d].
                        bool result = false;
                        if ((first_right_boundary_type && second_left_boundary_type) == closed) {
                            result = (result || util::is_zero(evaluate_det_form1(
                                                    static_cast<FirstUpperBound&&>(b))));
                        }
                        result = (result || has_parallel_point(static_cast<FirstLowerBound&&>(a),
                                                               first_left_boundary_type,
                                                               static_cast<SecondUpperBound&&>(d),
                                                               second_right_boundary_type));
                        return result;
                    }
                }
                else if (cyclic_order(a, b, c)) {
                    // [a,b,c] & [c,d,a] => [a,b,c,d].
                    if (cyclic_order(c, d, a)) {
                        // [a,b] cap [c,d] = empty.
                        return false;
                    }
                    // [b,d,c] & [c,a,b] => [a,b,d,c].
                    else if (cyclic_order(b, d, c)) {
                        // [a,b] cap [c,d] = [a,b].
                        return has_parallel_point(
                            static_cast<FirstLowerBound&&>(a), first_left_boundary_type,
                            static_cast<FirstUpperBound&&>(b), first_right_boundary_type);
                    }
                    // [d,b,c] & [c,a,d] => [a,d,b,c].
                    else {
                        // [a,b] cap [c,d] = [a,d].
                        return has_parallel_point(
                            static_cast<FirstLowerBound&&>(a), first_left_boundary_type,
                            static_cast<SecondUpperBound&&>(d), second_right_boundary_type);
                    }
                }
                else {
                    // [a,c,b] & [b,d,a] => [a,c,b,d].
                    if (cyclic_order(b, d, a)) {
                        // [a,b] cap [c,d] = [c,b].
                        return has_parallel_point(
                            static_cast<SecondLowerBound&&>(c), second_left_boundary_type,
                            static_cast<FirstUpperBound&&>(b), first_right_boundary_type);
                    }
                    // [a,c,d] & [d,b,a] => [a,c,d,b].
                    else if (cyclic_order(a, c, d)) {
                        // [a,b] cap [c,d] = [c,d].
                        return has_parallel_point(
                            static_cast<SecondLowerBound&&>(c), second_left_boundary_type,
                            static_cast<SecondUpperBound&&>(d), second_right_boundary_type);
                    }
                    // [a,d,c] & [c,b,a] => [a,d,c,b].
                    else {
                        // [a,b] cap [c,d] = [a,d] cup [c,b].
                        return has_parallel_point(static_cast<FirstLowerBound&&>(a),
                                                  first_left_boundary_type,
                                                  static_cast<SecondUpperBound&&>(d),
                                                  second_right_boundary_type) ||
                               has_parallel_point(
                                   static_cast<SecondLowerBound&&>(c), second_left_boundary_type,
                                   static_cast<SecondUpperBound&&>(d), second_right_boundary_type);
                    }
                }
            }

            enum class corner_index_t {
                bottom_left = 0,
                bottom_right = 1,
                top_right = 2,
                top_left = 3
            };
            template <class ValueTypeX, class ValueTypeY>
            struct corner_info_t {
                ValueTypeX const& horizontal_coord;
                ValueTypeX const& other_horizontal_coord;
                ValueTypeY const& vertical_coord;
                ValueTypeY const& other_vertical_coord;
                boundary_type_t horizontal_bdy_type;
                boundary_type_t other_horizontal_bdy_type;
                boundary_type_t vertical_bdy_type;
                boundary_type_t other_vertical_bdy_type;
            };

            template <class StaticShapeCyclicIntervalX, class StaticShapeCyclicIntervalY,
                      class CornerValue>
            constexpr bool
            should_remove_corner(StaticShapeCyclicIntervalX const&,
                                 StaticShapeCyclicIntervalY const& itv_y,
                                 corner_info_t<typename StaticShapeCyclicIntervalX::value_type,
                                               typename StaticShapeCyclicIntervalY::value_type> (
                                     &corner_info_arr)[4],
                                 CornerValue const& corner_value, corner_index_t corner_idx) const {
                auto a = corner_value.denominator * xnum_ynum_to_num() -
                         corner_value.numerator * xnum_ynum_to_den();
                auto minus_b = corner_value.numerator * xnum_yden_to_den() -
                               corner_value.denominator * xnum_yden_to_num();
                auto minus_c = corner_value.numerator * xden_ynum_to_den() -
                               corner_value.denominator * xden_ynum_to_num();
                auto d = corner_value.denominator * xden_yden_to_num() -
                         corner_value.numerator * xden_yden_to_den();

                auto determinant_sign = a * d <=> minus_b * minus_c;

                auto const& corner_info = corner_info_arr[std::size_t(corner_idx)];
                if (determinant_sign == 0) {
                    if ((corner_info.horizontal_bdy_type == boundary_type_t::closed &&
                         a * corner_info.horizontal_coord.numerator ==
                             minus_c * corner_info.horizontal_coord.denominator &&
                         d * corner_info.horizontal_coord.denominator ==
                             minus_b * corner_info.horizontal_coord.numerator) ||
                        (corner_info.vertical_bdy_type == boundary_type_t::closed &&
                         a * corner_info.vertical_coord.numerator ==
                             minus_b * corner_info.vertical_coord.denominator &&
                         d * corner_info.vertical_coord.denominator ==
                             minus_c * corner_info.vertical_coord.numerator)) {
                        return false;
                    }
                    else {
                        return true;
                    }
                }
                else if (((corner_idx == corner_index_t::bottom_left ||
                           corner_idx == corner_index_t::top_right) &&
                          determinant_sign > 0) ||
                         ((corner_idx == corner_index_t::bottom_right ||
                           corner_idx == corner_index_t::top_left) &&
                          determinant_sign < 0)) {
                    return false;
                }
                else {
                    auto endpoint_img = projective_rational<Int, Int>{
                        minus_b * corner_info.other_horizontal_coord.numerator -
                            d * corner_info.other_horizontal_coord.denominator,
                        a * corner_info.other_horizontal_coord.numerator -
                            minus_c * corner_info.other_horizontal_coord.denominator};

                    if (cyclic_order(itv_y.lower_bound(), endpoint_img, itv_y.upper_bound())) {
                        return false;
                    }
                    else if (endpoint_img == corner_info.other_vertical_coord) {
                        return corner_info.other_horizontal_bdy_type == boundary_type_t::open ||
                               corner_info.other_vertical_bdy_type == boundary_type_t::open;
                    }
                    else {
                        return true;
                    }
                }
            }

            template <bool away_from_indeterminacy_locus, class ReturnType,
                      class StaticShapeCyclicIntervalX, class StaticShapeCyclicIntervalY,
                      class StaticShapeCyclicIntervalEdgeUnion>
            constexpr ReturnType map_cyclic_rectangle_remove_corner_impl(
                StaticShapeCyclicIntervalX const& itv_x, StaticShapeCyclicIntervalY const& itv_y,
                StaticShapeCyclicIntervalEdgeUnion const& edge_union,
                corner_info_t<typename StaticShapeCyclicIntervalX::value_type,
                              typename StaticShapeCyclicIntervalY::value_type> (
                    &corner_info_arr)[4]) const {
                using enum cyclic_interval_type_t;
                using itv_x_type = StaticShapeCyclicIntervalX;
                using itv_y_type = StaticShapeCyclicIntervalY;
                using itv_type = StaticShapeCyclicIntervalEdgeUnion;
                using value_type = projective_rational<int_type, uint_type>;

                if constexpr (itv_type::interval_type() == single_point) {
                    return cyclic_interval<value_type, single_point>{
                        detail::sign_normalize(std::move(edge_union).lower_bound())};
                }
                else if constexpr (itv_type::interval_type() == entire) {
                    for (corner_index_t corner_idx = corner_index_t::bottom_left;
                         corner_idx != corner_index_t::top_left;
                         corner_idx = corner_index_t(std::size_t(corner_idx) + 1)) {
                        auto const& corner_info = corner_info_arr[std::size_t(corner_idx)];

                        if constexpr (!away_from_indeterminacy_locus) {
                            if (is_indeterminacy_locus(corner_info.horizontal_coord,
                                                       corner_info.vertical_coord)) {
                                continue;
                            }
                        }
                        auto corner_value = base_type::operator()(corner_info.horizontal_coord,
                                                                  corner_info.vertical_coord);

                        if (should_remove_corner(itv_x, itv_y, corner_info_arr, corner_value,
                                                 corner_idx)) {
                            return cyclic_interval<value_type, single_complement>{
                                detail::sign_normalize(std::move(corner_value))};
                        }
                    }
                    return cyclic_interval<value_type, entire>{};
                }
                else {
                    static_assert(itv_type::interval_type() == closed);

                    if constexpr (itv_x_type::interval_type() == open &&
                                  itv_y_type::interval_type() == open) {
                        return cyclic_interval<value_type, open>{
                            detail::sign_normalize(std::move(edge_union).lower_bound()),
                            detail::sign_normalize(std::move(edge_union).upper_bound())};
                    }
                    else {
                        bool remove_lower_bound = false;
                        bool remove_upper_bound = false;
                        for (corner_index_t corner_idx = corner_index_t::bottom_left;
                             corner_idx != corner_index_t::top_left;
                             corner_idx = corner_index_t(std::size_t(corner_idx) + 1)) {
                            auto const& corner_info = corner_info_arr[std::size_t(corner_idx)];

                            if constexpr (!away_from_indeterminacy_locus) {
                                if (is_indeterminacy_locus(corner_info.horizontal_coord,
                                                           corner_info.vertical_coord)) {
                                    continue;
                                }
                            }
                            auto corner_value = base_type::operator()(corner_info.horizontal_coord,
                                                                      corner_info.vertical_coord);

                            if (!remove_lower_bound && corner_value == edge_union.lower_bound()) {
                                if (should_remove_corner(itv_x, itv_y, corner_info_arr,
                                                         corner_value, corner_idx)) {
                                    remove_lower_bound = true;
                                }
                                continue;
                            }
                            if (!remove_upper_bound && corner_value == edge_union.upper_bound()) {
                                if (should_remove_corner(itv_x, itv_y, corner_info_arr,
                                                         corner_value, corner_idx)) {
                                    remove_upper_bound = true;
                                }
                                continue;
                            }
                        } // for loop.

                        if (remove_lower_bound && remove_upper_bound) {
                            return cyclic_interval<value_type, open>{
                                detail::sign_normalize(std::move(edge_union).lower_bound()),
                                detail::sign_normalize(std::move(edge_union).upper_bound())};
                        }
                        else if (remove_lower_bound) {
                            return cyclic_interval<value_type, left_open_right_closed>{
                                detail::sign_normalize(std::move(edge_union).lower_bound()),
                                detail::sign_normalize(std::move(edge_union).upper_bound())};
                        }
                        else {
                            // When both of the intervals are half-open,
                            // then the resulting interval can be closed
                            // only if the rectangle touches the
                            // indeterminacy locus at its boundary.
                            if constexpr (away_from_indeterminacy_locus &&
                                          itv_x_type::interval_type() != closed &&
                                          itv_y_type::interval_type() != closed) {
                                return cyclic_interval<value_type, left_closed_right_open>{
                                    detail::sign_normalize(std::move(edge_union).lower_bound()),
                                    detail::sign_normalize(std::move(edge_union).upper_bound())};
                            }
                            else {
                                if (remove_upper_bound) {
                                    return cyclic_interval<value_type, left_closed_right_open>{
                                        detail::sign_normalize(std::move(edge_union).lower_bound()),
                                        detail::sign_normalize(
                                            std::move(edge_union).upper_bound())};
                                }
                                else {
                                    return cyclic_interval<value_type, itv_type::interval_type()>{
                                        detail::sign_normalize(std::move(edge_union).lower_bound()),
                                        detail::sign_normalize(
                                            std::move(edge_union).upper_bound())};
                                }
                            }
                        }
                    }
                }
            }

            template <bool away_from_indeterminacy_locus, class ReturnType,
                      class StaticShapeCyclicIntervalX, class StaticShapeCyclicIntervalY>
            constexpr ReturnType
            map_cyclic_rectangle_impl(StaticShapeCyclicIntervalX const& itv_x,
                                      StaticShapeCyclicIntervalY const& itv_y) const {
                using enum cyclic_interval_type_t;
                using itv_x_type = StaticShapeCyclicIntervalX;
                using itv_y_type = StaticShapeCyclicIntervalY;
                using value_type = projective_rational<int_type, uint_type>;

                // Either is empty.
                if constexpr (itv_x_type::interval_type() == empty ||
                              itv_y_type::interval_type() == empty) {
                    return cyclic_interval<value_type, empty>{};
                }

                // Second interval is a point.
                else if constexpr (itv_y_type::interval_type() == single_point) {
                    auto specialized = extended_linear_fractional_mapping<Int>{
                        base_type::specialize_y(itv_y.lower_bound())};

                    if constexpr (away_from_indeterminacy_locus ||
                                  itv_x_type::interval_type() != single_point) {
                        return specialized.map_cyclic_interval(itv_x);
                    }
                    else {
                        if (specialized.determinant_sign() != 0) {
                            specialized.map_cyclic_interval(itv_x);
                        }
                        else {
                            auto constant_value =
                                detail::kernel_of_rank1_transform<Int>(specialized.coefficients());
                            auto const associated_indeterminacy_locus_x =
                                util::is_zero(constant_value.numerator)
                                    ? projective_rational{-Int{std::move(specialized).den_to_den()},
                                                          Int{std::move(specialized).num_to_den()}}
                                    : projective_rational{-Int{std::move(specialized).den_to_num()},
                                                          Int{std::move(specialized).num_to_num()}};

                            if (itv_x.lower_bound() == associated_indeterminacy_locus_x) {
                                return cyclic_interval<value_type, empty>{};
                            }
                            else {
                                return cyclic_interval<value_type, single_point>{
                                    detail::sign_normalize(std::move(constant_value))};
                            }
                        }
                    }
                }

                // First interval is a point and the second interval is not empty nor a
                // point.
                else if constexpr (itv_x_type::interval_type() == single_point) {
                    return extended_linear_fractional_mapping<Int>{
                        base_type::specialize_x(itv_x.lower_bound())}
                        .map_cyclic_interval(itv_y);
                }

                // Either is RP1 and both are not empty nor singleton.
                else if constexpr (itv_x_type::interval_type() == entire ||
                                   itv_y_type::interval_type() == entire) {
                    return cyclic_interval<value_type, entire>{};
                }

                // Second interval is the complemenet of a point, first interval is not
                // empty nor a point nor RP1.
                else if constexpr (itv_y_type::interval_type() == single_complement) {
                    if (!util::is_zero(evaluate_det_form1(itv_y.lower_bound()))) {
                        return cyclic_interval<value_type, entire>{};
                    }
                    else {
                        auto potentially_removed_point = detail::range_of_rank1_transform<Int>(
                            base_type::specialize_y(itv_y.lower_bound()));
                        auto inverse_image =
                            detail::kernel_of_rank1_transform<Int>(linear_fractional_mapping{
                                potentially_removed_point.denominator * xnum_ynum_to_num() -
                                    potentially_removed_point.numerator * xnum_ynum_to_den(),
                                potentially_removed_point.denominator * xden_ynum_to_num() -
                                    potentially_removed_point.numerator * xden_ynum_to_den(),
                                potentially_removed_point.denominator * xnum_yden_to_num() -
                                    potentially_removed_point.numerator * xnum_yden_to_den(),
                                potentially_removed_point.denominator * xden_yden_to_num() -
                                    potentially_removed_point.numerator * xden_yden_to_den()});

                        if (itv_x.contains(inverse_image)) {
                            return cyclic_interval<value_type, entire>{};
                        }
                        else {
                            return cyclic_interval<value_type, single_complement>{
                                detail::sign_normalize(std::move(potentially_removed_point))};
                        }
                    }
                }

                // First interval is the complemenet of a point, second interval is a
                // nontrivial bounded interval.
                else if constexpr (itv_x_type::interval_type() == single_complement) {
                    if (!util::is_zero(evaluate_det_form2(itv_x.lower_bound()))) {
                        return cyclic_interval<value_type, entire>{};
                    }
                    else {
                        auto potentially_removed_point = detail::range_of_rank1_transform<Int>(
                            base_type::specialize_x(itv_x.lower_bound()));
                        auto inverse_image =
                            detail::kernel_of_rank1_transform<Int>(linear_fractional_mapping{
                                potentially_removed_point.denominator * xnum_ynum_to_num() -
                                    potentially_removed_point.numerator * xnum_ynum_to_den(),
                                potentially_removed_point.denominator * xnum_yden_to_num() -
                                    potentially_removed_point.numerator * xnum_yden_to_den(),
                                potentially_removed_point.denominator * xden_ynum_to_num() -
                                    potentially_removed_point.numerator * xden_ynum_to_den(),
                                potentially_removed_point.denominator * xden_yden_to_num() -
                                    potentially_removed_point.numerator * xden_yden_to_den()});

                        if (itv_y.contains(inverse_image)) {
                            return cyclic_interval<value_type, entire>{};
                        }
                        else {
                            return cyclic_interval<value_type, single_complement>{
                                detail::sign_normalize(std::move(potentially_removed_point))};
                        }
                    }
                }

                // Both are nontrivial bounded intervals.
                else {
                    using closed_cyclic_interval =
                        variable_shape_cyclic_interval<projective_rational<Int, Int>, single_point,
                                                       closed, entire>;

                    // Compute itv union [from,to].
                    // from should be inside itv and should be distinct from to.
                    auto append_to_right = [](auto&& itv, auto&& from, auto&& to,
                                              auto const& continuation) -> closed_cyclic_interval {
                        if (itv.interval_type() == single_point ||
                            cyclic_order(itv.upper_bound(), to, itv.lower_bound())) {
                            return continuation(
                                cyclic_interval<projective_rational<Int, Int>, closed>{
                                    static_cast<decltype(itv)&&>(itv).lower_bound(),
                                    static_cast<decltype(to)&&>(to)});
                        }
                        else if (itv.upper_bound() == from ||
                                 cyclic_order(from, itv.upper_bound(), to)) {
                            return cyclic_interval<projective_rational<Int, Int>, entire>{};
                        }
                        else {
                            return continuation(static_cast<decltype(itv)&&>(itv));
                        }
                    };
                    // Compute itv union [to,from].
                    // from should be inside itv and should be distinct from to.
                    auto append_to_left = [](auto&& itv, auto&& from, auto&& to,
                                             auto const& continuation) -> closed_cyclic_interval {
                        if (itv.interval_type() == single_point ||
                            cyclic_order(itv.upper_bound(), to, itv.lower_bound())) {
                            return continuation(
                                cyclic_interval<projective_rational<Int, Int>, closed>{
                                    static_cast<decltype(to)&&>(to),
                                    static_cast<decltype(itv)&&>(itv).upper_bound()});
                        }
                        else if (itv.lower_bound() == from ||
                                 cyclic_order(to, itv.lower_bound(), from)) {
                            return cyclic_interval<projective_rational<Int, Int>, entire>{};
                        }
                        else {
                            return continuation(static_cast<decltype(itv)&&>(itv));
                        }
                    };

                    // Compute the union of the images of the four edges in
                    // continuation-passing style.
                    auto append_bottom_edge =
                        [&](auto const& append_right_edge, auto const& append_top_edge,
                            auto const& append_left_edge) -> closed_cyclic_interval {
                        auto specialized = extended_linear_fractional_mapping<Int>{
                            base_type::specialize_y(itv_y.lower_bound())};

                        if (!util::is_zero(specialized.determinant_sign())) {
                            projective_rational<Int, Int> start = specialized(itv_x.lower_bound());
                            projective_rational<Int, Int> end = specialized(itv_x.upper_bound());
                            auto bottom_edge =
                                specialized.determinant_sign() > 0
                                    ? cyclic_interval<projective_rational<Int, Int>, closed>{start,
                                                                                             end}
                                    : cyclic_interval<projective_rational<Int, Int>, closed>{end,
                                                                                             start};
                            return append_right_edge(std::move(start), std::move(end),
                                                     std::move(bottom_edge), append_top_edge,
                                                     append_left_edge);
                        }
                        else {
                            auto k2 =
                                detail::range_of_rank1_transform<Int>(specialized.coefficients());

                            if constexpr (away_from_indeterminacy_locus) {
                                return append_right_edge(
                                    k2, k2,
                                    cyclic_interval<projective_rational<Int, Int>, single_point>{
                                        k2},
                                    append_top_edge, append_left_edge);
                            }
                            else {
                                auto x0 = detail::kernel_of_rank1_transform<Int>(specialized);
                                auto k1 = detail::range_of_rank1_transform<Int>(
                                    base_type::specialize_x(x0));

                                if (cyclic_order(itv_x.lower_bound(), x0, itv_x.upper_bound())) {
                                    return cyclic_interval<projective_rational<Int, Int>, entire>{};
                                }
                                else if (cyclic_order(itv_x.upper_bound(), x0,
                                                      itv_x.lower_bound())) {
                                    return append_right_edge(
                                        k2, k2,
                                        cyclic_interval<projective_rational<Int, Int>,
                                                        single_point>{k2},
                                        append_top_edge, append_left_edge);
                                }
                                else if (k1 != k2) {
                                    projective_rational<Int, Int> start =
                                        x0 == itv_x.lower_bound() ? k1 : k2;
                                    projective_rational<Int, Int> end =
                                        x0 == itv_x.lower_bound() ? k2 : k1;
                                    auto bottom_edge =
                                        util::is_nonnegative(
                                            det_form1_.b * (itv_y.lower_bound().numerator *
                                                                itv_y.lower_bound().numerator -
                                                            itv_y.lower_bound().denominator *
                                                                itv_y.lower_bound().denominator) +
                                            (det_form1_.d - det_form1_.a) *
                                                itv_y.lower_bound().numerator *
                                                itv_y.lower_bound().denominator)
                                            ? cyclic_interval<projective_rational<Int, Int>,
                                                              closed>{end, start}
                                            : cyclic_interval<projective_rational<Int, Int>,
                                                              closed>{start, end};
                                    return append_right_edge(std::move(start), std::move(end),
                                                             std::move(bottom_edge),
                                                             append_top_edge, append_left_edge);
                                }
                                else {
                                    auto vertical_sign = util::is_nonnegative(evaluate_det_form1(
                                        projective_rational<Int const&, Int const&>(
                                            -itv_y.lower_bound().denominator,
                                            itv_y.lower_bound().numerator)));
                                    auto horizontal_sign = util::is_nonnegative(evaluate_det_form2(
                                        projective_rational<Int const&, Int const&>(-x0.denominator,
                                                                                    x0.numerator)));

                                    if ((vertical_sign != horizontal_sign &&
                                         x0 == itv_x.lower_bound()) ||
                                        (vertical_sign == horizontal_sign &&
                                         x0 == itv_x.upper_bound())) {
                                        return cyclic_interval<projective_rational<Int, Int>,
                                                               entire>{};
                                    }
                                    else {
                                        return append_right_edge(
                                            k2, k2,
                                            cyclic_interval<projective_rational<Int, Int>,
                                                            single_point>{k2},
                                            append_top_edge, append_left_edge);
                                    }
                                }
                            }
                        }
                    }; // append_bottom_edge.

                    auto append_right_edge =
                        [this, &itv_x, &itv_y, &append_to_right,
                         &append_to_left](auto&& start, auto&& end, auto&& edge_union,
                                          auto const& append_top_edge,
                                          auto const& append_left_edge) -> closed_cyclic_interval {
                        auto specialized = extended_linear_fractional_mapping<Int>{
                            base_type::specialize_x(itv_x.upper_bound())};

                        if (util::is_zero(specialized.determinant_sign())) {
                            if constexpr (!away_from_indeterminacy_locus) {
                                auto const y0 = detail::kernel_of_rank1_transform<Int>(specialized);

                                if (cyclic_order(itv_y.lower_bound(), y0, itv_y.upper_bound())) {
                                    return cyclic_interval<projective_rational<Int, Int>, entire>{};
                                }
                            }
                            return append_top_edge(static_cast<decltype(start)&&>(start),
                                                   static_cast<decltype(end)&&>(end),
                                                   static_cast<decltype(edge_union)&&>(edge_union),
                                                   append_left_edge);
                        }
                        else {
                            auto new_end = specialized(itv_y.upper_bound());
                            auto continuation = [&](auto&& new_interval) {
                                return append_top_edge(static_cast<decltype(start)&&>(start),
                                                       std::move(new_end), std::move(new_interval),
                                                       append_left_edge);
                            };
                            if (specialized.determinant_sign() > 0) {
                                return append_to_right(
                                    static_cast<decltype(edge_union)&&>(edge_union),
                                    static_cast<decltype(end)&&>(end), new_end, continuation);
                            }
                            else {
                                return append_to_left(
                                    static_cast<decltype(edge_union)&&>(edge_union),
                                    static_cast<decltype(end)&&>(end), new_end, continuation);
                            }
                        }
                    }; // append_right_edge.

                    auto append_top_edge =
                        [this, &itv_x, &itv_y, &append_to_right,
                         &append_to_left](auto&& start, auto&& end, auto&& edge_union,
                                          auto const& append_left_edge) -> closed_cyclic_interval {
                        auto specialized = extended_linear_fractional_mapping<Int>{
                            base_type::specialize_y(itv_y.upper_bound())};

                        int edge_direction = specialized.determinant_sign();
                        auto new_end = specialized(itv_x.lower_bound());

                        if (edge_direction == 0) {
                            auto& k2 = new_end;

                            if constexpr (away_from_indeterminacy_locus) {
                                return append_left_edge(
                                    static_cast<decltype(start)&&>(start),
                                    static_cast<decltype(end)&&>(end),
                                    static_cast<decltype(edge_union)&&>(edge_union));
                            }
                            else {
                                auto x0 = detail::kernel_of_rank1_transform<Int>(specialized);
                                auto k1 = detail::range_of_rank1_transform<Int>(
                                    base_type::specialize_x(x0));

                                if (cyclic_order(itv_x.lower_bound(), x0, itv_x.upper_bound())) {
                                    return cyclic_interval<projective_rational<Int, Int>, entire>{};
                                }
                                else if (cyclic_order(itv_x.upper_bound(), x0,
                                                      itv_x.lower_bound())) {
                                    return append_left_edge(
                                        static_cast<decltype(start)&&>(start),
                                        static_cast<decltype(end)&&>(end),
                                        static_cast<decltype(edge_union)&&>(edge_union));
                                }
                                else if (k1 != k2) {
                                    if (x0 == itv_x.lower_bound()) {
                                        new_end = k1;
                                    }
                                    edge_direction =
                                        util::is_nonnegative(
                                            det_form1_.b * (itv_y.upper_bound().numerator *
                                                                itv_y.upper_bound().numerator -
                                                            itv_y.upper_bound().denominator *
                                                                itv_y.upper_bound().denominator) +
                                            (det_form1_.d - det_form1_.a) *
                                                itv_y.upper_bound().numerator *
                                                itv_y.upper_bound().denominator)
                                            ? -1
                                            : 1;
                                }
                                else {
                                    auto vertical_sign = util::is_nonnegative(evaluate_det_form1(
                                        projective_rational<Int const&, Int const&>(
                                            -itv_y.upper_bound().denominator,
                                            itv_y.upper_bound().numerator)));
                                    auto horizontal_sign = util::is_nonnegative(evaluate_det_form2(
                                        projective_rational<Int const&, Int const&>(-x0.denominator,
                                                                                    x0.numerator)));

                                    if ((vertical_sign != horizontal_sign &&
                                         x0 == itv_x.upper_bound()) ||
                                        (vertical_sign == horizontal_sign &&
                                         x0 == itv_x.lower_bound())) {
                                        return cyclic_interval<projective_rational<Int, Int>,
                                                               entire>{};
                                    }
                                    else {
                                        return append_left_edge(
                                            static_cast<decltype(start)&&>(start),
                                            static_cast<decltype(end)&&>(end),
                                            static_cast<decltype(edge_union)&&>(edge_union));
                                    }
                                }
                            }
                        }

                        auto continuation = [&](auto&& new_interval) {
                            return append_left_edge(static_cast<decltype(start)&&>(start),
                                                    std::move(new_end), std::move(new_interval));
                        };
                        if (edge_direction < 0) {
                            return append_to_right(static_cast<decltype(edge_union)&&>(edge_union),
                                                   static_cast<decltype(end)&&>(end), new_end,
                                                   continuation);
                        }
                        else {
                            return append_to_left(static_cast<decltype(edge_union)&&>(edge_union),
                                                  static_cast<decltype(end)&&>(end), new_end,
                                                  continuation);
                        }
                    }; // append_top_edge.

                    auto append_left_edge =
                        [this, &itv_x, &itv_y, &append_to_right, &append_to_left](
                            auto&& start, auto&& end, auto&& edge_union) -> closed_cyclic_interval {
                        // To suppress warning.
                        static_cast<void>(itv_y);

                        auto specialized = extended_linear_fractional_mapping<Int>{
                            base_type::specialize_x(itv_x.lower_bound())};

                        if (util::is_zero(specialized.determinant_sign())) {
                            if constexpr (!away_from_indeterminacy_locus) {
                                auto const y0 = detail::kernel_of_rank1_transform<Int>(specialized);

                                if (cyclic_order(itv_y.lower_bound(), y0, itv_y.upper_bound())) {
                                    return cyclic_interval<projective_rational<Int, Int>, entire>{};
                                }
                            }
                            return static_cast<decltype(edge_union)&&>(edge_union);
                        }
                        else {
                            auto continuation = [&](auto&& new_interval) {
                                return std::move(new_interval);
                            };
                            if (specialized.determinant_sign() < 0) {
                                return append_to_right(
                                    static_cast<decltype(edge_union)&&>(edge_union),
                                    static_cast<decltype(end)&&>(end),
                                    static_cast<decltype(start)&&>(start), continuation);
                            }
                            else {
                                return append_to_left(
                                    static_cast<decltype(edge_union)&&>(edge_union),
                                    static_cast<decltype(end)&&>(end),
                                    static_cast<decltype(start)&&>(start), continuation);
                            }
                        }
                    }; // append_left_edge.

                    auto edge_union =
                        append_bottom_edge(append_right_edge, append_top_edge, append_left_edge);

                    // Remove boundary points if necessary.
                    if constexpr (itv_x_type::interval_type() != closed ||
                                  itv_y_type::interval_type() != closed) {

                        corner_info_t<typename itv_x_type::value_type,
                                      typename itv_y_type::value_type>
                            corner_info_arr[4] = {
                                {itv_x.lower_bound(), itv_x.upper_bound(), itv_y.lower_bound(),
                                 itv_y.upper_bound(), itv_x.left_boundary_type(),
                                 itv_x.right_boundary_type(), itv_y.left_boundary_type(),
                                 itv_y.right_boundary_type()},
                                {itv_x.upper_bound(), itv_x.lower_bound(), itv_y.lower_bound(),
                                 itv_y.upper_bound(), itv_x.right_boundary_type(),
                                 itv_x.left_boundary_type(), itv_y.left_boundary_type(),
                                 itv_y.right_boundary_type()},
                                {itv_x.upper_bound(), itv_x.lower_bound(), itv_y.upper_bound(),
                                 itv_y.lower_bound(), itv_x.right_boundary_type(),
                                 itv_x.left_boundary_type(), itv_y.right_boundary_type(),
                                 itv_y.left_boundary_type()},
                                {itv_x.lower_bound(), itv_x.upper_bound(), itv_y.upper_bound(),
                                 itv_y.lower_bound(), itv_x.left_boundary_type(),
                                 itv_x.right_boundary_type(), itv_y.right_boundary_type(),
                                 itv_y.left_boundary_type()}};

                        return std::move(edge_union).visit([&](auto&& edge_union_) -> ReturnType {
                            return map_cyclic_rectangle_remove_corner_impl<
                                away_from_indeterminacy_locus, ReturnType>(
                                itv_x, itv_y, edge_union_, corner_info_arr);
                        });
                    }
                    else {
                        return std::move(edge_union).visit([](auto&& itv) -> ReturnType {
                            using itv_type = std::remove_cvref_t<decltype(itv)>;
                            if constexpr (itv_type::interval_type() == single_point) {
                                return cyclic_interval<value_type, single_point>{
                                    detail::sign_normalize(std::move(itv).lower_bound())};
                            }
                            else if constexpr (itv_type::interval_type() == closed) {
                                return cyclic_interval<value_type, closed>{
                                    detail::sign_normalize(std::move(itv).lower_bound()),
                                    detail::sign_normalize(std::move(itv).upper_bound())};
                            }
                            else {
                                return cyclic_interval<value_type, entire>{};
                            }
                        });
                    }
                }
            }

        public:
            template <class XNumYNumNum, class XNumYDenNum, class XDenYNumNum, class XDenYDenNum,
                      class XNumYNumDen, class XNumYDenDen, class XDenYNumDen, class XDenYDenDen>
            constexpr extended_bilinear_fractional_mapping(XNumYNumNum&& a, XNumYDenNum&& b,
                                                           XDenYNumNum&& c, XDenYDenNum&& d,
                                                           XNumYNumDen&& e, XNumYDenDen&& f,
                                                           XDenYNumDen&& g, XDenYDenDen&& h)
                : base_type{static_cast<XNumYNumNum&&>(a), static_cast<XNumYDenNum&&>(b),
                            static_cast<XDenYNumNum&&>(c), static_cast<XDenYDenNum&&>(d),
                            static_cast<XNumYNumDen&&>(e), static_cast<XNumYDenDen&&>(f),
                            static_cast<XDenYNumDen&&>(g), static_cast<XDenYDenDen&&>(h)},
                  numerator_determinant_sign_{
                      linear_fractional_mapping{xnum_ynum_to_num(), xnum_yden_to_num(),
                                                xden_ynum_to_num(), xden_yden_to_num()}
                          .determinant_sign()},
                  denominator_determinant_sign_{
                      linear_fractional_mapping{xnum_ynum_to_den(), xnum_yden_to_den(),
                                                xden_ynum_to_den(), xden_yden_to_den()}
                          .determinant_sign()} {
                // Determine the kind and the number of points in the indeterminacy locus.

                // Step 1. Check for linear dependence of the numerator and the denominator.
                Int const* num_ptrs[] = {&xnum_ynum_to_num(), &xnum_yden_to_num(),
                                         &xden_ynum_to_num(), &xden_yden_to_num()};
                Int const* den_ptrs[] = {&xnum_ynum_to_den(), &xnum_yden_to_den(),
                                         &xden_ynum_to_den(), &xden_yden_to_den()};

                std::size_t first_nonzero_idx = 0;
                for (; first_nonzero_idx < 4; ++first_nonzero_idx) {
                    if (!util::is_zero(*num_ptrs[first_nonzero_idx]) ||
                        !util::is_zero(*den_ptrs[first_nonzero_idx])) {
                        break;
                    }
                }
                util::constexpr_assert(first_nonzero_idx != 4);

                bool is_linearly_dependent = true;
                for (std::size_t idx = first_nonzero_idx + 1; idx < 4; ++idx) {
                    if (*num_ptrs[first_nonzero_idx] * *den_ptrs[idx] !=
                        *num_ptrs[idx] * *den_ptrs[first_nonzero_idx]) {
                        is_linearly_dependent = false;
                        break;
                    }
                }

                // If they are linearly dependent, simplify the form of the function.
                if (is_linearly_dependent) {
                    kind_ = kind_t::constant;
                    auto reduced_fraction =
                        reduce_fraction(*num_ptrs[first_nonzero_idx], *den_ptrs[first_nonzero_idx]);
                    static_cast<base_type&>(*this) =
                        base_type{std::move(reduced_fraction.numerator),   zero{}, zero{}, zero{},
                                  std::move(reduced_fraction.denominator), zero{}, zero{}, zero{}};
                    return;
                }

                // Step 2. Represent the numerator and the denominator into simple tensors if
                // possible, and in that case check for independence on x or y.
                if (numerator_determinant_sign_ == 0 && denominator_determinant_sign_ == 0) {
                    // Both the numerator and the denominator must be of rank 1, because if one of
                    // them is of rank 0, then they are linearly dependent.
                    // Any rank-1 matrix (a b;c d) can be written as one of the forms
                    // (1/a) (a;c)(a b), (1/b) (b;d)(a b), (1/c) (a;c)(c d), (1/d) (b;d)(c d),
                    // where each of the above four is valid whenever the divisor is not zero.
                    struct simple_tensor_representation {
                        // (1/divisor) (first_row;second_row)(first_column second_column)
                        Int const& divisor;
                        Int const& first_row;
                        Int const& second_row;
                        Int const& first_column;
                        Int const& second_column;
                    };
                    auto decompose = [](Int const& a_, Int const& b_, Int const& c_,
                                        Int const& d_) -> simple_tensor_representation {
                        if (!util::is_zero(a_)) {
                            return {a_, a_, c_, a_, b_};
                        }
                        else if (!util::is_zero(b_)) {
                            return {b_, b_, d_, a_, b_};
                        }
                        else if (!util::is_zero(c_)) {
                            return {c_, a_, c_, c_, d_};
                        }
                        else {
                            util::constexpr_assert(!util::is_zero(d_));
                            return {d_, b_, d_, c_, d_};
                        }
                    };
                    auto num_decomposed =
                        decompose(*num_ptrs[0], *num_ptrs[1], *num_ptrs[2], *num_ptrs[3]);
                    auto den_decomposed =
                        decompose(*den_ptrs[0], *den_ptrs[1], *den_ptrs[2], *den_ptrs[3]);

                    // Independent to x if column vectors are linearly dependent.
                    if (num_decomposed.first_row * den_decomposed.second_row ==
                        num_decomposed.second_row * den_decomposed.first_row) {
                        kind_ = kind_t::nonconstant_unary_function_of_y;

                        auto reduced_fraction =
                            util::is_zero(num_decomposed.first_row) &&
                                    util::is_zero(den_decomposed.first_row)
                                ? reduce_fraction(
                                      num_decomposed.second_row * den_decomposed.divisor,
                                      den_decomposed.second_row * num_decomposed.divisor)
                                : reduce_fraction(num_decomposed.first_row * den_decomposed.divisor,
                                                  den_decomposed.first_row *
                                                      num_decomposed.divisor);

                        // Simplify the form of the function: we only take row vectors and the
                        // divisors from the decompositions.
                        static_cast<base_type&>(*this) = base_type{
                            Int{num_decomposed.first_column * reduced_fraction.numerator},
                            Int{num_decomposed.second_column * reduced_fraction.numerator},
                            zero{},
                            zero{},
                            Int{den_decomposed.first_column * reduced_fraction.denominator},
                            Int{den_decomposed.second_column * reduced_fraction.denominator},
                            zero{},
                            zero{}};
                        return;
                    }

                    // Independent to y if row vectors are linearly dependent.
                    if (num_decomposed.first_column * den_decomposed.second_column ==
                        num_decomposed.second_column * den_decomposed.first_column) {
                        kind_ = kind_t::nonconstant_unary_function_of_x;

                        auto reduced_fraction =
                            util::is_zero(num_decomposed.first_column) &&
                                    util::is_zero(den_decomposed.first_column)
                                ? reduce_fraction(
                                      num_decomposed.second_column * den_decomposed.divisor,
                                      den_decomposed.second_column * num_decomposed.divisor)
                                : reduce_fraction(
                                      num_decomposed.first_column * den_decomposed.divisor,
                                      den_decomposed.first_column * num_decomposed.divisor);

                        // Simplify the form of the function: we only take column vectors and the
                        // divisors from the decompositions.
                        static_cast<base_type&>(*this) = base_type{
                            Int{num_decomposed.first_row * reduced_fraction.numerator},    zero{},
                            Int{num_decomposed.second_row * reduced_fraction.numerator},   zero{},
                            Int{den_decomposed.first_row * reduced_fraction.denominator},  zero{},
                            Int{den_decomposed.second_row * reduced_fraction.denominator}, zero{}};
                        return;
                    }
                } // if (is_num_singular && is_den_singular)

                // Step 3. Compute the bilinear forms
                // L1 = B^T R A - A^T R B = 2sym(ag-ce  bg-de \\ ah-cf  bh-df) and
                // L2 = B R A^T - A R B^T = 2sym(af-be  cf-de \\ ah-bg  ch-dg).
                // At this point, the function does not admit any proper continuous extension, and
                // the number of points in the indeterminacy locus is determined by the determinant
                // of L1.
                det_form1_ = {((xnum_ynum_to_num() * xden_ynum_to_den() -
                                xden_ynum_to_num() * xnum_ynum_to_den())
                               << 1),
                              xnum_ynum_to_num() * xden_yden_to_den() +
                                  xnum_yden_to_num() * xden_ynum_to_den() -
                                  xden_ynum_to_num() * xnum_yden_to_den() -
                                  xden_yden_to_num() * xnum_ynum_to_den(),
                              ((xnum_yden_to_num() * xden_yden_to_den() -
                                xden_yden_to_num() * xnum_yden_to_den())
                               << 1)};
                det_form2_ = {((xnum_ynum_to_num() * xnum_yden_to_den() -
                                xnum_yden_to_num() * xnum_ynum_to_den())
                               << 1),
                              xnum_ynum_to_num() * xden_yden_to_den() -
                                  xnum_yden_to_num() * xden_ynum_to_den() +
                                  xden_ynum_to_num() * xnum_yden_to_den() -
                                  xden_yden_to_num() * xnum_ynum_to_den(),
                              ((xden_ynum_to_num() * xden_yden_to_den() -
                                xden_yden_to_num() * xden_ynum_to_den())
                               << 1)};

                auto det_L1_sign = (det_form1_.a * det_form1_.d) <=> (det_form1_.b * det_form1_.b);
                if (det_L1_sign == 0) {
                    number_of_points_in_indeterminacy_locus_ = 1;
                }
                else if (det_L1_sign < 0) {
                    number_of_points_in_indeterminacy_locus_ = 2;
                }
            }

            template <class XNumYNumNum, class XNumYDenNum, class XDenYNumNum, class XDenYDenNum,
                      class XNumYNumDen, class XNumYDenDen, class XDenYNumDen, class XDenYDenDen>
            constexpr extended_bilinear_fractional_mapping(
                bilinear_fractional_mapping<XNumYNumNum, XNumYDenNum, XDenYNumNum, XDenYDenNum,
                                            XNumYNumDen, XNumYDenDen, XDenYNumDen,
                                            XDenYDenDen> const& blft)
                : extended_bilinear_fractional_mapping(
                      blft.xnum_ynum_to_num(), blft.xnum_yden_to_num(), blft.xden_ynum_to_num(),
                      blft.xden_yden_to_num(), blft.xnum_ynum_to_den(), blft.xnum_yden_to_den(),
                      blft.xden_ynum_to_den(), blft.xden_yden_to_den()) {}

            template <class XNumYNumNum, class XNumYDenNum, class XDenYNumNum, class XDenYDenNum,
                      class XNumYNumDen, class XNumYDenDen, class XDenYNumDen, class XDenYDenDen>
            constexpr extended_bilinear_fractional_mapping(
                bilinear_fractional_mapping<XNumYNumNum, XNumYDenNum, XDenYNumNum, XDenYDenNum,
                                            XNumYNumDen, XNumYDenDen, XDenYNumDen, XDenYDenDen>&&
                    blft)
                : extended_bilinear_fractional_mapping(
                      std::move(blft).xnum_ynum_to_num(), std::move(blft).xnum_yden_to_num(),
                      std::move(blft).xden_ynum_to_num(), std::move(blft).xden_yden_to_num(),
                      std::move(blft).xnum_ynum_to_den(), std::move(blft).xnum_yden_to_den(),
                      std::move(blft).xden_ynum_to_den(), std::move(blft).xden_yden_to_den()) {}

            base_type const& coefficients() const noexcept {
                return static_cast<base_type const&>(*this);
            }
            using base_type::xnum_ynum_to_num;
            using base_type::xnum_yden_to_num;
            using base_type::xden_ynum_to_num;
            using base_type::xden_yden_to_num;
            using base_type::xnum_ynum_to_den;
            using base_type::xnum_yden_to_den;
            using base_type::xden_ynum_to_den;
            using base_type::xden_yden_to_den;

            template <class Num, class Den>
            constexpr projective_rational<Int, Int>
            operator()(projective_rational<Num, Den> const& x) const {
                switch (kind_) {
                    using enum kind_t;
                case nonconstant_unary_function_of_x:
                    return projective_rational<Int, Int>{
                        linear_fractional_mapping{xnum_ynum_to_num(), xden_ynum_to_num(),
                                                  xnum_ynum_to_den(), xden_yden_to_den()}(x)};

                case nonconstant_unary_function_of_y:
                    return projective_rational<Int, Int>{
                        linear_fractional_mapping{xnum_ynum_to_num(), xnum_yden_to_num(),
                                                  xnum_ynum_to_den(), xnum_yden_to_den()}(x)};

                case constant:
                    return projective_rational<Int, Int>{xnum_ynum_to_num(), xnum_ynum_to_den()};

                default:
                    util::constexpr_assert(kind_ == generic);
                    return projective_rational<Int, Int>{base_type::operator()(x)};
                }
            }

            constexpr kind_t kind() const noexcept { return kind_; }
            constexpr bool is_globally_well_defined() const noexcept {
                return number_of_points_in_indeterminacy_locus_ == 0;
            }
            constexpr int number_of_points_in_indeterminacy_locus() const noexcept {
                return number_of_points_in_indeterminacy_locus_;
            }

            // Check if the specified point is in the indeterminacy locus.
            template <class XNum, class XDen, class YNum, class YDen>
            constexpr bool is_indeterminacy_locus(projective_rational<XNum, XDen> const& x,
                                                  projective_rational<YNum, YDen> const& y) const {
                if (is_globally_well_defined()) {
                    return false;
                }

                return util::is_zero(xnum_ynum_to_num() * x.numerator * y.numerator +
                                     xnum_yden_to_num() * x.numerator * y.denominator +
                                     xden_ynum_to_num() * x.denominator * y.numerator +
                                     xden_yden_to_num() * x.denominator * y.denominator) &&
                       util::is_zero(xnum_ynum_to_den() * x.numerator * y.numerator +
                                     xnum_yden_to_den() * x.numerator * y.denominator +
                                     xden_ynum_to_den() * x.denominator * y.numerator +
                                     xden_yden_to_den() * x.denominator * y.denominator);
            }

            // Check if the specified rectangle intersects with the indeterminacy locus.
            template <class VariableShapeCyclicIntervalX, class VariableShapeCyclicIntervalY>
            constexpr bool
            intersects_indeterminacy_locus(VariableShapeCyclicIntervalX const& itv_x,
                                           VariableShapeCyclicIntervalY const& itv_y) const {
                if (is_globally_well_defined()) {
                    return false;
                }

                // When both the numerator and the denominator are of rank 1.
                if (numerator_determinant_sign_ == 0 && denominator_determinant_sign_ == 0) {
                    util::constexpr_assert(
                        !util::is_zero(xnum_ynum_to_num()) || !util::is_zero(xnum_yden_to_num()) ||
                        !util::is_zero(xden_ynum_to_num()) || !util::is_zero(xden_yden_to_num()));
                    util::constexpr_assert(
                        !util::is_zero(xnum_ynum_to_den()) || !util::is_zero(xnum_yden_to_den()) ||
                        !util::is_zero(xden_ynum_to_den()) || !util::is_zero(xden_yden_to_den()));

                    // The kernel of the numerator.
                    auto y1 =
                        util::is_zero(xnum_ynum_to_num()) && util::is_zero(xnum_yden_to_num())
                            ? projective_rational<Int, Int>{-xden_yden_to_num(), xden_ynum_to_num()}
                            : projective_rational<Int, Int>{-xnum_yden_to_num(),
                                                            xnum_ynum_to_num()};
                    auto RBy1 = projective_rational<Int, Int>{
                        -xden_ynum_to_den() * y1.numerator - xden_yden_to_den() * y1.denominator,
                        xnum_ynum_to_den() * y1.numerator + xnum_yden_to_den() * y1.denominator};

                    if (itv_x.contains(RBy1) && itv_y.contains(y1)) {
                        return true;
                    }

                    // The kernel of the denominator.
                    auto y2 =
                        util::is_zero(xnum_ynum_to_den()) && util::is_zero(xnum_yden_to_den())
                            ? projective_rational<Int, Int>{-xden_yden_to_den(), xden_ynum_to_den()}
                            : projective_rational<Int, Int>{-xnum_yden_to_den(),
                                                            xnum_ynum_to_den()};
                    auto RAy2 = projective_rational<Int, Int>{
                        -xden_ynum_to_num() * y2.numerator - xden_yden_to_num() * y2.denominator,
                        xnum_ynum_to_num() * y2.numerator + xnum_yden_to_num() * y2.denominator};

                    if (itv_x.contains(RAy2) && itv_y.contains(y2)) {
                        return true;
                    }

                    return false;
                } // if (numerator_determinant_sign_ == 0 && denominator_determinant_sign_ == 0)
                else {
                    // T_{A^-1 R}[I] if A is invertible, T_{B^-1 R}[I] otherwise.
                    auto transformed_itv_x =
                        extended_linear_fractional_mapping<Int>{
                            numerator_determinant_sign_ != 0
                                ? linear_fractional_mapping{-xnum_yden_to_num(),
                                                            -xden_yden_to_num(), xnum_ynum_to_num(),
                                                            xden_ynum_to_num()}
                                : linear_fractional_mapping{-xnum_yden_to_den(),
                                                            -xden_yden_to_den(), xnum_ynum_to_den(),
                                                            xden_ynum_to_den()}}
                            .map_cyclic_interval(itv_x);

                    // Compute the intersection, and check each component of the intersection.
                    return transformed_itv_x.visit([this, &itv_y](auto const& first_itv) {
                        using enum cyclic_interval_type_t;
                        using first_itv_type = std::remove_cvref_t<decltype(first_itv)>;

                        return itv_y.visit([this, &first_itv](auto const& second_itv) {
                            using second_itv_type = std::remove_cvref_t<decltype(second_itv)>;

                            // To suppress warning.
                            static_cast<void>(first_itv);

                            // Either of them is empty.
                            if constexpr (first_itv_type::interval_type() == empty ||
                                          second_itv_type::interval_type() == empty) {
                                // Intersection is empty.
                                return false;
                            }
                            // First interval is RP1, second interval is nonempty.
                            else if constexpr (first_itv_type::interval_type() == entire) {
                                if constexpr (second_itv_type::interval_type() == entire) {
                                    // Intersection is RP1.
                                    return true;
                                }
                                else if constexpr (second_itv_type::interval_type() ==
                                                   single_point) {
                                    // Intersection is a point.
                                    return util::is_zero(
                                        evaluate_det_form1(second_itv.lower_bound()));
                                }
                                else if constexpr (second_itv_type::interval_type() ==
                                                   single_complement) {
                                    // Intersection is the complement of a point.
                                    return number_of_points_in_indeterminacy_locus_ == 2 ||
                                           !util::is_zero(
                                               evaluate_det_form1(second_itv.lower_bound()));
                                }
                                else {
                                    // Intersection is a nontrivial bounded interval.
                                    return has_parallel_point(
                                        second_itv.lower_bound(), second_itv.left_boundary_type(),
                                        second_itv.upper_bound(), second_itv.right_boundary_type());
                                }
                            }
                            // First interval is not empty nor RP1, second interval is RP1.
                            else if constexpr (second_itv_type::interval_type() == entire) {
                                if constexpr (first_itv_type::interval_type() == single_point) {
                                    // Intersection is a single point.
                                    return util::is_zero(
                                        evaluate_det_form1(first_itv.lower_bound()));
                                }
                                else if constexpr (first_itv_type::interval_type() ==
                                                   single_complement) {
                                    // Intersection is the complement of a point.
                                    return number_of_points_in_indeterminacy_locus_ == 2 ||
                                           !util::is_zero(
                                               evaluate_det_form1(first_itv.lower_bound()));
                                }
                                else {
                                    // Intersection is a nontrivial bounded interval.
                                    return has_parallel_point(
                                        first_itv.lower_bound(), first_itv.left_boundary_type(),
                                        first_itv.upper_bound(), first_itv.right_boundary_type());
                                }
                            }
                            // First interval is a point, second interval is not empty nor RP1.
                            else if constexpr (first_itv_type::interval_type() == single_point) {
                                // Intersection is either empty or a point.
                                if (second_itv.contains(first_itv.lower_bound())) {
                                    // Point.
                                    return util::is_zero(
                                        evaluate_det_form1(first_itv.lower_bound()));
                                }
                                else {
                                    // Empty.
                                    return false;
                                }
                            }
                            // First interval is not empty nor RP1 nor a point, second interval is a
                            // point.
                            else if constexpr (second_itv_type::interval_type() == single_point) {
                                // Intersection is either empty or a point.
                                if (first_itv.contains(second_itv.lower_bound())) {
                                    // Point.
                                    return util::is_zero(
                                        evaluate_det_form1(second_itv.lower_bound()));
                                }
                                else {
                                    // Empty.
                                    return false;
                                }
                            }
                            // First interval is the complement of a point, second interval is not
                            // empty nor RP1 nor a point.
                            else if constexpr (first_itv_type::interval_type() ==
                                               single_complement) {
                                if constexpr (second_itv_type::interval_type() ==
                                              single_complement) {
                                    // Intersection is the complement of either one point or two
                                    // points.
                                    if (first_itv.lower_bound() == second_itv.lower_bound()) {
                                        // Complement of one point.
                                        return number_of_points_in_indeterminacy_locus_ == 2 ||
                                               !util::is_zero(
                                                   evaluate_det_form1(first_itv.lower_bound()));
                                    }
                                    else {
                                        // Complement of two points.
                                        auto count_in_complement =
                                            (util::is_zero(
                                                 evaluate_det_form1(first_itv.lower_bound()))
                                                 ? 1
                                                 : 0) +
                                            (util::is_zero(
                                                 evaluate_det_form1(second_itv.lower_bound()))
                                                 ? 1
                                                 : 0);
                                        return count_in_complement <
                                               number_of_points_in_indeterminacy_locus_;
                                    }
                                }
                                else {
                                    // Intersection is either a nontrivial bounded interval or a
                                    // nontrivial bounded interval minus a point.
                                    if (cyclic_order(second_itv.lower_bound(),
                                                     first_itv.lower_bound(),
                                                     second_itv.upper_bound())) {
                                        // Nontrivial interval minus a point in the interior.
                                        return has_parallel_point(second_itv.lower_bound(),
                                                                  second_itv.left_boundary_type(),
                                                                  first_itv.lower_bound(),
                                                                  boundary_type_t::open) ||
                                               has_parallel_point(first_itv.lower_bound(),
                                                                  boundary_type_t::open,
                                                                  second_itv.upper_bound(),
                                                                  second_itv.right_boundary_type());
                                    }
                                    else {
                                        // Nontrivial interval.
                                        return has_parallel_point(
                                            second_itv.lower_bound(),
                                            second_itv.lower_bound() == first_itv.lower_bound()
                                                ? boundary_type_t::open
                                                : second_itv.left_boundary_type(),
                                            second_itv.upper_bound(),
                                            second_itv.upper_bound() == first_itv.lower_bound()
                                                ? boundary_type_t::open
                                                : second_itv.right_boundary_type());
                                    }
                                }
                            }
                            // First interval is a nontrivial bounded interval, second interval is
                            // the complement of a point.
                            else if constexpr (second_itv_type::interval_type() ==
                                               single_complement) {
                                // Intersection is either a nontrivial bounded interval or a
                                // nontrivial bounded interval minus a point.
                                if (cyclic_order(first_itv.lower_bound(), second_itv.lower_bound(),
                                                 first_itv.upper_bound())) {
                                    // Nontrivial interval minus a point in the interior.
                                    return has_parallel_point(first_itv.lower_bound(),
                                                              first_itv.left_boundary_type(),
                                                              second_itv.lower_bound(),
                                                              boundary_type_t::open) ||
                                           has_parallel_point(second_itv.lower_bound(),
                                                              boundary_type_t::open,
                                                              first_itv.upper_bound(),
                                                              first_itv.right_boundary_type());
                                }
                                else {
                                    // Nontrivial bounded interval.
                                    return has_parallel_point(
                                        first_itv.lower_bound(),
                                        first_itv.lower_bound() == second_itv.lower_bound()
                                            ? boundary_type_t::open
                                            : first_itv.left_boundary_type(),
                                        first_itv.upper_bound(),
                                        first_itv.upper_bound() == second_itv.lower_bound()
                                            ? boundary_type_t::open
                                            : first_itv.right_boundary_type());
                                }
                            }
                            // Both are nontrivial bounded intervals.
                            else {
                                return intersects_indeterminacy_locus_interval_pair(
                                    first_itv.lower_bound(), first_itv.left_boundary_type(),
                                    first_itv.upper_bound(), first_itv.right_boundary_type(),
                                    second_itv.lower_bound(), second_itv.left_boundary_type(),
                                    second_itv.upper_bound(), second_itv.right_boundary_type());
                            }
                        });
                    });
                }
            }

            // Check if there is an indeterminacy locus in the closure of the specified rectangle.
            template <class VariableShapeCyclicIntervalX, class VariableShapeCyclicIntervalY>
            constexpr bool
            touches_indeterminacy_locus(VariableShapeCyclicIntervalX const& itv_x,
                                        VariableShapeCyclicIntervalY const& itv_y) const {
                return intersects_indeterminacy_locus(
                    itv_x.visit([](auto const& itv) {
                        return closure_impl<VariableShapeCyclicIntervalX>(itv);
                    }),
                    itv_y.visit([](auto const& itv) {
                        return closure_impl<VariableShapeCyclicIntervalY>(itv);
                    }));
            }

            // Compose the linear fractional transform given by the matrix (t s;0 t), which
            // corresponds to the translation map x |-> x + s/t. This map is undefined if
            // s/t = x = infinity, but when s/t = infinity, we can continuously extend the function
            // into the constant function x |-> infinity. This composition changes L1, L2 into
            // t^2 L1, t^2 L2, but unless t = 0 this does not affect any use of L1, L2.
            // !! NOTE !! This function is NOT exception-safe.
            template <class Num, class Den = unity>
            constexpr void translate(Num&& numerator, Den&& denominator = Den{unity{}}) {
                if (util::is_zero(denominator)) {
                    static_cast<base_type&>(*this) =
                        base_type{unity{}, zero{}, zero{}, zero{}, zero{}, zero{}, zero{}, zero{}};
                    kind_ = kind_t::constant;
                    number_of_points_in_indeterminacy_locus_ = 0;
                    det_form1_ = {Int{zero{}}, Int{zero{}}, Int{zero{}}};
                    det_form2_ = {Int{zero{}}, Int{zero{}}, Int{zero{}}};
                }
                else {
                    base_type::xnum_ynum_to_num_ *= denominator;
                    base_type::xnum_ynum_to_num_ += numerator * base_type::xnum_ynum_to_den();
                    base_type::xnum_yden_to_num_ *= denominator;
                    base_type::xnum_yden_to_num_ += numerator * base_type::xnum_yden_to_den();
                    base_type::xden_ynum_to_num_ *= denominator;
                    base_type::xden_ynum_to_num_ += numerator * base_type::xden_ynum_to_den();
                    base_type::xden_yden_to_num_ *= denominator;
                    base_type::xden_yden_to_num_ += numerator * base_type::xden_yden_to_den();
                    base_type::xnum_ynum_to_den_ *= denominator;
                    base_type::xnum_yden_to_den_ *= denominator;
                    base_type::xden_ynum_to_den_ *= denominator;
                    base_type::xden_yden_to_den_ *= denominator;
                }
            }

            // Compose the linear fractional transform given by the matrix (0 1;1 0), which
            // corresponds to the reciprocal map x |-> 1/x. This composition changes L1, L2 into
            // -L1, -L2, but this does not affect any use of L1, L2.
            // !! NOTE !! This function is exception-safe ONLY IF Int is no-throw swappable.
            constexpr void reflect() noexcept(std::is_nothrow_swappable_v<Int>) {
                using std::swap;
                swap(base_type::xnum_ynum_to_num_, base_type::xnum_ynum_to_den_);
                swap(base_type::xnum_yden_to_num_, base_type::xnum_yden_to_den_);
                swap(base_type::xden_ynum_to_num_, base_type::xden_ynum_to_den_);
                swap(base_type::xden_yden_to_num_, base_type::xden_yden_to_den_);
            }

            // Compute the image of a cyclic rectangle.
            template <bool away_from_indeterminacy_locus = false,
                      class VariableShapeCyclicIntervalX, class VariableShapeCyclicIntervalY>
            constexpr map_cyclic_rectangle_return_type<
                away_from_indeterminacy_locus,
                VariableShapeCyclicIntervalX::allowed_interval_types(),
                VariableShapeCyclicIntervalY::allowed_interval_types()>
            map_cyclic_rectangle(VariableShapeCyclicIntervalX const& itv_x,
                                 VariableShapeCyclicIntervalY const& itv_y) const {
                using return_type = map_cyclic_rectangle_return_type<
                    away_from_indeterminacy_locus,
                    VariableShapeCyclicIntervalX::allowed_interval_types(),
                    VariableShapeCyclicIntervalY::allowed_interval_types()>;

                using enum cyclic_interval_type_t;

                if (kind_ == kind_t::constant) {
                    if constexpr (return_type::is_allowed_interval_type(empty)) {
                        if (itv_x.interval_type() == empty || itv_y.interval_type() == empty) {
                            return cyclic_interval<projective_rational<int_type, uint_type>,
                                                   empty>{};
                        }
                        return cyclic_interval<projective_rational<int_type, uint_type>,
                                               single_point>{detail::sign_normalize(
                            projective_rational{xnum_ynum_to_num(), xnum_ynum_to_den()})};
                    }
                    else {
                        return cyclic_interval<projective_rational<int_type, uint_type>,
                                               single_point>{detail::sign_normalize(
                            projective_rational{xnum_ynum_to_num(), xnum_ynum_to_den()})};
                    }
                }
                else if (kind_ == kind_t::nonconstant_unary_function_of_x) {
                    return extended_linear_fractional_mapping<Int>{
                        xnum_ynum_to_num(), xden_ynum_to_num(), xnum_ynum_to_den(),
                        xden_ynum_to_den()}
                        .map_cyclic_interval(itv_x);
                }
                else if (kind_ == kind_t::nonconstant_unary_function_of_y) {
                    return extended_linear_fractional_mapping<Int>{
                        xnum_ynum_to_num(), xnum_yden_to_num(), xnum_ynum_to_den(),
                        xnum_yden_to_den()}
                        .map_cyclic_interval(itv_y);
                }
                else {
                    return itv_x.visit([this, &itv_y](auto const& itv_x_) -> return_type {
                        return itv_y.visit([this, &itv_x_](auto const& itv_y_) -> return_type {
                            return this->map_cyclic_rectangle_impl<away_from_indeterminacy_locus,
                                                                   return_type>(itv_x_, itv_y_);
                        });
                    });
                }
            }
        };
    }
}

#endif
