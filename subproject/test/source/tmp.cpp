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

#include <idiv/util.h>
#include <boost/ut.hpp>

void tmp_test() {
    using namespace boost::ut;
    using namespace jkj;
    
    "[tmp]"_test = [] {
        using typelist = tmp::typelist<int, double, float, char, float, double, unsigned int, long>;
        should("find_first_index") = [] {
            expect(tmp::find_first_index<int>(typelist{}) == 0);
            expect(tmp::find_first_index<double>(typelist{}) == 1);
            expect(tmp::find_first_index<float>(typelist{}) == 2);
            expect(tmp::find_first_index<char>(typelist{}) == 3);
            expect(tmp::find_first_index<unsigned int>(typelist{}) == 6);
            expect(tmp::find_first_index<long>(typelist{}) == 7);
            expect(tmp::find_first_index<short>(typelist{}) == 8);
        };

        should("get_type") = [] {
            expect(std::is_same_v<tmp::get_type<0, typelist>, int>);
            expect(std::is_same_v<tmp::get_type<1, typelist>, double>);
            expect(std::is_same_v<tmp::get_type<2, typelist>, float>);
            expect(std::is_same_v<tmp::get_type<3, typelist>, char>);
            expect(std::is_same_v<tmp::get_type<4, typelist>, float>);
            expect(std::is_same_v<tmp::get_type<5, typelist>, double>);
            expect(std::is_same_v<tmp::get_type<6, typelist>, unsigned int>);
            expect(std::is_same_v<tmp::get_type<7, typelist>, long>);
        };

        should("back_sublist") = [] {
            expect(std::is_same_v<tmp::back_sublist<0, typelist>, tmp::typelist<>>);
            expect(std::is_same_v<tmp::back_sublist<1, typelist>, tmp::typelist<long>>);
            expect(
                std::is_same_v<tmp::back_sublist<2, typelist>, tmp::typelist<unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<3, typelist>,
                                  tmp::typelist<double, unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<4, typelist>,
                                  tmp::typelist<float, double, unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<5, typelist>,
                                  tmp::typelist<char, float, double, unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<6, typelist>,
                                  tmp::typelist<float, char, float, double, unsigned int, long>>);
            expect(std::is_same_v<
                   tmp::back_sublist<7, typelist>,
                   tmp::typelist<double, float, char, float, double, unsigned int, long>>);
            expect(std::is_same_v<tmp::back_sublist<8, typelist>, typelist>);
        };

        should("remove_duplicate") = [] {
            expect(std::is_same_v<tmp::remove_duplicate<typelist>,
                                  tmp::typelist<int, double, float, char, unsigned int, long>>);
        };

        should("push_back") = [] {
            expect(std::is_same_v<tmp::push_back<typelist, unsigned short>,
                                  tmp::typelist<int, double, float, char, float, double,
                                                unsigned int, long, unsigned short>>);
        };

        should("join") = [] {
            expect(std::is_same_v<
                   tmp::join<typelist, tmp::typelist<unsigned short, unsigned short>, typelist>,
                   tmp::typelist<int, double, float, char, float, double, unsigned int, long,
                                 unsigned short, unsigned short, int, double, float, char, float,
                                 double, unsigned int, long>>);
        };

        should("filter") = [] {
            auto predicate = [](auto arg) {
                return std::is_integral_v<typename decltype(arg)::type>;
            };

            expect(std::is_same_v<tmp::filter<typelist, decltype(predicate)>,
                                  tmp::typelist<int, char, unsigned int, long>>);
        };
    };
}
