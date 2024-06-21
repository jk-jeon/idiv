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

void integer_util_test();
void tmp_test();
void wuint_test();
void bigint_test();
void continued_fraction_mixin_metaprogramming_test();
void rational_continued_fraction_test();
void gosper_algorithm_test();
void log_continued_fraction_test();
void best_rational_approx_test();
void idiv_test();

int main() {
    integer_util_test();
    tmp_test();
    wuint_test();
    bigint_test();
    continued_fraction_mixin_metaprogramming_test();
    rational_continued_fraction_test();
    gosper_algorithm_test();
    log_continued_fraction_test();
    best_rational_approx_test();
    idiv_test();
}
