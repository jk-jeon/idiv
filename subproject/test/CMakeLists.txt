cmake_minimum_required(VERSION 3.19)

project(idivTests LANGUAGES CXX)

include(../../cmake/project-is-top-level.cmake)
include(../../cmake/folders.cmake)
include(../../cmake/windows-set-path.cmake)

include(FetchContent)
if (NOT TARGET ut)
  FetchContent_Declare(
          ut
          GIT_REPOSITORY https://github.com/boost-ext/ut.git
          GIT_TAG        v1.1.9
  )
  FetchContent_MakeAvailable(ut)
endif()

if(PROJECT_IS_TOP_LEVEL)
  find_package(idiv REQUIRED)
  enable_testing()
endif()

add_executable(idiv_unit_test source/integer_util.cpp
                              source/tmp.cpp
                              source/wuint.cpp
                              source/bigint.cpp
                              source/continued_fraction_mixin_metaprogramming.cpp
                              source/rational_continued_fraction.cpp
                              source/gosper_algorithm.cpp
                              source/log_continued_fraction.cpp
                              source/best_rational_approx.cpp
                              source/idiv.cpp
                              source/main.cpp)

target_link_libraries(idiv_unit_test PRIVATE idiv::idiv)
target_link_libraries(idiv_unit_test PRIVATE Boost::ut)
target_compile_features(idiv_unit_test PRIVATE cxx_std_20)

add_test(NAME idiv_unit_test COMMAND idiv_unit_test)

add_folders(Test)
