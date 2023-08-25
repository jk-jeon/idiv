cmake_minimum_required(VERSION 3.19)

project(idivTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)
include(../cmake/windows-set-path.cmake)

if(PROJECT_IS_TOP_LEVEL)
  find_package(idiv REQUIRED)
  enable_testing()
endif()

add_executable(idiv_unit_test source/unit_test.cpp)
target_link_libraries(idiv_unit_test PRIVATE idiv::idiv)
target_compile_features(idiv_unit_test PRIVATE cxx_std_20)

add_test(NAME idiv_unit_test COMMAND idiv_unit_test)

add_folders(Test)
