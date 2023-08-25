cmake_minimum_required(VERSION 3.19)

project(idivExamples CXX)

include(../../cmake/project-is-top-level.cmake)
include(../../cmake/folders.cmake)

if(PROJECT_IS_TOP_LEVEL)
  find_package(idiv REQUIRED)
endif()

add_custom_target(run-examples)

function(add_example NAME)
  add_executable("${NAME}" "${NAME}.cpp")
  target_link_libraries("${NAME}" PRIVATE idiv::idiv)
  target_compile_features("${NAME}" PRIVATE cxx_std_20)
  add_custom_target("run_${NAME}" COMMAND "${NAME}" VERBATIM)
  add_dependencies("run_${NAME}" "${NAME}")
  add_dependencies(run-examples "run_${NAME}")
endfunction()

add_example(empty_example)

add_folders(Example)
