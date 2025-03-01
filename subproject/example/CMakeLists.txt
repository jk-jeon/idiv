cmake_minimum_required(VERSION 3.19)

project(idivExamples CXX)

include(../../cmake/project-is-top-level.cmake)
include(../../cmake/folders.cmake)

include(FetchContent)
if (NOT TARGET daw_json_link)
  FetchContent_Declare(
          daw_json_link
          GIT_REPOSITORY https://github.com/beached/daw_json_link.git
          GIT_TAG        v3.30.2
          GIT_SHALLOW    1
  )
  FetchContent_MakeAvailable(daw_json_link)
endif()

if(PROJECT_IS_TOP_LEVEL)
  find_package(idiv REQUIRED)
endif()

add_custom_target(run-examples)

function(add_example NAME)
  set(ADDITIONAL_LIBS ${ARGN})
  add_executable("${NAME}" "${NAME}.cpp")
  target_link_libraries("${NAME}" PRIVATE idiv::idiv "${ADDITIONAL_LIBS}")
  target_compile_features("${NAME}" PRIVATE cxx_std_20)
  add_custom_target("run_${NAME}" COMMAND "${NAME}" VERBATIM)
  add_dependencies("run_${NAME}" "${NAME}")
  add_dependencies(run-examples "run_${NAME}")
endfunction()

add_example(empty_example)
add_example(jeaiii_analysis)
add_example(xi_zeta_region_json daw::daw-json-link)

add_folders(Example)
