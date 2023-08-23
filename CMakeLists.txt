cmake_minimum_required(VERSION 3.19)

include(cmake/prelude.cmake)

project(
    idiv
    VERSION 0.1.0
    DESCRIPTION "A library for constant integer division and related tools"
    LANGUAGES NONE
)

docs_early_return()

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(idiv_idiv INTERFACE
	include/idiv/best_rational_approx.h
	include/idiv/big_int.h
	include/idiv/continued_fractions.h
	include/idiv/core.h
	include/idiv/frac.h
	include/idiv/idiv.h
	include/idiv/rational_continued_fractions.h
	include/idiv/wuint.h)
add_library(idiv::idiv ALIAS idiv_idiv)

set_property(
    TARGET idiv_idiv PROPERTY
    EXPORT_NAME idiv
)

target_include_directories(
    idiv_idiv ${warning_guard}
    INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_compile_features(idiv_idiv INTERFACE cxx_std_20)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Examples ----

if(PROJECT_IS_TOP_LEVEL)
  option(BUILD_EXAMPLES "Build examples tree.")
  if(BUILD_EXAMPLES)
    add_subdirectory(example)
  endif()
endif()