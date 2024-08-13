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
	include/idiv/continued_fraction/engine/gosper.h
	include/idiv/continued_fraction/engine/log.h
	include/idiv/continued_fraction/engine/rational.h
	include/idiv/continued_fraction/engine/type_erased.h
	include/idiv/continued_fraction/caching_generator.h
	include/idiv/continued_fraction/extended_linear_fractional_mapping.h
	include/idiv/continued_fraction/generator.h
	include/idiv/continued_fraction/projective_rational.h
	include/idiv/best_rational_approx.h
	include/idiv/bigint.h
	include/idiv/frac.h
	include/idiv/fractional_part_extremizer.h
	include/idiv/idiv.h
	include/idiv/interval.h
	include/idiv/optimal_multiply_shift.h
	include/idiv/simultaneous_floor.h
	include/idiv/tmp.h
	include/idiv/util.h
	include/idiv/variant_storage.h
	include/idiv/wuint.h
	include/idiv/xi_zeta_region.h)
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
  option(BUILD_TEST "Build test tree.")
  if(BUILD_EXAMPLES)
    add_subdirectory(subproject/example)
  endif()
  if(BUILD_TEST)
    add_subdirectory(subproject/test)
  endif()
endif()