cmake_minimum_required(VERSION 3.14)
project(MHDVSH
  VERSION   0.1.0
  LANGUAGES Fortran
)

include(GNUInstallDirs)
include(CTest)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(MHDVSH_FETCH_FORTVSH "Build FORTVSH from a local source checkout instead of find_package" OFF)
set(MHDVSH_FORTVSH_SOURCE_DIR "" CACHE PATH "Path to a FORTVSH source checkout (used when MHDVSH_FETCH_FORTVSH=ON)")
option(MHDVSH_ENABLE_OPENMP "Link against OpenMP for the radial-index-parallel transform/solve loops" ON)

# ── FORTVSH dependency ────────────────────────────────────────────────────────
if(MHDVSH_FETCH_FORTVSH)
  if(NOT MHDVSH_FORTVSH_SOURCE_DIR)
    message(FATAL_ERROR "MHDVSH_FETCH_FORTVSH=ON requires -DMHDVSH_FORTVSH_SOURCE_DIR=<path to fortran_VSH>")
  endif()
  include(FetchContent)
  FetchContent_Declare(fortvsh SOURCE_DIR "${MHDVSH_FORTVSH_SOURCE_DIR}")
  FetchContent_GetProperties(fortvsh)
  if(NOT fortvsh_POPULATED)
    FetchContent_Populate(fortvsh)
    add_subdirectory(${fortvsh_SOURCE_DIR} ${fortvsh_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
else()
  find_package(FORTVSH REQUIRED)
endif()

if(MHDVSH_ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
endif()

find_package(LAPACK REQUIRED)

add_compile_options(-O3 -Wall)

# ── Core library ───────────────────────────────────────────────────────────────
# NOTE: this file list currently covers Phase 0 (grid/field skeleton), Phase 1
# (transform_plan.f90, transforms.f90), the Phase 2 radial FD modules
# (fd_weights.f90, radial_operators.f90, implemented ahead of schedule since
# they have no dependency on the transform layer), the physics content of
# boundary_conditions.f90 (Phase 3's vacuum outer BC plus an inner Dirichlet
# BC), linear_solve.f90 (a generic dense per-system LAPACK wrapper with no
# knowledge of spherical harmonics or physics -- every regime assembles its
# own system matrix from these building blocks), regime_interface.f90/
# timestepper.f90 (the opaque per-regime ADVANCE contract and the shared
# fixed-step loop that drives it), and field_diagnostics.f90 (energy computed
# directly from Phi/Psi -- see its header for which parts are provisional,
# pending re-derivation of the underlying physics coefficients). The first
# real regime, src/regimes/diffusion/diffusion_regime.f90 (pure linear
# Ohmic diffusion), is wired up separately below -- see "Regimes". io_*.f90
# remains, added once there's a regime worth writing output for.
add_library(mhdvsh_core STATIC
  src/core/grid_angular.f90
  src/core/grid_radial.f90
  src/core/field_types.f90
  src/core/fd_weights.f90
  src/core/radial_operators.f90
  src/core/transform_plan.f90
  src/core/transforms.f90
  src/core/boundary_conditions.f90
  src/core/linear_solve.f90
  src/core/regime_interface.f90
  src/core/timestepper.f90
  src/core/field_diagnostics.f90
)
set_target_properties(mhdvsh_core PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_include_directories(mhdvsh_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
)
target_link_libraries(mhdvsh_core PUBLIC FORTVSH::vsh ${LAPACK_LIBRARIES})
if(MHDVSH_ENABLE_OPENMP AND OpenMP_Fortran_FOUND)
  target_link_libraries(mhdvsh_core PUBLIC OpenMP::OpenMP_Fortran)
endif()

# ── Regimes ──────────────────────────────────────────────────────────────────
# Regime sources are NOT part of mhdvsh_core (which stays regime-agnostic) --
# each regime compiles its own sources directly into its own driver
# executable(s) below, per-regime-executable, no runtime dispatch (see
# REGIME_INTERFACE/TIMESTEPPER).
set(MHDVSH_DIFFUSION_SOURCES src/regimes/diffusion/diffusion_regime.f90)

# ── App ──────────────────────────────────────────────────────────────────────
add_executable(mhdvsh_run app/mhdvsh_run.f90)
set_target_properties(mhdvsh_run PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(mhdvsh_run PRIVATE mhdvsh_core)

add_executable(mhdvsh_diffusion ${MHDVSH_DIFFUSION_SOURCES} app/mhdvsh_diffusion.f90)
set_target_properties(mhdvsh_diffusion PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(mhdvsh_diffusion PRIVATE mhdvsh_core)

# ── Artifacts (plots/data regenerated automatically with every build) ───────
# Directory convention: artifacts/<regime>/<diagnostic>/{data,plots}/ -- one
# data-producing driver + one scripts/plot_*.py + one add_custom_target(...
# ALL) per diagnostic, added here as each is built. Today: diffusion's
# energy_budget (total magnetic energy only -- see FIELD_DIAGNOSTICS and
# app/mhdvsh_diffusion.f90 for why the Joule/Poynting columns aren't there
# yet). CMake's OUTPUT-based custom-command dependency chain means each step
# only reruns when its inputs actually changed, same as any other build rule.
find_package(Python3 COMPONENTS Interpreter)

if(Python3_Interpreter_FOUND)
  set(MHDVSH_ARTIFACTS_DIR ${CMAKE_BINARY_DIR}/artifacts)
  set(MHDVSH_DIFFUSION_ENERGY_DATA
      ${MHDVSH_ARTIFACTS_DIR}/diffusion/energy_budget/data/energy_budget.dat)
  set(MHDVSH_DIFFUSION_ENERGY_PLOT
      ${MHDVSH_ARTIFACTS_DIR}/diffusion/energy_budget/plots/energy_budget.png)

  add_custom_command(
    OUTPUT ${MHDVSH_DIFFUSION_ENERGY_DATA}
    COMMAND ${CMAKE_COMMAND} -E make_directory
            ${MHDVSH_ARTIFACTS_DIR}/diffusion/energy_budget/data
    COMMAND ${CMAKE_COMMAND} -E make_directory
            ${MHDVSH_ARTIFACTS_DIR}/diffusion/energy_budget/plots
    COMMAND mhdvsh_diffusion ${MHDVSH_DIFFUSION_ENERGY_DATA}
    DEPENDS mhdvsh_diffusion
    COMMENT "Running diffusion regime to produce energy-budget time series"
    VERBATIM
  )

  add_custom_command(
    OUTPUT ${MHDVSH_DIFFUSION_ENERGY_PLOT}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/plot_energy_budget.py
            ${MHDVSH_DIFFUSION_ENERGY_DATA} ${MHDVSH_DIFFUSION_ENERGY_PLOT}
    DEPENDS ${MHDVSH_DIFFUSION_ENERGY_DATA} ${CMAKE_SOURCE_DIR}/scripts/plot_energy_budget.py
    COMMENT "Plotting diffusion regime energy budget"
    VERBATIM
  )

  add_custom_target(diffusion_energy_budget_plot ALL
    DEPENDS ${MHDVSH_DIFFUSION_ENERGY_PLOT}
  )
else()
  message(STATUS "Python3 interpreter not found -- skipping artifact plot targets")
endif()

# ── Tests ────────────────────────────────────────────────────────────────────
add_executable(test_transforms test/test_transforms.f90)
set_target_properties(test_transforms PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(test_transforms PRIVATE mhdvsh_core)
add_test(NAME transforms COMMAND test_transforms)

add_executable(test_boundary_conditions test/test_boundary_conditions.f90)
set_target_properties(test_boundary_conditions PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(test_boundary_conditions PRIVATE mhdvsh_core)
add_test(NAME boundary_conditions COMMAND test_boundary_conditions)

add_executable(test_linear_solve test/test_linear_solve.f90)
set_target_properties(test_linear_solve PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(test_linear_solve PRIVATE mhdvsh_core)
add_test(NAME linear_solve COMMAND test_linear_solve)

add_executable(test_radial_operators test/test_radial_operators.f90)
set_target_properties(test_radial_operators PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(test_radial_operators PRIVATE mhdvsh_core)
add_test(NAME radial_operators COMMAND test_radial_operators)

add_executable(test_timestepper test/test_timestepper.f90)
set_target_properties(test_timestepper PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(test_timestepper PRIVATE mhdvsh_core)
add_test(NAME timestepper COMMAND test_timestepper)

add_executable(test_field_diagnostics test/test_field_diagnostics.f90)
set_target_properties(test_field_diagnostics PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(test_field_diagnostics PRIVATE mhdvsh_core)
add_test(NAME field_diagnostics COMMAND test_field_diagnostics)

add_executable(test_diffusion_regime ${MHDVSH_DIFFUSION_SOURCES} test/test_diffusion_regime.f90)
set_target_properties(test_diffusion_regime PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)
target_link_libraries(test_diffusion_regime PRIVATE mhdvsh_core)
add_test(NAME diffusion_regime COMMAND test_diffusion_regime)
