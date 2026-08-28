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
  src/core/units.f90
  src/core/hall_induction.f90
  src/core/io_checkpoint.f90
  src/core/ode_integrator.f90
  src/core/eos_table.f90
  src/core/tov_solver.f90
  src/core/crust_conductivity.f90
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
set(MHDVSH_HALL_SOURCES
  src/regimes/diffusion/diffusion_regime.f90
  src/regimes/hall/hall_regime.f90)

# ── App ──────────────────────────────────────────────────────────────────────
add_executable(mhdvsh_run app/mhdvsh_run.f90)
set_target_properties(mhdvsh_run PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_run"
)
target_link_libraries(mhdvsh_run PRIVATE mhdvsh_core)

add_executable(mhdvsh_diffusion ${MHDVSH_DIFFUSION_SOURCES} app/mhdvsh_diffusion.f90)
set_target_properties(mhdvsh_diffusion PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_diffusion"
)
target_link_libraries(mhdvsh_diffusion PRIVATE mhdvsh_core)

# TOV stellar-structure driver (see src/core/tov_solver.f90,
# eos_table.f90, crust_conductivity.f90 -- ported from ~/Desktop/EOSNS/).
# Reads the EOS table from data/eos/ by default, path relative to the
# working directory the driver is run from (matching how mhdvsh_diffusion
# et al. take their data-output path as a plain CLI argument, not an
# install-time-resolved path).
add_executable(mhdvsh_tov app/mhdvsh_tov.f90)
set_target_properties(mhdvsh_tov PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_tov"
)
target_link_libraries(mhdvsh_tov PRIVATE mhdvsh_core)

# Combined resistive+Hall regime -- the "primary interest" case (see
# hall_regime.f90's own header for the splitting/BC design). Unlike the
# Stage 2 stability harness, this is a real regime now: in ALL/ctest,
# mirroring diffusion's wiring exactly.
add_executable(mhdvsh_hall ${MHDVSH_HALL_SOURCES} app/mhdvsh_hall.f90)
set_target_properties(mhdvsh_hall PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_hall"
)
target_link_libraries(mhdvsh_hall PRIVATE mhdvsh_core)

# Same regime, dynamic Hall-CFL-limited dt instead of mhdvsh_hall's
# fixed one (see hall_regime.f90's HALL_COMPUTE_DT/HALL_SET_DT and
# timestepper.f90's RUN_ADAPTIVE).
add_executable(mhdvsh_hall_adaptive ${MHDVSH_HALL_SOURCES} app/mhdvsh_hall_adaptive.f90)
set_target_properties(mhdvsh_hall_adaptive PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_hall_adaptive"
)
target_link_libraries(mhdvsh_hall_adaptive PRIVATE mhdvsh_core)

# Stage 2 Hall stability experiment (see .claude/plans/i-would-like-to-
# quizzical-comet.md) -- deliberately NOT added to the default ALL
# target or wired to any artifact/plot command: this is an exploratory
# harness whose parameters (dt, Lmax, F_Hall, ...) are expected to keep
# changing between runs, unlike diffusion's fixed regression-style demo.
# Build/run explicitly: `cmake --build . --target mhdvsh_hall_stability_experiment`.
add_executable(mhdvsh_hall_stability_experiment EXCLUDE_FROM_ALL
  app/mhdvsh_hall_stability_experiment.f90)
set_target_properties(mhdvsh_hall_stability_experiment PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_hall_stability_experiment"
)
target_link_libraries(mhdvsh_hall_stability_experiment PRIVATE mhdvsh_core)

# Combined-regime dt-stability sweep (see mhdvsh_hall_dt_sweep.f90's own
# header) -- same exploratory-harness treatment as
# mhdvsh_hall_stability_experiment above: not in ALL, parameters (dt,
# tmax) passed on the command line per invocation, not baked in.
add_executable(mhdvsh_hall_dt_sweep EXCLUDE_FROM_ALL
  ${MHDVSH_HALL_SOURCES} app/mhdvsh_hall_dt_sweep.f90)
set_target_properties(mhdvsh_hall_dt_sweep PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_hall_dt_sweep"
)
target_link_libraries(mhdvsh_hall_dt_sweep PRIVATE mhdvsh_core)

# Same regime/grid/timestepper as mhdvsh_hall, but with a physically-
# sourced (eta, f_H) pair read off the real crust profile instead of
# mhdvsh_hall's own toy constants -- see mhdvsh_hall_crust_sample.f90's
# own header for full input provenance. A deliberate, separate
# comparison run (user, 2026-08-24), not a modification of mhdvsh_hall
# itself -- same EXCLUDE_FROM_ALL treatment as the other one-off
# comparison drivers above (not part of the default ALL/ctest surface).
add_executable(mhdvsh_hall_crust_sample EXCLUDE_FROM_ALL
  ${MHDVSH_HALL_SOURCES} app/mhdvsh_hall_crust_sample.f90)
set_target_properties(mhdvsh_hall_crust_sample PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_hall_crust_sample"
)
target_link_libraries(mhdvsh_hall_crust_sample PRIVATE mhdvsh_core)

# Same regime/grid as mhdvsh_hall_crust_sample, but with a full,
# radially-varying eta(r)/f_H(r) (CRUST_CONDUCTIVITY::
# ETA_AND_F_HALL_ON_GRID, from a live TOV solve) instead of one
# physically-sourced constant pair -- see
# mhdvsh_hall_crust_profile.f90's own header. dt/n_steps/n_sub are CLI
# arguments (like mhdvsh_hall_dt_sweep) for the staged smoke-test/
# dt-check/full-1000yr-run workflow. Same EXCLUDE_FROM_ALL treatment as
# the other one-off comparison drivers.
add_executable(mhdvsh_hall_crust_profile EXCLUDE_FROM_ALL
  ${MHDVSH_HALL_SOURCES} app/mhdvsh_hall_crust_profile.f90)
set_target_properties(mhdvsh_hall_crust_profile PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/mhdvsh_hall_crust_profile"
)
target_link_libraries(mhdvsh_hall_crust_profile PRIVATE mhdvsh_core)

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

  set(MHDVSH_HALL_ENERGY_DATA
      ${MHDVSH_ARTIFACTS_DIR}/hall/energy_budget/data/energy_budget.dat)
  set(MHDVSH_HALL_ENERGY_PLOT
      ${MHDVSH_ARTIFACTS_DIR}/hall/energy_budget/plots/energy_budget.png)

  add_custom_command(
    OUTPUT ${MHDVSH_HALL_ENERGY_DATA}
    COMMAND ${CMAKE_COMMAND} -E make_directory
            ${MHDVSH_ARTIFACTS_DIR}/hall/energy_budget/data
    COMMAND ${CMAKE_COMMAND} -E make_directory
            ${MHDVSH_ARTIFACTS_DIR}/hall/energy_budget/plots
    COMMAND mhdvsh_hall ${MHDVSH_HALL_ENERGY_DATA}
    DEPENDS mhdvsh_hall
    COMMENT "Running combined resistive+Hall regime to produce energy-budget time series"
    VERBATIM
  )

  add_custom_command(
    OUTPUT ${MHDVSH_HALL_ENERGY_PLOT}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/plot_energy_budget.py
            ${MHDVSH_HALL_ENERGY_DATA} ${MHDVSH_HALL_ENERGY_PLOT}
    DEPENDS ${MHDVSH_HALL_ENERGY_DATA} ${CMAKE_SOURCE_DIR}/scripts/plot_energy_budget.py
    COMMENT "Plotting combined resistive+Hall regime energy budget"
    VERBATIM
  )

  # NOT in ALL, unlike diffusion's equivalent target: this driver's
  # LMAX=30/N_SUB=10 configuration costs ~5 minutes to run (the price of
  # a real Hall cascade at this resolution), vs. diffusion's sub-second
  # implicit solve -- forcing that into every `cmake --build .` would be
  # a surprising, unwelcome cost for routine development. Build/run
  # explicitly: `cmake --build . --target hall_energy_budget_plot`.
  add_custom_target(hall_energy_budget_plot
    DEPENDS ${MHDVSH_HALL_ENERGY_PLOT}
  )
else()
  message(STATUS "Python3 interpreter not found -- skipping artifact plot targets")
endif()

# ── Tests ────────────────────────────────────────────────────────────────────
add_executable(test_transforms test/test_transforms.f90)
set_target_properties(test_transforms PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_transforms"
)
target_link_libraries(test_transforms PRIVATE mhdvsh_core)
add_test(NAME transforms COMMAND test_transforms)

add_executable(test_boundary_conditions test/test_boundary_conditions.f90)
set_target_properties(test_boundary_conditions PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_boundary_conditions"
)
target_link_libraries(test_boundary_conditions PRIVATE mhdvsh_core)
add_test(NAME boundary_conditions COMMAND test_boundary_conditions)

add_executable(test_linear_solve test/test_linear_solve.f90)
set_target_properties(test_linear_solve PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_linear_solve"
)
target_link_libraries(test_linear_solve PRIVATE mhdvsh_core)
add_test(NAME linear_solve COMMAND test_linear_solve)

add_executable(test_radial_operators test/test_radial_operators.f90)
set_target_properties(test_radial_operators PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_radial_operators"
)
target_link_libraries(test_radial_operators PRIVATE mhdvsh_core)
add_test(NAME radial_operators COMMAND test_radial_operators)

add_executable(test_timestepper test/test_timestepper.f90)
set_target_properties(test_timestepper PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_timestepper"
)
target_link_libraries(test_timestepper PRIVATE mhdvsh_core)
add_test(NAME timestepper COMMAND test_timestepper)

add_executable(test_field_diagnostics test/test_field_diagnostics.f90)
set_target_properties(test_field_diagnostics PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_field_diagnostics"
)
target_link_libraries(test_field_diagnostics PRIVATE mhdvsh_core)
add_test(NAME field_diagnostics COMMAND test_field_diagnostics)

add_executable(test_units test/test_units.f90)
set_target_properties(test_units PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_units"
)
target_link_libraries(test_units PRIVATE mhdvsh_core)
add_test(NAME units COMMAND test_units)

add_executable(test_hall_induction test/test_hall_induction.f90)
set_target_properties(test_hall_induction PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_hall_induction"
)
target_link_libraries(test_hall_induction PRIVATE mhdvsh_core)
add_test(NAME hall_induction COMMAND test_hall_induction)

add_executable(test_io_checkpoint test/test_io_checkpoint.f90)
set_target_properties(test_io_checkpoint PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_io_checkpoint"
)
target_link_libraries(test_io_checkpoint PRIVATE mhdvsh_core)
add_test(NAME io_checkpoint COMMAND test_io_checkpoint)

add_executable(test_diffusion_regime ${MHDVSH_DIFFUSION_SOURCES} test/test_diffusion_regime.f90)
set_target_properties(test_diffusion_regime PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_diffusion_regime"
)
target_link_libraries(test_diffusion_regime PRIVATE mhdvsh_core)
add_test(NAME diffusion_regime COMMAND test_diffusion_regime)

add_executable(test_hall_regime ${MHDVSH_HALL_SOURCES} test/test_hall_regime.f90)
set_target_properties(test_hall_regime PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_hall_regime"
)
target_link_libraries(test_hall_regime PRIVATE mhdvsh_core)
add_test(NAME hall_regime COMMAND test_hall_regime)

# TOV/EOS port tests (see src/core/eos_table.f90, tov_solver.f90,
# crust_conductivity.f90) -- test_eos_table and test_tov_solver read
# data/eos/lowd-eos.ja.tab via a path relative to the source root, so
# their WORKING_DIRECTORY is set explicitly (unlike every other test
# here, none of which touch the filesystem).
add_executable(test_eos_table test/test_eos_table.f90)
set_target_properties(test_eos_table PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_eos_table"
)
target_link_libraries(test_eos_table PRIVATE mhdvsh_core)
add_test(NAME eos_table COMMAND test_eos_table WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_tov_solver test/test_tov_solver.f90)
set_target_properties(test_tov_solver PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_tov_solver"
)
target_link_libraries(test_tov_solver PRIVATE mhdvsh_core)
add_test(NAME tov_solver COMMAND test_tov_solver WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

# test_crust_conductivity's reference values were produced independently
# (see its own header) from the original, unmodified potekhinc.f -- it
# reads no files itself, no WORKING_DIRECTORY needed.
add_executable(test_crust_conductivity test/test_crust_conductivity.f90)
set_target_properties(test_crust_conductivity PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_crust_conductivity"
)
target_link_libraries(test_crust_conductivity PRIVATE mhdvsh_core)
add_test(NAME crust_conductivity COMMAND test_crust_conductivity)
