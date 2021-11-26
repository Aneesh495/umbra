cmake_minimum_required(VERSION 3.20)
project(umbra
    VERSION 1.0.0
    DESCRIPTION "High-Order Discontinuous Galerkin GRMHD and Geodesic Engine"
    LANGUAGES Fortran C
)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    add_compile_options(-O3 -Wall -Wextra -Wpedantic -Werror -fimplicit-none -fPIC)
endif()

set(CORE_SOURCES
    src/core/umbra_kinds.f90
    src/core/umbra_constants.f90
    src/core/umbra_error.f90
    src/core/umbra_logging.f90
    src/core/umbra_timer.f90
    src/core/umbra_math.f90
    src/core/umbra_polynomials.f90
    src/core/umbra_quadrature.f90
    src/spacetime/umbra_metric_base.f90
    src/spacetime/umbra_minkowski.f90
    src/spacetime/umbra_schwarzschild.f90
    src/spacetime/umbra_kerr_bl.f90
    src/spacetime/umbra_kerr_ks.f90
    src/spacetime/umbra_connection.f90
    src/spacetime/umbra_tetrad.f90
    src/eos/umbra_eos_base.f90
    src/eos/umbra_eos_ideal.f90
    src/eos/umbra_eos_polytropic.f90
    src/eos/umbra_eos_synge.f90
    src/mhd/umbra_state.f90
    src/mhd/umbra_con2prim.f90
    src/mhd/umbra_con2prim_1d.f90
    src/mhd/umbra_con2prim_2d.f90
    src/mhd/umbra_wavespeeds.f90
    src/mhd/umbra_fluxes.f90
    src/mhd/umbra_sources.f90
    src/mhd/umbra_riemann.f90
    src/mhd/umbra_riemann_llf.f90
    src/mhd/umbra_riemann_hll.f90
    src/mhd/umbra_riemann_hllc.f90
    src/mhd/umbra_riemann_hlld.f90
    src/divergence/umbra_ct.f90
    src/divergence/umbra_glm.f90
    src/mesh/umbra_grid1d.f90
    src/mesh/umbra_grid2d.f90
    src/mesh/umbra_grid3d.f90
    src/mesh/umbra_boundary.f90
    src/mesh/umbra_octree.f90
    src/reconstruction/umbra_limiters.f90
    src/reconstruction/umbra_reconstruct_tvd.f90
    src/reconstruction/umbra_reconstruct_ppm.f90
    src/reconstruction/umbra_reconstruct_weno.f90
    src/dg/umbra_dg_basis.f90
    src/dg/umbra_dg_element.f90
    src/dg/umbra_dg_volume.f90
    src/dg/umbra_dg_surface.f90
    src/dg/umbra_dg_solver.f90
    src/timestepping/umbra_cfl.f90
    src/timestepping/umbra_rk_ssprk.f90
    src/timestepping/umbra_rk_dopri5.f90
    src/raytracing/umbra_carter_constant.f90
    src/raytracing/umbra_geodesic_integrator.f90
    src/raytracing/umbra_synchrotron.f90
    src/raytracing/umbra_polarimetry.f90
    src/raytracing/umbra_radiative_transfer.f90
    src/raytracing/umbra_camera.f90
    src/benchmarks/umbra_bench_sod.f90
    src/benchmarks/umbra_bench_briowu.f90
    src/benchmarks/umbra_bench_blast.f90
    src/benchmarks/umbra_bench_kelvin_helmholtz.f90
    src/benchmarks/umbra_bench_alfven_wave.f90
    src/benchmarks/umbra_bench_orszag_tang.f90
    src/benchmarks/umbra_bench_fishbone_moncrief.f90
    src/benchmarks/umbra_bench_kerr_shadow.f90
    src/io/umbra_io_ascii.f90
    src/io/umbra_io_raw.f90
    src/io/umbra_io_vtk.f90
    src/io/umbra_checkpoint.f90
    src/c_api/umbra_c.f90
    src/umbra_mod.f90
)

# Shared Library
add_library(umbra_shared SHARED ${CORE_SOURCES})
set_target_properties(umbra_shared PROPERTIES
    OUTPUT_NAME "umbra"
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib
)

# Executable Driver
add_executable(umbra_bin ${CORE_SOURCES} src/driver/umbra_main.f90)
set_target_properties(umbra_bin PROPERTIES
    OUTPUT_NAME "umbra"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
)

# Test Executables
enable_testing()

set(UNIT_TEST_SOURCES
    tests/unit/test_framework.f90
    tests/unit/test_math_quadrature.f90
    tests/unit/test_spacetime_metric.f90
    tests/unit/test_eos.f90
    tests/unit/test_con2prim.f90
    tests/unit/test_riemann_solvers.f90
    tests/unit/test_reconstruction.f90
    tests/unit/test_dg_operators.f90
    tests/unit/test_geodesics.f90
    tests/unit/test_unit_runner.f90
)

add_executable(test_unit_runner ${CORE_SOURCES} ${UNIT_TEST_SOURCES})
set_target_properties(test_unit_runner PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
)
add_test(NAME UnitTests COMMAND test_unit_runner)

add_executable(test_acceptance_runner ${CORE_SOURCES} tests/acceptance/test_acceptance_runner.f90)
set_target_properties(test_acceptance_runner PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
)
add_test(NAME AcceptanceTests COMMAND test_acceptance_runner)

add_executable(test_negative_controls ${CORE_SOURCES} tests/negative/test_negative_controls.f90)
set_target_properties(test_negative_controls PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
)
add_test(NAME NegativeControls COMMAND test_negative_controls)
