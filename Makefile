# umbra: High-Order Discontinuous Galerkin GRMHD and Geodesic Engine
# Modern Fortran Build System

FC ?= gfortran
FFLAGS ?= -O3 -Wall -Wextra -Wpedantic -Werror -fimplicit-none -fPIC -Jinclude -Iinclude
PYTHON ?= /opt/homebrew/bin/python3

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHLIB_EXT = dylib
  SHLIB_FLAGS = -dynamiclib
else
  SHLIB_EXT = so
  SHLIB_FLAGS = -shared
endif

BIN_DIR = bin
LIB_DIR = lib
INC_DIR = include

CORE_SRCS = \
	src/core/umbra_kinds.f90 \
	src/core/umbra_constants.f90 \
	src/core/umbra_error.f90 \
	src/core/umbra_logging.f90 \
	src/core/umbra_timer.f90 \
	src/core/umbra_math.f90 \
	src/core/umbra_polynomials.f90 \
	src/core/umbra_quadrature.f90 \
	src/spacetime/umbra_metric_base.f90 \
	src/spacetime/umbra_minkowski.f90 \
	src/spacetime/umbra_schwarzschild.f90 \
	src/spacetime/umbra_kerr_bl.f90 \
	src/spacetime/umbra_kerr_ks.f90 \
	src/spacetime/umbra_connection.f90 \
	src/spacetime/umbra_tetrad.f90 \
	src/eos/umbra_eos_base.f90 \
	src/eos/umbra_eos_ideal.f90 \
	src/eos/umbra_eos_polytropic.f90 \
	src/eos/umbra_eos_synge.f90 \
	src/mhd/umbra_state.f90 \
	src/mhd/umbra_con2prim.f90 \
	src/mhd/umbra_con2prim_1d.f90 \
	src/mhd/umbra_con2prim_2d.f90 \
	src/mhd/umbra_wavespeeds.f90 \
	src/mhd/umbra_fluxes.f90 \
	src/mhd/umbra_sources.f90 \
	src/mhd/umbra_riemann.f90 \
	src/mhd/umbra_riemann_llf.f90 \
	src/mhd/umbra_riemann_hll.f90 \
	src/mhd/umbra_riemann_hllc.f90 \
	src/mhd/umbra_riemann_hlld.f90 \
	src/divergence/umbra_ct.f90 \
	src/divergence/umbra_glm.f90 \
	src/mesh/umbra_grid1d.f90 \
	src/mesh/umbra_grid2d.f90 \
	src/mesh/umbra_grid3d.f90 \
	src/mesh/umbra_boundary.f90 \
	src/mesh/umbra_octree.f90 \
	src/reconstruction/umbra_limiters.f90 \
	src/reconstruction/umbra_reconstruct_tvd.f90 \
	src/reconstruction/umbra_reconstruct_ppm.f90 \
	src/reconstruction/umbra_reconstruct_weno.f90 \
	src/dg/umbra_dg_basis.f90 \
	src/dg/umbra_dg_element.f90 \
	src/dg/umbra_dg_volume.f90 \
	src/dg/umbra_dg_surface.f90 \
	src/dg/umbra_dg_solver.f90 \
	src/timestepping/umbra_cfl.f90 \
	src/timestepping/umbra_rk_ssprk.f90 \
	src/timestepping/umbra_rk_dopri5.f90 \
	src/raytracing/umbra_carter_constant.f90 \
	src/raytracing/umbra_geodesic_integrator.f90 \
	src/raytracing/umbra_synchrotron.f90 \
	src/raytracing/umbra_polarimetry.f90 \
	src/raytracing/umbra_radiative_transfer.f90 \
	src/raytracing/umbra_camera.f90 \
	src/benchmarks/umbra_bench_sod.f90 \
	src/benchmarks/umbra_bench_briowu.f90 \
	src/benchmarks/umbra_bench_blast.f90 \
	src/benchmarks/umbra_bench_kelvin_helmholtz.f90 \
	src/benchmarks/umbra_bench_alfven_wave.f90 \
	src/benchmarks/umbra_bench_orszag_tang.f90 \
	src/benchmarks/umbra_bench_fishbone_moncrief.f90 \
	src/benchmarks/umbra_bench_kerr_shadow.f90 \
	src/io/umbra_io_ascii.f90 \
	src/io/umbra_io_raw.f90 \
	src/io/umbra_io_vtk.f90 \
	src/io/umbra_checkpoint.f90 \
	src/c_api/umbra_c.f90 \
	src/umbra_mod.f90

TARGET_CLI = $(BIN_DIR)/umbra
TARGET_LIB = $(LIB_DIR)/libumbra.$(SHLIB_EXT)

UNIT_TEST_SRCS = \
	tests/unit/test_framework.f90 \
	tests/unit/test_math_quadrature.f90 \
	tests/unit/test_spacetime_metric.f90 \
	tests/unit/test_eos.f90 \
	tests/unit/test_con2prim.f90 \
	tests/unit/test_riemann_solvers.f90 \
	tests/unit/test_reconstruction.f90 \
	tests/unit/test_dg_operators.f90 \
	tests/unit/test_geodesics.f90 \
	tests/unit/test_unit_runner.f90

ACCEPTANCE_TEST_SRCS = tests/acceptance/test_acceptance_runner.f90
NEGATIVE_TEST_SRCS = tests/negative/test_negative_controls.f90

TARGET_UNIT = $(BIN_DIR)/test_unit_runner
TARGET_ACCEPTANCE = $(BIN_DIR)/test_acceptance_runner
TARGET_NEGATIVE = $(BIN_DIR)/test_negative_controls

.PHONY: all clean test unit acceptance negative verify benchmark

all: dirs $(TARGET_LIB) $(TARGET_CLI)

dirs:
	@mkdir -p $(BIN_DIR) $(LIB_DIR) $(INC_DIR)

$(TARGET_LIB): $(CORE_SRCS)
	@mkdir -p $(LIB_DIR) $(INC_DIR)
	$(FC) $(FFLAGS) $(SHLIB_FLAGS) $(CORE_SRCS) -o $(TARGET_LIB)
	@if [ "$(UNAME_S)" = "Darwin" ]; then \
		cp $(TARGET_LIB) $(LIB_DIR)/libumbra.so 2>/dev/null || true; \
	fi

$(TARGET_CLI): $(CORE_SRCS) src/driver/umbra_main.f90
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) $(CORE_SRCS) src/driver/umbra_main.f90 -o $(TARGET_CLI)

$(TARGET_UNIT): $(CORE_SRCS) $(UNIT_TEST_SRCS)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) $(CORE_SRCS) $(UNIT_TEST_SRCS) -o $(TARGET_UNIT)

$(TARGET_ACCEPTANCE): $(CORE_SRCS) $(ACCEPTANCE_TEST_SRCS)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) $(CORE_SRCS) $(ACCEPTANCE_TEST_SRCS) -o $(TARGET_ACCEPTANCE)

$(TARGET_NEGATIVE): $(CORE_SRCS) $(NEGATIVE_TEST_SRCS)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) $(CORE_SRCS) $(NEGATIVE_TEST_SRCS) -o $(TARGET_NEGATIVE)

unit: $(TARGET_UNIT)
	@echo "=== Running unit test suite (70 assertions) ==="
	@./$(TARGET_UNIT)

acceptance: $(TARGET_ACCEPTANCE)
	@echo "=== Running acceptance integration test suite (5 benchmarks) ==="
	@./$(TARGET_ACCEPTANCE)

negative: $(TARGET_NEGATIVE)
	@echo "=== Running negative defensive test suite (6 traps) ==="
	@./$(TARGET_NEGATIVE)

test: unit acceptance negative verify

verify:
	@echo "=== Running 12 adversarial negative controls verification ==="
	@$(PYTHON) scripts/verify_evidence.py

benchmark:
	@echo "=== Executing benchmark evidence generator ==="
	@$(PYTHON) scripts/run_acceptance.py

clean:
	@rm -rf $(BIN_DIR)/* $(LIB_DIR)/*
	@rm -f $(INC_DIR)/*.mod
