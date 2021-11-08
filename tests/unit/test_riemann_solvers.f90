!> Unit tests for relativistic Riemann solvers (LLF, HLL, HLLD)
module test_riemann_solvers
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS
  use umbra_state, only: prim_state_t
  use umbra_metric_base, only: metric_3p1_t
  use umbra_minkowski, only: minkowski_metric_t
  use umbra_eos_ideal, only: ideal_gas_eos_t, make_ideal_gas_eos
  use umbra_riemann_llf, only: riemann_solve_llf
  use umbra_riemann_hll, only: riemann_solve_hll
  use umbra_riemann_hlld, only: riemann_solve_hlld
  use umbra_error, only: umbra_error_t
  use test_framework, only: test_suite_t, assert_true, assert_close_dp
  implicit none
  private

  public :: run_test_riemann_solvers

contains

  subroutine run_test_riemann_solvers(suite)
    type(test_suite_t), intent(inout) :: suite
    type(prim_state_t) :: pl, pr
    type(metric_3p1_t) :: m
    type(minkowski_metric_t) :: mink
    type(ideal_gas_eos_t) :: eos
    type(umbra_error_t) :: err
    real(dp), dimension(NCONS) :: f_llf, f_hll, f_hlld
    real(dp) :: cmax_llf, cmax_hll, cmax_hlld
    real(dp), dimension(0:3) :: x

    eos = make_ideal_gas_eos(1.4_dp)
    x = 0.0_dp
    call mink%eval(x, m)

    ! Setup isolated contact discontinuity along x (dir=1): equal pressure & velocity, different density
    pl%rho = 2.0_dp
    pl%v = [0.2_dp, 0.0_dp, 0.0_dp]
    pl%p = 1.0_dp
    pl%b = [0.5_dp, 0.0_dp, 0.0_dp]

    pr%rho = 1.0_dp
    pr%v = [0.2_dp, 0.0_dp, 0.0_dp]
    pr%p = 1.0_dp
    pr%b = [0.5_dp, 0.0_dp, 0.0_dp]

    ! Evaluate all three Riemann solvers
    call riemann_solve_llf(pl, pr, m, eos, 1_i4, f_llf, cmax_llf, err)
    call riemann_solve_hll(pl, pr, m, eos, 1_i4, f_hll, cmax_hll, err)
    call riemann_solve_hlld(pl, pr, m, eos, 1_i4, f_hlld, cmax_hlld, err)

    ! Test 1: Positive maximum wave speed across all solvers
    call assert_true(suite, cmax_llf > 0.0_dp, "LLF wave speed > 0")
    call assert_true(suite, cmax_hll > 0.0_dp, "HLL wave speed > 0")
    call assert_true(suite, cmax_hlld > 0.0_dp, "HLLD wave speed > 0")

    ! Test 2: Subluminal wave speeds (< 1 in geometric units)
    call assert_true(suite, cmax_llf < 1.0_dp, "LLF wave speed < 1")
    call assert_true(suite, cmax_hll < 1.0_dp, "HLL wave speed < 1")
    call assert_true(suite, cmax_hlld < 1.0_dp, "HLLD wave speed < 1")

    ! Test 3: Normal magnetic field flux along direction of Bx is identically 0
    call assert_close_dp(suite, f_llf(6), 0.0_dp, 1.0e-14_dp, "LLF normal B flux = 0")
    call assert_close_dp(suite, f_hll(6), 0.0_dp, 1.0e-14_dp, "HLL normal B flux = 0")
    call assert_close_dp(suite, f_hlld(6), 0.0_dp, 1.0e-14_dp, "HLLD normal B flux = 0")
  end subroutine run_test_riemann_solvers

end module test_riemann_solvers
