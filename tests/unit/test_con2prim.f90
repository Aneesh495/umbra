!> Unit tests for Conservative-to-Primitive inversion schemes
module test_con2prim
  use umbra_kinds, only: dp, i4
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  use umbra_metric_base, only: metric_3p1_t
  use umbra_minkowski, only: minkowski_metric_t
  use umbra_eos_ideal, only: ideal_gas_eos_t, make_ideal_gas_eos
  use umbra_con2prim_1d, only: con2prim_1d_solve
  use umbra_con2prim_2d, only: con2prim_2d_solve
  use umbra_error, only: umbra_error_t
  use test_framework, only: test_suite_t, assert_true, assert_close_dp, assert_ok
  implicit none
  private

  public :: run_test_con2prim

contains

  subroutine run_test_con2prim(suite)
    type(test_suite_t), intent(inout) :: suite
    type(prim_state_t) :: p_init, p_rec
    type(cons_state_t) :: u
    type(metric_3p1_t) :: m
    type(minkowski_metric_t) :: mink
    type(ideal_gas_eos_t) :: eos
    type(umbra_error_t) :: err
    real(dp), dimension(0:3) :: x

    eos = make_ideal_gas_eos(4.0_dp / 3.0_dp)
    x = 0.0_dp
    call mink%eval(x, m)

    ! Test 1: 1D Noble-Mignone unmagnetized static fluid inversion round-trip
    p_init%rho = 2.5_dp
    p_init%v = [0.1_dp, 0.0_dp, 0.0_dp]
    p_init%p = 1.8_dp
    p_init%b = 0.0_dp

    call prim_to_cons(p_init, m, eos, u)
    call con2prim_1d_solve(u, m, eos, p_rec, err)

    call assert_ok(suite, err, "1D con2prim solve status OK")
    call assert_close_dp(suite, p_rec%rho, p_init%rho, 1.0e-10_dp, "1D con2prim density recovery")
    call assert_close_dp(suite, p_rec%p, p_init%p, 1.0e-10_dp, "1D con2prim pressure recovery")
    call assert_close_dp(suite, p_rec%v(1), p_init%v(1), 1.0e-10_dp, "1D con2prim vx recovery")

    ! Test 2: 2D Newton-Raphson magnetized relativistic fluid inversion round-trip
    p_init%rho = 1.0_dp
    p_init%v = [0.3_dp, -0.2_dp, 0.1_dp]
    p_init%p = 0.5_dp
    p_init%b = [1.0_dp, 0.5_dp, -0.25_dp]

    call prim_to_cons(p_init, m, eos, u)
    call con2prim_2d_solve(u, m, eos, p_rec, err)

    call assert_ok(suite, err, "2D con2prim solve status OK")
    call assert_close_dp(suite, p_rec%rho, p_init%rho, 1.0e-8_dp, "2D con2prim density recovery")
    call assert_close_dp(suite, p_rec%p, p_init%p, 1.0e-8_dp, "2D con2prim pressure recovery")
    call assert_close_dp(suite, p_rec%v(1), p_init%v(1), 1.0e-8_dp, "2D con2prim vx recovery")
  end subroutine run_test_con2prim

end module test_con2prim
