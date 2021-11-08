!> Unit tests for Kerr null geodesic integration and Carter constant invariant
module test_geodesics
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: PI
  use umbra_kerr_bl, only: kerr_bl_metric_t
  use umbra_geodesic_integrator, only: ray_state_t, integrate_geodesic_step_rk4, geodesic_hamiltonian
  use umbra_carter_constant, only: compute_carter_constant
  use umbra_error, only: umbra_error_t
  use test_framework, only: test_suite_t, assert_true, assert_close_dp, assert_ok
  implicit none
  private

  public :: run_test_geodesics

contains

  subroutine run_test_geodesics(suite)
    type(test_suite_t), intent(inout) :: suite
    type(kerr_bl_metric_t) :: metric
    type(ray_state_t) :: ray
    type(umbra_error_t) :: err
    real(dp) :: h_val, q_init, q_now, d_lambda, r_horiz, r_escape
    integer(i4) :: step

    metric%mass = 1.0_dp
    metric%spin = 0.9_dp
    r_horiz = metric%horizon_outer() + 0.01_dp
    r_escape = 200.0_dp
    d_lambda = 0.05_dp

    ! Initial ray at r=20M, inclined at theta=pi/3 (60 deg)
    ray%x = [0.0_dp, 20.0_dp, PI / 3.0_dp, 0.0_dp]
    ray%p = [-1.0_dp, -0.98_dp, 0.15_dp, 3.5_dp]
    ray%affine_lambda = 0.0_dp
    ray%terminated = .false.
    ray%captured_by_horizon = .false.
    ray%escaped_to_infinity = .false.

    ! Evaluate initial Hamiltonian and Carter constant
    call geodesic_hamiltonian(metric, ray%x, ray%p, h_val, err)
    call assert_ok(suite, err, "Geodesic Hamiltonian evaluation status OK")

    q_init = compute_carter_constant(metric%mass, metric%spin, ray%x(1), ray%x(2), &
                                     ray%p(0), ray%p(1), ray%p(2), ray%p(3))

    ! Integrate 20 RK4 steps
    do step = 1, 20
      call integrate_geodesic_step_rk4(metric, ray, d_lambda, r_horiz, r_escape, err)
      if (ray%terminated) exit
    end do

    ! Evaluate final Carter constant
    q_now = compute_carter_constant(metric%mass, metric%spin, ray%x(1), ray%x(2), &
                                    ray%p(0), ray%p(1), ray%p(2), ray%p(3))

    ! Test 1: Carter constant conserved to within RK4 local discretization error (< 5e-4)
    call assert_close_dp(suite, q_now, q_init, 5.0e-4_dp, "Carter constant invariant conservation")

    ! Test 2: Ray propagates inward toward black hole (r decreases initially)
    call assert_true(suite, ray%x(1) < 20.0_dp, "Null ray propagates inward toward black hole")
  end subroutine run_test_geodesics

end module test_geodesics
