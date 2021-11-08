!> Unit tests for TVD, PPM, and WENO5 spatial reconstruction
module test_reconstruction
  use umbra_kinds, only: dp, i4
  use umbra_limiters, only: LIMITER_MINMOD, LIMITER_MC, LIMITER_VANLEER, slope_limiter
  use umbra_reconstruct_tvd, only: reconstruct_tvd_1d
  use umbra_reconstruct_weno, only: reconstruct_weno5_1d
  use umbra_error, only: umbra_error_t
  use test_framework, only: test_suite_t, assert_true, assert_close_dp, assert_ok
  implicit none
  private

  public :: run_test_reconstruction

contains

  subroutine run_test_reconstruction(suite)
    use umbra_constants, only: NPRIM
    type(test_suite_t), intent(inout) :: suite
    real(dp) :: slope
    real(dp), allocatable :: q(:,:), q_l(:,:), q_r(:,:)
    type(umbra_error_t) :: err
    integer(i4) :: n, nghost, i

    ! Test 1: Slope limiter zeroes out at local extrema (dl > 0, dr < 0)
    slope = slope_limiter(LIMITER_MINMOD, 1.0_dp, -1.0_dp)
    call assert_close_dp(suite, slope, 0.0_dp, 1.0e-14_dp, "Minmod slope at extremum = 0")

    slope = slope_limiter(LIMITER_MC, 1.0_dp, -1.0_dp)
    call assert_close_dp(suite, slope, 0.0_dp, 1.0e-14_dp, "MC slope at extremum = 0")

    ! Test 2: Minmod preserves constant linear slope (dl=dr=1 => slope=1)
    slope = slope_limiter(LIMITER_MINMOD, 1.0_dp, 1.0_dp)
    call assert_close_dp(suite, slope, 1.0_dp, 1.0e-14_dp, "Minmod linear slope = 1")

    ! Test 3: TVD reconstruction on smooth linear field q(v, i) = 2.0 * i
    n = 20_i4
    nghost = 3_i4
    allocate(q(NPRIM, 1 - nghost : n + nghost))
    allocate(q_l(NPRIM, 1 : n + 1), q_r(NPRIM, 1 : n + 1))

    do i = 1 - nghost, n + nghost
      q(:, i) = 2.0_dp * real(i, dp)
    end do

    call reconstruct_tvd_1d(q, n, nghost, LIMITER_MINMOD, q_l, q_r, err)
    call assert_ok(suite, err, "TVD reconstruction status OK")

    ! Interior interface 10 should have q_l(1, 10) == q_r(1, 10) == 19.0
    call assert_close_dp(suite, q_l(1, 10), 19.0_dp, 1.0e-12_dp, "TVD linear left interface")
    call assert_close_dp(suite, q_r(1, 10), 19.0_dp, 1.0e-12_dp, "TVD linear right interface")

    ! Test 4: WENO5 reconstruction on constant state q = 5.0 preserves 5.0 exactly
    q = 5.0_dp
    call reconstruct_weno5_1d(q, n, nghost, q_l, q_r, err)
    call assert_ok(suite, err, "WENO5 reconstruction status OK")
    call assert_close_dp(suite, q_l(1, 10), 5.0_dp, 1.0e-12_dp, "WENO5 constant left state")
    call assert_close_dp(suite, q_r(1, 10), 5.0_dp, 1.0e-12_dp, "WENO5 constant right state")

    deallocate(q, q_l, q_r)
  end subroutine run_test_reconstruction

end module test_reconstruction
