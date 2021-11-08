!> Unit tests for spacetime metrics and 3+1 ADM decompositions
module test_spacetime_metric
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: PI
  use umbra_metric_base, only: metric_3p1_t
  use umbra_minkowski, only: minkowski_metric_t
  use umbra_schwarzschild, only: schwarzschild_metric_t
  use umbra_kerr_bl, only: kerr_bl_metric_t
  use umbra_kerr_ks, only: kerr_ks_metric_t
  use test_framework, only: test_suite_t, assert_true, assert_close_dp
  implicit none
  private

  public :: run_test_spacetime_metric

contains

  subroutine run_test_spacetime_metric(suite)
    type(test_suite_t), intent(inout) :: suite
    type(minkowski_metric_t) :: mink
    type(schwarzschild_metric_t) :: schw
    type(kerr_bl_metric_t) :: kerr
    type(kerr_ks_metric_t) :: ks
    type(metric_3p1_t) :: m
    real(dp), dimension(0:3) :: x
    real(dp) :: r_plus, r_minus, r_ergo

    ! Test 1: Minkowski flat metric determinant and lapse
    x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call mink%eval(x, m)
    call assert_close_dp(suite, m%alpha, 1.0_dp, 1.0e-14_dp, "Minkowski lapse alpha = 1.0")
    call assert_close_dp(suite, m%det_gamma, 1.0_dp, 1.0e-14_dp, "Minkowski det(gamma) = 1.0")
    call assert_close_dp(suite, m%sqrt_det_gamma, 1.0_dp, 1.0e-14_dp, "Minkowski sqrt(det(gamma)) = 1.0")

    ! Test 2: Schwarzschild lapse at r = 4M (alpha = sqrt(1 - 2M/r) = sqrt(0.5) ~ 0.7071)
    schw%mass = 1.0_dp
    x = [0.0_dp, 4.0_dp, 0.5_dp * PI, 0.0_dp]
    call schw%eval(x, m)
    call assert_close_dp(suite, m%alpha, sqrt(0.5_dp), 1.0e-12_dp, "Schwarzschild lapse at r=4M")

    ! Test 3: Kerr BL horizons for a=0.8M (r+ = 1 + sqrt(1 - 0.64) = 1.6, r- = 0.4)
    kerr%mass = 1.0_dp
    kerr%spin = 0.8_dp
    r_plus = kerr%horizon_outer()
    r_minus = kerr%horizon_inner()
    call assert_close_dp(suite, r_plus, 1.6_dp, 1.0e-14_dp, "Kerr outer horizon r+ = 1.6")
    call assert_close_dp(suite, r_minus, 0.4_dp, 1.0e-14_dp, "Kerr inner horizon r- = 0.4")

    ! Test 4: Kerr ergosphere on equator theta=pi/2 (r_E = 2.0 for all spin)
    r_ergo = kerr%ergosphere_radius(0.5_dp * PI)
    call assert_close_dp(suite, r_ergo, 2.0_dp, 1.0e-14_dp, "Kerr equatorial ergosphere r_E = 2.0")

    ! Test 5: Kerr-Schild regularity at horizon r = r+
    ks%mass = 1.0_dp
    ks%spin = 0.8_dp
    x = [0.0_dp, 1.6_dp, 0.0_dp, 0.0_dp] ! On horizon
    call ks%eval(x, m)
    call assert_true(suite, m%alpha > 0.0_dp, "Kerr-Schild horizon lapse > 0")
    call assert_true(suite, m%sqrt_det_gamma > 0.0_dp, "Kerr-Schild horizon sqrt_det_gamma > 0")
  end subroutine run_test_spacetime_metric

end module test_spacetime_metric
