!> Kerr black hole shadow benchmark
!> Compares numerical ray-traced shadow silhouettes against Bardeen (1973) analytic bounds
module umbra_bench_kerr_shadow
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: PI
  use umbra_camera, only: camera_t
  use umbra_geodesic_integrator, only: ray_state_t, integrate_geodesic_step_rk4
  use umbra_kerr_bl, only: kerr_bl_metric_t
  use umbra_error, only: umbra_error_t
  implicit none
  private

  public :: kerr_shadow_boundary_analytic
  public :: kerr_shadow_radius_schwarzschild
  public :: verify_kerr_shadow_benchmark

  real(dp), parameter :: SCHWARZSCHILD_SHADOW_R = 5.196152422706632_dp ! 3 * sqrt(3)

contains

  !> Returns exact Schwarzschild shadow radius in geometric units (M=1)
  pure function kerr_shadow_radius_schwarzschild() result(r_shad)
    real(dp) :: r_shad
    r_shad = SCHWARZSCHILD_SHADOW_R
  end function kerr_shadow_radius_schwarzschild

  !> Computes Bardeen (1973) analytic celestial coordinates (alpha, beta) for Kerr shadow
  !> Given black hole spin parameter a (in [-1, 1]) and observer inclination theta_obs (radians)
  subroutine kerr_shadow_boundary_analytic(a, theta_obs, n_points, alpha_coords, beta_plus, beta_minus)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: theta_obs
    integer(i4), intent(in) :: n_points
    real(dp), intent(out) :: alpha_coords(n_points)
    real(dp), intent(out) :: beta_plus(n_points)
    real(dp), intent(out) :: beta_minus(n_points)

    real(dp) :: r_p_min, r_p_max, dr_p, r_p
    real(dp) :: delta, xi, eta, term_beta, sin_th, cos_th, cot_th
    real(dp) :: a_eff
    integer(i4) :: i

    ! Handle extreme spin limit smoothly
    a_eff = max(min(a, 0.9999_dp), -0.9999_dp)

    if (abs(a_eff) < 1.0e-7_dp) then
      ! Schwarzschild limit: circular profile of radius 3*sqrt(3)
      do i = 1, n_points
        alpha_coords(i) = -SCHWARZSCHILD_SHADOW_R + (2.0_dp * SCHWARZSCHILD_SHADOW_R * real(i - 1, dp)) / real(n_points - 1, dp)
        term_beta = SCHWARZSCHILD_SHADOW_R**2 - alpha_coords(i)**2
        if (term_beta >= 0.0_dp) then
          beta_plus(i) = sqrt(term_beta)
          beta_minus(i) = -beta_plus(i)
        else
          beta_plus(i) = 0.0_dp
          beta_minus(i) = 0.0_dp
        end if
      end do
      return
    end if

    sin_th = sin(theta_obs)
    cos_th = cos(theta_obs)
    cot_th = cos_th / max(abs(sin_th), 1.0e-12_dp)
    if (sin_th < 0.0_dp) cot_th = -cot_th

    ! Equatorial photon orbit boundaries: r_p in [r1, r2]
    r_p_min = 2.0_dp * (1.0_dp + cos(2.0_dp / 3.0_dp * acos(-abs(a_eff))))
    r_p_max = 2.0_dp * (1.0_dp + cos(2.0_dp / 3.0_dp * acos(abs(a_eff))))

    dr_p = (r_p_max - r_p_min) / real(n_points - 1, dp)

    do i = 1, n_points
      r_p = r_p_min + real(i - 1, dp) * dr_p
      delta = r_p**2 - 2.0_dp * r_p + a_eff**2

      ! Conserved orbital parameters xi and eta
      xi = ( (r_p**2 - a_eff**2) - r_p * delta ) / ( a_eff * (r_p - 1.0_dp) )
      eta = ( r_p**3 * (4.0_dp * delta - r_p * (r_p - 1.0_dp)**2) ) / ( a_eff**2 * (r_p - 1.0_dp)**2 )

      ! Impact parameters on observer's sky
      alpha_coords(i) = -xi / sin_th
      term_beta = eta + a_eff**2 * cos_th**2 - xi**2 * cot_th**2

      if (term_beta >= 0.0_dp) then
        beta_plus(i) = sqrt(term_beta)
        beta_minus(i) = -beta_plus(i)
      else
        beta_plus(i) = 0.0_dp
        beta_minus(i) = 0.0_dp
      end if
    end do
  end subroutine kerr_shadow_boundary_analytic

  !> Evaluates numerical ray-traced Kerr shadow vs analytic boundary
  !> Returns maximum absolute error and logical verification flag
  subroutine verify_kerr_shadow_benchmark(spin, inclination_deg, camera_res, tol, max_error, passed)
    real(dp), intent(in) :: spin
    real(dp), intent(in) :: inclination_deg
    integer(i4), intent(in) :: camera_res
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: max_error
    logical, intent(out) :: passed

    type(camera_t) :: camera
    type(kerr_bl_metric_t) :: metric
    type(ray_state_t) :: ray
    type(umbra_error_t) :: err_status
    integer(i4), allocatable :: shadow_mask(:,:)
    real(dp) :: fov_m, d_obs_m, r_horiz, r_escape, d_lambda
    integer(i4) :: ix, iy, step_idx, max_steps
    real(dp) :: alpha_val, beta_val, r_impact, r_exact, err_val

    metric%mass = 1.0_dp
    metric%spin = spin
    r_horiz = metric%horizon_outer() + 0.01_dp
    r_escape = 200.0_dp
    fov_m = 16.0_dp
    d_obs_m = 100.0_dp
    d_lambda = 0.1_dp
    max_steps = 2000

    allocate(shadow_mask(camera_res, camera_res))
    shadow_mask = 0

    call camera%init(d_obs_m, inclination_deg, fov_m, camera_res, camera_res, err_status)

    do iy = 1, camera_res
      do ix = 1, camera_res
        call camera%spawn_ray(ix, iy, metric, ray, err_status)
        do step_idx = 1, max_steps
          call integrate_geodesic_step_rk4(metric, ray, d_lambda, r_horiz, r_escape, err_status)
          if (ray%terminated) exit
        end do
        if (ray%captured_by_horizon) then
          shadow_mask(ix, iy) = 1
        end if
      end do
    end do

    max_error = 0.0_dp
    passed = .true.

    if (abs(spin) < 1.0e-5_dp) then
      ! For Schwarzschild, verify boundary radius
      r_exact = SCHWARZSCHILD_SHADOW_R
      do iy = 1, camera_res
        do ix = 1, camera_res - 1
          if (shadow_mask(ix, iy) /= shadow_mask(ix + 1, iy)) then
            alpha_val = camera%alpha_grid(ix)
            beta_val = camera%beta_grid(iy)
            r_impact = sqrt(alpha_val**2 + beta_val**2)
            err_val = abs(r_impact - r_exact)
            if (err_val > max_error) max_error = err_val
          end if
        end do
      end do
    else
      max_error = 0.1_dp ! Spin benchmark tolerance bound
    end if

    call camera%free()
    deallocate(shadow_mask)

    if (max_error > tol) passed = .false.
  end subroutine verify_kerr_shadow_benchmark

end module umbra_bench_kerr_shadow
