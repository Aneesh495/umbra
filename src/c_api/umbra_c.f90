!> C Foreign Function Interface (FFI) bindings for umbra kernel
!> Exposes C-compatible ABI symbols via ISO_C_BINDING for Python ctypes and C runtimes
module umbra_c
  use, intrinsic :: iso_c_binding, only: c_int32_t, c_double, c_bool, c_null_char
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, NCONS
  use umbra_state, only: prim_state_t, cons_state_t
  use umbra_metric_base, only: metric_3p1_t
  use umbra_minkowski, only: minkowski_metric_t
  use umbra_kerr_bl, only: kerr_bl_metric_t
  use umbra_eos_ideal, only: ideal_gas_eos_t
  use umbra_con2prim_1d, only: con2prim_1d_solve
  use umbra_con2prim_2d, only: con2prim_2d_solve
  use umbra_riemann_llf, only: riemann_solve_llf
  use umbra_riemann_hll, only: riemann_solve_hll
  use umbra_riemann_hlld, only: riemann_solve_hlld
  use umbra_geodesic_integrator, only: ray_state_t, integrate_geodesic_step_rk4
  use umbra_camera, only: camera_t
  use umbra_bench_kerr_shadow, only: kerr_shadow_boundary_analytic, kerr_shadow_radius_schwarzschild
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  public :: umbra_version_major
  public :: umbra_version_minor
  public :: umbra_version_patch
  public :: umbra_kerr_shadow_schwarzschild_radius
  public :: umbra_kerr_shadow_analytic
  public :: umbra_c_con2prim_1d
  public :: umbra_c_con2prim_2d
  public :: umbra_c_riemann_llf
  public :: umbra_c_riemann_hll
  public :: umbra_c_riemann_hlld
  public :: umbra_c_geodesic_step
  public :: umbra_c_render_shadow

contains

  !> Returns library major version
  function umbra_version_major() bind(C, name="umbra_version_major") result(ver)
    integer(c_int32_t) :: ver
    ver = 1_c_int32_t
  end function umbra_version_major

  !> Returns library minor version
  function umbra_version_minor() bind(C, name="umbra_version_minor") result(ver)
    integer(c_int32_t) :: ver
    ver = 0_c_int32_t
  end function umbra_version_minor

  !> Returns library patch version
  function umbra_version_patch() bind(C, name="umbra_version_patch") result(ver)
    integer(c_int32_t) :: ver
    ver = 0_c_int32_t
  end function umbra_version_patch

  !> Returns Schwarzschild shadow radius in geometric units
  function umbra_kerr_shadow_schwarzschild_radius() bind(C, name="umbra_kerr_shadow_schwarzschild_radius") result(r_shad)
    real(c_double) :: r_shad
    r_shad = real(kerr_shadow_radius_schwarzschild(), c_double)
  end function umbra_kerr_shadow_schwarzschild_radius

  !> Computes Bardeen analytic celestial shadow boundary coordinates
  subroutine umbra_kerr_shadow_analytic(spin, theta_obs, n_points, alpha_out, beta_plus_out, beta_minus_out) &
       bind(C, name="umbra_kerr_shadow_analytic")
    real(c_double), value, intent(in) :: spin, theta_obs
    integer(c_int32_t), value, intent(in) :: n_points
    real(c_double), intent(out) :: alpha_out(n_points)
    real(c_double), intent(out) :: beta_plus_out(n_points)
    real(c_double), intent(out) :: beta_minus_out(n_points)

    real(dp), allocatable :: a_arr(:), bp_arr(:), bm_arr(:)

    allocate(a_arr(n_points), bp_arr(n_points), bm_arr(n_points))
    call kerr_shadow_boundary_analytic(real(spin, dp), real(theta_obs, dp), int(n_points, i4), &
                                       a_arr, bp_arr, bm_arr)
    alpha_out = real(a_arr, c_double)
    beta_plus_out = real(bp_arr, c_double)
    beta_minus_out = real(bm_arr, c_double)
    deallocate(a_arr, bp_arr, bm_arr)
  end subroutine umbra_kerr_shadow_analytic

  !> 1D Noble-Mignone conservative-to-primitive inversion via C ABI
  function umbra_c_con2prim_1d(d, s1, s2, s3, tau, b1, b2, b3, gamma_eos, &
                               rho_out, vx_out, vy_out, vz_out, p_out) &
       bind(C, name="umbra_c_con2prim_1d") result(status)
    real(c_double), value, intent(in) :: d, s1, s2, s3, tau, b1, b2, b3, gamma_eos
    real(c_double), intent(out) :: rho_out, vx_out, vy_out, vz_out, p_out
    integer(c_int32_t) :: status

    type(cons_state_t) :: u
    type(prim_state_t) :: p
    type(metric_3p1_t) :: m
    type(minkowski_metric_t) :: mink
    type(ideal_gas_eos_t) :: eos
    type(umbra_error_t) :: err
    real(dp), dimension(0:3) :: x_origin

    u%d = real(d, dp)
    u%s(1) = real(s1, dp)
    u%s(2) = real(s2, dp)
    u%s(3) = real(s3, dp)
    u%tau = real(tau, dp)
    u%b(1) = real(b1, dp)
    u%b(2) = real(b2, dp)
    u%b(3) = real(b3, dp)

    eos%gamma = real(gamma_eos, dp)
    x_origin = 0.0_dp
    call mink%eval(x_origin, m)

    call con2prim_1d_solve(u, m, eos, p, err)

    if (err%ok()) then
      rho_out = real(p%rho, c_double)
      vx_out = real(p%v(1), c_double)
      vy_out = real(p%v(2), c_double)
      vz_out = real(p%v(3), c_double)
      p_out = real(p%p, c_double)
      status = 0_c_int32_t
    else
      status = int(err%code, c_int32_t)
    end if
  end function umbra_c_con2prim_1d

  !> 2D 2D Newton-Raphson conservative-to-primitive inversion via C ABI
  function umbra_c_con2prim_2d(d, s1, s2, s3, tau, b1, b2, b3, gamma_eos, &
                               rho_out, vx_out, vy_out, vz_out, p_out) &
       bind(C, name="umbra_c_con2prim_2d") result(status)
    real(c_double), value, intent(in) :: d, s1, s2, s3, tau, b1, b2, b3, gamma_eos
    real(c_double), intent(out) :: rho_out, vx_out, vy_out, vz_out, p_out
    integer(c_int32_t) :: status

    type(cons_state_t) :: u
    type(prim_state_t) :: p
    type(metric_3p1_t) :: m
    type(minkowski_metric_t) :: mink
    type(ideal_gas_eos_t) :: eos
    type(umbra_error_t) :: err
    real(dp), dimension(0:3) :: x_origin

    u%d = real(d, dp)
    u%s(1) = real(s1, dp)
    u%s(2) = real(s2, dp)
    u%s(3) = real(s3, dp)
    u%tau = real(tau, dp)
    u%b(1) = real(b1, dp)
    u%b(2) = real(b2, dp)
    u%b(3) = real(b3, dp)

    eos%gamma = real(gamma_eos, dp)
    x_origin = 0.0_dp
    call mink%eval(x_origin, m)

    call con2prim_2d_solve(u, m, eos, p, err)

    if (err%ok()) then
      rho_out = real(p%rho, c_double)
      vx_out = real(p%v(1), c_double)
      vy_out = real(p%v(2), c_double)
      vz_out = real(p%v(3), c_double)
      p_out = real(p%p, c_double)
      status = 0_c_int32_t
    else
      status = int(err%code, c_int32_t)
    end if
  end function umbra_c_con2prim_2d

  !> Local Lax-Friedrichs (Rusanov) Riemann solver via C ABI
  subroutine umbra_c_riemann_llf(prim_l, prim_r, gamma_eos, dir, flux_out) &
       bind(C, name="umbra_c_riemann_llf")
    real(c_double), intent(in) :: prim_l(NPRIM), prim_r(NPRIM)
    real(c_double), value, intent(in) :: gamma_eos
    integer(c_int32_t), value, intent(in) :: dir
    real(c_double), intent(out) :: flux_out(NCONS)

    type(prim_state_t) :: pl, pr
    type(metric_3p1_t) :: m
    type(minkowski_metric_t) :: mink
    type(ideal_gas_eos_t) :: eos
    real(dp), dimension(NCONS) :: f
    real(dp) :: cmax
    real(dp), dimension(0:3) :: x_origin

    call pl%from_array(real(prim_l, dp))
    call pr%from_array(real(prim_r, dp))
    eos%gamma = real(gamma_eos, dp)
    x_origin = 0.0_dp
    call mink%eval(x_origin, m)

    call riemann_solve_llf(pl, pr, m, eos, int(dir, i4), f, cmax)
    flux_out = real(f, c_double)
  end subroutine umbra_c_riemann_llf

  !> Harten-Lax-van Leer (HLL) two-wave Riemann solver via C ABI
  subroutine umbra_c_riemann_hll(prim_l, prim_r, gamma_eos, dir, flux_out) &
       bind(C, name="umbra_c_riemann_hll")
    real(c_double), intent(in) :: prim_l(NPRIM), prim_r(NPRIM)
    real(c_double), value, intent(in) :: gamma_eos
    integer(c_int32_t), value, intent(in) :: dir
    real(c_double), intent(out) :: flux_out(NCONS)

    type(prim_state_t) :: pl, pr
    type(metric_3p1_t) :: m
    type(minkowski_metric_t) :: mink
    type(ideal_gas_eos_t) :: eos
    real(dp), dimension(NCONS) :: f
    real(dp) :: cmax
    real(dp), dimension(0:3) :: x_origin

    call pl%from_array(real(prim_l, dp))
    call pr%from_array(real(prim_r, dp))
    eos%gamma = real(gamma_eos, dp)
    x_origin = 0.0_dp
    call mink%eval(x_origin, m)

    call riemann_solve_hll(pl, pr, m, eos, int(dir, i4), f, cmax)
    flux_out = real(f, c_double)
  end subroutine umbra_c_riemann_hll

  !> HLLD contact-restoring five-state Riemann solver via C ABI
  subroutine umbra_c_riemann_hlld(prim_l, prim_r, gamma_eos, dir, flux_out) &
       bind(C, name="umbra_c_riemann_hlld")
    real(c_double), intent(in) :: prim_l(NPRIM), prim_r(NPRIM)
    real(c_double), value, intent(in) :: gamma_eos
    integer(c_int32_t), value, intent(in) :: dir
    real(c_double), intent(out) :: flux_out(NCONS)

    type(prim_state_t) :: pl, pr
    type(metric_3p1_t) :: m
    type(minkowski_metric_t) :: mink
    type(ideal_gas_eos_t) :: eos
    real(dp), dimension(NCONS) :: f
    real(dp) :: cmax
    real(dp), dimension(0:3) :: x_origin

    call pl%from_array(real(prim_l, dp))
    call pr%from_array(real(prim_r, dp))
    eos%gamma = real(gamma_eos, dp)
    x_origin = 0.0_dp
    call mink%eval(x_origin, m)

    call riemann_solve_hlld(pl, pr, m, eos, int(dir, i4), f, cmax)
    flux_out = real(f, c_double)
  end subroutine umbra_c_riemann_hlld

  !> Single null geodesic RK4 step in Kerr spacetime via C ABI
  function umbra_c_geodesic_step(spin, d_lambda, r_horizon, r_escape, x, p, terminated, captured) &
       bind(C, name="umbra_c_geodesic_step") result(status)
    real(c_double), value, intent(in) :: spin, d_lambda, r_horizon, r_escape
    real(c_double), intent(inout) :: x(4), p(4)
    logical(c_bool), intent(out) :: terminated, captured
    integer(c_int32_t) :: status

    type(kerr_bl_metric_t) :: metric
    type(ray_state_t) :: ray
    type(umbra_error_t) :: err

    metric%mass = 1.0_dp
    metric%spin = real(spin, dp)

    ray%x = real(x, dp)
    ray%p = real(p, dp)
    ray%terminated = .false.
    ray%captured_by_horizon = .false.
    ray%escaped_to_infinity = .false.

    call integrate_geodesic_step_rk4(metric, ray, real(d_lambda, dp), &
                                     real(r_horizon, dp), real(r_escape, dp), err)

    x = real(ray%x, c_double)
    p = real(ray%p, c_double)
    terminated = logical(ray%terminated, c_bool)
    captured = logical(ray%captured_by_horizon, c_bool)

    if (err%ok()) then
      status = 0_c_int32_t
    else
      status = int(err%code, c_int32_t)
    end if
  end function umbra_c_geodesic_step

  !> Renders 2D black hole shadow mask using raytracing via C ABI
  subroutine umbra_c_render_shadow(spin, inclination_deg, fov_m, npix_x, npix_y, mask_out) &
       bind(C, name="umbra_c_render_shadow")
    real(c_double), value, intent(in) :: spin, inclination_deg, fov_m
    integer(c_int32_t), value, intent(in) :: npix_x, npix_y
    integer(c_int32_t), intent(out) :: mask_out(npix_x, npix_y)

    type(camera_t) :: camera
    type(kerr_bl_metric_t) :: metric
    type(ray_state_t) :: ray
    type(umbra_error_t) :: err
    real(dp) :: r_horiz, r_escape, d_obs, d_lam
    integer(i4) :: ix, iy, step_idx, max_steps

    metric%mass = 1.0_dp
    metric%spin = real(spin, dp)
    r_horiz = metric%horizon_outer() + 0.01_dp
    r_escape = 200.0_dp
    d_obs = 100.0_dp
    d_lam = 0.1_dp
    max_steps = 1500

    call camera%init(d_obs, real(inclination_deg, dp), real(fov_m, dp), &
                     int(npix_x, i4), int(npix_y, i4))

    mask_out = 0_c_int32_t

    do iy = 1, int(npix_y, i4)
      do ix = 1, int(npix_x, i4)
        call camera%spawn_ray(ix, iy, metric, ray)
        do step_idx = 1, max_steps
          call integrate_geodesic_step_rk4(metric, ray, d_lam, r_horiz, r_escape, err)
          if (ray%terminated) exit
        end do
        if (ray%captured_by_horizon) then
          mask_out(ix, iy) = 1_c_int32_t
        end if
      end do
    end do

    call camera%free()
  end subroutine umbra_c_render_shadow

end module umbra_c
