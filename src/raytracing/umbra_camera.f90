! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_camera
! Description: Virtual pinhole camera setup, impact parameter grids, and initial ray generation.

module umbra_camera
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: PI
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t
  use umbra_geodesic_integrator, only: ray_state_t
  implicit none
  private

  type, public :: camera_t
    real(dp) :: r_obs = 1000.0_dp              ! Observer radial distance in M
    real(dp) :: theta_obs = 1.0471975511965977_dp ! Inclination angle in radians (60 degrees)
    real(dp) :: phi_obs = 0.0_dp               ! Observer azimuth
    real(dp) :: fov = 30.0_dp                  ! Field of view in M
    integer(i4) :: npix_x = 64_i4              ! Horizontal screen resolution
    integer(i4) :: npix_y = 64_i4              ! Vertical screen resolution

    real(dp), allocatable, dimension(:) :: alpha_grid
    real(dp), allocatable, dimension(:) :: beta_grid
    real(dp), allocatable, dimension(:,:) :: intensity_map
  contains
    procedure, public :: init => camera_init
    procedure, public :: free => camera_free
    procedure, public :: spawn_ray => camera_spawn_ray
  end type camera_t

contains

  subroutine camera_init(this, r_dist, inc_deg, fov_m, nx, ny, err)
    class(camera_t), intent(inout) :: this
    real(dp), intent(in) :: r_dist, inc_deg, fov_m
    integer(i4), intent(in) :: nx, ny
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: d_alpha, d_beta
    integer(i4) :: i, j

    this%r_obs = r_dist
    this%theta_obs = inc_deg * (PI / 180.0_dp)
    this%fov = fov_m
    this%npix_x = nx
    this%npix_y = ny

    call this%free()

    allocate(this%alpha_grid(nx))
    allocate(this%beta_grid(ny))
    allocate(this%intensity_map(nx, ny))
    this%intensity_map = 0.0_dp

    d_alpha = fov_m / real(nx, dp)
    d_beta  = fov_m / real(ny, dp)

    do i = 1, nx
      this%alpha_grid(i) = -0.5_dp * fov_m + (real(i, dp) - 0.5_dp) * d_alpha
    end do

    do j = 1, ny
      this%beta_grid(j) = -0.5_dp * fov_m + (real(j, dp) - 0.5_dp) * d_beta
    end do

    if (present(err)) call err%clear()
  end subroutine camera_init

  subroutine camera_free(this)
    class(camera_t), intent(inout) :: this
    if (allocated(this%alpha_grid)) deallocate(this%alpha_grid)
    if (allocated(this%beta_grid)) deallocate(this%beta_grid)
    if (allocated(this%intensity_map)) deallocate(this%intensity_map)
  end subroutine camera_free

  subroutine camera_spawn_ray(this, ix, iy, metric, ray, err)
    ! Spawns a ray backward in time from pixel (ix, iy)
    class(camera_t), intent(in) :: this
    integer(i4), intent(in) :: ix, iy
    class(spacetime_metric_t), intent(in) :: metric
    type(ray_state_t), intent(out) :: ray
    type(umbra_error_t), intent(inout), optional :: err

    type(metric_3p1_t) :: m
    real(dp) :: alpha_val, beta_val, sin_obs
    real(dp) :: pt, pphi, ptheta, pr_sq

    alpha_val = this%alpha_grid(ix)
    beta_val  = this%beta_grid(iy)
    sin_obs = max(abs(sin(this%theta_obs)), 1.0e-12_dp)

    ! Initial 4-position x^mu = (0, r_obs, theta_obs, phi_obs)
    ray%x(0) = 0.0_dp
    ray%x(1) = this%r_obs
    ray%x(2) = this%theta_obs
    ray%x(3) = this%phi_obs
    ray%affine_lambda = 0.0_dp
    ray%terminated = .false.
    ray%captured_by_horizon = .false.
    ray%escaped_to_infinity = .false.

    call metric%eval(ray%x, m, err)

    ! Standard impact parameter momentum components:
    ! E = -p_t = 1
    ! L_z = p_phi = -alpha * sin(theta_obs)
    ! p_theta = beta
    pt = -1.0_dp
    pphi = -alpha_val * sin_obs
    ptheta = beta_val

    ! Compute p_r such that null condition holds: g^{mu nu} p_mu p_nu = 0
    ! g^{tt} pt^2 + 2 g^{tphi} pt pphi + g^{rr} pr^2 + g^{thth} ptheta^2 + g^{phiphi} pphi^2 = 0
    pr_sq = -(m%g_up(0,0)*pt*pt + 2.0_dp*m%g_up(0,3)*pt*pphi + &
              m%g_up(2,2)*ptheta*ptheta + m%g_up(3,3)*pphi*pphi) / max(m%g_up(1,1), 1.0e-30_dp)

    ray%p(0) = pt
    ray%p(1) = -sqrt(max(pr_sq, 0.0_dp)) ! Negative sign: ray propagates inward towards the black hole
    ray%p(2) = ptheta
    ray%p(3) = pphi

    if (present(err)) call err%clear()
  end subroutine camera_spawn_ray

end module umbra_camera
