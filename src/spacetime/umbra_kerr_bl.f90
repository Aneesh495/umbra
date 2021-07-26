! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_kerr_bl
! Description: Rotating Kerr black hole metric in standard Boyer-Lindquist coordinates.

module umbra_kerr_bl
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NDIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_HORIZON_PENETRATED, &
                         UMBRA_ERR_SINGULAR_METRIC
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t, compute_adm_4metric
  implicit none
  private

  type, extends(spacetime_metric_t), public :: kerr_bl_metric_t
  contains
    procedure, public :: eval => kerr_bl_eval
    procedure, public :: to_cartesian => kerr_bl_to_cartesian
    procedure, public :: horizon_outer => kerr_bl_horizon_outer
    procedure, public :: horizon_inner => kerr_bl_horizon_inner
    procedure, public :: ergosphere_radius => kerr_bl_ergosphere_radius
    procedure, public :: photon_orbit_radii => kerr_bl_photon_orbit_radii
  end type kerr_bl_metric_t

contains

  pure function kerr_bl_horizon_outer(this) result(r_plus)
    class(kerr_bl_metric_t), intent(in) :: this
    real(dp) :: r_plus
    r_plus = this%mass + sqrt(max(this%mass**2 - this%spin**2, 0.0_dp))
  end function kerr_bl_horizon_outer

  pure function kerr_bl_horizon_inner(this) result(r_minus)
    class(kerr_bl_metric_t), intent(in) :: this
    real(dp) :: r_minus
    r_minus = this%mass - sqrt(max(this%mass**2 - this%spin**2, 0.0_dp))
  end function kerr_bl_horizon_inner

  pure function kerr_bl_ergosphere_radius(this, theta) result(r_e)
    class(kerr_bl_metric_t), intent(in) :: this
    real(dp), intent(in) :: theta
    real(dp) :: r_e
    r_e = this%mass + sqrt(max(this%mass**2 - this%spin**2 * cos(theta)**2, 0.0_dp))
  end function kerr_bl_ergosphere_radius

  pure subroutine kerr_bl_photon_orbit_radii(this, r_pro, r_ret)
    class(kerr_bl_metric_t), intent(in) :: this
    real(dp), intent(out) :: r_pro, r_ret
    real(dp) :: a_star

    a_star = this%spin / max(this%mass, 1.0e-12_dp)
    a_star = max(min(a_star, 0.999999_dp), -0.999999_dp)

    r_pro = 2.0_dp * this%mass * (1.0_dp + cos((2.0_dp / 3.0_dp) * acos(-a_star)))
    r_ret = 2.0_dp * this%mass * (1.0_dp + cos((2.0_dp / 3.0_dp) * acos( a_star)))
  end subroutine kerr_bl_photon_orbit_radii

  subroutine kerr_bl_eval(this, x, m, err)
    class(kerr_bl_metric_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: x
    type(metric_3p1_t), intent(out) :: m
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: r, th, sin_th, cos_th, sin2_th
    real(dp) :: m_val, a_val, r_plus
    real(dp) :: sigma, delta, a_big

    m_val = this%mass
    a_val = this%spin
    r_plus = this%horizon_outer()

    r = x(1)
    th = x(2)
    sin_th = max(abs(sin(th)), 1.0e-12_dp)
    cos_th = cos(th)
    sin2_th = sin_th * sin_th

    if (r <= r_plus + 1.0e-5_dp) then
      if (present(err)) call err%set(UMBRA_ERR_HORIZON_PENETRATED, &
        "Point inside or touching outer event horizon in Boyer-Lindquist coordinates", "kerr_bl_eval")
    end if

    sigma = r * r + a_val * a_val * cos_th * cos_th
    delta = r * r - 2.0_dp * m_val * r + a_val * a_val
    if (delta < 1.0e-14_dp) delta = 1.0e-14_dp

    a_big = (r * r + a_val * a_val)**2 - a_val * a_val * delta * sin2_th

    ! 3+1 ADM lapse alpha = sqrt(sigma * delta / A)
    m%alpha = sqrt(sigma * delta / a_big)

    ! Shift vector: beta^phi = -omega = -2 M a r / A
    m%beta_up = 0.0_dp
    m%beta_up(3) = - (2.0_dp * m_val * a_val * r) / a_big

    ! Spatial 3-metric gamma_ij:
    ! gamma_rr = sigma / delta
    ! gamma_thth = sigma
    ! gamma_phiphi = (A / sigma) * sin^2(theta)
    m%gamma_down = 0.0_dp
    m%gamma_down(1,1) = sigma / delta
    m%gamma_down(2,2) = sigma
    m%gamma_down(3,3) = (a_big / sigma) * sin2_th

    ! Spatial determinant det(gamma) = A * sigma * sin^2(th) / delta
    m%det_gamma = (a_big * sigma * sin2_th) / delta
    m%sqrt_det_gamma = sqrt(m%det_gamma)

    ! Volume element: alpha * sqrt(det_gamma) = sigma * sin(theta)
    m%sqrt_det_g = sigma * sin_th

    ! Covariant shift beta_i = gamma_ij beta^j
    m%beta_down = 0.0_dp
    m%beta_down(3) = -(2.0_dp * m_val * a_val * r * sin2_th) / sigma

    ! Contravariant spatial metric gamma^ij
    m%gamma_up = 0.0_dp
    m%gamma_up(1,1) = delta / sigma
    m%gamma_up(2,2) = 1.0_dp / sigma
    m%gamma_up(3,3) = sigma / (a_big * sin2_th)

    ! Covariant 4-metric g_munu
    m%g_down = 0.0_dp
    m%g_down(0,0) = -(1.0_dp - (2.0_dp * m_val * r) / sigma)
    m%g_down(0,3) = -(2.0_dp * m_val * a_val * r * sin2_th) / sigma
    m%g_down(3,0) = m%g_down(0,3)
    m%g_down(1,1) = m%gamma_down(1,1)
    m%g_down(2,2) = m%gamma_down(2,2)
    m%g_down(3,3) = m%gamma_down(3,3)

    ! Contravariant 4-metric g^munu
    m%g_up = 0.0_dp
    m%g_up(0,0) = -a_big / (sigma * delta)
    m%g_up(0,3) = -(2.0_dp * m_val * a_val * r) / (sigma * delta)
    m%g_up(3,0) = m%g_up(0,3)
    m%g_up(1,1) = delta / sigma
    m%g_up(2,2) = 1.0_dp / sigma
    m%g_up(3,3) = (delta - a_val * a_val * sin2_th) / (sigma * delta * sin2_th)

    if (present(err)) call err%clear()
  end subroutine kerr_bl_eval

  pure subroutine kerr_bl_to_cartesian(this, x, x_cart)
    class(kerr_bl_metric_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: x
    real(dp), dimension(0:3), intent(out) :: x_cart

    real(dp) :: r, th, phi, a_val, rho_cyl

    a_val = this%spin
    r = x(1)
    th = x(2)
    phi = x(3)

    ! Kerr-Schild / Cartesian embedding of Boyer-Lindquist coords:
    ! x = sqrt(r^2 + a^2) * sin(th) * cos(phi)
    ! y = sqrt(r^2 + a^2) * sin(th) * sin(phi)
    ! z = r * cos(th)
    rho_cyl = sqrt(r * r + a_val * a_val)
    x_cart(0) = x(0)
    x_cart(1) = rho_cyl * sin(th) * cos(phi)
    x_cart(2) = rho_cyl * sin(th) * sin(phi)
    x_cart(3) = r * cos(th)
  end subroutine kerr_bl_to_cartesian

end module umbra_kerr_bl
