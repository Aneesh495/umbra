! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_kerr_ks
! Description: Horizon-penetrating Kerr-Schild spacetime metric for stable accretion simulation.

module umbra_kerr_ks
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NDIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_SINGULAR_METRIC
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t, compute_adm_4metric
  use umbra_math, only: mat3_inv
  implicit none
  private

  type, extends(spacetime_metric_t), public :: kerr_ks_metric_t
  contains
    procedure, public :: eval => kerr_ks_eval
    procedure, public :: to_cartesian => kerr_ks_to_cartesian
    procedure, public :: horizon_outer => kerr_ks_horizon_outer
    procedure, public :: horizon_inner => kerr_ks_horizon_inner
  end type kerr_ks_metric_t

contains

  pure function kerr_ks_horizon_outer(this) result(r_plus)
    class(kerr_ks_metric_t), intent(in) :: this
    real(dp) :: r_plus
    r_plus = this%mass + sqrt(max(this%mass**2 - this%spin**2, 0.0_dp))
  end function kerr_ks_horizon_outer

  pure function kerr_ks_horizon_inner(this) result(r_minus)
    class(kerr_ks_metric_t), intent(in) :: this
    real(dp) :: r_minus
    r_minus = this%mass - sqrt(max(this%mass**2 - this%spin**2, 0.0_dp))
  end function kerr_ks_horizon_inner

  subroutine kerr_ks_eval(this, x, m, err)
    class(kerr_ks_metric_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: x
    type(metric_3p1_t), intent(out) :: m
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: r, th, sin_th, cos_th, sin2_th
    real(dp) :: m_val, a_val
    real(dp) :: sigma, delta, a_big, two_m_r, z

    m_val = this%mass
    a_val = this%spin

    r = max(x(1), 1.0e-6_dp)
    th = x(2)
    sin_th = max(abs(sin(th)), 1.0e-12_dp)
    cos_th = cos(th)
    sin2_th = sin_th * sin_th

    sigma = r * r + a_val * a_val * cos_th * cos_th
    delta = r * r - 2.0_dp * m_val * r + a_val * a_val
    two_m_r = 2.0_dp * m_val * r
    z = two_m_r / sigma

    a_big = (r * r + a_val * a_val)**2 - a_val * a_val * delta * sin2_th

    ! Kerr-Schild 4-metric g_munu components in spherical coords:
    ! g_00 = -(1 - z)
    ! g_01 = z
    ! g_03 = -z * a * sin^2(th)
    ! g_11 = 1 + z
    ! g_13 = -a * sin^2(th) * (1 + z)
    ! g_22 = sigma
    ! g_33 = (A / sigma) * sin^2(th) = sin^2(th) * (sigma + a^2 sin^2(th) * (1 + z))
    m%g_down = 0.0_dp
    m%g_down(0,0) = -(1.0_dp - z)
    m%g_down(0,1) = z
    m%g_down(1,0) = z
    m%g_down(0,3) = -z * a_val * sin2_th
    m%g_down(3,0) = m%g_down(0,3)

    m%g_down(1,1) = 1.0_dp + z
    m%g_down(1,3) = -a_val * sin2_th * (1.0_dp + z)
    m%g_down(3,1) = m%g_down(1,3)

    m%g_down(2,2) = sigma
    m%g_down(3,3) = sin2_th * (sigma + a_val * a_val * sin2_th * (1.0_dp + z))

    ! Spatial 3-metric gamma_ij is the spatial submatrix of g_munu
    m%gamma_down(1,1) = m%g_down(1,1)
    m%gamma_down(1,2) = 0.0_dp
    m%gamma_down(1,3) = m%g_down(1,3)

    m%gamma_down(2,1) = 0.0_dp
    m%gamma_down(2,2) = m%g_down(2,2)
    m%gamma_down(2,3) = 0.0_dp

    m%gamma_down(3,1) = m%g_down(3,1)
    m%gamma_down(3,2) = 0.0_dp
    m%gamma_down(3,3) = m%g_down(3,3)

    ! Determinant det(gamma) = sigma * (1 + z) * [sigma * sin^2(th)] = (1 + z) * sigma^2 * sin^2(th)
    ! Volume factor sqrt(-g) = sigma * sin(theta)
    m%sqrt_det_g = sigma * sin_th

    ! Invert gamma_ij
    call mat3_inv(m%gamma_down, m%gamma_up, m%det_gamma, err)
    if (m%det_gamma > 0.0_dp) then
      m%sqrt_det_gamma = sqrt(m%det_gamma)
    else
      m%sqrt_det_gamma = 1.0_dp
    end if

    ! Lapse alpha and shift beta in Kerr-Schild:
    ! alpha = 1 / sqrt(1 + z)
    ! beta^1 = z / (1 + z) = 2 M r / (sigma + 2 M r)
    ! beta^2 = 0
    ! beta^3 = 0
    m%alpha = 1.0_dp / sqrt(1.0_dp + z)
    m%beta_up = 0.0_dp
    m%beta_up(1) = z / (1.0_dp + z)

    ! Covariant shift beta_i = gamma_ij beta^j = g_0i
    m%beta_down(1) = m%g_down(0,1)
    m%beta_down(2) = 0.0_dp
    m%beta_down(3) = m%g_down(0,3)

    ! Contravariant 4-metric g^munu:
    ! In Kerr-Schild, g^munu = eta^munu - 2 H l^mu l^nu
    ! For spherical coords, we compute from ADM relation
    m%g_up(0,0) = -1.0_dp / (m%alpha * m%alpha)
    m%g_up(0,1) = m%beta_up(1) / (m%alpha * m%alpha)
    m%g_up(1,0) = m%g_up(0,1)
    m%g_up(0,2) = 0.0_dp
    m%g_up(2,0) = 0.0_dp
    m%g_up(0,3) = 0.0_dp
    m%g_up(3,0) = 0.0_dp

    m%g_up(1,1) = m%gamma_up(1,1) - (m%beta_up(1)**2) / (m%alpha * m%alpha)
    m%g_up(1,2) = m%gamma_up(1,2)
    m%g_up(2,1) = m%g_up(1,2)
    m%g_up(1,3) = m%gamma_up(1,3)
    m%g_up(3,1) = m%g_up(1,3)

    m%g_up(2,2) = m%gamma_up(2,2)
    m%g_up(2,3) = 0.0_dp
    m%g_up(3,2) = 0.0_dp

    m%g_up(3,3) = m%gamma_up(3,3)

    if (present(err)) call err%clear()
  end subroutine kerr_ks_eval

  pure subroutine kerr_ks_to_cartesian(this, x, x_cart)
    class(kerr_ks_metric_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: x
    real(dp), dimension(0:3), intent(out) :: x_cart

    real(dp) :: r, th, phi, a_val

    a_val = this%spin
    r = x(1)
    th = x(2)
    phi = x(3)

    x_cart(0) = x(0)
    x_cart(1) = (r * cos(phi) - a_val * sin(phi)) * sin(th)
    x_cart(2) = (r * sin(phi) + a_val * cos(phi)) * sin(th)
    x_cart(3) = r * cos(th)
  end subroutine kerr_ks_to_cartesian

end module umbra_kerr_ks
