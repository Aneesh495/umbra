! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_carter_constant
! Description: Carter constant invariant evaluation and analytical Bardeen photon ring boundaries.

module umbra_carter_constant
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: EPSILON_DP, PI
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  public :: compute_carter_constant
  public :: compute_analytic_shadow_curve

contains

  pure function compute_carter_constant(mass, spin, r, theta, pt, pr, ptheta, pphi) result(q_val)
    ! Computes Carter constant Q = p_theta^2 + cos^2(theta) [ a^2 (-p_t^2) + p_phi^2 / sin^2(theta) ]
    real(dp), intent(in) :: mass, spin, r, theta, pt, pr, ptheta, pphi
    real(dp) :: q_val, sin_th, cos_th, e_val, lz_val, a_val

    a_val = spin
    sin_th = max(abs(sin(theta)), 1.0e-12_dp)
    cos_th = cos(theta)

    e_val = -pt
    lz_val = pphi

    ! Carter constant Q for null geodesics
    q_val = ptheta**2 + (cos_th**2) * (-a_val * a_val * (e_val**2) + (lz_val**2) / (sin_th**2)) + 0.0_dp * (mass + r + pr)
  end function compute_carter_constant

  subroutine compute_analytic_shadow_curve(mass, spin, theta_obs, n_points, alpha_pts, beta_pts, err)
    ! Computes the exact analytical Bardeen shadow curve (alpha, beta) on observer screen
    real(dp), intent(in) :: mass, spin, theta_obs
    integer(i4), intent(in) :: n_points
    real(dp), dimension(n_points), intent(out) :: alpha_pts, beta_pts
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: a_val, m_val, a_star, r_pro, r_ret, r, dr
    real(dp) :: delta, xi_c, eta_c, sin_obs, cos_obs, beta_sq
    integer(i4) :: k

    m_val = mass
    a_val = spin
    a_star = a_val / max(m_val, 1.0e-12_dp)
    a_star = max(min(a_star, 0.999999_dp), -0.999999_dp)

    sin_obs = max(abs(sin(theta_obs)), 1.0e-12_dp)
    cos_obs = cos(theta_obs)

    ! Photon orbit radial limits: r_prograde and r_retrograde
    r_pro = 2.0_dp * m_val * (1.0_dp + cos((2.0_dp / 3.0_dp) * acos(-a_star)))
    r_ret = 2.0_dp * m_val * (1.0_dp + cos((2.0_dp / 3.0_dp) * acos( a_star)))

    dr = (r_ret - r_pro) / real(n_points - 1, dp)

    do k = 1, n_points
      r = r_pro + real(k - 1, dp) * dr
      delta = r * r - 2.0_dp * m_val * r + a_val * a_val

      if (abs(a_val) < 1.0e-10_dp) then
        ! Schwarzschild limit: circular shadow of radius sqrt(27) M
        alpha_pts(k) = sqrt(27.0_dp) * m_val * cos(PI * real(k - 1, dp) / real(n_points - 1, dp))
        beta_pts(k)  = sqrt(27.0_dp) * m_val * sin(PI * real(k - 1, dp) / real(n_points - 1, dp))
      else
        ! Critical impact parameters xi_c and eta_c (Chandrasekhar 1983)
        xi_c = (m_val * (r * r - a_val * a_val) - r * delta) / (a_val * (r - m_val))
        eta_c = (r**3) * (4.0_dp * m_val * delta - r * (r - m_val)**2) / ((a_val * (r - m_val))**2)

        ! Celestial coordinates on screen (alpha, beta)
        alpha_pts(k) = -xi_c / sin_obs
        beta_sq = eta_c + a_val * a_val * cos_obs * cos_obs - (xi_c**2) * (cos_obs**2) / (sin_obs**2)
        beta_pts(k) = sqrt(max(beta_sq, 0.0_dp))
      end if
    end do

    if (present(err)) call err%clear()
  end subroutine compute_analytic_shadow_curve

end module umbra_carter_constant
