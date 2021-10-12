! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_synchrotron
! Description: Thermal relativistic synchrotron emissivity and self-absorption coefficients.

module umbra_synchrotron
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: PI, SI_E, SI_ME, SI_C, SI_KB
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  public :: compute_synchrotron_emissivity
  public :: compute_synchrotron_absorptivity

contains

  pure subroutine compute_synchrotron_emissivity(ne, b_mag, theta_e, nu, sin_pitch, j_nu, j_q, j_v)
    ! Evaluates relativistic thermal synchrotron emissivity j_nu for electrons in units of erg s^-1 cm^-3 Hz^-1 sr^-1
    ! Based on Leung, Gammie & Noble (2011) fitting formulas
    real(dp), intent(in) :: ne           ! Electron number density [cm^-3]
    real(dp), intent(in) :: b_mag        ! Magnetic field strength [Gauss]
    real(dp), intent(in) :: theta_e      ! Dimensionless electron temperature kB Te / (me c^2)
    real(dp), intent(in) :: nu           ! Observing frequency in fluid frame [Hz]
    real(dp), intent(in) :: sin_pitch    ! Sin of pitch angle between ray momentum and B-field
    real(dp), intent(out) :: j_nu        ! Total unpolarized emissivity j_I
    real(dp), intent(out) :: j_q         ! Linear polarization emissivity
    real(dp), intent(out) :: j_v         ! Circular polarization emissivity

    real(dp) :: nu_c, x_nu, exp_factor, pol_deg

    ! Critical synchrotron frequency: nu_c = (3 e B sin(theta_B) / 4 pi m_e c) * theta_e^2
    nu_c = (3.0_dp * SI_E * b_mag * max(sin_pitch, 1.0e-6_dp) / (4.0_dp * PI * SI_ME * SI_C)) * (theta_e * theta_e)
    nu_c = max(nu_c, 1.0e-10_dp)

    x_nu = nu / nu_c

    ! Fitting function for thermal synchrotron spectrum (Mahadevan et al. 1996)
    if (x_nu < 1.0e-3_dp) then
      exp_factor = 4.0505_dp * (x_nu**0.166666666666667_dp)
    else
      exp_factor = (sqrt(x_nu) + 2.0_dp**(11.0_dp / 12.0_dp) * (x_nu**(1.0_dp / 6.0_dp)))**2 &
                 * exp(- (x_nu**(1.0_dp / 3.0_dp)))
    end if

    j_nu = (sqrt(2.0_dp) * PI * (SI_E**2) * ne * nu_c / (3.0_dp * SI_C)) * exp_factor / max(theta_e * theta_e, 1.0e-12_dp)
    j_nu = max(j_nu, 0.0_dp)

    ! Fractional linear and circular polarizations
    pol_deg = (x_nu**(1.0_dp / 3.0_dp) + 0.5_dp) / (x_nu**(1.0_dp / 3.0_dp) + 1.8_dp)
    pol_deg = min(max(pol_deg, 0.0_dp), 0.75_dp)

    j_q = pol_deg * j_nu
    j_v = (2.0_dp / max(theta_e, 1.0e-4_dp)) * (1.0_dp / max(sqrt(x_nu), 1.0e-4_dp)) * j_nu * 0.01_dp
  end subroutine compute_synchrotron_emissivity

  pure subroutine compute_synchrotron_absorptivity(j_nu, nu, theta_e, alpha_nu)
    ! Evaluates absorption coefficient alpha_nu via Kirchhoff's law: alpha_nu = j_nu / B_nu(T)
    real(dp), intent(in) :: j_nu, nu, theta_e
    real(dp), intent(out) :: alpha_nu

    real(dp) :: b_nu, t_e

    ! Rayleigh-Jeans limit: B_nu(T) = 2 nu^2 kB Te / c^2 = 2 nu^2 m_e theta_e
    t_e = max(theta_e, 1.0e-6_dp)
    b_nu = (2.0_dp * (nu**2) * SI_ME * (SI_C**2) * t_e) / (SI_C**2)
    b_nu = max(b_nu, 1.0e-30_dp)

    alpha_nu = j_nu / b_nu
  end subroutine compute_synchrotron_absorptivity

end module umbra_synchrotron
