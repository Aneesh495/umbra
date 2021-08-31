! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_riemann_hllc
! Description: Relativistic HLLC Riemann solver restoring the contact wave discontinuity.

module umbra_riemann_hllc
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS, CONS_D, CONS_SX, CONS_SY, CONS_SZ, CONS_TAU, &
                             CONS_BX, CONS_BY, CONS_BZ, DHO_FLOOR, PRESS_FLOOR
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  use umbra_fluxes, only: compute_mhd_flux
  use umbra_wavespeeds, only: compute_signal_speed_bounds
  implicit none
  private

  public :: riemann_solve_hllc

contains

  subroutine riemann_solve_hllc(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
    type(prim_state_t), intent(in) :: p_left, p_right
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: dir
    real(dp), dimension(NCONS), intent(out) :: num_flux
    real(dp), intent(out), optional :: max_speed
    type(umbra_error_t), intent(inout), optional :: err

    type(cons_state_t) :: u_left, u_right
    real(dp), dimension(NCONS) :: f_left, f_right, u_star, u_l_arr, u_r_arr
    real(dp) :: sl, sr, s_star, p_star
    real(dp) :: gamma_l, gamma_r, eps_l, eps_r, h_l, h_r, rho_h_l, rho_h_r
    real(dp) :: a_l, a_r, vn_l, vn_r, p_tot_l, p_tot_r
    integer(i4) :: k

    call compute_mhd_flux(p_left, m, eos, dir, f_left, err)
    call compute_mhd_flux(p_right, m, eos, dir, f_right, err)

    call prim_to_cons(p_left, m, eos, u_left, err)
    call prim_to_cons(p_right, m, eos, u_right, err)

    call compute_signal_speed_bounds(p_left, p_right, m, eos, dir, sl, sr)
    if (present(max_speed)) max_speed = max(abs(sl), abs(sr))

    u_l_arr = u_left%to_array()
    u_r_arr = u_right%to_array()

    if (sl >= 0.0_dp) then
      num_flux = f_left
      return
    else if (sr <= 0.0_dp) then
      num_flux = f_right
      return
    end if

    gamma_l = p_left%lorentz_factor(m)
    gamma_r = p_right%lorentz_factor(m)
    eps_l = eos%specific_internal_energy(p_left%rho, p_left%p)
    eps_r = eos%specific_internal_energy(p_right%rho, p_right%p)
    h_l = eos%specific_enthalpy(p_left%rho, eps_l, p_left%p)
    h_r = eos%specific_enthalpy(p_right%rho, eps_r, p_right%p)

    rho_h_l = p_left%rho * h_l * (gamma_l**2)
    rho_h_r = p_right%rho * h_r * (gamma_r**2)

    vn_l = p_left%v(dir)
    vn_r = p_right%v(dir)
    p_tot_l = p_left%total_pressure(m)
    p_tot_r = p_right%total_pressure(m)

    a_l = (sl - vn_l) * rho_h_l
    a_r = (sr - vn_r) * rho_h_r

    ! Contact wave speed S_*
    s_star = (a_r * vn_r - a_l * vn_l + p_tot_l - p_tot_r) / max(a_r - a_l, 1.0e-30_dp)
    p_star = p_tot_l + a_l * (s_star - vn_l)

    if (s_star >= 0.0_dp) then
      ! Left star state U*_L
      u_star = u_l_arr
      u_star(CONS_D) = u_l_arr(CONS_D) * (sl - vn_l) / max(sl - s_star, 1.0e-30_dp)
      u_star(CONS_SX + dir - 1) = (u_l_arr(CONS_SX + dir - 1) * (sl - vn_l) + p_star - p_tot_l) / &
                                  max(sl - s_star, 1.0e-30_dp)
      u_star(CONS_TAU) = (u_l_arr(CONS_TAU) * (sl - vn_l) + p_star * s_star - p_tot_l * vn_l) / &
                         max(sl - s_star, 1.0e-30_dp)
      do k = 1, NCONS
        num_flux(k) = f_left(k) + sl * (u_star(k) - u_l_arr(k))
      end do
    else
      ! Right star state U*_R
      u_star = u_r_arr
      u_star(CONS_D) = u_r_arr(CONS_D) * (sr - vn_r) / max(sr - s_star, 1.0e-30_dp)
      u_star(CONS_SX + dir - 1) = (u_r_arr(CONS_SX + dir - 1) * (sr - vn_r) + p_star - p_tot_r) / &
                                  max(sr - s_star, 1.0e-30_dp)
      u_star(CONS_TAU) = (u_r_arr(CONS_TAU) * (sr - vn_r) + p_star * s_star - p_tot_r * vn_r) / &
                         max(sr - s_star, 1.0e-30_dp)
      do k = 1, NCONS
        num_flux(k) = f_right(k) + sr * (u_star(k) - u_r_arr(k))
      end do
    end if

    if (present(err)) call err%clear()
  end subroutine riemann_solve_hllc

end module umbra_riemann_hllc
