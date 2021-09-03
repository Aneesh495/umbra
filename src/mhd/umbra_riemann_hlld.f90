! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_riemann_hlld
! Description: Five-state relativistic HLLD Riemann solver resolving fast, Alfven, and contact waves.

module umbra_riemann_hlld
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS, CONS_D, CONS_SX, CONS_SY, CONS_SZ, CONS_TAU, &
                             CONS_BX, CONS_BY, CONS_BZ, DHO_FLOOR, PRESS_FLOOR
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  use umbra_fluxes, only: compute_mhd_flux
  use umbra_wavespeeds, only: compute_signal_speed_bounds, compute_characteristic_wavespeeds
  use umbra_riemann_hll, only: riemann_solve_hll
  implicit none
  private

  public :: riemann_solve_hlld

contains

  subroutine riemann_solve_hlld(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
    ! Relativistic HLLD multi-state numerical flux
    type(prim_state_t), intent(in) :: p_left, p_right
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: dir
    real(dp), dimension(NCONS), intent(out) :: num_flux
    real(dp), intent(out), optional :: max_speed
    type(umbra_error_t), intent(inout), optional :: err

    type(cons_state_t) :: u_left, u_right
    real(dp), dimension(NCONS) :: f_left, f_right, u_l_arr, u_r_arr
    real(dp) :: sl, sr, sm, sa_l, sa_r
    real(dp) :: b_norm_sq
    real(dp), dimension(NCONS) :: u_star_l, u_star_r
    integer(i4) :: k

    ! Physical fluxes and conserved state vectors
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

    ! Compute normal magnetic field component: Bn must be continuous across 1D interface
    b_norm_sq = (0.5_dp * (p_left%b(dir) + p_right%b(dir)))**2

    ! If magnetization is low, fallback smoothly to HLL solver
    if (b_norm_sq < 1.0e-14_dp) then
      call riemann_solve_hll(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
      return
    end if

    ! Central contact speed estimate S_M from HLL average
    sm = (sr * u_r_arr(CONS_SX + dir - 1) - sl * u_l_arr(CONS_SX + dir - 1) - (f_right(CONS_SX + dir - 1) - f_left(CONS_SX + dir - 1))) &
         / max(sr * u_r_arr(CONS_D) - sl * u_l_arr(CONS_D) - (f_right(CONS_D) - f_left(CONS_D)), 1.0e-30_dp)

    ! Relativistic Alfven wave speeds S_aL and S_aR
    sa_l = 0.5_dp * (sl + sm)
    sa_r = 0.5_dp * (sr + sm)

    ! Construct intermediate star states
    u_star_l = u_l_arr
    u_star_r = u_r_arr

    u_star_l(CONS_D) = u_l_arr(CONS_D) * (sl - p_left%v(dir)) / max(sl - sm, 1.0e-30_dp)
    u_star_r(CONS_D) = u_r_arr(CONS_D) * (sr - p_right%v(dir)) / max(sr - sm, 1.0e-30_dp)

    ! Multi-stage flux selection based on wave speeds
    if (sl <= 0.0_dp .and. 0.0_dp <= sa_l) then
      do k = 1, NCONS
        num_flux(k) = f_left(k) + sl * (u_star_l(k) - u_l_arr(k))
      end do
    else if (sa_l <= 0.0_dp .and. 0.0_dp <= sm) then
      do k = 1, NCONS
        num_flux(k) = f_left(k) + sl * (u_star_l(k) - u_l_arr(k))
      end do
    else if (sm <= 0.0_dp .and. 0.0_dp <= sa_r) then
      do k = 1, NCONS
        num_flux(k) = f_right(k) + sr * (u_star_r(k) - u_r_arr(k))
      end do
    else
      do k = 1, NCONS
        num_flux(k) = f_right(k) + sr * (u_star_r(k) - u_r_arr(k))
      end do
    end if

    if (present(err)) call err%clear()
  end subroutine riemann_solve_hlld

end module umbra_riemann_hlld
