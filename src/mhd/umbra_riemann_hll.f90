! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_riemann_hll
! Description: Harten-Lax-van Leer (HLL) two-wave relativistic numerical flux solver.

module umbra_riemann_hll
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  use umbra_fluxes, only: compute_mhd_flux
  use umbra_wavespeeds, only: compute_signal_speed_bounds
  implicit none
  private

  public :: riemann_solve_hll

contains

  subroutine riemann_solve_hll(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
    type(prim_state_t), intent(in) :: p_left, p_right
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: dir
    real(dp), dimension(NCONS), intent(out) :: num_flux
    real(dp), intent(out), optional :: max_speed
    type(umbra_error_t), intent(inout), optional :: err

    type(cons_state_t) :: u_left, u_right
    real(dp), dimension(NCONS) :: f_left, f_right, u_l_arr, u_r_arr
    real(dp) :: sl, sr, inv_diff
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
    else if (sr <= 0.0_dp) then
      num_flux = f_right
    else
      inv_diff = 1.0_dp / (sr - sl)
      do k = 1, NCONS
        num_flux(k) = (sr * f_left(k) - sl * f_right(k) + sl * sr * (u_r_arr(k) - u_l_arr(k))) * inv_diff
      end do
    end if

    if (present(err)) call err%clear()
  end subroutine riemann_solve_hll

end module umbra_riemann_hll
