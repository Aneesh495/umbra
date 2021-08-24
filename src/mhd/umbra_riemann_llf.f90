! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_riemann_llf
! Description: Local Lax-Friedrichs (Rusanov) relativistic numerical flux solver.

module umbra_riemann_llf
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  use umbra_fluxes, only: compute_mhd_flux
  use umbra_wavespeeds, only: compute_characteristic_wavespeeds
  implicit none
  private

  public :: riemann_solve_llf

contains

  subroutine riemann_solve_llf(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
    type(prim_state_t), intent(in) :: p_left, p_right
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: dir
    real(dp), dimension(NCONS), intent(out) :: num_flux
    real(dp), intent(out), optional :: max_speed
    type(umbra_error_t), intent(inout), optional :: err

    type(cons_state_t) :: u_left, u_right
    real(dp), dimension(NCONS) :: f_left, f_right, u_diff
    real(dp) :: sl_left, sr_left, sl_right, sr_right, c_max
    integer(i4) :: k

    ! Physical flux vectors
    call compute_mhd_flux(p_left, m, eos, dir, f_left, err)
    call compute_mhd_flux(p_right, m, eos, dir, f_right, err)

    ! Conserved states
    call prim_to_cons(p_left, m, eos, u_left, err)
    call prim_to_cons(p_right, m, eos, u_right, err)

    ! Wave speeds
    call compute_characteristic_wavespeeds(p_left, m, eos, dir, sl_left, sr_left)
    call compute_characteristic_wavespeeds(p_right, m, eos, dir, sl_right, sr_right)

    c_max = max(max(abs(sl_left), abs(sr_left)), max(abs(sl_right), abs(sr_right)))
    if (present(max_speed)) max_speed = c_max

    u_diff = u_right%to_array() - u_left%to_array()

    ! Rusanov dissipation formula: F_num = 0.5 * (F_L + F_R) - 0.5 * c_max * (U_R - U_L)
    do k = 1, NCONS
      num_flux(k) = 0.5_dp * (f_left(k) + f_right(k)) - 0.5_dp * c_max * u_diff(k)
    end do

    if (present(err)) call err%clear()
  end subroutine riemann_solve_llf

end module umbra_riemann_llf
