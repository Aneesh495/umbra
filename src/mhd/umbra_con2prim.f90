! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_con2prim
! Description: Unified conservative-to-primitive dispatcher with multi-stage fallbacks and physical floors.

module umbra_con2prim
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: DHO_FLOOR, PRESS_FLOOR, GAMMA_MAX, VEL_MAX, BSQ_OVER_RHO_MAX
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_CON2PRIM_FAILED
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t
  use umbra_con2prim_1d, only: con2prim_1d_solve
  use umbra_con2prim_2d, only: con2prim_2d_solve
  implicit none
  private

  public :: con2prim_solve
  public :: apply_physical_floors

contains

  subroutine con2prim_solve(u, m, eos, p, err)
    type(cons_state_t), intent(in) :: u
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    type(prim_state_t), intent(out) :: p
    type(umbra_error_t), intent(inout), optional :: err

    type(umbra_error_t) :: local_err
    logical :: success

    success = .false.

    ! Stage 1: Fast 1D Noble-Mignone solver
    call con2prim_1d_solve(u, m, eos, p, local_err)
    if (local_err%ok()) then
      if (p%rho > 0.0_dp .and. p%p > 0.0_dp) success = .true.
    end if

    ! Stage 2: Robust 2D Newton-Raphson fallback
    if (.not. success) then
      call local_err%clear()
      call con2prim_2d_solve(u, m, eos, p, local_err)
      if (local_err%ok()) then
        if (p%rho > 0.0_dp .and. p%p > 0.0_dp) success = .true.
      end if
    end if

    ! Stage 3: Cold atmospheric floor recovery if extreme shock or horizon breach
    if (.not. success) then
      p%rho = max(u%d, DHO_FLOOR)
      p%p = PRESS_FLOOR
      p%v = 0.0_dp
      p%b = u%b
    end if

    ! Enforce physical boundary floors and velocity caps
    call apply_physical_floors(p, m)

    if (present(err)) call err%clear()
  end subroutine con2prim_solve

  pure subroutine apply_physical_floors(p, m)
    type(prim_state_t), intent(inout) :: p
    type(metric_3p1_t), intent(in) :: m

    real(dp) :: v_sq, b2, sigma_m
    integer(i4) :: i, j

    ! Floor density and pressure
    p%rho = max(p%rho, DHO_FLOOR)
    p%p = max(p%p, PRESS_FLOOR)

    ! Velocity cap: ensure v^2 <= VEL_MAX^2
    v_sq = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        v_sq = v_sq + m%gamma_down(i, j) * p%v(i) * p%v(j)
      end do
    end do

    if (v_sq > VEL_MAX * VEL_MAX) then
      p%v = p%v * (VEL_MAX / sqrt(max(v_sq, 1.0e-30_dp)))
    end if

    ! Magnetization ceiling: b^2 / rho <= BSQ_OVER_RHO_MAX
    b2 = p%magnetic_energy(m)
    sigma_m = b2 / max(p%rho, 1.0e-30_dp)
    if (sigma_m > BSQ_OVER_RHO_MAX) then
      p%rho = b2 / BSQ_OVER_RHO_MAX
    end if
  end subroutine apply_physical_floors

end module umbra_con2prim
