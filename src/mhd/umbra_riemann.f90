! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_riemann
! Description: Unified polymorphic dispatcher for relativistic Riemann solvers.

module umbra_riemann
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t
  use umbra_riemann_llf, only: riemann_solve_llf
  use umbra_riemann_hll, only: riemann_solve_hll
  use umbra_riemann_hllc, only: riemann_solve_hllc
  use umbra_riemann_hlld, only: riemann_solve_hlld
  implicit none
  private

  integer(i4), parameter, public :: RIEMANN_LLF  = 1_i4
  integer(i4), parameter, public :: RIEMANN_HLL  = 2_i4
  integer(i4), parameter, public :: RIEMANN_HLLC = 3_i4
  integer(i4), parameter, public :: RIEMANN_HLLD = 4_i4

  public :: riemann_solve

contains

  subroutine riemann_solve(solver_type, p_left, p_right, m, eos, dir, num_flux, max_speed, err)
    integer(i4), intent(in) :: solver_type
    type(prim_state_t), intent(in) :: p_left, p_right
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: dir
    real(dp), dimension(NCONS), intent(out) :: num_flux
    real(dp), intent(out), optional :: max_speed
    type(umbra_error_t), intent(inout), optional :: err

    select case (solver_type)
      case (RIEMANN_LLF)
        call riemann_solve_llf(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
      case (RIEMANN_HLL)
        call riemann_solve_hll(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
      case (RIEMANN_HLLC)
        call riemann_solve_hllc(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
      case (RIEMANN_HLLD)
        call riemann_solve_hlld(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
      case default
        ! Default fallback is HLLD
        call riemann_solve_hlld(p_left, p_right, m, eos, dir, num_flux, max_speed, err)
    end select
  end subroutine riemann_solve

end module umbra_riemann
