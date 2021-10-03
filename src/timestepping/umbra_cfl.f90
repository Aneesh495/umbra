! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_cfl
! Description: Multi-dimensional Courant-Friedrichs-Lewy (CFL) time step calculation for relativistic flows.

module umbra_cfl
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t
  use umbra_wavespeeds, only: compute_characteristic_wavespeeds
  implicit none
  private

  public :: compute_cfl_timestep_2d

contains

  subroutine compute_cfl_timestep_2d(prim, nx, ny, nghost, dx, dy, m, eos, cfl_target, dt, err)
    ! Computes maximum stable timestep dt based on 2D relativistic spectral wave bounds:
    ! dt = CFL / max_ij [ (|lambda_x| / dx) + (|lambda_y| / dy) ]
    integer(i4), intent(in) :: nx, ny, nghost
    real(dp), dimension(NPRIM, 1-nghost : nx+nghost, 1-nghost : ny+nghost), intent(in) :: prim
    real(dp), intent(in) :: dx, dy, cfl_target
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    real(dp), intent(out) :: dt
    type(umbra_error_t), intent(inout), optional :: err

    type(prim_state_t) :: p
    real(dp) :: slx, srx, sly, sry, max_rate, cell_rate
    real(dp) :: inv_dx, inv_dy
    integer(i4) :: i, j

    inv_dx = 1.0_dp / dx
    inv_dy = 1.0_dp / dy
    max_rate = 1.0e-8_dp

    do j = 1, ny
      do i = 1, nx
        call p%from_array(prim(:, i, j))

        call compute_characteristic_wavespeeds(p, m, eos, 1_i4, slx, srx)
        call compute_characteristic_wavespeeds(p, m, eos, 2_i4, sly, sry)

        cell_rate = max(abs(slx), abs(srx)) * inv_dx + max(abs(sly), abs(sry)) * inv_dy
        max_rate = max(max_rate, cell_rate)
      end do
    end do

    dt = cfl_target / max_rate

    if (present(err)) call err%clear()
  end subroutine compute_cfl_timestep_2d

end module umbra_cfl
