! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_bench_sod
! Description: One-dimensional relativistic Sod shock tube benchmark (Marti & Mueller 1994).

module umbra_bench_sod
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_P, &
                             PRIM_BX, PRIM_BY, PRIM_BZ
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_grid1d, only: grid1d_t
  implicit none
  private

  public :: init_bench_sod_1d

contains

  subroutine init_bench_sod_1d(grid, x_discontinuity, err)
    ! Initializes relativistic Sod shock tube problem:
    ! Left (x < x_disc):  rho = 10.0, v = 0.0, p = 13.333333333333334
    ! Right (x >= x_disc): rho = 1.0,  v = 0.0, p = 1.0e-6
    type(grid1d_t), intent(inout) :: grid
    real(dp), intent(in), optional :: x_discontinuity
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: x_disc
    integer(i4) :: i

    x_disc = 0.5_dp
    if (present(x_discontinuity)) x_disc = x_discontinuity

    do i = 1 - grid%nghost, grid%nx + grid%nghost
      grid%prim(:, i) = 0.0_dp
      if (grid%x_cell(i) < x_disc) then
        grid%prim(PRIM_RHO, i) = 10.0_dp
        grid%prim(PRIM_VX,  i) = 0.0_dp
        grid%prim(PRIM_VY,  i) = 0.0_dp
        grid%prim(PRIM_VZ,  i) = 0.0_dp
        grid%prim(PRIM_P,   i) = 40.0_dp / 3.0_dp
      else
        grid%prim(PRIM_RHO, i) = 1.0_dp
        grid%prim(PRIM_VX,  i) = 0.0_dp
        grid%prim(PRIM_VY,  i) = 0.0_dp
        grid%prim(PRIM_VZ,  i) = 0.0_dp
        grid%prim(PRIM_P,   i) = 1.0e-6_dp
      end if
    end do

    if (present(err)) call err%clear()
  end subroutine init_bench_sod_1d

end module umbra_bench_sod
