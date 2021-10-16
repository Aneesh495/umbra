! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_bench_briowu
! Description: One-dimensional relativistic magnetized Brio-Wu shock tube benchmark (Balsara 2001).

module umbra_bench_briowu
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_P, &
                             PRIM_BX, PRIM_BY, PRIM_BZ
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_grid1d, only: grid1d_t
  implicit none
  private

  public :: init_bench_briowu_1d

contains

  subroutine init_bench_briowu_1d(grid, x_discontinuity, err)
    ! Initializes relativistic Brio-Wu magnetized shock tube:
    ! Left (x < x_disc):  rho = 1.0,   v = 0, p = 1.0, Bx = 0.5, By =  1.0, Bz = 0
    ! Right (x >= x_disc): rho = 0.125, v = 0, p = 0.1, Bx = 0.5, By = -1.0, Bz = 0
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
        grid%prim(PRIM_RHO, i) = 1.0_dp
        grid%prim(PRIM_VX,  i) = 0.0_dp
        grid%prim(PRIM_VY,  i) = 0.0_dp
        grid%prim(PRIM_VZ,  i) = 0.0_dp
        grid%prim(PRIM_P,   i) = 1.0_dp
        grid%prim(PRIM_BX,  i) = 0.5_dp
        grid%prim(PRIM_BY,  i) = 1.0_dp
        grid%prim(PRIM_BZ,  i) = 0.0_dp
      else
        grid%prim(PRIM_RHO, i) = 0.125_dp
        grid%prim(PRIM_VX,  i) = 0.0_dp
        grid%prim(PRIM_VY,  i) = 0.0_dp
        grid%prim(PRIM_VZ,  i) = 0.0_dp
        grid%prim(PRIM_P,   i) = 0.1_dp
        grid%prim(PRIM_BX,  i) = 0.5_dp
        grid%prim(PRIM_BY,  i) = -1.0_dp
        grid%prim(PRIM_BZ,  i) = 0.0_dp
      end if
    end do

    if (present(err)) call err%clear()
  end subroutine init_bench_briowu_1d

end module umbra_bench_briowu
