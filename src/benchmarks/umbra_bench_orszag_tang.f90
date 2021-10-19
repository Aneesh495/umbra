! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_bench_orszag_tang
! Description: Two-dimensional relativistic Orszag-Tang vortex benchmark for magnetized turbulence.

module umbra_bench_orszag_tang
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_P, &
                             PRIM_BX, PRIM_BY, PRIM_BZ, TWO_PI, FOUR_PI
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_grid2d, only: grid2d_t
  implicit none
  private

  public :: init_bench_orszag_tang_2d

contains

  subroutine init_bench_orszag_tang_2d(grid, err)
    ! Relativistic Orszag-Tang vortex configuration on periodic domain [0, 1] x [0, 1]:
    ! rho = 25 / (36 pi), p = 5 / (12 pi)
    ! vx = -v0 * sin(2 pi y), vy = v0 * sin(2 pi x) with v0 = 0.5 / sqrt(2)
    ! Bx = -B0 * sin(2 pi y), By = B0 * sin(4 pi x) with B0 = 1 / sqrt(4 pi)
    type(grid2d_t), intent(inout) :: grid
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: v0, b0, x, y
    integer(i4) :: i, j

    v0 = 0.5_dp / sqrt(2.0_dp)
    b0 = 1.0_dp / sqrt(FOUR_PI)

    do j = 1 - grid%nghost, grid%ny + grid%nghost
      y = grid%y_cell(j)
      do i = 1 - grid%nghost, grid%nx + grid%nghost
        x = grid%x_cell(i)

        grid%prim(:, i, j) = 0.0_dp
        grid%prim(PRIM_RHO, i, j) = 25.0_dp / (36.0_dp * (TWO_PI / 2.0_dp))
        grid%prim(PRIM_P,   i, j) = 5.0_dp / (12.0_dp * (TWO_PI / 2.0_dp))
        grid%prim(PRIM_VX,  i, j) = -v0 * sin(TWO_PI * y)
        grid%prim(PRIM_VY,  i, j) =  v0 * sin(TWO_PI * x)
        grid%prim(PRIM_VZ,  i, j) =  0.0_dp
        grid%prim(PRIM_BX,  i, j) = -b0 * sin(TWO_PI * y)
        grid%prim(PRIM_BY,  i, j) =  b0 * sin(2.0_dp * TWO_PI * x)
        grid%prim(PRIM_BZ,  i, j) =  0.0_dp
      end do
    end do

    ! Set face-centered staggered fields for Constrained Transport
    do j = 1, grid%ny
      y = grid%y_cell(j)
      do i = 1, grid%nx + 1
        grid%bx_face(i, j) = -b0 * sin(TWO_PI * y)
      end do
    end do

    do j = 1, grid%ny + 1
      do i = 1, grid%nx
        x = grid%x_cell(i)
        grid%by_face(i, j) = b0 * sin(2.0_dp * TWO_PI * x)
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine init_bench_orszag_tang_2d

end module umbra_bench_orszag_tang
