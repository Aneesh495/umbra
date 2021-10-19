! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_bench_kelvin_helmholtz
! Description: Two-dimensional relativistic Kelvin-Helmholtz shear layer instability benchmark.

module umbra_bench_kelvin_helmholtz
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_P, &
                             PRIM_BX, PRIM_BY, PRIM_BZ, TWO_PI
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_grid2d, only: grid2d_t
  implicit none
  private

  public :: init_bench_kelvin_helmholtz_2d

contains

  subroutine init_bench_kelvin_helmholtz_2d(grid, err)
    ! Relativistic Kelvin-Helmholtz instability on domain [-0.5, 0.5] x [-0.5, 0.5]:
    ! Upper and lower fluid ( |y| > 0.25 ): rho = 1.0, vx = -0.5, p = 1.0
    ! Central shear layer  ( |y| <= 0.25 ): rho = 2.0, vx =  0.5, p = 1.0
    ! Transverse velocity perturbation: vy = 0.01 * sin(2 pi x) * exp(- (y / 0.05)^2)
    ! Uniform longitudinal magnetic field Bx = 0.1
    type(grid2d_t), intent(inout) :: grid
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: x, y, vy_pert
    integer(i4) :: i, j

    do j = 1 - grid%nghost, grid%ny + grid%nghost
      y = grid%y_cell(j)
      do i = 1 - grid%nghost, grid%nx + grid%nghost
        x = grid%x_cell(i)

        vy_pert = 0.01_dp * sin(TWO_PI * x) * exp(-((abs(y) - 0.25_dp) / 0.05_dp)**2)

        grid%prim(:, i, j) = 0.0_dp
        if (abs(y) <= 0.25_dp) then
          grid%prim(PRIM_RHO, i, j) = 2.0_dp
          grid%prim(PRIM_VX,  i, j) = 0.5_dp
        else
          grid%prim(PRIM_RHO, i, j) = 1.0_dp
          grid%prim(PRIM_VX,  i, j) = -0.5_dp
        end if

        grid%prim(PRIM_VY, i, j) = vy_pert
        grid%prim(PRIM_VZ, i, j) = 0.0_dp
        grid%prim(PRIM_P,  i, j) = 1.0_dp
        grid%prim(PRIM_BX, i, j) = 0.1_dp
        grid%prim(PRIM_BY, i, j) = 0.0_dp
        grid%prim(PRIM_BZ, i, j) = 0.0_dp
      end do
    end do

    do j = 1, grid%ny
      do i = 1, grid%nx + 1
        grid%bx_face(i, j) = 0.1_dp
      end do
    end do
    grid%by_face = 0.0_dp

    if (present(err)) call err%clear()
  end subroutine init_bench_kelvin_helmholtz_2d

end module umbra_bench_kelvin_helmholtz
