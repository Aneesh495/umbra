! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_bench_blast
! Description: Two-dimensional cylindrical magnetized relativistic blast wave benchmark (Komissarov 1999).

module umbra_bench_blast
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_P, &
                             PRIM_BX, PRIM_BY, PRIM_BZ
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_grid2d, only: grid2d_t
  implicit none
  private

  public :: init_bench_blast_2d

contains

  subroutine init_bench_blast_2d(grid, r_blast, p_in, p_out, bx_val, err)
    ! Cylindrical relativistic blast wave in uniform magnetic field B = (bx_val, 0, 0):
    ! Inside cylinder (r <= r_blast): rho = 1.0e-2, p = 1.0 (or p_in), v = 0
    ! Outside cylinder (r > r_blast): rho = 1.0e-4, p = 1.0e-3 (or p_out), v = 0
    type(grid2d_t), intent(inout) :: grid
    real(dp), intent(in), optional :: r_blast, p_in, p_out, bx_val
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: rb, pin, pout, bx0, x, y, r
    integer(i4) :: i, j

    rb = 0.8_dp
    pin = 1.0_dp
    pout = 1.0e-3_dp
    bx0 = 0.5_dp

    if (present(r_blast)) rb = r_blast
    if (present(p_in)) pin = p_in
    if (present(p_out)) pout = p_out
    if (present(bx_val)) bx0 = bx_val

    do j = 1 - grid%nghost, grid%ny + grid%nghost
      y = grid%y_cell(j)
      do i = 1 - grid%nghost, grid%nx + grid%nghost
        x = grid%x_cell(i)
        r = sqrt(x * x + y * y)

        grid%prim(:, i, j) = 0.0_dp
        if (r <= rb) then
          grid%prim(PRIM_RHO, i, j) = 1.0e-2_dp
          grid%prim(PRIM_P,   i, j) = pin
        else
          grid%prim(PRIM_RHO, i, j) = 1.0e-4_dp
          grid%prim(PRIM_P,   i, j) = pout
        end if

        grid%prim(PRIM_VX, i, j) = 0.0_dp
        grid%prim(PRIM_VY, i, j) = 0.0_dp
        grid%prim(PRIM_VZ, i, j) = 0.0_dp
        grid%prim(PRIM_BX, i, j) = bx0
        grid%prim(PRIM_BY, i, j) = 0.0_dp
        grid%prim(PRIM_BZ, i, j) = 0.0_dp
      end do
    end do

    do j = 1, grid%ny
      do i = 1, grid%nx + 1
        grid%bx_face(i, j) = bx0
      end do
    end do
    grid%by_face = 0.0_dp

    if (present(err)) call err%clear()
  end subroutine init_bench_blast_2d

end module umbra_bench_blast
