! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_dg_surface
! Description: Discontinuous Galerkin surface flux integration across element interfaces.

module umbra_dg_surface
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_dg_basis, only: dg_basis_t
  implicit none
  private

  public :: compute_dg_surface_operator_1d
  public :: interpolate_to_element_boundaries_1d

contains

  pure subroutine interpolate_to_element_boundaries_1d(basis, u_nodal, u_left_boundary, u_right_boundary)
    ! Evaluates the polynomial solution at the element physical boundaries xi = -1 and xi = +1
    type(dg_basis_t), intent(in) :: basis
    real(dp), dimension(NCONS, basis%num_nodes), intent(in) :: u_nodal
    real(dp), dimension(NCONS), intent(out) :: u_left_boundary
    real(dp), dimension(NCONS), intent(out) :: u_right_boundary

    integer(i4) :: j, v, n

    n = basis%num_nodes
    u_left_boundary = 0.0_dp
    u_right_boundary = 0.0_dp

    do j = 1, n
      do v = 1, NCONS
        u_left_boundary(v)  = u_left_boundary(v)  + basis%b_left(j)  * u_nodal(v, j)
        u_right_boundary(v) = u_right_boundary(v) + basis%b_right(j) * u_nodal(v, j)
      end do
    end do
  end subroutine interpolate_to_element_boundaries_1d

  subroutine compute_dg_surface_operator_1d(basis, f_num_left, f_int_left, f_num_right, f_int_right, &
                                           dx, surf_residual, err)
    ! Evaluates the surface penalty:
    ! surf_residual(:, i) = (2 / (dx * w_i)) * [ b_right(i) * (f_num_right - f_int_right) &
    !                                          - b_left(i)  * (f_num_left  - f_int_left) ]
    type(dg_basis_t), intent(in) :: basis
    real(dp), dimension(NCONS), intent(in) :: f_num_left, f_int_left
    real(dp), dimension(NCONS), intent(in) :: f_num_right, f_int_right
    real(dp), intent(in) :: dx
    real(dp), dimension(NCONS, basis%num_nodes), intent(out) :: surf_residual
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: inv_dx, factor
    real(dp), dimension(NCONS) :: jump_left, jump_right
    integer(i4) :: i, v, n

    n = basis%num_nodes
    inv_dx = 2.0_dp / dx
    jump_left  = f_num_left  - f_int_left
    jump_right = f_num_right - f_int_right
    surf_residual = 0.0_dp

    do i = 1, n
      factor = inv_dx / max(basis%weights(i), 1.0e-30_dp)
      do v = 1, NCONS
        surf_residual(v, i) = factor * (basis%b_right(i) * jump_right(v) - basis%b_left(i) * jump_left(v))
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine compute_dg_surface_operator_1d

end module umbra_dg_surface
