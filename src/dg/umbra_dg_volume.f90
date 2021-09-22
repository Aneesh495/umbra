! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_dg_volume
! Description: Weak and strong form Discontinuous Galerkin volume integral operators.

module umbra_dg_volume
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_dg_basis, only: dg_basis_t
  implicit none
  private

  public :: compute_dg_volume_operator_1d
  public :: compute_dg_volume_operator_2d

contains

  subroutine compute_dg_volume_operator_1d(basis, f_nodal, dx, vol_residual, err)
    ! Evaluates strong-form volume derivative term: vol_residual(:, i) = (2 / dx) * sum_j D_{ij} F_j
    type(dg_basis_t), intent(in) :: basis
    real(dp), dimension(NCONS, basis%num_nodes), intent(in) :: f_nodal
    real(dp), intent(in) :: dx
    real(dp), dimension(NCONS, basis%num_nodes), intent(out) :: vol_residual
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: inv_jacobian
    integer(i4) :: i, j, v, n

    n = basis%num_nodes
    inv_jacobian = 2.0_dp / dx
    vol_residual = 0.0_dp

    do i = 1, n
      do j = 1, n
        do v = 1, NCONS
          vol_residual(v, i) = vol_residual(v, i) + basis%d_mat(i, j) * f_nodal(v, j)
        end do
      end do
      vol_residual(:, i) = vol_residual(:, i) * inv_jacobian
    end do

    if (present(err)) call err%clear()
  end subroutine compute_dg_volume_operator_1d

  subroutine compute_dg_volume_operator_2d(basis, fx_nodal, fy_nodal, dx, dy, vol_residual, err)
    ! Evaluates 2D tensor-product volume derivative: (2/dx) dFx/dxi + (2/dy) dFy/deta
    type(dg_basis_t), intent(in) :: basis
    real(dp), dimension(NCONS, basis%num_nodes, basis%num_nodes), intent(in) :: fx_nodal
    real(dp), dimension(NCONS, basis%num_nodes, basis%num_nodes), intent(in) :: fy_nodal
    real(dp), intent(in) :: dx, dy
    real(dp), dimension(NCONS, basis%num_nodes, basis%num_nodes), intent(out) :: vol_residual
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: inv_jx, inv_jy
    integer(i4) :: i, j, k, v, n

    n = basis%num_nodes
    inv_jx = 2.0_dp / dx
    inv_jy = 2.0_dp / dy
    vol_residual = 0.0_dp

    do j = 1, n
      do i = 1, n
        ! Derivative along xi (x-direction): sum_k D_{ik} Fx_{kj}
        do k = 1, n
          do v = 1, NCONS
            vol_residual(v, i, j) = vol_residual(v, i, j) + inv_jx * basis%d_mat(i, k) * fx_nodal(v, k, j)
          end do
        end do
        ! Derivative along eta (y-direction): sum_k D_{jk} Fy_{ik}
        do k = 1, n
          do v = 1, NCONS
            vol_residual(v, i, j) = vol_residual(v, i, j) + inv_jy * basis%d_mat(j, k) * fy_nodal(v, i, k)
          end do
        end do
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine compute_dg_volume_operator_2d

end module umbra_dg_volume
