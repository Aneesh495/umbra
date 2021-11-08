!> Unit tests for Discontinuous Galerkin basis and differentiation operators
module test_dg_operators
  use umbra_kinds, only: dp, i4
  use umbra_dg_basis, only: dg_basis_t
  use test_framework, only: test_suite_t, assert_true, assert_close_dp
  implicit none
  private

  public :: run_test_dg_operators

contains

  subroutine run_test_dg_operators(suite)
    type(test_suite_t), intent(inout) :: suite
    type(dg_basis_t) :: basis
    real(dp), allocatable :: u_nodes(:), du_nodes(:)
    integer(i4) :: i, j, n

    ! Test 1: DG basis initialization at degree N=3 (4 nodes)
    call basis%init(3_i4, use_lobatto=.true.)
    call assert_true(suite, basis%num_nodes == 4_i4, "DG degree 3 has 4 nodal points")
    call assert_close_dp(suite, basis%b_left(1), 1.0_dp, 1.0e-12_dp, "DG LGL left boundary evaluation")
    call assert_close_dp(suite, basis%b_right(4), 1.0_dp, 1.0e-12_dp, "DG LGL right boundary evaluation")

    ! Test 2: Differentiation of constant field u(xi) = 3.0 gives exactly du/dxi = 0.0
    n = basis%num_nodes
    allocate(u_nodes(n), du_nodes(n))
    u_nodes = 3.0_dp
    du_nodes = 0.0_dp

    do i = 1, n
      do j = 1, n
        du_nodes(i) = du_nodes(i) + basis%d_mat(i, j) * u_nodes(j)
      end do
      call assert_close_dp(suite, du_nodes(i), 0.0_dp, 1.0e-12_dp, "DG constant differentiation = 0")
    end do

    ! Test 3: Differentiation of linear field u(xi) = 2.0 * xi gives exactly du/dxi = 2.0
    do i = 1, n
      u_nodes(i) = 2.0_dp * basis%nodes(i)
    end do
    du_nodes = 0.0_dp

    do i = 1, n
      do j = 1, n
        du_nodes(i) = du_nodes(i) + basis%d_mat(i, j) * u_nodes(j)
      end do
      call assert_close_dp(suite, du_nodes(i), 2.0_dp, 1.0e-12_dp, "DG linear differentiation = 2.0")
    end do

    call basis%free()
    deallocate(u_nodes, du_nodes)
  end subroutine run_test_dg_operators

end module test_dg_operators
