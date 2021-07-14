! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_polynomials
! Description: Legendre polynomials, nodal Vandermonde matrices, and spectral differentiation operators.

module umbra_polynomials
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: EPSILON_DP
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_QUADRATURE_ORDER
  implicit none
  private

  public :: legendre_poly_and_deriv
  public :: lagrange_interpolant
  public :: nodal_differentiation_matrix
  public :: modal_to_nodal_vandermonde

contains

  pure subroutine legendre_poly_and_deriv(n, x, p, dp_dx)
    ! Evaluates Legendre polynomial P_n(x) and its derivative P'_n(x)
    integer(i4), intent(in) :: n
    real(dp), intent(in) :: x
    real(dp), intent(out) :: p
    real(dp), intent(out) :: dp_dx

    integer(i4) :: k
    real(dp) :: p_prev, p_curr, p_next

    if (n == 0) then
      p = 1.0_dp
      dp_dx = 0.0_dp
      return
    else if (n == 1) then
      p = x
      dp_dx = 1.0_dp
      return
    end if

    p_prev = 1.0_dp
    p_curr = x

    do k = 2, n
      p_next = ((2.0_dp * real(k, dp) - 1.0_dp) * x * p_curr - (real(k, dp) - 1.0_dp) * p_prev) / real(k, dp)
      p_prev = p_curr
      p_curr = p_next
    end do

    p = p_curr
    if (abs(x * x - 1.0_dp) > 1.0e-14_dp) then
      dp_dx = real(n, dp) * (p_prev - x * p) / (1.0_dp - x * x)
    else if (x > 0.0_dp) then
      dp_dx = 0.5_dp * real(n * (n + 1), dp)
    else
      dp_dx = 0.5_dp * real(n * (n + 1), dp) * (-1.0_dp)**(n - 1)
    end if
  end subroutine legendre_poly_and_deriv

  pure function lagrange_interpolant(i, x, nodes) result(val)
    ! Evaluates i-th Lagrange cardinal polynomial ell_i(x) through given nodes
    integer(i4), intent(in) :: i
    real(dp), intent(in) :: x
    real(dp), dimension(:), intent(in) :: nodes
    real(dp) :: val

    integer(i4) :: j, n

    n = size(nodes)
    val = 1.0_dp

    do j = 1, n
      if (j /= i) then
        val = val * (x - nodes(j)) / (nodes(i) - nodes(j))
      end if
    end do
  end function lagrange_interpolant

  subroutine nodal_differentiation_matrix(nodes, d_mat, err)
    ! Computes nodal differentiation matrix D_ij = ell'_j(xi_i)
    real(dp), dimension(:), intent(in) :: nodes
    real(dp), dimension(:,:), intent(out) :: d_mat
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j, k, n
    real(dp) :: c_i, c_j, sum_val

    n = size(nodes)
    if (size(d_mat, 1) /= n .or. size(d_mat, 2) /= n) then
      if (present(err)) call err%set(UMBRA_ERR_QUADRATURE_ORDER, &
        "Differentiation matrix dimensions must match node count", "nodal_differentiation_matrix")
      return
    end if

    do i = 1, n
      sum_val = 0.0_dp
      do j = 1, n
        if (i /= j) then
          ! Evaluate weight product c_i / c_j
          c_i = 1.0_dp
          c_j = 1.0_dp
          do k = 1, n
            if (k /= i) c_i = c_i * (nodes(i) - nodes(k))
            if (k /= j) c_j = c_j * (nodes(j) - nodes(k))
          end do
          d_mat(i, j) = (c_i / c_j) / (nodes(i) - nodes(j))
          sum_val = sum_val + d_mat(i, j)
        end if
      end do
      ! Negative sum condition ensures constant polynomials have zero derivative
      d_mat(i, i) = -sum_val
    end do

    if (present(err)) call err%clear()
  end subroutine nodal_differentiation_matrix

  subroutine modal_to_nodal_vandermonde(nodes, v_mat, err)
    ! Computes generalized Vandermonde matrix V_ij = P_{j-1}(xi_i)
    real(dp), dimension(:), intent(in) :: nodes
    real(dp), dimension(:,:), intent(out) :: v_mat
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j, n
    real(dp) :: p_val, dp_val

    n = size(nodes)
    if (size(v_mat, 1) /= n .or. size(v_mat, 2) /= n) then
      if (present(err)) call err%set(UMBRA_ERR_QUADRATURE_ORDER, &
        "Vandermonde matrix dimensions must match node count", "modal_to_nodal_vandermonde")
      return
    end if

    do i = 1, n
      do j = 1, n
        call legendre_poly_and_deriv(j - 1, nodes(i), p_val, dp_val)
        v_mat(i, j) = p_val
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine modal_to_nodal_vandermonde

end module umbra_polynomials
