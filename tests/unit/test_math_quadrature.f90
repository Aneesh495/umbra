!> Unit tests for numerical quadrature and Legendre polynomials
module test_math_quadrature
  use umbra_kinds, only: dp, i4
  use umbra_quadrature, only: gauss_legendre_nodes_weights, gauss_lobatto_nodes_weights
  use umbra_polynomials, only: legendre_poly_and_deriv
  use test_framework, only: test_suite_t, assert_true, assert_close_dp
  implicit none
  private

  public :: run_test_math_quadrature

contains

  subroutine run_test_math_quadrature(suite)
    type(test_suite_t), intent(inout) :: suite
    real(dp), allocatable :: x_nodes(:), w_weights(:)
    real(dp) :: sum_weights, p_val, dp_val
    integer(i4) :: n

    ! Test 1: Gauss-Legendre weights sum to 2.0 exactly
    do n = 2, 8
      allocate(x_nodes(n), w_weights(n))
      call gauss_legendre_nodes_weights(n, x_nodes, w_weights)
      sum_weights = sum(w_weights)
      call assert_close_dp(suite, sum_weights, 2.0_dp, 1.0e-13_dp, "Gauss-Legendre weight sum = 2.0")
      deallocate(x_nodes, w_weights)
    end do

    ! Test 2: Gauss-Lobatto endpoints at -1 and +1
    allocate(x_nodes(4), w_weights(4))
    call gauss_lobatto_nodes_weights(4_i4, x_nodes, w_weights)
    call assert_close_dp(suite, x_nodes(1), -1.0_dp, 1.0e-14_dp, "Gauss-Lobatto left endpoint -1")
    call assert_close_dp(suite, x_nodes(4), 1.0_dp, 1.0e-14_dp, "Gauss-Lobatto right endpoint +1")
    deallocate(x_nodes, w_weights)

    ! Test 3: Legendre polynomial values at x=1 (P_n(1) = 1 for all n)
    do n = 0, 5
      call legendre_poly_and_deriv(n, 1.0_dp, p_val, dp_val)
      call assert_close_dp(suite, p_val, 1.0_dp, 1.0e-14_dp, "Legendre P_n(1) = 1")
    end do

    ! Test 4: Legendre polynomial derivative at x=0 for P_1(x)=x => P'_1(0) = 1
    call legendre_poly_and_deriv(1_i4, 0.0_dp, p_val, dp_val)
    call assert_close_dp(suite, dp_val, 1.0_dp, 1.0e-14_dp, "Legendre P'_1(0) = 1")

    ! Test 5: Legendre polynomial derivative at x=0 for P_2(x)=(3x^2-1)/2 => P'_2(0) = 0
    call legendre_poly_and_deriv(2_i4, 0.0_dp, p_val, dp_val)
    call assert_close_dp(suite, dp_val, 0.0_dp, 1.0e-14_dp, "Legendre P'_2(0) = 0")
  end subroutine run_test_math_quadrature

end module test_math_quadrature
