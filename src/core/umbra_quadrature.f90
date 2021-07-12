! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_quadrature
! Description: High-order Gauss-Legendre and Gauss-Lobatto-Legendre quadrature rules.

module umbra_quadrature
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: PI, EPSILON_DP
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_QUADRATURE_ORDER
  implicit none
  private

  integer(i4), parameter, public :: MAX_QUAD_ORDER = 8_i4

  public :: gauss_legendre_nodes_weights
  public :: gauss_lobatto_nodes_weights

contains

  subroutine gauss_legendre_nodes_weights(n, nodes, weights, err)
    integer(i4), intent(in) :: n
    real(dp), dimension(:), intent(out) :: nodes
    real(dp), dimension(:), intent(out) :: weights
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j, k, m
    real(dp) :: z, z1, p1, p2, pp
    real(dp), parameter :: tol = 1.0e-15_dp
    integer(i4), parameter :: max_iter = 100_i4

    if (n < 1 .or. n > size(nodes) .or. n > size(weights)) then
      if (present(err)) call err%set(UMBRA_ERR_QUADRATURE_ORDER, &
        "Invalid quadrature order or array dimension", "gauss_legendre_nodes_weights")
      return
    end if

    ! Pre-computed exact values for small N for maximum precision and speed
    select case (n)
      case (1)
        nodes(1) = 0.0_dp
        weights(1) = 2.0_dp
        return

      case (2)
        nodes(1) = -1.0_dp / sqrt(3.0_dp)
        nodes(2) =  1.0_dp / sqrt(3.0_dp)
        weights(1) = 1.0_dp
        weights(2) = 1.0_dp
        return

      case (3)
        nodes(1) = -sqrt(0.6_dp)
        nodes(2) =  0.0_dp
        nodes(3) =  sqrt(0.6_dp)
        weights(1) = 5.0_dp / 9.0_dp
        weights(2) = 8.0_dp / 9.0_dp
        weights(3) = 5.0_dp / 9.0_dp
        return

      case (4)
        nodes(1) = -sqrt((3.0_dp + 2.0_dp * sqrt(6.0_dp / 5.0_dp)) / 7.0_dp)
        nodes(2) = -sqrt((3.0_dp - 2.0_dp * sqrt(6.0_dp / 5.0_dp)) / 7.0_dp)
        nodes(3) = -nodes(2)
        nodes(4) = -nodes(1)
        weights(1) = (18.0_dp - sqrt(30.0_dp)) / 36.0_dp
        weights(2) = (18.0_dp + sqrt(30.0_dp)) / 36.0_dp
        weights(3) = weights(2)
        weights(4) = weights(1)
        return

      case default
        ! High-order Newton-Raphson calculation via Legendre polynomial recurrence
        m = (n + 1) / 2
        do i = 1, m
          ! Initial Chebyshev-Gauss approximation
          z = cos(PI * (real(i, dp) - 0.25_dp) / (real(n, dp) + 0.5_dp))
          do j = 1, max_iter
            p1 = 1.0_dp
            p2 = 0.0_dp
            do k = 1, n
              ! Bonnet's recurrence relation
              pp = p2
              p2 = p1
              p1 = ((2.0_dp * real(k, dp) - 1.0_dp) * z * p2 - (real(k, dp) - 1.0_dp) * pp) / real(k, dp)
            end do
            ! Derivative of P_n(z)
            pp = real(n, dp) * (z * p1 - p2) / (z * z - 1.0_dp)
            z1 = z
            z = z1 - p1 / pp
            if (abs(z - z1) <= tol) exit
          end do

          nodes(i) = -z
          nodes(n + 1 - i) = z
          weights(i) = 2.0_dp / ((1.0_dp - z * z) * pp * pp)
          weights(n + 1 - i) = weights(i)
        end do
    end select

    if (present(err)) call err%clear()
  end subroutine gauss_legendre_nodes_weights

  subroutine gauss_lobatto_nodes_weights(n, nodes, weights, err)
    ! Computes n Gauss-Lobatto-Legendre (LGL) nodes including endpoints -1 and +1
    integer(i4), intent(in) :: n
    real(dp), dimension(:), intent(out) :: nodes
    real(dp), dimension(:), intent(out) :: weights
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j, k
    real(dp) :: z, z1, p, p_prev, p_next, p_deriv2
    real(dp), parameter :: tol = 1.0e-15_dp
    integer(i4), parameter :: max_iter = 100_i4

    if (n < 2 .or. n > size(nodes) .or. n > size(weights)) then
      if (present(err)) call err%set(UMBRA_ERR_QUADRATURE_ORDER, &
        "Gauss-Lobatto requires order N >= 2", "gauss_lobatto_nodes_weights")
      return
    end if

    select case (n)
      case (2)
        nodes(1) = -1.0_dp
        nodes(2) =  1.0_dp
        weights(1) = 1.0_dp
        weights(2) = 1.0_dp
        return

      case (3)
        nodes(1) = -1.0_dp
        nodes(2) =  0.0_dp
        nodes(3) =  1.0_dp
        weights(1) = 1.0_dp / 3.0_dp
        weights(2) = 4.0_dp / 3.0_dp
        weights(3) = 1.0_dp / 3.0_dp
        return

      case (4)
        nodes(1) = -1.0_dp
        nodes(2) = -1.0_dp / sqrt(5.0_dp)
        nodes(3) =  1.0_dp / sqrt(5.0_dp)
        nodes(4) =  1.0_dp
        weights(1) = 1.0_dp / 6.0_dp
        weights(2) = 5.0_dp / 6.0_dp
        weights(3) = 5.0_dp / 6.0_dp
        weights(4) = 1.0_dp / 6.0_dp
        return

      case default
        nodes(1) = -1.0_dp
        nodes(n) =  1.0_dp
        weights(1) = 2.0_dp / (real(n * (n - 1), dp))
        weights(n) = weights(1)

        ! Internal nodes are the zeros of P'_{n-1}(x)
        do i = 2, n - 1
          ! Initial guess based on Chebyshev-Lobatto points
          z = -cos(PI * real(i - 1, dp) / real(n - 1, dp))
          do j = 1, max_iter
            ! Evaluate P_{n-1}(z), P'_{n-1}(z), and P''_{n-1}(z)
            p_prev = 1.0_dp
            p = z
            do k = 2, n - 1
              p_next = ((2.0_dp * real(k, dp) - 1.0_dp) * z * p - (real(k, dp) - 1.0_dp) * p_prev) / real(k, dp)
              p_prev = p
              p = p_next
            end do
            ! p is now P_{n-1}(z), p_prev is P_{n-2}(z)
            ! First derivative: (1-z^2) P'_{n-1} = (n-1)(P_{n-2} - z P_{n-1})
            ! Second derivative: (1-z^2) P''_{n-1} = 2 z P'_{n-1} - (n-1) n P_{n-1}
            p_deriv2 = (2.0_dp * z * real(n - 1, dp) * (p_prev - z * p) / (1.0_dp - z * z) &
                       - real(n * (n - 1), dp) * p) / (1.0_dp - z * z)

            z1 = z
            z = z1 - (real(n - 1, dp) * (p_prev - z * p) / (1.0_dp - z * z)) / p_deriv2
            if (abs(z - z1) <= tol) exit
          end do
          nodes(i) = z
          ! Weight formula: w_i = 2 / (n(n-1) [P_{n-1}(x_i)]^2)
          weights(i) = 2.0_dp / (real(n * (n - 1), dp) * p * p)
        end do
    end select

    if (present(err)) call err%clear()
  end subroutine gauss_lobatto_nodes_weights

end module umbra_quadrature
