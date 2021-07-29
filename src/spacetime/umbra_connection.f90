! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_connection
! Description: Christoffel connection coefficients, metric derivatives, and geodesic accelerations.

module umbra_connection
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NDIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t
  implicit none
  private

  public :: compute_christoffel_symbols
  public :: compute_geodesic_acceleration
  public :: compute_metric_derivatives

contains

  subroutine compute_metric_derivatives(metric, x, dg, err)
    ! Computes dg(mu, nu, lambda) = partial_lambda g_{mu nu} via 4th-order central difference
    class(spacetime_metric_t), intent(in) :: metric
    real(dp), dimension(0:3), intent(in) :: x
    real(dp), dimension(0:3, 0:3, 0:3), intent(out) :: dg
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), parameter :: h = 1.0e-5_dp
    real(dp), parameter :: inv_12h = 1.0_dp / (12.0_dp * h)
    real(dp), dimension(0:3) :: xp2, xp1, xm1, xm2
    type(metric_3p1_t) :: mp2, mp1, mm1, mm2
    integer(i4) :: lam, mu, nu

    dg = 0.0_dp

    ! Perturb along spatial coordinates lambda = 1, 2, 3 (metrics are stationary: lambda=0 derivative is zero)
    do lam = 1, 3
      xp2 = x
      xp1 = x
      xm1 = x
      xm2 = x

      xp2(lam) = x(lam) + 2.0_dp * h
      xp1(lam) = x(lam) + 1.0_dp * h
      xm1(lam) = x(lam) - 1.0_dp * h
      xm2(lam) = x(lam) - 2.0_dp * h

      call metric%eval(xp2, mp2, err)
      call metric%eval(xp1, mp1, err)
      call metric%eval(xm1, mm1, err)
      call metric%eval(xm2, mm2, err)

      do nu = 0, 3
        do mu = 0, 3
          dg(mu, nu, lam) = (-mp2%g_down(mu, nu) + 8.0_dp * mp1%g_down(mu, nu) &
                            - 8.0_dp * mm1%g_down(mu, nu) + mm2%g_down(mu, nu)) * inv_12h
        end do
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine compute_metric_derivatives

  subroutine compute_christoffel_symbols(metric, x, gamma_conn, err)
    ! Computes Gamma^mu_{alpha, beta} = 0.5 * g^{mu, lambda} (d_alpha g_{beta, lambda} + d_beta g_{alpha, lambda} - d_lambda g_{alpha, beta})
    class(spacetime_metric_t), intent(in) :: metric
    real(dp), dimension(0:3), intent(in) :: x
    real(dp), dimension(0:3, 0:3, 0:3), intent(out) :: gamma_conn
    type(umbra_error_t), intent(inout), optional :: err

    type(metric_3p1_t) :: m
    real(dp), dimension(0:3, 0:3, 0:3) :: dg
    integer(i4) :: mu, alpha, beta, lambda
    real(dp) :: term

    call metric%eval(x, m, err)
    call compute_metric_derivatives(metric, x, dg, err)

    gamma_conn = 0.0_dp

    do mu = 0, 3
      do alpha = 0, 3
        do beta = alpha, 3
          term = 0.0_dp
          do lambda = 0, 3
            term = term + 0.5_dp * m%g_up(mu, lambda) * &
              (dg(beta, lambda, alpha) + dg(alpha, lambda, beta) - dg(alpha, beta, lambda))
          end do
          gamma_conn(mu, alpha, beta) = term
          gamma_conn(mu, beta, alpha) = term ! Torsion-free Levi-Civita connection symmetry
        end do
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine compute_christoffel_symbols

  pure subroutine compute_geodesic_acceleration(gamma_conn, u_up, a_up)
    ! Computes a^mu = -Gamma^mu_{alpha, beta} u^alpha u^beta
    real(dp), dimension(0:3, 0:3, 0:3), intent(in) :: gamma_conn
    real(dp), dimension(0:3), intent(in) :: u_up
    real(dp), dimension(0:3), intent(out) :: a_up

    integer(i4) :: mu, alpha, beta

    a_up = 0.0_dp
    do mu = 0, 3
      do alpha = 0, 3
        do beta = 0, 3
          a_up(mu) = a_up(mu) - gamma_conn(mu, alpha, beta) * u_up(alpha) * u_up(beta)
        end do
      end do
    end do
  end subroutine compute_geodesic_acceleration

end module umbra_connection
