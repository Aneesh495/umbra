! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_reconstruct_weno
! Description: Fifth-order Weighted Essentially Non-Oscillatory (WENO-5) spatial reconstruction.

module umbra_reconstruct_weno
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  public :: reconstruct_weno5_1d

contains

  subroutine reconstruct_weno5_1d(q, n, nghost, q_l, q_r, err)
    ! Reconstructs left and right interface states at interfaces i = 1 .. n+1
    integer(i4), intent(in) :: n, nghost
    real(dp), dimension(NPRIM, 1-nghost : n+nghost), intent(in) :: q
    real(dp), dimension(NPRIM, 1 : n+1), intent(out) :: q_l
    real(dp), dimension(NPRIM, 1 : n+1), intent(out) :: q_r
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), parameter :: eps_weno = 1.0e-6_dp
    real(dp), parameter :: c13_12 = 13.0_dp / 12.0_dp
    real(dp) :: v0, v1, v2, v3, v4
    real(dp) :: beta0, beta1, beta2
    real(dp) :: w0, w1, w2, alpha0, alpha1, alpha2, sum_alpha
    real(dp) :: p0, p1, p2
    integer(i4) :: i, v

    ! Left state q_l at interface i-1/2 from stencil [i-3, i-2, i-1, i, i+1]
    do i = 1, n + 1
      do v = 1, NPRIM
        ! Reconstruction from left cell i-1
        v0 = q(v, i - 3)
        v1 = q(v, i - 2)
        v2 = q(v, i - 1)
        v3 = q(v, i)
        v4 = q(v, i + 1)

        ! Smoothness indicators beta_0, beta_1, beta_2
        beta0 = c13_12 * (v0 - 2.0_dp * v1 + v2)**2 + 0.25_dp * (v0 - 4.0_dp * v1 + 3.0_dp * v2)**2
        beta1 = c13_12 * (v1 - 2.0_dp * v2 + v3)**2 + 0.25_dp * (v1 - v3)**2
        beta2 = c13_12 * (v2 - 2.0_dp * v3 + v4)**2 + 0.25_dp * (3.0_dp * v2 - 4.0_dp * v3 + v4)**2

        ! Unnormalized nonlinear weights
        alpha0 = 0.1_dp / (eps_weno + beta0)**2
        alpha1 = 0.6_dp / (eps_weno + beta1)**2
        alpha2 = 0.3_dp / (eps_weno + beta2)**2
        sum_alpha = alpha0 + alpha1 + alpha2

        w0 = alpha0 / sum_alpha
        w1 = alpha1 / sum_alpha
        w2 = alpha2 / sum_alpha

        ! Stencil candidate polynomials
        p0 = (1.0_dp / 3.0_dp) * v0 - (7.0_dp / 6.0_dp) * v1 + (11.0_dp / 6.0_dp) * v2
        p1 = -(1.0_dp / 6.0_dp) * v1 + (5.0_dp / 6.0_dp) * v2 + (1.0_dp / 3.0_dp) * v3
        p2 = (1.0_dp / 3.0_dp) * v2 + (5.0_dp / 6.0_dp) * v3 - (1.0_dp / 6.0_dp) * v4

        q_l(v, i) = w0 * p0 + w1 * p1 + w2 * p2

        ! Reconstruction from right cell i: mirrored stencil [i+2, i+1, i, i-1, i-2]
        v0 = q(v, i + 2)
        v1 = q(v, i + 1)
        v2 = q(v, i)
        v3 = q(v, i - 1)
        v4 = q(v, i - 2)

        beta0 = c13_12 * (v0 - 2.0_dp * v1 + v2)**2 + 0.25_dp * (v0 - 4.0_dp * v1 + 3.0_dp * v2)**2
        beta1 = c13_12 * (v1 - 2.0_dp * v2 + v3)**2 + 0.25_dp * (v1 - v3)**2
        beta2 = c13_12 * (v2 - 2.0_dp * v3 + v4)**2 + 0.25_dp * (3.0_dp * v2 - 4.0_dp * v3 + v4)**2

        alpha0 = 0.1_dp / (eps_weno + beta0)**2
        alpha1 = 0.6_dp / (eps_weno + beta1)**2
        alpha2 = 0.3_dp / (eps_weno + beta2)**2
        sum_alpha = alpha0 + alpha1 + alpha2

        w0 = alpha0 / sum_alpha
        w1 = alpha1 / sum_alpha
        w2 = alpha2 / sum_alpha

        p0 = (1.0_dp / 3.0_dp) * v0 - (7.0_dp / 6.0_dp) * v1 + (11.0_dp / 6.0_dp) * v2
        p1 = -(1.0_dp / 6.0_dp) * v1 + (5.0_dp / 6.0_dp) * v2 + (1.0_dp / 3.0_dp) * v3
        p2 = (1.0_dp / 3.0_dp) * v2 + (5.0_dp / 6.0_dp) * v3 - (1.0_dp / 6.0_dp) * v4

        q_r(v, i) = w0 * p0 + w1 * p1 + w2 * p2
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine reconstruct_weno5_1d

end module umbra_reconstruct_weno
