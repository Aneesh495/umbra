! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_reconstruct_tvd
! Description: Second-order TVD MUSCL interface state reconstruction.

module umbra_reconstruct_tvd
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_limiters, only: slope_limiter, LIMITER_MC
  implicit none
  private

  public :: reconstruct_tvd_1d

contains

  subroutine reconstruct_tvd_1d(q, n, nghost, limiter_type, q_l, q_r, err)
    ! Reconstructs left and right interface states at interfaces i = 1 .. n+1
    ! q_l(v, i) is q^-(i - 1/2) coming from left cell i-1
    ! q_r(v, i) is q^+(i - 1/2) coming from right cell i
    integer(i4), intent(in) :: n, nghost, limiter_type
    real(dp), dimension(NPRIM, 1-nghost : n+nghost), intent(in) :: q
    real(dp), dimension(NPRIM, 1 : n+1), intent(out) :: q_l
    real(dp), dimension(NPRIM, 1 : n+1), intent(out) :: q_r
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), dimension(NPRIM, 1-nghost : n+nghost) :: dq
    real(dp) :: dl, dr
    integer(i4) :: i, v

    ! Compute limited slopes dq(:, i) inside cell i
    do i = 1 - nghost + 1, n + nghost - 1
      do v = 1, NPRIM
        dl = q(v, i) - q(v, i - 1)
        dr = q(v, i + 1) - q(v, i)
        dq(v, i) = slope_limiter(limiter_type, dl, dr)
      end do
    end do

    ! Interface reconstruction at face i-1/2:
    ! Left state from cell i-1: q_l(:, i) = q(:, i-1) + 0.5 * dq(:, i-1)
    ! Right state from cell i:  q_r(:, i) = q(:, i)   - 0.5 * dq(:, i)
    do i = 1, n + 1
      do v = 1, NPRIM
        q_l(v, i) = q(v, i - 1) + 0.5_dp * dq(v, i - 1)
        q_r(v, i) = q(v, i)     - 0.5_dp * dq(v, i)
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine reconstruct_tvd_1d

end module umbra_reconstruct_tvd
