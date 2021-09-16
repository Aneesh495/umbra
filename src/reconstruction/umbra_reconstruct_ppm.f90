! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_reconstruct_ppm
! Description: High-order Piecewise Parabolic Method (PPM) spatial reconstruction with monotonicity constraints.

module umbra_reconstruct_ppm
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  public :: reconstruct_ppm_1d

contains

  subroutine reconstruct_ppm_1d(q, n, nghost, q_l, q_r, err)
    ! PPM reconstruction returning left and right face states at interface i-1/2 for i = 1 .. n+1
    integer(i4), intent(in) :: n, nghost
    real(dp), dimension(NPRIM, 1-nghost : n+nghost), intent(in) :: q
    real(dp), dimension(NPRIM, 1 : n+1), intent(out) :: q_l
    real(dp), dimension(NPRIM, 1 : n+1), intent(out) :: q_r
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), dimension(NPRIM, 1-nghost : n+nghost) :: q_minus, q_plus
    real(dp) :: q_face, dql, dqr, diff, term
    integer(i4) :: i, v

    ! 4th-order polynomial interpolation to cell edges:
    ! q_{i+1/2} = 7/12 (q_i + q_{i+1}) - 1/12 (q_{i+2} + q_{i-1})
    do i = 1 - nghost + 2, n + nghost - 2
      do v = 1, NPRIM
        q_face = (7.0_dp / 12.0_dp) * (q(v, i) + q(v, i + 1)) - (1.0_dp / 12.0_dp) * (q(v, i + 2) + q(v, i - 1))
        q_minus(v, i + 1) = q_face ! Left boundary of cell i+1
        q_plus(v, i)      = q_face ! Right boundary of cell i
      end do
    end do

    ! Apply Colella-Woodward monotonicity filter within each cell i
    do i = 1 - nghost + 2, n + nghost - 2
      do v = 1, NPRIM
        dql = q_minus(v, i) - q(v, i)
        dqr = q_plus(v, i)  - q(v, i)

        ! If local extremum, flatten to cell average
        if (dql * dqr >= 0.0_dp) then
          q_minus(v, i) = q(v, i)
          q_plus(v, i)  = q(v, i)
          cycle
        end if

        diff = q_plus(v, i) - q_minus(v, i)
        term = 6.0_dp * (q(v, i) - 0.5_dp * (q_minus(v, i) + q_plus(v, i)))

        if (diff * term > diff * diff) then
          q_minus(v, i) = 3.0_dp * q(v, i) - 2.0_dp * q_plus(v, i)
        else if (diff * term < -diff * diff) then
          q_plus(v, i) = 3.0_dp * q(v, i) - 2.0_dp * q_minus(v, i)
        end if
      end do
    end do

    ! Assign to interface arrays:
    ! Interface i-1/2: left state is q_plus from cell i-1, right state is q_minus from cell i
    do i = 1, n + 1
      do v = 1, NPRIM
        q_l(v, i) = q_plus(v, i - 1)
        q_r(v, i) = q_minus(v, i)
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine reconstruct_ppm_1d

end module umbra_reconstruct_ppm
