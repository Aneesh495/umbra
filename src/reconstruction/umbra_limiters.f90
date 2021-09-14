! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_limiters
! Description: Total Variation Diminishing (TVD) slope limiters for shock capturing.

module umbra_limiters
  use umbra_kinds, only: dp, i4
  implicit none
  private

  integer(i4), parameter, public :: LIMITER_MINMOD   = 1_i4
  integer(i4), parameter, public :: LIMITER_MC       = 2_i4
  integer(i4), parameter, public :: LIMITER_SUPERBEE = 3_i4
  integer(i4), parameter, public :: LIMITER_VANLEER  = 4_i4

  public :: slope_limiter

contains

  pure function slope_limiter(limiter_type, dl, dr) result(slope)
    integer(i4), intent(in) :: limiter_type
    real(dp), intent(in) :: dl, dr
    real(dp) :: slope, r, phi

    ! If slopes have opposite signs, extrema present: zero slope
    if (dl * dr <= 0.0_dp) then
      slope = 0.0_dp
      return
    end if

    r = dl / dr

    select case (limiter_type)
      case (LIMITER_MINMOD)
        phi = max(0.0_dp, min(1.0_dp, r))

      case (LIMITER_MC)
        ! Monotonized Central (MC) limiter
        phi = max(0.0_dp, min(2.0_dp, min(2.0_dp * r, 0.5_dp * (1.0_dp + r))))

      case (LIMITER_SUPERBEE)
        phi = max(0.0_dp, max(min(1.0_dp, 2.0_dp * r), min(2.0_dp, r)))

      case (LIMITER_VANLEER)
        phi = (r + abs(r)) / (1.0_dp + abs(r))

      case default
        phi = max(0.0_dp, min(1.0_dp, r))
    end select

    slope = phi * dr
  end function slope_limiter

end module umbra_limiters
