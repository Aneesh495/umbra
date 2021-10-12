! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_polarimetry
! Description: Polarized Stokes vector (I, Q, U, V) transport with Faraday rotation and conversion.

module umbra_polarimetry
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NSTOKES, STOKES_I, STOKES_Q, STOKES_U, STOKES_V
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  type, public :: stokes_vector_t
    real(dp) :: i = 0.0_dp
    real(dp) :: q = 0.0_dp
    real(dp) :: u = 0.0_dp
    real(dp) :: v = 0.0_dp
  contains
    procedure, public :: to_array => stokes_to_array
    procedure, public :: from_array => stokes_from_array
    procedure, public :: linear_polarization_fraction => stokes_m_pol
    procedure, public :: polarization_angle => stokes_chi
  end type stokes_vector_t

  public :: integrate_polarized_transfer_step

contains

  pure function stokes_to_array(this) result(arr)
    class(stokes_vector_t), intent(in) :: this
    real(dp), dimension(NSTOKES) :: arr
    arr(STOKES_I) = this%i
    arr(STOKES_Q) = this%q
    arr(STOKES_U) = this%u
    arr(STOKES_V) = this%v
  end function stokes_to_array

  pure subroutine stokes_from_array(this, arr)
    class(stokes_vector_t), intent(inout) :: this
    real(dp), dimension(NSTOKES), intent(in) :: arr
    this%i = arr(STOKES_I)
    this%q = arr(STOKES_Q)
    this%u = arr(STOKES_U)
    this%v = arr(STOKES_V)
  end subroutine stokes_from_array

  pure function stokes_m_pol(this) result(m_frac)
    class(stokes_vector_t), intent(in) :: this
    real(dp) :: m_frac
    m_frac = sqrt(this%q**2 + this%u**2) / max(this%i, 1.0e-30_dp)
    m_frac = min(max(m_frac, 0.0_dp), 1.0_dp)
  end function stokes_m_pol

  pure function stokes_chi(this) result(angle)
    class(stokes_vector_t), intent(in) :: this
    real(dp) :: angle
    ! Electric vector polarization angle (EVPA): chi = 0.5 * atan2(U, Q)
    angle = 0.5_dp * atan2(this%u, this%q)
  end function stokes_chi

  subroutine integrate_polarized_transfer_step(s, j_i, j_q, j_u, j_v, alpha_i, rho_v, d_lambda, err)
    ! Formal solution step for polarized transfer equation: dS/dlambda = j - K S
    type(stokes_vector_t), intent(inout) :: s
    real(dp), intent(in) :: j_i, j_q, j_u, j_v, alpha_i, rho_v, d_lambda
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: d_tau, trans, cos_rot, sin_rot, q_rot, u_rot

    d_tau = max(alpha_i * d_lambda, 0.0_dp)
    trans = exp(-d_tau)

    ! Source integration
    if (d_tau > 1.0e-6_dp) then
      s%i = s%i * trans + (j_i / max(alpha_i, 1.0e-30_dp)) * (1.0_dp - trans)
    else
      s%i = s%i + j_i * d_lambda
    end if

    ! Faraday rotation in Q-U plane: angle = rho_v * d_lambda
    cos_rot = cos(rho_v * d_lambda)
    sin_rot = sin(rho_v * d_lambda)

    q_rot = s%q * cos_rot - s%u * sin_rot + j_q * d_lambda
    u_rot = s%q * sin_rot + s%u * cos_rot + j_u * d_lambda

    s%q = q_rot * trans
    s%u = u_rot * trans
    s%v = s%v * trans + j_v * d_lambda

    if (present(err)) call err%clear()
  end subroutine integrate_polarized_transfer_step

end module umbra_polarimetry
