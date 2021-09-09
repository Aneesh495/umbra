! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_glm
! Description: Hyperbolic divergence cleaning (GLM-MHD) Dedner formulation with mixed damping.

module umbra_glm
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: GLM_ALPHA_P, GLM_CR
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  type, public :: glm_state_t
    real(dp) :: ch = 1.0_dp          ! Hyperbolic divergence propagation speed
    real(dp) :: cp = 0.18_dp         ! Parabolic damping parameter
    real(dp) :: damping_rate = 1.0_dp
  contains
    procedure, public :: init => glm_init
    procedure, public :: damping_factor => glm_damping_factor
  end type glm_state_t

  public :: glm_source_step_2d
  public :: glm_flux_gradient_2d

contains

  subroutine glm_init(this, max_signal_speed, dx_min)
    class(glm_state_t), intent(inout) :: this
    real(dp), intent(in) :: max_signal_speed, dx_min

    this%ch = max(max_signal_speed, 1.0e-3_dp)
    this%cp = sqrt(dx_min * this%ch * GLM_ALPHA_P)
    this%damping_rate = (this%ch * this%ch) / max(this%cp * this%cp, 1.0e-30_dp)
  end subroutine glm_init

  pure function glm_damping_factor(this, dt) result(factor)
    class(glm_state_t), intent(in) :: this
    real(dp), intent(in) :: dt
    real(dp) :: factor

    ! Exact exponential analytical decay: exp(-dt * ch^2 / cp^2)
    factor = exp(-dt * this%damping_rate)
  end function glm_damping_factor

  subroutine glm_source_step_2d(glm, psi, dt, nx, ny, err)
    type(glm_state_t), intent(in) :: glm
    real(dp), dimension(:,:), intent(inout) :: psi
    real(dp), intent(in) :: dt
    integer(i4), intent(in) :: nx, ny
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: factor
    integer(i4) :: i, j

    factor = glm%damping_factor(dt)
    do j = 1, ny
      do i = 1, nx
        psi(i, j) = psi(i, j) * factor
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine glm_source_step_2d

  subroutine glm_flux_gradient_2d(psi, div_b, ch, dx, dy, nx, ny, d_bx, d_by, d_psi, err)
    ! Evaluates Dedner GLM coupling terms:
    ! d_Bx/dt += -d_psi/dx
    ! d_By/dt += -d_psi/dy
    ! d_psi/dt += -ch^2 * div(B)
    real(dp), dimension(:,:), intent(in) :: psi
    real(dp), dimension(:,:), intent(in) :: div_b
    real(dp), intent(in) :: ch, dx, dy
    integer(i4), intent(in) :: nx, ny
    real(dp), dimension(:,:), intent(out) :: d_bx
    real(dp), dimension(:,:), intent(out) :: d_by
    real(dp), dimension(:,:), intent(out) :: d_psi
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: inv_2dx, inv_2dy, ch_sq
    integer(i4) :: i, j

    inv_2dx = 0.5_dp / dx
    inv_2dy = 0.5_dp / dy
    ch_sq = ch * ch

    do j = 1, ny
      do i = 1, nx
        d_bx(i, j) = -(psi(min(i + 1, nx), j) - psi(max(i - 1, 1), j)) * inv_2dx
        d_by(i, j) = -(psi(i, min(j + 1, ny)) - psi(i, max(j - 1, 1))) * inv_2dy
        d_psi(i, j) = -ch_sq * div_b(i, j)
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine glm_flux_gradient_2d

end module umbra_glm
