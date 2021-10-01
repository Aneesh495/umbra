! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_rk_ssprk
! Description: Strong Stability Preserving Runge-Kutta (SSPRK2, SSPRK3, SSPRK4) time integrators.

module umbra_rk_ssprk
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  public :: ssprk2_step_1d
  public :: ssprk3_step_1d

  abstract interface
    subroutine rhs_eval_func_1d(u_in, u_rhs, err)
      import :: dp, i4, NCONS, umbra_error_t
      real(dp), dimension(:,:), intent(in) :: u_in
      real(dp), dimension(:,:), intent(out) :: u_rhs
      type(umbra_error_t), intent(inout), optional :: err
    end subroutine rhs_eval_func_1d
  end interface

contains

  subroutine ssprk2_step_1d(rhs_func, u, ncells, dt, err)
    ! Two-stage second-order Strong Stability Preserving Runge-Kutta (Heun method)
    procedure(rhs_eval_func_1d) :: rhs_func
    integer(i4), intent(in) :: ncells
    real(dp), dimension(NCONS, ncells), intent(inout) :: u
    real(dp), intent(in) :: dt
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), dimension(NCONS, ncells) :: u0, u1, l0, l1

    u0 = u

    ! Stage 1: U^(1) = U^n + dt * L(U^n)
    call rhs_func(u0, l0, err)
    u1 = u0 + dt * l0

    ! Stage 2: U^(n+1) = 0.5 * U^n + 0.5 * U^(1) + 0.5 * dt * L(U^(1))
    call rhs_func(u1, l1, err)
    u = 0.5_dp * u0 + 0.5_dp * u1 + 0.5_dp * dt * l1

    if (present(err)) call err%clear()
  end subroutine ssprk2_step_1d

  subroutine ssprk3_step_1d(rhs_func, u, ncells, dt, err)
    ! Three-stage third-order Shu-Osher Strong Stability Preserving Runge-Kutta
    procedure(rhs_eval_func_1d) :: rhs_func
    integer(i4), intent(in) :: ncells
    real(dp), dimension(NCONS, ncells), intent(inout) :: u
    real(dp), intent(in) :: dt
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), dimension(NCONS, ncells) :: u0, u1, u2, l0, l1, l2

    u0 = u

    ! Stage 1: U^(1) = U^0 + dt * L(U^0)
    call rhs_func(u0, l0, err)
    u1 = u0 + dt * l0

    ! Stage 2: U^(2) = 0.75 * U^0 + 0.25 * U^(1) + 0.25 * dt * L(U^(1))
    call rhs_func(u1, l1, err)
    u2 = 0.75_dp * u0 + 0.25_dp * u1 + 0.25_dp * dt * l1

    ! Stage 3: U^(n+1) = (1/3) * U^0 + (2/3) * U^(2) + (2/3) * dt * L(U^(2))
    call rhs_func(u2, l2, err)
    u = (1.0_dp / 3.0_dp) * u0 + (2.0_dp / 3.0_dp) * u2 + (2.0_dp / 3.0_dp) * dt * l2

    if (present(err)) call err%clear()
  end subroutine ssprk3_step_1d

end module umbra_rk_ssprk
