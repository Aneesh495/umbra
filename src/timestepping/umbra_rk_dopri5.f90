! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_rk_dopri5
! Description: Embedded Dormand-Prince 5(4) adaptive Runge-Kutta integrator with error control.

module umbra_rk_dopri5
  use umbra_kinds, only: dp, i4
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  public :: dopri5_adaptive_step

  abstract interface
    subroutine ode_deriv_func(t, y, dydt, err)
      import :: dp, umbra_error_t
      real(dp), intent(in) :: t
      real(dp), dimension(:), intent(in) :: y
      real(dp), dimension(:), intent(out) :: dydt
      type(umbra_error_t), intent(inout), optional :: err
    end subroutine ode_deriv_func
  end interface

  ! Dormand-Prince 5(4) Butcher Tableau Coefficients
  real(dp), parameter :: c2 = 1.0_dp / 5.0_dp
  real(dp), parameter :: c3 = 3.0_dp / 10.0_dp
  real(dp), parameter :: c4 = 4.0_dp / 5.0_dp
  real(dp), parameter :: c5 = 8.0_dp / 9.0_dp

  real(dp), parameter :: a21 = 1.0_dp / 5.0_dp
  real(dp), parameter :: a31 = 3.0_dp / 40.0_dp,        a32 = 9.0_dp / 40.0_dp
  real(dp), parameter :: a41 = 44.0_dp / 45.0_dp,       a42 = -56.0_dp / 15.0_dp,     a43 = 32.0_dp / 9.0_dp
  real(dp), parameter :: a51 = 19372.0_dp / 6561.0_dp,  a52 = -25360.0_dp / 2187.0_dp, a53 = 64448.0_dp / 6561.0_dp, a54 = -212.0_dp / 729.0_dp
  real(dp), parameter :: a61 = 9017.0_dp / 3168.0_dp,   a62 = -355.0_dp / 33.0_dp,    a63 = 46732.0_dp / 5247.0_dp, a64 = 49.0_dp / 176.0_dp, a65 = -5103.0_dp / 18656.0_dp
  real(dp), parameter :: a71 = 35.0_dp / 384.0_dp,      a73 = 500.0_dp / 1113.0_dp,   a74 = 125.0_dp / 192.0_dp,    a75 = -2187.0_dp / 6784.0_dp, a76 = 11.0_dp / 84.0_dp

  ! 5th-order solution weights
  real(dp), parameter :: b1 = 35.0_dp / 384.0_dp,       b3 = 500.0_dp / 1113.0_dp,    b4 = 125.0_dp / 192.0_dp,     b5 = -2187.0_dp / 6784.0_dp, b6 = 11.0_dp / 84.0_dp

  ! 4th-order solution weights for error estimation: e_i = b_i - bhat_i
  real(dp), parameter :: e1 = 71.0_dp / 57600.0_dp,     e3 = -71.0_dp / 16695.0_dp,   e4 = 71.0_dp / 1920.0_dp,     e5 = -17253.0_dp / 339200.0_dp, e6 = 22.0_dp / 525.0_dp, e7 = -1.0_dp / 40.0_dp

contains

  subroutine dopri5_adaptive_step(derivs, t, y, dt, dt_next, tol, err)
    procedure(ode_deriv_func) :: derivs
    real(dp), intent(inout) :: t
    real(dp), dimension(:), intent(inout) :: y
    real(dp), intent(inout) :: dt
    real(dp), intent(out) :: dt_next
    real(dp), intent(in) :: tol
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), dimension(size(y)) :: k1, k2, k3, k4, k5, k6, k7, y_stage, err_est
    real(dp) :: err_max, scale, fac
    integer(i4) :: n, i
    integer(i4) :: iter

    n = size(y)

    do iter = 1, 10
      ! Stage 1
      call derivs(t, y, k1, err)

      ! Stage 2
      y_stage = y + dt * (a21 * k1)
      call derivs(t + c2 * dt, y_stage, k2, err)

      ! Stage 3
      y_stage = y + dt * (a31 * k1 + a32 * k2)
      call derivs(t + c3 * dt, y_stage, k3, err)

      ! Stage 4
      y_stage = y + dt * (a41 * k1 + a42 * k2 + a43 * k3)
      call derivs(t + c4 * dt, y_stage, k4, err)

      ! Stage 5
      y_stage = y + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
      call derivs(t + c5 * dt, y_stage, k5, err)

      ! Stage 6
      y_stage = y + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
      call derivs(t + dt, y_stage, k6, err)

      ! 5th-order new state
      y_stage = y + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)

      ! Stage 7 (FSAL property: k7 is k1 for next step)
      call derivs(t + dt, y_stage, k7, err)

      ! Error estimate vector: E = dt * sum e_i * k_i
      err_est = dt * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7)

      err_max = 0.0_dp
      do i = 1, n
        scale = tol * (1.0_dp + max(abs(y(i)), abs(y_stage(i))))
        err_max = max(err_max, abs(err_est(i)) / scale)
      end do

      ! Step adjustment factor
      if (err_max <= 1.0_dp) then
        ! Step accepted
        t = t + dt
        y = y_stage
        fac = min(2.0_dp, max(0.2_dp, 0.9_dp * (1.0_dp / max(err_max, 1.0e-10_dp))**0.2_dp))
        dt_next = dt * fac
        if (present(err)) call err%clear()
        return
      else
        ! Step rejected: reduce dt and retry
        fac = min(1.0_dp, max(0.1_dp, 0.9_dp * (1.0_dp / err_max)**0.25_dp))
        dt = dt * fac
      end if
    end do

    dt_next = dt
    if (present(err)) call err%clear()
  end subroutine dopri5_adaptive_step

end module umbra_rk_dopri5
