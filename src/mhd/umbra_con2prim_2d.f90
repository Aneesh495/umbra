! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_con2prim_2d
! Description: Two-dimensional multidimensional Newton-Raphson conservative-to-primitive solver for extreme magnetization.

module umbra_con2prim_2d
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: DHO_FLOOR, PRESS_FLOOR, GAMMA_MAX, VEL_MAX, TOL_CONVERGENCE
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_CON2PRIM_FAILED
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t
  implicit none
  private

  public :: con2prim_2d_solve

contains

  subroutine con2prim_2d_solve(u, m, eos, p, err)
    type(cons_state_t), intent(in) :: u
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    type(prim_state_t), intent(out) :: p
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: b2, s2, s_dot_b, s_dot_b_sq
    real(dp) :: w, v2, gamma, rho, press, eps, h_val
    real(dp) :: f1, f2, df1_dw, df1_dv2, df2_dw, df2_dv2, det_j, dw, dv2
    real(dp) :: w_plus_b2, res_norm, step_size
    real(dp), dimension(3) :: s_up
    integer(i4) :: iter, i, j
    integer(i4), parameter :: max_iter = 50_i4
    real(dp), parameter :: tol = 1.0e-11_dp

    if (u%d <= 0.0_dp) then
      if (present(err)) call err%set(UMBRA_ERR_CON2PRIM_FAILED, "Non-positive D in con2prim_2d", "con2prim_2d_solve")
      return
    end if

    ! Compute metric contractions
    b2 = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        b2 = b2 + m%gamma_down(i, j) * u%b(i) * u%b(j)
      end do
    end do

    s2 = 0.0_dp
    s_up = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        s_up(i) = s_up(i) + m%gamma_up(i, j) * u%s(j)
        s2 = s2 + m%gamma_up(i, j) * u%s(i) * u%s(j)
      end do
    end do

    s_dot_b = u%s(1)*u%b(1) + u%s(2)*u%b(2) + u%s(3)*u%b(3)
    s_dot_b_sq = s_dot_b * s_dot_b

    ! Initial estimates for W and v^2
    w = max(u%d + u%tau, sqrt(s2) + 1.0e-6_dp)
    v2 = min(s2 / (w * w + 1.0e-12_dp), VEL_MAX * VEL_MAX)

    do iter = 1, max_iter
      w_plus_b2 = w + b2
      gamma = 1.0_dp / sqrt(max(1.0_dp - v2, 1.0e-12_dp))
      gamma = min(gamma, GAMMA_MAX)
      rho = max(u%d / gamma, DHO_FLOOR)

      ! Equation of state evaluations
      press = ((eos%gamma - 1.0_dp) / eos%gamma) * (w - gamma * u%d) / (gamma * gamma)
      press = max(press, PRESS_FLOOR)
      eps = eos%specific_internal_energy(rho, press)
      h_val = eos%specific_enthalpy(rho, eps, press)

      ! F1(W, v^2) = W^2 (W + B^2)^2 v^2 - [S^2 W^2 + 2 (S.B)^2 W + (S.B)^2 B^2] = 0
      f1 = (w**2) * (w_plus_b2**2) * v2 - (s2 * (w**2) + 2.0_dp * s_dot_b_sq * w + s_dot_b_sq * b2)

      ! F2(W, v^2) = W - p + 0.5 * B^2 * (1 + v^2) - 0.5 * (S.B)^2 / W^2 - D - tau = 0
      f2 = w - press + 0.5_dp * b2 * (1.0_dp + v2) - 0.5_dp * s_dot_b_sq / max(w * w, 1.0e-30_dp) - u%d - u%tau

      res_norm = sqrt(f1 * f1 / max(w**4, 1.0e-20_dp) + f2 * f2 / max(w * w, 1.0e-20_dp))
      if (res_norm < tol) exit

      ! 2x2 Jacobian evaluation
      df1_dw = 2.0_dp * w * (w_plus_b2**2) * v2 + 2.0_dp * (w**2) * w_plus_b2 * v2 - (2.0_dp * s2 * w + 2.0_dp * s_dot_b_sq)
      df1_dv2 = (w**2) * (w_plus_b2**2)

      df2_dw = 1.0_dp - ((eos%gamma - 1.0_dp) / (eos%gamma * gamma * gamma)) + s_dot_b_sq / max(w**3, 1.0e-30_dp)
      df2_dv2 = 0.5_dp * b2 - 0.5_dp * (press / (1.0_dp - v2))

      det_j = df1_dw * df2_dv2 - df1_dv2 * df2_dw
      if (abs(det_j) < 1.0e-30_dp) then
        ! Ill-conditioned Jacobian: apply safe damping step
        dw = -0.5_dp * f2
        dv2 = -0.1_dp * f1 / max(df1_dv2, 1.0e-12_dp)
      else
        dw  = (df1_dv2 * f2 - df2_dv2 * f1) / det_j
        dv2 = (df2_dw  * f1 - df1_dw  * f2) / det_j
      end if

      ! Backtracking line search / step limiter
      step_size = 1.0_dp
      if (abs(dw) > 0.5_dp * w) step_size = min(step_size, 0.5_dp * w / abs(dw))
      if (abs(dv2) > 0.2_dp) step_size = min(step_size, 0.2_dp / abs(dv2))

      w = max(w + step_size * dw, u%d)
      v2 = min(max(v2 + step_size * dv2, 0.0_dp), VEL_MAX * VEL_MAX)
    end do

    w_plus_b2 = w + b2
    gamma = 1.0_dp / sqrt(max(1.0_dp - v2, 1.0e-12_dp))
    gamma = min(gamma, GAMMA_MAX)

    rho = max(u%d / gamma, DHO_FLOOR)
    press = ((eos%gamma - 1.0_dp) / eos%gamma) * (w - gamma * u%d) / (gamma * gamma)
    press = max(press, PRESS_FLOOR)

    p%rho = rho
    p%p = press
    p%b = u%b
    do i = 1, 3
      p%v(i) = (s_up(i) + (s_dot_b / w) * u%b(i)) / w_plus_b2
    end do

    if (present(err)) call err%clear()
  end subroutine con2prim_2d_solve

end module umbra_con2prim_2d
