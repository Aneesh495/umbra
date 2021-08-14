! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_con2prim_1d
! Description: One-dimensional master variable Newton-Raphson conservative-to-primitive solver.

module umbra_con2prim_1d
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: DHO_FLOOR, PRESS_FLOOR, GAMMA_MAX, VEL_MAX, TOL_CONVERGENCE
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_CON2PRIM_FAILED
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t
  use umbra_math, only: solve_newton_bisection
  implicit none
  private

  public :: con2prim_1d_solve

  ! Thread-local parameters for the 1D scalar function residual
  real(dp) :: g_d, g_tau, g_b2, g_s2, g_s_dot_b, g_gamma_eos
  type(metric_3p1_t) :: g_metric

contains

  subroutine con2prim_1d_residual_and_deriv(w, f, df)
    real(dp), intent(in) :: w
    real(dp), intent(out) :: f, df

    real(dp) :: w_sq, w_plus_b2, v_sq, gamma, rho, p
    real(dp) :: term_num, term_den, dv2_dw, dgamma_dw, dp_dw
    real(dp) :: s_dot_b_sq

    w_sq = w * w
    w_plus_b2 = w + g_b2
    s_dot_b_sq = g_s_dot_b * g_s_dot_b

    ! Velocity squared: v^2 = [S^2 W^2 + 2 (S.B)^2 W + (S.B)^2 B^2] / [W^2 (W + B^2)^2]
    term_num = g_s2 * w_sq + 2.0_dp * s_dot_b_sq * w + s_dot_b_sq * g_b2
    term_den = w_sq * (w_plus_b2 * w_plus_b2)

    v_sq = term_num / max(term_den, 1.0e-30_dp)
    v_sq = min(max(v_sq, 0.0_dp), VEL_MAX * VEL_MAX)

    gamma = 1.0_dp / sqrt(max(1.0_dp - v_sq, 1.0e-12_dp))
    gamma = min(gamma, GAMMA_MAX)

    rho = max(g_d / gamma, DHO_FLOOR)

    ! Gas pressure for ideal gas: p = (gamma_eos - 1) * (W - gamma * D) / (gamma_eos * gamma^2)
    p = ((g_gamma_eos - 1.0_dp) / g_gamma_eos) * (w - gamma * g_d) / (gamma * gamma)
    p = max(p, PRESS_FLOOR)

    ! Residual: f(W) = W - p + 0.5 * B^2 * (1 + v^2) - 0.5 * (S.B)^2 / W^2 - D - tau
    f = w - p + 0.5_dp * g_b2 * (1.0_dp + v_sq) - 0.5_dp * s_dot_b_sq / max(w_sq, 1.0e-30_dp) - g_d - g_tau

    ! Derivative df/dW via finite difference step
    dv2_dw = -2.0_dp * (g_s2 * w_sq * w_plus_b2 + s_dot_b_sq * (3.0_dp * w * w_plus_b2 + g_b2 * w_plus_b2 + w * g_b2)) &
             / max(w_sq * w * (w_plus_b2**3), 1.0e-30_dp)
    dgamma_dw = 0.5_dp * (gamma**3) * dv2_dw
    dp_dw = ((g_gamma_eos - 1.0_dp) / g_gamma_eos) * (1.0_dp / (gamma * gamma) - 2.0_dp * (w - gamma * g_d) * dgamma_dw / (gamma**3) &
            - g_d * dgamma_dw / (gamma * gamma))

    df = 1.0_dp - dp_dw + 0.5_dp * g_b2 * dv2_dw + s_dot_b_sq / max(w_sq * w, 1.0e-30_dp)
  end subroutine con2prim_1d_residual_and_deriv

  subroutine con2prim_1d_solve(u, m, eos, p, err)
    type(cons_state_t), intent(in) :: u
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    type(prim_state_t), intent(out) :: p
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: w_min, w_max, w_root
    real(dp) :: gamma, rho, press, v_sq, w_plus_b2
    real(dp), dimension(3) :: s_up
    integer(i4) :: i, j

    ! Check positivity of density
    if (u%d <= 0.0_dp) then
      if (present(err)) call err%set(UMBRA_ERR_CON2PRIM_FAILED, "Non-positive D in con2prim_1d", "con2prim_1d_solve")
      return
    end if

    ! Compute metric contractions
    g_d = u%d
    g_tau = u%tau
    g_gamma_eos = eos%gamma
    g_metric = m

    ! B^2 = gamma_ij B^i B^j
    g_b2 = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        g_b2 = g_b2 + m%gamma_down(i, j) * u%b(i) * u%b(j)
      end do
    end do

    ! S^2 = gamma^ij S_i S_j
    g_s2 = 0.0_dp
    s_up = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        s_up(i) = s_up(i) + m%gamma_up(i, j) * u%s(j)
        g_s2 = g_s2 + m%gamma_up(i, j) * u%s(i) * u%s(j)
      end do
    end do

    ! S . B = S_i B^i
    g_s_dot_b = u%s(1)*u%b(1) + u%s(2)*u%b(2) + u%s(3)*u%b(3)

    ! Initial bracket bounds for W = gamma^2 * rho * h
    ! W >= D is a strict mathematical requirement
    w_min = max(g_d, sqrt(g_s2) + 1.0e-8_dp)
    w_max = max(100.0_dp * (g_tau + g_d + g_b2), 10.0_dp * w_min)

    call solve_newton_bisection(con2prim_1d_residual_and_deriv, w_min, w_max, w_root, 60_i4, TOL_CONVERGENCE, err)

    if (present(err)) then
      if (err%failed()) return
    end if

    w_plus_b2 = w_root + g_b2
    v_sq = (g_s2 * (w_root**2) + 2.0_dp * (g_s_dot_b**2) * w_root + (g_s_dot_b**2) * g_b2) / &
           max((w_root**2) * (w_plus_b2**2), 1.0e-30_dp)
    v_sq = min(max(v_sq, 0.0_dp), VEL_MAX * VEL_MAX)

    gamma = 1.0_dp / sqrt(max(1.0_dp - v_sq, 1.0e-12_dp))
    gamma = min(gamma, GAMMA_MAX)

    rho = max(g_d / gamma, DHO_FLOOR)
    press = ((g_gamma_eos - 1.0_dp) / g_gamma_eos) * (w_root - gamma * g_d) / (gamma * gamma)
    press = max(press, PRESS_FLOOR)

    ! Recover spatial 3-velocity v^i = (S^i + (S.B / W) B^i) / (W + B^2)
    p%rho = rho
    p%p = press
    p%b = u%b
    do i = 1, 3
      p%v(i) = (s_up(i) + (g_s_dot_b / w_root) * u%b(i)) / w_plus_b2
    end do

    if (present(err)) call err%clear()
  end subroutine con2prim_1d_solve

end module umbra_con2prim_1d
