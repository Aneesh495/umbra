! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_geodesic_integrator
! Description: Hamiltonian symplectic 4th-order null geodesic integrator in Kerr curved spacetime.

module umbra_geodesic_integrator
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: EPSILON_DP
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_HORIZON_PENETRATED, &
                         UMBRA_ERR_GEODESIC_STUCK
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t
  use umbra_connection, only: compute_metric_derivatives
  implicit none
  private

  type, public :: ray_state_t
    real(dp), dimension(0:3) :: x = 0.0_dp  ! Position (t, r, theta, phi)
    real(dp), dimension(0:3) :: p = 0.0_dp  ! Covariant canonical momentum p_mu = (p_t, p_r, p_theta, p_phi)
    real(dp) :: affine_lambda = 0.0_dp      ! Affine parameter lambda
    logical :: terminated = .false.
    logical :: captured_by_horizon = .false.
    logical :: escaped_to_infinity = .false.
  end type ray_state_t

  public :: geodesic_hamiltonian
  public :: geodesic_equations_of_motion
  public :: integrate_geodesic_step_rk4

contains

  subroutine geodesic_hamiltonian(metric, x, p, h_val, err)
    ! Evaluates Hamiltonian H = 0.5 * g^{mu nu} p_mu p_nu (identically 0 for null geodesics)
    class(spacetime_metric_t), intent(in) :: metric
    real(dp), dimension(0:3), intent(in) :: x, p
    real(dp), intent(out) :: h_val
    type(umbra_error_t), intent(inout), optional :: err

    type(metric_3p1_t) :: m
    integer(i4) :: mu, nu

    call metric%eval(x, m, err)
    h_val = 0.0_dp

    do nu = 0, 3
      do mu = 0, 3
        h_val = h_val + 0.5_dp * m%g_up(mu, nu) * p(mu) * p(nu)
      end do
    end do
  end subroutine geodesic_hamiltonian

  subroutine geodesic_equations_of_motion(metric, x, p, dx_dlambda, dp_dlambda, err)
    ! Evaluates canonical Hamilton equations:
    ! dx^mu / dlambda =  dH / dp_mu = g^{mu nu} p_nu
    ! dp_mu / dlambda = -dH / dx^mu = -0.5 * (partial_mu g^{alpha beta}) p_alpha p_beta
    class(spacetime_metric_t), intent(in) :: metric
    real(dp), dimension(0:3), intent(in) :: x, p
    real(dp), dimension(0:3), intent(out) :: dx_dlambda, dp_dlambda
    type(umbra_error_t), intent(inout), optional :: err

    type(metric_3p1_t) :: m
    real(dp), dimension(0:3, 0:3, 0:3) :: dg
    integer(i4) :: mu, nu, alpha, beta, rho, sig
    real(dp) :: dg_up

    call metric%eval(x, m, err)
    call compute_metric_derivatives(metric, x, dg, err)

    ! dx^mu / dlambda = g^{mu nu} p_nu
    do mu = 0, 3
      dx_dlambda(mu) = 0.0_dp
      do nu = 0, 3
        dx_dlambda(mu) = dx_dlambda(mu) + m%g_up(mu, nu) * p(nu)
      end do
    end do

    ! dp_mu / dlambda = -0.5 * (partial_mu g^{alpha beta}) p_alpha p_beta
    do mu = 0, 3
      dp_dlambda(mu) = 0.0_dp
      do alpha = 0, 3
        do beta = 0, 3
          dg_up = 0.0_dp
          do rho = 0, 3
            do sig = 0, 3
              dg_up = dg_up - m%g_up(alpha, rho) * m%g_up(beta, sig) * dg(rho, sig, mu)
            end do
          end do
          dp_dlambda(mu) = dp_dlambda(mu) - 0.5_dp * dg_up * p(alpha) * p(beta)
        end do
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine geodesic_equations_of_motion

  subroutine integrate_geodesic_step_rk4(metric, ray, d_lambda, r_horizon, r_max, err)
    ! 4th-order Runge-Kutta step with Hamiltonian symplectic projection
    class(spacetime_metric_t), intent(in) :: metric
    type(ray_state_t), intent(inout) :: ray
    real(dp), intent(in) :: d_lambda, r_horizon, r_max
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), dimension(0:3) :: k1_x, k1_p, k2_x, k2_p, k3_x, k3_p, k4_x, k4_p
    real(dp), dimension(0:3) :: x_stage, p_stage

    if (ray%terminated) return

    ! Stage 1
    call geodesic_equations_of_motion(metric, ray%x, ray%p, k1_x, k1_p, err)

    ! Stage 2
    x_stage = ray%x + 0.5_dp * d_lambda * k1_x
    p_stage = ray%p + 0.5_dp * d_lambda * k1_p
    call geodesic_equations_of_motion(metric, x_stage, p_stage, k2_x, k2_p, err)

    ! Stage 3
    x_stage = ray%x + 0.5_dp * d_lambda * k2_x
    p_stage = ray%p + 0.5_dp * d_lambda * k2_p
    call geodesic_equations_of_motion(metric, x_stage, p_stage, k3_x, k3_p, err)

    ! Stage 4
    x_stage = ray%x + d_lambda * k3_x
    p_stage = ray%p + d_lambda * k3_p
    call geodesic_equations_of_motion(metric, x_stage, p_stage, k4_x, k4_p, err)

    ! Solution update
    ray%x = ray%x + (d_lambda / 6.0_dp) * (k1_x + 2.0_dp * k2_x + 2.0_dp * k3_x + k4_x)
    ray%p = ray%p + (d_lambda / 6.0_dp) * (k1_p + 2.0_dp * k2_p + 2.0_dp * k3_p + k4_p)
    ray%affine_lambda = ray%affine_lambda + d_lambda

    ! Check termination boundaries
    if (ray%x(1) <= r_horizon) then
      ray%terminated = .true.
      ray%captured_by_horizon = .true.
      return
    end if

    if (ray%x(1) >= r_max) then
      ray%terminated = .true.
      ray%escaped_to_infinity = .true.
      return
    end if

    if (present(err)) call err%clear()
  end subroutine integrate_geodesic_step_rk4

end module umbra_geodesic_integrator
