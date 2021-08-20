! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_sources
! Description: Curved spacetime geometric source terms for 3+1 ADM relativistic MHD evolution.

module umbra_sources
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS, CONS_D, CONS_SX, CONS_SY, CONS_SZ, CONS_TAU, &
                             CONS_BX, CONS_BY, CONS_BZ
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  use umbra_fluxes, only: compute_stress_energy_tensor
  implicit none
  private

  public :: compute_geometric_sources

contains

  subroutine compute_geometric_sources(metric, x, p, m, eos, source, err)
    ! Evaluates 3+1 ADM geometric source vector S_geom
    class(spacetime_metric_t), intent(in) :: metric
    real(dp), dimension(0:3), intent(in) :: x
    type(prim_state_t), intent(in) :: p
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    real(dp), dimension(NCONS), intent(out) :: source
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), parameter :: h = 1.0e-5_dp
    real(dp), parameter :: inv_2h = 0.5_dp / h
    type(cons_state_t) :: u
    real(dp), dimension(0:3, 0:3) :: t_down, t_up
    real(dp), dimension(3) :: grad_alpha
    real(dp), dimension(3, 3) :: grad_beta
    real(dp), dimension(3, 3, 3) :: dgamma
    real(dp), dimension(0:3) :: xp, xm
    type(metric_3p1_t) :: mp, mm
    integer(i4) :: j, i, k

    source = 0.0_dp
    call prim_to_cons(p, m, eos, u, err)
    call compute_stress_energy_tensor(p, m, eos, t_down, t_up)

    ! Compute gradients of lapse alpha, shift beta^i, and spatial metric gamma_ik
    do j = 1, 3
      xp = x
      xm = x
      xp(j) = x(j) + h
      xm(j) = x(j) - h
      call metric%eval(xp, mp, err)
      call metric%eval(xm, mm, err)

      grad_alpha(j) = (mp%alpha - mm%alpha) * inv_2h
      do i = 1, 3
        grad_beta(i, j) = (mp%beta_up(i) - mm%beta_up(i)) * inv_2h
        do k = 1, 3
          dgamma(i, k, j) = (mp%gamma_down(i, k) - mm%gamma_down(i, k)) * inv_2h
        end do
      end do
    end do

    ! Source term for mass density: S_D = 0
    source(CONS_D) = 0.0_dp

    ! Source term for momentum density: S_{S_j}
    do j = 1, 3
      source(CONS_SX + j - 1) = -(u%tau + u%d) * grad_alpha(j)
      do i = 1, 3
        source(CONS_SX + j - 1) = source(CONS_SX + j - 1) + u%s(i) * grad_beta(i, j)
        do k = 1, 3
          source(CONS_SX + j - 1) = source(CONS_SX + j - 1) + 0.5_dp * m%alpha * t_up(i, k) * dgamma(i, k, j)
        end do
      end do
    end do

    ! Source term for energy density: S_tau
    do i = 1, 3
      source(CONS_TAU) = source(CONS_TAU) + (m%gamma_up(i, 1)*u%s(1) + m%gamma_up(i, 2)*u%s(2) + m%gamma_up(i, 3)*u%s(3)) &
                         * grad_alpha(i)
      do k = 1, 3
        do j = 1, 3
          source(CONS_TAU) = source(CONS_TAU) + 0.5_dp * t_up(i, k) * m%beta_up(j) * dgamma(i, k, j)
        end do
        source(CONS_TAU) = source(CONS_TAU) + t_up(i, k) * m%gamma_down(k, 1) * grad_beta(1, i)
        source(CONS_TAU) = source(CONS_TAU) + t_up(i, k) * m%gamma_down(k, 2) * grad_beta(2, i)
        source(CONS_TAU) = source(CONS_TAU) + t_up(i, k) * m%gamma_down(k, 3) * grad_beta(3, i)
      end do
    end do

    ! Source term for magnetic field: S_B = 0
    source(CONS_BX:CONS_BZ) = 0.0_dp

    if (present(err)) call err%clear()
  end subroutine compute_geometric_sources

end module umbra_sources
