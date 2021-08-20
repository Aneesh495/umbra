! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_fluxes
! Description: Relativistic magnetohydrodynamic flux tensors in 3+1 ADM conservative form.

module umbra_fluxes
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS, CONS_D, CONS_SX, CONS_SY, CONS_SZ, CONS_TAU, &
                             CONS_BX, CONS_BY, CONS_BZ
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  implicit none
  private

  public :: compute_mhd_flux
  public :: compute_stress_energy_tensor

contains

  subroutine compute_mhd_flux(p, m, eos, dir, flux, err)
    ! Computes flux vector F^k along spatial coordinate direction dir in {1, 2, 3}
    type(prim_state_t), intent(in) :: p
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: dir
    real(dp), dimension(NCONS), intent(out) :: flux
    type(umbra_error_t), intent(inout), optional :: err

    type(cons_state_t) :: u
    real(dp) :: gamma, p_tot, b_0, v_rel
    real(dp), dimension(3) :: b_i, b_down
    integer(i4) :: j, k

    flux = 0.0_dp
    call prim_to_cons(p, m, eos, u, err)

    gamma = p%lorentz_factor(m)
    p_tot = p%total_pressure(m)
    call p%magnetic_field_4vec(m, b_0, b_i)

    ! Relative transport velocity: v^k - beta^k / alpha
    v_rel = p%v(dir) - m%beta_up(dir) / m%alpha

    ! Covariant b_j = gamma_jk b^k
    b_down = 0.0_dp
    do k = 1, 3
      do j = 1, 3
        b_down(j) = b_down(j) + m%gamma_down(j, k) * b_i(k)
      end do
    end do

    ! F^k_D = D (v^k - beta^k / alpha)
    flux(CONS_D) = u%d * v_rel

    ! F^k_{S_j} = S_j (v^k - beta^k / alpha) + p^* delta_j^k - (b_j B^k) / gamma
    do j = 1, 3
      flux(CONS_SX + j - 1) = u%s(j) * v_rel - (b_down(j) * p%b(dir)) / gamma
      if (j == dir) then
        flux(CONS_SX + j - 1) = flux(CONS_SX + j - 1) + p_tot
      end if
    end do

    ! F^k_tau = tau (v^k - beta^k / alpha) + p^* v^k - (alpha b^0 B^k) / gamma
    flux(CONS_TAU) = u%tau * v_rel + p_tot * p%v(dir) - (m%alpha * b_0 * p%b(dir)) / gamma

    ! F^k_{B^j} = B^j (v^k - beta^k / alpha) - B^k (v^j - beta^j / alpha)
    do j = 1, 3
      flux(CONS_BX + j - 1) = p%b(j) * v_rel - p%b(dir) * (p%v(j) - m%beta_up(j) / m%alpha)
    end do

    if (present(err)) call err%clear()
  end subroutine compute_mhd_flux

  pure subroutine compute_stress_energy_tensor(p, m, eos, t_down, t_up)
    ! Computes full 4D stress-energy tensor T_{mu nu} = (rho h + b^2) u_mu u_nu + (p + 0.5 b^2) g_{mu nu} - b_mu b_nu
    type(prim_state_t), intent(in) :: p
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    real(dp), dimension(0:3, 0:3), intent(out) :: t_down
    real(dp), dimension(0:3, 0:3), intent(out), optional :: t_up

    real(dp) :: gamma, eps, h_val, b2, b_0, p_tot, rho_h_tot
    real(dp), dimension(3) :: b_spatial
    real(dp), dimension(0:3) :: u_up, u_down, b_4up, b_4down
    integer(i4) :: mu, nu, i

    gamma = p%lorentz_factor(m)
    eps = eos%specific_internal_energy(p%rho, p%p)
    h_val = eos%specific_enthalpy(p%rho, eps, p%p)
    b2 = p%magnetic_energy(m)
    p_tot = p%p + 0.5_dp * b2
    rho_h_tot = p%rho * h_val + b2

    ! 4-velocity contravariant u^mu = (gamma / alpha, gamma (v^i - beta^i / alpha))
    u_up(0) = gamma / m%alpha
    do i = 1, 3
      u_up(i) = gamma * (p%v(i) - m%beta_up(i) / m%alpha)
    end do

    ! 4-velocity covariant u_mu = g_{mu nu} u^nu
    do mu = 0, 3
      u_down(mu) = 0.0_dp
      do nu = 0, 3
        u_down(mu) = u_down(mu) + m%g_down(mu, nu) * u_up(nu)
      end do
    end do

    ! Magnetic 4-vector b^mu
    call p%magnetic_field_4vec(m, b_0, b_spatial)
    b_4up(0) = b_0
    b_4up(1:3) = b_spatial

    ! b_mu = g_{mu nu} b^nu
    do mu = 0, 3
      b_4down(mu) = 0.0_dp
      do nu = 0, 3
        b_4down(mu) = b_4down(mu) + m%g_down(mu, nu) * b_4up(nu)
      end do
    end do

    ! T_{mu nu} = (rho h + b^2) u_mu u_nu + p^* g_{mu nu} - b_mu b_nu
    do nu = 0, 3
      do mu = 0, 3
        t_down(mu, nu) = rho_h_tot * u_down(mu) * u_down(nu) + p_tot * m%g_down(mu, nu) - b_4down(mu) * b_4down(nu)
      end do
    end do

    ! Contravariant T^{mu nu} if requested
    if (present(t_up)) then
      do nu = 0, 3
        do mu = 0, 3
          t_up(mu, nu) = rho_h_tot * u_up(mu) * u_up(nu) + p_tot * m%g_up(mu, nu) - b_4up(mu) * b_4up(nu)
        end do
      end do
    end if
  end subroutine compute_stress_energy_tensor

end module umbra_fluxes
