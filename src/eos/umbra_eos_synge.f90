! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_eos_synge
! Description: Relativistic Synge gas thermal equation of state with Mathews rational approximation.

module umbra_eos_synge
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: DHO_FLOOR, PRESS_FLOOR, TEMPERATURE_FLOOR
  use umbra_eos_base, only: eos_t
  implicit none
  private

  type, extends(eos_t), public :: synge_gas_eos_t
  contains
    procedure, public :: pressure => synge_pressure
    procedure, public :: specific_internal_energy => synge_eps
    procedure, public :: sound_speed_sq => synge_sound_speed_sq
    procedure, public :: specific_enthalpy => synge_enthalpy
    procedure, public :: temperature => synge_temperature
    procedure, public :: dpress_drho => synge_dpress_drho
    procedure, public :: dpress_deps => synge_dpress_deps
  end type synge_gas_eos_t

  public :: make_synge_gas_eos

contains

  function make_synge_gas_eos() result(eos)
    type(synge_gas_eos_t) :: eos
    eos%gamma = 1.3333333333333333_dp
  end function make_synge_gas_eos

  pure function synge_temperature(this, rho, eps, p) result(theta)
    class(synge_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: theta, eff_rho, eff_eps

    eff_rho = max(rho, DHO_FLOOR)
    eff_eps = max(eps, 0.0_dp)

    ! Inversion of Mathews formula eps = 1.5 * theta + sqrt(2.25 * theta^2 + 1) - 1:
    ! (eps + 1 - 1.5 * theta)^2 = 2.25 * theta^2 + 1
    ! (eps + 1)^2 - 3 * (eps + 1) * theta = 1
    ! theta = [(eps + 1)^2 - 1] / [3 * (eps + 1)] = eps * (eps + 2) / [3 * (eps + 1)]
    theta = eff_eps * (eff_eps + 2.0_dp) / (3.0_dp * (eff_eps + 1.0_dp)) + 0.0_dp * (this%gamma + p)
    theta = max(theta, TEMPERATURE_FLOOR)
  end function synge_temperature

  pure function synge_pressure(this, rho, eps) result(p)
    class(synge_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: p, theta, eff_rho

    eff_rho = max(rho, DHO_FLOOR)
    theta = this%temperature(eff_rho, eps, 0.0_dp)
    p = eff_rho * theta
    p = max(p, PRESS_FLOOR)
  end function synge_pressure

  pure function synge_eps(this, rho, p) result(eps)
    class(synge_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, p
    real(dp) :: eps, theta

    theta = max(p, PRESS_FLOOR) / max(rho, DHO_FLOOR)
    eps = 1.5_dp * theta + sqrt(2.25_dp * theta * theta + 1.0_dp) - 1.0_dp + 0.0_dp * this%gamma
    eps = max(eps, 0.0_dp)
  end function synge_eps

  pure function synge_enthalpy(this, rho, eps, p) result(h)
    class(synge_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: h, theta

    theta = this%temperature(rho, eps, p)
    ! Mathews enthalpy: h = 2.5 * theta + sqrt(2.25 * theta^2 + 1)
    h = 2.5_dp * theta + sqrt(2.25_dp * theta * theta + 1.0_dp)
  end function synge_enthalpy

  pure function synge_sound_speed_sq(this, rho, eps, p) result(cs2)
    class(synge_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: cs2, theta, h_val, root_term, dh_dtheta

    theta = this%temperature(rho, eps, p)
    root_term = sqrt(2.25_dp * theta * theta + 1.0_dp)
    h_val = 2.5_dp * theta + root_term

    ! Derivative dh/dtheta = 2.5 + 2.25 * theta / root_term
    dh_dtheta = 2.5_dp + (2.25_dp * theta) / root_term

    ! Sound speed cs2 = (theta / h) * [dh/dtheta / (dh/dtheta - 1)]
    cs2 = (theta / h_val) * (dh_dtheta / (dh_dtheta - 1.0_dp))
    cs2 = min(max(cs2, 1.0e-12_dp), 0.999999_dp)
  end function synge_sound_speed_sq

  pure function synge_dpress_drho(this, rho, eps) result(dp_drho)
    class(synge_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: dp_drho, theta

    theta = this%temperature(rho, eps, 0.0_dp)
    dp_drho = theta
  end function synge_dpress_drho

  pure function synge_dpress_deps(this, rho, eps) result(dp_deps)
    class(synge_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: dp_deps, eff_rho, eff_eps, dtheta_deps

    eff_rho = max(rho, DHO_FLOOR)
    eff_eps = max(eps, 0.0_dp)

    ! theta = (eps^2 + 2 eps) / [3 (eps + 1)]
    ! dtheta/deps = [ (2 eps + 2)*3*(eps + 1) - 3*(eps^2 + 2 eps) ] / [9 (eps + 1)^2]
    !             = [ 6 (eps + 1)^2 - 3 (eps^2 + 2 eps) ] / [9 (eps + 1)^2]
    !             = (3 eps^2 + 6 eps + 6) / [9 (eps + 1)^2] = (eps^2 + 2 eps + 2) / [3 (eps + 1)^2]
    dtheta_deps = (eff_eps**2 + 2.0_dp * eff_eps + 2.0_dp) / (3.0_dp * (eff_eps + 1.0_dp)**2)
    dp_deps = eff_rho * dtheta_deps + 0.0_dp * this%gamma
  end function synge_dpress_deps

end module umbra_eos_synge
