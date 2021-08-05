! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_eos_ideal
! Description: Standard Gamma-law ideal gas equation of state for relativistic gas dynamics.

module umbra_eos_ideal
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: DHO_FLOOR, PRESS_FLOOR, TEMPERATURE_FLOOR
  use umbra_eos_base, only: eos_t
  implicit none
  private

  type, extends(eos_t), public :: ideal_gas_eos_t
  contains
    procedure, public :: pressure => ideal_pressure
    procedure, public :: specific_internal_energy => ideal_eps
    procedure, public :: sound_speed_sq => ideal_sound_speed_sq
    procedure, public :: specific_enthalpy => ideal_enthalpy
    procedure, public :: temperature => ideal_temperature
    procedure, public :: dpress_drho => ideal_dpress_drho
    procedure, public :: dpress_deps => ideal_dpress_deps
  end type ideal_gas_eos_t

  public :: make_ideal_gas_eos

contains

  function make_ideal_gas_eos(gamma_val) result(eos)
    real(dp), intent(in) :: gamma_val
    type(ideal_gas_eos_t) :: eos
    eos%gamma = gamma_val
  end function make_ideal_gas_eos

  pure function ideal_pressure(this, rho, eps) result(p)
    class(ideal_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: p

    p = (this%gamma - 1.0_dp) * max(rho, DHO_FLOOR) * max(eps, 0.0_dp)
    p = max(p, PRESS_FLOOR)
  end function ideal_pressure

  pure function ideal_eps(this, rho, p) result(eps)
    class(ideal_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, p
    real(dp) :: eps

    eps = max(p, PRESS_FLOOR) / ((this%gamma - 1.0_dp) * max(rho, DHO_FLOOR))
  end function ideal_eps

  pure function ideal_enthalpy(this, rho, eps, p) result(h)
    class(ideal_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: h, eff_p, eff_rho

    eff_p = max(p, PRESS_FLOOR)
    eff_rho = max(rho, DHO_FLOOR)
    h = 1.0_dp + (this%gamma / (this%gamma - 1.0_dp)) * (eff_p / eff_rho) + 0.0_dp * eps
  end function ideal_enthalpy

  pure function ideal_sound_speed_sq(this, rho, eps, p) result(cs2)
    class(ideal_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: cs2, h_val, eff_p, eff_rho

    eff_p = max(p, PRESS_FLOOR)
    eff_rho = max(rho, DHO_FLOOR)
    h_val = this%specific_enthalpy(eff_rho, eps, eff_p)

    cs2 = (this%gamma * eff_p) / (eff_rho * h_val)
    ! Sound speed must not exceed speed of light (cs2 < 1 in geometric units)
    cs2 = min(max(cs2, 1.0e-12_dp), 0.999999_dp)
  end function ideal_sound_speed_sq

  pure function ideal_temperature(this, rho, eps, p) result(theta)
    class(ideal_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: theta

    theta = max(p, PRESS_FLOOR) / max(rho, DHO_FLOOR) + 0.0_dp * eps
    theta = max(theta, TEMPERATURE_FLOOR)
  end function ideal_temperature

  pure function ideal_dpress_drho(this, rho, eps) result(dp_drho)
    class(ideal_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: dp_drho

    dp_drho = (this%gamma - 1.0_dp) * max(eps, 0.0_dp) + 0.0_dp * rho
  end function ideal_dpress_drho

  pure function ideal_dpress_deps(this, rho, eps) result(dp_deps)
    class(ideal_gas_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: dp_deps

    dp_deps = (this%gamma - 1.0_dp) * max(rho, DHO_FLOOR) + 0.0_dp * eps
  end function ideal_dpress_deps

end module umbra_eos_ideal
