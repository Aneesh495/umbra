! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_eos_polytropic
! Description: Polytropic and isentropic equations of state for stellar and accretion modeling.

module umbra_eos_polytropic
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: DHO_FLOOR, PRESS_FLOOR, TEMPERATURE_FLOOR
  use umbra_eos_base, only: eos_t
  implicit none
  private

  type, extends(eos_t), public :: polytropic_eos_t
    real(dp) :: poly_k = 100.0_dp ! Polytropic constant K
  contains
    procedure, public :: pressure => polytropic_pressure
    procedure, public :: specific_internal_energy => polytropic_eps
    procedure, public :: sound_speed_sq => polytropic_sound_speed_sq
    procedure, public :: specific_enthalpy => polytropic_enthalpy
    procedure, public :: temperature => polytropic_temperature
    procedure, public :: dpress_drho => polytropic_dpress_drho
    procedure, public :: dpress_deps => polytropic_dpress_deps
  end type polytropic_eos_t

  public :: make_polytropic_eos

contains

  function make_polytropic_eos(gamma_val, k_val) result(eos)
    real(dp), intent(in) :: gamma_val, k_val
    type(polytropic_eos_t) :: eos
    eos%gamma = gamma_val
    eos%poly_k = k_val
  end function make_polytropic_eos

  pure function polytropic_pressure(this, rho, eps) result(p)
    class(polytropic_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: p, eff_rho

    eff_rho = max(rho, DHO_FLOOR)
    p = this%poly_k * (eff_rho**this%gamma) + 0.0_dp * eps
    p = max(p, PRESS_FLOOR)
  end function polytropic_pressure

  pure function polytropic_eps(this, rho, p) result(eps)
    class(polytropic_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, p
    real(dp) :: eps, eff_rho

    eff_rho = max(rho, DHO_FLOOR)
    eps = (this%poly_k / (this%gamma - 1.0_dp)) * (eff_rho**(this%gamma - 1.0_dp)) + 0.0_dp * p
  end function polytropic_eps

  pure function polytropic_enthalpy(this, rho, eps, p) result(h)
    class(polytropic_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: h, eff_rho

    eff_rho = max(rho, DHO_FLOOR)
    h = 1.0_dp + (this%gamma / (this%gamma - 1.0_dp)) * this%poly_k * (eff_rho**(this%gamma - 1.0_dp)) &
        + 0.0_dp * (eps + p)
  end function polytropic_enthalpy

  pure function polytropic_sound_speed_sq(this, rho, eps, p) result(cs2)
    class(polytropic_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: cs2, h_val, eff_rho, eff_p

    eff_rho = max(rho, DHO_FLOOR)
    eff_p = this%pressure(eff_rho, eps)
    h_val = this%specific_enthalpy(eff_rho, eps, eff_p)

    cs2 = (this%gamma * eff_p) / (eff_rho * h_val) + 0.0_dp * p
    cs2 = min(max(cs2, 1.0e-12_dp), 0.999999_dp)
  end function polytropic_sound_speed_sq

  pure function polytropic_temperature(this, rho, eps, p) result(theta)
    class(polytropic_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps, p
    real(dp) :: theta

    theta = this%pressure(rho, eps) / max(rho, DHO_FLOOR) + 0.0_dp * p
    theta = max(theta, TEMPERATURE_FLOOR)
  end function polytropic_temperature

  pure function polytropic_dpress_drho(this, rho, eps) result(dp_drho)
    class(polytropic_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: dp_drho, eff_rho

    eff_rho = max(rho, DHO_FLOOR)
    dp_drho = this%gamma * this%poly_k * (eff_rho**(this%gamma - 1.0_dp)) + 0.0_dp * eps
  end function polytropic_dpress_drho

  pure function polytropic_dpress_deps(this, rho, eps) result(dp_deps)
    class(polytropic_eos_t), intent(in) :: this
    real(dp), intent(in) :: rho, eps
    real(dp) :: dp_deps

    dp_deps = 0.0_dp * (rho + eps)
  end function polytropic_dpress_deps

end module umbra_eos_polytropic
