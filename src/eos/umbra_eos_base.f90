! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_eos_base
! Description: Polymorphic base class for fluid equations of state.

module umbra_eos_base
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: DHO_FLOOR, PRESS_FLOOR, TEMPERATURE_FLOOR
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  type, abstract, public :: eos_t
    real(dp) :: gamma = 1.3333333333333333_dp ! Default 4/3 for relativistic fluid
  contains
    procedure(eos_press_iface), deferred, public :: pressure
    procedure(eos_eps_iface), deferred, public :: specific_internal_energy
    procedure(eos_cs2_iface), deferred, public :: sound_speed_sq
    procedure(eos_enthalpy_iface), deferred, public :: specific_enthalpy
    procedure(eos_temp_iface), deferred, public :: temperature
    procedure(eos_dp_drho_iface), deferred, public :: dpress_drho
    procedure(eos_dp_deps_iface), deferred, public :: dpress_deps
  end type eos_t

  abstract interface
    pure function eos_press_iface(this, rho, eps) result(p)
      import :: eos_t, dp
      class(eos_t), intent(in) :: this
      real(dp), intent(in) :: rho, eps
      real(dp) :: p
    end function eos_press_iface

    pure function eos_eps_iface(this, rho, p) result(eps)
      import :: eos_t, dp
      class(eos_t), intent(in) :: this
      real(dp), intent(in) :: rho, p
      real(dp) :: eps
    end function eos_eps_iface

    pure function eos_cs2_iface(this, rho, eps, p) result(cs2)
      import :: eos_t, dp
      class(eos_t), intent(in) :: this
      real(dp), intent(in) :: rho, eps, p
      real(dp) :: cs2
    end function eos_cs2_iface

    pure function eos_enthalpy_iface(this, rho, eps, p) result(h)
      import :: eos_t, dp
      class(eos_t), intent(in) :: this
      real(dp), intent(in) :: rho, eps, p
      real(dp) :: h
    end function eos_enthalpy_iface

    pure function eos_temp_iface(this, rho, eps, p) result(theta)
      import :: eos_t, dp
      class(eos_t), intent(in) :: this
      real(dp), intent(in) :: rho, eps, p
      real(dp) :: theta
    end function eos_temp_iface

    pure function eos_dp_drho_iface(this, rho, eps) result(dp_drho)
      import :: eos_t, dp
      class(eos_t), intent(in) :: this
      real(dp), intent(in) :: rho, eps
      real(dp) :: dp_drho
    end function eos_dp_drho_iface

    pure function eos_dp_deps_iface(this, rho, eps) result(dp_deps)
      import :: eos_t, dp
      class(eos_t), intent(in) :: this
      real(dp), intent(in) :: rho, eps
      real(dp) :: dp_deps
    end function eos_dp_deps_iface
  end interface

end module umbra_eos_base
