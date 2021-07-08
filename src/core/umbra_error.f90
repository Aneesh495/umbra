! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_error
! Description: Structured error codes, diagnostics, and defensive assertion primitives.

module umbra_error
  use umbra_kinds, only: i4, dp
  implicit none
  private

  ! Return codes
  integer(i4), parameter, public :: UMBRA_SUCCESS                  = 0_i4
  integer(i4), parameter, public :: UMBRA_ERR_GENERIC              = 1_i4
  integer(i4), parameter, public :: UMBRA_ERR_CON2PRIM_FAILED      = 2_i4
  integer(i4), parameter, public :: UMBRA_ERR_NON_PHYSICAL_STATE   = 3_i4
  integer(i4), parameter, public :: UMBRA_ERR_SUPERLUMINAL_VEL     = 4_i4
  integer(i4), parameter, public :: UMBRA_ERR_NEGATIVE_PRESSURE    = 5_i4
  integer(i4), parameter, public :: UMBRA_ERR_NEGATIVE_DENSITY     = 6_i4
  integer(i4), parameter, public :: UMBRA_ERR_SINGULAR_METRIC      = 7_i4
  integer(i4), parameter, public :: UMBRA_ERR_HORIZON_PENETRATED   = 8_i4
  integer(i4), parameter, public :: UMBRA_ERR_CFL_VIOLATION        = 9_i4
  integer(i4), parameter, public :: UMBRA_ERR_DIVERGENCE_EXCEEDED  = 10_i4
  integer(i4), parameter, public :: UMBRA_ERR_NAN_DETECTED         = 11_i4
  integer(i4), parameter, public :: UMBRA_ERR_QUADRATURE_ORDER     = 12_i4
  integer(i4), parameter, public :: UMBRA_ERR_IO_FAILURE           = 13_i4
  integer(i4), parameter, public :: UMBRA_ERR_BOUNDARY_CONDITION   = 14_i4
  integer(i4), parameter, public :: UMBRA_ERR_GEODESIC_STUCK       = 15_i4

  type, public :: umbra_error_t
    integer(i4) :: code = UMBRA_SUCCESS
    character(len=256) :: message = ""
    character(len=128) :: location = ""
  contains
    procedure, public :: ok => error_is_ok
    procedure, public :: failed => error_is_failed
    procedure, public :: set => error_set
    procedure, public :: clear => error_clear
    procedure, public :: print => error_print
  end type umbra_error_t

  public :: umbra_error_string
  public :: umbra_check_finite
  public :: umbra_assert

contains

  pure function error_is_ok(this) result(res)
    class(umbra_error_t), intent(in) :: this
    logical :: res
    res = (this%code == UMBRA_SUCCESS)
  end function error_is_ok

  pure function error_is_failed(this) result(res)
    class(umbra_error_t), intent(in) :: this
    logical :: res
    res = (this%code /= UMBRA_SUCCESS)
  end function error_is_failed

  subroutine error_set(this, code, msg, loc)
    class(umbra_error_t), intent(inout) :: this
    integer(i4), intent(in) :: code
    character(len=*), intent(in) :: msg
    character(len=*), intent(in), optional :: loc

    this%code = code
    this%message = trim(msg)
    if (present(loc)) then
      this%location = trim(loc)
    else
      this%location = ""
    end if
  end subroutine error_set

  pure subroutine error_clear(this)
    class(umbra_error_t), intent(inout) :: this
    this%code = UMBRA_SUCCESS
    this%message = ""
    this%location = ""
  end subroutine error_clear

  subroutine error_print(this, unit)
    class(umbra_error_t), intent(in) :: this
    integer(i4), intent(in), optional :: unit
    integer(i4) :: u

    u = 0
    if (present(unit)) u = unit

    if (this%failed()) then
      if (len_trim(this%location) > 0) then
        write(u, '(A,I0,A,A,A,A,A)') "[umbra::error] Code ", this%code, &
          " (", trim(umbra_error_string(this%code)), ") at ", trim(this%location), ": ", trim(this%message)
      else
        write(u, '(A,I0,A,A,A,A)') "[umbra::error] Code ", this%code, &
          " (", trim(umbra_error_string(this%code)), "): ", trim(this%message)
      end if
    end if
  end subroutine error_print

  pure function umbra_error_string(code) result(str)
    integer(i4), intent(in) :: code
    character(len=64) :: str

    select case (code)
      case (UMBRA_SUCCESS)
        str = "SUCCESS"
      case (UMBRA_ERR_GENERIC)
        str = "GENERIC_ERROR"
      case (UMBRA_ERR_CON2PRIM_FAILED)
        str = "CON2PRIM_FAILED"
      case (UMBRA_ERR_NON_PHYSICAL_STATE)
        str = "NON_PHYSICAL_STATE"
      case (UMBRA_ERR_SUPERLUMINAL_VEL)
        str = "SUPERLUMINAL_VELOCITY"
      case (UMBRA_ERR_NEGATIVE_PRESSURE)
        str = "NEGATIVE_PRESSURE"
      case (UMBRA_ERR_NEGATIVE_DENSITY)
        str = "NEGATIVE_DENSITY"
      case (UMBRA_ERR_SINGULAR_METRIC)
        str = "SINGULAR_METRIC"
      case (UMBRA_ERR_HORIZON_PENETRATED)
        str = "HORIZON_PENETRATED"
      case (UMBRA_ERR_CFL_VIOLATION)
        str = "CFL_VIOLATION"
      case (UMBRA_ERR_DIVERGENCE_EXCEEDED)
        str = "DIVERGENCE_EXCEEDED"
      case (UMBRA_ERR_NAN_DETECTED)
        str = "NAN_OR_INF_DETECTED"
      case (UMBRA_ERR_QUADRATURE_ORDER)
        str = "INVALID_QUADRATURE_ORDER"
      case (UMBRA_ERR_IO_FAILURE)
        str = "IO_FAILURE"
      case (UMBRA_ERR_BOUNDARY_CONDITION)
        str = "INVALID_BOUNDARY_CONDITION"
      case (UMBRA_ERR_GEODESIC_STUCK)
        str = "GEODESIC_INTEGRATION_STUCK"
      case default
        str = "UNKNOWN_ERROR_CODE"
    end select
  end function umbra_error_string

  pure function umbra_check_finite(val) result(res)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(dp), intent(in) :: val
    logical :: res
    res = ieee_is_finite(val)
  end function umbra_check_finite

  subroutine umbra_assert(condition, message, err)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    type(umbra_error_t), intent(inout), optional :: err

    if (.not. condition) then
      if (present(err)) then
        call err%set(UMBRA_ERR_GENERIC, message)
      else
        write(0, '(A,A)') "[FATAL ASSERTION FAILURE] ", trim(message)
        error stop 1
      end if
    end if
  end subroutine umbra_assert

end module umbra_error
