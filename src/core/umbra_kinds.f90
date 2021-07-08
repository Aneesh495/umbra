! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_kinds
! Description: Architecture-neutral numeric precision kinds and interop standards.

module umbra_kinds
  use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64, real32, real64, real128
  use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_float, c_double, c_bool, c_ptr, c_null_ptr
  implicit none
  private

  ! Integer precision kinds
  integer, parameter, public :: i1 = int8
  integer, parameter, public :: i2 = int16
  integer, parameter, public :: i4 = int32
  integer, parameter, public :: i8 = int64
  integer, parameter, public :: ip = int64

  ! Floating-point precision kinds
  integer, parameter, public :: sp = real32
  integer, parameter, public :: dp = real64
  integer, parameter, public :: qp = real128
  integer, parameter, public :: wp = real64

  ! Logical kinds
  integer, parameter, public :: lk = kind(.true.)

  ! C interoperability types re-exported for foreign function boundaries
  public :: c_int, c_int64_t, c_float, c_double, c_bool, c_ptr, c_null_ptr

end module umbra_kinds
