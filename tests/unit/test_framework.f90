!> umbra: Deterministic unit testing framework
!> Lightweight assertion engine and TAP (Test Anything Protocol) reporter
module test_framework
  use umbra_kinds, only: dp, i4
  use umbra_error, only: umbra_error_t
  implicit none
  private

  public :: test_suite_t
  public :: assert_true
  public :: assert_false
  public :: assert_equal_i4
  public :: assert_equal_dp
  public :: assert_close_dp
  public :: assert_ok

  type, public :: test_suite_t
    character(len=64) :: name = "unnamed"
    integer(i4) :: num_tests = 0_i4
    integer(i4) :: num_passed = 0_i4
    integer(i4) :: num_failed = 0_i4
  contains
    procedure, public :: init => suite_init
    procedure, public :: report => suite_report
    procedure, public :: is_success => suite_is_success
  end type test_suite_t

contains

  subroutine suite_init(this, suite_name)
    class(test_suite_t), intent(inout) :: this
    character(len=*), intent(in) :: suite_name
    this%name = trim(suite_name)
    this%num_tests = 0_i4
    this%num_passed = 0_i4
    this%num_failed = 0_i4
  end subroutine suite_init

  subroutine suite_report(this)
    class(test_suite_t), intent(in) :: this
    write(*, '(A, A, A, I0, A, I0, A, I0, A)') "Suite '", trim(this%name), "': ", &
      this%num_tests, " tests (", this%num_passed, " passed, ", this%num_failed, " failed)"
  end subroutine suite_report

  pure function suite_is_success(this) result(res)
    class(test_suite_t), intent(in) :: this
    logical :: res
    res = (this%num_failed == 0_i4)
  end function suite_is_success

  subroutine assert_true(suite, condition, test_name)
    type(test_suite_t), intent(inout) :: suite
    logical, intent(in) :: condition
    character(len=*), intent(in) :: test_name

    suite%num_tests = suite%num_tests + 1_i4
    if (condition) then
      suite%num_passed = suite%num_passed + 1_i4
    else
      suite%num_failed = suite%num_failed + 1_i4
      write(*, '(A, A, A, A)') "  FAIL: [", trim(suite%name), "] ", trim(test_name)
    end if
  end subroutine assert_true

  subroutine assert_false(suite, condition, test_name)
    type(test_suite_t), intent(inout) :: suite
    logical, intent(in) :: condition
    character(len=*), intent(in) :: test_name

    call assert_true(suite, .not. condition, test_name)
  end subroutine assert_false

  subroutine assert_equal_i4(suite, actual, expected, test_name)
    type(test_suite_t), intent(inout) :: suite
    integer(i4), intent(in) :: actual, expected
    character(len=*), intent(in) :: test_name

    call assert_true(suite, actual == expected, test_name)
  end subroutine assert_equal_i4

  subroutine assert_equal_dp(suite, actual, expected, test_name)
    type(test_suite_t), intent(inout) :: suite
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: test_name

    real(dp) :: diff
    diff = abs(actual - expected)
    call assert_true(suite, diff <= tiny(1.0_dp), test_name)
  end subroutine assert_equal_dp

  subroutine assert_close_dp(suite, actual, expected, tol, test_name)
    type(test_suite_t), intent(inout) :: suite
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: test_name

    real(dp) :: diff
    diff = abs(actual - expected)
    call assert_true(suite, diff <= tol, test_name)
  end subroutine assert_close_dp

  subroutine assert_ok(suite, err, test_name)
    type(test_suite_t), intent(inout) :: suite
    type(umbra_error_t), intent(in) :: err
    character(len=*), intent(in) :: test_name

    call assert_true(suite, err%ok(), test_name)
  end subroutine assert_ok

end module test_framework
