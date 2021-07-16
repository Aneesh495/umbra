! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_timer
! Description: High-resolution monotonic timers for profiling kernel execution phases.

module umbra_timer
  use umbra_kinds, only: dp, i8
  implicit none
  private

  type, public :: umbra_timer_t
    integer(i8) :: t_start = 0_i8
    integer(i8) :: t_stop = 0_i8
    integer(i8) :: rate = 1_i8
    real(dp)    :: elapsed_sec = 0.0_dp
    logical     :: is_running = .false.
  contains
    procedure, public :: start => timer_start
    procedure, public :: stop => timer_stop
    procedure, public :: reset => timer_reset
    procedure, public :: elapsed => timer_elapsed
  end type umbra_timer_t

contains

  subroutine timer_start(this)
    class(umbra_timer_t), intent(inout) :: this
    integer(i8) :: count, rate, count_max
    call system_clock(count, rate, count_max)
    this%t_start = count
    this%rate = rate
    this%is_running = .true.
  end subroutine timer_start

  subroutine timer_stop(this)
    class(umbra_timer_t), intent(inout) :: this
    integer(i8) :: count, rate, count_max
    if (.not. this%is_running) return
    call system_clock(count, rate, count_max)
    this%t_stop = count
    this%elapsed_sec = this%elapsed_sec + real(this%t_stop - this%t_start, dp) / real(this%rate, dp)
    this%is_running = .false.
  end subroutine timer_stop

  pure subroutine timer_reset(this)
    class(umbra_timer_t), intent(inout) :: this
    this%t_start = 0_i8
    this%t_stop = 0_i8
    this%elapsed_sec = 0.0_dp
    this%is_running = .false.
  end subroutine timer_reset

  pure function timer_elapsed(this) result(sec)
    class(umbra_timer_t), intent(in) :: this
    real(dp) :: sec
    sec = this%elapsed_sec
  end function timer_elapsed

end module umbra_timer
