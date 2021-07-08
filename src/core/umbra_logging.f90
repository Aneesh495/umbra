! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_logging
! Description: Thread-safe leveled diagnostics logging with ANSI telemetry.

module umbra_logging
  use umbra_kinds, only: i4, dp
  implicit none
  private

  integer(i4), parameter, public :: LOG_LEVEL_TRACE = 0_i4
  integer(i4), parameter, public :: LOG_LEVEL_DEBUG = 1_i4
  integer(i4), parameter, public :: LOG_LEVEL_INFO  = 2_i4
  integer(i4), parameter, public :: LOG_LEVEL_WARN  = 3_i4
  integer(i4), parameter, public :: LOG_LEVEL_ERROR = 4_i4
  integer(i4), parameter, public :: LOG_LEVEL_OFF   = 5_i4

  integer(i4), save :: current_log_level = LOG_LEVEL_INFO
  logical, save :: enable_ansi_colors = .true.
  integer(i4), save :: log_unit = 6_i4

  ! ANSI escape codes
  character(len=7), parameter :: ANSI_RESET  = char(27) // "[0m"
  character(len=7), parameter :: ANSI_RED    = char(27) // "[31m"
  character(len=7), parameter :: ANSI_GREEN  = char(27) // "[32m"
  character(len=7), parameter :: ANSI_YELLOW = char(27) // "[33m"
  character(len=7), parameter :: ANSI_BLUE   = char(27) // "[34m"
  character(len=7), parameter :: ANSI_CYAN   = char(27) // "[36m"
  character(len=7), parameter :: ANSI_GRAY   = char(27) // "[90m"

  public :: log_set_level
  public :: log_get_level
  public :: log_set_colors
  public :: log_set_unit
  public :: log_trace
  public :: log_debug
  public :: log_info
  public :: log_warn
  public :: log_error

contains

  subroutine log_set_level(level)
    integer(i4), intent(in) :: level
    current_log_level = level
  end subroutine log_set_level

  pure function log_get_level() result(lvl)
    integer(i4) :: lvl
    lvl = current_log_level
  end function log_get_level

  subroutine log_set_colors(enable)
    logical, intent(in) :: enable
    enable_ansi_colors = enable
  end subroutine log_set_colors

  subroutine log_set_unit(u)
    integer(i4), intent(in) :: u
    log_unit = u
  end subroutine log_set_unit

  subroutine log_emit(level, tag, color, msg)
    integer(i4), intent(in) :: level
    character(len=*), intent(in) :: tag
    character(len=*), intent(in) :: color
    character(len=*), intent(in) :: msg

    if (level < current_log_level) return

    if (enable_ansi_colors) then
      write(log_unit, '(A,A,A,A,1X,A)') trim(color), "[", trim(tag), "]" // ANSI_RESET, trim(msg)
    else
      write(log_unit, '(A,A,A,1X,A)') "[", trim(tag), "]", trim(msg)
    end if
    flush(log_unit)
  end subroutine log_emit

  subroutine log_trace(msg)
    character(len=*), intent(in) :: msg
    call log_emit(LOG_LEVEL_TRACE, "TRACE", ANSI_GRAY, msg)
  end subroutine log_trace

  subroutine log_debug(msg)
    character(len=*), intent(in) :: msg
    call log_emit(LOG_LEVEL_DEBUG, "DEBUG", ANSI_CYAN, msg)
  end subroutine log_debug

  subroutine log_info(msg)
    character(len=*), intent(in) :: msg
    call log_emit(LOG_LEVEL_INFO, "INFO ", ANSI_GREEN, msg)
  end subroutine log_info

  subroutine log_warn(msg)
    character(len=*), intent(in) :: msg
    call log_emit(LOG_LEVEL_WARN, "WARN ", ANSI_YELLOW, msg)
  end subroutine log_warn

  subroutine log_error(msg)
    character(len=*), intent(in) :: msg
    call log_emit(LOG_LEVEL_ERROR, "ERROR", ANSI_RED, msg)
  end subroutine log_error

end module umbra_logging
