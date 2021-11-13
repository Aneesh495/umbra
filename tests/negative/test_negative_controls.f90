!> umbra: Rigorous negative control and adversarial stress test suite
!> Validates defensive barriers, unphysical state rejections, and error handling
program test_negative_controls
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use umbra_kinds, only: dp, i4

  use umbra_error, only: umbra_error_t, UMBRA_ERR_NEGATIVE_DENSITY, &
                         UMBRA_ERR_NEGATIVE_PRESSURE, UMBRA_ERR_CON2PRIM_FAILED, &
                         UMBRA_ERR_SINGULAR_METRIC, UMBRA_ERR_IO_FAILURE, &
                         umbra_check_finite
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  use umbra_metric_base, only: metric_3p1_t
  use umbra_minkowski, only: minkowski_metric_t
  use umbra_eos_ideal, only: ideal_gas_eos_t, make_ideal_gas_eos
  use umbra_con2prim_1d, only: con2prim_1d_solve
  use umbra_math, only: mat3_inv
  use umbra_checkpoint, only: checkpoint_meta_t, write_checkpoint_2d, read_checkpoint_2d
  implicit none

  integer(i4) :: passed_count, total_count
  type(minkowski_metric_t) :: mink
  type(metric_3p1_t) :: m
  type(ideal_gas_eos_t) :: eos
  type(prim_state_t) :: p
  type(cons_state_t) :: u
  type(umbra_error_t) :: err
  real(dp), dimension(3,3) :: a_sing, a_inv
  real(dp) :: det_val, nan_val
  type(checkpoint_meta_t) :: meta
  real(dp), allocatable :: prim_arr(:,:,:), cons_arr(:,:,:), psi_arr(:,:)
  real(dp), allocatable :: prim_read(:,:,:), cons_read(:,:,:), psi_read(:,:)
  integer(i4) :: u_file

  passed_count = 0
  total_count = 0
  eos = make_ideal_gas_eos(1.4_dp)
  call mink%eval([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], m)

  write(*, '(A)') "============================================================"
  write(*, '(A)') "umbra: Kernel Adversarial Negative Control Suite"
  write(*, '(A)') "============================================================"

  ! NC-1: Negative density injection into prim_to_cons
  total_count = total_count + 1
  p%rho = -1.0_dp
  p%p = 1.0_dp
  p%v = 0.0_dp
  p%b = 0.0_dp
  call prim_to_cons(p, m, eos, u, err)
  if (err%failed() .and. err%code == UMBRA_ERR_NEGATIVE_DENSITY) then
    passed_count = passed_count + 1
    write(*, '(A)') "  [PASS] NC-01: Rejection of negative rest-mass density (rho < 0)"
  else
    write(*, '(A)') "  [FAIL] NC-01: Did not reject negative density"
  end if

  ! NC-2: Negative pressure injection into prim_to_cons
  total_count = total_count + 1
  p%rho = 1.0_dp
  p%p = -0.5_dp
  p%v = 0.0_dp
  p%b = 0.0_dp
  call prim_to_cons(p, m, eos, u, err)
  if (err%failed() .and. err%code == UMBRA_ERR_NEGATIVE_PRESSURE) then
    passed_count = passed_count + 1
    write(*, '(A)') "  [PASS] NC-02: Rejection of negative gas pressure (p < 0)"
  else
    write(*, '(A)') "  [FAIL] NC-02: Did not reject negative pressure"
  end if

  ! NC-3: Non-positive conserved D in con2prim_1d
  total_count = total_count + 1
  u%d = -2.0_dp
  u%tau = 1.0_dp
  u%s = 0.0_dp
  u%b = 0.0_dp
  call con2prim_1d_solve(u, m, eos, p, err)
  if (err%failed() .and. err%code == UMBRA_ERR_CON2PRIM_FAILED) then
    passed_count = passed_count + 1
    write(*, '(A)') "  [PASS] NC-03: Rejection of non-positive conserved D in 1D inversion"
  else
    write(*, '(A)') "  [FAIL] NC-03: Did not reject non-positive D"
  end if

  ! NC-4: Singular matrix inversion defense
  total_count = total_count + 1
  a_sing = 0.0_dp ! Determinant identically 0
  call mat3_inv(a_sing, a_inv, det_val, err)
  if (err%failed() .and. err%code == UMBRA_ERR_SINGULAR_METRIC) then
    passed_count = passed_count + 1
    write(*, '(A)') "  [PASS] NC-04: Detection and rejection of singular 3x3 metric matrix"
  else
    write(*, '(A)') "  [FAIL] NC-04: Did not reject singular matrix"
  end if

  ! NC-5: IEEE-754 NaN detection
  total_count = total_count + 1
  nan_val = ieee_value(1.0_dp, ieee_quiet_nan)
  if (.not. umbra_check_finite(nan_val)) then
    passed_count = passed_count + 1
    write(*, '(A)') "  [PASS] NC-05: Trap and detection of unphysical IEEE-754 NaN scalar"
  else
    write(*, '(A)') "  [FAIL] NC-05: Failed to trap NaN scalar"
  end if

  ! NC-6: Checkpoint corrupted magic header detection
  total_count = total_count + 1
  meta%nx = 4
  meta%ny = 4
  allocate(prim_arr(8, 4, 4), cons_arr(8, 4, 4), psi_arr(4, 4))
  prim_arr = 1.0_dp
  cons_arr = 1.0_dp
  psi_arr = 0.0_dp
  call write_checkpoint_2d("corrupt_test.chk", meta, prim_arr, cons_arr, psi_arr, err)

  ! Tamper magic bytes
  open(newunit=u_file, file="corrupt_test.chk", access='stream', form='unformatted', status='old')
  write(u_file) "BADMAGIC"
  close(u_file)

  call read_checkpoint_2d("corrupt_test.chk", meta, prim_read, cons_read, psi_read, err)
  if (err%failed() .and. err%code == UMBRA_ERR_IO_FAILURE) then
    passed_count = passed_count + 1
    write(*, '(A)') "  [PASS] NC-06: Detection and rejection of corrupted checkpoint magic header"
  else
    write(*, '(A)') "  [FAIL] NC-06: Failed to catch corrupt checkpoint magic"
  end if

  ! Clean up
  open(newunit=u_file, file="corrupt_test.chk", status='old')
  close(u_file, status='delete')
  deallocate(prim_arr, cons_arr, psi_arr)

  write(*, '(A)') "============================================================"
  write(*, '(A, I0, A, I0, A)') "NEGATIVE CONTROLS RESULT: ", passed_count, " / ", total_count, " PASSED"
  if (passed_count == total_count) then
    write(*, '(A)') ">>> ALL KERNEL NEGATIVE CONTROLS PASSED <<<"
    write(*, '(A)') "============================================================"
  else
    write(*, '(A)') ">>> ERROR: ONE OR MORE NEGATIVE CONTROLS FAILED <<<"
    write(*, '(A)') "============================================================"
    stop 1
  end if

end program test_negative_controls
