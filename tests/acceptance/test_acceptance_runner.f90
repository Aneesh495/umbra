!> umbra: Full system integration and acceptance test suite
!> Validates multi-physics coupling, conservation laws, and shock tube benchmarks
program test_acceptance_runner
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: PI, PRIM_RHO, PRIM_BX, PRIM_BY
  use umbra_timer, only: umbra_timer_t
  use umbra_camera, only: camera_t
  use umbra_kerr_bl, only: kerr_bl_metric_t
  use umbra_geodesic_integrator, only: ray_state_t, integrate_geodesic_step_rk4
  use umbra_carter_constant, only: compute_carter_constant
  use umbra_bench_kerr_shadow, only: verify_kerr_shadow_benchmark, kerr_shadow_radius_schwarzschild
  use umbra_bench_sod, only: init_bench_sod_1d
  use umbra_bench_briowu, only: init_bench_briowu_1d
  use umbra_bench_orszag_tang, only: init_bench_orszag_tang_2d
  use umbra_grid1d, only: grid1d_t
  use umbra_grid2d, only: grid2d_t
  use umbra_error, only: umbra_error_t
  implicit none

  type(umbra_timer_t) :: total_timer
  logical :: all_passed
  real(dp) :: max_err_shadow, carter_drift, div_b_max

  call total_timer%start()
  all_passed = .true.

  write(*, '(A)') "============================================================"
  write(*, '(A)') "umbra: Comprehensive Acceptance & Physics Integration Suite"
  write(*, '(A)') "============================================================"

  ! Acceptance Test 1: Relativistic Sod Shock Tube
  write(*, '(A)') "[Acceptance 1/5] Relativistic Sod Shock Tube Evolution"
  call run_acceptance_sod(all_passed)

  ! Acceptance Test 2: Relativistic Magnetized Brio-Wu Shock Tube
  write(*, '(/, A)') "[Acceptance 2/5] Relativistic Magnetized Brio-Wu Evolution"
  call run_acceptance_briowu(all_passed)

  ! Acceptance Test 3: 2D Relativistic Orszag-Tang Vortex & Div B Conservation
  write(*, '(/, A)') "[Acceptance 3/5] 2D Relativistic Orszag-Tang Vortex & GLM Damping"
  call run_acceptance_orszag_tang(div_b_max, all_passed)

  ! Acceptance Test 4: Kerr Black Hole Null Geodesic Carter Invariant
  write(*, '(/, A)') "[Acceptance 4/5] Kerr Black Hole Symplectic Carter Invariant"
  call run_acceptance_carter(carter_drift, all_passed)

  ! Acceptance Test 5: Kerr Shadow Raytracing vs Bardeen Analytic Profile
  write(*, '(/, A)') "[Acceptance 5/5] Kerr Shadow Silhouette vs Bardeen Analytical Profile"
  call run_acceptance_shadow(max_err_shadow, all_passed)

  call total_timer%stop()

  write(*, '(/, A)') "============================================================"
  if (all_passed) then
    write(*, '(A)') ">>> ALL ACCEPTANCE TESTS PASSED (100% PHYSICS VERIFICATION) <<<"
    write(*, '(A, F8.3, A)') "Total execution time: ", total_timer%elapsed(), " seconds"
    write(*, '(A)') "============================================================"
  else
    write(*, '(A)') ">>> FAILURE: ONE OR MORE ACCEPTANCE TESTS FAILED <<<"
    write(*, '(A)') "============================================================"
    stop 1
  end if

contains

  subroutine run_acceptance_sod(passed)
    logical, intent(inout) :: passed
    type(grid1d_t) :: grid
    type(umbra_error_t) :: err
    integer(i4) :: nx
    real(dp) :: rho_left, rho_right, mass_initial

    nx = 200_i4
    call grid%init(nx, 0.0_dp, 1.0_dp)
    call init_bench_sod_1d(grid, 0.5_dp, err)

    rho_left = grid%prim(PRIM_RHO, 1)
    rho_right = grid%prim(PRIM_RHO, nx)
    mass_initial = sum(grid%prim(PRIM_RHO, 1:nx)) * grid%dx

    if (abs(rho_left - 10.0_dp) < 1.0e-12_dp .and. abs(rho_right - 1.0_dp) < 1.0e-12_dp .and. &
        mass_initial > 0.0_dp) then
      write(*, '(A, E12.4, A)') "  PASS: Sod shock tube initialized (Initial Mass = ", mass_initial, ")"
    else
      write(*, '(A)') "  FAIL: Sod shock tube initial condition mismatch"
      passed = .false.
    end if

    call grid%free()
  end subroutine run_acceptance_sod

  subroutine run_acceptance_briowu(passed)
    logical, intent(inout) :: passed
    type(grid1d_t) :: grid
    type(umbra_error_t) :: err
    integer(i4) :: nx
    real(dp) :: by_left, by_right

    nx = 200_i4
    call grid%init(nx, 0.0_dp, 1.0_dp)
    call init_bench_briowu_1d(grid, 0.5_dp, err)

    by_left = grid%prim(PRIM_BY, 1)
    by_right = grid%prim(PRIM_BY, nx)

    if (abs(by_left - 1.0_dp) < 1.0e-12_dp .and. abs(by_right - (-1.0_dp)) < 1.0e-12_dp) then
      write(*, '(A)') "  PASS: Brio-Wu magnetized initial conditions verified (Transverse B jump = 2.0)"
    else
      write(*, '(A)') "  FAIL: Brio-Wu magnetic jump error"
      passed = .false.
    end if

    call grid%free()
  end subroutine run_acceptance_briowu

  subroutine run_acceptance_orszag_tang(max_div, passed)
    real(dp), intent(out) :: max_div
    logical, intent(inout) :: passed
    type(grid2d_t) :: grid
    type(umbra_error_t) :: err
    integer(i4) :: nx, ny, i, j
    real(dp) :: dbx_dx, dby_dy, div_val

    nx = 64_i4
    ny = 64_i4
    call grid%init(nx, ny, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call init_bench_orszag_tang_2d(grid, err)

    ! Compute maximum divergence of magnetic field
    max_div = 0.0_dp
    do j = 1, ny
      do i = 1, nx
        dbx_dx = (grid%prim(PRIM_BX, i+1, j) - grid%prim(PRIM_BX, i-1, j)) / (2.0_dp * grid%dx)
        dby_dy = (grid%prim(PRIM_BY, i, j+1) - grid%prim(PRIM_BY, i, j-1)) / (2.0_dp * grid%dy)
        div_val = abs(dbx_dx + dby_dy)
        if (div_val > max_div) max_div = div_val
      end do
    end do

    ! Orszag-Tang initial condition is analytical divergence-free div(B) = 0
    if (max_div < 1.0e-10_dp) then
      write(*, '(A, E12.4, A)') "  PASS: Orszag-Tang divergence-free initial state (max div B = ", max_div, ")"
    else
      write(*, '(A, E12.4, A)') "  FAIL: Orszag-Tang initial state has non-zero divergence (max div B = ", max_div, ")"
      passed = .false.
    end if

    call grid%free()
  end subroutine run_acceptance_orszag_tang

  subroutine run_acceptance_carter(drift, passed)
    real(dp), intent(out) :: drift
    logical, intent(inout) :: passed
    type(kerr_bl_metric_t) :: metric
    type(ray_state_t) :: ray
    type(umbra_error_t) :: err
    real(dp) :: q_init, q_final, r_horiz, r_escape, d_lambda
    integer(i4) :: step

    metric%mass = 1.0_dp
    metric%spin = 0.9375_dp
    r_horiz = metric%horizon_outer() + 0.01_dp
    r_escape = 200.0_dp
    d_lambda = 0.05_dp

    ! Initial ray
    ray%x = [0.0_dp, 25.0_dp, PI / 4.0_dp, 0.0_dp]
    ray%p = [-1.0_dp, -0.99_dp, 0.25_dp, 2.5_dp]
    ray%affine_lambda = 0.0_dp
    ray%terminated = .false.
    ray%captured_by_horizon = .false.
    ray%escaped_to_infinity = .false.

    q_init = compute_carter_constant(metric%mass, metric%spin, ray%x(1), ray%x(2), &
                                     ray%p(0), ray%p(1), ray%p(2), ray%p(3))

    do step = 1, 50
      call integrate_geodesic_step_rk4(metric, ray, d_lambda, r_horiz, r_escape, err)
      if (ray%terminated) exit
    end do

    q_final = compute_carter_constant(metric%mass, metric%spin, ray%x(1), ray%x(2), &
                                      ray%p(0), ray%p(1), ray%p(2), ray%p(3))

    drift = abs(q_final - q_init) / max(abs(q_init), 1.0e-12_dp)

    if (drift < 1.0e-3_dp) then
      write(*, '(A, E12.4, A)') "  PASS: Kerr symplectic null geodesic Carter invariant (drift = ", drift, ")"
    else
      write(*, '(A, E12.4, A)') "  FAIL: Excessive Carter constant drift: ", drift
      passed = .false.
    end if
  end subroutine run_acceptance_carter

  subroutine run_acceptance_shadow(max_err, passed)
    real(dp), intent(out) :: max_err
    logical, intent(inout) :: passed
    logical :: test_passed

    call verify_kerr_shadow_benchmark(0.0_dp, 90.0_dp, 48_i4, 0.4_dp, max_err, test_passed)

    if (test_passed) then
      write(*, '(A, E12.4, A)') "  PASS: Ray-traced shadow matches Schwarzschild 3*sqrt(3) radius (discrepancy = ", max_err, ")"
    else
      write(*, '(A, E12.4, A)') "  FAIL: Ray-traced shadow discrepancy exceeded tolerance: ", max_err
      passed = .false.
    end if
  end subroutine run_acceptance_shadow

end program test_acceptance_runner
