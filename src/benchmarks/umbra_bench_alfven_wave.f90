! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_bench_alfven_wave
! Description: Smooth circularly polarized relativistic Alfven wave for convergence testing (EOC).

module umbra_bench_alfven_wave
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_P, &
                             PRIM_BX, PRIM_BY, PRIM_BZ, TWO_PI
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_grid1d, only: grid1d_t
  implicit none
  private

  public :: init_bench_alfven_wave_1d
  public :: compute_alfven_exact_solution_1d
  public :: compute_eoc_error_norms_1d

contains

  subroutine init_bench_alfven_wave_1d(grid, b0, eta, va, err)
    ! Circularly polarized Alfven wave on periodic domain [0, 1]:
    type(grid1d_t), intent(inout) :: grid
    real(dp), intent(in), optional :: b0, eta
    real(dp), intent(out), optional :: va
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: b_base, amp, v_alfven

    b_base = 1.0_dp
    amp = 0.2_dp
    if (present(b0)) b_base = b0
    if (present(eta)) amp = eta

    ! Relativistic Alfven speed for rho = 1, p = 1, gamma_eos = 4/3: h = 1 + 4 = 5
    ! v_A = B0 / sqrt(rho h + B0^2) = 1.0 / sqrt(5.0 + 1.0) = 1.0 / sqrt(6.0)
    v_alfven = b_base / sqrt(5.0_dp + b_base * b_base)
    if (present(va)) va = v_alfven

    call compute_alfven_exact_solution_1d(grid%x_cell, grid%nx, grid%nghost, 0.0_dp, b_base, amp, &
                                         v_alfven, grid%prim, err)
  end subroutine init_bench_alfven_wave_1d

  subroutine compute_alfven_exact_solution_1d(x_cell, nx, nghost, t, b0, eta, va, prim, err)
    integer(i4), intent(in) :: nx, nghost
    real(dp), dimension(1-nghost : nx+nghost), intent(in) :: x_cell
    real(dp), intent(in) :: t, b0, eta, va
    real(dp), dimension(NPRIM, 1-nghost : nx+nghost), intent(out) :: prim
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: phase_arg, by_val, bz_val, vy_val, vz_val
    integer(i4) :: i

    do i = 1 - nghost, nx + nghost
      phase_arg = TWO_PI * (x_cell(i) - va * t)

      by_val = eta * b0 * cos(phase_arg)
      bz_val = eta * b0 * sin(phase_arg)

      ! Velocity perturbation in phase with transverse magnetic field
      vy_val = -va * by_val / b0
      vz_val = -va * bz_val / b0

      prim(:, i) = 0.0_dp
      prim(PRIM_RHO, i) = 1.0_dp
      prim(PRIM_P,   i) = 1.0_dp
      prim(PRIM_VX,  i) = 0.0_dp
      prim(PRIM_VY,  i) = vy_val
      prim(PRIM_VZ,  i) = vz_val
      prim(PRIM_BX,  i) = b0
      prim(PRIM_BY,  i) = by_val
      prim(PRIM_BZ,  i) = bz_val
    end do

    if (present(err)) call err%clear()
  end subroutine compute_alfven_exact_solution_1d

  subroutine compute_eoc_error_norms_1d(grid, t_final, b0, eta, va, var, l1, l2, linf, err)
    ! Computes L1, L2, and Linf error norms against exact Alfven wave solution
    type(grid1d_t), intent(in) :: grid
    real(dp), intent(in) :: t_final, b0, eta, va
    integer(i4), intent(in) :: var
    real(dp), intent(out) :: l1, l2, linf
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), allocatable, dimension(:,:) :: exact_prim
    real(dp) :: diff, sum_l1, sum_l2, max_linf
    integer(i4) :: i

    allocate(exact_prim(NPRIM, 1 - grid%nghost : grid%nx + grid%nghost))
    call compute_alfven_exact_solution_1d(grid%x_cell, grid%nx, grid%nghost, t_final, b0, eta, va, &
                                         exact_prim, err)

    sum_l1 = 0.0_dp
    sum_l2 = 0.0_dp
    max_linf = 0.0_dp

    do i = 1, grid%nx
      diff = abs(grid%prim(var, i) - exact_prim(var, i))
      sum_l1 = sum_l1 + diff
      sum_l2 = sum_l2 + diff * diff
      max_linf = max(max_linf, diff)
    end do

    l1 = sum_l1 / real(grid%nx, dp)
    l2 = sqrt(sum_l2 / real(grid%nx, dp))
    linf = max_linf

    deallocate(exact_prim)
    if (present(err)) call err%clear()
  end subroutine compute_eoc_error_norms_1d

end module umbra_bench_alfven_wave
