! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_bench_fishbone_moncrief
! Description: Stationary relativistic accretion torus equilibrium in Kerr geometry (Fishbone & Moncrief 1976).

module umbra_bench_fishbone_moncrief
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_P, &
                             PRIM_BX, PRIM_BY, PRIM_BZ, DHO_FLOOR, PRESS_FLOOR, PI
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_kerr_ks, only: kerr_ks_metric_t
  use umbra_eos_ideal, only: ideal_gas_eos_t
  use umbra_grid2d, only: grid2d_t
  implicit none
  private

  public :: init_bench_fishbone_moncrief_2d

contains

  subroutine init_bench_fishbone_moncrief_2d(grid, metric, eos, r_in, r_max_density, err)
    ! Initializes Fishbone-Moncrief relativistic thick torus around spinning Kerr black hole
    type(grid2d_t), intent(inout) :: grid
    type(kerr_ks_metric_t), intent(in) :: metric
    type(ideal_gas_eos_t), intent(in) :: eos
    real(dp), intent(in), optional :: r_in, r_max_density
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: rin, rmax, m_val, a_val, l_spec
    real(dp) :: r, th, sin_th, cos_th, rho_max
    real(dp) :: ut_sq_inv, ut_in_sq_inv, h_val, rho_val, p_val, u_phi
    real(dp), dimension(0:3) :: x_pt
    type(metric_3p1_t) :: m
    integer(i4) :: i, j

    m_val = metric%mass
    a_val = metric%spin

    rin = 6.0_dp * m_val
    rmax = 12.0_dp * m_val
    if (present(r_in)) rin = r_in
    if (present(r_max_density)) rmax = r_max_density

    ! Specific angular momentum l = -u_phi / u_t evaluated at pressure maximum r_max
    ! For geodesic circular orbit at r_max:
    l_spec = (rmax**1.5_dp + a_val * sqrt(m_val)) * &
             (rmax**2 - 2.0_dp * a_val * sqrt(m_val * rmax) + a_val**2) / &
             (rmax**1.5_dp * sqrt(rmax - 2.0_dp * m_val) + a_val * sqrt(m_val * rmax))

    ! Potential at inner cusp r_in: ut_in_sq_inv
    x_pt = 0.0_dp
    x_pt(1) = rin
    x_pt(2) = 0.5_dp * PI
    call metric%eval(x_pt, m, err)

    ut_in_sq_inv = -(m%g_up(0,0) - 2.0_dp * l_spec * m%g_up(0,3) + (l_spec**2) * m%g_up(3,3))

    rho_max = 1.0_dp

    do j = 1 - grid%nghost, grid%ny + grid%nghost
      th = grid%y_cell(j)
      sin_th = max(abs(sin(th)), 1.0e-6_dp)
      cos_th = cos(th)

      do i = 1 - grid%nghost, grid%nx + grid%nghost
        r = grid%x_cell(i)
        x_pt(1) = r
        x_pt(2) = th

        call metric%eval(x_pt, m, err)

        grid%prim(:, i, j) = 0.0_dp

        ! Evaluate relativistic effective potential W(r, theta) = ln(h)
        ut_sq_inv = -(m%g_up(0,0) - 2.0_dp * l_spec * m%g_up(0,3) + (l_spec**2) * m%g_up(3,3))

        if (ut_sq_inv > 0.0_dp .and. ut_in_sq_inv > 0.0_dp .and. r >= rin) then
          h_val = sqrt(ut_in_sq_inv / ut_sq_inv)
        else
          h_val = 1.0_dp
        end if

        if (h_val > 1.000001_dp) then
          ! Inside torus
          p_val = ((eos%gamma - 1.0_dp) / eos%gamma) * (h_val - 1.0_dp)
          rho_val = (p_val / 0.01_dp)**(1.0_dp / eos%gamma)

          ! Angular velocity Omega = v^phi = - (g_{t phi} + l g_{tt}) / (g_{phi phi} + l g_{t phi})
          u_phi = (m%g_down(0,3) + l_spec * m%g_down(0,0)) / &
                  max(m%g_down(3,3) + l_spec * m%g_down(0,3), 1.0e-30_dp)

          grid%prim(PRIM_RHO, i, j) = max(rho_val, DHO_FLOOR)
          grid%prim(PRIM_P,   i, j) = max(p_val, PRESS_FLOOR)
          grid%prim(PRIM_VX,  i, j) = 0.0_dp
          grid%prim(PRIM_VY,  i, j) = 0.0_dp
          grid%prim(PRIM_VZ,  i, j) = -u_phi ! Azimuthal rotation v^phi
        else
          ! Ambient magnetospheric atmosphere
          grid%prim(PRIM_RHO, i, j) = DHO_FLOOR * (r / rin)**(-1.5_dp)
          grid%prim(PRIM_P,   i, j) = PRESS_FLOOR * (r / rin)**(-2.5_dp)
          grid%prim(PRIM_VX:PRIM_VZ, i, j) = 0.0_dp
        end if

        ! Seed magnetic loop in torus: A_phi = max(rho - 0.2 * rho_max, 0)
        if (grid%prim(PRIM_RHO, i, j) > 0.2_dp * rho_max) then
          grid%prim(PRIM_BX, i, j) = 0.01_dp * (grid%prim(PRIM_RHO, i, j) / rho_max) * cos_th
          grid%prim(PRIM_BY, i, j) = -0.01_dp * (grid%prim(PRIM_RHO, i, j) / rho_max) * sin_th
        else
          grid%prim(PRIM_BX:PRIM_BZ, i, j) = 0.0_dp
        end if
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine init_bench_fishbone_moncrief_2d

end module umbra_bench_fishbone_moncrief
