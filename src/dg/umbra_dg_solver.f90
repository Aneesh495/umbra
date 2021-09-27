! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_dg_solver
! Description: Full Discontinuous Galerkin spectral element solver with Zhang-Shu positivity limiters.

module umbra_dg_solver
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS, NPRIM, DHO_FLOOR, PRESS_FLOOR
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t, cons_state_t, prim_to_cons
  use umbra_con2prim, only: con2prim_solve
  use umbra_fluxes, only: compute_mhd_flux
  use umbra_riemann, only: riemann_solve, RIEMANN_HLLD
  use umbra_dg_basis, only: dg_basis_t
  use umbra_dg_element, only: dg_element_1d_t
  use umbra_dg_volume, only: compute_dg_volume_operator_1d
  use umbra_dg_surface, only: compute_dg_surface_operator_1d, interpolate_to_element_boundaries_1d
  implicit none
  private

  public :: compute_dg_rhs_1d
  public :: apply_positivity_limiter_1d

contains

  subroutine compute_dg_rhs_1d(elements, num_elements, basis, metric, eos, riemann_type, dt_cfl, err)
    ! Evaluates the complete spatial Discontinuous Galerkin RHS: dU/dt = - Volume_Op - Surface_Penalty
    integer(i4), intent(in) :: num_elements
    type(dg_element_1d_t), dimension(num_elements), intent(inout) :: elements
    type(dg_basis_t), intent(in) :: basis
    class(spacetime_metric_t), intent(in) :: metric
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: riemann_type
    real(dp), intent(out), optional :: dt_cfl
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), dimension(NCONS, basis%num_nodes) :: f_nodal, vol_res, surf_res
    real(dp), dimension(NCONS) :: u_l_bnd, u_r_bnd, f_l_int, f_r_int
    real(dp), dimension(NCONS) :: num_flux_l, num_flux_r
    type(prim_state_t) :: p_node, p_l_bnd, p_r_bnd, p_neighbor
    type(cons_state_t) :: u_node, u_bnd
    type(metric_3p1_t) :: m
    real(dp), dimension(0:3) :: x_pt
    real(dp) :: elem_speed, max_speed_global
    integer(i4) :: e, i, n, left_e, right_e

    n = basis%num_nodes
    max_speed_global = 1.0e-6_dp
    x_pt = 0.0_dp

    ! Step 1: Volume contribution within each element
    do e = 1, num_elements
      do i = 1, n
        x_pt(1) = elements(e)%x_nodes(i)
        call metric%eval(x_pt, m, err)

        call u_node%from_array(elements(e)%u_cons(:, i))
        call con2prim_solve(u_node, m, eos, p_node, err)
        elements(e)%p_prim(:, i) = p_node%to_array()

        call compute_mhd_flux(p_node, m, eos, 1_i4, f_nodal(:, i), err)
      end do

      call compute_dg_volume_operator_1d(basis, f_nodal, elements(e)%dx, vol_res, err)
      elements(e)%residual = -vol_res
    end do

    ! Step 2: Surface penalty contributions from element interfaces
    do e = 1, num_elements
      ! Left boundary of element e
      x_pt(1) = elements(e)%x_min
      call metric%eval(x_pt, m, err)

      call interpolate_to_element_boundaries_1d(basis, elements(e)%u_cons, u_l_bnd, u_r_bnd)
      call u_bnd%from_array(u_l_bnd)
      call con2prim_solve(u_bnd, m, eos, p_l_bnd, err)
      call compute_mhd_flux(p_l_bnd, m, eos, 1_i4, f_l_int, err)

      left_e = e - 1
      if (left_e < 1) left_e = num_elements ! Periodic boundary by default
      ! Neighbor state on left side of interface
      call u_bnd%from_array(elements(left_e)%u_cons(:, n))
      call con2prim_solve(u_bnd, m, eos, p_neighbor, err)

      call riemann_solve(riemann_type, p_neighbor, p_l_bnd, m, eos, 1_i4, num_flux_l, elem_speed, err)
      max_speed_global = max(max_speed_global, elem_speed)

      ! Right boundary of element e
      x_pt(1) = elements(e)%x_max
      call metric%eval(x_pt, m, err)

      call u_bnd%from_array(u_r_bnd)
      call con2prim_solve(u_bnd, m, eos, p_r_bnd, err)
      call compute_mhd_flux(p_r_bnd, m, eos, 1_i4, f_r_int, err)

      right_e = e + 1
      if (right_e > num_elements) right_e = 1
      call u_bnd%from_array(elements(right_e)%u_cons(:, 1))
      call con2prim_solve(u_bnd, m, eos, p_neighbor, err)

      call riemann_solve(riemann_type, p_r_bnd, p_neighbor, m, eos, 1_i4, num_flux_r, elem_speed, err)
      max_speed_global = max(max_speed_global, elem_speed)

      ! Surface operator
      call compute_dg_surface_operator_1d(basis, num_flux_l, f_l_int, num_flux_r, f_r_int, &
                                          elements(e)%dx, surf_res, err)

      elements(e)%residual = elements(e)%residual - surf_res
    end do

    if (present(dt_cfl)) then
      ! DG CFL limit: dt <= CFL * dx / ( (2 N + 1) * max_speed )
      dt_cfl = elements(1)%dx / (real(2 * basis%order + 1, dp) * max_speed_global)
    end if

    if (present(err)) call err%clear()
  end subroutine compute_dg_rhs_1d

  subroutine apply_positivity_limiter_1d(elem, basis, metric, eos, err)
    ! Zhang-Shu scaling limiter ensuring rho_i >= rho_floor and p_i >= p_floor without violating conservation
    type(dg_element_1d_t), intent(inout) :: elem
    type(dg_basis_t), intent(in) :: basis
    class(spacetime_metric_t), intent(in) :: metric
    class(eos_t), intent(in) :: eos
    type(umbra_error_t), intent(inout), optional :: err

    real(dp), dimension(NCONS) :: u_avg
    real(dp) :: theta_rho, theta_p, theta, rho_avg, p_avg
    real(dp) :: rho_i, p_i
    type(cons_state_t) :: u_test
    type(prim_state_t) :: p_test
    type(metric_3p1_t) :: m
    real(dp), dimension(0:3) :: x_pt
    integer(i4) :: i, v, n

    n = basis%num_nodes
    x_pt = 0.0_dp
    x_pt(1) = 0.5_dp * (elem%x_min + elem%x_max)
    call metric%eval(x_pt, m, err)

    ! Compute cell average U_avg = sum_i 0.5 * w_i * U_i
    do v = 1, NCONS
      u_avg(v) = elem%cell_average(basis, v)
    end do

    call u_test%from_array(u_avg)
    call con2prim_solve(u_test, m, eos, p_test, err)
    rho_avg = p_test%rho
    p_avg = p_test%p

    theta_rho = 1.0_dp
    theta_p = 1.0_dp

    ! Find minimum scaling factor theta in [0, 1]
    do i = 1, n
      call u_test%from_array(elem%u_cons(:, i))
      call con2prim_solve(u_test, m, eos, p_test, err)
      rho_i = p_test%rho
      p_i = p_test%p

      if (rho_i < DHO_FLOOR) then
        theta_rho = min(theta_rho, abs(rho_avg - DHO_FLOOR) / max(abs(rho_avg - rho_i), 1.0e-30_dp))
      end if

      if (p_i < PRESS_FLOOR) then
        theta_p = min(theta_p, abs(p_avg - PRESS_FLOOR) / max(abs(p_avg - p_i), 1.0e-30_dp))
      end if
    end do

    theta = min(theta_rho, theta_p)
    theta = min(max(theta, 0.0_dp), 1.0_dp)

    ! Rescale nodes towards average: U_i = U_avg + theta * (U_i - U_avg)
    if (theta < 0.999999_dp) then
      do i = 1, n
        do v = 1, NCONS
          elem%u_cons(v, i) = u_avg(v) + theta * (elem%u_cons(v, i) - u_avg(v))
        end do
      end do
    end if

    if (present(err)) call err%clear()
  end subroutine apply_positivity_limiter_1d

end module umbra_dg_solver
