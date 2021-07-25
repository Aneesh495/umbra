! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_schwarzschild
! Description: Schwarzschild spacetime metric in standard and isotropic coordinates.

module umbra_schwarzschild
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NDIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_HORIZON_PENETRATED, &
                         UMBRA_ERR_SINGULAR_METRIC
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t, compute_adm_4metric
  implicit none
  private

  type, extends(spacetime_metric_t), public :: schwarzschild_metric_t
    logical :: isotropic_coords = .false.
  contains
    procedure, public :: eval => schwarzschild_eval
    procedure, public :: to_cartesian => schwarzschild_to_cartesian
    procedure, public :: horizon_radius => schwarzschild_horizon_radius
  end type schwarzschild_metric_t

contains

  pure function schwarzschild_horizon_radius(this) result(rh)
    class(schwarzschild_metric_t), intent(in) :: this
    real(dp) :: rh
    if (.not. this%isotropic_coords) then
      rh = 2.0_dp * this%mass
    else
      rh = 0.5_dp * this%mass
    end if
  end function schwarzschild_horizon_radius

  subroutine schwarzschild_eval(this, x, m, err)
    class(schwarzschild_metric_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: x
    type(metric_3p1_t), intent(out) :: m
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: r, th, sin_th, f_r, psi, rh

    rh = this%horizon_radius()

    if (.not. this%isotropic_coords) then
      ! Standard Schwarzschild coordinates (t, r, theta, phi)
      r = x(1)
      th = x(2)
      sin_th = max(abs(sin(th)), 1.0e-12_dp)

      if (r <= rh + 1.0e-6_dp) then
        if (present(err)) call err%set(UMBRA_ERR_HORIZON_PENETRATED, &
          "Point inside or at Schwarzschild horizon", "schwarzschild_eval")
      end if

      f_r = 1.0_dp - (2.0_dp * this%mass) / max(r, 1.0e-12_dp)
      if (f_r <= 1.0e-14_dp) f_r = 1.0e-14_dp

      m%alpha = sqrt(f_r)
      m%beta_up = 0.0_dp
      m%gamma_down = 0.0_dp

      m%gamma_down(1,1) = 1.0_dp / f_r
      m%gamma_down(2,2) = r * r
      m%gamma_down(3,3) = r * r * sin_th * sin_th

      m%det_gamma = (r**4 * sin_th**2) / f_r
      m%sqrt_det_gamma = (r * r * sin_th) / sqrt(f_r)
      m%sqrt_det_g = m%alpha * m%sqrt_det_gamma

      m%beta_down = 0.0_dp
      m%gamma_up = 0.0_dp
      m%gamma_up(1,1) = f_r
      m%gamma_up(2,2) = 1.0_dp / (r * r)
      m%gamma_up(3,3) = 1.0_dp / (r * r * sin_th * sin_th)

      m%g_down = 0.0_dp
      m%g_down(0,0) = -f_r
      m%g_down(1,1) = m%gamma_down(1,1)
      m%g_down(2,2) = m%gamma_down(2,2)
      m%g_down(3,3) = m%gamma_down(3,3)

      m%g_up = 0.0_dp
      m%g_up(0,0) = -1.0_dp / f_r
      m%g_up(1,1) = m%gamma_up(1,1)
      m%g_up(2,2) = m%gamma_up(2,2)
      m%g_up(3,3) = m%gamma_up(3,3)

    else
      ! Isotropic coordinates (t, r_bar, theta, phi)
      r = max(x(1), 1.0e-12_dp)
      th = x(2)
      sin_th = max(abs(sin(th)), 1.0e-12_dp)

      psi = 1.0_dp + (0.5_dp * this%mass) / r
      m%alpha = abs(1.0_dp - (0.5_dp * this%mass) / r) / psi
      if (m%alpha < 1.0e-14_dp) m%alpha = 1.0e-14_dp

      m%beta_up = 0.0_dp
      m%gamma_down = 0.0_dp
      m%gamma_down(1,1) = psi**4
      m%gamma_down(2,2) = psi**4 * r * r
      m%gamma_down(3,3) = psi**4 * r * r * sin_th * sin_th

      m%det_gamma = psi**12 * r**4 * sin_th**2
      m%sqrt_det_gamma = psi**6 * r * r * sin_th
      m%sqrt_det_g = m%alpha * m%sqrt_det_gamma

      m%beta_down = 0.0_dp
      m%gamma_up = 0.0_dp
      m%gamma_up(1,1) = 1.0_dp / psi**4
      m%gamma_up(2,2) = 1.0_dp / (psi**4 * r * r)
      m%gamma_up(3,3) = 1.0_dp / (psi**4 * r * r * sin_th * sin_th)

      m%g_down = 0.0_dp
      m%g_down(0,0) = -m%alpha**2
      m%g_down(1,1) = m%gamma_down(1,1)
      m%g_down(2,2) = m%gamma_down(2,2)
      m%g_down(3,3) = m%gamma_down(3,3)

      m%g_up = 0.0_dp
      m%g_up(0,0) = -1.0_dp / (m%alpha**2)
      m%g_up(1,1) = m%gamma_up(1,1)
      m%g_up(2,2) = m%gamma_up(2,2)
      m%g_up(3,3) = m%gamma_up(3,3)
    end if

    if (present(err)) call err%clear()
  end subroutine schwarzschild_eval

  pure subroutine schwarzschild_to_cartesian(this, x, x_cart)
    class(schwarzschild_metric_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: x
    real(dp), dimension(0:3), intent(out) :: x_cart

    real(dp) :: r, r_iso

    if (.not. this%isotropic_coords) then
      r = x(1)
    else
      r_iso = x(1)
      r = r_iso * (1.0_dp + (0.5_dp * this%mass) / max(r_iso, 1.0e-12_dp))**2
    end if

    x_cart(0) = x(0)
    x_cart(1) = r * sin(x(2)) * cos(x(3))
    x_cart(2) = r * sin(x(2)) * sin(x(3))
    x_cart(3) = r * cos(x(2))
  end subroutine schwarzschild_to_cartesian

end module umbra_schwarzschild
