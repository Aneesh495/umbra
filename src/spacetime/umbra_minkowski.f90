! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_minkowski
! Description: Flat Minkowski spacetime metric implementation with coordinate mapping.

module umbra_minkowski
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NDIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t
  implicit none
  private

  type, extends(spacetime_metric_t), public :: minkowski_metric_t
    logical :: spherical_coords = .false.
  contains
    procedure, public :: eval => minkowski_eval
    procedure, public :: to_cartesian => minkowski_to_cartesian
  end type minkowski_metric_t

contains

  subroutine minkowski_eval(this, x, m, err)
    class(minkowski_metric_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: x
    type(metric_3p1_t), intent(out) :: m
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: r, th, sin_th

    m%alpha = 1.0_dp
    m%beta_up = 0.0_dp
    m%gamma_down = 0.0_dp

    if (.not. this%spherical_coords) then
      ! Cartesian coordinates (t, x, y, z)
      m%gamma_down(1,1) = 1.0_dp
      m%gamma_down(2,2) = 1.0_dp
      m%gamma_down(3,3) = 1.0_dp
      m%det_gamma = 1.0_dp
      m%sqrt_det_gamma = 1.0_dp
      m%sqrt_det_g = 1.0_dp

      m%beta_down = 0.0_dp
      m%gamma_up = m%gamma_down

      m%g_down = 0.0_dp
      m%g_down(0,0) = -1.0_dp
      m%g_down(1,1) = 1.0_dp
      m%g_down(2,2) = 1.0_dp
      m%g_down(3,3) = 1.0_dp

      m%g_up = m%g_down
    else
      ! Spherical polar coordinates (t, r, theta, phi)
      r = max(x(1), 1.0e-12_dp)
      th = x(2)
      sin_th = max(abs(sin(th)), 1.0e-12_dp)

      m%gamma_down(1,1) = 1.0_dp
      m%gamma_down(2,2) = r * r
      m%gamma_down(3,3) = r * r * sin_th * sin_th

      m%det_gamma = r**4 * sin_th**2
      m%sqrt_det_gamma = r * r * sin_th
      m%sqrt_det_g = m%sqrt_det_gamma

      m%beta_down = 0.0_dp
      m%gamma_up = 0.0_dp
      m%gamma_up(1,1) = 1.0_dp
      m%gamma_up(2,2) = 1.0_dp / (r * r)
      m%gamma_up(3,3) = 1.0_dp / (r * r * sin_th * sin_th)

      m%g_down = 0.0_dp
      m%g_down(0,0) = -1.0_dp
      m%g_down(1,1) = m%gamma_down(1,1)
      m%g_down(2,2) = m%gamma_down(2,2)
      m%g_down(3,3) = m%gamma_down(3,3)

      m%g_up = 0.0_dp
      m%g_up(0,0) = -1.0_dp
      m%g_up(1,1) = m%gamma_up(1,1)
      m%g_up(2,2) = m%gamma_up(2,2)
      m%g_up(3,3) = m%gamma_up(3,3)
    end if

    if (present(err)) call err%clear()
  end subroutine minkowski_eval

  pure subroutine minkowski_to_cartesian(this, x, x_cart)
    class(minkowski_metric_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: x
    real(dp), dimension(0:3), intent(out) :: x_cart

    if (.not. this%spherical_coords) then
      x_cart = x
    else
      ! Spherical to Cartesian: x = r sin(th) cos(phi), y = r sin(th) sin(phi), z = r cos(th)
      x_cart(0) = x(0)
      x_cart(1) = x(1) * sin(x(2)) * cos(x(3))
      x_cart(2) = x(1) * sin(x(2)) * sin(x(3))
      x_cart(3) = x(1) * cos(x(2))
    end if
  end subroutine minkowski_to_cartesian

end module umbra_minkowski
