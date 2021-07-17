! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_metric_base
! Description: Abstract base interface for 4D spacetime metrics and 3+1 ADM decomposition.

module umbra_metric_base
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NDIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_math, only: mat3_inv, mat4_inv
  implicit none
  private

  ! 3+1 ADM metric evaluation packet at a spacetime point x^mu = (t, x1, x2, x3)
  type, public :: metric_3p1_t
    real(dp) :: alpha                       ! Lapse function alpha
    real(dp), dimension(3) :: beta_up       ! Shift vector beta^i (contravariant)
    real(dp), dimension(3) :: beta_down     ! Shift vector beta_i = gamma_ij beta^j (covariant)
    real(dp), dimension(3,3) :: gamma_down  ! Spatial 3-metric gamma_ij
    real(dp), dimension(3,3) :: gamma_up    ! Inverse spatial 3-metric gamma^ij
    real(dp) :: det_gamma                   ! Determinant gamma = det(gamma_ij)
    real(dp) :: sqrt_det_gamma              ! sqrt(gamma)
    real(dp) :: sqrt_det_g                  ! 4D volume element alpha * sqrt(gamma)
    real(dp), dimension(0:3, 0:3) :: g_down ! Covariant 4-metric g_munu
    real(dp), dimension(0:3, 0:3) :: g_up   ! Contravariant 4-metric g^munu
  contains
    procedure, public :: compute_4metric => metric_3p1_compute_4metric
  end type metric_3p1_t

  ! Abstract metric class
  type, abstract, public :: spacetime_metric_t
    real(dp) :: mass = 1.0_dp
    real(dp) :: spin = 0.0_dp
  contains
    procedure(metric_eval_iface), deferred, public :: eval
    procedure(metric_coord_iface), deferred, public :: to_cartesian
  end type spacetime_metric_t

  abstract interface
    subroutine metric_eval_iface(this, x, m, err)
      import :: spacetime_metric_t, metric_3p1_t, umbra_error_t, dp
      class(spacetime_metric_t), intent(in) :: this
      real(dp), dimension(0:3), intent(in) :: x
      type(metric_3p1_t), intent(out) :: m
      type(umbra_error_t), intent(inout), optional :: err
    end subroutine metric_eval_iface

    pure subroutine metric_coord_iface(this, x, x_cart)
      import :: spacetime_metric_t, dp
      class(spacetime_metric_t), intent(in) :: this
      real(dp), dimension(0:3), intent(in) :: x
      real(dp), dimension(0:3), intent(out) :: x_cart
    end subroutine metric_coord_iface
  end interface

  public :: compute_adm_4metric

contains

  subroutine metric_3p1_compute_4metric(this, err)
    class(metric_3p1_t), intent(inout) :: this
    type(umbra_error_t), intent(inout), optional :: err

    call compute_adm_4metric(this%alpha, this%beta_up, this%gamma_down, &
                             this%beta_down, this%gamma_up, this%det_gamma, &
                             this%sqrt_det_gamma, this%sqrt_det_g, &
                             this%g_down, this%g_up, err)
  end subroutine metric_3p1_compute_4metric

  subroutine compute_adm_4metric(alpha, beta_up, gamma_down, beta_down, gamma_up, &
                                det_gamma, sqrt_det_gamma, sqrt_det_g, g_down, g_up, err)
    real(dp), intent(in) :: alpha
    real(dp), dimension(3), intent(in) :: beta_up
    real(dp), dimension(3,3), intent(in) :: gamma_down
    real(dp), dimension(3), intent(out) :: beta_down
    real(dp), dimension(3,3), intent(out) :: gamma_up
    real(dp), intent(out) :: det_gamma
    real(dp), intent(out) :: sqrt_det_gamma
    real(dp), intent(out) :: sqrt_det_g
    real(dp), dimension(0:3, 0:3), intent(out) :: g_down
    real(dp), dimension(0:3, 0:3), intent(out) :: g_up
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j
    real(dp) :: beta_sq

    ! Invert 3-metric gamma_ij to get gamma^ij and det(gamma)
    call mat3_inv(gamma_down, gamma_up, det_gamma, err)
    if (det_gamma <= 0.0_dp) then
      sqrt_det_gamma = 1.0_dp
      sqrt_det_g = alpha
      return
    end if

    sqrt_det_gamma = sqrt(det_gamma)
    sqrt_det_g = alpha * sqrt_det_gamma

    ! Compute covariant shift beta_i = gamma_ij beta^j
    beta_down = 0.0_dp
    beta_sq = 0.0_dp
    do i = 1, 3
      do j = 1, 3
        beta_down(i) = beta_down(i) + gamma_down(i, j) * beta_up(j)
      end do
      beta_sq = beta_sq + beta_down(i) * beta_up(i)
    end do

    ! Construct covariant 4-metric g_munu:
    ! g_00 = -alpha^2 + beta_i beta^i
    ! g_0i = beta_i
    ! g_ij = gamma_ij
    g_down(0, 0) = -alpha * alpha + beta_sq
    do i = 1, 3
      g_down(0, i) = beta_down(i)
      g_down(i, 0) = beta_down(i)
      do j = 1, 3
        g_down(i, j) = gamma_down(i, j)
      end do
    end do

    ! Construct contravariant 4-metric g^munu:
    ! g^00 = -1 / alpha^2
    ! g^0i = beta^i / alpha^2
    ! g^ij = gamma^ij - (beta^i beta^j) / alpha^2
    g_up(0, 0) = -1.0_dp / (alpha * alpha)
    do i = 1, 3
      g_up(0, i) = beta_up(i) / (alpha * alpha)
      g_up(i, 0) = g_up(0, i)
      do j = 1, 3
        g_up(i, j) = gamma_up(i, j) - (beta_up(i) * beta_up(j)) / (alpha * alpha)
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine compute_adm_4metric

end module umbra_metric_base
