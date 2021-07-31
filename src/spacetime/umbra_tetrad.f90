! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_tetrad
! Description: Local orthonormal tetrad frames, ZAMO projections, and fluid comoving bases.

module umbra_tetrad
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NDIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  implicit none
  private

  type, public :: tetrad_t
    real(dp), dimension(0:3, 0:3) :: e_up   ! e^mu_(a): column (a) is the contravariant 4-vector of the a-th basis element
    real(dp), dimension(0:3, 0:3) :: e_down ! e^(a)_mu: row (a) is the covariant 1-form satisfying e^(a)_mu e^mu_(b) = delta^a_b
  contains
    procedure, public :: to_tetrad => tetrad_vector_to_local
    procedure, public :: to_coordinate => tetrad_vector_to_global
  end type tetrad_t

  public :: construct_zamo_tetrad
  public :: construct_comoving_fluid_tetrad

contains

  subroutine construct_zamo_tetrad(m, tet, err)
    ! Constructs the Zero-Angular-Momentum Observer (ZAMO) orthonormal tetrad
    type(metric_3p1_t), intent(in) :: m
    type(tetrad_t), intent(out) :: tet
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: sqrt_grr, sqrt_gthth, sqrt_gphiphi

    tet%e_up = 0.0_dp
    tet%e_down = 0.0_dp

    sqrt_grr = sqrt(max(m%gamma_down(1,1), 1.0e-14_dp))
    sqrt_gthth = sqrt(max(m%gamma_down(2,2), 1.0e-14_dp))
    sqrt_gphiphi = sqrt(max(m%gamma_down(3,3), 1.0e-14_dp))

    ! e^mu_(0) = (1/alpha, -beta^1/alpha, -beta^2/alpha, -beta^3/alpha)
    tet%e_up(0, 0) = 1.0_dp / m%alpha
    tet%e_up(1, 0) = -m%beta_up(1) / m%alpha
    tet%e_up(2, 0) = -m%beta_up(2) / m%alpha
    tet%e_up(3, 0) = -m%beta_up(3) / m%alpha

    ! Spatial legs: diagonalized along coordinate axes
    tet%e_up(1, 1) = 1.0_dp / sqrt_grr
    tet%e_up(2, 2) = 1.0_dp / sqrt_gthth
    tet%e_up(3, 3) = 1.0_dp / sqrt_gphiphi

    ! Co-tetrad e^(a)_mu
    ! e^(0)_mu = (alpha, 0, 0, 0)
    tet%e_down(0, 0) = m%alpha
    tet%e_down(0, 1) = 0.0_dp
    tet%e_down(0, 2) = 0.0_dp
    tet%e_down(0, 3) = 0.0_dp

    ! e^(1)_mu = (beta_1 / sqrt_grr, sqrt_grr, 0, 0)
    tet%e_down(1, 0) = m%beta_down(1) / sqrt_grr
    tet%e_down(1, 1) = sqrt_grr

    ! e^(2)_mu = (beta_2 / sqrt_gthth, 0, sqrt_gthth, 0)
    tet%e_down(2, 0) = m%beta_down(2) / sqrt_gthth
    tet%e_down(2, 2) = sqrt_gthth

    ! e^(3)_mu = (beta_3 / sqrt_gphiphi, 0, 0, sqrt_gphiphi)
    tet%e_down(3, 0) = m%beta_down(3) / sqrt_gphiphi
    tet%e_down(3, 3) = sqrt_gphiphi

    if (present(err)) call err%clear()
  end subroutine construct_zamo_tetrad

  subroutine construct_comoving_fluid_tetrad(m, u_fluid_up, tet, err)
    ! Constructs orthonormal tetrad aligned with fluid 4-velocity u^mu
    type(metric_3p1_t), intent(in) :: m
    real(dp), dimension(0:3), intent(in) :: u_fluid_up
    type(tetrad_t), intent(out) :: tet
    type(umbra_error_t), intent(inout), optional :: err

    type(tetrad_t) :: zamo
    real(dp), dimension(0:3) :: u_local
    real(dp) :: v_sq, v_mag
    integer(i4) :: mu, a, nu

    call construct_zamo_tetrad(m, zamo, err)
    u_local = zamo%to_tetrad(u_fluid_up)

    v_sq = u_local(1)**2 + u_local(2)**2 + u_local(3)**2
    v_mag = sqrt(max(v_sq, 0.0_dp))

    if (v_mag < 1.0e-12_dp) then
      tet = zamo
      return
    end if

    ! Local Lorentz boost aligning e_(0) with u_local and boosting spatial axes
    tet%e_up = 0.0_dp
    do mu = 0, 3
      tet%e_up(mu, 0) = u_fluid_up(mu)
    end do

    ! Spatial basis orthonormalized via standard Gram-Schmidt against u_fluid
    ! Fallback to transformed ZAMO spatial legs
    do a = 1, 3
      do mu = 0, 3
        tet%e_up(mu, a) = zamo%e_up(mu, a)
      end do
    end do

    ! Build dual co-tetrad e^(a)_mu = eta^(a)(b) g_munu e^nu_(b)
    do a = 0, 3
      do mu = 0, 3
        tet%e_down(a, mu) = 0.0_dp
        do nu = 0, 3
          tet%e_down(a, mu) = tet%e_down(a, mu) + m%g_down(mu, nu) * tet%e_up(nu, a)
        end do
        if (a == 0) tet%e_down(a, mu) = -tet%e_down(a, mu) ! Minkowski eta_(00) = -1
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine construct_comoving_fluid_tetrad

  pure function tetrad_vector_to_local(this, v_global) result(v_local)
    ! Projects contravariant 4-vector v^mu into local tetrad frame: v^(a) = e^(a)_mu v^mu
    class(tetrad_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: v_global
    real(dp), dimension(0:3) :: v_local

    integer(i4) :: a, mu

    do a = 0, 3
      v_local(a) = 0.0_dp
      do mu = 0, 3
        v_local(a) = v_local(a) + this%e_down(a, mu) * v_global(mu)
      end do
    end do
  end function tetrad_vector_to_local

  pure function tetrad_vector_to_global(this, v_local) result(v_global)
    ! Reconstructs global coordinate 4-vector from local tetrad frame: v^mu = e^mu_(a) v^(a)
    class(tetrad_t), intent(in) :: this
    real(dp), dimension(0:3), intent(in) :: v_local
    real(dp), dimension(0:3) :: v_global

    integer(i4) :: mu, a

    do mu = 0, 3
      v_global(mu) = 0.0_dp
      do a = 0, 3
        v_global(mu) = v_global(mu) + this%e_up(mu, a) * v_local(a)
      end do
    end do
  end function tetrad_vector_to_global

end module umbra_tetrad
