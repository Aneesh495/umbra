! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_state
! Description: Relativistic magnetohydrodynamic primitive and conserved state vectors.

module umbra_state
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, NCONS, DHO_FLOOR, PRESS_FLOOR, GAMMA_MAX, VEL_MAX, &
                             PRIM_RHO, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_P, &
                             PRIM_BX, PRIM_BY, PRIM_BZ, &
                             CONS_D, CONS_SX, CONS_SY, CONS_SZ, CONS_TAU, &
                             CONS_BX, CONS_BY, CONS_BZ
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_SUPERLUMINAL_VEL, &
                         UMBRA_ERR_NEGATIVE_PRESSURE, UMBRA_ERR_NEGATIVE_DENSITY
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  implicit none
  private

  type, public :: prim_state_t
    real(dp) :: rho = 1.0_dp                ! Rest-mass density rho
    real(dp), dimension(3) :: v = 0.0_dp    ! Spatial 3-velocity v^i (relative to normal observer)
    real(dp) :: p = 1.0_dp                  ! Thermal gas pressure p
    real(dp), dimension(3) :: b = 0.0_dp    ! Spatial magnetic field B^i
  contains
    procedure, public :: to_array => prim_to_array
    procedure, public :: from_array => prim_from_array
    procedure, public :: lorentz_factor => prim_lorentz_factor
    procedure, public :: magnetic_field_4vec => prim_b_4vec
    procedure, public :: magnetic_energy => prim_b_sq
    procedure, public :: total_pressure => prim_total_pressure
  end type prim_state_t

  type, public :: cons_state_t
    real(dp) :: d = 1.0_dp                  ! Relativistic mass density D = gamma * rho
    real(dp), dimension(3) :: s = 0.0_dp    ! Relativistic momentum density S_i
    real(dp) :: tau = 1.0_dp                ! Relativistic energy density minus rest mass tau
    real(dp), dimension(3) :: b = 0.0_dp    ! Spatial magnetic field B^i
  contains
    procedure, public :: to_array => cons_to_array
    procedure, public :: from_array => cons_from_array
  end type cons_state_t

  public :: prim_to_cons

contains

  pure function prim_to_array(this) result(arr)
    class(prim_state_t), intent(in) :: this
    real(dp), dimension(NPRIM) :: arr
    arr(PRIM_RHO) = this%rho
    arr(PRIM_VX)  = this%v(1)
    arr(PRIM_VY)  = this%v(2)
    arr(PRIM_VZ)  = this%v(3)
    arr(PRIM_P)   = this%p
    arr(PRIM_BX)  = this%b(1)
    arr(PRIM_BY)  = this%b(2)
    arr(PRIM_BZ)  = this%b(3)
  end function prim_to_array

  pure subroutine prim_from_array(this, arr)
    class(prim_state_t), intent(inout) :: this
    real(dp), dimension(NPRIM), intent(in) :: arr
    this%rho  = arr(PRIM_RHO)
    this%v(1) = arr(PRIM_VX)
    this%v(2) = arr(PRIM_VY)
    this%v(3) = arr(PRIM_VZ)
    this%p    = arr(PRIM_P)
    this%b(1) = arr(PRIM_BX)
    this%b(2) = arr(PRIM_BY)
    this%b(3) = arr(PRIM_BZ)
  end subroutine prim_from_array

  pure function cons_to_array(this) result(arr)
    class(cons_state_t), intent(in) :: this
    real(dp), dimension(NCONS) :: arr
    arr(CONS_D)   = this%d
    arr(CONS_SX)  = this%s(1)
    arr(CONS_SY)  = this%s(2)
    arr(CONS_SZ)  = this%s(3)
    arr(CONS_TAU) = this%tau
    arr(CONS_BX)  = this%b(1)
    arr(CONS_BY)  = this%b(2)
    arr(CONS_BZ)  = this%b(3)
  end function cons_to_array

  pure subroutine cons_from_array(this, arr)
    class(cons_state_t), intent(inout) :: this
    real(dp), dimension(NCONS), intent(in) :: arr
    this%d    = arr(CONS_D)
    this%s(1) = arr(CONS_SX)
    this%s(2) = arr(CONS_SY)
    this%s(3) = arr(CONS_SZ)
    this%tau  = arr(CONS_TAU)
    this%b(1) = arr(CONS_BX)
    this%b(2) = arr(CONS_BY)
    this%b(3) = arr(CONS_BZ)
  end subroutine cons_from_array

  pure function prim_lorentz_factor(this, m) result(gamma)
    class(prim_state_t), intent(in) :: this
    type(metric_3p1_t), intent(in) :: m
    real(dp) :: gamma
    real(dp) :: v_sq
    integer(i4) :: i, j

    v_sq = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        v_sq = v_sq + m%gamma_down(i, j) * this%v(i) * this%v(j)
      end do
    end do

    if (v_sq >= 1.0_dp) then
      gamma = GAMMA_MAX
    else
      gamma = 1.0_dp / sqrt(max(1.0_dp - v_sq, 1.0e-12_dp))
      gamma = min(gamma, GAMMA_MAX)
    end if
  end function prim_lorentz_factor

  pure subroutine prim_b_4vec(this, m, b_0, b_i)
    class(prim_state_t), intent(in) :: this
    type(metric_3p1_t), intent(in) :: m
    real(dp), intent(out) :: b_0
    real(dp), dimension(3), intent(out) :: b_i

    real(dp) :: gamma, b_dot_v
    integer(i4) :: i, j

    gamma = this%lorentz_factor(m)

    ! B . v = gamma_ij B^i v^j
    b_dot_v = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        b_dot_v = b_dot_v + m%gamma_down(i, j) * this%b(i) * this%v(j)
      end do
    end do

    ! Comoving magnetic 4-vector components:
    ! b^0 = gamma (B . v) / alpha
    ! b^i = B^i / gamma + alpha b^0 (v^i - beta^i / alpha)
    b_0 = (gamma * b_dot_v) / m%alpha
    do i = 1, 3
      b_i(i) = this%b(i) / gamma + m%alpha * b_0 * (this%v(i) - m%beta_up(i) / m%alpha)
    end do
  end subroutine prim_b_4vec

  pure function prim_b_sq(this, m) result(b2)
    class(prim_state_t), intent(in) :: this
    type(metric_3p1_t), intent(in) :: m
    real(dp) :: b2

    real(dp) :: gamma, b_dot_v, b_mag_sq
    integer(i4) :: i, j

    gamma = this%lorentz_factor(m)

    b_mag_sq = 0.0_dp
    b_dot_v = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        b_mag_sq = b_mag_sq + m%gamma_down(i, j) * this%b(i) * this%b(j)
        b_dot_v  = b_dot_v  + m%gamma_down(i, j) * this%b(i) * this%v(j)
      end do
    end do

    ! b^2 = B^2 / gamma^2 + (B . v)^2
    b2 = b_mag_sq / (gamma * gamma) + b_dot_v * b_dot_v
  end function prim_b_sq

  pure function prim_total_pressure(this, m) result(p_tot)
    class(prim_state_t), intent(in) :: this
    type(metric_3p1_t), intent(in) :: m
    real(dp) :: p_tot

    p_tot = this%p + 0.5_dp * this%magnetic_energy(m)
  end function prim_total_pressure

  subroutine prim_to_cons(p, m, eos, u, err)
    type(prim_state_t), intent(in) :: p
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    type(cons_state_t), intent(out) :: u
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: gamma, eps, h, b2, w, b_dot_v, b_mag_sq, v_sq
    real(dp), dimension(3) :: v_down, b_down
    integer(i4) :: i, j

    if (p%rho <= 0.0_dp) then
      if (present(err)) call err%set(UMBRA_ERR_NEGATIVE_DENSITY, "Non-positive density in prim_to_cons", "prim_to_cons")
      return
    end if

    if (p%p <= 0.0_dp) then
      if (present(err)) call err%set(UMBRA_ERR_NEGATIVE_PRESSURE, "Non-positive pressure in prim_to_cons", "prim_to_cons")
      return
    end if

    gamma = p%lorentz_factor(m)
    eps = eos%specific_internal_energy(p%rho, p%p)
    h = eos%specific_enthalpy(p%rho, eps, p%p)
    b2 = p%magnetic_energy(m)

    ! W = gamma^2 * rho * h
    w = gamma * gamma * p%rho * h

    ! Compute down-index vectors: v_i = gamma_ij v^j, B_i = gamma_ij B^j
    v_down = 0.0_dp
    b_down = 0.0_dp
    b_dot_v = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        v_down(i) = v_down(i) + m%gamma_down(i, j) * p%v(j)
        b_down(i) = b_down(i) + m%gamma_down(i, j) * p%b(j)
        b_dot_v = b_dot_v + m%gamma_down(i, j) * p%b(i) * p%v(j)
      end do
    end do

    ! Conserved variables:
    ! D = gamma * rho
    u%d = gamma * p%rho

    b_mag_sq = 0.0_dp
    do j = 1, 3
      do i = 1, 3
        b_mag_sq = b_mag_sq + m%gamma_down(i, j) * p%b(i) * p%b(j)
      end do
    end do

    ! S_i = (W + B^2) v_i - (B . v) B_i
    do i = 1, 3
      u%s(i) = (w + b_mag_sq) * v_down(i) - b_dot_v * b_down(i)
    end do

    ! tau = (W - p) + b^2 * gamma^2 - 0.5 * b^2 - (alpha b^0)^2 - D
    !     = W - p + 0.5 * b^2 * (2 gamma^2 - 1) - (gamma * b_dot_v)^2 - D
    v_sq = max(1.0_dp - 1.0_dp / (gamma * gamma), 0.0_dp)
    u%tau = (w - p%p) + 0.5_dp * b_mag_sq * (1.0_dp + v_sq) - 0.5_dp * (b_dot_v**2) - u%d

    ! B^i remains invariant
    u%b = p%b

    if (present(err)) call err%clear()
  end subroutine prim_to_cons

end module umbra_state
