! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_wavespeeds
! Description: Relativistic magnetosonic characteristic speeds and spectral eigenvalue bounds.

module umbra_wavespeeds
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: DHO_FLOOR, PRESS_FLOOR, VEL_MAX
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: metric_3p1_t
  use umbra_eos_base, only: eos_t
  use umbra_state, only: prim_state_t
  implicit none
  private

  public :: compute_characteristic_wavespeeds
  public :: compute_signal_speed_bounds

contains

  pure subroutine compute_characteristic_wavespeeds(p, m, eos, dir, lambda_minus, lambda_plus)
    ! Computes the fast magnetosonic wave speeds lambda_- and lambda_+ along coordinate direction dir
    type(prim_state_t), intent(in) :: p
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: dir
    real(dp), intent(out) :: lambda_minus, lambda_plus

    real(dp) :: eps, h_val, b2, rho_h, va2, cs2, vms2, vms
    real(dp) :: vn, lambda_minus_local, lambda_plus_local, tmp

    eps = eos%specific_internal_energy(p%rho, p%p)
    h_val = eos%specific_enthalpy(p%rho, eps, p%p)
    b2 = p%magnetic_energy(m)

    rho_h = max(p%rho * h_val, DHO_FLOOR)
    cs2 = eos%sound_speed_sq(p%rho, eps, p%p)

    ! Relativistic Alfven speed squared: v_A^2 = b^2 / (rho h + b^2)
    va2 = b2 / max(rho_h + b2, 1.0e-30_dp)
    va2 = min(max(va2, 0.0_dp), 0.999999_dp)

    ! Fast magnetosonic speed in fluid comoving frame: v_ms^2 = c_s^2 + v_A^2 - c_s^2 * v_A^2
    vms2 = cs2 + va2 - cs2 * va2
    vms2 = min(max(vms2, 1.0e-12_dp), 0.999999_dp)
    vms = sqrt(vms2)

    ! Normal velocity component
    vn = p%v(dir)

    ! Relativistic velocity addition: lambda_local = (v_n +/- v_ms) / (1 +/- v_n * v_ms)
    lambda_minus_local = (vn - vms) / (1.0_dp - vn * vms)
    lambda_plus_local  = (vn + vms) / (1.0_dp + vn * vms)

    ! Transform to coordinate frame with shift and lapse: lambda = alpha * lambda_local - beta^k
    lambda_minus = m%alpha * lambda_minus_local - m%beta_up(dir)
    lambda_plus  = m%alpha * lambda_plus_local  - m%beta_up(dir)

    ! Ensure ordering
    if (lambda_minus > lambda_plus) then
      tmp = lambda_minus
      lambda_minus = lambda_plus
      lambda_plus = tmp
    end if
  end subroutine compute_characteristic_wavespeeds

  subroutine compute_signal_speed_bounds(p_left, p_right, m, eos, dir, sl, sr)
    ! Computes Davis/Einfeldt bounding signal speeds across an interface
    type(prim_state_t), intent(in) :: p_left, p_right
    type(metric_3p1_t), intent(in) :: m
    class(eos_t), intent(in) :: eos
    integer(i4), intent(in) :: dir
    real(dp), intent(out) :: sl, sr

    real(dp) :: sl_left, sr_left, sl_right, sr_right

    call compute_characteristic_wavespeeds(p_left, m, eos, dir, sl_left, sr_left)
    call compute_characteristic_wavespeeds(p_right, m, eos, dir, sl_right, sr_right)

    sl = min(0.0_dp, min(sl_left, sl_right))
    sr = max(0.0_dp, max(sr_left, sr_right))
  end subroutine compute_signal_speed_bounds

end module umbra_wavespeeds
