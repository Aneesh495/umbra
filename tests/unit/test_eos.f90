!> Unit tests for Equations of State (Ideal, Polytropic, Synge)
module test_eos
  use umbra_kinds, only: dp, i4
  use umbra_eos_ideal, only: ideal_gas_eos_t, make_ideal_gas_eos
  use umbra_eos_polytropic, only: polytropic_eos_t, make_polytropic_eos
  use umbra_eos_synge, only: synge_gas_eos_t, make_synge_gas_eos
  use test_framework, only: test_suite_t, assert_true, assert_close_dp
  implicit none
  private

  public :: run_test_eos

contains

  subroutine run_test_eos(suite)
    type(test_suite_t), intent(inout) :: suite
    type(ideal_gas_eos_t) :: ideal
    type(polytropic_eos_t) :: poly
    type(synge_gas_eos_t) :: synge
    real(dp) :: p, cs_sq, h

    ! Test 1: Ideal gas pressure p = (gamma - 1) * rho * eps
    ideal = make_ideal_gas_eos(1.4_dp)
    p = ideal%pressure(1.0_dp, 2.5_dp)
    call assert_close_dp(suite, p, 1.0_dp, 1.0e-12_dp, "Ideal gas p = (1.4-1)*1*2.5 = 1.0")

    ! Test 2: Ideal gas sound speed squared cs^2 = gamma * p / (rho * h)
    cs_sq = ideal%sound_speed_sq(1.0_dp, 2.5_dp, p)
    h = 1.0_dp + 2.5_dp + p / 1.0_dp ! h = 1 + eps + p/rho = 4.5
    call assert_close_dp(suite, cs_sq, 1.4_dp * 1.0_dp / 4.5_dp, 1.0e-12_dp, "Ideal gas sound speed squared")

    ! Test 3: Polytropic EOS p = K * rho^Gamma (K=100, Gamma=2)
    poly = make_polytropic_eos(2.0_dp, 100.0_dp)
    p = poly%pressure(0.1_dp, 0.0_dp)
    call assert_close_dp(suite, p, 1.0_dp, 1.0e-12_dp, "Polytropic p = 100 * (0.1)^2 = 1.0")

    ! Test 4: Synge relativistic gas specific enthalpy h > 1
    synge = make_synge_gas_eos()
    h = synge%specific_enthalpy(1.0_dp, 1.0_dp, 1.0_dp)
    call assert_true(suite, h > 1.0_dp, "Synge specific enthalpy h > 1")
  end subroutine run_test_eos

end module test_eos
