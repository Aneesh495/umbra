!> umbra: Comprehensive unit test runner executable
!> Executes all subsystem test suites and reports deterministic TAP results
program test_unit_runner
  use umbra_kinds, only: i4
  use test_framework, only: test_suite_t
  use test_math_quadrature, only: run_test_math_quadrature
  use test_spacetime_metric, only: run_test_spacetime_metric
  use test_eos, only: run_test_eos
  use test_con2prim, only: run_test_con2prim
  use test_riemann_solvers, only: run_test_riemann_solvers
  use test_reconstruction, only: run_test_reconstruction
  use test_dg_operators, only: run_test_dg_operators
  use test_geodesics, only: run_test_geodesics
  implicit none

  type(test_suite_t) :: suite_math, suite_spacetime, suite_eos, suite_con2prim
  type(test_suite_t) :: suite_riemann, suite_reconstruct, suite_dg, suite_geodesic
  integer(i4) :: total_tests, total_passed, total_failed
  logical :: all_ok

  write(*, '(A)') "============================================================"
  write(*, '(A)') "umbra: Comprehensive Unit Test Suite Execution"
  write(*, '(A)') "============================================================"

  call suite_math%init("math_quadrature")
  call run_test_math_quadrature(suite_math)
  call suite_math%report()

  call suite_spacetime%init("spacetime_metric")
  call run_test_spacetime_metric(suite_spacetime)
  call suite_spacetime%report()

  call suite_eos%init("equations_of_state")
  call run_test_eos(suite_eos)
  call suite_eos%report()

  call suite_con2prim%init("con2prim_inversion")
  call run_test_con2prim(suite_con2prim)
  call suite_con2prim%report()

  call suite_riemann%init("riemann_solvers")
  call run_test_riemann_solvers(suite_riemann)
  call suite_riemann%report()

  call suite_reconstruct%init("reconstruction_limiters")
  call run_test_reconstruction(suite_reconstruct)
  call suite_reconstruct%report()

  call suite_dg%init("dg_operators")
  call run_test_dg_operators(suite_dg)
  call suite_dg%report()

  call suite_geodesic%init("kerr_geodesics")
  call run_test_geodesics(suite_geodesic)
  call suite_geodesic%report()

  total_tests = suite_math%num_tests + suite_spacetime%num_tests + suite_eos%num_tests + &
                suite_con2prim%num_tests + suite_riemann%num_tests + suite_reconstruct%num_tests + &
                suite_dg%num_tests + suite_geodesic%num_tests

  total_passed = suite_math%num_passed + suite_spacetime%num_passed + suite_eos%num_passed + &
                 suite_con2prim%num_passed + suite_riemann%num_passed + suite_reconstruct%num_passed + &
                 suite_dg%num_passed + suite_geodesic%num_passed

  total_failed = suite_math%num_failed + suite_spacetime%num_failed + suite_eos%num_failed + &
                 suite_con2prim%num_failed + suite_riemann%num_failed + suite_reconstruct%num_failed + &
                 suite_dg%num_failed + suite_geodesic%num_failed

  all_ok = (total_failed == 0_i4)

  write(*, '(A)') "============================================================"
  write(*, '(A, I0, A, I0, A, I0, A)') "OVERALL RESULT: ", total_tests, " tests (", &
                                       total_passed, " passed, ", total_failed, " failed)"
  if (all_ok) then
    write(*, '(A)') ">>> ALL UNIT TESTS PASSED SUCCESSFULLY <<<"
    write(*, '(A)') "============================================================"
  else
    write(*, '(A)') ">>> ERROR: ONE OR MORE TESTS FAILED <<<"
    write(*, '(A)') "============================================================"
    stop 1
  end if

end program test_unit_runner
