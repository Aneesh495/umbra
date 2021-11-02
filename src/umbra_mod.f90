!> umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
!> Main facade module re-exporting the complete public API of the engine
module umbra
  ! Core foundations
  use umbra_kinds
  use umbra_constants
  use umbra_error
  use umbra_logging
  use umbra_quadrature
  use umbra_polynomials
  use umbra_math
  use umbra_timer

  ! Spacetime geometry
  use umbra_metric_base
  use umbra_minkowski
  use umbra_schwarzschild
  use umbra_kerr_bl
  use umbra_kerr_ks
  use umbra_connection
  use umbra_tetrad

  ! Equations of state
  use umbra_eos_base
  use umbra_eos_ideal
  use umbra_eos_polytropic
  use umbra_eos_synge

  ! Relativistic MHD
  use umbra_state
  use umbra_wavespeeds
  use umbra_fluxes
  use umbra_sources
  use umbra_con2prim
  use umbra_con2prim_1d
  use umbra_con2prim_2d
  use umbra_riemann
  use umbra_riemann_llf
  use umbra_riemann_hll
  use umbra_riemann_hllc
  use umbra_riemann_hlld

  ! Divergence cleaning
  use umbra_ct
  use umbra_glm

  ! Mesh hierarchy
  use umbra_boundary
  use umbra_grid1d
  use umbra_grid2d
  use umbra_grid3d
  use umbra_octree

  ! High-order reconstruction
  use umbra_limiters
  use umbra_reconstruct_tvd
  use umbra_reconstruct_ppm
  use umbra_reconstruct_weno

  ! Discontinuous Galerkin formulation
  use umbra_dg_basis
  use umbra_dg_volume
  use umbra_dg_surface
  use umbra_dg_element
  use umbra_dg_solver

  ! Time integration
  use umbra_rk_ssprk
  use umbra_rk_dopri5
  use umbra_cfl

  ! Kerr null geodesics & radiative transfer
  use umbra_geodesic_integrator
  use umbra_carter_constant
  use umbra_camera
  use umbra_synchrotron
  use umbra_polarimetry
  use umbra_radiative_transfer

  ! Standard benchmarks
  use umbra_bench_sod
  use umbra_bench_briowu
  use umbra_bench_orszag_tang
  use umbra_bench_kelvin_helmholtz
  use umbra_bench_blast
  use umbra_bench_alfven_wave
  use umbra_bench_fishbone_moncrief
  use umbra_bench_kerr_shadow

  ! Input / Output
  use umbra_io_raw
  use umbra_io_vtk
  use umbra_io_ascii
  use umbra_checkpoint

  ! C ABI interop
  use umbra_c

  implicit none
  public

end module umbra
