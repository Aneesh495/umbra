! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_constants
! Description: Fundamental physical, mathematical, and numerical threshold constants.

module umbra_constants
  use umbra_kinds, only: dp, i4
  implicit none
  private

  ! Mathematical constants
  real(dp), parameter, public :: PI = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter, public :: TWO_PI = 2.0_dp * PI
  real(dp), parameter, public :: HALF_PI = 0.5_dp * PI
  real(dp), parameter, public :: FOUR_PI = 4.0_dp * PI
  real(dp), parameter, public :: INV_PI = 1.0_dp / PI
  real(dp), parameter, public :: INV_FOUR_PI = 1.0_dp / FOUR_PI
  real(dp), parameter, public :: SQRT_TWO = 1.414213562373095048801688724209698078570_dp
  real(dp), parameter, public :: SQRT_THREE = 1.732050807568877293527446341505872366943_dp
  real(dp), parameter, public :: EULER_E = 2.718281828459045235360287471352662497757_dp
  real(dp), parameter, public :: LN_TWO = 0.693147180559945309417232121458176568075_dp

  ! Machine epsilon and numerical tolerances
  real(dp), parameter, public :: EPSILON_DP = epsilon(1.0_dp)
  real(dp), parameter, public :: TINY_DP = tiny(1.0_dp)
  real(dp), parameter, public :: HUGE_DP = huge(1.0_dp)
  real(dp), parameter, public :: TOL_CONVERGENCE = 1.0e-12_dp
  real(dp), parameter, public :: TOL_LOOSE = 1.0e-8_dp

  ! Physical constants in SI units
  real(dp), parameter, public :: SI_C = 299792458.0_dp                     ! Speed of light [m/s]
  real(dp), parameter, public :: SI_G = 6.67430e-11_dp                     ! Gravitational constant [m^3 kg^-1 s^-2]
  real(dp), parameter, public :: SI_M_SUN = 1.98847e30_dp                  ! Solar mass [kg]
  real(dp), parameter, public :: SI_HBAR = 1.054571817e-34_dp              ! Reduced Planck constant [J s]
  real(dp), parameter, public :: SI_KB = 1.380649e-23_dp                   ! Boltzmann constant [J/K]
  real(dp), parameter, public :: SI_E = 1.602176634e-19_dp                 ! Elementary charge [C]
  real(dp), parameter, public :: SI_ME = 9.1093837015e-31_dp               ! Electron mass [kg]
  real(dp), parameter, public :: SI_MP = 1.67262192369e-27_dp              ! Proton mass [kg]
  real(dp), parameter, public :: SI_THOMSON = 6.6524587321e-29_dp          ! Thomson cross-section [m^2]

  ! Geometric unit system (G = c = M = 1)
  real(dp), parameter, public :: GEOM_C = 1.0_dp
  real(dp), parameter, public :: GEOM_G = 1.0_dp
  real(dp), parameter, public :: GEOM_M = 1.0_dp

  ! Physical limits and atmosphere thresholds for relativistic hydrodynamics
  real(dp), parameter, public :: DHO_FLOOR = 1.0e-10_dp                    ! Rest-mass density floor
  real(dp), parameter, public :: PRESS_FLOOR = 1.0e-10_dp                  ! Thermal gas pressure floor
  real(dp), parameter, public :: GAMMA_MAX = 50.0_dp                       ! Maximum Lorentz factor
  real(dp), parameter, public :: VEL_MAX = 0.99979997999_dp                ! Maximum 3-velocity magnitude
  real(dp), parameter, public :: BSQ_OVER_RHO_MAX = 100.0_dp               ! Maximum magnetization limit
  real(dp), parameter, public :: TEMPERATURE_FLOOR = 1.0e-6_dp             ! Minimum dimensionless temperature

  ! Divergence cleaning parameters (Dedner hyperbolic GLM-MHD)
  real(dp), parameter, public :: GLM_ALPHA_P = 0.18_dp                     ! Damping coefficient ratio
  real(dp), parameter, public :: GLM_CR = 0.18_dp                          ! Damping rate scale

  ! Spatial dimensions and tensor rank constants
  integer(i4), parameter, public :: NDIM = 3_i4
  integer(i4), parameter, public :: NPRIM = 8_i4                           ! rho, vx, vy, vz, p, Bx, By, Bz
  integer(i4), parameter, public :: NCONS = 8_i4                           ! D, Sx, Sy, Sz, tau, Bx, By, Bz
  integer(i4), parameter, public :: NSTOKES = 4_i4                         ! I, Q, U, V polarimetric Stokes vector

  ! Primitive variable indices
  integer(i4), parameter, public :: PRIM_RHO = 1_i4
  integer(i4), parameter, public :: PRIM_VX  = 2_i4
  integer(i4), parameter, public :: PRIM_VY  = 3_i4
  integer(i4), parameter, public :: PRIM_VZ  = 4_i4
  integer(i4), parameter, public :: PRIM_P   = 5_i4
  integer(i4), parameter, public :: PRIM_BX  = 6_i4
  integer(i4), parameter, public :: PRIM_BY  = 7_i4
  integer(i4), parameter, public :: PRIM_BZ  = 8_i4

  ! Conserved variable indices
  integer(i4), parameter, public :: CONS_D   = 1_i4
  integer(i4), parameter, public :: CONS_SX  = 2_i4
  integer(i4), parameter, public :: CONS_SY  = 3_i4
  integer(i4), parameter, public :: CONS_SZ  = 4_i4
  integer(i4), parameter, public :: CONS_TAU = 5_i4
  integer(i4), parameter, public :: CONS_BX  = 6_i4
  integer(i4), parameter, public :: CONS_BY  = 7_i4
  integer(i4), parameter, public :: CONS_BZ  = 8_i4

  ! Stokes vector indices
  integer(i4), parameter, public :: STOKES_I = 1_i4
  integer(i4), parameter, public :: STOKES_Q = 2_i4
  integer(i4), parameter, public :: STOKES_U = 3_i4
  integer(i4), parameter, public :: STOKES_V = 4_i4

end module umbra_constants
