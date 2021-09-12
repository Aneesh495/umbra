! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_grid3d
! Description: Three-dimensional structured computational grid for volumetric GRMHD domains.

module umbra_grid3d
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  type, public :: grid3d_t
    integer(i4) :: nx = 32_i4
    integer(i4) :: ny = 32_i4
    integer(i4) :: nz = 32_i4
    integer(i4) :: nghost = 2_i4
    real(dp) :: x_min = -1.0_dp, x_max = 1.0_dp
    real(dp) :: y_min = -1.0_dp, y_max = 1.0_dp
    real(dp) :: z_min = -1.0_dp, z_max = 1.0_dp
    real(dp) :: dx = 0.0_dp, dy = 0.0_dp, dz = 0.0_dp

    real(dp), allocatable, dimension(:,:,:,:) :: prim
    real(dp), allocatable, dimension(:,:,:,:) :: cons
  contains
    procedure, public :: init => grid3d_init
    procedure, public :: free => grid3d_free
  end type grid3d_t

contains

  subroutine grid3d_init(this, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, nghost, err)
    class(grid3d_t), intent(inout) :: this
    integer(i4), intent(in) :: nx, ny, nz
    real(dp), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer(i4), intent(in), optional :: nghost
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: ng

    ng = 2_i4
    if (present(nghost)) ng = nghost

    this%nx = nx
    this%ny = ny
    this%nz = nz
    this%nghost = ng
    this%x_min = xmin
    this%x_max = xmax
    this%y_min = ymin
    this%y_max = ymax
    this%z_min = zmin
    this%z_max = zmax

    this%dx = (xmax - xmin) / real(nx, dp)
    this%dy = (ymax - ymin) / real(ny, dp)
    this%dz = (zmax - zmin) / real(nz, dp)

    call this%free()

    allocate(this%prim(NPRIM, -ng + 1 : nx + ng, -ng + 1 : ny + ng, -ng + 1 : nz + ng))
    allocate(this%cons(NCONS, -ng + 1 : nx + ng, -ng + 1 : ny + ng, -ng + 1 : nz + ng))
    this%prim = 0.0_dp
    this%cons = 0.0_dp

    if (present(err)) call err%clear()
  end subroutine grid3d_init

  subroutine grid3d_free(this)
    class(grid3d_t), intent(inout) :: this
    if (allocated(this%prim)) deallocate(this%prim)
    if (allocated(this%cons)) deallocate(this%cons)
  end subroutine grid3d_free

end module umbra_grid3d
