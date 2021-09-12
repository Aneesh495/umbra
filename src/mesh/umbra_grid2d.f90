! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_grid2d
! Description: Two-dimensional curvilinear structured computational grid with metric caching.

module umbra_grid2d
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_metric_base, only: spacetime_metric_t, metric_3p1_t
  use umbra_boundary, only: apply_boundary_conditions_2d, BC_PERIODIC, BC_OUTFLOW
  implicit none
  private

  type, public :: grid2d_t
    integer(i4) :: nx = 64_i4
    integer(i4) :: ny = 64_i4
    integer(i4) :: nghost = 3_i4
    real(dp) :: x_min = -1.0_dp
    real(dp) :: x_max =  1.0_dp
    real(dp) :: y_min = -1.0_dp
    real(dp) :: y_max =  1.0_dp
    real(dp) :: dx = 0.0_dp
    real(dp) :: dy = 0.0_dp

    ! Coordinates
    real(dp), allocatable, dimension(:) :: x_cell
    real(dp), allocatable, dimension(:) :: y_cell
    real(dp), allocatable, dimension(:) :: x_face
    real(dp), allocatable, dimension(:) :: y_face

    ! State fields: primitive P(NPRIM, i, j) and conserved U(NCONS, i, j)
    real(dp), allocatable, dimension(:,:,:) :: prim
    real(dp), allocatable, dimension(:,:,:) :: cons

    ! Staggered face magnetic fields for Constrained Transport:
    ! bx_face(nx+1, ny), by_face(nx, ny+1)
    real(dp), allocatable, dimension(:,:) :: bx_face
    real(dp), allocatable, dimension(:,:) :: by_face

    ! Divergence cleaning scalar potential
    real(dp), allocatable, dimension(:,:) :: psi

    ! Boundary condition types
    integer(i4) :: bc_x_min = BC_PERIODIC
    integer(i4) :: bc_x_max = BC_PERIODIC
    integer(i4) :: bc_y_min = BC_PERIODIC
    integer(i4) :: bc_y_max = BC_PERIODIC
  contains
    procedure, public :: init => grid2d_init
    procedure, public :: free => grid2d_free
    procedure, public :: apply_boundaries => grid2d_apply_boundaries
  end type grid2d_t

contains

  subroutine grid2d_init(this, nx, ny, xmin, xmax, ymin, ymax, nghost, err)
    class(grid2d_t), intent(inout) :: this
    integer(i4), intent(in) :: nx, ny
    real(dp), intent(in) :: xmin, xmax, ymin, ymax
    integer(i4), intent(in), optional :: nghost
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: ng, i, j

    ng = 3_i4
    if (present(nghost)) ng = nghost

    this%nx = nx
    this%ny = ny
    this%nghost = ng
    this%x_min = xmin
    this%x_max = xmax
    this%y_min = ymin
    this%y_max = ymax

    this%dx = (xmax - xmin) / real(nx, dp)
    this%dy = (ymax - ymin) / real(ny, dp)

    call this%free()

    allocate(this%x_cell(-ng + 1 : nx + ng))
    allocate(this%y_cell(-ng + 1 : ny + ng))
    allocate(this%x_face(-ng + 1 : nx + ng + 1))
    allocate(this%y_face(-ng + 1 : ny + ng + 1))

    allocate(this%prim(NPRIM, -ng + 1 : nx + ng, -ng + 1 : ny + ng))
    allocate(this%cons(NCONS, -ng + 1 : nx + ng, -ng + 1 : ny + ng))

    allocate(this%bx_face(1 : nx + 1, 1 : ny))
    allocate(this%by_face(1 : nx, 1 : ny + 1))
    allocate(this%psi(1 : nx, 1 : ny))

    ! Coordinate generation
    do i = -ng + 1, nx + ng
      this%x_cell(i) = this%x_min + (real(i, dp) - 0.5_dp) * this%dx
      this%x_face(i) = this%x_min + real(i - 1, dp) * this%dx
    end do
    this%x_face(nx + ng + 1) = this%x_min + real(nx + ng, dp) * this%dx

    do j = -ng + 1, ny + ng
      this%y_cell(j) = this%y_min + (real(j, dp) - 0.5_dp) * this%dy
      this%y_face(j) = this%y_min + real(j - 1, dp) * this%dy
    end do
    this%y_face(ny + ng + 1) = this%y_min + real(ny + ng, dp) * this%dy

    this%prim = 0.0_dp
    this%cons = 0.0_dp
    this%bx_face = 0.0_dp
    this%by_face = 0.0_dp
    this%psi = 0.0_dp

    if (present(err)) call err%clear()
  end subroutine grid2d_init

  subroutine grid2d_free(this)
    class(grid2d_t), intent(inout) :: this

    if (allocated(this%x_cell)) deallocate(this%x_cell)
    if (allocated(this%y_cell)) deallocate(this%y_cell)
    if (allocated(this%x_face)) deallocate(this%x_face)
    if (allocated(this%y_face)) deallocate(this%y_face)

    if (allocated(this%prim)) deallocate(this%prim)
    if (allocated(this%cons)) deallocate(this%cons)
    if (allocated(this%bx_face)) deallocate(this%bx_face)
    if (allocated(this%by_face)) deallocate(this%by_face)
    if (allocated(this%psi)) deallocate(this%psi)
  end subroutine grid2d_free

  subroutine grid2d_apply_boundaries(this, err)
    class(grid2d_t), intent(inout) :: this
    type(umbra_error_t), intent(inout), optional :: err

    call apply_boundary_conditions_2d(this%prim, this%nx, this%ny, this%nghost, &
                                     this%bc_x_min, this%bc_x_max, this%bc_y_min, this%bc_y_max, err)
  end subroutine grid2d_apply_boundaries

end module umbra_grid2d
