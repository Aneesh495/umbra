! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_grid1d
! Description: One-dimensional structured computational grid for shock tube and convergence benchmarks.

module umbra_grid1d
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_boundary, only: apply_boundary_conditions_1d, BC_OUTFLOW
  implicit none
  private

  type, public :: grid1d_t
    integer(i4) :: nx = 128_i4
    integer(i4) :: nghost = 3_i4
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: dx = 0.0_dp

    real(dp), allocatable, dimension(:) :: x_cell
    real(dp), allocatable, dimension(:) :: x_face

    real(dp), allocatable, dimension(:,:) :: prim
    real(dp), allocatable, dimension(:,:) :: cons

    integer(i4) :: bc_left = BC_OUTFLOW
    integer(i4) :: bc_right = BC_OUTFLOW
  contains
    procedure, public :: init => grid1d_init
    procedure, public :: free => grid1d_free
    procedure, public :: apply_boundaries => grid1d_apply_boundaries
  end type grid1d_t

contains

  subroutine grid1d_init(this, nx, xmin, xmax, nghost, err)
    class(grid1d_t), intent(inout) :: this
    integer(i4), intent(in) :: nx
    real(dp), intent(in) :: xmin, xmax
    integer(i4), intent(in), optional :: nghost
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: ng, i

    ng = 3_i4
    if (present(nghost)) ng = nghost

    this%nx = nx
    this%nghost = ng
    this%x_min = xmin
    this%x_max = xmax
    this%dx = (xmax - xmin) / real(nx, dp)

    call this%free()

    allocate(this%x_cell(-ng + 1 : nx + ng))
    allocate(this%x_face(-ng + 1 : nx + ng + 1))
    allocate(this%prim(NPRIM, -ng + 1 : nx + ng))
    allocate(this%cons(NCONS, -ng + 1 : nx + ng))

    do i = -ng + 1, nx + ng
      this%x_cell(i) = this%x_min + (real(i, dp) - 0.5_dp) * this%dx
      this%x_face(i) = this%x_min + real(i - 1, dp) * this%dx
    end do
    this%x_face(nx + ng + 1) = this%x_min + real(nx + ng, dp) * this%dx

    this%prim = 0.0_dp
    this%cons = 0.0_dp

    if (present(err)) call err%clear()
  end subroutine grid1d_init

  subroutine grid1d_free(this)
    class(grid1d_t), intent(inout) :: this

    if (allocated(this%x_cell)) deallocate(this%x_cell)
    if (allocated(this%x_face)) deallocate(this%x_face)
    if (allocated(this%prim)) deallocate(this%prim)
    if (allocated(this%cons)) deallocate(this%cons)
  end subroutine grid1d_free

  subroutine grid1d_apply_boundaries(this, err)
    class(grid1d_t), intent(inout) :: this
    type(umbra_error_t), intent(inout), optional :: err

    call apply_boundary_conditions_1d(this%prim, this%nx, this%nghost, &
                                     this%bc_left, this%bc_right, err)
  end subroutine grid1d_apply_boundaries

end module umbra_grid1d
