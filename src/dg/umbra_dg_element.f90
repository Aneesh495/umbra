! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_dg_element
! Description: Spectral element geometric mappings and volume polynomial representations.

module umbra_dg_element
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS, NPRIM
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_dg_basis, only: dg_basis_t
  implicit none
  private

  type, public :: dg_element_1d_t
    integer(i4) :: id = 0_i4
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: dx = 1.0_dp

    ! Physical coordinates of nodes x_i = 0.5 * (x_min + x_max) + 0.5 * dx * xi_i
    real(dp), allocatable, dimension(:) :: x_nodes

    ! State at collocation nodes
    real(dp), allocatable, dimension(:,:) :: u_cons ! (NCONS, num_nodes)
    real(dp), allocatable, dimension(:,:) :: p_prim ! (NPRIM, num_nodes)
    real(dp), allocatable, dimension(:,:) :: residual ! dU/dt (NCONS, num_nodes)
  contains
    procedure, public :: init => dg_element_1d_init
    procedure, public :: free => dg_element_1d_free
    procedure, public :: cell_average => dg_element_1d_cell_average
  end type dg_element_1d_t

contains

  subroutine dg_element_1d_init(this, elem_id, xmin, xmax, basis, err)
    class(dg_element_1d_t), intent(inout) :: this
    integer(i4), intent(in) :: elem_id
    real(dp), intent(in) :: xmin, xmax
    type(dg_basis_t), intent(in) :: basis
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: n, i

    this%id = elem_id
    this%x_min = xmin
    this%x_max = xmax
    this%dx = xmax - xmin
    n = basis%num_nodes

    call this%free()

    allocate(this%x_nodes(n))
    allocate(this%u_cons(NCONS, n))
    allocate(this%p_prim(NPRIM, n))
    allocate(this%residual(NCONS, n))

    do i = 1, n
      this%x_nodes(i) = 0.5_dp * (xmin + xmax) + 0.5_dp * this%dx * basis%nodes(i)
    end do

    this%u_cons = 0.0_dp
    this%p_prim = 0.0_dp
    this%residual = 0.0_dp

    if (present(err)) call err%clear()
  end subroutine dg_element_1d_init

  subroutine dg_element_1d_free(this)
    class(dg_element_1d_t), intent(inout) :: this
    if (allocated(this%x_nodes)) deallocate(this%x_nodes)
    if (allocated(this%u_cons)) deallocate(this%u_cons)
    if (allocated(this%p_prim)) deallocate(this%p_prim)
    if (allocated(this%residual)) deallocate(this%residual)
  end subroutine dg_element_1d_free

  pure function dg_element_1d_cell_average(this, basis, var) result(avg)
    class(dg_element_1d_t), intent(in) :: this
    type(dg_basis_t), intent(in) :: basis
    integer(i4), intent(in) :: var
    real(dp) :: avg

    integer(i4) :: i, n

    n = basis%num_nodes
    avg = 0.0_dp
    do i = 1, n
      avg = avg + 0.5_dp * basis%weights(i) * this%u_cons(var, i)
    end do
  end function dg_element_1d_cell_average

end module umbra_dg_element
