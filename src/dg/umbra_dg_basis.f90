! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_dg_basis
! Description: Discontinuous Galerkin spectral element polynomial basis and differentiation operators.

module umbra_dg_basis
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: EPSILON_DP
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  use umbra_quadrature, only: gauss_legendre_nodes_weights, gauss_lobatto_nodes_weights
  use umbra_polynomials, only: nodal_differentiation_matrix, lagrange_interpolant
  implicit none
  private

  type, public :: dg_basis_t
    integer(i4) :: order = 3_i4                ! Polynomial degree N
    integer(i4) :: num_nodes = 4_i4            ! Number of nodal collocation points N + 1
    logical     :: use_lobatto = .true.        ! True for LGL, False for GL nodes

    real(dp), allocatable, dimension(:) :: nodes   ! xi_i in [-1, 1]
    real(dp), allocatable, dimension(:) :: weights ! w_i
    real(dp), allocatable, dimension(:,:) :: d_mat ! D_ij = ell'_j(xi_i)

    ! Boundary evaluations: b_left(j) = ell_j(-1), b_right(j) = ell_j(+1)
    real(dp), allocatable, dimension(:) :: b_left
    real(dp), allocatable, dimension(:) :: b_right
  contains
    procedure, public :: init => dg_basis_init
    procedure, public :: free => dg_basis_free
  end type dg_basis_t

contains

  subroutine dg_basis_init(this, order, use_lobatto, err)
    class(dg_basis_t), intent(inout) :: this
    integer(i4), intent(in) :: order
    logical, intent(in), optional :: use_lobatto
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: n, j

    this%order = order
    this%num_nodes = order + 1_i4
    this%use_lobatto = .true.
    if (present(use_lobatto)) this%use_lobatto = use_lobatto

    call this%free()

    n = this%num_nodes
    allocate(this%nodes(n))
    allocate(this%weights(n))
    allocate(this%d_mat(n, n))
    allocate(this%b_left(n))
    allocate(this%b_right(n))

    if (this%use_lobatto) then
      call gauss_lobatto_nodes_weights(n, this%nodes, this%weights, err)
    else
      call gauss_legendre_nodes_weights(n, this%nodes, this%weights, err)
    end if

    call nodal_differentiation_matrix(this%nodes, this%d_mat, err)

    ! Evaluate boundary Lagrange interpolants at xi = -1 and xi = +1
    do j = 1, n
      this%b_left(j)  = lagrange_interpolant(j, -1.0_dp, this%nodes)
      this%b_right(j) = lagrange_interpolant(j,  1.0_dp, this%nodes)
    end do

    if (present(err)) call err%clear()
  end subroutine dg_basis_init

  subroutine dg_basis_free(this)
    class(dg_basis_t), intent(inout) :: this
    if (allocated(this%nodes)) deallocate(this%nodes)
    if (allocated(this%weights)) deallocate(this%weights)
    if (allocated(this%d_mat)) deallocate(this%d_mat)
    if (allocated(this%b_left)) deallocate(this%b_left)
    if (allocated(this%b_right)) deallocate(this%b_right)
  end subroutine dg_basis_free

end module umbra_dg_basis
