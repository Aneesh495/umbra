! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_boundary
! Description: Boundary condition handlers for relativistic fluid and electromagnetic fields.

module umbra_boundary
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, PRIM_VX, PRIM_VY, PRIM_VZ, PRIM_BX, PRIM_BY, PRIM_BZ
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  integer(i4), parameter, public :: BC_PERIODIC           = 1_i4
  integer(i4), parameter, public :: BC_OUTFLOW            = 2_i4
  integer(i4), parameter, public :: BC_REFLECTIVE         = 3_i4
  integer(i4), parameter, public :: BC_DIRICHLET          = 4_i4
  integer(i4), parameter, public :: BC_HORIZON_ABSORBING  = 5_i4

  public :: apply_boundary_conditions_1d
  public :: apply_boundary_conditions_2d

contains

  subroutine apply_boundary_conditions_1d(u, nx, nghost, bc_left, bc_right, err)
    integer(i4), intent(in) :: nx, nghost
    real(dp), dimension(NPRIM, 1-nghost : nx+nghost), intent(inout) :: u
    integer(i4), intent(in) :: bc_left, bc_right
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: g

    ! Left boundary: g in 1..nghost -> ghost cell 1-g
    select case (bc_left)
      case (BC_PERIODIC)
        do g = 1, nghost
          u(:, 1 - g) = u(:, nx + 1 - g)
        end do

      case (BC_OUTFLOW, BC_HORIZON_ABSORBING)
        do g = 1, nghost
          u(:, 1 - g) = u(:, 1)
        end do

      case (BC_REFLECTIVE)
        do g = 1, nghost
          u(:, 1 - g) = u(:, g)
          ! Invert normal velocity and normal magnetic field
          u(PRIM_VX, 1 - g) = -u(PRIM_VX, g)
          u(PRIM_BX, 1 - g) = -u(PRIM_BX, g)
        end do
    end select

    ! Right boundary: g in 1..nghost -> ghost cell nx+g
    select case (bc_right)
      case (BC_PERIODIC)
        do g = 1, nghost
          u(:, nx + g) = u(:, g)
        end do

      case (BC_OUTFLOW, BC_HORIZON_ABSORBING)
        do g = 1, nghost
          u(:, nx + g) = u(:, nx)
        end do

      case (BC_REFLECTIVE)
        do g = 1, nghost
          u(:, nx + g) = u(:, nx + 1 - g)
          u(PRIM_VX, nx + g) = -u(PRIM_VX, nx + 1 - g)
          u(PRIM_BX, nx + g) = -u(PRIM_BX, nx + 1 - g)
        end do
    end select

    if (present(err)) call err%clear()
  end subroutine apply_boundary_conditions_1d

  subroutine apply_boundary_conditions_2d(u, nx, ny, nghost, bc_x_min, bc_x_max, bc_y_min, bc_y_max, err)
    integer(i4), intent(in) :: nx, ny, nghost
    real(dp), dimension(NPRIM, 1-nghost : nx+nghost, 1-nghost : ny+nghost), intent(inout) :: u
    integer(i4), intent(in) :: bc_x_min, bc_x_max, bc_y_min, bc_y_max
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j, g

    ! X-min boundary
    select case (bc_x_min)
      case (BC_PERIODIC)
        do j = 1, ny
          do g = 1, nghost
            u(:, 1 - g, j) = u(:, nx + 1 - g, j)
          end do
        end do
      case (BC_OUTFLOW, BC_HORIZON_ABSORBING)
        do j = 1, ny
          do g = 1, nghost
            u(:, 1 - g, j) = u(:, 1, j)
          end do
        end do
      case (BC_REFLECTIVE)
        do j = 1, ny
          do g = 1, nghost
            u(:, 1 - g, j) = u(:, g, j)
            u(PRIM_VX, 1 - g, j) = -u(PRIM_VX, g, j)
            u(PRIM_BX, 1 - g, j) = -u(PRIM_BX, g, j)
          end do
        end do
    end select

    ! X-max boundary
    select case (bc_x_max)
      case (BC_PERIODIC)
        do j = 1, ny
          do g = 1, nghost
            u(:, nx + g, j) = u(:, g, j)
          end do
        end do
      case (BC_OUTFLOW, BC_HORIZON_ABSORBING)
        do j = 1, ny
          do g = 1, nghost
            u(:, nx + g, j) = u(:, nx, j)
          end do
        end do
      case (BC_REFLECTIVE)
        do j = 1, ny
          do g = 1, nghost
            u(:, nx + g, j) = u(:, nx + 1 - g, j)
            u(PRIM_VX, nx + g, j) = -u(PRIM_VX, nx + 1 - g, j)
            u(PRIM_BX, nx + g, j) = -u(PRIM_BX, nx + 1 - g, j)
          end do
        end do
    end select

    ! Y-min boundary
    select case (bc_y_min)
      case (BC_PERIODIC)
        do i = -nghost + 1, nx + nghost
          do g = 1, nghost
            u(:, i, 1 - g) = u(:, i, ny + 1 - g)
          end do
        end do
      case (BC_OUTFLOW, BC_HORIZON_ABSORBING)
        do i = -nghost + 1, nx + nghost
          do g = 1, nghost
            u(:, i, 1 - g) = u(:, i, 1)
          end do
        end do
      case (BC_REFLECTIVE)
        do i = -nghost + 1, nx + nghost
          do g = 1, nghost
            u(:, i, 1 - g) = u(:, i, g)
            u(PRIM_VY, i, 1 - g) = -u(PRIM_VY, i, g)
            u(PRIM_BY, i, 1 - g) = -u(PRIM_BY, i, g)
          end do
        end do
    end select

    ! Y-max boundary
    select case (bc_y_max)
      case (BC_PERIODIC)
        do i = -nghost + 1, nx + nghost
          do g = 1, nghost
            u(:, i, ny + g) = u(:, i, g)
          end do
        end do
      case (BC_OUTFLOW, BC_HORIZON_ABSORBING)
        do i = -nghost + 1, nx + nghost
          do g = 1, nghost
            u(:, i, ny + g) = u(:, i, ny)
          end do
        end do
      case (BC_REFLECTIVE)
        do i = -nghost + 1, nx + nghost
          do g = 1, nghost
            u(:, i, ny + g) = u(:, i, ny + 1 - g)
            u(PRIM_VY, i, ny + g) = -u(PRIM_VY, i, ny + 1 - g)
            u(PRIM_BY, i, ny + g) = -u(PRIM_BY, i, ny + 1 - g)
          end do
        end do
    end select

    if (present(err)) call err%clear()
  end subroutine apply_boundary_conditions_2d

end module umbra_boundary
