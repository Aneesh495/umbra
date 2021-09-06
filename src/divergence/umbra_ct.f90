! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_ct
! Description: Staggered-grid Yee mesh Constrained Transport (CT) operator maintaining div(B) = 0.

module umbra_ct
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: EPSILON_DP
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  public :: compute_corner_electric_fields_2d
  public :: update_staggered_magnetic_fields_2d
  public :: compute_magnetic_divergence_2d
  public :: face_to_cell_magnetic_field_2d

contains

  subroutine compute_corner_electric_fields_2d(fx_by, fy_bx, nx, ny, ez, err)
    ! Computes corner-centered electric field Ez at (i-1/2, j-1/2) from face fluxes
    ! fx_by(i, j) is the flux of By in x-direction at interface (i-1/2, j)
    ! fy_bx(i, j) is the flux of Bx in y-direction at interface (i, j-1/2)
    real(dp), dimension(:,:), intent(in) :: fx_by
    real(dp), dimension(:,:), intent(in) :: fy_bx
    integer(i4), intent(in) :: nx, ny
    real(dp), dimension(:,:), intent(out) :: ez
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j

    ! 4-way upwind corner averaging (Gardiner & Stone 2005)
    do j = 1, ny + 1
      do i = 1, nx + 1
        ez(i, j) = 0.25_dp * ( &
          fx_by(i, max(j-1, 1)) + fx_by(i, min(j, ny)) &
        - fy_bx(max(i-1, 1), j) - fy_bx(min(i, nx), j) )
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine compute_corner_electric_fields_2d

  subroutine update_staggered_magnetic_fields_2d(ez, dt, dx, dy, nx, ny, bx_face, by_face, err)
    ! Updates face-centered magnetic fields Bx(i-1/2, j) and By(i, j-1/2) via discrete curl:
    ! dBx/dt = -dEz/dy
    ! dBy/dt =  dEz/dx
    real(dp), dimension(:,:), intent(in) :: ez
    real(dp), intent(in) :: dt, dx, dy
    integer(i4), intent(in) :: nx, ny
    real(dp), dimension(:,:), intent(inout) :: bx_face
    real(dp), dimension(:,:), intent(inout) :: by_face
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j
    real(dp) :: dt_over_dy, dt_over_dx

    dt_over_dy = dt / dy
    dt_over_dx = dt / dx

    ! Update Bx at face (i-1/2, j): dBx = -dt * (Ez(i, j+1) - Ez(i, j)) / dy
    do j = 1, ny
      do i = 1, nx + 1
        bx_face(i, j) = bx_face(i, j) - dt_over_dy * (ez(i, j + 1) - ez(i, j))
      end do
    end do

    ! Update By at face (i, j-1/2): dBy =  dt * (Ez(i+1, j) - Ez(i, j)) / dx
    do j = 1, ny + 1
      do i = 1, nx
        by_face(i, j) = by_face(i, j) + dt_over_dx * (ez(i + 1, j) - ez(i, j))
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine update_staggered_magnetic_fields_2d

  subroutine compute_magnetic_divergence_2d(bx_face, by_face, dx, dy, nx, ny, div_b, max_div, err)
    ! Evaluates discrete divergence: div B = (Bx(i+1, j) - Bx(i, j))/dx + (By(i, j+1) - By(i, j))/dy
    real(dp), dimension(:,:), intent(in) :: bx_face
    real(dp), dimension(:,:), intent(in) :: by_face
    real(dp), intent(in) :: dx, dy
    integer(i4), intent(in) :: nx, ny
    real(dp), dimension(:,:), intent(out) :: div_b
    real(dp), intent(out), optional :: max_div
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j
    real(dp) :: inv_dx, inv_dy, max_val

    inv_dx = 1.0_dp / dx
    inv_dy = 1.0_dp / dy
    max_val = 0.0_dp

    do j = 1, ny
      do i = 1, nx
        div_b(i, j) = (bx_face(i + 1, j) - bx_face(i, j)) * inv_dx &
                    + (by_face(i, j + 1) - by_face(i, j)) * inv_dy
        max_val = max(max_val, abs(div_b(i, j)))
      end do
    end do

    if (present(max_div)) max_div = max_val
    if (present(err)) call err%clear()
  end subroutine compute_magnetic_divergence_2d

  subroutine face_to_cell_magnetic_field_2d(bx_face, by_face, nx, ny, bx_cell, by_cell, err)
    ! Interpolates face-centered magnetic fields to cell centers via second-order averaging
    real(dp), dimension(:,:), intent(in) :: bx_face
    real(dp), dimension(:,:), intent(in) :: by_face
    integer(i4), intent(in) :: nx, ny
    real(dp), dimension(:,:), intent(out) :: bx_cell
    real(dp), dimension(:,:), intent(out) :: by_cell
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: i, j

    do j = 1, ny
      do i = 1, nx
        bx_cell(i, j) = 0.5_dp * (bx_face(i, j) + bx_face(i + 1, j))
        by_cell(i, j) = 0.5_dp * (by_face(i, j) + by_face(i, j + 1))
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine face_to_cell_magnetic_field_2d

end module umbra_ct
