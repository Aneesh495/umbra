!> VTK structured rectilinear mesh output writer for scientific visualization
!> Produces standard legacy VTK datasets for ParaView and VisIt visualization
module umbra_io_vtk
  use umbra_kinds, only: dp, i4
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_IO_FAILURE
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_P, PRIM_VX, PRIM_VY, PRIM_VZ, &
                         PRIM_BX, PRIM_BY, PRIM_BZ
  implicit none
  private

  public :: write_vtk_rectilinear_2d
  public :: write_vtk_rectilinear_3d

contains

  !> Writes 2D GRMHD state variables to standard legacy VTK rectilinear file
  subroutine write_vtk_rectilinear_2d(filename, title, nx, ny, x_coords, y_coords, prim_arr, err)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: title
    integer(i4), intent(in) :: nx, ny
    real(dp), intent(in) :: x_coords(nx)
    real(dp), intent(in) :: y_coords(ny)
    real(dp), intent(in) :: prim_arr(NPRIM, nx, ny)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat, i, j
    integer(i4) :: n_total

    n_total = nx * ny

    open(newunit=u, file=trim(filename), status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Unable to create VTK file")
      return
    end if

    ! VTK Header
    write(u, '(A)') "# vtk DataFile Version 3.0"
    write(u, '(A)') trim(title)
    write(u, '(A)') "ASCII"
    write(u, '(A)') "DATASET RECTILINEAR_GRID"
    write(u, '(A, 3(1X, I0))') "DIMENSIONS", nx, ny, 1

    ! Coordinate vectors
    write(u, '(A, 1X, I0, 1X, A)') "X_COORDINATES", nx, "double"
    do i = 1, nx
      write(u, '(E16.8)') x_coords(i)
    end do

    write(u, '(A, 1X, I0, 1X, A)') "Y_COORDINATES", ny, "double"
    do j = 1, ny
      write(u, '(E16.8)') y_coords(j)
    end do

    write(u, '(A, 1X, I0, 1X, A)') "Z_COORDINATES", 1, "double"
    write(u, '(E16.8)') 0.0_dp

    ! Point data fields
    write(u, '(A, 1X, I0)') "POINT_DATA", n_total

    ! Rest-mass density scalar
    write(u, '(A)') "SCALARS density double 1"
    write(u, '(A)') "LOOKUP_TABLE default"
    do j = 1, ny
      do i = 1, nx
        write(u, '(E16.8)') prim_arr(PRIM_RHO, i, j)
      end do
    end do

    ! Gas pressure scalar
    write(u, '(A)') "SCALARS pressure double 1"
    write(u, '(A)') "LOOKUP_TABLE default"
    do j = 1, ny
      do i = 1, nx
        write(u, '(E16.8)') prim_arr(PRIM_P, i, j)
      end do
    end do

    ! Velocity vector
    write(u, '(A)') "VECTORS velocity double"
    do j = 1, ny
      do i = 1, nx
        write(u, '(3(1X, E16.8))') prim_arr(PRIM_VX, i, j), &
                                   prim_arr(PRIM_VY, i, j), &
                                   prim_arr(PRIM_VZ, i, j)
      end do
    end do

    ! Magnetic field vector
    write(u, '(A)') "VECTORS magnetic_field double"
    do j = 1, ny
      do i = 1, nx
        write(u, '(3(1X, E16.8))') prim_arr(PRIM_BX, i, j), &
                                   prim_arr(PRIM_BY, i, j), &
                                   prim_arr(PRIM_BZ, i, j)
      end do
    end do

    close(u)
    if (present(err)) call err%clear()
  end subroutine write_vtk_rectilinear_2d

  !> Writes 3D GRMHD state variables to standard legacy VTK rectilinear file
  subroutine write_vtk_rectilinear_3d(filename, title, nx, ny, nz, x_coords, y_coords, z_coords, prim_arr, err)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: title
    integer(i4), intent(in) :: nx, ny, nz
    real(dp), intent(in) :: x_coords(nx)
    real(dp), intent(in) :: y_coords(ny)
    real(dp), intent(in) :: z_coords(nz)
    real(dp), intent(in) :: prim_arr(NPRIM, nx, ny, nz)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat, i, j, k
    integer(i4) :: n_total

    n_total = nx * ny * nz

    open(newunit=u, file=trim(filename), status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Unable to create 3D VTK file")
      return
    end if

    ! VTK Header
    write(u, '(A)') "# vtk DataFile Version 3.0"
    write(u, '(A)') trim(title)
    write(u, '(A)') "ASCII"
    write(u, '(A)') "DATASET RECTILINEAR_GRID"
    write(u, '(A, 3(1X, I0))') "DIMENSIONS", nx, ny, nz

    ! Coordinates
    write(u, '(A, 1X, I0, 1X, A)') "X_COORDINATES", nx, "double"
    do i = 1, nx
      write(u, '(E16.8)') x_coords(i)
    end do

    write(u, '(A, 1X, I0, 1X, A)') "Y_COORDINATES", ny, "double"
    do j = 1, ny
      write(u, '(E16.8)') y_coords(j)
    end do

    write(u, '(A, 1X, I0, 1X, A)') "Z_COORDINATES", nz, "double"
    do k = 1, nz
      write(u, '(E16.8)') z_coords(k)
    end do

    ! Point data fields
    write(u, '(A, 1X, I0)') "POINT_DATA", n_total

    ! Density
    write(u, '(A)') "SCALARS density double 1"
    write(u, '(A)') "LOOKUP_TABLE default"
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(u, '(E16.8)') prim_arr(PRIM_RHO, i, j, k)
        end do
      end do
    end do

    ! Gas pressure
    write(u, '(A)') "SCALARS pressure double 1"
    write(u, '(A)') "LOOKUP_TABLE default"
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(u, '(E16.8)') prim_arr(PRIM_P, i, j, k)
        end do
      end do
    end do

    ! Velocity vector
    write(u, '(A)') "VECTORS velocity double"
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(u, '(3(1X, E16.8))') prim_arr(PRIM_VX, i, j, k), &
                                     prim_arr(PRIM_VY, i, j, k), &
                                     prim_arr(PRIM_VZ, i, j, k)
        end do
      end do
    end do

    ! Magnetic field vector
    write(u, '(A)') "VECTORS magnetic_field double"
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          write(u, '(3(1X, E16.8))') prim_arr(PRIM_BX, i, j, k), &
                                     prim_arr(PRIM_BY, i, j, k), &
                                     prim_arr(PRIM_BZ, i, j, k)
        end do
      end do
    end do

    close(u)
    if (present(err)) call err%clear()
  end subroutine write_vtk_rectilinear_3d

end module umbra_io_vtk
