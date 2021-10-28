!> Formatted human-readable ASCII table and diagnostics writer
!> Generates column-oriented tables for 1D slices, time histories, and ASCII PPM images
module umbra_io_ascii
  use umbra_kinds, only: dp, i4
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_IO_FAILURE
  use umbra_constants, only: NPRIM, PRIM_RHO, PRIM_P, PRIM_VX, PRIM_VY, PRIM_VZ, &
                         PRIM_BX, PRIM_BY, PRIM_BZ
  implicit none
  private

  public :: write_slice_1d_ascii
  public :: write_diagnostics_row
  public :: write_diagnostics_header
  public :: write_ppm_image

contains

  !> Writes 1D spatial profile along x-axis to columnated text file
  subroutine write_slice_1d_ascii(filename, nx, x_coords, prim_arr, time_val, err)
    character(len=*), intent(in) :: filename
    integer(i4), intent(in) :: nx
    real(dp), intent(in) :: x_coords(nx)
    real(dp), intent(in) :: prim_arr(NPRIM, nx)
    real(dp), intent(in) :: time_val
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat, i

    open(newunit=u, file=trim(filename), status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open 1D slice ASCII file")
      return
    end if

    write(u, '(A, 1X, E16.8)') "# time =", time_val
    write(u, '(A)') "# 1:x 2:rho 3:p 4:vx 5:vy 6:vz 7:Bx 8:By 9:Bz"

    do i = 1, nx
      write(u, '(9(1X, E16.8))') x_coords(i), &
                                 prim_arr(PRIM_RHO, i), &
                                 prim_arr(PRIM_P, i), &
                                 prim_arr(PRIM_VX, i), &
                                 prim_arr(PRIM_VY, i), &
                                 prim_arr(PRIM_VZ, i), &
                                 prim_arr(PRIM_BX, i), &
                                 prim_arr(PRIM_BY, i), &
                                 prim_arr(PRIM_BZ, i)
    end do

    close(u)
    if (present(err)) call err%clear()
  end subroutine write_slice_1d_ascii

  !> Writes header for time-series diagnostics log
  subroutine write_diagnostics_header(unit_num)
    integer(i4), intent(in) :: unit_num
    write(unit_num, '(A)') "# cycle time dt mass_total energy_total max_div_b"
  end subroutine write_diagnostics_header

  !> Writes a single time-series row to diagnostics file
  subroutine write_diagnostics_row(unit_num, cycle, time_val, dt, mass, energy, div_b)
    integer(i4), intent(in) :: unit_num
    integer(i4), intent(in) :: cycle
    real(dp), intent(in) :: time_val, dt, mass, energy, div_b

    write(unit_num, '(I8, 5(1X, E16.8))') cycle, time_val, dt, mass, energy, div_b
  end subroutine write_diagnostics_row

  !> Writes 2D ray-traced intensity map to standard ASCII PPM image format (P3)
  subroutine write_ppm_image(filename, npix_x, npix_y, intensity_map, err)
    character(len=*), intent(in) :: filename
    integer(i4), intent(in) :: npix_x, npix_y
    real(dp), intent(in) :: intensity_map(npix_x, npix_y)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat, i, j, val_int
    real(dp) :: max_val, min_val, range_val, normalized

    open(newunit=u, file=trim(filename), status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open PPM image file")
      return
    end if

    max_val = maxval(intensity_map)
    min_val = minval(intensity_map)
    range_val = max(max_val - min_val, 1.0e-12_dp)

    ! Netpbm ASCII PPM header
    write(u, '(A)') "P3"
    write(u, '(I0, 1X, I0)') npix_x, npix_y
    write(u, '(A)') "255"

    do j = npix_y, 1, -1 ! Invert y for screen raster orientation
      do i = 1, npix_x
        normalized = (intensity_map(i, j) - min_val) / range_val
        val_int = max(0_i4, min(255_i4, int(255.0_dp * normalized, i4)))
        ! Grayscale RGB
        write(u, '(3(1X, I3))') val_int, val_int, val_int
      end do
    end do

    close(u)
    if (present(err)) call err%clear()
  end subroutine write_ppm_image

end module umbra_io_ascii
