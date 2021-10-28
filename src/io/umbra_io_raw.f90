!> High-throughput raw binary stream I/O for GRMHD simulations
!> Provides portable stream read/write routines for grid fields and camera frames
module umbra_io_raw
  use umbra_kinds, only: dp, i4
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_IO_FAILURE
  implicit none
  private

  public :: write_raw_field_2d
  public :: read_raw_field_2d
  public :: write_raw_field_3d
  public :: read_raw_field_3d
  public :: write_raw_image
  public :: read_raw_image

  character(len=4), parameter :: MAGIC_HEADER = "UMBR"
  integer(i4), parameter :: FILE_VERSION = 1_i4

contains

  !> Writes a 2D scalar or vector field to unformatted stream binary file
  subroutine write_raw_field_2d(filename, nvars, nx, ny, time_val, data_arr, err)
    character(len=*), intent(in) :: filename
    integer(i4), intent(in) :: nvars, nx, ny
    real(dp), intent(in) :: time_val
    real(dp), intent(in) :: data_arr(nvars, nx, ny)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open raw binary file for writing")
      return
    end if

    write(u, iostat=io_stat) MAGIC_HEADER
    write(u, iostat=io_stat) FILE_VERSION
    write(u, iostat=io_stat) nvars, nx, ny
    write(u, iostat=io_stat) time_val
    write(u, iostat=io_stat) data_arr

    close(u, iostat=io_stat)
    if (present(err)) call err%clear()
  end subroutine write_raw_field_2d

  !> Reads a 2D scalar or vector field from unformatted stream binary file
  subroutine read_raw_field_2d(filename, nvars, nx, ny, time_val, data_arr, err)
    character(len=*), intent(in) :: filename
    integer(i4), intent(out) :: nvars, nx, ny
    real(dp), intent(out) :: time_val
    real(dp), allocatable, intent(out) :: data_arr(:,:,:)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat, version
    character(len=4) :: magic

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='old', action='read', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open raw binary file for reading")
      return
    end if

    read(u, iostat=io_stat) magic
    if (magic /= MAGIC_HEADER) then
      close(u)
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Invalid magic header in raw binary file")
      return
    end if

    read(u, iostat=io_stat) version
    read(u, iostat=io_stat) nvars, nx, ny
    read(u, iostat=io_stat) time_val

    allocate(data_arr(nvars, nx, ny))
    read(u, iostat=io_stat) data_arr

    close(u, iostat=io_stat)
    if (present(err)) call err%clear()
  end subroutine read_raw_field_2d

  !> Writes a 3D scalar or vector field to unformatted stream binary file
  subroutine write_raw_field_3d(filename, nvars, nx, ny, nz, time_val, data_arr, err)
    character(len=*), intent(in) :: filename
    integer(i4), intent(in) :: nvars, nx, ny, nz
    real(dp), intent(in) :: time_val
    real(dp), intent(in) :: data_arr(nvars, nx, ny, nz)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open 3D raw binary file for writing")
      return
    end if

    write(u, iostat=io_stat) MAGIC_HEADER
    write(u, iostat=io_stat) FILE_VERSION
    write(u, iostat=io_stat) nvars, nx, ny, nz
    write(u, iostat=io_stat) time_val
    write(u, iostat=io_stat) data_arr

    close(u, iostat=io_stat)
    if (present(err)) call err%clear()
  end subroutine write_raw_field_3d

  !> Reads a 3D scalar or vector field from unformatted stream binary file
  subroutine read_raw_field_3d(filename, nvars, nx, ny, nz, time_val, data_arr, err)
    character(len=*), intent(in) :: filename
    integer(i4), intent(out) :: nvars, nx, ny, nz
    real(dp), intent(out) :: time_val
    real(dp), allocatable, intent(out) :: data_arr(:,:,:,:)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat, version
    character(len=4) :: magic

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='old', action='read', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open 3D raw binary file for reading")
      return
    end if

    read(u, iostat=io_stat) magic
    if (magic /= MAGIC_HEADER) then
      close(u)
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Invalid magic header in raw binary file")
      return
    end if

    read(u, iostat=io_stat) version
    read(u, iostat=io_stat) nvars, nx, ny, nz
    read(u, iostat=io_stat) time_val

    allocate(data_arr(nvars, nx, ny, nz))
    read(u, iostat=io_stat) data_arr

    close(u, iostat=io_stat)
    if (present(err)) call err%clear()
  end subroutine read_raw_field_3d

  !> Writes a ray-traced intensity image to unformatted stream binary file
  subroutine write_raw_image(filename, npix_x, npix_y, fov_m, image_data, err)
    character(len=*), intent(in) :: filename
    integer(i4), intent(in) :: npix_x, npix_y
    real(dp), intent(in) :: fov_m
    real(dp), intent(in) :: image_data(npix_x, npix_y)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open image file for writing")
      return
    end if

    write(u, iostat=io_stat) "UIMG"
    write(u, iostat=io_stat) FILE_VERSION
    write(u, iostat=io_stat) npix_x, npix_y
    write(u, iostat=io_stat) fov_m
    write(u, iostat=io_stat) image_data

    close(u, iostat=io_stat)
    if (present(err)) call err%clear()
  end subroutine write_raw_image

  !> Reads a ray-traced intensity image from unformatted stream binary file
  subroutine read_raw_image(filename, npix_x, npix_y, fov_m, image_data, err)
    character(len=*), intent(in) :: filename
    integer(i4), intent(out) :: npix_x, npix_y
    real(dp), intent(out) :: fov_m
    real(dp), allocatable, intent(out) :: image_data(:,:)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat, version
    character(len=4) :: magic

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='old', action='read', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open image file for reading")
      return
    end if

    read(u, iostat=io_stat) magic
    if (magic /= "UIMG") then
      close(u)
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Invalid image magic header")
      return
    end if

    read(u, iostat=io_stat) version
    read(u, iostat=io_stat) npix_x, npix_y
    read(u, iostat=io_stat) fov_m

    allocate(image_data(npix_x, npix_y))
    read(u, iostat=io_stat) image_data

    close(u, iostat=io_stat)
    if (present(err)) call err%clear()
  end subroutine read_raw_image

end module umbra_io_raw
