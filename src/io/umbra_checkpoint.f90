!> Checkpoint and restart manager with checksum integrity verification
!> Serializes complete simulation state to stream binary files for disaster recovery
module umbra_checkpoint
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NPRIM, NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_IO_FAILURE, UMBRA_ERR_GENERIC
  implicit none
  private

  public :: checkpoint_meta_t
  public :: write_checkpoint_2d
  public :: read_checkpoint_2d
  public :: compute_buffer_checksum

  character(len=8), parameter :: CHECKPOINT_MAGIC = "UMBRACHK"
  integer(i4), parameter :: CHECKPOINT_VERSION = 2_i4

  type, public :: checkpoint_meta_t
    integer(i4) :: cycle = 0
    real(dp) :: time_val = 0.0_dp
    real(dp) :: dt = 0.0_dp
    integer(i4) :: nx = 0
    integer(i4) :: ny = 0
    integer(i4) :: nz = 1
    real(dp) :: x_min = 0.0_dp
    real(dp) :: x_max = 1.0_dp
    real(dp) :: y_min = 0.0_dp
    real(dp) :: y_max = 1.0_dp
    real(dp) :: gamma_eos = 1.3333333333333333_dp
    integer(i4) :: checksum = 0
  end type checkpoint_meta_t

contains

  !> Computes 32-bit bitwise polynomial checksum over 1D real buffer
  pure function compute_buffer_checksum(buffer, n_elems) result(cs)
    integer(i4), intent(in) :: n_elems
    real(dp), intent(in) :: buffer(n_elems)
    integer(i4) :: cs
    integer(i4) :: i
    integer(i4) :: byte_rep(2)

    cs = 123456789_i4
    do i = 1, n_elems
      byte_rep = transfer(buffer(i), byte_rep)
      cs = ieor(cs, byte_rep(1)) + byte_rep(2)
      cs = ishftc(cs, 5)
    end do
  end function compute_buffer_checksum

  !> Writes complete 2D simulation checkpoint
  subroutine write_checkpoint_2d(filename, meta, prim_arr, cons_arr, psi_arr, err)
    character(len=*), intent(in) :: filename
    type(checkpoint_meta_t), intent(inout) :: meta
    real(dp), intent(in) :: prim_arr(NPRIM, meta%nx, meta%ny)
    real(dp), intent(in) :: cons_arr(NCONS, meta%nx, meta%ny)
    real(dp), intent(in) :: psi_arr(meta%nx, meta%ny)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat
    integer(i4) :: n_total_prim

    n_total_prim = NPRIM * meta%nx * meta%ny
    meta%checksum = compute_buffer_checksum(reshape(prim_arr, [n_total_prim]), n_total_prim)

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open checkpoint file for writing")
      return
    end if

    write(u, iostat=io_stat) CHECKPOINT_MAGIC
    write(u, iostat=io_stat) CHECKPOINT_VERSION
    write(u, iostat=io_stat) meta%cycle
    write(u, iostat=io_stat) meta%time_val
    write(u, iostat=io_stat) meta%dt
    write(u, iostat=io_stat) meta%nx, meta%ny, meta%nz
    write(u, iostat=io_stat) meta%x_min, meta%x_max, meta%y_min, meta%y_max
    write(u, iostat=io_stat) meta%gamma_eos
    write(u, iostat=io_stat) meta%checksum

    write(u, iostat=io_stat) prim_arr
    write(u, iostat=io_stat) cons_arr
    write(u, iostat=io_stat) psi_arr

    close(u, iostat=io_stat)
    if (present(err)) call err%clear()
  end subroutine write_checkpoint_2d

  !> Reads and validates 2D simulation checkpoint
  subroutine read_checkpoint_2d(filename, meta, prim_arr, cons_arr, psi_arr, err)
    character(len=*), intent(in) :: filename
    type(checkpoint_meta_t), intent(out) :: meta
    real(dp), allocatable, intent(out) :: prim_arr(:,:,:)
    real(dp), allocatable, intent(out) :: cons_arr(:,:,:)
    real(dp), allocatable, intent(out) :: psi_arr(:,:)
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: u, io_stat, version, stored_checksum, computed_checksum
    character(len=8) :: magic
    integer(i4) :: n_total_prim

    open(newunit=u, file=trim(filename), access='stream', form='unformatted', &
         status='old', action='read', iostat=io_stat)
    if (io_stat /= 0) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Failed to open checkpoint file for reading")
      return
    end if

    read(u, iostat=io_stat) magic
    if (magic /= CHECKPOINT_MAGIC) then
      close(u)
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Corrupt checkpoint magic header")
      return
    end if

    read(u, iostat=io_stat) version
    if (version /= CHECKPOINT_VERSION) then
      close(u)
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Unsupported checkpoint file version")
      return
    end if

    read(u, iostat=io_stat) meta%cycle
    read(u, iostat=io_stat) meta%time_val
    read(u, iostat=io_stat) meta%dt
    read(u, iostat=io_stat) meta%nx, meta%ny, meta%nz
    read(u, iostat=io_stat) meta%x_min, meta%x_max, meta%y_min, meta%y_max
    read(u, iostat=io_stat) meta%gamma_eos
    read(u, iostat=io_stat) stored_checksum
    meta%checksum = stored_checksum

    allocate(prim_arr(NPRIM, meta%nx, meta%ny))
    allocate(cons_arr(NCONS, meta%nx, meta%ny))
    allocate(psi_arr(meta%nx, meta%ny))

    read(u, iostat=io_stat) prim_arr
    read(u, iostat=io_stat) cons_arr
    read(u, iostat=io_stat) psi_arr

    close(u)

    n_total_prim = NPRIM * meta%nx * meta%ny
    computed_checksum = compute_buffer_checksum(reshape(prim_arr, [n_total_prim]), n_total_prim)

    if (computed_checksum /= stored_checksum) then
      if (present(err)) call err%set(UMBRA_ERR_IO_FAILURE, "Checksum mismatch: checkpoint data corrupted")
      return
    end if

    if (present(err)) call err%clear()
  end subroutine read_checkpoint_2d

end module umbra_checkpoint
