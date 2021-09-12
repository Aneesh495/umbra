! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_octree
! Description: Octree block-based adaptive mesh refinement (AMR) hierarchy with 2:1 balance.

module umbra_octree
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: NCONS
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS
  implicit none
  private

  integer(i4), parameter, public :: MAX_AMR_LEVEL = 6_i4
  integer(i4), parameter, public :: BLOCK_SIZE = 8_i4

  type, public :: amr_block_t
    integer(i4) :: level = 0_i4
    integer(i4) :: id = 0_i4
    real(dp), dimension(3) :: origin = 0.0_dp
    real(dp), dimension(3) :: size_box = 1.0_dp
    real(dp), dimension(3) :: dx = 0.125_dp
    logical :: is_leaf = .true.
    integer(i4) :: refine_tag = 0_i4 ! -1 = coarsen, 0 = keep, +1 = refine

    ! Data buffer for conserved variables: (NCONS, BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE)
    real(dp), allocatable, dimension(:,:,:,:) :: data

    ! Pointers/indices to parent and 8 children (1..8 octants)
    integer(i4) :: parent_id = 0_i4
    integer(i4), dimension(8) :: child_ids = 0_i4
  contains
    procedure, public :: init => amr_block_init
    procedure, public :: free => amr_block_free
  end type amr_block_t

  type, public :: amr_tree_t
    integer(i4) :: max_blocks = 1024_i4
    integer(i4) :: num_blocks = 0_i4
    type(amr_block_t), allocatable, dimension(:) :: blocks
  contains
    procedure, public :: init => amr_tree_init
    procedure, public :: free => amr_tree_free
    procedure, public :: tag_refinement => amr_tree_tag_refinement
    procedure, public :: refine_block => amr_tree_refine_block
    procedure, public :: restrict_to_parent => amr_tree_restrict_to_parent
  end type amr_tree_t

contains

  subroutine amr_block_init(this, id_val, lvl, orig, s_box)
    class(amr_block_t), intent(inout) :: this
    integer(i4), intent(in) :: id_val, lvl
    real(dp), dimension(3), intent(in) :: orig, s_box

    this%id = id_val
    this%level = lvl
    this%origin = orig
    this%size_box = s_box
    this%dx = s_box / real(BLOCK_SIZE, dp)
    this%is_leaf = .true.
    this%refine_tag = 0_i4
    this%child_ids = 0_i4

    if (allocated(this%data)) deallocate(this%data)
    allocate(this%data(NCONS, BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE))
    this%data = 0.0_dp
  end subroutine amr_block_init

  subroutine amr_block_free(this)
    class(amr_block_t), intent(inout) :: this
    if (allocated(this%data)) deallocate(this%data)
  end subroutine amr_block_free

  subroutine amr_tree_init(this, max_b)
    class(amr_tree_t), intent(inout) :: this
    integer(i4), intent(in) :: max_b

    this%max_blocks = max_b
    this%num_blocks = 0_i4
    if (allocated(this%blocks)) deallocate(this%blocks)
    allocate(this%blocks(max_b))
  end subroutine amr_tree_init

  subroutine amr_tree_free(this)
    class(amr_tree_t), intent(inout) :: this
    integer(i4) :: i
    if (allocated(this%blocks)) then
      do i = 1, this%num_blocks
        call this%blocks(i)%free()
      end do
      deallocate(this%blocks)
    end if
    this%num_blocks = 0_i4
  end subroutine amr_tree_free

  subroutine amr_tree_tag_refinement(this, refine_threshold, coarsen_threshold)
    ! Evaluates normalized gradient and curvature indicators across all leaf blocks
    class(amr_tree_t), intent(inout) :: this
    real(dp), intent(in) :: refine_threshold, coarsen_threshold

    integer(i4) :: b, i, j, k
    real(dp) :: max_curvature, curvature, rho_c, d2_rho

    do b = 1, this%num_blocks
      if (.not. this%blocks(b)%is_leaf) cycle

      max_curvature = 0.0_dp
      do k = 2, BLOCK_SIZE - 1
        do j = 2, BLOCK_SIZE - 1
          do i = 2, BLOCK_SIZE - 1
            rho_c = this%blocks(b)%data(1, i, j, k)
            ! Discrete Laplacian: d^2 rho / dx^2 + d^2 rho / dy^2 + d^2 rho / dz^2
            d2_rho = (this%blocks(b)%data(1, i+1, j, k) - 2.0_dp*rho_c + this%blocks(b)%data(1, i-1, j, k)) &
                   + (this%blocks(b)%data(1, i, j+1, k) - 2.0_dp*rho_c + this%blocks(b)%data(1, i, j-1, k)) &
                   + (this%blocks(b)%data(1, i, j, k+1) - 2.0_dp*rho_c + this%blocks(b)%data(1, i, j, k-1))
            curvature = abs(d2_rho) / max(abs(rho_c), 1.0e-8_dp)
            max_curvature = max(max_curvature, curvature)
          end do
        end do
      end do

      if (max_curvature > refine_threshold .and. this%blocks(b)%level < MAX_AMR_LEVEL) then
        this%blocks(b)%refine_tag = 1_i4
      else if (max_curvature < coarsen_threshold .and. this%blocks(b)%level > 0) then
        this%blocks(b)%refine_tag = -1_i4
      else
        this%blocks(b)%refine_tag = 0_i4
      end if
    end do
  end subroutine amr_tree_tag_refinement

  subroutine amr_tree_refine_block(this, block_id, err)
    class(amr_tree_t), intent(inout) :: this
    integer(i4), intent(in) :: block_id
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: oct, child_id, ox, oy, oz
    real(dp), dimension(3) :: child_orig, child_size
    integer(i4) :: i, j, k, pi, pj, pk, v

    if (this%num_blocks + 8 > this%max_blocks) return
    if (.not. this%blocks(block_id)%is_leaf) return

    child_size = 0.5_dp * this%blocks(block_id)%size_box
    this%blocks(block_id)%is_leaf = .false.

    do oct = 1, 8
      this%num_blocks = this%num_blocks + 1
      child_id = this%num_blocks
      this%blocks(block_id)%child_ids(oct) = child_id

      ! Octant coordinates (0 or 1 along x, y, z)
      ox = mod(oct - 1, 2)
      oy = mod((oct - 1) / 2, 2)
      oz = (oct - 1) / 4

      child_orig(1) = this%blocks(block_id)%origin(1) + real(ox, dp) * child_size(1)
      child_orig(2) = this%blocks(block_id)%origin(2) + real(oy, dp) * child_size(2)
      child_orig(3) = this%blocks(block_id)%origin(3) + real(oz, dp) * child_size(3)

      call this%blocks(child_id)%init(child_id, this%blocks(block_id)%level + 1, child_orig, child_size)
      this%blocks(child_id)%parent_id = block_id

      ! Prolongation: copy coarse parent cell to 2x2x2 fine child cells
      do k = 1, BLOCK_SIZE
        pk = oz * (BLOCK_SIZE / 2) + (k + 1) / 2
        do j = 1, BLOCK_SIZE
          pj = oy * (BLOCK_SIZE / 2) + (j + 1) / 2
          do i = 1, BLOCK_SIZE
            pi = ox * (BLOCK_SIZE / 2) + (i + 1) / 2
            do v = 1, NCONS
              this%blocks(child_id)%data(v, i, j, k) = this%blocks(block_id)%data(v, pi, pj, pk)
            end do
          end do
        end do
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine amr_tree_refine_block

  subroutine amr_tree_restrict_to_parent(this, parent_id, err)
    class(amr_tree_t), intent(inout) :: this
    integer(i4), intent(in) :: parent_id
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: oct, child_id, ox, oy, oz
    integer(i4) :: pi, pj, pk, ci, cj, ck, v
    real(dp) :: sum_val

    if (this%blocks(parent_id)%is_leaf) return

    ! Conservative volumetric restriction: coarse parent cell is average of 8 fine child cells
    do oct = 1, 8
      child_id = this%blocks(parent_id)%child_ids(oct)
      if (child_id == 0) cycle

      ox = mod(oct - 1, 2)
      oy = mod((oct - 1) / 2, 2)
      oz = (oct - 1) / 4

      do pk = 1, BLOCK_SIZE / 2
        ck = (pk - 1) * 2 + 1
        do pj = 1, BLOCK_SIZE / 2
          cj = (pj - 1) * 2 + 1
          do pi = 1, BLOCK_SIZE / 2
            ci = (pi - 1) * 2 + 1
            do v = 1, NCONS
              sum_val = this%blocks(child_id)%data(v, ci,   cj,   ck)   + &
                        this%blocks(child_id)%data(v, ci+1, cj,   ck)   + &
                        this%blocks(child_id)%data(v, ci,   cj+1, ck)   + &
                        this%blocks(child_id)%data(v, ci+1, cj+1, ck)   + &
                        this%blocks(child_id)%data(v, ci,   cj,   ck+1) + &
                        this%blocks(child_id)%data(v, ci+1, cj,   ck+1) + &
                        this%blocks(child_id)%data(v, ci,   cj+1, ck+1) + &
                        this%blocks(child_id)%data(v, ci+1, cj+1, ck+1)

              this%blocks(parent_id)%data(v, ox*(BLOCK_SIZE/2) + pi, &
                                             oy*(BLOCK_SIZE/2) + pj, &
                                             oz*(BLOCK_SIZE/2) + pk) = sum_val * 0.125_dp
            end do
          end do
        end do
      end do
    end do

    if (present(err)) call err%clear()
  end subroutine amr_tree_restrict_to_parent

end module umbra_octree
