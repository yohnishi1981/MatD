  private

  public :: matd_matrix
  type matd_matrix
    private
    integer :: nprocs
    integer :: myrank
    integer :: win
    integer :: dim1
    integer :: dim2
    integer :: nblocks1
    integer :: nblocks2
    integer :: nblocks
    integer, pointer :: map1(:)
    integer, pointer :: map2(:)
    MATD_ELEMTYPE, pointer :: storage(:)
    logical :: is_scalapack
    integer :: procgrid1
    integer :: procgrid2
    integer :: nprocgrids1
    integer :: nprocgrids2
    integer :: nprocgrids
  end type matd_matrix


  public :: matd_create_reg
  interface matd_create_reg
    module procedure matd_create_reg_
  end interface matd_create_reg

  public :: matd_create_irreg
  interface matd_create_irreg
    module procedure matd_create_irreg_
  end interface matd_create_irreg

  public :: matd_create_scalapack
  interface matd_create_scalapack
    module procedure matd_create_scalapack_
  end interface matd_create_scalapack

  public :: matd_put
  interface matd_put
    module procedure matd_put_
  end interface matd_put

  public :: matd_get
  interface matd_get
    module procedure matd_get_
  end interface matd_get

  public :: matd_fence
  interface matd_fence
    module procedure matd_fence_
  end interface matd_fence

  public :: matd_destroy
  interface matd_destroy
    module procedure matd_destroy_
  end interface matd_destroy

  public :: matd_copy
  interface matd_copy
    module procedure matd_copy_
  end interface matd_copy

  public :: matd_print
  interface matd_print
    module procedure matd_print_
  end interface matd_print

  public :: matd_print_info
  interface matd_print_info
    module procedure matd_print_info_
  end interface matd_print_info

  public :: matd_locate
  interface matd_locate
    module procedure matd_locate_
  end interface matd_locate

  public :: matd_data
  interface matd_data
    module procedure matd_data_
  end interface matd_data

  public :: matd_distribution
  interface matd_distribution
    module procedure matd_distribution_
  end interface matd_distribution

  public :: matd_get_nblocks_owned
  interface matd_get_nblocks_owned
    module procedure matd_get_nblocks_owned_
  end interface matd_get_nblocks_owned

  ! =========================

  public :: matd_get_nprocs
  interface matd_get_nprocs
    module procedure matd_get_nprocs_
  end interface matd_get_nprocs

  public :: matd_get_myrank
  interface matd_get_myrank
    module procedure matd_get_myrank_
  end interface matd_get_myrank

  public :: matd_get_win
  interface matd_get_win
    module procedure matd_get_win_
  end interface matd_get_win

  public :: matd_get_dim1
  interface matd_get_dim1
    module procedure matd_get_dim1_
  end interface matd_get_dim1

  public :: matd_get_dim2
  interface matd_get_dim2
    module procedure matd_get_dim2_
  end interface matd_get_dim2

  public :: matd_get_nblocks1
  interface matd_get_nblocks1
    module procedure matd_get_nblocks1_
  end interface matd_get_nblocks1

  public :: matd_get_nblocks2
  interface matd_get_nblocks2
    module procedure matd_get_nblocks2_
  end interface matd_get_nblocks2

  public :: matd_get_nblocks
  interface matd_get_nblocks
    module procedure matd_get_nblocks_
  end interface matd_get_nblocks

  public :: matd_get_procgrid1
  interface matd_get_procgrid1
    module procedure matd_get_procgrid1_
  end interface matd_get_procgrid1

  public :: matd_get_procgrid2
  interface matd_get_procgrid2
    module procedure matd_get_procgrid2_
  end interface matd_get_procgrid2

  public :: matd_get_nprocgrids1
  interface matd_get_nprocgrids1
    module procedure matd_get_nprocgrids1_
  end interface matd_get_nprocgrids1

  public :: matd_get_nprocgrids2
  interface matd_get_nprocgrids2
    module procedure matd_get_nprocgrids2_
  end interface matd_get_nprocgrids2

  public :: matd_get_nprocgrids
  interface matd_get_nprocgrids
    module procedure matd_get_nprocgrids_
  end interface matd_get_nprocgrids

  public :: matd_get_map1
  interface matd_get_map1
    module procedure matd_get_map1_
  end interface matd_get_map1

  public :: matd_get_map2
  interface matd_get_map2
    module procedure matd_get_map2_
  end interface matd_get_map2

  public :: matd_get_is_scalapack
  interface matd_get_is_scalapack
    module procedure matd_get_is_scalapack_
  end interface matd_get_is_scalapack

  ! ----------------------------

  public :: matd_create_reg_gellan
  interface matd_create_reg_gellan
    module procedure matd_create_reg_gellan_
  end interface matd_create_reg_gellan

  public :: matd_create_irreg_gellan
  interface matd_create_irreg_gellan
    module procedure matd_create_irreg_gellan_
  end interface matd_create_irreg_gellan

  public :: matd_create_scalapack_gellan
  interface matd_create_scalapack_gellan
    module procedure matd_create_scalapack_gellan_
  end interface matd_create_scalapack_gellan

  public :: matd_destroy_gellan
  interface matd_destroy_gellan
    module procedure matd_destroy_gellan_
  end interface matd_destroy_gellan


contains

  subroutine matd_create_window_(self, nelems, comm)
!!
!! Subroutine to allocate memory for my block and make the window for RMA. 
!!
    use iso_c_binding
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: nelems, comm
    integer :: sizeoftype, ierr
    integer(kind=mpi_address_kind) :: bytesize
    type(c_ptr) :: baseptr

    call mpi_type_size(MATD_MPI_TYPE, sizeoftype, ierr)
    bytesize = sizeoftype * nelems

!! mpi_alloc_mem for c_ptr does not work in some cases.
!! So, allocate is employed.
!! When you use mpi_alloc_mem, you have to change matd_destroy_
!    call mpi_alloc_mem(bytesize, mpi_info_null, baseptr, ierr)
!    call c_f_pointer(baseptr, self%storage, [nelems])

    allocate(self%storage(nelems))
    call mpi_win_create(self%storage(1), bytesize, sizeoftype, &
                        mpi_info_null, comm, self%win, ierr)
    return
  end subroutine matd_create_window_


  subroutine matd_get_block_index_(self, b, block_index1, block_index2)
!!
!!  Subroutine to get indices of the block of a matrix.
!!
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: b
    integer, intent(out) :: block_index1, block_index2

    block_index1 = mod(b, self%nblocks1) + 1
    block_index2 = b / self%nblocks1 + 1
    return
  end subroutine matd_get_block_index_


  subroutine matd_get_block_index_by_index_(self, index1, index2, block_index1, block_index2)
!!
!!  Subroutine to get the indices of the block from the indices of an element.
!!
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: index1, index2
    integer, intent(out) :: block_index1, block_index2
    integer :: itr

    if (self%map1(self%nblocks1) <= index1) then
      block_index1 = self%nblocks1
    else
      do itr = 1, self%nblocks1
        if (self%map1(itr) - 1 >= index1) then
          block_index1 = itr - 1
          exit
        endif
      enddo
    endif

    if (self%map2(self%nblocks2) <= index2) then
      block_index2 = self%nblocks2
    else
      do itr = 1, self%nblocks2
        if (self%map2(itr) - 1 >= index2) then
          block_index2 = itr - 1
          exit
        endif
      enddo
    endif
    return
  end subroutine matd_get_block_index_by_index_

  subroutine matd_get_block_number_by_index_(self, index1, index2, block_number)
!!
!!  Subroutine to get the block number from the indices of an element.
!!
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: index1, index2
    integer, intent(out) :: block_number
    integer :: block_index1, block_index2

    call matd_get_block_index_by_index_(self, index1, index2, block_index1, block_index2)
    call matd_get_block_number_by_block_index_(self, block_index1, block_index2, block_number)
    return
  end subroutine matd_get_block_number_by_index_


  subroutine matd_get_block_number_by_block_index_(self, block_index1, block_index2, block_number)
!!
!!  Subroutine to get the block number from the indices of a matrix
!!
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: block_index1, block_index2
    integer, intent(out) :: block_number

    block_number = (block_index2 - 1) * self%nblocks1 + (block_index1 - 1)
    return
  end subroutine matd_get_block_number_by_block_index_

  subroutine matd_get_block_range_(self, b, low1, high1, low2, high2)
!!
!!  Subroutine to return the index by calculating the range of a block.
!!
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: b
    integer, intent(out) :: low1, high1, low2, high2
    integer :: block_index1, block_index2

    call matd_get_block_index_(self, b, block_index1, block_index2)

    low1 = self%map1(block_index1)
    low2 = self%map2(block_index2)

    if (block_index1 == self%nblocks1) then
      high1 = self%dim1
    else
      high1 = self%map1(block_index1 + 1) - 1
    endif

    if (block_index2 == self%nblocks2) then
      high2 = self%dim2
    else
      high2 = self%map2(block_index2 + 1) - 1
    endif
    return
  end subroutine matd_get_block_range_

  subroutine matd_get_block_size_(self, b, block_size)
!!
!!  Subroutine to get the number of element in a block
!!
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: b
    integer, intent(out) :: block_size
    integer :: low1, high1, low2, high2

    call matd_get_block_range_(self, b, low1, high1, low2, high2)
    block_size = (high1 - low1 + 1) * (high2 - low2 + 1)
    return
  end subroutine matd_get_block_size_

  subroutine matd_print_info_irreg_blockcyclic_(self)
!!
!!  Subroutine to print the information of irregular block cyclic distribution
!!
    type(matd_matrix), intent(in) :: self
    integer :: i

    call matd_fence_(self)
    if (self%myrank == 3) then
      write (*, '("(dim1=", i0, ", dim2=", i0")")') self%dim1, self%dim2
      write (*, '("(nblocks1=", i0, ", nblocks2=", i0, ", nblocks=", i0")")') &
        self%nblocks1, self%nblocks2, self%nblocks
      write (*, '(a)', advance='no') "map1( "
      do i = 1, self%nblocks1
        write (*, '(i0, " ")', advance='no') self%map1(i)
      end do
      write (*, '(a)') ")"
      write (*, '(a)', advance='no') "map2( "
      do i = 1, self%nblocks2
        write (*, '(i0, " ")', advance='no') self%map2(i)
      end do
      write (*, '(a)') ")"
    end if
    call matd_fence_(self)
    write (*, '("myrank=", i0, " ,storage size=", i0)') &
      self%myrank, size(self%storage)
    call matd_fence_(self)
    return
  end subroutine matd_print_info_irreg_blockcyclic_

  subroutine matd_fence_(self)
!!
!!  Subroutine to fence 
!!
    type(matd_matrix), intent(in) :: self
    integer :: ierr

    call mpi_win_fence(0, self%win, ierr)
    return
  end subroutine matd_fence_

  subroutine matd_create_irreg_blockcyclic_(self, dim1, dim2, map1, map2, comm)
!!
!!  Subroutine to generate the irregular cyclic block distribution.
!!
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, map1(:), map2(:), comm
    integer :: ierr, itr_b, nelems, block_size

    call mpi_comm_size(comm, self%nprocs, ierr)
    call mpi_comm_rank(comm, self%myrank, ierr)
    self%dim1 = dim1
    self%dim2 = dim2
    self%nblocks1 = size(map1)
    self%nblocks2 = size(map2)
    self%nblocks = self%nblocks1 * self%nblocks2
    allocate(self%map1(self%nblocks1), self%map2(self%nblocks2))
    self%map1 = map1
    self%map2 = map2
    self%is_scalapack = .false.

    nelems = 0 ; itr_b = self%myrank
    do while (itr_b < self%nblocks)
      call matd_get_block_size_(self, itr_b, block_size)
      nelems = nelems + block_size
      itr_b = itr_b + self%nprocs
    enddo
    call matd_create_window_(self, nelems, comm)
    return
  end subroutine matd_create_irreg_blockcyclic_

  subroutine matd_get_block_index_of_rank_in_procgrid_(self, r, block_index1, block_index2)
!!
!!  Subroutine to calculate the indices of the block related to the rank in a process grid.
!!
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: r
    integer, intent(out) :: block_index1, block_index2

    block_index1 = mod(r, self%procgrid1) + 1
    block_index2 = r / self%procgrid1 + 1
    return
  end subroutine matd_get_block_index_of_rank_in_procgrid_

  subroutine matd_destroy_(self)
!!
!!  Subroutine to deallocate the memory and destroy the window
!!
    type(matd_matrix), intent(out) :: self
    integer :: ierr

    call mpi_win_free(self%win, ierr)
    deallocate(self%map1, self%map2)

!!  When you used mpi_alloc_mem, activate the following line and
!!  comment out "deallocate"
!    call mpi_free_mem(self%storage, ierr)
    deallocate(self%storage)

    return
  end subroutine matd_destroy_


  !
  ! プロセッサグリッドの範囲を求めブロックに関するインデックスで返すサブルーチン
  !
  subroutine matd_get_procgrid_range_(self, pg, blow1, bhigh1, blow2, bhigh2)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: pg
    integer, intent(out) :: blow1, bhigh1, blow2, bhigh2

    blow1 = mod(pg, self%nprocgrids1) * self%procgrid1 + 1
    blow2 = pg / self%nprocgrids1 * self%procgrid2 + 1

    if ((blow1 + self%procgrid1 - 1) > (self%nblocks1 - 1)) then
      bhigh1 = self%nblocks1 - 1
    else
      bhigh1 = blow1 + self%procgrid1 - 1
    endif

    if ((blow2 + self%procgrid2 - 1) > (self%nblocks2 - 1)) then
      bhigh2 = self%nblocks2 - 1
    else
      bhigh2 = blow2 + self%procgrid2 - 1
    endif
    return
  end subroutine matd_get_procgrid_range_


  !
  ! 指定されたプロセッサグリッド内の指定されたランクが担当するブロックの
  ! ブロック番号を求めるサブルーチン
  ! ブロックが存在しない場合は-1を返すこととする
  !
  subroutine matd_get_block_number_of_rank_in_procgrid_(self, pg, r, block_number)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: pg, r
    integer, intent(out) :: block_number
    integer :: block_index1_disp, block_index2_disp, blow1, bhigh1, blow2, bhigh2

    call matd_get_block_index_of_rank_in_procgrid_(self, r, block_index1_disp, block_index2_disp)
    call matd_get_procgrid_range_(self, pg, blow1, bhigh1, blow2, bhigh2)
    if (((blow1 + block_index1_disp - 1) > self%nblocks1) .or. &
      ((blow2 + block_index2_disp - 1) > self%nblocks2)) then
      block_number = -1
    else
      call matd_get_block_number_by_block_index_(self, blow1 + block_index1_disp - 1, &
        blow2 + block_index2_disp - 1, block_number)
    endif
    return
  end subroutine matd_get_block_number_of_rank_in_procgrid_


  !
  ! プロセッサグリッド付きイレギュラーブロックサイクリック分散行列
  ! の情報を出力するサブルーチン
  !
  subroutine matd_print_info_irreg_scalapack_(self)
    type(matd_matrix), intent(in) :: self
    integer :: i

    call matd_fence_(self)
    if (self%myrank == 3) then
      write (*, '("(dim1=", i0, ", dim2=", i0")")') self%dim1, self%dim2
      write (*, '("(nblocks1=", i0, ", nblocks2=", i0, ", nblocks=", i0")")') &
        self%nblocks1, self%nblocks2, self%nblocks
      write (*, '(a)', advance='no') "map1( "
      do i = 1, self%nblocks1
        write (*, '(i0, " ")', advance='no') self%map1(i)
      end do
      write (*, '(a)') ")"
      write (*, '(a)', advance='no') "map2( "
      do i = 1, self%nblocks2
        write (*, '(i0, " ")', advance='no') self%map2(i)
      end do
      write (*, '(a)') ")"
      write (*, '("procgrid1=", i0, ", procgrid2=", i0)') self%procgrid1, self%procgrid2
      write (*, '("nprocgrids1=", i0, ", nprocgrids2=", i0, ", nprocgrids=", i0)') &
        self%nprocgrids1, self%nprocgrids2, self%nprocgrids
    end if
    call matd_fence_(self)
    write (*, '("myrank=", i0, " ,storage size=", i0)') &
      self%myrank, size(self%storage)
    call matd_fence_(self)
    return
  end subroutine matd_print_info_irreg_scalapack_


  !
  ! プロセッサグリッド付きイレギュラーブロックサイクリック分散行列
  ! を生成するサブルーチン
  !
  subroutine matd_create_irreg_scalapack_(self, dim1, dim2, map1, map2, procgrid1, procgrid2, comm)
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, map1(:), map2(:), procgrid1, procgrid2, comm
    integer :: ierr, nelems, blow1, bhigh1, blow2, bhigh2, block_index1_disp, block_index2_disp, &
      block_number, block_size, itr_pg

    call mpi_comm_size(comm, self%nprocs, ierr)
    call mpi_comm_rank(comm, self%myrank, ierr)
    self%dim1 = dim1
    self%dim2 = dim2
    self%nblocks1 = size(map1)
    self%nblocks2 = size(map2)
    self%nblocks = self%nblocks1 * self%nblocks2
    allocate(self%map1(self%nblocks1), self%map2(self%nblocks2))
    self%map1 = map1
    self%map2 = map2
    self%procgrid1 = procgrid1
    self%procgrid2 = procgrid2
    self%is_scalapack = .true.

    if (mod(self%nblocks1, self%procgrid1) == 0) then
      self%nprocgrids1 = self%nblocks1 / self%procgrid1
    else
      self%nprocgrids1 = self%nblocks1 / self%procgrid1 + 1
    endif
    if (mod(self%nblocks2, self%procgrid2) == 0) then
      self%nprocgrids2 = self%nblocks2 / self%procgrid2
    else
      self%nprocgrids2 = self%nblocks2 / self%procgrid2 + 1
    endif
    self%nprocgrids = self%nprocgrids1 * self%nprocgrids2

    nelems = 0
    do itr_pg = 0, self%nprocgrids - 1
      call matd_get_block_number_of_rank_in_procgrid_(self, itr_pg, self%myrank, block_number)
      if (block_number >= 0) then
        call matd_get_block_size_(self, block_number, block_size)
        nelems = nelems + block_size
      endif
    enddo
    call matd_create_window_(self, nelems, comm)
    return
  end subroutine matd_create_irreg_scalapack_


  !
  ! イレギュラーブロックサイクリック分散行列のLocate操作を行うサブルーチン
  !
  subroutine matd_locate_irreg_blockcyclic_(self, index1, index2, owner_rank)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: index1, index2
    integer, intent(out) :: owner_rank
    integer :: block_number

    call matd_get_block_number_by_index_(self, index1, index2, block_number)
    owner_rank = mod(block_number, self%nprocs)
    return
  end subroutine matd_locate_irreg_blockcyclic_


  !
  ! イレギュラーブロックサイクリック分散行列へのPut/Get操作を行うサブルーチン
  !
  subroutine matd_put_get_common_irreg_blockcyclic_(self, low1, high1, low2, high2, buf, is_put)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: low1, high1, low2, high2
    MATD_ELEMTYPE, intent(inout) :: buf(:)
    logical, intent(in) :: is_put
    integer :: cur_index2, cur_index1, cur_low1, cur_high1, &
      cur_low2, cur_high2, cur_block_number, cur_rank, cur_fst1, cur_lst1, nelems, &
      cur_block_index1, cur_block_index2, ierr, itr_b, sub_block_index1, sub_block_index2, &
      sub_low1, sub_high1, sub_low2, sub_high2, sub_block_size
    integer(kind=mpi_address_kind) :: disp
    integer(8) :: buf_index

    cur_index2 = low2 ;  buf_index = 1
    do while (cur_index2 <= high2)
      cur_index1 = low1
      do while (cur_index1 <= high1)
        call matd_get_block_index_by_index_(self, cur_index1, cur_index2, &
          cur_block_index1, cur_block_index2)
        call matd_get_block_number_by_block_index_(self, cur_block_index1, cur_block_index2, &
          cur_block_number)
        call matd_locate_irreg_blockcyclic_(self, cur_index1, cur_index2, cur_rank)
        call matd_get_block_range_(self, cur_block_number, cur_low1, cur_high1, cur_low2, cur_high2)

        cur_fst1 = max(low1, cur_low1) ; cur_lst1 = min(high1, cur_high1)
        nelems = cur_lst1 - cur_fst1 + 1

        ! ブロックを辿り、displacementを計算する
        disp = 0 ; itr_b = cur_rank
        do while (itr_b < self%nblocks)
          call matd_get_block_index_(self, itr_b, sub_block_index1, sub_block_index2)
          if (sub_block_index2 > cur_block_index2) exit

          if (cur_block_index2 /= sub_block_index2) then
            call matd_get_block_size_(self, itr_b, sub_block_size)
            disp = disp + sub_block_size
          else
            call matd_get_block_range_(self, itr_b, sub_low1, sub_high1, sub_low2, sub_high2)
            if (sub_high1 < cur_index1) then
              disp = disp + (cur_index2 - sub_low2 + 1) * (sub_high1 - sub_low1 + 1)
            else
              disp = disp + (cur_index2 - sub_low2) * (sub_high1 - sub_low1 + 1) + &
                max(0, cur_index1 - sub_low1)
            endif
          endif
          itr_b = itr_b + self%nprocs
        enddo

        call mpi_win_lock(mpi_lock_exclusive, cur_rank, mpi_mode_nocheck, self%win, ierr)
        if (is_put) then
          call mpi_put(buf(buf_index), nelems, MATD_MPI_TYPE, &
            cur_rank, disp, nelems, MATD_MPI_TYPE, self%win, ierr)
        else
          call mpi_get(buf(buf_index), nelems, MATD_MPI_TYPE, &
            cur_rank, disp, nelems, MATD_MPI_TYPE, self%win, ierr)
        endif
        call mpi_win_unlock(cur_rank, self%win, ierr)

        buf_index = buf_index + nelems
        cur_index1 = cur_lst1 + 1
      enddo
      cur_index2 = cur_index2 + 1
    enddo
    return
  end subroutine matd_put_get_common_irreg_blockcyclic_


  !
  ! イレギュラーブロックサイクリック分散行列へのPut操作を行うサブルーチン
  !
  subroutine matd_put_irreg_blockcyclic_(self, low1, high1, low2, high2, buf)
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: low1, high1, low2, high2
    MATD_ELEMTYPE, intent(inout) :: buf(:)

    call matd_put_get_common_irreg_blockcyclic_(self, low1, high1, low2, high2, buf, .true.)
    return
  end subroutine matd_put_irreg_blockcyclic_


  !
  ! イレギュラーブロックサイクリック分散行列へのGet操作を行うサブルーチン
  !
  subroutine matd_get_irreg_blockcyclic_(self, low1, high1, low2, high2, buf)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: low1, high1, low2, high2
    MATD_ELEMTYPE, intent(inout) :: buf(:)

    call matd_put_get_common_irreg_blockcyclic_(self, low1, high1, low2, high2, buf, .false.)
    return
  end subroutine matd_get_irreg_blockcyclic_


  !
  ! プロセッサグリッド付きイレギュラーブロックサイクリック分散行列
  ! へのLocate操作を行うサブルーチン
  !
  subroutine matd_locate_irreg_scalapack_(self, index1, index2, owner_rank)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: index1, index2
    integer, intent(out) :: owner_rank
    integer :: block_index1, block_index2

    call matd_get_block_index_by_index_(self, index1, index2, block_index1, block_index2)
    owner_rank = (mod(block_index2 - 1, self%procgrid2) * self%procgrid1) + &
      mod(block_index1 - 1, self%procgrid1)
    return
  end subroutine matd_locate_irreg_scalapack_


  !
  ! プロセッサグリッド付きイレギュラーブロックサイクリック分散行列
  ! へのPut/Get操作を行うサブルーチン
  !
  subroutine matd_put_get_common_irreg_scalapack_(self, low1, high1, low2, high2, buf, is_put)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: low1, high1, low2, high2
    MATD_ELEMTYPE, intent(inout) :: buf(:)
    logical, intent(in) :: is_put
    integer :: cur_index2, cur_index1, cur_block_index1, cur_block_index2, &
      cur_block_number, cur_rank, cur_fst1, cur_lst1, cur_low1, cur_high1, cur_low2, cur_high2, &
      itr_pg, sub_block_number, sub_block_index1, sub_block_index2, sub_block_size, &
      sub_low1, sub_high1, sub_low2, sub_high2, nelems, ierr
    integer(kind=mpi_address_kind) :: disp
    integer(8) :: buf_index ! これはかなり大きくなる可能性があるの８バイト

    cur_index2 = low2 ; buf_index = 1
    do while (cur_index2 <= high2)
      cur_index1 = low1
      do while (cur_index1 <= high1)
        call matd_get_block_index_by_index_(self, cur_index1, cur_index2, &
          cur_block_index1, cur_block_index2)
        call matd_get_block_number_by_block_index_(self, cur_block_index1, cur_block_index2, &
          cur_block_number)
        call matd_locate_irreg_scalapack_(self, cur_index1, cur_index2, cur_rank)
        call matd_get_block_range_(self, cur_block_number, cur_low1, cur_high1, &
          cur_low2, cur_high2)

        cur_fst1 = max(low1, cur_low1) ; cur_lst1 = min(high1, cur_high1)
        nelems = cur_lst1 - cur_fst1 + 1

        disp = 0
        do itr_pg = 0, self%nprocgrids
          call matd_get_block_number_of_rank_in_procgrid_(self, itr_pg, cur_rank, sub_block_number)
          if (sub_block_number >= 0) then
            call matd_get_block_index_(self, sub_block_number, sub_block_index1, sub_block_index2)
            if (sub_block_index2 > cur_block_index2) exit

            if (cur_block_index2 /= sub_block_index2) then
              call matd_get_block_size_(self, sub_block_number, sub_block_size)
              disp = disp + sub_block_size
            else
              call matd_get_block_range_(self, sub_block_number, &
                sub_low1, sub_high1, sub_low2, sub_high2)
              if (sub_high1 < cur_index1) then
                disp = disp + (cur_index2 - sub_low2 + 1) * (sub_high1 - sub_low1 + 1)
              else
                disp = disp + (cur_index2 - sub_low2) * (sub_high1 - sub_low1 + 1) + &
                  max(0, cur_index1 - sub_low1)
              endif
            endif
          endif
        enddo

        call mpi_win_lock(mpi_lock_exclusive, cur_rank, mpi_mode_nocheck, self%win, ierr)
        if (is_put) then
          call mpi_put(buf(buf_index), nelems, MATD_MPI_TYPE, &
            cur_rank, disp, nelems, MATD_MPI_TYPE, self%win, ierr)
        else
          call mpi_get(buf(buf_index), nelems, MATD_MPI_TYPE, &
            cur_rank, disp, nelems, MATD_MPI_TYPE, self%win, ierr)
        endif
        call mpi_win_unlock(cur_rank, self%win, ierr)

        buf_index = buf_index + nelems
        cur_index1 = cur_lst1 + 1
      enddo
      cur_index2 = cur_index2 + 1
    enddo
    return
  end subroutine matd_put_get_common_irreg_scalapack_


  !
  ! プロセッサグリッド付きイレギュラーブロックサイクリック分散行列
  ! へのPut操作を行うサブルーチン
  !
  subroutine matd_put_irreg_scalapack_(self, low1, high1, low2, high2, buf)
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: low1, high1, low2, high2
    MATD_ELEMTYPE, intent(inout) :: buf(:)

    call matd_put_get_common_irreg_scalapack_(self, low1, high1, low2, high2, buf, .true.)
    return
  end subroutine matd_put_irreg_scalapack_


  !
  ! プロセッサグリッド付きイレギュラーブロックサイクリック分散行列
  ! へのGet操作を行うサブルーチン
  !
  subroutine matd_get_irreg_scalapack_(self, low1, high1, low2, high2, buf)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: low1, high1, low2, high2
    MATD_ELEMTYPE, intent(inout) :: buf(:)

    call matd_put_get_common_irreg_scalapack_(self, low1, high1, low2, high2, buf, .false.)
    return
  end subroutine matd_get_irreg_scalapack_


  !
  ! ブロック(サイクリック)分散行列を生成するサブルーチン
  !
  subroutine matd_create_reg_(self, dim1, dim2, block_size1, block_size2, comm, is_blockcyclic)
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, block_size1, block_size2, comm
    logical, intent(in), optional :: is_blockcyclic
    integer, allocatable :: map1(:), map2(:)
    integer :: nblocks1, nblocks2, itr, nprocs, ierr
    logical :: is_bc

    if (mod(dim1, block_size1) == 0) then
      nblocks1 = dim1 / block_size1
    else
      nblocks1 = dim1 / block_size1 + 1
    endif

    if (mod(dim2, block_size2) == 0) then
      nblocks2 = dim2 / block_size2
    else
      nblocks2 = dim2 / block_size2 + 1
    endif

    if (present(is_blockcyclic)) then
      is_bc = is_blockcyclic
    else
      is_bc = .false.
    endif
    if (.not. is_bc) then
      call mpi_comm_size(comm, nprocs, ierr)
      if (nprocs < nblocks1 * nblocks2) then
        write (*, *) "ERROR : Creating a matrix of regular block distribution"
        stop
      endif
    endif

    allocate(map1(nblocks1), map2(nblocks2))
    do itr = 1, nblocks1
      map1(itr) = (itr - 1) * block_size1 + 1
    enddo
    do itr = 1, nblocks2
      map2(itr) = (itr - 1) * block_size2 + 1
    enddo

    call matd_create_irreg_blockcyclic_(self, dim1, dim2, map1, map2, comm)
    deallocate(map1, map2)
    return
  end subroutine matd_create_reg_


  !
  ! イレギュラーブロック(サイクリック)分散行列を生成するサブルーチン
  !
  subroutine matd_create_irreg_(self, dim1, dim2, map1, map2, comm, is_blockcyclic)
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, map1(:), map2(:), comm
    logical, intent(in), optional :: is_blockcyclic
    integer :: ierr, nprocs
    logical :: is_bc

    if (present(is_blockcyclic)) then
      is_bc = is_blockcyclic
    else
      is_bc = .false.
    endif
    if (.not. is_bc) then
      call mpi_comm_size(comm, nprocs, ierr)
      if (nprocs < size(map1) * size(map2)) then
        write (*, *) "ERROR : Creating a matrix of irregular block distribution"
        stop
      endif
    endif

    call matd_create_irreg_blockcyclic_(self, dim1, dim2, map1, map2, comm)
    return
  end subroutine matd_create_irreg_


  !
  ! プロセッサグリッド付きレギュラーブロックサイクリック分散行列
  ! を生成するサブルーチン
  !
  subroutine matd_create_scalapack_(self, dim1, dim2, block_size1, block_size2, &
      procgrid1, procgrid2, comm)
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, block_size1, block_size2, procgrid1, procgrid2, comm
    integer, allocatable :: map1(:), map2(:)
    integer :: nblocks1, nblocks2, itr

    if (mod(dim1, block_size1) == 0) then
      nblocks1 = dim1 / block_size1
    else
      nblocks1 = dim1 / block_size1 + 1
    endif

    if (mod(dim2, block_size2) == 0) then
      nblocks2 = dim2 / block_size2
    else
      nblocks2 = dim2 / block_size2 + 1
    endif

    allocate(map1(nblocks1), map2(nblocks2))
    do itr = 1, nblocks1
      map1(itr) = (itr - 1) * block_size1 + 1
    enddo
    do itr = 1, nblocks2
      map2(itr) = (itr - 1) * block_size2 + 1
    enddo

    call matd_create_irreg_scalapack_(self, dim1, dim2, map1, map2, procgrid1, procgrid2, comm)
    deallocate(map1, map2)
    return
  end subroutine matd_create_scalapack_


  !
  ! 任意の分散行列へPut操作を行うサブルーチン
  !
  subroutine matd_put_(self, low1, high1, low2, high2, buf)
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: low1, high1, low2, high2
    MATD_ELEMTYPE, intent(inout) :: buf(:)

    if (self%is_scalapack) then
      call matd_put_irreg_scalapack_(self, low1, high1, low2, high2, buf)
    else
      call matd_put_irreg_blockcyclic_(self, low1, high1, low2, high2, buf)
    endif
    return
  end subroutine matd_put_


  !
  ! 任意の分散行列へGet操作を行うサブルーチン
  !
  subroutine matd_get_(self, low1, high1, low2, high2, buf)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: low1, high1, low2, high2
    MATD_ELEMTYPE, intent(inout) :: buf(:)

    if (self%is_scalapack) then
      call matd_get_irreg_scalapack_(self, low1, high1, low2, high2, buf)
    else
      call matd_get_irreg_blockcyclic_(self, low1, high1, low2, high2, buf)
    endif
    return
  end subroutine matd_get_


  !
  ! 任意の分散行列へのLocate操作を行うサブルーチン
  !
  subroutine matd_locate_(self, index1, index2, owner_rank)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: index1, index2
    integer, intent(out) :: owner_rank

    if (self%is_scalapack) then
      call matd_locate_irreg_scalapack_(self, index1, index2, owner_rank)
    else
      call matd_locate_irreg_blockcyclic_(self, index1, index2, owner_rank)
    endif
    return
  end subroutine matd_locate_


  !
  ! 行列全体を同じサイズの行列にコピーするサブルーチン
  !
  subroutine matd_copy_(self, source)
    type(matd_matrix), intent(out) :: self
    type(matd_matrix), intent(in) :: source
    integer :: itr_b, low1, high1, low2, high2
    MATD_ELEMTYPE, allocatable :: buf(:)

    itr_b = self%myrank

    do while(itr_b < self%nblocks)
      call matd_get_block_range_(self, itr_b, low1, high1, low2, high2)
      allocate(buf((high1 - low1 + 1) * (high2 - low2 + 1)))
      call matd_get_(source, low1, high1, low2, high2, buf)
      call matd_put_(self, low1, high1, low2, high2, buf)
      deallocate(buf)

      itr_b = itr_b + self%nprocs
    enddo

    return
  end subroutine matd_copy_


  !
  ! 行列全体を出力するサブルーチン(最悪の効率 デバッグ用)
  !
  subroutine matd_print_(self)
    type(matd_matrix), intent(in) :: self
    MATD_ELEMTYPE :: buf(1)
    integer :: itr1, itr2
    if (self%myrank == 0) then
      write (*, *) "["
      do itr1 = 1, self%dim1
        do itr2 = 1, self%dim2
          call matd_get_(self, itr1, itr1, itr2, itr2, buf)
          ! それっぽく動いていることさえ分かればいいので出力の精度はかなり低く設定
          write (*, '(f10.2, " ")', advance='no') buf(1)
        enddo
        write (*, *) " "
      enddo
      write (*, *) "]"
    endif
    return
  end subroutine matd_print_


  !
  ! 行列の分散状況を出力するサブルーチン
  !
  subroutine matd_print_info_(self)
    type(matd_matrix), intent(in) :: self

    if (self%is_scalapack) then
      call matd_print_info_irreg_scalapack_(self)
    else
      call matd_print_info_irreg_blockcyclic_(self)
    endif
    return
  end subroutine matd_print_info_


  !
  ! ストレージ領域のアドレスを返すサブルーチン
  !
  subroutine matd_data_(self, ptr)
    type(matd_matrix), intent(in) :: self
    MATD_ELEMTYPE, pointer :: ptr(:)

    ptr => self%storage
    return
  end subroutine matd_data_


  !
  ! イレギュラーブロックサイクリック分散行列のディストリビューション処理
  !
  subroutine matd_distribution_irreg_blockcyclic_(self, r, lows1, highs1, lows2, highs2)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: r
    integer, intent(out) :: lows1(:), highs1(:), lows2(:), highs2(:)
    integer :: itr_b, i

    i = 1;  itr_b = r
    do while (itr_b < self%nblocks)
      call matd_get_block_range_(self, itr_b, lows1(i), highs1(i), lows2(i), highs2(i))
      i = i + 1
      itr_b = itr_b + self%nprocs
    enddo
    return
  end subroutine matd_distribution_irreg_blockcyclic_


  !
  ! ScaLAPACK型イレギュラーブロックサイクリック分散行列の
  ! ディストリビューション処理
  !
  subroutine matd_distribution_irreg_scalapack_(self, r, lows1, highs1, lows2, highs2)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: r
    integer, intent(out) :: lows1(:), highs1(:), lows2(:), highs2(:)
    integer :: itr_pg, i, b

    i = 1
    do itr_pg = 0, self%nprocgrids - 1
      call matd_get_block_number_of_rank_in_procgrid_(self, itr_pg, r, b)
      if (b >= 0) then
        call matd_get_block_range_(self, b, lows1(i), highs1(i), lows2(i), highs2(i))
        i = i + 1
      endif
    enddo
    return
  end subroutine matd_distribution_irreg_scalapack_


  !
  ! ディストリビューション処理
  !
  subroutine matd_distribution_(self, r, lows1, highs1, lows2, highs2)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: r
    integer, intent(out) :: lows1(:), highs1(:), lows2(:), highs2(:)

    if (self%is_scalapack) then
      call matd_distribution_irreg_scalapack_(self, r, lows1, highs1, lows2, highs2)
    else
      call matd_distribution_irreg_blockcyclic_(self, r, lows1, highs1, lows2, highs2)
    endif
    return
  end subroutine matd_distribution_


  !
  ! イレギュラーブロックサイクリック分散行列の
  ! あるランクのプロセスが保持しているブロックの数を返す
  !
  subroutine matd_get_nblocks_owned_irreg_blockcyclic_(self, r, nblocks)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: r
    integer, intent(out) :: nblocks
    integer :: itr_b

    nblocks = 0;  itr_b = r
    do while (itr_b < self%nblocks)
      nblocks = nblocks + 1
      itr_b = itr_b + self%nprocs
    enddo
    return
  end subroutine matd_get_nblocks_owned_irreg_blockcyclic_


  !
  ! ScaLAPAKC型イレギュラーブロックサイクリック分散行列の
  ! あるランクのプロセスが保持しているブロックの数を返す
  !
  subroutine matd_get_nblocks_owned_irreg_scalapack_(self, r, nblocks)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: r
    integer, intent(out) :: nblocks
    integer :: itr_pg, b

    nblocks = 0;
    do itr_pg = 0, self%nprocgrids - 1
      call matd_get_block_number_of_rank_in_procgrid_(self, itr_pg, r, b)
      if (b >= 0) then
        nblocks = nblocks + 1
      endif
    enddo
    return
  end subroutine matd_get_nblocks_owned_irreg_scalapack_


  !
  ! あるランクのプロセスが保持しているブロックの数を返す
  !
  subroutine matd_get_nblocks_owned_(self, r, nblocks)
    type(matd_matrix), intent(in) :: self
    integer, intent(in) :: r
    integer, intent(out) :: nblocks

    if (self%is_scalapack) then
      call matd_get_nblocks_owned_irreg_scalapack_(self, r, nblocks)
    else
      call matd_get_nblocks_owned_irreg_blockcyclic_(self, r, nblocks)
    endif
    return
  end subroutine matd_get_nblocks_owned_


  !
  ! メンバ変数のゲッタ
  !
  subroutine matd_get_nprocs_(self, nprocs)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: nprocs

    nprocs = self%nprocs
    return
  end subroutine matd_get_nprocs_

  subroutine matd_get_myrank_(self, myrank)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: myrank

    myrank = self%myrank
    return
  end subroutine matd_get_myrank_

  subroutine matd_get_win_(self, win)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: win

    win = self%win
    return
  end subroutine matd_get_win_

  subroutine matd_get_dim1_(self, dim1)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: dim1

    dim1 = self%dim1
    return
  end subroutine matd_get_dim1_

  subroutine matd_get_dim2_(self, dim2)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: dim2

    dim2 = self%dim2
    return
  end subroutine matd_get_dim2_

  subroutine matd_get_nblocks1_(self, nblocks1)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: nblocks1

    nblocks1 = self%nblocks1
    return
  end subroutine matd_get_nblocks1_

  subroutine matd_get_nblocks2_(self, nblocks2)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: nblocks2

    nblocks2 = self%nblocks2
    return
  end subroutine matd_get_nblocks2_

  subroutine matd_get_nblocks_(self, nblocks)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: nblocks

    nblocks = self%nblocks
    return
  end subroutine matd_get_nblocks_

  subroutine matd_get_procgrid1_(self, procgrid1)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: procgrid1

    procgrid1 = self%procgrid1
    return
  end subroutine matd_get_procgrid1_

  subroutine matd_get_procgrid2_(self, procgrid2)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: procgrid2

    procgrid2 = self%procgrid2
    return
  end subroutine matd_get_procgrid2_

  subroutine matd_get_nprocgrids1_(self, nprocgrids1)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: nprocgrids1

    nprocgrids1 = self%nprocgrids1
    return
  end subroutine matd_get_nprocgrids1_

  subroutine matd_get_nprocgrids2_(self, nprocgrids2)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: nprocgrids2

    nprocgrids2 = self%nprocgrids2
    return
  end subroutine matd_get_nprocgrids2_

  subroutine matd_get_nprocgrids_(self, nprocgrids)
    type(matd_matrix), intent(in) :: self
    integer, intent(out) :: nprocgrids

    nprocgrids = self%nprocgrids
    return
  end subroutine matd_get_nprocgrids_

  subroutine matd_get_is_scalapack_(self, is_scalapack)
    type(matd_matrix), intent(in) :: self
    logical, intent(out) :: is_scalapack

    is_scalapack = self%is_scalapack
    return
  end subroutine matd_get_is_scalapack_

  subroutine matd_get_map1_(self, map1)
    type(matd_matrix), intent(in) :: self
    MATD_ELEMTYPE :: map1(:)

    map1 = self%map1
    return
  end subroutine matd_get_map1_

  subroutine matd_get_map2_(self, map2)
    type(matd_matrix), intent(in) :: self
    MATD_ELEMTYPE :: map2(:)

    map2 = self%map2
    return
  end subroutine matd_get_map2_


  ! ========================================================
  ! GELLAN用のサブルーチン群
  ! ========================================================

  !
  ! [GELLAN用]
  ! イレギュラーブロックサイクリック分散行列を生成するサブルーチン
  !
  subroutine matd_create_irreg_blockcyclic_gellan_(&
      self, dim1, dim2, map1, map2, comm, &
      mem_pool, & ! X自体を渡す
      ibegin,   & ! Xの何番目のインデックスから使うかを渡す
      iend,     & ! どこまで使ったかを返す（使った次のインデックス）。
      map_flag  & ! map用の領域をアロケートする（true）かそのまま使うか（false）
      )
    use iso_c_binding

    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, comm
    integer, target, intent(in) :: map1(:), map2(:)
    double precision, target, intent(in) :: mem_pool(:)
    integer, intent(in) :: ibegin
    integer, intent(out) :: iend
    logical, intent(in) :: map_flag

    integer :: ierr, itr_b, nelems, block_size
    type(c_ptr) :: cptr
    integer :: pool_index, sizeoftype
    integer(kind=mpi_address_kind) :: bytesize

    pool_index = ibegin

    call mpi_comm_size(comm, self%nprocs, ierr)
    call mpi_comm_rank(comm, self%myrank, ierr)
    self%dim1 = dim1
    self%dim2 = dim2
    self%nblocks1 = size(map1)
    self%nblocks2 = size(map2)
    self%nblocks = self%nblocks1 * self%nblocks2
    self%is_scalapack = .false.

    ! マップ情報の登録
    if (map_flag) then ! Xから取る場合
      cptr = c_loc(mem_pool(pool_index))
      call c_f_pointer(cptr, self%map1, shape=[self%nblocks1])

      ! integerは4バイトなので、8バイトのXの領域を取っていく際、
      ! 要素数が奇数の場合+1しなければならない。
      if (mod(self%nblocks1, 2) == 0) then
        pool_index = pool_index + (self%nblocks1 / 2)
      else
        pool_index = pool_index + (self%nblocks1 / 2) + 1
      endif

      cptr = c_loc(mem_pool(pool_index))
      call c_f_pointer(cptr, self%map2, shape=[self%nblocks2])

      if (mod(self%nblocks2, 2) == 0) then
        pool_index = pool_index + (self%nblocks1 / 2)
      else
        pool_index = pool_index + (self%nblocks1 / 2) + 1
      endif

      self%map1 = map1
      self%map2 = map2

    else
      self%map1 => map1
      self%map2 => map2
    endif

    nelems = 0 ; itr_b = self%myrank
    do while (itr_b < self%nblocks)
      call matd_get_block_size_(self, itr_b, block_size)
      nelems = nelems + block_size
      itr_b = itr_b + self%nprocs
    enddo

    ! 領域の割り当て
    call mpi_type_size(MATD_MPI_TYPE, sizeoftype, ierr)
    cptr = c_loc(mem_pool(pool_index))
    call c_f_pointer(cptr, self%storage, shape=[nelems])
    bytesize = sizeoftype * nelems

    if (mod(bytesize, 8) == 0) then
      pool_index = pool_index + (bytesize / 8)
    else
      pool_index = pool_index + (bytesize / 8) + 1
    endif

    call mpi_win_create(self%storage(1), bytesize, sizeoftype, &
      mpi_info_null, comm, self%win, ierr)
    iend = pool_index
    return
  end subroutine matd_create_irreg_blockcyclic_gellan_


  !
  ! [GELLAN用]
  ! ScaLAPACK型イレギュラーブロックサイクリック分散行列を生成するサブルーチン
  !
  subroutine matd_create_irreg_scalapack_gellan_(&
      self, dim1, dim2, map1, map2, procgrid1, procgrid2,&
      comm, mem_pool, ibegin, iend)
    use iso_c_binding
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, procgrid1, procgrid2, comm
    integer, target, intent(in) :: map1(:), map2(:)
    double precision, target, intent(in) :: mem_pool(:)
    integer, intent(in) :: ibegin
    integer, intent(out) :: iend
    integer :: ierr, nelems, blow1, bhigh1, blow2, bhigh2, block_index1_disp, block_index2_disp, &
      block_number, block_size, itr_pg
    type(c_ptr) :: cptr
    integer :: pool_index, sizeoftype
    integer(kind=mpi_address_kind) :: bytesize

    call mpi_comm_size(comm, self%nprocs, ierr)
    call mpi_comm_rank(comm, self%myrank, ierr)
    self%dim1 = dim1
    self%dim2 = dim2
    self%nblocks1 = size(map1)
    self%nblocks2 = size(map2)
    self%nblocks = self%nblocks1 * self%nblocks2

    self%map1 => map1
    self%map2 => map2
    self%procgrid1 = procgrid1
    self%procgrid2 = procgrid2
    self%is_scalapack = .true.

    if (mod(self%nblocks1, self%procgrid1) == 0) then
      self%nprocgrids1 = self%nblocks1 / self%procgrid1
    else
      self%nprocgrids1 = self%nblocks1 / self%procgrid1 + 1
    endif
    if (mod(self%nblocks2, self%procgrid2) == 0) then
      self%nprocgrids2 = self%nblocks2 / self%procgrid2
    else
      self%nprocgrids2 = self%nblocks2 / self%procgrid2 + 1
    endif
    self%nprocgrids = self%nprocgrids1 * self%nprocgrids2

    nelems = 0
    do itr_pg = 0, self%nprocgrids - 1
      call matd_get_block_number_of_rank_in_procgrid_(self, itr_pg, self%myrank, block_number)
      if (block_number >= 0) then
        call matd_get_block_size_(self, block_number, block_size)
        nelems = nelems + block_size
      endif
    enddo

    call mpi_type_size(MATD_MPI_TYPE, sizeoftype, ierr)
    cptr = c_loc(mem_pool(pool_index))
    call c_f_pointer(cptr, self%storage, shape=[nelems])
    bytesize = sizeoftype * nelems

    pool_index = ibegin
    if (mod(bytesize, 8) == 0) then
      pool_index = pool_index + (bytesize / 8)
    else
      pool_index = pool_index + (bytesize / 8) + 1
    endif

    call mpi_win_create(self%storage(1), bytesize, sizeoftype, &
      mpi_info_null, comm, self%win, ierr)
    iend = pool_index
    return
  end subroutine matd_create_irreg_scalapack_gellan_


  !
  ! [GELLAN用]
  ! ブロック(サイクリック)分散行列を生成するサブルーチン
  !
  subroutine matd_create_reg_gellan_(self, dim1, dim2, block_size1, block_size2, &
      comm, mem_pool, ibegin, iend, is_blockcyclic)
    use iso_c_binding
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, block_size1, block_size2, comm
    double precision, target, intent(in) :: mem_pool(:)
    integer, intent(in) :: ibegin
    integer, intent(out) :: iend
    logical, intent(in), optional :: is_blockcyclic
    integer, pointer :: map1(:), map2(:)
    integer :: nblocks1, nblocks2, itr, nprocs, ierr
    logical :: is_bc
    type(c_ptr) :: cptr
    integer :: pool_index

    if (mod(dim1, block_size1) == 0) then
      nblocks1 = dim1 / block_size1
    else
      nblocks1 = dim1 / block_size1 + 1
    endif

    if (mod(dim2, block_size2) == 0) then
      nblocks2 = dim2 / block_size2
    else
      nblocks2 = dim2 / block_size2 + 1
    endif

    if (present(is_blockcyclic)) then
      is_bc = is_blockcyclic
    else
      is_bc = .false.
    endif
    if (.not. is_bc) then
      call mpi_comm_size(comm, nprocs, ierr)
      if (nprocs < nblocks1 * nblocks2) then
        write (*, *) "ERROR : Creating a matrix of regular block distribution"
        stop
      endif
    endif

    pool_index = ibegin
    ! マップ情報の領域もプールから取る
    cptr = c_loc(mem_pool(pool_index))
    call c_f_pointer(cptr, map1, shape=[nblocks1])
    ! integerは4バイトなのでインデックスは要素数の半分加算
    if (mod(nblocks1, 2) == 0) then
      pool_index = pool_index + (nblocks1 / 2)
    else
      pool_index = pool_index + (nblocks1 / 2) + 1
    endif

    cptr = c_loc(mem_pool(pool_index))
    call c_f_pointer(cptr, map2, shape=[nblocks2])
    ! integerは4バイトなのでインデックスは要素数の半分加算
    if (mod(nblocks2, 2) == 0) then
      pool_index = pool_index + (nblocks2 / 2)
    else
      pool_index = pool_index + (nblocks2 / 2) + 1
    endif

    do itr = 1, nblocks1
      map1(itr) = (itr - 1) * block_size1 + 1
    enddo
    do itr = 1, nblocks2
      map2(itr) = (itr - 1) * block_size2 + 1
    enddo

    call matd_create_irreg_blockcyclic_gellan_(&
      self, dim1, dim2, map1, map2, comm, mem_pool, pool_index, iend, .false.)
  end subroutine matd_create_reg_gellan_


  !
  ! [GELLAN用]
  ! イレギュラーブロック(サイクリック)分散行列を生成するサブルーチン
  !
  subroutine matd_create_irreg_gellan_(self, dim1, dim2, map1, map2, comm, &
      mem_pool, ibegin, iend, is_blockcyclic)
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, comm
    integer, intent(in) :: map1(:), map2(:)
    double precision, target, intent(in) :: mem_pool(:)
    integer, intent(in) :: ibegin
    integer, intent(out) :: iend

    logical, intent(in), optional :: is_blockcyclic
    integer :: ierr, nprocs
    logical :: is_bc

    if (present(is_blockcyclic)) then
      is_bc = is_blockcyclic
    else
      is_bc = .false.
    endif
    if (.not. is_bc) then
      call mpi_comm_size(comm, nprocs, ierr)
      if (nprocs < size(map1) * size(map2)) then
        write (*, *) "ERROR : Creating a matrix of irregular block distribution"
        stop
      endif
    endif

    call matd_create_irreg_blockcyclic_gellan_(&
      self, dim1, dim2, map1, map2, comm, mem_pool, ibegin, iend, .true.)
    return
  end subroutine matd_create_irreg_gellan_


  !
  ! [GELLAN用]
  ! ScaLAPACK型ブロックサイクリック分散行列を生成するサブルーチン
  !
  subroutine matd_create_scalapack_gellan_(self, dim1, dim2, block_size1, block_size2, &
      procgrid1, procgrid2, comm, mem_pool, ibegin, iend)
    use iso_c_binding
    type(matd_matrix), intent(out) :: self
    integer, intent(in) :: dim1, dim2, block_size1, block_size2, procgrid1, procgrid2, comm
    integer, pointer :: map1(:), map2(:)
    double precision, target, intent(in) :: mem_pool(:)
    integer, intent(in) :: ibegin
    integer, intent(out) :: iend

    integer :: nblocks1, nblocks2, itr
    type(c_ptr) :: cptr
    integer :: pool_index

    if (mod(dim1, block_size1) == 0) then
      nblocks1 = dim1 / block_size1
    else
      nblocks1 = dim1 / block_size1 + 1
    endif

    if (mod(dim2, block_size2) == 0) then
      nblocks2 = dim2 / block_size2
    else
      nblocks2 = dim2 / block_size2 + 1
    endif

    pool_index = ibegin
    ! マップ情報の領域もプールから取る
    cptr = c_loc(mem_pool(pool_index))
    call c_f_pointer(cptr, map1, shape=[nblocks1])
    ! integerは4バイトなのでインデックスは要素数の半分加算
    if (mod(nblocks1, 2) == 0) then
      pool_index = pool_index + (nblocks1 / 2)
    else
      pool_index = pool_index + (nblocks1 / 2) + 1
    endif

    cptr = c_loc(mem_pool(pool_index))
    call c_f_pointer(cptr, map2, shape=[nblocks2])
    ! integerは4バイトなのでインデックスは要素数の半分加算
    if (mod(nblocks2, 2) == 0) then
      pool_index = pool_index + (nblocks2 / 2)
    else
      pool_index = pool_index + (nblocks2 / 2) + 1
    endif

    do itr = 1, nblocks1
      map1(itr) = (itr - 1) * block_size1 + 1
    enddo
    do itr = 1, nblocks2
      map2(itr) = (itr - 1) * block_size2 + 1
    enddo

    call matd_create_irreg_scalapack_gellan_(&
      self, dim1, dim2, map1, map2, procgrid1, procgrid2, comm, mem_pool, pool_index, iend)
  end subroutine matd_create_scalapack_gellan_


  !
  ! [GELLAN用]
  ! 行列を開放するサブルーチン
  !
  subroutine matd_destroy_gellan_(self)
    type(matd_matrix), intent(out) :: self
    integer :: ierr

    call mpi_win_free(self%win, ierr)
    return
  end subroutine matd_destroy_gellan_
