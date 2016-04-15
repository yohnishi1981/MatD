program main
  use mpi
  use matd
  implicit none

  type(matd_int_matrix) :: m
  integer :: ierr, i, myrank, j
  integer :: buf(110)
  integer, pointer :: ptr(:)
  integer :: nblocks_owned
  integer, allocatable :: lows1(:), highs1(:), lows2(:), highs2(:)

  integer :: map1(5) = (/ 1, 3, 6, 7, 10 /)
  integer :: map2(4) = (/ 1, 3, 6, 10 /)

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)

  call matd_create_irreg(m, 10, 11, map1, map2, mpi_comm_world, .true.)
  call matd_fence(m)

  do i = 0, 5
    if (myrank == i) then
      print *, "[RANK]=", i
      call matd_get_nblocks_owned(m, myrank, nblocks_owned)
      print *, "NBlocks ", nblocks_owned
      allocate(lows1(nblocks_owned), highs1(nblocks_owned), &
               lows2(nblocks_owned), highs2(nblocks_owned))
      call matd_distribution(m, myrank, lows1, highs1, lows2, highs2)
      do j = 1, nblocks_owned
        print *, "[", j, "] = ", lows1(j), highs1(j), lows2(j), highs2(j)
      enddo
      deallocate(lows1, highs1, lows2, highs2)
    endif
    call matd_fence(m)
  enddo

  if (myrank == 0) then
    do i = 1, 110
      buf(i) = i
    enddo
    call matd_put(m, 1, 10, 1, 11, buf)
  endif
  call matd_data(m, ptr)

  call matd_fence(m)

  do i = 0, 5
    if (myrank == i) then
      print *, "[Rank] : <<", i, ">>"
      print *, ptr
    endif
    call matd_fence(m)
  enddo

  call matd_destroy(m)
  call mpi_finalize(ierr)
end program main
