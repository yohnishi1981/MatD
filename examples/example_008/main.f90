program main
  use mpi
  use matd
  implicit none

  type(matd_int_matrix) :: m1, m2
  integer :: ierr, i, myrank, j, nblocks_owned
  integer :: map1(1) = (/1/)
  integer :: map2(6) = (/1, 4, 6, 8, 9, 11/)
  integer, pointer :: ptr(:)
  integer, allocatable :: lows1(:), highs1(:), lows2(:), highs2(:)


  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)

  call matd_create_irreg(m1, 10, 11, map1, map2, mpi_comm_world)
  call matd_data(m1, ptr)

  ptr = myrank

  call matd_create_scalapack(m2, 10, 11, 2, 2, 2, 3, mpi_comm_world)
  call matd_fence(m1)
  call matd_fence(m2)
  call matd_copy(m2, m1)
  call matd_fence(m1)
  call matd_fence(m2)

  do i = 0, 5
    if (myrank == i) then
      print *, "[RANK]=", i
      call matd_get_nblocks_owned(m1, myrank, nblocks_owned)
      print *, "NBlocks ", nblocks_owned
      allocate(lows1(nblocks_owned), highs1(nblocks_owned), &
               lows2(nblocks_owned), highs2(nblocks_owned))
      call matd_distribution(m1, myrank, lows1, highs1, lows2, highs2)
      do j = 1, nblocks_owned
        print *, "[", j, "] = ", lows1(j), highs1(j), lows2(j), highs2(j)
      enddo
      deallocate(lows1, highs1, lows2, highs2)
    endif
    call matd_fence(m1)
  enddo

  do i = 0, 5
    if (myrank == i) then
      print *, "[RANK]=", i
      call matd_get_nblocks_owned(m2, myrank, nblocks_owned)
      print *, "NBlocks ", nblocks_owned
      allocate(lows1(nblocks_owned), highs1(nblocks_owned), &
               lows2(nblocks_owned), highs2(nblocks_owned))
      call matd_distribution(m2, myrank, lows1, highs1, lows2, highs2)
      do j = 1, nblocks_owned
        print *, "[", j, "] = ", lows1(j), highs1(j), lows2(j), highs2(j)
      enddo
      deallocate(lows1, highs1, lows2, highs2)
    endif
    call matd_fence(m2)
  enddo

  call matd_data(m2, ptr)

  do i = 0, 5
    if (myrank == i) then
      print *, "== Rank ", i
      print *, ptr
    endif
    call matd_fence(m2)
  enddo


  call matd_destroy(m1)
  call matd_destroy(m2)

  call mpi_finalize(ierr)
end program main
