!!  1. Create 10x11 block distribution matrix, whose block size is 5x4.
!!  2. Insert the number into it from 1 to 110.
!!  3. Print out the data in each rank.
!!  * NProcs must be 6.
program main
  use mpi
  use matd
  implicit none

  type(matd_int_matrix) :: m
  integer,parameter :: nelm = 10*11
  integer(4) :: ierr, myrank, nprocs
  integer :: i
  integer :: buf(nelm)
  integer, pointer :: ptr(:)

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  call mpi_comm_size(mpi_comm_world, nprocs, ierr)
  if (nprocs /= 6) then
    if (myrank == 0) write(6,'(a)') " NProcs must be 6 in this program. Program stop."
    call mpi_finalize(ierr)
    stop
  endif

!!  1. Create 10x11 block distribution matrix, whose block size is 5x4.
!!  * NProcs must be 6.
  call matd_create_reg(m, 10, 11, 5, 4, mpi_comm_world)

!!  Fence required
  call matd_fence(m)

!!  Generate data to put distributed matrix only on rank 0.
  buf(1:nelm) = 0
  if (myrank == 0) then
    do i = 1, nelm
      buf(i) = i
    enddo
!!  The size of distributed matrix is 10x11.
    call matd_put(m, 1, 10, 1, 11, buf)
  endif
!!  Mapping the data in m to ptr.
  call matd_data(m, ptr)

!!  Fence required to wait for put operation.
  call matd_fence(m)

!!  Print out the elements on each rank.
  do i = 0, 5
    if (myrank == i) then
      print *, "== Rank ", i
      print *, ptr
    endif
    call matd_fence(m)
  enddo

!!  Destroy the matrix
  call matd_destroy(m)

  call mpi_finalize(ierr)
end program main
