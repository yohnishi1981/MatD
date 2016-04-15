! example_001の分散方法の行列から
! example_005の分散方法の行列へ要素をコピーする例
program main
  use mpi
  use matd
  implicit none

  type(matd_int_matrix) :: m1, m2
  integer :: ierr, i, myrank
  integer :: buf(110)
  integer, pointer :: ptr(:)

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)

  call matd_create_reg(m1, 10, 11, 5, 4, mpi_comm_world)
  call matd_create_scalapack(m2, 10, 11, 2, 2, 2, 3, mpi_comm_world)
  call matd_fence(m1)

  if (myrank == 0) then
    do i = 1, 110
      buf(i) = i
    enddo
    call matd_put(m1, 1, 10, 1, 11, buf)
  endif
  call matd_data(m1, ptr)

  call matd_fence(m1)

  ! m1について各プロセスが保持している要素を出力
  if (myrank == 0) print *, "==<< REGULAR BLOCK"
  do i = 0, 5
    if (myrank == i) then
      print *, "== Rank ", i
      print *, ptr
    endif
    call matd_fence(m1)
  enddo

  ! コピー処理
  call matd_fence(m1)
  call matd_fence(m2)
  call matd_copy(m2, m1)
  call matd_fence(m1)
  call matd_fence(m2)

  call matd_data(m2, ptr)
  ! m1について各プロセスが保持している要素を出力
  if (myrank == 0) print *, "==<< SCALAPACK BLOCKCYCLIC"
  do i = 0, 5
    if (myrank == i) then
      print *, "== Rank ", i
      print *, ptr
    endif
    call matd_fence(m2)
  enddo
  
  ! 行列の破壊処理
  call matd_destroy(m1)
  call matd_destroy(m2)

  call mpi_finalize(ierr)
end program main

