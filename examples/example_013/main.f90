program main
!  use mpi
  use matd
  implicit none
  include 'mpif.h'
  type(matd_int_matrix) :: m

  integer, parameter :: block_size1 = 5, block_size2 = 4
  integer, parameter :: dim1 = 10, dim2 = 11

  integer :: buf(110)
  integer :: ierr, myrank, i, ibegin, iend
  integer, pointer :: ptr(:)

  double precision, allocatable :: X(:) ! GELLANの領域と考える
  allocate(X(10000)) ! メモリプールの確保

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)

  ibegin = 200 ! index=200からX領域を使い始める

  ! X領域を使ってイレギュラーブロックサイクリック分散行列を生成
  call matd_create_reg_gellan( &
    m, dim1, dim2, block_size1, block_size2, &
    mpi_comm_world, X, ibegin, iend, .true. &
  )
  call matd_fence(m)

  ! 各プロセスがX領域のどの部分を利用したかを示す。
  ! 例えばこの分散方法の場合ランク0のプロセスは
  ! 20個の4バイト整数の要素を持つ。
  ! したがって要素分としては20 / 2 = 10
  ! インデックスが進む。
  ! またマップ情報分(4バイト整数*(2+3)) 5 / 2 + 1 = 3
  ! インデックスが進み、計13インデックスが進む。
  do i = 0, 5
    if (myrank == i) then
      print *, iend - ibegin
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
      print *, "== Rank ==", i, ", SIZE=", size(ptr)
      print *, ptr
    endif
    call matd_fence(m)
  enddo

  call matd_destroy_gellan(m)
  call mpi_finalize(ierr)
end program main
