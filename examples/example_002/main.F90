! ブロックサイズが2x3の10x11のブロックサイクリック分散行列
! を生成し、1から110までの数値を順番に代入し、
! 各プロセスが保持している要素を出力して確認するプログラム
program main
!  use mpi
  use matd
  implicit none
  include 'mpif.h'
#if defined (INT8)
  type(matd_int8_matrix) :: m
#else
  type(matd_int_matrix) :: m
#endif
  integer :: ierr, i, myrank
  integer :: buf(110)
  integer, pointer :: ptr(:)

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)

  ! ブロックサイズが2x3の10x11の行列を生成
  ! 末尾のオプショナル引数を.true.とする必要がある
  call matd_create_reg(m, 10, 11, 2, 3, &
    mpi_comm_world, .true.)

  ! フェンスを置く
  call matd_fence(m)
  
  ! 分散行列へ代入するデータをランク0で生成し、
  ! Put処理で行列全体へ代入する
  if (myrank == 0) then
    do i = 1, 110
      buf(i) = i
    enddo
    ! 列方向で1-10、行方向で1-11の範囲に代入
    call matd_put(m, 1, 10, 1, 11, buf)
  endif
  call matd_data(m, ptr)
  
  ! フェンスを置いてput操作の待ち合わせ
  call matd_fence(m)

  ! 各プロセスが保持している要素を出力
  do i = 0, 5
    if (myrank == i) then
      print *, "== Rank ", i
      print *, ptr
    endif
    call matd_fence(m)
  enddo

  ! 行列の破壊処理
  call matd_destroy(m)
  
  call mpi_finalize(ierr)
end program main
