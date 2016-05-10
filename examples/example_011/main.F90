PROGRAM MAIN
  USE MatD
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  TYPE(MatD_Real8_matrix) :: M ! 8-byte double precision
  INTEGER,PARAMETER :: MaxMem = 100000
  INTEGER :: Map1(5) = (/1,3,6,7,10/)
  INTEGER :: Map2(4) = (/1,3,6,10/)
  INTEGER,PARAMETER :: Dim1 = 10, Dim2 = 11
  INTEGER(4) :: IErr,MyRank,NProcs
  INTEGER :: I, IBegin, IEnd
  DOUBLE PRECISION :: Buf(110)
  DOUBLE PRECISION,POINTER :: Ptr(:)

  DOUBLE PRECISION,ALLOCATABLE :: X(:) ! As a memory pool of Gellan
  ALLOCATE(X(MaxMem))

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,MyRank,IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,NProcs,IErr)
  IF (NProcs /= 6) THEN
    IF (MyRank == 0) WRITE(6,'(A)') " ERROR! NProcs must be 6"
    STOP
  ENDIF

  IBegin = 200 + 1 ! Kind of "CALL MemTop(IBegin)" in Gellan

!!  Make an irregular block cyclic distributed matrix using X.
  CALL MatD_Create_irreg_gellan( &
    M,Dim1,Dim2,Map1,Map2,MPI_COMM_WORLD,X,IBegin,IEnd,.TRUE. &
  )
  CALL MatD_Fence(M)
  CALL MatD_Print_info(M)
  CALL MatD_Fence(M)

!!  Print out the memory space of each rank.
!!  For example, in this case, the rank 0 owns 23 8-byte
!!  double precision numbers.  Therefore, the M%Storage
!!  consumes 23 words.  In addition, since the last argument
!!  of the above subroutine is .TRUE., the map information 
!!  also consumes X.  Since the Map information is 4-byte integer
!!  and has 9 elements, 9/2 + 1 = 5 words are consumed.
!!  In total, 23 + 5 = 28 words are consumed on the rank 0.
  DO I = 0, NProcs-1
    IF (MyRank == I) THEN
      WRITE(6,'("Rank = ",I3,", X(",I3,":",I3,")")') I,IBegin,IEnd-1
    ENDIF
    CALL MatD_Fence(M)
  ENDDO

  IF (MyRank == 0) THEN
    DO I = 1, 110
      Buf(I) = DBLE(I)
    ENDDO
    CALL MatD_Put(M,1,10,1,11,Buf)
  ENDIF
  CALL MatD_Data(M,Ptr)

  CALL MatD_Fence(M)

  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I3,A,I5)') "== Rank ==", I, ", SIZE=", SIZE(Ptr)
      WRITE(6,'(10F6.1)') Ptr
    ENDIF
    CALL MatD_Fence(M)
  enddo

  CALL MatD_Destroy_gellan(M)
  CALL MPI_Finalize(IErr)

END PROGRAM MAIN
