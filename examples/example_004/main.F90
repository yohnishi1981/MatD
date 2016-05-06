PROGRAM MAIN
!!
!!  1. Make an irregular block cyclic distribution matrix
!!     whose size is 10 x 11 with the map information of
!!     (1,3,6,7,10) and (1,3,6,10).
!!  2. Insert the data from 1 to 110.
!!  3. Dump the data on each rank.
!!
  USE MatD
  IMPLICIT NONE
  INCLUDE 'mpif.h'
#if defined (INT8)
  TYPE(MatD_Int8_matrix) :: M
#else
  TYPE(MatD_Int_matrix) :: M
#endif
  INTEGER(4) :: IErr, MyRank, NProcs
  INTEGER :: I,K,Buf(110)
  INTEGER,POINTER :: Ptr(:)
!! Map Information
  INTEGER :: Map1(5) = (/1,3,6,7,10/)
  INTEGER :: Map2(4) = (/1,3,6,10/)

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,MyRank,IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,NProcs,IErr)
  IF (NProcs /= 6) THEN
    IF (MyRank == 0) WRITE(6,'(A)') " Error! NProcs must be 6."
    CALL MPI_Finalize(IErr)
    STOP
  ENDIF

!!  Make a matrix whose size is 10 x 11 with Map information of Map1, Map2.
!!  The last argument must be ".TRUE."
  CALL MatD_Create_irreg(M,10,11,Map1,Map2,MPI_COMM_WORLD,.TRUE.)

!!  Fence required
  CALL MatD_Fence(M)

!!  Generate the data on rank 0 and insert by Put.
  IF (MyRank == 0) THEN
    DO I = 1, 110
      Buf(I) = I
    ENDDO
!!  1 - 10 to column direction. 1 - 11 to row direction.
    CALL MatD_Put(M,1,10,1,11,Buf)
  ENDIF
  CALL MatD_Data(M,Ptr)

!!  Fence required to wait for put.
  CALL MatD_Fence(M)

!!  Dump the data.
  IF (MyRank == 0) THEN
    WRITE(6,'(A)') " Global data "
    DO I = 1, 10
      WRITE(6,'(11I5)') (Buf(I+(K-1)*10), K = 1,11)
    ENDDO
    WRITE(6,'(A,10I5)') " Map1 ", Map1
    WRITE(6,'(A,10I5)') " Map2 ", Map2
  ENDIF

  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5)') " == Rank ", I
      WRITE(6,'(10I5)') Ptr
    ENDIF
    CALL MatD_Fence(M)
  ENDDO

!!  Destroy before quit
  CALL MatD_Destroy(M)

  CALL MPI_Finalize(IErr)
END PROGRAM MAIN
