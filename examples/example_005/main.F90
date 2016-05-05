PROGRAM MAIN
!!
!!  1. Make a ScaLAPACK-type distribution whose size is 10 x 11
!!     with the 2x2 block size.
!!  2. Insert from 1 to 110 into 1-D Buf array.
!!  3. Output the data on each rank.
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

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,MyRank,IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,NProcs,IErr)
  IF (NProcs /= 6) THEN
    IF (MyRank == 0) WRITE(6,'(A)') " Error! NProcs must be 6."
    CALL MPI_Finalize(IErr)
    STOP
  ENDIF

!!  Make a matrix whose size is 10 x 11 and with Map information of Map1 and 2.
  CALL MatD_Create_scalapack(M,10,11,2,2,2,3,MPI_COMM_WORLD)

!!  Fence required
  CALL MatD_Fence(M)
  
!!  Generate the data on rank 0 and insert by Put.
  IF (MyRank == 0) THEN
    DO I = 1, 110
      Buf(I) = I
    ENDDO
!!  1 - 10 to column direction.  1 - 11 to row direction.
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
