PROGRAM MAIN
!!
!!  1. Make a block-cyclic distributed matrix whose size is 10 x 11 with 
!!     the block size of 2 x 3.
!!  2. Insert 1 to 110.
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
  INTEGER(4) :: MyRank, NProcs, IErr
  INTEGER :: I
  INTEGER :: Buf(110)
  INTEGER,POINTER :: Ptr(:)

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, MyRank, IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, NProcs, IErr)
  IF (NProcs /= 6) THEN
    IF (MyRank == 0) WRITE(6,'(A)') " ERROR! NProcs must be 6."
    CALL MPI_Finalize(IErr)
    STOP
  ENDIF

!!  Make a matrix whose size is 10 x 11 with the block size of 2 x 3.
!!  The last argument must be ".TRUE."
  CALL MatD_Create_reg(M,10,11,2,3,MPI_COMM_WORLD,.TRUE.)

!!  Fence required.
  CALL MatD_Fence(M)
  CALL MatD_Print_info(M)
  CALL MatD_Fence(M)

!!  Generate the data on rank 0 and insert by Put.  
  IF (MyRank == 0) THEN
    DO I = 1, 110
      Buf(I) = I
    ENDDO
!!  1 - 10 to Row and 1 - 11 to Col
    CALL MatD_Put(M,1,10,1,11,Buf)
  ENDIF
  CALL MatD_Data(M,Ptr)

!!  Fence required to wait Put  
  CALL MatD_Fence(M)

!!  Dump the data
  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5)') " == Rank ", I
      WRITE(6,'(10I5)') Ptr
    ENDIF
    CALL MatD_Fence(M)
  ENDDO

!!  Destroy before quit.
  CALL MatD_Destroy(M)
  
  CALL MPI_Finalize(IErr)
END PROGRAM MAIN
