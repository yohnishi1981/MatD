PROGRAM MAIN
!!  Copy from the matrix of example_005 to one of example_001
  USE MatD
  IMPLICIT NONE
  INCLUDE 'mpif.h'
#if defined (INT8)
  TYPE(MatD_Int8_matrix) :: M1, M2
#else
  TYPE(MatD_Int_matrix) :: M1, M2
#endif
  INTEGER(4) :: IErr, MyRank, NProcs
  INTEGER :: I
  INTEGER :: Buf(110)
  INTEGER,POINTER :: Ptr(:)

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,MyRank,IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,NProcs,IErr)
  IF (NProcs /= 6) THEN
    IF (MyRank == 0) WRITE(6,'(A)') " Error! NProcs must be 6."
    CALL MPI_Finalize(IErr)
    STOP
  ENDIF

  CALL MatD_Create_reg(M1,10,11,5,4,MPI_COMM_WORLD)
  CALL MatD_Create_scalapack(M2,10,11,2,2,2,3,MPI_COMM_WORLD)
  CALL MatD_Fence(M1)

  IF (MyRank == 0) THEN
    DO I = 1, 110
      Buf(I) = I
    ENDDO
    CALL MatD_Put(M1,1,10,1,11,Buf)
  ENDIF
  CALL MatD_Data(M1,Ptr)

  CALL MatD_Fence(M1)

!!  Output the data on rank 0 of matrix M1.
  IF (MyRank == 0) WRITE(6,'(A)') "==<< REGULAR BLOCK "
  CALL MatD_Print_Info(M1)
  CALL MatD_Fence(M1)
  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5)') "== Rank ", I
      WRITE(6,'(10I5)') Ptr
    ENDIF
    CALL MatD_Fence(M1)
  ENDDO

!!  Copy
  CALL MatD_Fence(M1)
  CALL MatD_Fence(M2)
  CALL MatD_Copy(M2,M1)   !! MatD_Copy(B,A) is A => B 
  CALL MatD_Fence(M1)
  CALL MatD_Fence(M2)

  CALL MatD_Data(M2,Ptr)
!!  Output the data on rank 0 of matrix M2.
  IF (MyRank == 0) WRITE(6,'(A)') "==<< SCALAPACK BLOCKCYCLIC "
  CALL MatD_Print_Info(M2)
  CALL MatD_Fence(M2)
  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5)') "== Rank ", I
      WRITE(6,'(10I5)') Ptr
    ENDIF
    CALL MatD_Fence(M2)
  ENDDO
  
!!  Destroy before quit.
  CALL MatD_Destroy(M1)
  CALL MatD_Destroy(M2)

  CALL MPI_Finalize(IErr)
END PROGRAM MAIN

