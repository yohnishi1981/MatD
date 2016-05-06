!!  1. Create 10x11 block distribution matrix, whose block size is 5x4.
!!  2. Insert the number into it from 1 to 110.
!!  3. Print out the data in each rank.
!!  * NProcs must be 6.
PROGRAM MAIN
!  use mpi
  USE MatD
  IMPLICIT NONE
  INCLUDE "mpif.h"
#if defined (INT8)
  TYPE(MatD_Int8_matrix) :: M
#else
  TYPE(MatD_Int_matrix) :: M
#endif
  INTEGER,PARAMETER :: NElm = 10*11
  INTEGER(4) :: IErr, MyRank, NProcs
  INTEGER :: I
  INTEGER :: Buf(NElm)
  INTEGER,POINTER :: Ptr(:)

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,MyRank,IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,NProcs,IErr)
  IF (NProcs /= 6) THEN
    IF (MyRank == 0) WRITE(6,'(A)') " NProcs must be 6 in this program. Program stop."
    CALL MPI_Finalize(IErr)
    STOP
  ENDIF

!!  1. Create 10x11 block distribution matrix, whose block size is 5x4.
!!  * NProcs must be 6.
  CALL MatD_Create_reg(M,10,11,5,4,MPI_COMM_WORLD)

!!  Fence required
  CALL MatD_Fence(M)
  CALL MatD_Print_info(M)
  CALL MatD_Fence(M)
  
!!  Generate data to put distributed matrix only on rank 0.
  Buf(1:NElm) = 0
  IF (MyRank == 0) then
    DO I = 1, NElm
      Buf(I) = I
    ENDDO
!!  The size of distributed matrix is 10x11.
    CALL Matd_Put(M,1,10,1,11,Buf)
  ENDIF
!!  Mapping the data in m to ptr.
  CALL MatD_Data(M,Ptr)

!!  Fence required to wait for put operation.
  CALL MatD_Fence(M)

!!  Print out the elements on each rank.
  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5)') " == Rank ", I
      WRITE(6,'(10I5)') Ptr
    ENDIF
    CALL MatD_Fence(M)
  ENDDO

!!  Destroy the matrix
  CALL MatD_Destroy(M)
  CALL MPI_Finalize(IErr)
END PROGRAM MAIN
