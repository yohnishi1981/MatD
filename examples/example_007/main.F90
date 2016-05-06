PROGRAM MAIN
!!
!!  1. Make the matrix of example_004.
!!  2. Dump the number of blocks and block positions of each rank.
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
  INTEGER :: I, J, Buf(110)
  INTEGER,POINTER :: Ptr(:)
  INTEGER :: NBlk_Own
  INTEGER,ALLOCATABLE :: Lows1(:),Highs1(:),Lows2(:),Highs2(:)

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

!!  This is the same as the matrix of example_004
  CALL MatD_Create_irreg(M,10,11,Map1,Map2,MPI_COMM_WORLD,.TRUE.)
  CALL MatD_Fence(M)
  CALL MatD_Print_info(M)
  CALL MatD_Fence(M)

  DO I = 0, 5
    IF (MyRank == I) THEN
      CALL MatD_Get_nblocks_owned(M,MyRank,NBlk_Own)
      WRITE(6,'("Rank [",I2,"] owns "I2," blocks")') I, NBlk_Own
!!  Allocate the arrays to store the position of blocks.
      ALLOCATE(Lows1(NBlk_Own), Highs1(NBlk_Own), &
               Lows2(NBlk_Own), Highs2(NBlk_Own))
      CALL MatD_Distribution(M,MyRank,Lows1,Highs1,Lows2,Highs2)
      DO J = 1, NBlk_Own
        WRITE(6,'(" Block ",I2,", Positions ",5I5)') J,Lows1(J), Highs1(J), Lows2(J), Highs2(J)
      ENDDO
      DEALLOCATE(Lows1, Highs1, Lows2, Highs2)
    ENDIF
    CALL MatD_Fence(M)
  ENDDO

  IF (MyRank == 0) THEN
    DO I = 1, 110
      Buf(I) = I
    ENDDO
    CALL MatD_Put(M,1,10,1,11,Buf)
  ENDIF
  CALL MatD_Data(M,Ptr)

  CALL MatD_Fence(M)

  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5)') "Rank : ", I
      WRITE(6,'(10I5)') Ptr
    ENDIF
    CALL MatD_Fence(M)
  ENDDO

  CALL MatD_Destroy(M)
  CALL MPI_finalize(IErr)
END PROGRAM MAIN
