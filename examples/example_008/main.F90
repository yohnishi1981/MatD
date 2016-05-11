PROGRAM MAIN
  USE MatD
  IMPLICIT NONE
  INCLUDE 'mpif.h'
#if defined (INT8)
  TYPE(MatD_Int8_matrix) :: M1, M2
#else
  TYPE(MatD_Int_matrix) :: M1, M2
#endif
  INTEGER(4) :: IErr, MyRank, NProcs
  INTEGER :: I, J, NBlk_Own
  INTEGER :: Map1(1) = (/1/)
  INTEGER :: Map2(6) = (/1, 4, 6, 8, 9, 11/)
  INTEGER,POINTER :: Ptr(:)
  INTEGER,ALLOCATABLE :: Lows1(:), Highs1(:), Lows2(:), Highs2(:)

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,MyRank,IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,NProcs,IErr)
  IF (NProcs /= 6) THEN
    IF (MyRank == 0) WRITE(6,'(A)') " Error! NProcs must be 6."
    CALL MPI_Finalize(IErr)
    STOP
  ENDIF

  CALL MatD_Create_irreg(M1,10,11,Map1,Map2,MPI_COMM_WORLD)
  CALL MatD_Data(M1,Ptr)

  Ptr = MyRank

  CALL MatD_Create_scalapack(M2,10,11,2,2,2,3,MPI_COMM_WORLD)
  CALL MatD_Fence(M1)
  CALL MatD_Fence(M2)
  CALL MatD_Copy(M2, M1)    !! MatD_Copy(B,A) is A => B
  CALL MatD_Fence(M1)
  CALL MatD_Fence(M2)

  IF (MyRank == 0) WRITE(6,'(A)') " === Matrix M1 (Irregular Block) === "
  CALL MatD_Print_info(M1)
  DO I = 0, 5
    IF (MyRank == I) THEN
      CALL MatD_get_nblocks_owned(M1,MyRank,NBlk_Own)
      WRITE(6,'("Rank [",I2,"] owns "I2," blocks")') I, NBlk_Own
      ALLOCATE(Lows1(NBlk_Own), Highs1(NBlk_Own), &
               Lows2(NBlk_Own), Highs2(NBlk_Own))
      CALL MatD_Distribution(M1,MyRank,Lows1,Highs1,Lows2,Highs2)
      DO J = 1, NBlk_Own
        WRITE(6,'(" Block ",I2,", Positions ",5I5)') J,Lows1(J), Highs1(J), Lows2(J), Highs2(J)
      ENDDO
      DEALLOCATE(Lows1,Highs1,Lows2,Highs2)
    ENDIF
    CALL MatD_Fence(M1)
  ENDDO
!!  Dump Ptr of M1
  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5)') "Rank : ", I
      WRITE(6,'(10I5)') Ptr
    ENDIF
    CALL MatD_Fence(M1)
  ENDDO


  IF (MyRank == 0) WRITE(6,'(/A)') " === Matrix M2 (ScaLAPACK) === "
  CALL MatD_Print_info(M2)
  DO I = 0, 5
    IF (MyRank == I) THEN
      CALL MatD_get_nblocks_owned(M2,MyRank,NBlk_Own)
      WRITE(6,'("Rank [",I2,"] owns "I2," blocks")') I, NBlk_Own
      ALLOCATE(Lows1(NBlk_Own), Highs1(NBlk_Own), &
               Lows2(NBlk_Own), Highs2(NBlk_Own))
      CALL MatD_Distribution(M2,MyRank,Lows1,Highs1,Lows2,Highs2)
      DO J = 1, NBlk_Own
        WRITE(6,'(" Block ",I2,", Positions ",5I5)') J,Lows1(J), Highs1(J), Lows2(J), Highs2(J)
      ENDDO
      DEALLOCATE(Lows1,Highs1,Lows2,Highs2)
    ENDIF
    CALL MatD_Fence(M2)
  ENDDO

  CALL MatD_Data(M2, ptr)
!!  Dump Ptr of M2 after copy
  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5)') "Rank : ", I
      WRITE(6,'(10I5)') Ptr
    ENDIF
    CALL MatD_Fence(M2)
  ENDDO

  CALL MatD_Destroy(M1)
  CALL MatD_Destroy(M2)

  CALL MPI_Finalize(IErr)
END PROGRAM MAIN
