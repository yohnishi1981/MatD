PROGRAM MAIN
  USE MatD
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  TYPE(MatD_Int_matrix) :: M
  INTEGER :: Map1(5) = (/1,3,6,7,10/)
  INTEGER :: Map2(4) = (/1,3,6,10/)
  INTEGER,PARAMETER :: Dim1 = 10, Dim2 = 11
  INTEGER(4) :: MyRank, NProcs, IErr
  INTEGER :: I, IBegin, IEnd
  INTEGER :: Buf(110)
  INTEGER,POINTER :: Ptr(:)

  DOUBLE PRECISION, ALLOCATABLE :: X(:)   !! Kind of memory pool in Gellan.
  ALLOCATE(X(10000))

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,MyRank,IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,NProcs,IErr)
  IF (NProcs /= 6) THEN
    IF (MyRank == 0) WRITE(6,'(A)') " ERROR! NProcs must be 6."
    STOP
  ENDIF

  IBegin = 200 + 1 ! Kind of "CALL MemTop(IBegin)"

!!  Make an irregular block cyclic distribution using X.
  CALL MatD_Create_irreg_gellan( &
    M,Dim1,Dim2,Map1,Map2,MPI_COMM_WORLD,X,IBegin,IEnd,.TRUE.)
  CALL MatD_Fence(M)
!!  In Gellan, MemAdv(IEnd-IBegin) is required after Create. 
  CALL MatD_Print_info(M)
  CALL MatD_Fence(M)

!!  Print out the memory space on each rank.
!!  ex.) Rank 0 owns 23 integer(4) numbers and thus 
!!  the index advances by 23/2 + 1 = 12.
!!  Also, the memory for the map information 
!!  is added (9/2 + 1 = 5) because the map information is 4-byte
!!  integer.
!!  Finally, on rank 0, the index advances by 17.
  DO I = 0, 5
    IF (MyRank == I) THEN
      WRITE(6,'(A,I5,A,I5,A)') "Rank ", I, " Memory = ", IEnd-IBegin, " Words."
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
      WRITE(6,'(A,I5,A,I5)')  "== Rank ==", I, ", NElems =", SIZE(Ptr)
      WRITE(6, '(10I5)') Ptr
    endif
    CALL MatD_Fence(M)
  enddo

  CALL MatD_Destroy_gellan(M)
  CALL MPI_Finalize(IErr)
END PROGRAM MAIN
