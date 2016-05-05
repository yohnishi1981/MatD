PROGRAM MAIN
  USE matd
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  TYPE(MATD_REAL8_MATRIX) :: M1,M2,M3
  INTEGER,PARAMETER :: NRow = 5, NCol = 5
  INTEGER :: NBRow,NBCol
  INTEGER(4) :: MyRank4, NProcs4, IErr, ICTxt4, MyPRow4, MyPCol4, Info
  INTEGER :: MyRank, NProcs, NPRow, NPCol, MyNBlk, MyNRow, MyNCol, ICTxt
  INTEGER :: MyPRow, MyPCol
  INTEGER :: I, J, IJ, M, N, K
  INTEGER,ALLOCATABLE :: LoR(:),LoC(:),HiR(:),HiC(:)
  DOUBLE PRECISION :: A_G(NRow*NCol),B_G(NRow*NCol),C_G(NRow*NCol)
  DOUBLE PRECISION,POINTER :: A_L(:), B_L(:)
  DOUBLE PRECISION,ALLOCATABLE :: C_L(:)
  DOUBLE PRECISION :: WTS,WTE,CTS,CTE
  INTEGER :: DescA(9), DescB(9), DescC(9)
  INTEGER,EXTERNAL :: NUMROC

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, MyRank4, IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, NProcs4, IErr)
  MyRank = MyRank4
  NProcs = NProcs4

  CALL Calc_NPRC(NProcs,NPRow,NPCol)  
  IF (MyRank == 0) WRITE(6,'(A,I5,A,I5)') " Process Grid ", NPRow, ' X ', NPCol

  A_G(:) = 0.0D0
  B_G(:) = 0.0D0
  IF (MyRank == 0) THEN
    IJ = 0
    DO I = 1, NCol
      DO J = 1, NRow
        IJ = IJ + 1
        A_G(IJ) = DBLE(IJ)
        B_G(IJ) = -DBLE(IJ)
      ENDDO
    ENDDO
  ENDIF
 
  NBRow = INT(CEILING(DBLE(NRow)/DBLE(NPRow)))
  NBCol = INT(CEILING(DBLE(NCol)/DBLE(NPCol)))
 
! Create Global Array M1 and distribute to A_L
  CALL MatD_Create_scalapack(M1,NRow,NCol,NBRow,NBCol,NPRow,NPCol,MPI_COMM_WORLD)
  CALL MatD_Fence(M1)
  IF (MyRank == 0) THEN
    CALL MatD_Put(M1,1,NRow,1,NCol,A_G)
  ENDIF
  CALL MatD_Data(M1,A_L)
  CALL MatD_Fence(M1)

! Create Global Array M2 and distribute to B_L
  CALL Matd_Create_scalapack(M2,NRow,NCol,NBRow,NBCol,NPRow,NPCol,MPI_COMM_WORLD)
  CALL MatD_Fence(M2)
  IF (MyRank == 0) THEN
    CALL MatD_Put(M2,1,NRow,1,NCol,B_G)
  ENDIF
  CALL MatD_Data(M2,B_L)
  CALL MatD_Fence(M2)

! BLACS ROUTINES
  CALL BLACS_Get(0,0,ICTxt4)
  ICTxt = ICTxt4
  CALL BLACS_GridInit(ICTxt,'COL',NPRow,NPCol)
  CALL BLACS_GridInfo(ICTxt,NPRow,NPCol,MyPRow4,MyPCol4)
  MyPRow = MyPRow4
  MyPCol = MyPCol4
  MyNRow = NUMROC(NRow,NBRow,MyPRow,0,NPRow)
  MyNCol = NUMROC(NCol,NBCol,MyPCol,0,NPCol)

  ALLOCATE(C_L(1:MyNRow*MyNCol))

!  WRITE(6,'(A,I5,A,2I5)') " MyRank = ", MyRank, " MyNRow,MyNCol = ", MyNRow, MyNCol

  CALL DescInit(DescA,NRow,NCol,NBRow,NBCol,0,0,ICTxt,MAX(1,MyNRow),Info)
  CALL DescInit(DescB,NRow,NCol,NBRow,NBCol,0,0,ICTxt,MAX(1,MyNRow),Info)
  CALL DescInit(DescC,NRow,NCol,NBRow,NBCol,0,0,ICTxt,MAX(1,MyNRow),Info)

!! Dump the data of M1 on each rank
!  IF (MyRank == 0) WRITE(6,'(A)') "==<< ScaLAPACK IBC >>=="
!  DO I = 0, NProcs-1
!    IF (MyRank == I) THEN
!      WRITE(6,'(/A,I5)') "== Rank ", I
!      CALL DumpDPMat(A_L,MyNRow,MyNCol,6)
!    ENDIF
!    CALL MatD_Fence(m1)
!  ENDDO
!! Dump the data of M2 on each rank
!  IF (MyRank == 0) WRITE(6,'(A)') "==<< ScaLAPACK IBC >>=="
!  DO I = 0, NProcs-1
!    IF (MyRank == I) THEN
!      WRITE(6,'(/A,I5)') "== Rank ", I
!      CALL DumpDPMat(B_L,MyNRow,MyNCol,6)
!    ENDIF
!    CALL MatD_Fence(M2)
!  ENDDO
  
  M = NRow
  N = NCol
  K = NCol
  CALL DGEMM('N','N',M,N,K,1.0D0,A_G,M,B_G,N,0.0D0,C_G,M)
  IF (MyRank == 0) THEN
    WRITE(6,'(A)') " Reference Result "
    CALL DumpDPMat(C_G,NRow,NCol,6) 
  ENDIF

  CALL PDGEMM('N','N',M,N,K,1.0d0,A_L,1,1,DescA,B_L,1,1,DescB,0.0D0,C_L,1,1,DescC)
! Dump the data of C_L on each rank
  IF (MyRank == 0) WRITE(6,'(/A)') "==<< PDGEMM Result >>=="
  DO I = 0, NProcs-1
    IF (MyRank == I) THEN
      WRITE(6,'(/A,I5)') "== Rank ", I
      CALL DumpDPMat(C_L,MyNRow,MyNCol,6)
    ENDIF
    CALL MPI_Barrier(MPI_COMM_WORLD,IErr)
  ENDDO
 
! Destroy matrices
  CALL MatD_Destroy(M1)
  CALL MatD_Destroy(M2)
  CALL BLACS_GridExit(ICTxt)
  CALL MPI_finalize(IErr)
END PROGRAM MAIN


SUBROUTINE DumpIMat(A,N,M,ILU)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,M,ILU
  INTEGER,INTENT(IN) :: A(N,M)
  INTEGER :: I, J, K, L

  DO I = 1, M, 10
    L = MIN(I+9,M)
    WRITE(6,'(4X,A,10I12)')' |',(K, K=I,L)
    WRITE(6,*) REPEAT('-',6+12*(L-I+1))
    DO J = I, N
      WRITE(ILU,'(I4,A,10I12)') J, ' |', (A(J,K), K=I,L)
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE DumpIMat

SUBROUTINE DumpDPMat(A,N,M,ILU)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,M,ILU
  DOUBLE PRECISION,INTENT(IN) :: A(N,M)
  INTEGER :: I, J, K, L

  DO I = 1, M, 10
    L = MIN(I+9,M)
    WRITE(6,'(4X,A,10I12)')' |',(K, K=I,L)
    WRITE(6,*) REPEAT('-',6+12*(L-I+1))
    DO J = I, N
      WRITE(ILU,'(I4,A,10F12.5)') J, ' |', (A(J,K), K=I,L)
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE DumpDPMat

SUBROUTINE Calc_NPRC(NProcs,NPRow,NPCol)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NProcs
    INTEGER, INTENT(OUT) :: NPRow, NPCol
    INTEGER :: I, IPRow
    LOGICAL :: LgEN  !   LOGICAL FOR ELEMENTARY NUMBER

    IF (NProcs == 1) THEN
      NPRow = 1
      NPCol = 1
      RETURN
    ELSE IF (NProcs == 2) THEN
      NPRow = 1
      NPCol = 2
      RETURN
    ELSE
      DO I = 2, NProcs-1
        IF (MOD(NProcs,I) == 0) THEN
          LgEN = .FALSE.
          EXIT
        ENDIF
        LgEN = .TRUE.
      ENDDO
    ENDIF

    IF (LgEN) THEN
      NPRow = 1
      NPCol = NProcs
    ELSE
      NPRow = INT(DSQRT(DBLE(NProcs)))
      DO I = 0, NPRow-1
        IPRow = NPRow - I
        NPCol = NProcs/IPRow
        IF (IPRow*NPCol == NProcs) THEN
          NPRow = IPRow
          RETURN
        ELSE
          CYCLE
        ENDIF
      ENDDO
    ENDIF

    RETURN
END SUBROUTINE Calc_NPRC


SUBROUTINE Wall_and_CPU_Time(WALL,CPU)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: WALL,CPU
    INTEGER(8) :: T, T_RATE, T_MAX

    CALL CPU_TIME(CPU)
    CALL SYSTEM_CLOCK(T, T_RATE, T_MAX)
    WALL = DBLE(T)/DBLE(T_RATE)

END SUBROUTINE

