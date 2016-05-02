PROGRAM MAIN
!  use MPI
  USE matd
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  TYPE(MATD_REAL8_MATRIX) :: M1
  INTEGER,PARAMETER :: NRow = 20000, NCol = 20000
  INTEGER :: NBRow,NBCol
  INTEGER(4) :: MyRank4,NProcs4,IErr
  INTEGER :: MyRank, NProcs
  INTEGER :: NPRow, NPCol, MyNBlk, MyNRow, MyNCol
  INTEGER :: I, J, IJ
  INTEGER,ALLOCATABLE :: LoR(:),LoC(:),HiR(:),HiC(:)
  DOUBLE PRECISION :: Buf(NRow*NCol)
  DOUBLE PRECISION,POINTER :: Ptr(:)
  DOUBLE PRECISION :: WTS,WTE,CTS,CTE

  CALL MPI_Init(IErr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, MyRank4, IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, NProcs4, IErr)
  MyRank = MyRank4
  NProcs = NProcs4

  CALL Calc_NPRC(NProcs,NPRow,NPCol)  
!  NPRow = 1
!  NPCol = NProcs
  NBRow = NRow/NPRow
  NBCol = NCol/NPCol
!  NBRow = 2500
!  NBCol = 2500
  IF (MyRank == 0) THEN
    WRITE(6,'(A,I5,A,I5)') " Process Grid for BLACS = ", NPRow, " x ",NPCol
    WRITE(6,'(A,I5,A,I5)') " Size of Global Matrix = ", NRow, " X", NCol 
    WRITE(6,'(A,I5,A,I5)') " Size of Block Size: Row = ", NBRow, ", Col = ", NBCol 
  ENDIF
  CALL Wall_and_CPU_Time(WTS,CTS)
  CALL MatD_Create_scalapack(M1,NRow,NCol,NBRow,NBCol,NPRow,NPCol,MPI_COMM_WORLD)
  CALL Wall_and_CPU_Time(WTE,CTE)
  IF (MyRank == 0) WRITE(6,'(A,F10.5)') " Wall time for creating M1 (sec.)", WTE-WTS
  CALL MatD_Fence(m1)
  CALL Wall_and_CPU_Time(WTS,CTS)

  IF (MyRank == 0) THEN
    IJ = 0
    DO I = 1, NCol
      DO J = 1, NRow
        IJ = IJ + 1
        Buf(IJ) = DBLE(IJ)
      ENDDO
    ENDDO
    CALL Wall_and_CPU_Time(WTS,CTS)
    CALL matd_put(M1,1,NRow,1,NCol,Buf)
    CALL Wall_and_CPU_Time(WTE,CTE)
    WRITE(6,'(A,F10.5)') " Wall time for putting data to M1 (sec.)", WTE-WTS
  ENDIF
  CALL Wall_and_CPU_Time(WTS,CTS)
  CALL matd_data(M1, Ptr)
  CALL matd_fence(M1)
  CALL Wall_and_CPU_Time(WTE,CTE)
  IF (MyRank == 0) WRITE(6,'(A,F10.5)') " Wall time for getting data from M1 (sec.)", WTE-WTS

! Calculate the size of Local Matrix
  CALL matd_get_nblocks_owned(M1,MyRank,MyNBlk)
  ALLOCATE(LoR(1:MyNBlk),LoC(1:MyNBlk),HiR(1:MyNBlk),HiC(1:MyNBlk))
  CALL matd_distribution(M1,MyRank,LoR,HiR,LoC,HiC)
!  DO I = 0, NProcs-1
!    IF (MyRank == I) THEN
!      WRITE(6,'(A,I5)') " Rank = ", MyRank
!      DO J = 1, MyNBlk
!        WRITE(6,'(A,I5,A,4I5)') "  MyBlock ",J," Range",LoR(J),HiR(J),LoC(J),HiC(J)
!      ENDDO
!    ENDIF
!  ENDDO 

!! Dump the data of M1 on each rank
!  IF (MyRank == 0) WRITE(6,'(A)') "==<< ScaLAPACK IBC >>=="
!  DO I = 0, NProcs-1
!    IF (MyRank == I) THEN
!      WRITE(6,'(/A,I5)') "== Rank ", I
!!      CALL DumpDPMat(Ptr,MyNRow,MyNCol,6)
!      WRITE(6,'(10F12.5)') Ptr
!    ENDIF
!    CALL matd_fence(m1)
!  ENDDO


! Destroy matrices
  CALL matd_destroy(M1)
  CALL MPI_finalize(IErr)
END PROGRAM MAIN


subroutine dumpimat(a,n,m,ilu)
  implicit none
  INTEGER,intent(in) :: n,m,ilu
  INTEGER,intent(in) :: a(n,m)
  INTEGER :: i, j, k, l

  do i = 1, m, 10
    l = min(i+9,m)
    write(6,'(4X,A,10i12)')' |',(k, k=i,l)
    write(6,*) REPEAT('-',6+12*(l-i+1))
    do j = i, n
      write(ilu,'(i4,a,10i12)') j, ' |', (a(j,k), k=i,l)
    enddo
  enddo

  return
end subroutine dumpimat

SUBROUTINE DumpDPMat(A,N,M,ILU)
  implicit none
  INTEGER,intent(in) :: n,m,ilu
  DOUBLE PRECISION,intent(in) :: a(n,m)
  INTEGER :: i, j, k, l

  do i = 1, m, 10
    l = min(i+9,m)
    write(6,'(4X,A,10i12)')' |',(k, k=i,l)
    write(6,*) REPEAT('-',6+12*(l-i+1))
    do j = i, n
      write(ilu,'(i4,a,10f12.5)') j, ' |', (a(j,k), k=i,l)
    enddo
  enddo

  return
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

