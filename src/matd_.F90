MODULE MatD_Common_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: MatD_SP, MatD_DP

  INTEGER,PARAMETER :: MatD_SP = SELECTED_REAL_KIND(6,37)
  INTEGER,PARAMETER :: MatD_DP = SELECTED_REAL_KIND(15,307)
END MODULE MatD_Common_module


MODULE MatD_Matrix_int_module
!!  4-byte integer
  USE MatD_Common_module
  IMPLICIT NONE
  INCLUDE "mpif.h"
#define MATD_ELEMTYPE integer
#define MATD_MPI_TYPE MPI_INTEGER4
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
END MODULE MatD_Matrix_int_module


MODULE MatD_Matrix_int8_module
!!  8-byte integer
  USE MatD_Common_module
  IMPLICIT NONE
  INCLUDE "mpif.h"
#define MATD_ELEMTYPE integer(8)
#define MATD_MPI_TYPE MPI_INTEGER8
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
END MODULE MatD_Matrix_int8_module


MODULE MatD_Matrix_real4_module
  USE MatD_Common_module
  IMPLICIT NONE
  INCLUDE "mpif.h"
#define MATD_ELEMTYPE real(matd_sp)
#define MATD_MPI_TYPE MPI_REAL4
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
END MODULE MatD_Matrix_real4_module


MODULE MatD_matrix_real8_module
  USE MatD_Common_module
  IMPLICIT NONE
  INCLUDE "mpif.h"
#define MATD_ELEMTYPE real(matd_dp)
#define MATD_MPI_TYPE MPI_REAL8
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
END MODULE MatD_Matrix_real8_module


MODULE MatD
  USE MatD_Common_MODULE
  USE MatD_Matrix_int_MODULE,   MatD_Int_matrix   => MatD_Matrix
  USE MatD_Matrix_int8_MODULE,  MatD_Int8_matrix  => MatD_Matrix
  USE MatD_Matrix_real4_MODULE, MatD_Real4_matrix => MatD_Matrix
  USE MatD_Matrix_real8_MODULE, MatD_Real8_matrix => MatD_Matrix
  IMPLICIT NONE
END MODULE MatD
