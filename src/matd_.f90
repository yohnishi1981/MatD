module matd_common_module
  implicit none
  private
  public :: matd_sp, matd_dp

  integer, parameter :: matd_sp = selected_real_kind(6, 37)
  integer, parameter :: matd_dp = selected_real_kind(15, 307)
end module matd_common_module


module matd_matrix_int_module
!!  4-byte integer
!  use mpi
  use matd_common_module
  implicit none
  include "mpif.h"
#define MATD_ELEMTYPE integer
#define MATD_MPI_TYPE MPI_INTEGER4
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
end module matd_matrix_int_module


module matd_matrix_int8_module
!!  8-byte integer
!  use mpi
  use matd_common_module
  implicit none
  include "mpif.h"
#define MATD_ELEMTYPE integer(8)
#define MATD_MPI_TYPE MPI_INTEGER8
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
end module matd_matrix_int8_module


module matd_matrix_real4_module
!  use mpi
  use matd_common_module
  implicit none
  include "mpif.h"
#define MATD_ELEMTYPE real(matd_sp)
#define MATD_MPI_TYPE MPI_REAL4
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
end module matd_matrix_real4_module


module matd_matrix_real8_module
!  use mpi
  use matd_common_module
  implicit none
  include "mpif.h"
#define MATD_ELEMTYPE real(matd_dp)
#define MATD_MPI_TYPE MPI_REAL8
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
end module matd_matrix_real8_module


MODULE MatD
  USE MatD_Common_module
  USE MatD_Matrix_int_module, MatD_Int_matrix => MatD_Matrix
  USE MatD_Matrix_int8_module, MatD_Int8_matrix => MatD_Matrix
  USE MatD_Matrix_real4_module, MatD_Real4_matrix => MatD_Matrix
  USE MatD_Matrix_real8_module, MatD_Real8_matrix => MatD_Matrix
  IMPLICIT NONE
END MODULE MatD
