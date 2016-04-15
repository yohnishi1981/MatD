
module matd_common_module
  implicit none
  private
  public :: matd_sp, matd_dp

  integer, parameter :: matd_sp = selected_real_kind(6, 37)
  integer, parameter :: matd_dp = selected_real_kind(15, 307)
end module matd_common_module


module matd_matrix_int_module
  use mpi
  use matd_common_module
  implicit none
#define MATD_ELEMTYPE integer
#define MATD_MPI_TYPE MPI_INT
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
end module matd_matrix_int_module


module matd_matrix_int8_module
  use mpi
  use matd_common_module
  implicit none
#define MATD_ELEMTYPE integer(8)
#define MATD_MPI_TYPE MPI_LONG_LONG_INT
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
end module matd_matrix_int8_module


module matd_matrix_real4_module
  use mpi
  use matd_common_module
  implicit none
#define MATD_ELEMTYPE real(matd_sp)
#define MATD_MPI_TYPE MPI_REAL4
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
end module matd_matrix_real4_module


module matd_matrix_real8_module
  use mpi
  use matd_common_module
  implicit none
#define MATD_ELEMTYPE real(matd_dp)
#define MATD_MPI_TYPE MPI_REAL8
#include "matd.h"
#undef MATD_ELEMTYPE
#undef MATD_MPI_TYPE
end module matd_matrix_real8_module


module matd
  use matd_common_module
  use matd_matrix_int_module, matd_int_matrix => matd_matrix
  use matd_matrix_int8_module, matd_int8_matrix => matd_matrix
  use matd_matrix_real4_module, matd_real4_matrix => matd_matrix
  use matd_matrix_real8_module, matd_real8_matrix => matd_matrix
  implicit none
end module matd
