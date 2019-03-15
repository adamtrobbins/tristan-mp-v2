#include "../defs.F90"

module m_communications
  use m_globalnamespace
  implicit none
  #ifdef MPI
    include "mpif.h"
  #endif

  ! mpi variables
  integer           :: mpi_rank, mpi_size, mpi_statsize

  integer           :: sizex, sizey, sizez
end module m_communications
