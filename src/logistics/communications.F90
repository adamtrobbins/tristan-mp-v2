#include "../defs.F90"

module m_communications
  use m_globalnamespace
  implicit none

  ! mpi variables
  integer           :: my_rank, size0, statsize
  integer, private  :: ierr

  integer           :: sizex, sizey, sizez
end module m_communications
