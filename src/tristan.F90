#include "defs.F90"

program tristan
  use m_globalnamespace
  use m_aux
  use m_initialize
  use m_writeoutput
  use m_finalize
  use m_userfile
  use m_mainloop

  use m_fldsolver
  use m_fields
  implicit none
  !----- main code --------------------------
  integer :: status, ierr, ind1, ind2, ind3, s, p
  integer :: mpi_sendto, mpi_recvfrom

  call initializeAll()
  call mainloop()
  call finalizeAll()

  !..... main code ..........................
end program tristan
