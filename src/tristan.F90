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
  integer :: ierr

  call initializeAll()
  ! call fillGhosts()
  ! if (mpi_rank .eq. 3) print *, ex(-1, 0, 0)
  call mainloop()
  call finalizeAll()

  !..... main code ..........................

  ! stop
end program tristan
