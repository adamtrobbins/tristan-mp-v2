#include "defs.F90"

program tristan
  use m_globalnamespace
  use m_aux
  use m_initialize
  use m_writeoutput
  use m_finalize
  use m_userfile
  use m_mainloop
  implicit none
  !----- main code --------------------------
  integer :: i

  call initializeAll()
  call mainloop()
  call finalizeAll()

  !..... main code ..........................

  ! stop
end program tristan
