#include "defs.F90"

program tristan
  use m_globalnamespace
  use m_communications
  use m_initialize
  use m_userfile
  use m_mainloop
  implicit none
  !----- main code --------------------------

  call initializeAll()

  !..... main code ..........................

  stop
end program tristan
