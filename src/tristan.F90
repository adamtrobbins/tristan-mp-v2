#include "defs.F90"

program tristan
  use m_globalnamespace
  use m_initialize
  use m_finalize
  use m_userfile
  use m_mainloop
  implicit none
  !----- main code --------------------------

  call initializeAll()

  print *, 'RNK', my_rank, this_meshblock%sx, this_meshblock%sy, this_meshblock%sz,&
         & this_meshblock%x0, this_meshblock%y0, this_meshblock%z0

  call finalizeAll()

  !..... main code ..........................

  ! stop
end program tristan
