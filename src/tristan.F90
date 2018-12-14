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

  if (my_rank .eq. 2) then
    print *, this_meshblock%ptr%neighbor(-1, 0, 0)%ptr%rnk
    print *, this_meshblock%ptr%neighbor(+1, 0, 0)%ptr%rnk
    print *, this_meshblock%ptr%neighbor(0, -1, 0)%ptr%rnk
    print *, this_meshblock%ptr%neighbor(0, +1, 0)%ptr%rnk
  end if

  call finalizeAll()

  !..... main code ..........................

  ! stop
end program tristan
