#include "defs.F90"

program tristan
  use m_globalnamespace
  use m_communications
  use m_initialize
  use m_userfile
  use m_mainloop
  implicit none
  !----- main code --------------------------
  real(mprec) :: myvar

  call initializeAll()
  call getInput('particles', 'delgam2', myvar, toMPREC(123))
  print *, myvar, kind(myvar)

  !..... main code ..........................

  stop
end program tristan
