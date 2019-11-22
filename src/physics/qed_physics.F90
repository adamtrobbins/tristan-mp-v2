#include "../defs.F90"

#ifdef QED

module m_qedphysics
  use m_globalnamespace
  use m_bwpairproduction
  implicit none

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine QEDstep(timestep)
    implicit none
    integer, intent(in) :: timestep
    if (modulo(timestep, BW_interval) .eq. 0) then
      call bwPairProduction()
    end if
  end subroutine QEDstep

end module m_qedphysics

#endif
