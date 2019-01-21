#include "defs.F90"

module m_mainloop
  use m_globalnamespace
  use m_writeoutput
  use m_mover
  implicit none

  integer       :: timestep
contains
  subroutine mainloop()
    implicit none
    ! ADD needs to be changed for restart
    do timestep = 1, final_timestep
      call moveParticles()
      if (modulo(timestep, output_interval) .eq. 0) then
        call writeOutput(INT(timestep / output_interval))
      end if
    end do
  end subroutine mainloop
end module m_mainloop
