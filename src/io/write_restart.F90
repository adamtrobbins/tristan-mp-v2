#include "../defs.F90"

module m_writerestart
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  implicit none
contains
  subroutine writeRestart(timestep)
    implicit none
    integer, intent(in) :: timestep
  end subroutine
end module m_writerestart
