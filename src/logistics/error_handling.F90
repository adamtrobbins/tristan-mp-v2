#include "../defs.F90"

module m_errors
  use m_globalnamespace
  use m_finalize
  use m_communications
  implicit none
contains
  subroutine throwError(msg)
    character(len=*), intent(in)  :: msg
    print *, msg
    call finalizeAll()
    stop 'TERMINATING EXECUTION'
  end subroutine
end module m_errors
