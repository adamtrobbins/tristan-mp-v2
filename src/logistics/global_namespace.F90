#include "../defs.F90"

!--- GLOBAL_NAMESPACE ------------------------------------------!
! To store predefined(!) parameters and functions/interfaces
!   and share them between modules
!...............................................................!

module m_globalnamespace
  implicit none
  integer, parameter :: dprec = kind(1.0d0)
  integer, parameter :: sprec = kind(1.0e0)
  ! main precision: `mprec`
  #ifdef DPREC
    integer, parameter :: mprec = dprec
  #else
    integer, parameter :: mprec = sprec
  #endif

  character(len=STR_MAX) :: input_file_name = 'input',&
                          & output_dir_name = 'output',&
                          & restart_dir_name = 'restart'

  integer, parameter :: UNIT_input = 10

  interface toMPREC
    module procedure intToMprec
    module procedure sprecToMprec
    module procedure dprecToMprec
  end interface toMPREC

  private :: intToMprec, sprecToMprec, dprecToMprec

contains
  real(mprec) function intToMprec(inp)
    implicit none
    integer, intent(in) :: inp
    if (mprec .eq. dprec) then
      intToMprec = DBLE(inp)
    else
      intToMprec = REAL(inp)
    end if
  end function intToMprec

  real(mprec) function sprecToMprec(inp)
    implicit none
    real(sprec), intent(in) :: inp
    if (mprec .eq. dprec) then
      sprecToMprec = DBLE(inp)
    else
      sprecToMprec = REAL(inp)
    end if
  end function sprecToMprec

  real(mprec) function dprecToMprec(inp)
    implicit none
    real(dprec), intent(in) :: inp
    if (mprec .eq. dprec) then
      dprecToMprec = DBLE(inp)
    else
      dprecToMprec = REAL(inp)
    end if
  end function dprecToMprec

  subroutine throwError(msg)
    character(len=*), intent(in)  :: msg
    print *, msg
    stop 'TERMINATING EXECUTION'
  end subroutine
end module m_globalnamespace
