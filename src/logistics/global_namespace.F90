#include "../defs.F90"

!--- GLOBAL_NAMESPACE ------------------------------------------!
! To store predefined(!) parameters and functions/interfaces
!   and share them between modules
!...............................................................!

module m_globalnamespace
  implicit none
  logical            :: mpi_initialized = .false.
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
    module procedure int4ToMprec
    module procedure int8ToMprec
    module procedure real4ToMprec
    module procedure real8ToMprec
  end interface toMPREC

  private :: int4ToMprec, int8ToMprec, real4ToMprec, real8ToMprec
contains
  real(mprec) function int4ToMprec(inp)
    implicit none
    integer, intent(in) :: inp
    if (mprec .eq. dprec) then
      int4ToMprec = DBLE(inp)
    else
      int4ToMprec = REAL(inp)
    end if
  end function int4ToMprec

  real(mprec) function int8ToMprec(inp)
    implicit none
    integer*8, intent(in) :: inp
    if (mprec .eq. dprec) then
      int8ToMprec = DBLE(inp)
    else
      int8ToMprec = REAL(inp)
    end if
  end function int8ToMprec

  real(mprec) function real4ToMprec(inp)
    implicit none
    real(sprec), intent(in) :: inp
    if (mprec .eq. dprec) then
      real4ToMprec = DBLE(inp)
    else
      real4ToMprec = REAL(inp)
    end if
  end function real4ToMprec

  real(mprec) function real8ToMprec(inp)
    implicit none
    real(dprec), intent(in) :: inp
    if (mprec .eq. dprec) then
      real8ToMprec = DBLE(inp)
    else
      real8ToMprec = REAL(inp)
    end if
  end function real8ToMprec
end module m_globalnamespace
