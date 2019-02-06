#include "../defs.F90"

!--- GLOBAL_NAMESPACE ------------------------------------------!
! To store predefined(!) parameters and functions/interfaces
!   and share them between modules
!...............................................................!

module m_globalnamespace
  implicit none
  logical                :: mpi_initialized = .false.
  integer, parameter     :: dprec = kind(1.0d0)
  integer, parameter     :: sprec = kind(1.0e0)
  integer, parameter     :: UNIT_input = 10

  ! simulation parameters
  integer                :: final_timestep
  character(len=STR_MAX) :: input_file_name = 'input',&
                          & output_dir_name = 'output',&
                          & restart_dir_name = 'restart'
end module m_globalnamespace
