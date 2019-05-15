#include "defs.F90"

!--- GLOBAL_NAMESPACE ------------------------------------------!
! To store predefined(!) parameters and functions/interfaces
!   and share them between modules
!...............................................................!

module m_globalnamespace
  #ifdef MPI
    use mpi_f08
  #endif
  implicit none
  integer, parameter     :: dprec = kind(1.0d0)
  integer, parameter     :: sprec = kind(1.0e0)
  integer, parameter     :: UNIT_input = 10, UNIT_output = 20, UNIT_history = 30

  ! plasma parameters
  real :: ppc0, c_omp, sigma, B_norm, unit_ch

  ! simulation parameters
  integer                :: final_timestep
  character(len=STR_MAX) :: input_file_name = 'input',&
                          & output_dir_name = 'output',&
                          & restart_dir_name = 'restart'

  integer       :: nfilter, spec_num
  real          :: spec_min, spec_max

  ! mpi variables
  integer       :: mpi_rank, mpi_size, mpi_statsize
  integer       :: sizex, sizey, sizez
end module m_globalnamespace
