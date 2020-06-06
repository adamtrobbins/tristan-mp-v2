#include "../defs.F90"

module m_qedphysics
#ifdef QED

  use m_aux
  use m_globalnamespace
  #ifdef BWPAIRPRODUCTION
    use m_bwpairproduction
  #endif
  #ifdef COMPTONSCATTERING
    use m_compton
  #endif
  implicit none

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine QEDstep(timestep)
    implicit none
    integer, intent(in) :: timestep

    #ifdef COMPTONSCATTERING
      if (modulo(timestep, Compton_interval) .eq. 0) then
        call comptonScattering()
      end if
    #endif

    #ifdef BWPAIRPRODUCTION
      if (modulo(timestep, BW_interval) .eq. 0) then
        call bwPairProduction()
      end if
    #endif

    call printDiag((mpi_rank .eq. 0), "QEDstep()", .true.)
  end subroutine QEDstep

#endif
end module m_qedphysics
