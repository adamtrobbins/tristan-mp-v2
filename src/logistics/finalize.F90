#include "../defs.F90"

module m_finalize
  use m_globalnamespace
  use m_communications
  use m_domain
  use m_particles
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: finalizeCommunications, deallocateArrays
  !...............................................................!
contains
  ! initialize all the necessary functions
  subroutine finalizeAll()
    implicit none
    call deallocateArrays()
    call finalizeCommunications()
  end subroutine finalizeAll

  subroutine deallocateArrays()
    implicit none
    if (allocated(meshblocks)) deallocate(meshblocks)
    this_meshblock%ptr => null()
    if (allocated(sp_)) deallocate(sp_)
    if (allocated(spp_)) deallocate(spp_)
  end subroutine deallocateArrays

  subroutine finalizeCommunications()
    implicit none
    integer :: ierr
    if (mpi_initialized) call MPI_Finalize(ierr)
  end subroutine finalizeCommunications
end module m_finalize
