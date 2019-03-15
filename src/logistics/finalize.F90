#include "../defs.F90"

module m_finalize
  use m_globalnamespace
  use m_aux
  use m_communications
  use m_domain
  use m_particles
  use m_fields
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: finalizeCommunications, deallocateArrays
  !...............................................................!
contains
  subroutine finalizeAll()
    implicit none
    integer :: ierr
    call deallocateArrays()
    call finalizeCommunications()
    call printReport((mpi_rank .eq. 0), "finalizeAll()" // TAB // TAB // TAB // TAB // "[OK]")
  end subroutine finalizeAll

  subroutine deallocateArrays()
    implicit none
    integer :: i, j, k, ierr

    ! dealloc meshblocks
    if (allocated(meshblocks)) deallocate(meshblocks)
    nullify(this_meshblock%ptr)

    ! dealloc particle arrays
    if (allocated(sp_)) deallocate(sp_)
    if (allocated(spp_)) deallocate(spp_)

    ! dealloc exchange arrays
    do i = -1, 1
      do j = -1, 1
        do k = -1, 1
          if (allocated(enroute_bot%get(i,j,k)%send_enroute)) deallocate(enroute_bot%get(i,j,k)%send_enroute)
        end do
      end do
    end do
    if (allocated(recv_enroute)) deallocate(recv_enroute)
    call MPI_TYPE_FREE(myMPI_ENROUTE, ierr)

    ! dealloc field arrays
    if (allocated(ex)) deallocate(ex)
    if (allocated(ey)) deallocate(ey)
    if (allocated(ez)) deallocate(ez)
    if (allocated(bx)) deallocate(bx)
    if (allocated(by)) deallocate(by)
    if (allocated(bz)) deallocate(bz)
  end subroutine deallocateArrays

  subroutine finalizeCommunications()
    implicit none
    integer :: ierr
    if (mpi_initialized) call MPI_FINALIZE(ierr)
  end subroutine finalizeCommunications
end module m_finalize
