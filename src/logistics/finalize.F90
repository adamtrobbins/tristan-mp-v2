#include "../defs.F90"

module m_finalize
  use m_globalnamespace
  use m_aux
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

    ! dealloc particle species
    if (allocated(species)) deallocate(species)

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

    ! dealloc field exchange
    if (allocated(send_fld)) deallocate(send_fld)
    if (allocated(recv_fld)) deallocate(recv_fld)

    ! dealloc field output
    if (allocated(scalar_array)) deallocate(scalar_array)
  end subroutine deallocateArrays

  subroutine finalizeCommunications()
    implicit none
    integer :: ierr
    #ifdef MPI
      call MPI_FINALIZE(ierr)
    #endif
  end subroutine finalizeCommunications
end module m_finalize
