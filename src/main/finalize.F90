#include "../defs.F90"

module m_finalize
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_particles
  use m_fieldlogistics, only: deallocateFields, deallocateFieldBackups
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: finalizeCommunications, deallocateArrays
  !...............................................................!
contains
  subroutine finalizeAll()
    implicit none
    integer :: ierr
    call deallocateArrays()
      call printDiag("deallocateArrays()", 1)
      call finalizeCommunications()
      call printDiag("finalizeCommunications()", 1)

    call printDiag("FinalizeAll()", 0)
  end subroutine finalizeAll

  subroutine deallocateArrays()
    implicit none
    integer :: i, j, k, ierr

    ! dealloc meshblocks
    if (allocated(meshblocks)) deallocate(meshblocks)
    if (allocated(new_meshblocks)) deallocate(new_meshblocks)
    nullify(this_meshblock%ptr)
    ! SLB/ALB arrays
    if (allocated(lb_load_glob)) deallocate(lb_load_glob)
    if (allocated(lb_group_x0)) deallocate(lb_group_x0)
    if (allocated(lb_group_x1)) deallocate(lb_group_x1)
    if (allocated(lb_group_y0)) deallocate(lb_group_y0)
    if (allocated(lb_group_y1)) deallocate(lb_group_y1)
    if (allocated(lb_group_z0)) deallocate(lb_group_z0)
    if (allocated(lb_group_z1)) deallocate(lb_group_z1)

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

    call deallocateFields()
    call deallocateFieldBackups()
  end subroutine deallocateArrays

  subroutine finalizeCommunications()
    implicit none
    integer :: ierr
    call MPI_FINALIZE(ierr)
  end subroutine finalizeCommunications
end module m_finalize
