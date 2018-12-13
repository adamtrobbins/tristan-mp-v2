#include "../defs.F90"

module m_communications
  use m_globalnamespace
  use m_readinput
  use m_domain

  implicit none

  include "mpif.h"

  ! mpi variables
  integer           :: my_rank, size0, statsize
  integer, private  :: ierr

  integer           :: sizex, sizey, sizez
contains
  subroutine initializeCommunications()
    implicit none
    integer ierr

    call getInput('node_configuration', 'sizex', sizex, 1)
    call getInput('node_configuration', 'sizey', sizey, 1)
    #ifdef threeD
      sizez = getInput('node_configuration', 'sizez', sizez, 1)
    #else
      sizez = 1
    #endif

    allocate(meshblocks(size0))
    global_mesh%x0 = 1
    global_mesh%y0 = 1
    global_mesh%z0 = 1
    call getInput('grid', 'mx0', global_mesh%sx, 1)
    call getInput('grid', 'my0', global_mesh%sy, 1)
    #ifdef threeD
      call getInput('grid', 'mz0', global_mesh%sz, 1)
    #else
      global_mesh%sz = 1
    #endif

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_Comm_world, my_rank, ierr)
    call MPI_Comm_size(MPI_Comm_world, size0, ierr)
    statsize = MPI_STATUS_SIZE

    if (size0 .ne. sizex * sizey * sizez) then

    end if
    print *, my_rank, size0, ierr, statsize
  end subroutine initializeCommunications
end module m_communications
