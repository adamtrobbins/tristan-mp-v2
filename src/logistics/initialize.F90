#include "../defs.F90"

module m_initialize
  use m_globalnamespace
  use m_readinput
  use m_communications
  use m_domain
  use m_particles
  implicit none
  include "mpif.h"

  !--- PRIVATE functions -----------------------------------------!
  private :: readCommandlineArgs, initializeCommunications,&
           & firstRankInitialize, initializeParticles,&
           & distributeMeshblocks
  !...............................................................!
contains
  ! initialize all the necessary functions
  subroutine initializeAll()
    implicit none
    call readCommandlineArgs()
    call initializeCommunications()
    call initializeParticles()
    if (my_rank .eq. 0) call firstRankInitialize()
  end subroutine initializeAll

  subroutine initializeParticles()
    implicit none
    integer                 :: i
    character(len=STR_MAX)  :: var_name

    call getInput('particles', 'nspec', nspec, 2)

    allocate(spp_(nspec))
    allocate(sp_(nspec))

    do i = 1, nspec
      write (var_name, "(A6,I1)") "maxptl", i
      call getInput('particles', var_name, spp_(i)%maxptl_sp)
    end do
  end subroutine initializeParticles

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
      call throwError('ERROR: # of processors is not equal to the number of processors from input')
    end if
    ! ADD possibility to define meshblock distribution in userfile
    call distributeMeshblocks()
  end subroutine initializeCommunications

  subroutine distributeMeshblocks()
    implicit none
    integer :: ir, jr, kr, xr, yr, zr, mx, my, mz
    integer :: rnk
    mx = global_mesh%sx / sizex
    my = global_mesh%sy / sizey
    mz = global_mesh%sz / sizez
    allocate(meshblocks(size0))
    do rnk = 0, size0 - 1
      kr = rnk / (sizex * sizey)
      jr = (rnk - sizex * sizey * kr) / sizex
      ir = rnk - sizex * sizey * kr - sizex * jr
      meshblocks(rnk + 1)%sx = mx
      meshblocks(rnk + 1)%sy = my
      meshblocks(rnk + 1)%sz = mz
      meshblocks(rnk + 1)%x0 = ir * mx + global_mesh%x0
      meshblocks(rnk + 1)%y0 = jr * my + global_mesh%y0
      meshblocks(rnk + 1)%z0 = kr * mz + global_mesh%z0
    end do
    this_meshblock => meshblocks(my_rank + 1)
  end subroutine distributeMeshblocks

  subroutine firstRankInitialize()
    ! create output/restart directories
    !   if does not already exist
    call system('mkdir -p '//trim(output_dir_name))
    call system('mkdir -p '//trim(restart_dir_name))
  end subroutine firstRankInitialize
end module m_initialize
