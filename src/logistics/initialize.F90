#include "../defs.F90"

module m_initialize
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_communications
  use m_domain
  use m_particles
  use m_userfile
  implicit none
  include "mpif.h"

  !--- PRIVATE functions -----------------------------------------!
  private :: initializeCommunications,&
           & firstRankInitialize, initializeParticles,&
           & distributeMeshblocks, initializeDomain,&
           & rnkToInd, indToRnk, assignNeighbor,&
           & allocateParticles
  !...............................................................!
contains
  ! initialize all the necessary functions
  subroutine initializeAll()
    implicit none
    call readCommandlineArgs()
    call initializeDomain()
    call initializeCommunications()
    ! ADD possibility to define meshblock distribution in userfile
    call distributeMeshblocks()
    call initializeParticles()

    call initializeRandomSeed(my_rank)

    if (my_rank .eq. 0) call firstRankInitialize()

    call userInitialize()
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
      write (var_name, "(A1,I1)") "m", i
      call getInput('particles', var_name, spp_(i)%m_sp)
      write (var_name, "(A2,I1)") "ch", i
      call getInput('particles', var_name, spp_(i)%ch_sp)
      call allocateParticles(sp_(i), spp_(i)%maxptl_sp)
      spp_(i)%npart_sp = 0
    end do
  end subroutine initializeParticles

  subroutine allocateParticles(prt, sz)
    implicit none
    type(species), intent(inout)    :: prt
    integer, intent(in)             :: sz
    if (allocated(prt%x)) deallocate(prt%x)
    if (allocated(prt%y)) deallocate(prt%y)
    if (allocated(prt%z)) deallocate(prt%z)
    if (allocated(prt%u)) deallocate(prt%u)
    if (allocated(prt%v)) deallocate(prt%v)
    if (allocated(prt%w)) deallocate(prt%w)
    allocate(prt%x(sz)); allocate(prt%y(sz)); allocate(prt%z(sz))
    allocate(prt%u(sz)); allocate(prt%v(sz)); allocate(prt%w(sz))
    allocate(prt%ind(sz)); allocate(prt%proc(sz))
  end subroutine allocateParticles

  subroutine initializeCommunications()
    implicit none
    integer ierr
    call MPI_Init(ierr)
    ! ADD if statement here
    mpi_initialized = .true.
    call MPI_Comm_rank(MPI_Comm_world, my_rank, ierr)
    call MPI_Comm_size(MPI_Comm_world, size0, ierr)
    statsize = MPI_STATUS_SIZE
    if (size0 .ne. sizex * sizey * sizez) then
      call throwError('ERROR: # of processors is not equal to the number of processors from input')
    end if
  end subroutine initializeCommunications

  subroutine initializeDomain()
    implicit none
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
    call getInput('grid', 'boundary_x', boundary_x, 1)
    call getInput('grid', 'boundary_y', boundary_y, 1)
    #ifdef threeD
      call getInput('grid', 'boundary_z', boundary_z, 1)
    #else
      boundary_z = 0
    #endif
  end subroutine initializeDomain

  subroutine distributeMeshblocks()
    implicit none
    integer, dimension(3) :: ind, m
    integer               :: rnk, rnk2
    m(1) = global_mesh%sx / sizex
    m(2) = global_mesh%sy / sizey
    m(3) = global_mesh%sz / sizez
    allocate(meshblocks(size0))
    this_meshblock%ptr => meshblocks(my_rank + 1)
    do rnk = 0, size0 - 1
      ind = rnkToInd(rnk)
      meshblocks(rnk + 1)%rnk = rnk
      ! find sizes and corner coords
      meshblocks(rnk + 1)%sx = m(1)
      meshblocks(rnk + 1)%sy = m(2)
      meshblocks(rnk + 1)%sz = m(3)
      meshblocks(rnk + 1)%x0 = ind(1) * m(1) + global_mesh%x0
      meshblocks(rnk + 1)%y0 = ind(2) * m(2) + global_mesh%y0
      meshblocks(rnk + 1)%z0 = ind(3) * m(3) + global_mesh%z0

      ! assign neighbors
      call assignNeighbor(rnk, (/ 0, 0, 0/))
      call assignNeighbor(rnk, (/ 0, 0,-1/))
      call assignNeighbor(rnk, (/ 0,-1, 0/))
      call assignNeighbor(rnk, (/-1, 0, 0/))
      call assignNeighbor(rnk, (/ 0, 0,+1/))
      call assignNeighbor(rnk, (/ 0,+1, 0/))
      call assignNeighbor(rnk, (/+1, 0, 0/))
      call assignNeighbor(rnk, (/ 0,-1,-1/))
      call assignNeighbor(rnk, (/ 0,-1,+1/))
      call assignNeighbor(rnk, (/ 0,+1,-1/))
      call assignNeighbor(rnk, (/ 0,+1,+1/))
      call assignNeighbor(rnk, (/-1, 0,-1/))
      call assignNeighbor(rnk, (/-1, 0,+1/))
      call assignNeighbor(rnk, (/+1, 0,-1/))
      call assignNeighbor(rnk, (/+1, 0,+1/))
      call assignNeighbor(rnk, (/-1,-1, 0/))
      call assignNeighbor(rnk, (/-1,+1, 0/))
      call assignNeighbor(rnk, (/+1,-1, 0/))
      call assignNeighbor(rnk, (/+1,+1, 0/))
      call assignNeighbor(rnk, (/-1,-1,-1/))
      call assignNeighbor(rnk, (/-1,-1,+1/))
      call assignNeighbor(rnk, (/-1,+1,-1/))
      call assignNeighbor(rnk, (/-1,+1,+1/))
      call assignNeighbor(rnk, (/+1,-1,-1/))
      call assignNeighbor(rnk, (/+1,-1,+1/))
      call assignNeighbor(rnk, (/+1,+1,-1/))
      call assignNeighbor(rnk, (/+1,+1,+1/))
    end do
  end subroutine distributeMeshblocks

  subroutine assignNeighbor(rnk, inds1)
    implicit none
    integer, intent(in)       :: rnk, inds1(3)
    integer                   :: rnk2, inds0(3)
    inds0 = rnkToInd(rnk)
    rnk2 = indToRnk([inds0(1) + inds1(1), inds0(2) + inds1(2), inds0(3) + inds1(3)])
    if (rnk2 .eq. -1) then
      meshblocks(rnk + 1)%neighbor(inds1(1), inds1(2), inds1(3))%ptr => null()
    else
      meshblocks(rnk + 1)%neighbor(inds1(1), inds1(2), inds1(3))%ptr => meshblocks(rnk2 + 1)
    end if
  end subroutine

  function rnkToInd(rnk)
    implicit none
    integer, intent(in)   :: rnk
    integer, dimension(3) :: rnkToInd
    if ((rnk .lt. 0) .or. (rnk .ge. size0)) then
      rnkToInd = (/-1, -1, -1/)
    else
      rnkToInd(3) = rnk / (sizex * sizey)
      rnkToInd(2) = (rnk - sizex * sizey * rnkToInd(3)) / sizex
      rnkToInd(1) = rnk - sizex * sizey * rnkToInd(3) - sizex * rnkToInd(2)
    end if
  end function rnkToInd

  function indToRnk(ind)
    implicit none
    integer, intent(in)                 :: ind(3)
    integer                             :: ind_(3), indToRnk
    ind_ = ind
    if (boundary_x .eq. 1) ind_(1) = modulo(ind_(1), sizex)
    if (boundary_y .eq. 1) ind_(2) = modulo(ind_(2), sizey)
    if (boundary_z .eq. 1) ind_(3) = modulo(ind_(3), sizez)
    if ((ind_(1) .lt. 0) .or. (ind_(1) .ge. sizex) .or.&
      & (ind_(2) .lt. 0) .or. (ind_(2) .ge. sizey) .or.&
      & (ind_(3) .lt. 0) .or. (ind_(3) .ge. sizez)) then
      indToRnk = -1
    else
      indToRnk = ind_(3) * sizex * sizey + ind_(2) * sizex + ind_(1)
    end if
  end function indToRnk

  subroutine firstRankInitialize()
    ! create output/restart directories
    !   if does not already exist
    call system('mkdir -p '//trim(output_dir_name))
    call system('mkdir -p '//trim(restart_dir_name))
  end subroutine firstRankInitialize
end module m_initialize
