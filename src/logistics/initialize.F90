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

  !--- PRIVATE functions -----------------------------------------!
  private :: initializeCommunications, initializeOutput, &
           & firstRankInitialize, initializeParticles,&
           & distributeMeshblocks, initializeDomain,&
           & rnkToInd, indToRnk, assignNeighbor,&
           & allocateParticles, initializeSimulation,&
           & initializePrtlExchange
  !...............................................................!
contains
  ! initialize all the necessary functions
  subroutine initializeAll()
    implicit none
    call readCommandlineArgs()
    call initializeOutput()
    ! ADD possibility to define output function in userfile
    ! ADD hst file?
    call initializeDomain()
    call initializeCommunications()
    ! ADD possibility to define meshblock distribution in userfile
    call distributeMeshblocks()
    call initializeParticles()
    call initializeSimulation()
    call initializePrtlExchange()

    call initializeRandomSeed(mpi_rank)

    if (mpi_rank .eq. 0) call firstRankInitialize()
    call userInitialize()

    call print_diag((mpi_rank .eq. 0), "initializeAll()" // TAB // TAB // TAB // "[OK]")
  end subroutine initializeAll

  subroutine initializeOutput()
    implicit none
    call getInput('output', 'stride', output_stride, 10)
    call getInput('output', 'interval', output_interval, 10)
  end subroutine initializeOutput

  subroutine initializePrtlExchange()
    implicit none
    integer            :: buffsize, buffsize_x, buffsize_y
    integer            :: buffsize_xy
    integer            :: buffsize_z, buffsize_xz, buffsize_yz
    integer            :: buffsize_xyz
    integer            :: multiplier, ierr, ind1, ind2, ind3

    integer, dimension(0:2) :: oldtypes, blockcounts, offsets
    integer                 :: extent_int2, extent_real

    call getInput('particles', 'ppc0', multiplier)
    multiplier = max(multiplier, 1) * 1000
    ! FIX this might change over time (due to load balancing)
    buffsize = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * multiplier
    buffsize_x = this_meshblock%ptr%sy * this_meshblock%ptr%sz * multiplier
    buffsize_y = this_meshblock%ptr%sx * this_meshblock%ptr%sz * multiplier
    buffsize_xy = this_meshblock%ptr%sz * multiplier
    #ifdef threeD
      buffsize_z = this_meshblock%ptr%sx * this_meshblock%ptr%sy * multiplier
      buffsize_xz = this_meshblock%ptr%sz * multiplier
      buffsize_yz = this_meshblock%ptr%sx * multiplier
      buffsize_xyz = multiplier
    #else
      buffsize_z = 0; buffsize_xz = 0; buffsize_yz = 0; buffsize_xyz = 0
    #endif

    allocate(recv_enroute(buffsize))

    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifndef threeD
            if (ind3 .ne. 0) cycle
          #endif
          if ((ind2 .eq. 0) .and. (ind3 .eq. 0)) then
            buffsize = buffsize_x
          else if ((ind1 .eq. 0) .and. (ind3 .eq. 0)) then
            buffsize = buffsize_y
          else if ((ind1 .eq. 0) .and. (ind2 .eq. 0)) then
            buffsize = buffsize_z
          else if (ind3 .eq. 0) then
            buffsize = buffsize_xy
          else if (ind2 .eq. 0) then
            buffsize = buffsize_xz
          else if (ind1 .eq. 0) then
            buffsize = buffsize_yz
          else
            buffsize = buffsize_xyz
          end if
          allocate(enroute_bot%get(ind1, ind2, ind3)%send_enroute(buffsize))
        end do
      end do
    end do

    ! new type for myMPI_ENROUTE
    !   BY DEFAULT:
    !     # of blockcounts = 3:
    !       3 x integer2  [xi, yi, zi]
    !       6 x real      [dx, dy, dz, u, v, w]
    !       2 x integer   [ind, proc]
    call MPI_TYPE_EXTENT(MPI_INTEGER2, extent_int2, ierr)
    call MPI_TYPE_EXTENT(MPI_REAL, extent_real, ierr)
    blockcounts(0) = 3
    oldtypes(0) = MPI_INTEGER2
    blockcounts(1) = 6
    oldtypes(1) = MPI_REAL
    blockcounts(2) = 2
    oldtypes(2) = MPI_INTEGER
    offsets(0) = 0
    offsets(1) = blockcounts(0) * extent_int2 + offsets(0)
    offsets(2) = blockcounts(1) * extent_real + offsets(1)
  	call MPI_TYPE_STRUCT(3, blockcounts, offsets, oldtypes, myMPI_ENROUTE, ierr)
  	call MPI_TYPE_COMMIT(myMPI_ENROUTE, ierr)
  end subroutine initializePrtlExchange

  subroutine initializeSimulation()
    implicit none
    call getInput('time', 'last', final_timestep, 1000)
  end subroutine initializeSimulation

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
      spp_(i)%cntr_sp = 0
    end do
  end subroutine initializeParticles

  subroutine allocateParticles(prt, sz)
    implicit none
    type(species), intent(inout)    :: prt
    integer, intent(in)             :: sz
    if (allocated(prt%xi)) deallocate(prt%xi)
    if (allocated(prt%yi)) deallocate(prt%yi)
    if (allocated(prt%zi)) deallocate(prt%zi)
    if (allocated(prt%dx)) deallocate(prt%dx)
    if (allocated(prt%dy)) deallocate(prt%dy)
    if (allocated(prt%dz)) deallocate(prt%dz)
    if (allocated(prt%u)) deallocate(prt%u)
    if (allocated(prt%v)) deallocate(prt%v)
    if (allocated(prt%w)) deallocate(prt%w)
    allocate(prt%xi(sz)); allocate(prt%yi(sz)); allocate(prt%zi(sz))
    allocate(prt%dx(sz)); allocate(prt%dy(sz)); allocate(prt%dz(sz))
    allocate(prt%u(sz)); allocate(prt%v(sz)); allocate(prt%w(sz))
    allocate(prt%ind(sz)); allocate(prt%proc(sz))
  end subroutine allocateParticles

  subroutine initializeCommunications()
    implicit none
    integer :: ierr
    call MPI_Init(ierr)
    ! ADD if statement here
    mpi_initialized = .true.
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
    mpi_statsize = MPI_STATUS_SIZE
    if (mpi_size .ne. sizex * sizey * sizez) then
      call throwError('ERROR: # of processors is not equal to the number of processors from input')
    end if
  end subroutine initializeCommunications

  subroutine initializeDomain()
    implicit none
    call getInput('node_configuration', 'sizex', sizex, 1)
    call getInput('node_configuration', 'sizey', sizey, 1)
    #ifdef threeD
      call getInput('node_configuration', 'sizez', sizez, 1)
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
    integer               :: rnk, ind1, ind2, ind3
    m(1) = global_mesh%sx / sizex
    m(2) = global_mesh%sy / sizey
    m(3) = global_mesh%sz / sizez
    allocate(meshblocks(mpi_size))
    do rnk = 0, mpi_size - 1
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
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            call assignNeighbor(rnk, (/ ind1, ind2, ind3/))
          end do
        end do
      end do
    end do
    this_meshblock%ptr => meshblocks(mpi_rank + 1)
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
    if ((rnk .lt. 0) .or. (rnk .ge. mpi_size)) then
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
    call system('mkdir -p ' // trim(output_dir_name))
    call system('mkdir -p ' // trim(restart_dir_name))
  end subroutine firstRankInitialize
end module m_initialize
