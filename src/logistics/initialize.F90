#include "../defs.F90"

module m_initialize
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_particlelogistics
  use m_fields
  use m_userfile
  use m_helpers
  use m_errors
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: initializeCommunications, initializeOutput, &
           & firstRankInitialize, initializeParticles,&
           & distributeMeshblocks, initializeDomain,&
           & initializePrtlExchange, initializeFields,&
           & assignNeighbor, allocateParticles,&
           & initializeSimulation, checkEverything
  !...............................................................!
contains
  ! initialize all the necessary arrays and variables
  subroutine initializeAll()
    implicit none
    call readCommandlineArgs()
    call initializeSimulation()
    call initializeOutput()
    ! ADD possibility to define output function in userfile
    ! ADD hst file?
    call initializeDomain()
    call initializeCommunications()
    ! ADD possibility to define meshblock distribution in userfile
    call distributeMeshblocks()

    call initializeFields()
    call initializeParticles()
    call initializePrtlExchange()

    call initializeRandomSeed(mpi_rank)

    if (mpi_rank .eq. 0) call firstRankInitialize()
    call userInitialize()

    ! check everything before moving forward
    call checkEverything()

    call printReport((mpi_rank .eq. 0), "initializeAll()" // TAB // TAB // TAB // "[OK]")
  end subroutine initializeAll

  subroutine initializeOutput()
    implicit none
    call getInput('output', 'stride', output_stride, 10)
    call getInput('output', 'interval', output_interval, 10)
    call getInput('output', 'istep', output_istep, 4)

    call getInput('output', 'spec_min', spec_min, 1e-2)
    call getInput('output', 'spec_max', spec_max, 1e2)
    call getInput('output', 'spec_num', spec_num, 100)
    spec_min = log(spec_min)
    spec_max = log(spec_max)
  end subroutine initializeOutput

  subroutine initializeSimulation()
    implicit none
    call getInput('time', 'last', final_timestep, 1000)
    call getInput('algorithm', 'nfilter', nfilter, 2)
  end subroutine initializeSimulation

  subroutine initializeParticles()
    implicit none
    integer                 :: s, ti, tj, tk
    character(len=STR_MAX)  :: var_name
		integer        					:: maxptl_

    call getInput('plasma', 'ppc0', ppc0)
    call getInput('plasma', 'sigma', sigma)
    call getInput('plasma', 'c_omp', c_omp)
    B_norm = CC**2 * sqrt(sigma) / c_omp
    unit_ch = CC**2 / (ppc0 * c_omp**2)

    call getInput('particles', 'nspec', nspec, 2)

    allocate(species(nspec))
		do s = 1, nspec
			call getInput('grid', 'tileX', species(s)%tile_sx)
			call getInput('grid', 'tileY', species(s)%tile_sy)
			call getInput('grid', 'tileZ', species(s)%tile_sz)
			#ifndef threeD
				species(s)%tile_sz = 1
			#endif
			species(s)%tile_nx = ceiling(real(this_meshblock%ptr%sx) / real(species(s)%tile_sx))
			species(s)%tile_ny = ceiling(real(this_meshblock%ptr%sy) / real(species(s)%tile_sy))
			species(s)%tile_nz = ceiling(real(this_meshblock%ptr%sz) / real(species(s)%tile_sz))
			allocate(species(s)%prtl_tile(species(s)%tile_nx,&
			                            & species(s)%tile_ny,&
																	& species(s)%tile_nz))
		end do

    do s = 1, nspec
      write (var_name, "(A6,I1)") "maxptl", s
      call getInput('particles', var_name, maxptl_)
      write (var_name, "(A1,I1)") "m", s
      call getInput('particles', var_name, species(s)%m_sp)
      write (var_name, "(A2,I1)") "ch", s
      call getInput('particles', var_name, species(s)%ch_sp)
			do ti = 1, species(s)%tile_nx
				do tj = 1, species(s)%tile_ny
					do tk = 1, species(s)%tile_nz
						species(s)%prtl_tile(ti, tj, tk)%maxptl_sp = max(maxptl_ / &
															& (species(s)%tile_nx * species(s)%tile_ny * species(s)%tile_nz),&
															& max(INT(ppc0), 1) * species(s)%tile_sx * species(s)%tile_sy * species(s)%tile_sz)
						species(s)%prtl_tile(ti, tj, tk)%npart_sp = 0
						call allocateParticles(species(s)%prtl_tile(ti, tj, tk),&
																 & species(s)%prtl_tile(ti, tj, tk)%maxptl_sp)
					end do
				end do
			end do
      species(s)%cntr_sp = 0
    end do
  end subroutine initializeParticles

  subroutine initializePrtlExchange()
    implicit none
    integer            :: buffsize, buffsize_x, buffsize_y
    integer            :: buffsize_xy
    integer            :: buffsize_z, buffsize_xz, buffsize_yz
    integer            :: buffsize_xyz
    integer            :: multiplier, ierr, ind1, ind2, ind3

    type(MPI_DATATYPE), dimension(0:2)              :: oldtypes
    integer, dimension(0:2)                         :: blockcounts
    integer(kind=MPI_ADDRESS_KIND), dimension(0:2)  :: offsets
    integer(kind=MPI_COUNT_KIND)                    :: extent_int2, extent_real, lb

    multiplier = max(INT(ppc0), 1) * 1000
    ! FIX this might change over time (due to load balancing)
    buffsize_x = this_meshblock%ptr%sy * this_meshblock%ptr%sz * multiplier
    buffsize_y = this_meshblock%ptr%sx * this_meshblock%ptr%sz * multiplier
    buffsize_xy = this_meshblock%ptr%sz * multiplier
    #ifdef threeD
			buffsize = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * multiplier

      buffsize_z = this_meshblock%ptr%sx * this_meshblock%ptr%sy * multiplier
      buffsize_xz = this_meshblock%ptr%sz * multiplier
      buffsize_yz = this_meshblock%ptr%sx * multiplier
      buffsize_xyz = multiplier
    #else
			buffsize = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz) * multiplier

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
    call MPI_TYPE_GET_EXTENT(MPI_INTEGER2, lb, extent_int2, ierr)
    call MPI_TYPE_GET_EXTENT(MPI_REAL, lb, extent_real, ierr)
    blockcounts(0) = 3
    oldtypes(0) = MPI_INTEGER2
    blockcounts(1) = 6
    oldtypes(1) = MPI_REAL
    blockcounts(2) = 2
    oldtypes(2) = MPI_INTEGER
    offsets(0) = 0
    offsets(1) = blockcounts(0) * extent_int2 + offsets(0)
    offsets(2) = blockcounts(1) * extent_real + offsets(1)
  	call MPI_TYPE_CREATE_STRUCT(3, blockcounts, offsets, oldtypes, myMPI_ENROUTE, ierr)
  	call MPI_TYPE_COMMIT(myMPI_ENROUTE, ierr)
  end subroutine initializePrtlExchange

  subroutine initializeFields()
    implicit none
    if (allocated(ex)) deallocate(ex)
    if (allocated(ey)) deallocate(ey)
    if (allocated(ez)) deallocate(ez)
    if (allocated(bx)) deallocate(bx)
    if (allocated(by)) deallocate(by)
    if (allocated(bz)) deallocate(bz)
    if (allocated(jx)) deallocate(jx)
    if (allocated(jy)) deallocate(jy)
    if (allocated(jz)) deallocate(jz)
    allocate(ex(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(ey(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(ez(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(by(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(bz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(jx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(jy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))
    allocate(jz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
              & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
              & fldBoundZ))

    ! exchange fields
    ! 20 = max # of fields sent in each direction
    #ifndef threeD
      sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz) * NGHOST * 20
      sendrecv_buffsz = sendrecv_offsetsz * 10
      ! 8 (~10) directions to send/recv in 2D
    #else
      sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * NGHOST * 20
      sendrecv_buffsz = sendrecv_offsetsz * 30
      ! 26 (~30) directions to send/recv in 3D
    #endif

    if (allocated(send_fld)) deallocate(send_fld)
    allocate(send_fld(sendrecv_buffsz))
    if (allocated(recv_fld)) deallocate(recv_fld)
    allocate(recv_fld(sendrecv_offsetsz))

    ! output fields
    if (allocated(scalar_array)) deallocate(scalar_array)
    allocate(scalar_array(0:this_meshblock%ptr%sx - 1,&
                        & 0:this_meshblock%ptr%sy - 1,&
                        & 0:this_meshblock%ptr%sz - 1))
  end subroutine initializeFields

  subroutine initializeCommunications()
    implicit none
    integer :: ierr
    call MPI_Init(ierr)
    ! ADD `ifdef MPI`-statement here
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
    global_mesh%x0 = 0
    global_mesh%y0 = 0
    global_mesh%z0 = 0
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
    integer               :: rnk, ind1, ind2, ind3, cntr
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

    cntr = 0
    ! find the number of neighbors
    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifndef threeD
            if (ind3 .ne. 0) cycle
          #endif
          if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
          cntr = cntr + 1
        end do
      end do
    end do
    sendrecv_neighbors = cntr
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
  end subroutine assignNeighbor

  subroutine firstRankInitialize()
    ! create output/restart directories
    !   if does not already exist
    call system('mkdir -p ' // trim(output_dir_name))
    call system('mkdir -p ' // trim(restart_dir_name))
  end subroutine firstRankInitialize

  subroutine checkEverything()
    implicit none
    ! check that the domain size is larger than the number of ghost zones
    #ifndef threeD
      if ((this_meshblock%ptr%sx .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sy .lt. NGHOST)) then
        call throwError('ERROR: ghost zones overflow the domain size in ' // trim(STR(mpi_rank)))
      end if
    #else
      if ((this_meshblock%ptr%sx .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sy .lt. NGHOST) .or.&
        & (this_meshblock%ptr%sz .lt. NGHOST)) then
        call throwError('ERROR: ghost zones overflow the domain size in ' // trim(STR(mpi_rank)))
      end if
    #endif
  end subroutine checkEverything
end module m_initialize
