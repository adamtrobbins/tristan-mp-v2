#include "../defs.F90"

module m_particles
  use m_globalnamespace
  implicit none

  type :: particle_tile
    ! DEP_PRT [particle-dependent]
    integer                                     :: npart_sp, maxptl_sp
    ! species index for a given tile
    integer                                     :: spec
    ! tile boundaries in local coordinates
    integer                                     :: x1, x2, y1, y2, z1, z2
    integer(kind=2), allocatable, dimension(:)  :: xi, yi, zi
    real, allocatable, dimension(:)             :: dx, dy, dz
    real, allocatable, dimension(:)             :: u, v, w
    real, allocatable, dimension(:)             :: weight
    integer, allocatable, dimension(:)          :: ind, proc
    !dir$ attributes align: 64 :: xi, yi, zi, dx, dy, dz, u, v, w, weight, ind, proc
    ! > `proc < 0` means the particle will be deleted once the `clearGhostParticles()` is called
  end type particle_tile

  type :: particle_species
    integer     :: cntr_sp
    real        :: m_sp, ch_sp
    ! sizes and boundaries of the tiles in each direction
    integer     :: tile_sx, tile_sy, tile_sz
    ! numbers of the tiles in each direction
    integer     :: tile_nx, tile_ny, tile_nz
    type (particle_tile), allocatable, dimension(:,:,:) :: prtl_tile

    ! extra physics properties
    #ifdef RADIATION
      ! `true/false` - either apply cooling to species or not
      logical     :: cool_sp
    #endif

    #ifdef BWPAIRPRODUCTION
      ! `0` means species does not participate in BW process
      ! `1` and `2` would be separate BW groups
      integer     :: bw_sp
    #endif

    #ifdef DOWNSAMPLING
    ! `true/false` - either downsample species or not
      logical     :: dwn_sp
    #endif
  end type particle_species

  ! particle types for exchange between processors />
  type :: prtl_enroute
    ! DEP_PRT [particle-dependent]
    integer(kind=2)   :: xi, yi, zi
    real              :: dx, dy, dz
    real              :: u, v, w
    real              :: weight
    integer           :: ind, proc
  end type prtl_enroute

  type :: enroute_array
    type(prtl_enroute), allocatable    :: send_enroute(:)
    integer                            :: cnt_send
  end type enroute_array

  type :: enroute_handler
    type(enroute_array), dimension(-1:1,-1:1,-1:1)   :: get
  end type enroute_handler
  ! </ particle types for exchange between processors

  ! main container for particles
  type(particle_species), target, allocatable   :: species(:)
  ! number of species
  integer                                       :: nspec

  type(prtl_enroute), allocatable, dimension(:)    :: recv_enroute
  type(enroute_handler)                            :: enroute_bot

  #ifdef MPI08
    type(MPI_DATATYPE)                             :: myMPI_ENROUTE
  #endif

  #ifdef MPI
    integer                                        :: myMPI_ENROUTE
  #endif

end module m_particles
