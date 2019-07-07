#include "../defs.F90"

module m_particles
  use m_globalnamespace
  implicit none

  type :: particle_tile
    integer                                     :: npart_sp, maxptl_sp
    integer(kind=2), allocatable, dimension(:)  :: xi, yi, zi
    real, allocatable, dimension(:)             :: dx, dy, dz
    real, allocatable, dimension(:)             :: u, v, w
    integer, allocatable, dimension(:)          :: ind, proc
    !dir$ attributes align: 64 :: xi, yi, zi, dx, dy, dz, u, v, w, ind, proc
  end type particle_tile

  type :: particle_species
    integer     :: cntr_sp
    real        :: m_sp, ch_sp
    ! sizes of the tiles in each direction
    integer     :: tile_sx, tile_sy, tile_sz
    ! numbers of the tiles in each direction
    integer     :: tile_nx, tile_ny, tile_nz
    type (particle_tile), allocatable, dimension(:,:,:) :: prtl_tile
  end type particle_species

  ! particle types for exchange between processors />
  type :: prtl_enroute
    integer(kind=2)   :: xi, yi, zi
    real              :: dx, dy, dz
    real              :: u, v, w
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

  type(particle_species), target, allocatable :: species(:)
  integer                                      :: nspec

  type(prtl_enroute), allocatable, dimension(:)    :: recv_enroute
  type(enroute_handler)                            :: enroute_bot

  #ifdef MPI08
    type(MPI_DATATYPE)                             :: myMPI_ENROUTE
  #endif

  #ifdef MPI
    integer                                        :: myMPI_ENROUTE
  #endif
  
end module m_particles
