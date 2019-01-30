#include "../defs.F90"

module m_particles
  use m_globalnamespace
  implicit none

  type :: species
    real(mprec), allocatable, dimension(:)  :: x, y, z, u, v, w
    integer, allocatable, dimension(:)      :: ind, proc
    !dir$ attributes align: 64 :: x, y, z, u, v, w, ind, proc
  end type species

  type :: species_param
    integer     :: npart_sp
    integer     :: maxptl_sp
    real        :: m_sp, ch_sp
  end type species_param

  type(species), target, allocatable :: sp_(:)
  type(species_param), allocatable   :: spp_(:)
  integer                            :: nspec

  type(sometype), allocatable        :: part_send_(:), part_recv_(:)
end module m_particles
