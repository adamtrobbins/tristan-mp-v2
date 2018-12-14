#include "../defs.F90"

module m_particles
  use m_globalnamespace
  implicit none

  type :: species
    real(mprec), allocatable, dimension(:)  :: x, y, z, u, v, w
    integer*4, allocatable, dimension(:)    :: ind, proc
    !dir$ attributes align: 64 :: x, y, z, u, v, w, ind, proc
  end type species

  type :: species_param
    integer*8   :: npart_sp
    integer*8   :: maxptl_sp
    real        :: m_sp, ch_sp
  end type species_param

  type(species), allocatable        :: sp_(:)
  type(species_param), allocatable  :: spp_(:)
  integer                           :: nspec
end module m_particles
