#include "../defs.F90"

module m_particles
  use m_globalnamespace
  implicit none

  type :: species
    integer(kind=2), allocatable, dimension(:)  :: xi, yi, zi
    real, allocatable, dimension(:)             :: dx, dy, dz
    real, allocatable, dimension(:)             :: u, v, w
    integer, allocatable, dimension(:)          :: ind, proc
    !dir$ attributes align: 64 :: xi, yi, zi, dx, dy, dz, u, v, w, ind, proc
  end type species

  type :: species_param
    integer     :: npart_sp, maxptl_sp, cntr_sp
    real        :: m_sp, ch_sp
  end type species_param

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

  type(species), target, allocatable :: sp_(:)
  type(species_param), allocatable   :: spp_(:)
  integer                            :: nspec

  type(prtl_enroute), allocatable, dimension(:)    :: recv_enroute
  type(enroute_handler)                            :: enroute_bot
  type(MPI_DATATYPE)                               :: myMPI_ENROUTE
contains
  subroutine copyParticleFromTo(s, p_from, p_to)
    implicit none
    integer, intent(in)   :: s, p_from, p_to
    sp_(s)%xi(p_to) = sp_(s)%xi(p_from); sp_(s)%yi(p_to) = sp_(s)%yi(p_from); sp_(s)%zi(p_to) = sp_(s)%zi(p_from)
    sp_(s)%dx(p_to) = sp_(s)%dx(p_from); sp_(s)%dy(p_to) = sp_(s)%dy(p_from); sp_(s)%dz(p_to) = sp_(s)%dz(p_from)
    sp_(s)%u(p_to) = sp_(s)%u(p_from); sp_(s)%v(p_to) = sp_(s)%v(p_from); sp_(s)%w(p_to) = sp_(s)%w(p_from)
    sp_(s)%ind(p_to) = sp_(s)%ind(p_from); sp_(s)%proc(p_to) = sp_(s)%proc(p_from)
  end subroutine copyParticleFromTo

  subroutine removeParticle(s, p)
    implicit none
    integer, intent(in)   :: s, p
    if (p .ne. spp_(s)%npart_sp) call copyParticleFromTo(s, spp_(s)%npart_sp, p)
    spp_(s)%npart_sp = spp_(s)%npart_sp - 1
  end subroutine removeParticle

  subroutine createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w)
    implicit none
    integer, intent(in)           :: s
    integer(kind=2), intent(in)   :: xi, yi, zi
    real, intent(in)              :: dx, dy, dz, u, v, w
    integer                       :: p
    spp_(s)%npart_sp = spp_(s)%npart_sp + 1
    p = spp_(s)%npart_sp
    sp_(s)%xi(p) = xi; sp_(s)%dx(p) = dx
    sp_(s)%yi(p) = yi; sp_(s)%dy(p) = dy
    sp_(s)%zi(p) = zi; sp_(s)%dz(p) = dz
    sp_(s)%u(p) = u; sp_(s)%v(p) = v; sp_(s)%w(p) = w
    sp_(s)%ind(p) = spp_(s)%cntr_sp; sp_(s)%proc(p) = mpi_rank
    spp_(s)%cntr_sp = spp_(s)%cntr_sp + 1
  end subroutine createParticle

end module m_particles
