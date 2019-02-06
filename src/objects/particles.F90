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
  integer                                          :: myMPI_ENROUTE
contains
  subroutine removeParticle(s, p)
    implicit none
    integer, intent(in)   :: s, p
    integer               :: nprt
    nprt = spp_(s)%npart_sp
    sp_(s)%xi(p) = sp_(s)%xi(nprt); sp_(s)%yi(p) = sp_(s)%yi(nprt); sp_(s)%zi(p) = sp_(s)%zi(nprt)
    sp_(s)%dx(p) = sp_(s)%dx(nprt); sp_(s)%dy(p) = sp_(s)%dy(nprt); sp_(s)%dz(p) = sp_(s)%dz(nprt)
    sp_(s)%u(p) = sp_(s)%u(nprt); sp_(s)%v(p) = sp_(s)%v(nprt); sp_(s)%w(p) = sp_(s)%w(nprt)
    sp_(s)%ind(p) = sp_(s)%ind(nprt); sp_(s)%proc(p) = sp_(s)%proc(nprt)
    spp_(s)%npart_sp = spp_(s)%npart_sp - 1
  end subroutine

  subroutine copyToEnroute(spec_id, prtl_id, enroute)
    implicit none
    integer, intent(in)               :: spec_id, prtl_id
    type(prtl_enroute), intent(inout) :: enroute
    enroute%xi = sp_(spec_id)%xi(prtl_id)
    enroute%yi = sp_(spec_id)%yi(prtl_id)
    enroute%zi = sp_(spec_id)%zi(prtl_id)
    enroute%dx = sp_(spec_id)%dx(prtl_id)
    enroute%dy = sp_(spec_id)%dy(prtl_id)
    enroute%dz = sp_(spec_id)%dz(prtl_id)
    enroute%u = sp_(spec_id)%u(prtl_id)
    enroute%v = sp_(spec_id)%v(prtl_id)
    enroute%w = sp_(spec_id)%w(prtl_id)
    enroute%ind = sp_(spec_id)%ind(prtl_id)
    enroute%proc = sp_(spec_id)%proc(prtl_id)
  end subroutine

  subroutine copyFromEnroute(enroute, spec_id, prtl_id)
    implicit none
    type(prtl_enroute), intent(inout) :: enroute
    integer, intent(in)               :: spec_id, prtl_id
    sp_(spec_id)%xi(prtl_id) = enroute%xi
    sp_(spec_id)%yi(prtl_id) = enroute%yi
    sp_(spec_id)%zi(prtl_id) = enroute%zi
    sp_(spec_id)%dx(prtl_id) = enroute%dx
    sp_(spec_id)%dy(prtl_id) = enroute%dy
    sp_(spec_id)%dz(prtl_id) = enroute%dz
    sp_(spec_id)%u(prtl_id) = enroute%u
    sp_(spec_id)%v(prtl_id) = enroute%v
    sp_(spec_id)%w(prtl_id) = enroute%w
    sp_(spec_id)%ind(prtl_id) = enroute%ind
    sp_(spec_id)%proc(prtl_id) = enroute%proc
  end subroutine

  subroutine extractParticlesFromEnroute(cnt, spec_id)
    implicit none
    integer, intent(in)            :: cnt, spec_id
    integer                        :: p
    do p = 1, cnt
      spp_(spec_id)%npart_sp = spp_(spec_id)%npart_sp + 1
      call copyFromEnroute(recv_enroute(p), spec_id, spp_(spec_id)%npart_sp)
    end do
  end subroutine

end module m_particles
