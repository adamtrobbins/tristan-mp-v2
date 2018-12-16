#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_communications
  use m_domain
  use m_particles
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitParticles
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    integer        :: npart
    integer        :: ppc0
    call getInput('particles', 'ppc0', ppc0)
    npart = this_meshblock%ptr%sx * this_meshblock%ptr%sy * this_meshblock%ptr%sz * ppc0
    call userInitParticles(npart)
  end subroutine userInitialize

  subroutine userInitParticles(npart)
    ! FIX ghost zones to be included here
    implicit none
    integer, intent(in) :: npart
    integer             :: p
    real(mprec)         :: x1, y1, z1
    z1 = 0.5
    do p = 1, npart
      x1 = random(dseed) * this_meshblock%ptr%sx
      y1 = random(dseed) * this_meshblock%ptr%sy
      #ifdef threeD
        z1 = random(dseed) * this_meshblock%ptr%sz
      #endif
      spp_(1)%npart_sp = spp_(1)%npart_sp + 1
      spp_(2)%npart_sp = spp_(2)%npart_sp + 1
      sp_(1)%x(p) = x1; sp_(2)%x(p) = x1
      sp_(1)%y(p) = y1; sp_(2)%y(p) = y1
      sp_(1)%z(p) = z1; sp_(2)%z(p) = z1
      sp_(1)%u(p) = 0; sp_(2)%u(p) = 0
      sp_(1)%v(p) = 0; sp_(2)%v(p) = 0
      sp_(1)%w(p) = 0; sp_(2)%w(p) = 0
      sp_(1)%ind(p) = spp_(1)%npart_sp; sp_(2)%ind(p) = spp_(2)%npart_sp
      sp_(1)%proc(p) = my_rank; sp_(2)%proc(p) = my_rank
    end do
  end subroutine userInitParticles
end module m_userfile
