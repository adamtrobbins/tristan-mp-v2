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
    real           :: ppc0
    call getInput('particles', 'ppc0', ppc0)
    npart = INT(this_meshblock%ptr%sx * this_meshblock%ptr%sy * this_meshblock%ptr%sz * ppc0)
    call userInitParticles(npart)
  end subroutine userInitialize

  subroutine userInitParticles(npart)
    implicit none
    integer, intent(in) :: npart
    integer(kind=2)     :: xi_, yi_, zi_
    integer             :: p
    real                :: dx_, dy_, dz_

    #ifndef threeD
      dz_ = 0.5; zi_ = 0
    #else
      dz_ = 10; zi_ = 0.43
    #endif

    if (mpi_rank .eq. 0) then
      p = 1
      sp_(1)%xi(p) = 10; sp_(1)%dx(p) = 0.2
      sp_(1)%yi(p) = 10; sp_(1)%dy(p) = 0.2
      sp_(1)%zi(p) = zi_; sp_(1)%dz(p) = dz_
      sp_(1)%u = 1.5; sp_(1)%v = 1.5; sp_(1)%w = 0

      sp_(1)%ind(p) = spp_(1)%cntr_sp; sp_(1)%proc(p) = mpi_rank

      spp_(1)%npart_sp = spp_(1)%npart_sp + 1; spp_(1)%cntr_sp = spp_(1)%cntr_sp + 1
    endif

    ! if (mpi_rank .eq. 3) then
    !   p = 1
    !   sp_(1)%xi(p) = 40; sp_(1)%dx(p) = 0.356
    !   sp_(1)%yi(p) = 60; sp_(1)%dy(p) = 0.23
    !   sp_(1)%zi(p) = zi_; sp_(1)%dz(p) = dz_
    !   sp_(1)%u = 0; sp_(1)%v = 2.6; sp_(1)%w = 0
    !
    !   sp_(1)%ind(p) = spp_(1)%cntr_sp; sp_(1)%proc(p) = mpi_rank
    !
    !   spp_(1)%npart_sp = spp_(1)%npart_sp + 1; spp_(1)%cntr_sp = spp_(1)%cntr_sp + 1
    ! endif

    ! do p = 1, npart
    !   dx_ = random(dseed); dy_ = random(dseed);
    !   xi_ = INT(random(dseed) * this_meshblock%ptr%sx, kind(2))
    !   yi_ = INT(random(dseed) * this_meshblock%ptr%sy, kind(2))
    !   #ifdef threeD
    !     dz_ = random(dseed); zi_ = INT(random(dseed) * this_meshblock%ptr%sz)
    !   #endif
    !   sp_(1)%xi(p) = xi_; sp_(1)%dx(p) = dx_
    !   sp_(1)%yi(p) = yi_; sp_(1)%dy(p) = dy_
    !   sp_(1)%zi(p) = zi_; sp_(1)%dz(p) = dz_
    !   sp_(2)%xi(p) = xi_; sp_(2)%dx(p) = dx_
    !   sp_(2)%yi(p) = yi_; sp_(2)%dy(p) = dy_
    !   sp_(2)%zi(p) = zi_; sp_(2)%dz(p) = dz_
    !
    !   sp_(1)%u(p) = random(dseed); sp_(2)%u(p) = random(dseed)
    !   sp_(1)%v(p) = random(dseed); sp_(2)%v(p) = random(dseed)
    !   sp_(1)%w(p) = random(dseed); sp_(2)%w(p) = random(dseed)
    !
    !   sp_(1)%ind(p) = spp_(1)%cntr_sp; sp_(2)%ind(p) = spp_(2)%cntr_sp
    !   sp_(1)%proc(p) = mpi_rank; sp_(2)%proc(p) = mpi_rank
    !
    !   spp_(1)%npart_sp = spp_(1)%npart_sp + 1; spp_(2)%npart_sp = spp_(2)%npart_sp + 1
    !   spp_(1)%cntr_sp = spp_(1)%cntr_sp + 1; spp_(2)%cntr_sp = spp_(2)%cntr_sp + 1
    ! end do
  end subroutine userInitParticles
end module m_userfile
