#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
	use m_particlelogistics
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitParticles, userInitFields
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    integer        :: npart
    npart = INT(this_meshblock%ptr%sx * this_meshblock%ptr%sy * this_meshblock%ptr%sz * ppc0)
    call userInitParticles(npart)
    call userInitFields()
  end subroutine userInitialize

  subroutine userInitParticles(npart)
    implicit none
    integer, intent(in) :: npart
    integer(kind=2)     :: xi_, yi_, zi_
    integer(kind=2)     :: i, j
    integer             :: p, s
    real                :: dx_, dy_, dz_, u_, v_

    #ifndef threeD
      dz_ = 0.5; zi_ = 0
    #else
      dz_ = 10; zi_ = 0.43
    #endif

    ! if (mpi_rank .eq. 1) then
    u_ = 2 * (random(dseed) - 0.5)
    v_ = 2 * (random(dseed) - 0.5)
    do i = 0, this_meshblock%ptr%sx - 1
      do j = 0, this_meshblock%ptr%sy - 1
        if (((i - this_meshblock%ptr%sx / 2.0) / (this_meshblock%ptr%sx / 3.0))**2 +&
          & ((j - this_meshblock%ptr%sy / 2.0) / (this_meshblock%ptr%sy / 3.0))**2 < 1) then
          call createParticle(1, i, j, zi_, 0.5, 0.5, dz_, u_, v_, 0.0)
        end if
      end do
    end do
    ! endif

  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: ind1, ind2, ind3
    integer :: i_glob, j_glob, k_glob
    do ind1 = 0, this_meshblock%ptr%sx - 1
      do ind2 = 0, this_meshblock%ptr%sy - 1
        do ind3 = 0, this_meshblock%ptr%sz - 1
          ex(ind1, ind2, ind3) = mpi_rank
        end do
      end do
    end do
  end subroutine userInitFields
end module m_userfile
