#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_communications
  use m_domain
  use m_particles
  use m_fields
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitParticles, userInitFields
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    integer        :: npart
    call userInitParticles()
    call userInitFields()
  end subroutine userInitialize

  subroutine userInitParticles()
    implicit none
    integer(kind=2)     :: xi_, yi_, zi_
    if (mpi_rank .eq. 0) then
      xi_ = 12; yi_ = 5; zi_ = 0
      call createParticle(1, xi_, yi_, zi_, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0)
      xi_ = 10; yi_ = 2; zi_ = 0
      call createParticle(2, xi_, yi_, zi_, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0)
    end if
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = -1; ey(:,:,:) = -1; ez(:,:,:) = -1
    bx(:,:,:) = -1; by(:,:,:) = -1; bz(:,:,:) = -1
    do i = 0, this_meshblock%ptr%sx - 1
      do j = 0, this_meshblock%ptr%sy - 1
        do k = 0, this_meshblock%ptr%sz - 1
          ex(i, j, k) = 0
          ey(i, j, k) = 0.1
          ez(i, j, k) = 0
          bx(i, j, k) = 0
          by(i, j, k) = 0
          bz(i, j, k) = 1
        end do
      end do
    end do
  end subroutine userInitFields
end module m_userfile
