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
    call userInitFields()
  end subroutine userInitialize

  subroutine userDriveParticles()
    implicit none
    ! do nothing
  end subroutine userDriveParticles

  subroutine userInitParticles()
    implicit none
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real :: kx, ky, ex_norm, ey_norm, exy_norm

    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0
    jx(20, 2, 0) = 1000
    jy(20, 2, 0) = 1000
    jz(20, 2, 0) = 1000
  end subroutine userInitFields
end module m_userfile
