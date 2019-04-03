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
    real           :: ppc0
    call getInput('particles', 'ppc0', ppc0)
    ! npart = INT(this_meshblock%ptr%sx * this_meshblock%ptr%sy * this_meshblock%ptr%sz * ppc0)
    ! call userInitParticles(npart)
    call userInitFields()
  end subroutine userInitialize

  subroutine userInitParticles()
    implicit none
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = -1
    do i = 0, this_meshblock%ptr%sx - 1
      do j = 0, this_meshblock%ptr%sy - 1
        do k = 0, this_meshblock%ptr%sz - 1
          ex(i, j, k) = i + this_meshblock%ptr%x0
          ey(i, j, k) = j + this_meshblock%ptr%y0
          ez(i, j, k) = -(this_meshblock%ptr%y0 + j)**2
          bx(i, j, k) = mpi_rank * 100 + 10 * i + j
          by(i, j, k) = -mpi_rank
          bz(i, j, k) = mpi_rank + 123.5
        end do
      end do
    end do
  end subroutine userInitFields
end module m_userfile
