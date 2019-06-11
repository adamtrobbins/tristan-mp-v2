#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
	use m_particlelogistics
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real :: plasma_temp
  private :: plasma_temp
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitParticles, userInitFields, userReadInput
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    call userReadInput()
      call printDiag((mpi_rank .eq. 0), "...userReadInput()", .true.)

    call userInitParticles()
      call printDiag((mpi_rank .eq. 0), "...userInitParticles()", .true.)

    call userInitFields()
      call printDiag((mpi_rank .eq. 0), "...userInitFields()", .true.)
  end subroutine userInitialize

  subroutine userReadInput()
    implicit none
    call getInput('problem', 'temperature', plasma_temp)
  end subroutine userReadInput

  subroutine userInitParticles()
    implicit none
    type(region)        :: user_region
    integer             :: npart

    user_region%x_min = 0
    user_region%x_max = this_meshblock%ptr%sx - 1e-6
    user_region%y_min = 0
    user_region%y_max = this_meshblock%ptr%sy - 1e-6

    npart = this_meshblock%ptr%sx * this_meshblock%ptr%sy * ppc0

    call fillRegionWithThermalPlasma(user_region, (/1, 2/), 2, npart, plasma_temp)
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0
    ! do i = 0, this_meshblock%ptr%sx - 1
    !   i_glob = i + this_meshblock%ptr%x0
    !   do j = 0, this_meshblock%ptr%sy - 1
    !     j_glob = j + this_meshblock%ptr%y0
    !     do k = 0, this_meshblock%ptr%sz - 1
    !       k_glob = k + this_meshblock%ptr%z0
    !       ...
    !     end do
    !   end do
    ! end do
  end subroutine userInitFields

  subroutine userDriveParticles()
    implicit none
    ! integer :: s, ti, tj, tk, p
    ! do s = 1, nspec
		! 	do ti = 1, species(s)%tile_nx
		! 		do tj = 1, species(s)%tile_ny
		! 			do tk = 1, species(s)%tile_nz
    !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
    !           ...
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
  end subroutine userDriveParticles
end module m_userfile
