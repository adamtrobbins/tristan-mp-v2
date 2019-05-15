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
  real :: nCS_over_nUP, current_width, upstream_T, cs_x1, cs_x2
  private :: nCS_over_nUP, current_width, upstream_T, cs_x1, cs_x2
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitParticles, userInitFields, userReadInput
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    call userReadInput()
    call userInitParticles()
    call userInitFields()
  end subroutine userInitialize

  subroutine userReadInput()
    implicit none
    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'nCS_nUP', nCS_over_nUP)
    call getInput('problem', 'current_width', current_width)
    cs_x1 = 1.0 / 3.0; cs_x2 = 2.0 / 3.0
  end subroutine userReadInput

  subroutine userInitParticles()
    implicit none
    real                :: nUP
    integer             :: nUP_tot
    type(region)        :: back_region, cs_region

    nUP = ppc0
    nUP_tot = INT(nUP * this_meshblock%ptr%sx * this_meshblock%ptr%sy)

    back_region%x_min = 0
    back_region%x_max = this_meshblock%ptr%sx
    back_region%y_min = 0
    back_region%y_max = this_meshblock%ptr%sy

    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP_tot, upstream_T)
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0
    k = 0
    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        by(i, j, k) = tanh((i_glob + 0.5 - cs_x1 * global_mesh%sx) / current_width) -&
                    & tanh((i_glob + 0.5 - cs_x2 * global_mesh%sx) / current_width) - 1.0
      end do
    end do
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
