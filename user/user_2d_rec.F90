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
  real    :: nCS_over_nUP, current_width, upstream_T, cs_x
  real    :: injector_x1, injector_x2, injector_sx, injector_betax

  private :: nCS_over_nUP, current_width, upstream_T, cs_x
  private :: injector_x1, injector_x2, injector_sx, injector_betax
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userInitParticles, userInitFields, userReadInput,&
           & userSpatialDistribution
  !...............................................................!
contains
  subroutine userInitialize()
    implicit none
    call userReadInput()
    call userInitParticles()
    call userInitFields()
  end subroutine userInitialize

  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'nCS_nUP', nCS_over_nUP)
    call getInput('problem', 'current_width', current_width)
    call getInput('problem', 'injector_sx', injector_sx)
    call getInput('problem', 'injector_betax', injector_betax)
    cs_x = 0.5
    injector_x1 = ...
    injector_x2 = ...
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    if (present(x_glob) .and. present(dummy1) .and. present(dummy2)) then
      userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2
    else
      call throwError("ERROR: variable not present in `userSpatialDistribution()`")
    end if
    return
    return
  end function

  subroutine userInitParticles()
    implicit none
    real                :: nUP
    integer             :: nUP_tot, nCS_tot
    type(region)        :: back_region
    real                :: sx_glob, shift_gamma, shift_beta, current_sheet_T
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    nUP = ppc0
    nUP_tot = INT(0.5 * nUP * this_meshblock%ptr%sx * this_meshblock%ptr%sy)
    nCS_tot = INT(0.5 * nUP * nCS_over_nUP * this_meshblock%ptr%sx * this_meshblock%ptr%sy)

    back_region%x_min = 0
    back_region%x_max = this_meshblock%ptr%sx
    back_region%y_min = 0
    back_region%y_max = this_meshblock%ptr%sy

    sx_glob = REAL(global_mesh%sx)

    shift_beta = sqrt(sigma) * c_omp / (current_width * nCS_over_nUP)
    if (shift_beta .ge. 1) then
      call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
    end if
    shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
    current_sheet_T = 0.5 * sigma / nCS_over_nUP

    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP_tot, upstream_T)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nCS_tot, current_sheet_T,&
                                   & shift_gamma = shift_gamma, shift_dir = -3,&
                                   & spat_distr_ptr = spat_distr_ptr,&
                                   & dummy1 = cs_x * sx_glob, dummy2 = current_width)
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob
    real    :: x_glob, sx_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0
    k = 0
    sx_glob = REAL(global_mesh%sx)
    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      x_glob = REAL(i_glob)
      do j = 0, this_meshblock%ptr%sy - 1
        by(i, j, k) = tanh((x_glob - cs_x * sx_glob) / current_width)
      end do
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userDriveParticles()
    implicit none
    ! ... dummy loop ...
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
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions()
    implicit none
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions()
    implicit none
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
