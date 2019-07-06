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
  end function userSpatialDistribution

  subroutine userInitParticles()
    implicit none
    real                :: nUP
    integer             :: nUP_tot, nCS_tot
    type(region)        :: back_region
    real                :: sx_glob, sy_glob, shift_gamma, shift_beta, current_sheet_T
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    nUP = ppc0
    nUP_tot = INT(0.5 * nUP * this_meshblock%ptr%sx * this_meshblock%ptr%sy)
    nCS_tot = INT(0.5 * nUP * nCS_over_nUP * this_meshblock%ptr%sx * this_meshblock%ptr%sy)

    back_region%x_min = REAL(0)
    back_region%x_max = REAL(this_meshblock%ptr%sx)
    back_region%y_min = REAL(0)
    back_region%y_max = REAL(this_meshblock%ptr%sy)

    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)

    shift_beta = sqrt(sigma) * c_omp / (current_width * nCS_over_nUP)
    if (shift_beta .ge. 1) then
      call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
    end if
    shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
    current_sheet_T = 0.5 * sigma / nCS_over_nUP

    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP_tot, upstream_T)
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nCS_tot, current_sheet_T,&
                                   & shift_gamma = shift_gamma, shift_dir = 3,&
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

    cs_x = 0.5
    injector_x1 = injector_sx - 1.0e-3
    injector_x2 = REAL(global_mesh%sx - 1) - injector_sx + 1.0e-3

    k = 0
    sx_glob = REAL(global_mesh%sx)
    do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
      i_glob = i + this_meshblock%ptr%x0
      x_glob = REAL(i_glob)
      do j = -NGHOST, this_meshblock%ptr%sy - 1 + NGHOST
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
    !   do ti = 1, species(s)%tile_nx
    !     do tj = 1, species(s)%tile_ny
    !       do tk = 1, species(s)%tile_nz
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
  function userSpatialDistributionInjector(x_glob, y_glob, z_glob,&
                                         & dummy1, dummy2, dummy3)
    real :: userSpatialDistributionInjector
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    integer                     :: i_glob, injector_i1_glob, injector_i2_glob
    ! only inject in cell columns along the injectors
    if (present(x_glob) .and. present(dummy1) .and. present(dummy2)) then
      i_glob = INT(x_glob)
      injector_i1_glob = INT(dummy1)
      injector_i2_glob = INT(dummy2)
      if ((i_glob .eq. injector_i1_glob) .or. (i_glob .eq. injector_i2_glob)) then
        userSpatialDistributionInjector = 1.0
      else
        userSpatialDistributionInjector = 0.0
      end if
    else
      call throwError("ERROR: variable not present in `userSpatialDistribution()`")
    end if
    return
    return
  end function userSpatialDistributionInjector

  subroutine userParticleBoundaryConditions()
    implicit none
    real          :: nUP
    integer       :: s, ti, tj, tk, p, nUP_tot
    integer       :: i_glob, injector_i1_glob, injector_i2_glob
    type(region)  :: back_region
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistributionInjector

    ! move the injectors
    injector_x1 = injector_x1 - injector_betax * CC
    injector_x2 = injector_x2 + injector_betax * CC
    ! reset the injector positions if necessary
    if (injector_x1 .le. 0.0) then
      injector_x1 = injector_x1 + injector_sx
    end if
    if (injector_x2 .ge. REAL(global_mesh%sx - 1)) then
      injector_x2 = injector_x2 - injector_sx
    end if

    injector_i1_glob = INT(injector_x1)
    injector_i2_glob = INT(injector_x2)

    if (((injector_i1_glob .ge. this_meshblock%ptr%x0) .and. (injector_i1_glob .lt. this_meshblock%ptr%x0 + this_meshblock%ptr%sx)) .or.&
      & ((injector_i2_glob .ge. this_meshblock%ptr%x0) .and. (injector_i2_glob .lt. this_meshblock%ptr%x0 + this_meshblock%ptr%sx))) then

      ! remove particles left and right from the injectors
      do s = 1, nspec
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                i_glob = species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0
                if ((i_glob .le. injector_i1_glob) .or. (i_glob .ge. injector_i2_glob)) then
                  call removeParticleFromTile(s, ti, tj, tk, p)
                end if
              end do
            end do
          end do
        end do
      end do

      ! inject background particles at the injectors' positions
      nUP = ppc0
      nUP_tot = INT(0.5 * nUP * this_meshblock%ptr%sy)

      ! left injector
      back_region%x_min = REAL(injector_i1_glob - this_meshblock%ptr%x0)
      back_region%x_max = REAL(injector_i1_glob - this_meshblock%ptr%x0 + 1)
      back_region%y_min = REAL(0)
      back_region%y_max = REAL(this_meshblock%ptr%sy)

      call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP_tot, upstream_T)
      ! ,&
      !                                & spat_distr_ptr = spat_distr_ptr,&
      !                                & dummy1 = injector_x1, dummy2 = injector_x2)

      ! right injector
      back_region%x_min = REAL(injector_i2_glob - this_meshblock%ptr%x0)
      back_region%x_max = REAL(injector_i2_glob - this_meshblock%ptr%x0 + 1)
      back_region%y_min = REAL(0)
      back_region%y_max = REAL(this_meshblock%ptr%sy)

      call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, nUP_tot, upstream_T)
      ! ,&
      !                               & spat_distr_ptr = spat_distr_ptr,&
      !                               & dummy1 = injector_x1, dummy2 = injector_x2)
    end if
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions()
    implicit none
    real    :: sx_glob, x_glob
    integer :: i, j, k
    integer :: i_glob, injector_i1_glob, injector_i2_glob
    injector_i1_glob = INT(injector_x1)
    injector_i2_glob = INT(injector_x2)

    if ((injector_i1_glob .le. this_meshblock%ptr%x0 + this_meshblock%ptr%sx) .or.&
      & (injector_i2_glob .ge. this_meshblock%ptr%x0)) then

      ! reset fields left and right from the injectors
      ! FIX0: do I need to do anything with the currents?
      k = 0
      sx_glob = REAL(global_mesh%sx)
      do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
        i_glob = i + this_meshblock%ptr%x0
        x_glob = REAL(i_glob)
        do j = -NGHOST, this_meshblock%ptr%sy - 1 + NGHOST
          if ((i_glob .le. injector_i1_glob) .or. (i_glob .ge. injector_i2_glob)) then
            ex(i, j, k) = 0; ey(i, j, k) = 0; ez(i, j, k) = 0
            bx(i, j, k) = 0; bz(i, j, k) = 0
            by(i, j, k) = tanh((x_glob - cs_x * sx_glob) / current_width)
          end if
        end do
      end do

    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
