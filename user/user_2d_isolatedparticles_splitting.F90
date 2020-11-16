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

  procedure (spatialDistribution), pointer :: user_slb_load_ptr => userSLBload

  !--- PRIVATE variables -----------------------------------------!

  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    return
  end function

  function userSLBload(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    return
  end function

  function userSplitting(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    logical :: userSplitting
    ! global coordinates
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    real :: splitregion_xmin = 24.5
    real :: splitregion_xmax = 25.5
    userSplitting = ((x_glob.ge.splitregion_xmin).and.(x_glob.le.splitregion_xmax))
    return 
  end function

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field

    integer               :: n_split, n_split_frac, split_ind
    integer, allocatable  :: ind_split(:)
    integer               :: s, ti, tj, tk, p, i
    real                  :: x_glob, y_glob, z_glob
    real                  :: weight_D
    real                  :: x0, y0, z0, vx, vy, vz, dxyz
   
    integer               :: p_ind
    real                  :: weight_min, wiggler, splitfactor
    real                  :: alpha, beta, gamma, xA, yA, zA, xB, yB, zB

    ! MINIMUM WEIGHT. SPLITTING IS DONE BY FACTORS OF 2.
    weight_min = 0.25
    ! FRACTION OF PARTICLES IN A SPLITTING BIN THAT WILL BE RANDOMLY SELECTED.
    splitfactor = 1.0
    ! PERTURBATION TO THE POSITION OF DAUGHTER PARTICLES IN UNITS OF THE CELL WIDTH.
    wiggler = 0.25

    ! DECIDE WHICH SPECIES IS AFFECTED BY SPLITTING
    do s = 1, 1

      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            
            allocate(ind_split(species(s)%prtl_tile(ti, tj, tk)%npart_sp))
            n_split = 0

            ! IDENTIFY SPLITTABLE PARTICLES THAT ARE LOCATED IN THE SPLITTING REGION 
            ! AND THROW THEIR INDICES INTO A BIN 
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              
              x_glob = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                     & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
              y_glob = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
                     & + species(s)%prtl_tile(ti, tj, tk)%dy(p)
              z_glob = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
                     & + species(s)%prtl_tile(ti, tj, tk)%dz(p)

              ! CHECK EACH PARTICLE FOR SPLITTING CRITERION AND WEIGHT (userSplitting)
              if(userSplitting(x_glob, y_glob, z_glob).and.((0.5*species(s)%prtl_tile(ti, tj, tk)%weight(p)).gt.(weight_min*0.99))) then
                  n_split = n_split + 1
                  ind_split(n_split) = p
              endif

            end do

            ! DETERMINE THE NUMBER OF PARTICLES THAT SHOULD BE SPLIT BY A FRACTION
            ! OF THE TOTAL NUMBER OF PARTICLES IN THE BIN
            n_split_frac = floor(n_split * splitfactor)
            do i = 1, n_split_frac

              ! RANDOMLY SELECT PARTICLES FROM THE BIN THAT HAVE NOT YET BEEN SELECTED
              split_ind = -1
              do while (split_ind.lt.0)
                 p_ind = INT((random(dseed) * n_split + 1)) 
                 split_ind = ind_split(p_ind)
              end do

              ! LOAD POSITION AND VELOCITY OF THE PARENT PARTICLE
              x0 = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(split_ind)) + species(s)%prtl_tile(ti, tj, tk)%dx(split_ind)
              y0 = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(split_ind)) + species(s)%prtl_tile(ti, tj, tk)%dy(split_ind)
              z0 = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(split_ind)) + species(s)%prtl_tile(ti, tj, tk)%dz(split_ind)
              vx = species(s)%prtl_tile(ti, tj, tk)%u(split_ind)
              vy = species(s)%prtl_tile(ti, tj, tk)%v(split_ind)
              vz = species(s)%prtl_tile(ti, tj, tk)%w(split_ind)

              ! DETERMINE WEIGHT OF DAUGHTER PARTICLES (HALF THE WEIGHT OF THE PARENT)
              weight_D = species(s)%prtl_tile(ti, tj, tk)%weight(split_ind) / 2.0
              
              ! ADD A POSITION PERTURBATION (THAT DEPENDS ON THE REFINEMENT LEVEL)
              dxyz = wiggler * weight_D * 2.0
              xA = dxyz
              xB = - dxyz
              yA = 0.0
              yB = 0.0
              zA = 0.0
              zB = 0.0

              ! ROTATE THE TWO PARTICLES IN A PLANE PERPENDICULAR TO THE 
              ! MOTION OF THE PARENT PARTICLE. IN 3D ADD A RANDOM ORBITAL ROTATION.
              #ifdef twoD

                alpha = acos(vx / sqrt(vx**2 + vy**2))
                call rotateIn3D(xA, yA, zA, 0.0, 0.0, 1.0, alpha - M_PI/2.0) 
                call rotateIn3D(xB, yB, zB, 0.0, 0.0, 1.0, alpha - M_PI/2.0) 

              #elif threeD

                alpha = acos(vx / sqrt(vx**2 + vy**2))
                beta  = acos(vz / sqrt(vx**2 + vy**2 + vz**2))
                gamma = 2.0 * random(dseed) * M_PI
                call rotateIn3D(xA, yA, zA, 0.0, 1.0, 0.0, beta) 
                call rotateIn3D(xB, yB, zB, 0.0, 1.0, 0.0, beta)     
                call rotateIn3D(xA, yA, zA, 0.0, 0.0, 1.0, alpha) 
                call rotateIn3D(xB, yB, zB, 0.0, 0.0, 1.0, alpha) 
                call rotateIn3D(xA, yA, zA, vx, vy, vz, gamma) 
                call rotateIn3D(xB, yB, zB, vx, vy, vz, gamma) 

              #endif

                xA = xA + x0
                xB = xB + x0
                yA = yA + y0
                yB = yB + y0
                zA = zA + z0
                zB = zB + z0

              ! CREATE NEW PARTICLES AT NEW POSITIONS WITH NEW WEIGHT
              call injectParticleLocally(s, xA, yA, zA, vx, vy, vz, weight=weight_D)
              call injectParticleLocally(s, xB, yB, zB, vx, vy, vz, weight=weight_D)

              ! DEPOSIT A CURRENT CORRECTION TO ACCOUNT FOR INSTANTANEOUS POSITION CHANGE
              call depositCurrentsFromSingleParticle(s, species(s)%prtl_tile(ti, tj, tk), species(s)%prtl_tile(ti, tj, tk)%npart_sp - 1, x0, y0, z0, xA, yA, zA)
              call depositCurrentsFromSingleParticle(s, species(s)%prtl_tile(ti, tj, tk), species(s)%prtl_tile(ti, tj, tk)%npart_sp, x0, y0, z0, xB, yB, zB)

              ! FLAGG PARENT PARTICLE IN THE TILE ARRAY AND THE SPLITTING BIN
              species(s)%prtl_tile(ti, tj, tk)%proc(split_ind) = -1
              ind_split(p_ind) = -1

            end do

            deallocate(ind_split)

          end do
        end do
      end do
    end do 

  end subroutine userCurrentDeposit

  subroutine userInitParticles()
    implicit none
    real            :: u, v, w, vx, vy, vz, gamma
    real            :: xg, yg, zg
    real            :: vpm
    integer         :: ntot, n, i, j, k
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution


    do n = 1, ppc0
      xg = 30.0 + random(dseed) * 15.0
      yg = 7.5 + random(dseed) * 5.0
      zg = 0.5

      vx = (25.5 - xg) * 0.01
      vy = (25.5 - yg) * 0.01
      vz = 0.0
      gamma = 1.0 / sqrt(1.0 - vx**2.0 - vy**2.0 - vz**2.0)

      u = gamma*vx
      v = gamma*vy
      w = gamma*vz
      ! particles at rest
      call injectParticleGlobally(1, xg, yg, zg, u, v, w)
      call injectParticleGlobally(2, xg, yg, zg, u, v, w)
    end do


  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0
    ! ... dummy loop ...
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
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
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

  subroutine userExternalFields(xp, yp, zp,&
                              & ex_ext, ey_ext, ez_ext,&
                              & bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in)  :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step

  
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical                       :: updateE_, updateB_

    if (present(updateE)) then
      updateE_ = updateE
    else
      updateE_ = .true.
    end if

    if (present(updateB)) then
      updateB_ = updateB
    else
      updateB_ = .true.
    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
