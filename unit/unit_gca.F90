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
  use m_helpers
  implicit none

  procedure (spatialDistribution), pointer :: user_slb_load_ptr => userSLBload

  !--- PRIVATE variables -----------------------------------------!
  real, private     :: velocity, angle
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    ! B-field geometry: 1 = monopole, 2 = dipole
    call getInput('problem', 'angle', angle)
    call getInput('problem', 'velocity', velocity)
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
    real                        :: radius2
    radius2 = (dummy1 * 0.5 - x_glob)**2 + (dummy2 * 0.5 - y_glob)**2 + (dummy3 * 0.5 - z_glob)**2 + 1.0
    userSLBload = 40**2 / radius2
    if (radius2 .lt. 40**2) then
      userSLBload = 1.0 / exp((40**2 - radius2) / 40**2)
    end if
    return
  end function

  subroutine userInitParticles()
    implicit none
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real    :: bx0, by0, bz0
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      bx(i,:,:) = (REAL(i_glob) - REAL(global_mesh%sx) * 0.5) / REAL(global_mesh%sx)
    end do
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
    integer                       :: s, ti, tj, tk, p
    real                          :: x_g, y_g, z_g, r_g
    real                          :: Ux, Uy, Uz, posx, posy, posz

    posx = REAL(global_mesh%sx) * 0.05
    posz = 0.5

    Ux = velocity * cos(angle * M_PI / 180.0)
    Uy = velocity * sin(angle * M_PI / 180.0)
    Uz = 0.0
    posy = REAL(global_mesh%sy) * 0.25
    if (step .eq. 0) call injectParticleGlobally(1, posx, posy, posz, Ux, Uy, Uz)
    posy = REAL(global_mesh%sy) * 0.75
    if (step .eq. 0) call injectParticleGlobally(2, posx, posy, posz, Ux, Uy, Uz)

    ! remove particles near Y-boundaries
    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                  & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
              if ((x_g .lt. 5.0) .or. (x_g .gt. REAL(global_mesh%sx) - 6.0)) then
                species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
              end if
            end do
          end do
        end do
      end do
    end do
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical                       :: updateE_, updateB_
    integer                       :: i, j, k, i_glob, j_glob, k_glob, mx, my, mz
    real                          :: supersph_radius_sq
    real, allocatable             :: ex_new(:,:,:), ey_new(:,:,:), ez_new(:,:,:),&
                                   & bx_new(:,:,:), by_new(:,:,:), bz_new(:,:,:)
    real                          :: rx, ry, rz, rr_sqr, shift_B, s
    real                          :: bx_dip, by_dip, bz_dip, b_dip_dot_r, b_int_dot_r
    real                          :: scaleEpar, scaleEperp, scaleBperp, scaleBpar, scale
    real                          :: vx, vy, vz, ex_dip, ey_dip, ez_dip
    real                          :: shift_E, e_int_dot_r, e_dip_dot_r

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

  !--- auxiliary functions ------------------------------------!
  !............................................................!
end module m_userfile
