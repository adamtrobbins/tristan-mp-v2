#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  implicit none

  procedure (spatialDistribution), pointer :: user_slb_load_ptr => null()

  !--- PRIVATE variables -----------------------------------------!
  integer   :: ph_ndot
  real      :: ph_energy

  private   :: ph_ndot, ph_energy
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
    user_slb_load_ptr => userSLBload
  end subroutine userInitialize

  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'ndot', ph_ndot)
    call getInput('problem', 'energy', ph_energy)
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

  subroutine userInitParticles()
    implicit none
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
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
    integer                       :: i, thet, rnd
    real                          :: u_, v_, w_
    real                          :: x1_g, y1_g, x2_g, y2_g
    real                          :: x1_l, y1_l, x2_l, y2_l
    real                          :: dx_, dy_, dz_, x_, y_
    integer                       :: xi_, yi_, zi_

    x1_g = global_mesh%sx * 0.25; y1_g = global_mesh%sy * 0.5
    x2_g = global_mesh%sx * 0.75; y2_g = global_mesh%sy * 0.5

    dz_ = 0.5; zi = 0

    call globalToLocalCoords(x1_g, y1_g, 0.0,&
                           & x1_l, y1_l, rnd)
    call globalToLocalCoords(x2_g, y2_g, 0.0,&
                           & x2_l, y2_l, rnd)

    if (((x1_l .ge. 0) .and. (x1_l .lt. this_meshblock%ptr%sx)) .and.&
      & ((y1_l .ge. 0) .and. (y1_l .lt. this_meshblock%ptr%sy))) then

      do i = 1, ph_ndot
        rnd = random(dseed)
        thet = 2 * M_PI * rnd
        u_ = cos(thet) * ph_energy
        v_ = sin(thet) * ph_energy
        w_ = 0.0

        x_ = x1_l + u_ * (random(dseed) - 0.5) * 2
        y_ = y1_l + v_ * (random(dseed) - 0.5) * 2

        xi_ = INT(FLOOR(x_), 2); dx_ = x_ - FLOOR(x_)
        if (xi_ .eq. this_meshblock%ptr%sx) then
          xi_ = xi_ - 1; dx_ = dx_ + 1.0
        end if
        yi_ = INT(FLOOR(y_), 2); dy_ = y_ - FLOOR(y_)
        if (yi_ .eq. this_meshblock%ptr%sy) then
          yi_ = yi_ - 1; dy_ = dy_ + 1.0
        end if

        call createParticle(1, xi_, yi_, yi, dx_, dy_, dz_, u_, v_, w_)
      end do
    end if

    if (((x2_l .ge. 0) .and. (x2_l .lt. this_meshblock%ptr%sx)) .and.&
      & ((y2_l .ge. 0) .and. (y2_l .lt. this_meshblock%ptr%sy))) then
      do i = 1, ph_ndot
        rnd = random(dseed)
        thet = 2 * M_PI * rnd
        u_ = cos(thet) * ph_energy
        v_ = sin(thet) * ph_energy
        w_ = 0.0

        x_ = x2_l + u_ * (random(dseed) - 0.5) * 2
        y_ = y2_l + v_ * (random(dseed) - 0.5) * 2

        xi_ = INT(FLOOR(x_), 2); dx_ = x_ - FLOOR(x_)
        if (xi_ .eq. this_meshblock%ptr%sx) then
          xi_ = xi_ - 1; dx_ = dx_ + 1.0
        end if
        yi_ = INT(FLOOR(y_), 2); dy_ = y_ - FLOOR(y_)
        if (yi_ .eq. this_meshblock%ptr%sy) then
          yi_ = yi_ - 1; dy_ = dy_ + 1.0
        end if

        call createParticle(2, xi_, yi_, yi, dx_, dy_, dz_, u_, v_, w_)
      end do
    end if

  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
