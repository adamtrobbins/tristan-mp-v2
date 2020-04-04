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

  procedure (spatialDistribution), pointer :: user_slb_load_ptr => null()

  !--- PRIVATE variables -----------------------------------------!
  real :: xc_g, yc_g, zc_g
  real :: psr_angle, psr_period, psr_omega, psr_radius

  private :: xc_g, yc_g, zc_g, psr_angle, psr_period, psr_omega, psr_radius
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
    call getInput('problem', 'psr_radius', psr_radius)
    call getInput('problem', 'psr_angle', psr_angle)
    call getInput('problem', 'psr_period', psr_period)

    psr_angle = psr_angle * M_PI / 180.0
    psr_omega = 2.0 * M_PI / psr_period

    xc_g = 0.5 * global_mesh%sx
    yc_g = 0.5 * global_mesh%sy
    zc_g = 0.5 * global_mesh%sz
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
    real    :: bx0, by0, bz0
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        do k = 0, this_meshblock%ptr%sz - 1
          k_glob = k + this_meshblock%ptr%z0
          call getDipole(0, 0.0, REAL(i_g), REAL(j_g) + 0.5, REAL(k_g) + 0.5, bx0, by0, bz0)
          bx(i, j, k) = bx0
          call getDipole(0, 0.0, REAL(i_g) + 0.5, REAL(j_g), REAL(k_g) + 0.5, bx0, by0, bz0)
          by(i, j, k) = by0
          call getDipole(0, 0.0, REAL(i_g) + 0.5, REAL(j_g) + 0.5, REAL(k_g), bx0, by0, bz0)
          bz(i, j, k) = bz0
        end do
      end do
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

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                  & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
              y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
                  & + species(s)%prtl_tile(ti, tj, tk)%dy(p)
              z_g = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
                  & + species(s)%prtl_tile(ti, tj, tk)%dz(p)
              r_g = sqrt(x_g**2 + y_g**2 + z_g**2)
              if (r_g .le. (psr_radius - 0.5)) then
                species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
              end if
            end do
          end do
        end do
      end do
    end do

  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    integer                       :: i, j, k, i_glob, j_glob, k_glob
    real                          :: supersph_radius_sq
    real, allocatable             :: ex_new(:,:,:), ey_new(:,:,:), ez_new(:,:,:),&
                                   & bx_new(:,:,:), by_new(:,:,:), bz_new(:,:,:)

    allocate(ex_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                  & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                  & fldBoundZ))
    allocate(ey_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                  & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                  & fldBoundZ))
    allocate(ez_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                  & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                  & fldBoundZ))
    allocate(bx_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                  & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                  & fldBoundZ))
    allocate(by_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                  & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                  & fldBoundZ))
    allocate(bz_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                  & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                  & fldBoundZ))

    supersph_radius_sq = (radius * 2)**2

    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        do k = 0, this_meshblock%ptr%sz - 1
          k_glob = k + this_meshblock%ptr%z0
          if ((i_glog - xc_g)**2 + (j_glob - yc_g)**2 + (k_glob - zc_g)**2 .lt. supersph_radius_sq) then

          end if
        end do
      end do
    end do

    deallocate(ex_new, ey_new, ez_new)
    deallocate(bx_new, by_new, bz_new)

  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- auxiliary functions ------------------------------------!
  subroutine getDipole(step, offset, x_g, y_g, z_g,&
                     & obx, oby, obz)
    implicit none
    integer, intent(in) :: step
    real, intent(in)    :: x_g, y_g, z_g, offset
    real, intent(out)   :: obx, oby, obz
    real                :: phase, nx, ny, nz, rr, mux, muy, muz, mu_dot_n

    phase = rotation_omega * step
    nx = x_g - xc_g
    ny = y_g - yc_g
    nz = z_g - zc_g

    rr = sqrt(nx**2 + ny**2 + nz**2)
    nx = nx / rr
    ny = ny / rr
    nz = nz / rr

    rr = 1.0 / rr**3

    mux = psr_radius**3 * sin(psr_angle) * cos(phase + offset)
    muy = psr_radius**3 * sin(psr_angle) * sin(phase + offset)
    muz = psr_radius**3 * cos(psr_angle)

    mu_dot_n = mux * nx + muy * ny + muz * nz

    obx = (3.0 * mux * mu_dot_n - mux) * rr
    oby = (3.0 * muy * mu_dot_n - muy) * rr
    obz = (3.0 * muz * mu_dot_n - muz) * rr
  end subroutine getDipole
  !............................................................!
end module m_userfile
