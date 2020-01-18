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
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
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
    real :: kx, ky, kz
    real :: ex_norm, ey_norm, ez_norm, exyz_norm
    real :: bx_norm, by_norm, bz_norm, bxyz_norm

    ex(:,:,:) = -1; ey(:,:,:) = -1; ez(:,:,:) = -1
    bx(:,:,:) = -1; by(:,:,:) = -1; bz(:,:,:) = -1

    kx = 5; ky = 2; kz = 2
    kx = kx * 2 * M_PI / global_mesh%sx
    ky = ky * 2 * M_PI / global_mesh%sy
    kz = kz * 2 * M_PI / global_mesh%sz

    ex_norm = 0; ey_norm = 2; ez_norm = -2;
    bx_norm = -8; by_norm = 10; bz_norm = 10;

    exyz_norm = sqrt(ex_norm**2 + ey_norm**2 + ez_norm**2)
    ex_norm = ex_norm / exyz_norm
    ey_norm = ey_norm / exyz_norm
    ez_norm = ez_norm / exyz_norm

    bxyz_norm = sqrt(bx_norm**2 + by_norm**2 + bz_norm**2)
    bx_norm = bx_norm / bxyz_norm
    by_norm = by_norm / bxyz_norm
    bz_norm = bz_norm / bxyz_norm

    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        do k = 0, this_meshblock%ptr%sz - 1
          k_glob = k + this_meshblock%ptr%z0
          ex(i, j, k) = ex_norm * sin((i_glob) * kx       + (j_glob - 0.5) * ky + (k_glob - 0.5) * kz)
          ey(i, j, k) = ey_norm * sin((i_glob - 0.5) * kx + (j_glob) * ky       + (k_glob - 0.5) * kz)
          ez(i, j, k) = ez_norm * sin((i_glob - 0.5) * kx + (j_glob - 0.5) * ky + (k_glob) * kz)
          bx(i, j, k) = bx_norm * sin((i_glob - 0.5) * kx + (j_glob) * ky       + (k_glob) * kz)
          by(i, j, k) = by_norm * sin((i_glob) * kx       + (j_glob - 0.5) * ky + (k_glob) * kz)
          bz(i, j, k) = bz_norm * sin((i_glob) * kx       + (j_glob) * ky       + (k_glob - 0.5) * kz)
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
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userFieldBoundaryConditions
  !............................................................!
end module m_userfile
