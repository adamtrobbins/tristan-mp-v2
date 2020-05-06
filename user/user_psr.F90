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
  integer :: fld_geometry
  real :: xc_g, yc_g, zc_g
  real :: psr_angle, psr_period, psr_omega, psr_radius
  real :: weight_mult, e_dot_b_thr
  real :: shell_width

  private :: fld_geometry
  private :: xc_g, yc_g, zc_g, psr_angle, psr_period, psr_omega, psr_radius
  private :: weight_mult, e_dot_b_thr
  private :: shell_width
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
    call getInput('problem', 'e_dot_b_thr', e_dot_b_thr, 0.01)
    call getInput('problem', 'inj_shell', shell_width, 0.75)
    call getInput('problem', 'weight_mult', weight_mult, 1.0)
    call getInput('problem', 'fld_geometry', fld_geometry, 2)

    ! safety check
    if ((psr_angle .le. 1e-2) .or. (psr_angle .ge. 1.0)) then
      psr_angle = psr_angle * M_PI / 180.0
    end if

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
    real                        :: radius
    radius = (dummy1 * 0.5 - x_glob)**2 + (dummy2 * 0.5 - y_glob)**2 + (dummy3 * 0.5 - z_glob)**2 + 1.0
    userSLBload = 10.0 / radius
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
          call getBfield(0, 0.0, REAL(i_glob), REAL(j_glob) + 0.5, REAL(k_glob) + 0.5, bx0, by0, bz0)
          bx(i, j, k) = bx0
          call getBfield(0, 0.0, REAL(i_glob) + 0.5, REAL(j_glob), REAL(k_glob) + 0.5, bx0, by0, bz0)
          by(i, j, k) = by0
          call getBfield(0, 0.0, REAL(i_glob) + 0.5, REAL(j_glob) + 0.5, REAL(k_glob), bx0, by0, bz0)
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
    integer                       :: n_part, n
    real                          :: x_loc, y_loc, z_loc, dx, dy, dz
    integer(kind=2)               :: xi, yi, zi
    real                          :: x_glob, y_glob, z_glob, weight, ppc
    real                          :: ex0, ey0, ez0, bx0, by0, bz0, e_dot_b, b_sqr

    ! inject new particles in the spherical shell
    ppc = 0.5 * ppc0
    n_part = INT((4.0 * M_PI / 3.0) * ((psr_radius + shell_width)**3 - psr_radius**3) * ppc)
    do n = 1, n_part
      call randomPointInSphericalShell(psr_radius, psr_radius + shell_width, x_glob, y_glob, z_glob)
      x_glob = x_glob + xc_g
      y_glob = y_glob + yc_g
      z_glob = z_glob + zc_g

      ! convert global coordinates to local
      call globalToLocalCoords(x_glob, y_glob, z_glob,&
                             & x_loc, y_loc, z_loc, .true.)
      call localToCellBasedCoords(x_loc, y_loc, z_loc,&
                                & xi, yi, zi, dx, dy, dz)
      ! interpolate fields on particle position
      call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
      call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
      b_sqr = bx0**2 + by0**2 + bz0**2
      e_dot_b = abs(ex0 * bx0 + ey0 * by0 + ez0 * bz0)
      weight = weight_mult * B_norm * (e_dot_b / sqrt(b_sqr)) / (abs(unit_ch) * ppc * shell_width)

      if (e_dot_b / b_sqr .gt. e_dot_b_thr) then
        call injectParticleGlobally(1, x_glob, y_glob, z_glob, 0.0, 0.0, 0.0, weight)
        call injectParticleGlobally(2, x_glob, y_glob, z_glob, 0.0, 0.0, 0.0, weight)
      end if
    end do

    ! remove particles falling into the star
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
              r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
              if (r_g .le. (psr_radius - 1.0)) then
                species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
              end if
            end do
          end do
        end do
      end do
    end do
  end subroutine userParticleBoundaryConditions

  subroutine randomPointInSphericalShell(rmin, rmax, x, y, z)
    implicit none
    real, intent(in)    :: rmin, rmax
    real, intent(out)   :: x, y, z
    real                :: TH, X0, Y0, Z0, T, R
    TH = random(dseed) * 2.0 * M_PI
    Z0 = 2.0 * (random(dseed) - 0.5)
    X0 = sqrt(1.0 - Z0**2) * cos(TH)
    Y0 = sqrt(1.0 - Z0**2) * sin(TH)
    T = random(dseed) * (rmax**3 - rmin**3) + rmin**3
    R = T**(1.0 / 3.0)
    x = X0 * R
    y = Y0 * R
    z = Z0 * R
  end subroutine randomPointInSphericalShell

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

    if (updateE_) then
      allocate(ex_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                    & fldBoundZ))
      allocate(ey_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                    & fldBoundZ))
      allocate(ez_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                    & fldBoundZ))
    end if
    if (updateB_) then
      allocate(bx_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                    & fldBoundZ))
      allocate(by_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                    & fldBoundZ))
      allocate(bz_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                    & fldBoundZ))
    end if

    ! update only within this "supersphere"
    supersph_radius_sq = (psr_radius * 2)**2
    ! B fields are set few cells below the E fields
    shift_B = 4.0
    shift_E = 0.0

    scaleEpar = 0.5
    scaleEperp = 0.25
    scaleBperp = 0.5
    scaleBpar = 0.25

    if (updateE_ .or. updateB_) then
      ! saving boundaries into `e/b_new` arrays
      do i = 0, this_meshblock%ptr%sx - 1
        i_glob = i + this_meshblock%ptr%x0
        do j = 0, this_meshblock%ptr%sy - 1
          j_glob = j + this_meshblock%ptr%y0
          do k = 0, this_meshblock%ptr%sz - 1
            k_glob = k + this_meshblock%ptr%z0
            if ((i_glob - xc_g)**2 + (j_glob - yc_g)**2 + (k_glob - zc_g)**2 .lt. supersph_radius_sq) then
              ! ... setting B-field
              if (updateB_) then
                ! setting `B_x`
                rx = REAL(i_glob) - xc_g
                ry = REAL(j_glob) + 0.5 - xc_g
                rz = REAL(k_glob) + 0.5 - xc_g
                rr_sqr = rx**2 + ry**2 + rz**2
                b_int_dot_r = bx(i,j,k) * rx +&
                            & 0.25 * (by(i,j,k) + by(i,j+1,k) + by(i-1,j,k) + by(i-1,j+1,k)) * ry +&
                            & 0.25 * (bz(i,j,k) + bz(i,j,k+1) + bz(i-1,j,k) + bz(i-1,j,k+1)) * rz
                call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

                ! scale = scaleBperp
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                ! bx_new(i, j, k) = (b_dip_dot_r * rx / rr_sqr) +&
                !                 & (b_int_dot_r - b_dip_dot_r) * (rx / rr_sqr) * (1.0 - s)
                ! scale = scaleBpar
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                ! bx_new(i, j, k) = bx_new(i, j, k) + bx_dip - (b_dip_dot_r * rx / rr_sqr) +&
                !                 & ((bx(i, j, k) - b_int_dot_r * rx / rr_sqr) -&
                !                 & (bx_dip - b_dip_dot_r * rx / rr_sqr)) * (1.0 - s)
                scale = scaleBperp
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                bx_new(i, j, k) = (b_int_dot_r - b_dip_dot_r) * (rx / rr_sqr) * (1.0 - s)
                scale = scaleBpar
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                bx_new(i, j, k) = bx_new(i, j, k) + bx_dip +&
                                & ((bx(i, j, k) - b_int_dot_r * rx / rr_sqr) -&
                                & (bx_dip - b_dip_dot_r * rx / rr_sqr)) * (1.0 - s)


                ! setting `B_y`
                rx = REAL(i_glob) + 0.5 - xc_g
                ry = REAL(j_glob) - xc_g
                rz = REAL(k_glob) + 0.5 - xc_g
                rr_sqr = rx**2 + ry**2 + rz**2
                b_int_dot_r = 0.25 * (bx(i,j,k) + bx(i+1,j,k) + bx(i,j-1,k) + bx(i+1,j-1,k)) * rx +&
                            & by(i,j,k) * ry +&
                            & 0.25 * (bz(i,j,k) + bz(i,j,k+1) + bz(i,j-1,k) + bz(i,j-1,k+1)) * rz
                call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

                ! scale = scaleBperp
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                ! by_new(i, j, k) = (b_dip_dot_r * ry / rr_sqr) +&
                !                 & (b_int_dot_r - b_dip_dot_r) * (ry / rr_sqr) * (1.0 - s)
                ! scale = scaleBpar
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                ! by_new(i, j, k) = by_new(i, j, k) + by_dip - (b_dip_dot_r * ry / rr_sqr) +&
                !                 & ((by(i, j, k) - b_int_dot_r * ry / rr_sqr) -&
                !                 & (by_dip - b_dip_dot_r * ry / rr_sqr)) * (1.0 - s)
                scale = scaleBperp
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                by_new(i, j, k) = (b_int_dot_r - b_dip_dot_r) * (ry / rr_sqr) * (1.0 - s)
                scale = scaleBpar
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                by_new(i, j, k) = by_new(i, j, k) + by_dip +&
                                & ((by(i, j, k) - b_int_dot_r * ry / rr_sqr) -&
                                & (by_dip - b_dip_dot_r * ry / rr_sqr)) * (1.0 - s)

                ! setting `B_z`
                rx = REAL(i_glob) + 0.5 - xc_g
                ry = REAL(j_glob) + 0.5 - xc_g
                rz = REAL(k_glob) - xc_g
                rr_sqr = rx**2 + ry**2 + rz**2
                b_int_dot_r = 0.25 * (bx(i,j,k) + bx(i+1,j,k) + bx(i,j,k-1) + bx(i+1,j,k-1)) * rx +&
                            & 0.25 * (by(i,j,k) + by(i,j+1,k) + by(i,j,k-1) + by(i,j+1,k-1)) * ry +&
                            & bz(i,j,k) * rz
                call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

                ! scale = scaleBperp
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                ! bz_new(i, j, k) = (b_dip_dot_r * rz / rr_sqr) +&
                !                 & (b_int_dot_r - b_dip_dot_r) * (rz / rr_sqr) * (1.0 - s)
                ! scale = scaleBpar
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                ! bz_new(i, j, k) = bz_new(i, j, k) + bz_dip - (b_dip_dot_r * rz / rr_sqr) +&
                !                 & ((bz(i, j, k) - b_int_dot_r * rz / rr_sqr) -&
                !                 & (bz_dip - b_dip_dot_r * rz / rr_sqr)) * (1.0 - s)
                scale = scaleBperp
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                bz_new(i, j, k) = (b_int_dot_r - b_dip_dot_r) * (rz / rr_sqr) * (1.0 - s)
                scale = scaleBpar
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_B) / scale)
                bz_new(i, j, k) = bz_new(i, j, k) + bz_dip +&
                                & ((bz(i, j, k) - b_int_dot_r * rz / rr_sqr) -&
                                & (bz_dip - b_dip_dot_r * rz / rr_sqr)) * (1.0 - s)
              end if
              ! ... setting E-field
              if (updateE_) then
                ! setting `E_x`
                rx = REAL(i_glob) + 0.5 - xc_g
                ry = REAL(j_glob) - xc_g
                rz = REAL(k_glob) - xc_g
                rr_sqr = rx**2 + ry**2 + rz**2
                e_int_dot_r = ex(i,j,k) * rx +&
                            & 0.25 * (ey(i,j,k) + ey(i+1,j,k) + ey(i,j-1,k) + ey(i+1,j-1,k)) * ry +&
                            & 0.25 * (ez(i,j,k) + ez(i+1,j,k) + ez(i,j,k-1) + ez(i+1,j,k-1)) * rz
                call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                vx = -psr_omega * ry
                vy = psr_omega * rx
                vz = 0.0
                ex_dip = -(vy * bz_dip - vz * by_dip) * CCINV
                ey_dip = (vx * bz_dip - vz * bx_dip) * CCINV
                ez_dip = -(vx * by_dip - vy * bx_dip) * CCINV
                e_dip_dot_r = ex_dip * rx + ey_dip * ry + ez_dip * rz

                ! scale = scaleEpar
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ! ex_new(i, j, k) = (ex_dip - e_dip_dot_r * rx / rr_sqr) +&
                !                 & ((ex(i, j, k) - e_int_dot_r * rx / rr_sqr) -&
                !                 & (ex_dip - e_dip_dot_r * rx / rr_sqr)) * (1.0 - s)
                ! scale = scaleEperp
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ! ex_new(i, j, k) = ex_new(i, j, k) + (e_dip_dot_r * rx / rr_sqr) +&
                !                 & (e_int_dot_r - e_dip_dot_r) * (rx / rr_sqr) * (1.0 - s)
                scale = scaleEpar
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ex_new(i, j, k) = ex_dip +&
                                & ((ex(i, j, k) - e_int_dot_r * rx / rr_sqr) -&
                                & (ex_dip - e_dip_dot_r * rx / rr_sqr)) * (1.0 - s)
                scale = scaleEperp
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ex_new(i, j, k) = ex_new(i, j, k) +&
                                & (e_int_dot_r - e_dip_dot_r) * (rx / rr_sqr) * (1.0 - s)

                ! setting `E_y`
                rx = REAL(i_glob) - xc_g
                ry = REAL(j_glob) + 0.5 - xc_g
                rz = REAL(k_glob) - xc_g
                rr_sqr = rx**2 + ry**2 + rz**2
                e_int_dot_r = 0.25 * (ex(i,j,k) + ex(i,j+1,k) + ex(i-1,j,k) + ex(i-1,j+1,k)) * rx +&
                            & ey(i,j,k) * ry +&
                            & 0.25 * (ez(i,j,k) + ez(i,j+1,k) + ez(i,j,k-1) + ez(i,j+1,k-1)) * rz
                call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                vx = -psr_omega * ry
                vy = psr_omega * rx
                vz = 0.0
                ex_dip = -(vy * bz_dip - vz * by_dip) * CCINV
                ey_dip = (vx * bz_dip - vz * bx_dip) * CCINV
                ez_dip = -(vx * by_dip - vy * bx_dip) * CCINV
                e_dip_dot_r = ex_dip * rx + ey_dip * ry + ez_dip * rz

                ! scale = scaleEpar
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ! ey_new(i, j, k) = (ey_dip - e_dip_dot_r * ry / rr_sqr) +&
                !                 & ((ey(i, j, k) - e_int_dot_r * ry / rr_sqr) -&
                !                 & (ey_dip - e_dip_dot_r * ry / rr_sqr)) * (1.0 - s)
                ! scale = scaleEperp
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ! ey_new(i, j, k) = ey_new(i, j, k) + (e_dip_dot_r * ry / rr_sqr) +&
                !                 & (e_int_dot_r - e_dip_dot_r) * (ry / rr_sqr) * (1.0 - s)
                scale = scaleEpar
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ey_new(i, j, k) = ey_dip +&
                                & ((ey(i, j, k) - e_int_dot_r * ry / rr_sqr) -&
                                & (ey_dip - e_dip_dot_r * ry / rr_sqr)) * (1.0 - s)
                scale = scaleEperp
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ey_new(i, j, k) = ey_new(i, j, k) +&
                                & (e_int_dot_r - e_dip_dot_r) * (ry / rr_sqr) * (1.0 - s)

                ! setting `E_z`
                rx = REAL(i_glob) - xc_g
                ry = REAL(j_glob) - xc_g
                rz = REAL(k_glob) + 0.5 - xc_g
                rr_sqr = rx**2 + ry**2 + rz**2
                e_int_dot_r = 0.25 * (ex(i,j,k) + ex(i,j,k+1) + ex(i-1,j,k) + ex(i-1,j,k+1)) * rx +&
                            & 0.25 * (ey(i,j,k) + ey(i,j-1,k) + ey(i,j,k+1) + ey(i,j-1,k+1)) * ry +&
                            & ez(i,j,k) * rz
                call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                vx = -psr_omega * ry
                vy = psr_omega * rx
                vz = 0.0
                ex_dip = -(vy * bz_dip - vz * by_dip) * CCINV
                ey_dip = (vx * bz_dip - vz * bx_dip) * CCINV
                ez_dip = -(vx * by_dip - vy * bx_dip) * CCINV
                e_dip_dot_r = ex_dip * rx + ey_dip * ry + ez_dip * rz

                ! scale = scaleEpar
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ! ez_new(i, j, k) = (ez_dip - e_dip_dot_r * rz / rr_sqr) +&
                !                 & ((ez(i, j, k) - e_int_dot_r * rz / rr_sqr) -&
                !                 & (ez_dip - e_dip_dot_r * rz / rr_sqr)) * (1.0 - s)
                ! scale = scaleEperp
                ! s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ! ez_new(i, j, k) = ez_new(i, j, k) + (e_dip_dot_r * rz / rr_sqr) +&
                !                 & (e_int_dot_r - e_dip_dot_r) * (rz / rr_sqr) * (1.0 - s)
                scale = scaleEpar
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ez_new(i, j, k) = ez_dip +&
                                & ((ez(i, j, k) - e_int_dot_r * rz / rr_sqr) -&
                                & (ez_dip - e_dip_dot_r * rz / rr_sqr)) * (1.0 - s)
                scale = scaleEperp
                s = shape(sqrt(rr_sqr) / scale, (psr_radius - shift_E) / scale)
                ez_new(i, j, k) = ez_new(i, j, k) +&
                                & (e_int_dot_r - e_dip_dot_r) * (rz / rr_sqr) * (1.0 - s)
              end if
            end if
          end do
        end do
      end do
      ! copying to real array
      do i = 0, this_meshblock%ptr%sx - 1
        i_glob = i + this_meshblock%ptr%x0
        do j = 0, this_meshblock%ptr%sy - 1
          j_glob = j + this_meshblock%ptr%y0
          do k = 0, this_meshblock%ptr%sz - 1
            k_glob = k + this_meshblock%ptr%z0
            if ((i_glob - xc_g)**2 + (j_glob - yc_g)**2 + (k_glob - zc_g)**2 .lt. supersph_radius_sq) then
              if (updateE_) then
                ex(i, j, k) = ex_new(i, j, k)
                ey(i, j, k) = ey_new(i, j, k)
                ez(i, j, k) = ez_new(i, j, k)
              end if
              if (updateB_) then
                bx(i, j, k) = bx_new(i, j, k)
                by(i, j, k) = by_new(i, j, k)
                bz(i, j, k) = bz_new(i, j, k)
              end if
            end if
          end do
        end do
      end do
    end if

    if (updateE_) deallocate(ex_new, ey_new, ez_new)
    if (updateB_) deallocate(bx_new, by_new, bz_new)
  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- auxiliary functions ------------------------------------!
  subroutine getBfield(step, offset, x_g, y_g, z_g,&
                     & obx, oby, obz)
    integer, intent(in) :: step
    real, intent(in)    :: x_g, y_g, z_g, offset
    real, intent(out)   :: obx, oby, obz
    if (fld_geometry .eq. 1) then
      call getMonopole(step, offset, x_g, y_g, z_g, obx, oby, obz)
    else if (fld_geometry .eq. 2) then
      call getDipole(step, offset, x_g, y_g, z_g, obx, oby, obz)
    else
      print *, "Something went wrong in `usr_psr`."
      stop
    end if
  end subroutine getBfield

  subroutine getDipole(step, offset, x_g, y_g, z_g,&
                     & obx, oby, obz)
    implicit none
    integer, intent(in) :: step
    real, intent(in)    :: x_g, y_g, z_g, offset
    real, intent(out)   :: obx, oby, obz
    real                :: phase, nx, ny, nz, rr, mux, muy, muz, mu_dot_n

    phase = psr_omega * step
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

    obx = (3.0 * nx * mu_dot_n - mux) * rr
    oby = (3.0 * ny * mu_dot_n - muy) * rr
    obz = (3.0 * nz * mu_dot_n - muz) * rr
  end subroutine getDipole

  subroutine getMonopole(step, offset, x_g, y_g, z_g,&
                       & obx, oby, obz)
    implicit none
    integer, intent(in) :: step
    real, intent(in)    :: x_g, y_g, z_g, offset
    real, intent(out)   :: obx, oby, obz
    real                :: nx, ny, nz, rr
    nx = x_g - xc_g
    ny = y_g - yc_g
    nz = z_g - zc_g

    rr = sqrt(nx**2 + ny**2 + nz**2)
    rr = 1.0 / rr**3

    obx = psr_radius**2 * nx * rr
    oby = psr_radius**2 * ny * rr
    obz = psr_radius**2 * nz * rr
  end subroutine getMonopole

  real function shape(rad, rad0)
    implicit none
    real, intent(in)  :: rad, rad0
    real              :: del
    del = 1.0
    shape = 0.5 * (1.0 - tanh((rad - rad0) / del))
  end function shape

  !............................................................!
end module m_userfile
