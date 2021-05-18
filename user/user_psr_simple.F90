#include "../src/defs.F90"

! Configuration for this userfile:
! ```
!   $ python configure.py -3d --user=user_psr_simple -slb --gca=2 --radiation=sync -absorb
! ```

module m_userfile
  use m_globalnamespace
  use m_qednamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  use m_helpers
  #ifdef USROUTPUT
    use m_writeusroutput
  #endif
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private     :: xc_g, yc_g, zc_g, psr_bstar, cooling_on
  real, private     :: psr_angle, psr_period, psr_omega, psr_radius
  real, private     :: inj_mult
  real, private     :: shell_width, prtl_kick, rmin_dr
  real, private     :: sigma_nGJ, nGJ, inj_dr
  real, private     :: sigGJ_limiter, edotb_thr_closed
  real, private     :: omegaB0, psr_rlc
  #ifdef GCA
    real, private     :: psr_gca_enforce_rad
  #endif

  real, private     :: gammarad_over_sigma_LC, gammarad_dummy
  real, private     :: eph_at_LC, gammac_dummy
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution, randomPointInSphericalShell, distance3D_sq
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none

    call getInput('problem', 'psr_radius', psr_radius)
    call getInput('problem', 'psr_angle', psr_angle)
    psr_angle = psr_angle * M_PI / 180.0
    call getInput('problem', 'psr_period', psr_period)
    psr_omega = 2.0 * M_PI / psr_period
    call getInput('problem', 'psr_bstar', psr_bstar)
    call getInput('problem', 'prtl_kick', prtl_kick, 0.0)

    call getInput('problem', 'inj_mult', inj_mult, 1.0)
    call getInput('problem', 'sigGJ_limiter', sigGJ_limiter, 1000.0)
    call getInput('problem', 'edotb_thr_closed', edotb_thr_closed, 0.002)

    call getInput('problem', 'cooling_on', cooling_on, 0.0)

    #ifdef GCA
      call getInput('problem', 'gca_radius', psr_gca_enforce_rad)
    #endif

    psr_rlc = CC / psr_omega
    omegaB0 = sqrt(sigma) * CC / c_omp

    #ifdef RADIATION
      ! gamma_rad / sigma_LC for the field near LC
      call getInput('problem', 'grad_sigma_LC', gammarad_over_sigma_LC, 0.5)
      gammarad_dummy = 0.5 * gammarad_over_sigma_LC * (omegaB0 * psr_radius / CC) * sqrt(psr_radius / psr_rlc) * psr_bstar

      ! energy of the photon (in me c^2) radiated in the field of LC by a particle with gamma ~ sigmaLC
      call getInput('problem', 'eph_at_LC', eph_at_LC, 100.0)
      gammac_dummy = 0.5 * eph_at_LC**(-0.5) * psr_bstar**1.5 * (omegaB0 * psr_radius / CC) * (psr_radius / psr_rlc)**(3.5)
    #endif

    call getInput('problem', 'rmin_dr', rmin_dr, 1.0)
    call getInput('problem', 'shell_width', shell_width, 2.0)
    call getInput('problem', 'inj_dr', inj_dr, 1.0)

    xc_g = 0.5 * global_mesh%sx
    yc_g = 0.5 * global_mesh%sy
    zc_g = 0.5 * global_mesh%sz

    global_usr_variable_1 = psr_radius
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
    real                        :: radius2, psrrad
    psrrad = global_usr_variable_1
    radius2 = (dummy1 * 0.5 - x_glob)**2 + (dummy2 * 0.5 - y_glob)**2 + (dummy3 * 0.5 - z_glob)**2 + 1.0
    userSLBload = psrrad**2 / radius2 + exp(-(dummy3 * 0.5 - z_glob)**2 / (psrrad * 0.5)**2)
    if (radius2 .lt. psrrad**2) then
      userSLBload = 1.0 / exp((psrrad**2 - radius2) / psrrad**2)
    end if
    return
  end function

  subroutine userInitParticles()
    implicit none
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution
    #ifdef RADIATION
      ! redefine `gamma_syn`, `emit_gamma_syn`
      cool_gamma_syn = gammarad_dummy
      emit_gamma_syn = gammac_dummy
    #endif
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real    :: bx0, by0, bz0, x_, y_, z_
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    do i = 0, this_meshblock%ptr%sx - 1
      i_glob = i + this_meshblock%ptr%x0
      do j = 0, this_meshblock%ptr%sy - 1
        j_glob = j + this_meshblock%ptr%y0
        do k = 0, this_meshblock%ptr%sz - 1
          k_glob = k + this_meshblock%ptr%z0

          x_ = REAL(i_glob);  y_ = REAL(j_glob) + 0.5;  z_ = REAL(k_glob) + 0.5
          call getBfield(0, 0.0, x_, y_, z_, bx0, by0, bz0)
          bx(i, j, k) = bx0

          x_ = REAL(i_glob) + 0.5;  y_ = REAL(j_glob);  z_ = REAL(k_glob) + 0.5
          call getBfield(0, 0.0, x_, y_, z_, bx0, by0, bz0)
          by(i, j, k) = by0

          x_ = REAL(i_glob) + 0.5;  y_ = REAL(j_glob) + 0.5;  z_ = REAL(k_glob)
          call getBfield(0, 0.0, x_, y_, z_, bx0, by0, bz0)
          bz(i, j, k) = bz0
        end do
      end do
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
  end subroutine userCurrentDeposit

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

  #ifdef GCA
    logical function userEnforceGCA(xi, yi, zi, dx, dy, dz, u, v, w, weight)
      implicit none
      integer(kind=2), intent(in), optional   :: xi, yi, zi
      real, intent(in), optional              :: dx, dy, dz, u, v, w
      real, intent(in), optional              :: weight
      real :: rr
      rr = sqrt((real(xi + this_meshblock%ptr%x0) + dx - xc_g)**2 +&
              & (real(yi + this_meshblock%ptr%y0) + dy - yc_g)**2 +&
              & (real(zi + this_meshblock%ptr%z0) + dz - zc_g)**2)
      userEnforceGCA = ((rr .lt. psr_radius + psr_gca_enforce_rad) .and. (psr_gca_enforce_rad .ne. 0))
    end function userEnforceGCA
  #endif
  !............................................................!

  !--- boundaries ---------------------------------------------!
  #include "user_psr_simple.F08"

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

  real function distance3D_sq(xA, yA, zA, xB, yB, zB)
    real, intent(in) :: xA, yA, zA, xB, yB, zB
    distance3D_sq = (xA - xB)**2 + (yA - yB)**2 + (zA - zB)**2
  end function distance3D_sq

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
    real                          :: rr, x_, y_, z_, rlimit
    logical                       :: dummy_log

    ! update only within this "supersphere"
    supersph_radius_sq = (psr_radius * 2)**2
    ! test if meshblock "touches" the supersphere
    x_ = max(REAL(this_meshblock%ptr%x0 - NGHOST), min(REAL(this_meshblock%ptr%x0 + this_meshblock%ptr%sx - 1 + NGHOST), xc_g))
    y_ = max(REAL(this_meshblock%ptr%y0 - NGHOST), min(REAL(this_meshblock%ptr%y0 + this_meshblock%ptr%sy - 1 + NGHOST), yc_g))
    z_ = max(REAL(this_meshblock%ptr%z0 - NGHOST), min(REAL(this_meshblock%ptr%z0 + this_meshblock%ptr%sz - 1 + NGHOST), zc_g))
    dummy_log = (distance3D_sq(x_, y_, z_, xc_g, yc_g, zc_g) .lt. supersph_radius_sq)

    if (dummy_log) then

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
                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
        allocate(ey_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
        allocate(ez_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      end if
      if (updateB_) then
        allocate(bx_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
        allocate(by_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
        allocate(bz_new(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                      & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      end if

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
                  ry = REAL(j_glob) + 0.5 - yc_g
                  rz = REAL(k_glob) + 0.5 - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  b_int_dot_r = bx(i,j,k) * rx +&
                              & 0.25 * (by(i,j,k) + by(i,j+1,k) + by(i-1,j,k) + by(i-1,j+1,k)) * ry +&
                              & 0.25 * (bz(i,j,k) + bz(i,j,k+1) + bz(i-1,j,k) + bz(i-1,j,k+1)) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

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
                  ry = REAL(j_glob) - yc_g
                  rz = REAL(k_glob) + 0.5 - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  b_int_dot_r = 0.25 * (bx(i,j,k) + bx(i+1,j,k) + bx(i,j-1,k) + bx(i+1,j-1,k)) * rx +&
                              & by(i,j,k) * ry +&
                              & 0.25 * (bz(i,j,k) + bz(i,j,k+1) + bz(i,j-1,k) + bz(i,j-1,k+1)) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

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
                  ry = REAL(j_glob) + 0.5 - yc_g
                  rz = REAL(k_glob) - zc_g
                  rr_sqr = rx**2 + ry**2 + rz**2
                  b_int_dot_r = 0.25 * (bx(i,j,k) + bx(i+1,j,k) + bx(i,j,k-1) + bx(i+1,j,k-1)) * rx +&
                              & 0.25 * (by(i,j,k) + by(i,j+1,k) + by(i,j,k-1) + by(i,j+1,k-1)) * ry +&
                              & bz(i,j,k) * rz
                  call getBfield(step, 0.0, rx + xc_g, ry + yc_g, rz + zc_g, bx_dip, by_dip, bz_dip)
                  b_dip_dot_r = bx_dip * rx + by_dip * ry + bz_dip * rz

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
                  ry = REAL(j_glob) - yc_g
                  rz = REAL(k_glob) - zc_g
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
                  ry = REAL(j_glob) + 0.5 - yc_g
                  rz = REAL(k_glob) - zc_g
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
                  ry = REAL(j_glob) - yc_g
                  rz = REAL(k_glob) + 0.5 - zc_g
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

    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- auxiliary functions ------------------------------------!
  subroutine getBfield(step, offset, x_g, y_g, z_g,&
                     & obx, oby, obz)
    integer, intent(in) :: step
    real, intent(in)    :: x_g, y_g, z_g, offset
    real, intent(out)   :: obx, oby, obz
    call getDipole(step, offset, x_g, y_g, z_g, obx, oby, obz)
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

    rr = 1.0 / sqrt(nx**2 + ny**2 + nz**2)
    nx = nx * rr
    ny = ny * rr
    nz = nz * rr
    rr = rr**3

    mux = psr_radius**3 * sin(psr_angle) * cos(phase + offset)
    muy = psr_radius**3 * sin(psr_angle) * sin(phase + offset)
    muz = psr_radius**3 * cos(psr_angle)

    mu_dot_n = mux * nx + muy * ny + muz * nz

    obx = psr_bstar * (3.0 * nx * mu_dot_n - mux) * rr
    oby = psr_bstar * (3.0 * ny * mu_dot_n - muy) * rr
    obz = psr_bstar * (3.0 * nz * mu_dot_n - muz) * rr
  end subroutine getDipole

  real function shape(rad, rad0)
    implicit none
    real, intent(in)  :: rad, rad0
    real              :: del
    del = 1.0
    shape = 0.5 * (1.0 - tanh((rad - rad0) / del))
  end function shape

  !............................................................!

  !--- user-specific output -----------------------------------!
  #ifdef USROUTPUT
    subroutine userOutput(step)
      implicit none
      integer, optional, intent(in) :: step
      integer                       :: root_rank = 0
      real, allocatable             :: r_bins(:)
      real                          :: dr, x_glob, y_glob, z_glob, r_glob, fr_factor
      real                          :: dummy_x, dummy_y, dummy_z, dummy1, dummy2
      real, allocatable             :: sum_ExBr_f(:), sum_f(:), sum_ExBr_f_global(:), sum_f_global(:)
      integer                       :: ri, rnum = 50, i, j, k, ierr

      allocate(r_bins(rnum))
      allocate(sum_ExBr_f(rnum), sum_f(rnum))
      allocate(sum_ExBr_f_global(rnum), sum_f_global(rnum))
      sum_ExBr_f(:) = 0.0
      sum_f(:) = 0.0

      dr = (REAL(global_mesh%sx) * 0.5 - psr_radius) / REAL(rnum)
      do ri = 1, rnum
        r_bins(ri) = psr_radius + ri * dr
      end do

      do i = 0, this_meshblock%ptr%sx - 1
        x_glob = REAL(i + this_meshblock%ptr%x0)
        do j = 0, this_meshblock%ptr%sy - 1
          y_glob = REAL(j + this_meshblock%ptr%y0)
          do k = 0, this_meshblock%ptr%sz - 1
            z_glob = REAL(k + this_meshblock%ptr%z0)
            r_glob = sqrt((x_glob - xc_g)**2 + (y_glob - yc_g)**2 + (z_glob - zc_g)**2)
            ! compute ExB_r
            dummy_x = -(ez(i,j,k) * by(i,j,k)) + ey(i,j,k) * bz(i,j,k)
            dummy_y = ez(i,j,k) * bx(i,j,k) - ex(i,j,k) * bz(i,j,k)
            dummy_z = -(ey(i,j,k) * bx(i,j,k)) + ex(i,j,k) * by(i,j,k)
            if (r_glob .gt. psr_radius / 2.0) then
              dummy1 = (dummy_x * (x_glob - xc_g) +&
                     & dummy_y * (y_glob - yc_g) +&
                     & dummy_z * (z_glob - zc_g)) / r_glob
            else
              dummy1 = 0.0
            end if
            do ri = 1, rnum
              fr_factor = exp(-(r_glob - r_bins(ri))**2 / (dr * 0.5)**2)
              sum_ExBr_f(ri) = sum_ExBr_f(ri) + dummy1 * fr_factor 
              sum_f(ri) = sum_f(ri) + fr_factor
            end do
          end do
        end do
      end do

      call MPI_REDUCE(sum_ExBr_f, sum_ExBr_f_global, rnum, MPI_REAL, MPI_SUM, root_rank, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(sum_f, sum_f_global, rnum, MPI_REAL, MPI_SUM, root_rank, MPI_COMM_WORLD, ierr)

      if (mpi_rank .eq. root_rank) then
        ! normalizations
        sum_ExBr_f_global(:) = sum_ExBr_f_global(:) * r_bins(:)**2 * CC * B_norm**2 / sum_f_global(:)
        call writeUsrOutputTimestep(step)
        call writeUsrOutputArray('r_bins', r_bins)
        call writeUsrOutputArray('ExB_flux', sum_ExBr_f_global)
        call writeUsrOutputEnd()
      end if
    end subroutine userOutput

    logical function userExcludeParticles(s, ti, tj, tk, p)
      implicit none
      integer, intent(in)       :: s, ti, tj, tk, p
      real                      :: xx, yy, zz, rr
      
      xx = REAL(this_meshblock%ptr%x0 + species(s)%prtl_tile(ti, tj, tk)%xi(p))
      xx = xx + species(s)%prtl_tile(ti, tj, tk)%dx(p)
      xx = xx - REAL(global_mesh%sx) * 0.5

      yy = REAL(this_meshblock%ptr%y0 + species(s)%prtl_tile(ti, tj, tk)%yi(p))
      yy = yy + species(s)%prtl_tile(ti, tj, tk)%dy(p)
      yy = yy - REAL(global_mesh%sy) * 0.5

      zz = REAL(this_meshblock%ptr%z0 + species(s)%prtl_tile(ti, tj, tk)%zi(p))
      zz = zz + species(s)%prtl_tile(ti, tj, tk)%dz(p)
      zz = zz - REAL(global_mesh%sz) * 0.5

      rr = sqrt(xx**2 + yy**2 + zz**2)
      
      userExcludeParticles = ((rr .gt. global_usr_variable_1 + 10) .and. (rr .lt. REAL(global_mesh%sx) * 0.5 - 100))
    end function userExcludeParticles
  #endif
  !............................................................!
end module m_userfile
