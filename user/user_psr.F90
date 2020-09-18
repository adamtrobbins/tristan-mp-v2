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

  !--- PRIVATE variables -----------------------------------------!
  integer, private  :: fld_geometry, inj_method, e_par_method
  real, private     :: xc_g, yc_g, zc_g, psr_spinupT
  real, private     :: psr_angle, psr_period, psr_omega, psr_omega0, psr_radius
  real, private     :: inj_mult, e_thr
  real, private     :: shell_width, prtl_kick, rmin_dr, e_dr
  real, private     :: sigma_nGJ, nGJ, inj_dr
  real, private     :: nGJ_limiter, sigGJ_limiter, jdotb_limiter
  #ifdef GCA
    real, private     :: psr_gca_enforce_rad
  #endif
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution, getEparAt, randomPointInSphericalShell,&
           & getDeltaErAt
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    ! B-field geometry: 1 = monopole, 2 = dipole
    call getInput('problem', 'fld_geometry', fld_geometry)
    call getInput('problem', 'psr_radius', psr_radius)
    call getInput('problem', 'psr_angle', psr_angle)
    call getInput('problem', 'psr_period', psr_period)
    call getInput('problem', 'psr_spinupT', psr_spinupT, 0.0)

    ! remove particles which fall below `radius - rmin_dr`
    call getInput('problem', 'rmin_dr', rmin_dr)

    ! injection method:
    !     0 = no injection at all
    !     1 = inject with 0 velocity in the shell with weight ~ E.B
    !     2 = inject with a kick and constant weight
    call getInput('problem', 'inj_method', inj_method)
    call getInput('problem', 'inj_shell', shell_width)

    ! for method #1/#2
    if (inj_method .ge. 1) then
      ! probe slightl above the injection point
      call getInput('problem', 'e_dr', e_dr)
      ! inject slightly above the star
      call getInput('problem', 'inj_dr', inj_dr)
      call getInput('problem', 'inj_mult', inj_mult)

      call getInput('problem', 'nGJ_limiter', nGJ_limiter)
    end if

    ! for method #1
    if (inj_method .eq. 1) then
      call getInput('problem', 'e_thr', e_thr)
      call getInput('problem', 'e_par_method', e_par_method)
    end if

    ! for method #2
    if (inj_method .eq. 2) then
      call getInput('problem', 'sigGJ_limiter', sigGJ_limiter)
      call getInput('problem', 'jdotb_limiter', jdotb_limiter)
      call getInput('problem', 'prtl_kick', prtl_kick)
    end if

    #ifdef GCA
      call getInput('problem', 'gca_radius', psr_gca_enforce_rad)
    #endif

    ! safety check
    if ((psr_angle .le. 1e-2) .or. (psr_angle .ge. 1.0)) then
      psr_angle = psr_angle * M_PI / 180.0
    end if

    psr_omega0 = 2.0 * M_PI / psr_period

    if (psr_spinupT .eq. 0.0) then
      psr_omega = psr_omega0
    else
      psr_omega = 0.0
    end if

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

  #define PSRRADIUS 40

  function userSLBload(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    real                        :: radius2
    radius2 = (dummy1 * 0.5 - x_glob)**2 + (dummy2 * 0.5 - y_glob)**2 + (dummy3 * 0.5 - z_glob)**2 + 1.0
    userSLBload = PSRRADIUS**2 / radius2
    if (radius2 .lt. PSRRADIUS**2) then
      userSLBload = 1.0 / exp((PSRRADIUS**2 - radius2) / PSRRADIUS**2)
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
      userEnforceGCA = (rr .lt. psr_radius + psr_gca_enforce_rad)
    end function userEnforceGCA
  #endif
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
    integer                       :: s, ti, tj, tk, p
    real                          :: x_g, y_g, z_g, r_g
    integer                       :: n_part, n, sign
    real                          :: x_loc, y_loc, z_loc, dx, dy, dz
    integer(kind=2)               :: xi, yi, zi, xi_eb, yi_eb, zi_eb
    real                          :: x_glob, y_glob, z_glob, weight, ppc, dens, sig
    real                          :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
    real                          :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma
    real                          :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
    logical                       :: dummy_flag
    #ifdef GCA
      real                          :: dummy_, vE_x, vE_y, vE_z
      real                          :: e0_SQR, b0_SQR
    #endif

    nGJ = 2 * psr_omega0 * B_norm / (CC * abs(unit_ch))
    sigma_nGJ = sigma * ppc0 / nGJ

    if (inj_method .eq. 1) then
      ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ! inject particles with `w ~ E.B` at rest
      ppc = 0.5 * ppc0
      if (nGJ_limiter .ne. 0) then
        ! compute charge density and write to `lg_arr`
        call computeDensity(1, reset=.true., ds=0, charge=.true.)
        call computeDensity(2, reset=.false., ds=0, charge=.true.)
      end if
      n_part = INT((4.0 * M_PI / 3.0) * ((psr_radius + shell_width + inj_dr)**3 - (psr_radius + inj_dr)**3) * ppc)
      do n = 1, n_part
        ! inject slightly above the star
        call randomPointInSphericalShell(psr_radius + inj_dr, psr_radius + inj_dr + shell_width, x_glob, y_glob, z_glob)
        rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
        nx = x_glob / rr
        ny = y_glob / rr
        nz = z_glob / rr
        x_glob = x_glob + xc_g
        y_glob = y_glob + yc_g
        z_glob = z_glob + zc_g
        ! screen E-parallel
        if (e_par_method .eq. 1) then
          ! ... measured `e_dr` cells above the injection point
          call getEparAt(e_dot_b, b_sqr, x_glob, y_glob, z_glob, dummy_flag)
          b_sqr = b_sqr + TINYFLD
          e_b_scale = abs(e_dot_b) / b_sqr
          weight = inj_mult * (B_norm * abs(e_dot_b) / sqrt(b_sqr)) / (unit_ch * ppc * shell_width)
        else if (e_par_method .eq. 2) then
          ! ... measured `e_dr` cells above the injection point
          call getDeltaErAt(delta_er, b_sqr, x_glob, y_glob, z_glob, dummy_flag)
          e_b_scale = abs(delta_er) / sqrt(b_sqr)
          weight = inj_mult * (B_norm * abs(delta_er)) / (unit_ch * ppc * shell_width)
        end if

        ! if point is contained in MPI block
        if (dummy_flag) then
          if (e_b_scale .gt. e_thr) then
            if (nGJ_limiter .ne. 0) then
              call globalToLocalCoords(x_glob, y_glob, z_glob,&
                                     & x_loc, y_loc, z_loc, .true.)
              call localToCellBasedCoords(x_loc, y_loc, z_loc,&
                                        & xi, yi, zi, dx, dy, dz)
              dens = lg_arr(xi, yi, zi)
              dens_GJ = 2.0 * psr_omega * bz(xi, yi, zi) * B_norm / (CC * unit_ch)
              u_ = 0.0; v_ = 0.0; w_ = 0.0
            end if
            if ((nGJ_limiter .eq. 0) .or.&
              & ((dens_GJ .lt. 0) .and. (dens .gt. dens_GJ)) .or.&
              & ((dens_GJ .gt. 0) .and. (dens .lt. dens_GJ))) then
              call injectParticleGlobally(1, x_glob, y_glob, z_glob, u_, v_, w_, weight)
              call injectParticleGlobally(2, x_glob, y_glob, z_glob, u_, v_, w_, weight)
            end if
          end if
        end if
      end do
    else if (inj_method .eq. 2) then
      ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ! second injection option: fraction of polar GJ
      ppc = 0.5 * ppc0
      if ((nGJ_limiter .ne. 0) .or. (sigGJ_limiter .ne. 0)) then
        ! compute number density and write to `lg_arr`
        call computeDensity(1, reset=.true., ds=0, charge=.false.)
        call computeDensity(2, reset=.false., ds=0, charge=.false.)
      end if
      n_part = INT((4.0 * M_PI / 3.0) * ((psr_radius + shell_width + inj_dr)**3 - (psr_radius + inj_dr)**3) * ppc)
      do n = 1, n_part
        call randomPointInSphericalShell(psr_radius + inj_dr, psr_radius + inj_dr + shell_width, x_glob, y_glob, z_glob)
        rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
        nx = x_glob / rr
        ny = y_glob / rr
        nz = z_glob / rr
        x_glob = x_glob + xc_g
        y_glob = y_glob + yc_g
        z_glob = z_glob + zc_g

        call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
        ! if particle is within the current MPI meshblock
        if (dummy_flag) then
          call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
          dummy_flag = .true.
          if ((sigGJ_limiter .ne. 0) .and. dummy_flag) then
            call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
            density = lg_arr(xi, yi, zi)
            b_sqr = bx0**2 + by0**2 + bz0**2
            if (density .gt. 0) then
              sig = b_sqr * sigma * ppc0 / density
            else
              sig = sigma
            end if
            dummy_flag = (sig .gt. sigma_nGJ * sigGJ_limiter)
          end if

          if ((nGJ_limiter .ne. 0) .and. dummy_flag) then
            density = lg_arr(xi, yi, zi)
            dummy_flag = (density .lt. nGJ * nGJ_limiter)
          end if

          if ((jdotb_limiter .ne. 0) .and. dummy_flag) then
            call interpFromEdges(dx, dy, dz, xi, yi, zi, jx, jy, jz, jx0, jy0, jz0)
            call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
            j_dot_b = (jx0 * bx0 + jy0 * by0 + jz0 * bz0) / sqrt(bx0**2 + by0**2 + bz0**2)
            dummy_flag = ((abs(j_dot_b) * B_norm .gt. jdotb_limiter * nGJ * CC * unit_ch) .or. (step .lt. 0.1 * psr_period))
          end if

          if (dummy_flag) then
            ! kick along local b-field:
            call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
            b_sqr = sqrt(bx0**2 + by0**2 + bz0**2)
            sign = 1
            if (bx0 * nx + by0 * ny + bz0 * nz .lt. 0) then
              sign = -1
              bx0 = -bx0; by0 = -by0; bz0 = -bz0
            end if
            nx = bx0 / b_sqr
            ny = by0 / b_sqr
            nz = bz0 / b_sqr

            weight = inj_mult * nGJ / ppc

            #ifndef GCA
              u_ = nx * prtl_kick
              v_ = ny * prtl_kick
              w_ = nz * prtl_kick
              call createParticle(1, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
              call createParticle(2, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
            #else
              ! kick along ExB:
              call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
              call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
              b0_SQR = bx0**2 + by0**2 + bz0**2
              e0_SQR = ex0**2 + ey0**2 + ez0**2

              dummy_ = 1.0 / (b0_SQR + TINYFLD)

              vE_x = (bz0 * ey0 - by0 * ez0) * dummy_
              vE_y = (-bz0 * ex0 + bx0 * ez0) * dummy_
              vE_z = (by0 * ex0 - bx0 * ey0) * dummy_

              dummy_ = 1.0 / sqrt(abs(1.0 - vE_x**2 - vE_y**2 - vE_z**2) + TINYFLD)

              u_ = vE_x * dummy_
              v_ = vE_y * dummy_
              w_ = vE_z * dummy_

              dummy_ = sqrt(abs(prtl_kick**2 - dummy_**2))

              u_ = u_ + nx * dummy_
              v_ = v_ + ny * dummy_
              w_ = w_ + nz * dummy_

              call createParticleFromAttributes(1, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
                                                & xi_past=xi, yi_past=yi, zi_past=zi,&
                                                & dx_past=dx, dy_past=dy, dz_past=dz,&
                                                & u=u_, v=v_, w=w_,&
                                                & u_eff=u_, v_eff=v_, w_eff=w_, u_par=sign*dummy_, u_perp=0.0,&
                                                & ind=species(1)%cntr_sp, proc=mpi_rank + 2 * mpi_size, weight=weight)
              species(1)%cntr_sp = species(1)%cntr_sp + 1
              call createParticleFromAttributes(2, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
                                                & xi_past=xi, yi_past=yi, zi_past=zi,&
                                                & dx_past=dx, dy_past=dy, dz_past=dz,&
                                                & u=u_, v=v_, w=w_,&
                                                & u_eff=u_, v_eff=v_, w_eff=w_, u_par=sign*dummy_, u_perp=0.0,&
                                                & ind=species(2)%cntr_sp, proc=mpi_rank + 2 * mpi_size, weight=weight)
              species(2)%cntr_sp = species(2)%cntr_sp + 1
            #endif
          end if
        end if
      end do
    end if
    if (inj_method .gt. 0) then
      rr = 0.5 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz) - ds_abs / 2.0
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
                if ((r_g .lt. (psr_radius - rmin_dr)) .or. (r_g .gt. rr)) then
                  species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
                end if
              end do
            end do
          end do
        end do
      end do
    end if
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

  subroutine getEparAt(E_dot_B, B_sqr, x0, y0, z0, contained_flag)
    implicit none
    real, intent(out)     :: E_dot_B, B_sqr
    real, intent(in)      :: x0, y0, z0
    logical, intent(out)  :: contained_flag
    real                  :: x_loc, y_loc, z_loc, dx, dy, dz
    integer(kind=2)       :: xi, yi, zi
    real                  :: ex0, ey0, ez0, bx0, by0, bz0
    real                  :: nx, ny, nz, rr
    nx = x0 - xc_g
    ny = y0 - yc_g
    nz = z0 - zc_g
    rr = sqrt(nx**2 + ny**2 + nz**2)
    nx = nx / rr; ny = ny / rr; nz = nz / rr
    ! ... measure `e_dr` cells above the injection point
    call globalToLocalCoords(x0 + nx * e_dr,&
                           & y0 + ny * e_dr,&
                           & z0 + nz * e_dr,&
                           & x_loc, y_loc, z_loc, containedQ=contained_flag)
    if (contained_flag) then
      call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
      ! interpolate fields on particle position + dr
      call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
      call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
      B_sqr = bx0**2 + by0**2 + bz0**2
      E_dot_B = (ex0 * bx0 + ey0 * by0 + ez0 * bz0)
    end if
  end subroutine getEparAt

  subroutine getDeltaErAt(delta_Er, B_sqr, x0, y0, z0, contained_flag)
    implicit none
    real, intent(out)     :: delta_Er, B_sqr
    real, intent(in)      :: x0, y0, z0
    logical, intent(out)  :: contained_flag
    real                  :: x_loc, y_loc, z_loc, dx, dy, dz
    integer(kind=2)       :: xi, yi, zi
    real                  :: ex0, ey0, ez0, bx0, by0, bz0, vx, vy, vz
    real                  :: nx, ny, nz, rr, ex_cor, ey_cor, ez_cor, er_cor, er
    nx = x0 - xc_g
    ny = y0 - yc_g
    nz = z0 - zc_g
    rr = sqrt(nx**2 + ny**2 + nz**2)
    nx = nx / rr; ny = ny / rr; nz = nz / rr
    ! ... measure `e_dr` cells above the injection point
    call globalToLocalCoords(x0 + nx * e_dr,&
                           & y0 + ny * e_dr,&
                           & z0 + nz * e_dr,&
                           & x_loc, y_loc, z_loc, containedQ=contained_flag)
    if (contained_flag) then
      call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
      ! interpolate fields on particle position + dr
      call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
      call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
      vx = -psr_omega * ny * rr
      vy = psr_omega * nx * rr
      ex_cor = -vy * bz0 * CCINV
      ey_cor = vx * bz0 * CCINV
      ez_cor = -(vx * by0 - vy * bx0) * CCINV
      er_cor = ex_cor * nx + ey_cor * ny + ez_cor * nz
      er = ex0 * nx + ey0 * ny + ez0 * nz
      delta_Er = er - er_cor
      B_sqr = bx0**2 + by0**2 + bz0**2
    end if
  end subroutine getDeltaErAt

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

    rlimit = step * CC + psr_radius + shell_width + inj_dr
    if (rlimit .lt. 0.6 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz)) then
      ! check if the sphere with radius `rlimit` intersects the current meshblock ...
      ! ... this additional step should speed things up a bit
      x_ = max(REAL(this_meshblock%ptr%x0), min(xc_g, REAL(this_meshblock%ptr%x0 + this_meshblock%ptr%sx - 1)))
      y_ = max(REAL(this_meshblock%ptr%y0), min(yc_g, REAL(this_meshblock%ptr%y0 + this_meshblock%ptr%sy - 1)))
      z_ = max(REAL(this_meshblock%ptr%z0), min(zc_g, REAL(this_meshblock%ptr%z0 + this_meshblock%ptr%sz - 1)))
      rr = sqrt(REAL(x_ - xc_g)**2 + REAL(y_ - yc_g)**2 + REAL(z_ - zc_g)**2)
      if (rr .le. rlimit) then
        ! damp E-field inside a sphere
        do i = 0, this_meshblock%ptr%sx - 1
          i_glob = i + this_meshblock%ptr%x0
          do j = 0, this_meshblock%ptr%sy - 1
            j_glob = j + this_meshblock%ptr%y0
            do k = 0, this_meshblock%ptr%sz - 1
              k_glob = k + this_meshblock%ptr%z0
              x_ = REAL(i_glob);  y_ = REAL(j_glob);  z_ = REAL(k_glob)
              rr = sqrt(REAL(x_ - xc_g)**2 + REAL(y_ - yc_g)**2 + REAL(z_ - zc_g)**2)
              if (rr .gt. rlimit) then
                ex(i, j, k) = 0;  ey(i, j, k) = 0;  ez(i, j, k) = 0
              end if
            end do
          end do
        end do
      end if
    end if

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

    if ((psr_spinupT .eq. 0.0) .or. (step .gt. 4 * psr_spinupT)) then
      psr_omega = psr_omega0
    else
      psr_omega = psr_omega0 * tanh(REAL(step) / psr_spinupT)
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
