subroutine dumpParticlesHere(xi, yi, zi, dx, dy, dz, weight)
  implicit none
  integer(kind=2), intent(in) :: xi, yi, zi
  real, intent(in)            :: dx, dy, dz, weight
  call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
  call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
end subroutine dumpParticlesHere

subroutine userParticleBoundaryConditions(step)
  implicit none
  integer, optional, intent(in) :: step
  integer                       :: s, ti, tj, tk, p
  real                          :: x_g, y_g, z_g, r_g
  integer                       :: n_part, n, sign
  real                          :: x_loc, y_loc, z_loc, dx, dy, dz
  integer(kind=2)               :: xi, yi, zi, xi_eb, yi_eb, zi_eb
  real                          :: x_glob, y_glob, z_glob, weight, ppc, dens, sig, rmax, rmin
  real                          :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
  real                          :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma, rmax_sph, smax_car
  real                          :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
  real                          :: bz_at_surface, mu_x, mu_y, mu_z, mu_dot_n, phase, dummy_, cos_theta
  logical                       :: dummy_flag
  #ifdef GCA
    real                          :: vE_x, vE_y, vE_z
    real                          :: e0_SQR, b0_SQR
  #endif
  real        :: dr_i, rmin_inj, rmax_inj
  integer     :: ri, factor = 10, fi
  dr_i = 2.0

  #ifdef RADIATION
    ! disable cooling until certain point
    if (step .lt. cooling_on * psr_period) then
      species(1)%cool_sp = .false.
      species(2)%cool_sp = .false.
    else
      species(1)%cool_sp = .true.
      species(2)%cool_sp = .true.
    end if
  #endif

  if (sigGJ_limiter .ne. 0) then
    ! compute number density and write to `lg_arr`
    call computeDensity(1, reset=.true., ds=0, charge=.false.)
    call computeDensity(2, reset=.false., ds=0, charge=.false.)
  end if

  nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))
  ppc = 0.5 * ppc0
  weight = inj_mult * nGJ / ppc
  do ri = 0, INT(shell_width / dr_i)
    ! iterate over spherical shells
    rmin_inj = psr_radius + inj_dr + REAL(ri) * dr_i
    rmax_inj = psr_radius + inj_dr + REAL(ri + 1) * dr_i
    n_part = INT(2.0 * (4.0 * M_PI / 3.0) * (rmax_inj**3 - rmin_inj**3) * ppc / factor)
    do n = 1, n_part
      ! for each spherical shell iterate over the given number of particles to inject
      call randomPointInSphericalShell(rmin_inj, rmax_inj, x_glob, y_glob, z_glob)
      rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
      nx = x_glob / rr
      ny = y_glob / rr
      nz = z_glob / rr
      cos_theta = (nz)
      x_glob = x_glob + xc_g
      y_glob = y_glob + yc_g
      z_glob = z_glob + zc_g
      call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
      ! if particle is within the current MPI meshblock
      !if (dummy_flag .and. (rr .gt. 0.6 * psr_rlc * (1.0 - cos_theta**2))) then
      if (dummy_flag) then
        ! open zone injection
        call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
        dummy_flag = .true.
        ! limiter on sigma
        if (sigGJ_limiter .ne. 0) then
          density = lg_arr(xi, yi, zi)
          b_sqr = bx(xi, yi, zi)**2 + by(xi, yi, zi)**2 + bz(xi, yi, zi)**2
          if (density .gt. 0) then
            sig = b_sqr * sigma * ppc0 / density
          else
            sig = sigma
          end if
          dummy_flag = (sig .gt. sigGJ_limiter * (psr_radius / rr)**2)
        end if
        if (dummy_flag) then
          !weight = inj_mult * nGJ * (psr_radius / rr)**3 / ppc
          do fi = 1, factor
            dx = random(dseed)
            dy = random(dseed)
            dz = random(dseed)
            call dumpParticlesHere(xi, yi, zi, dx, dy, dz, weight)
          end do
        end if
      !else if (dummy_flag) then
        !! closed zone injection
        !call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
        !call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
        !call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
        !dummy_ = abs(ex0 * bx0 + ey0 * by0 + ez0 * bz0) / (bx0**2 + by0**2 + bz0**2)
        !if (dummy_ .gt. edotb_thr_closed) then
          !weight = dummy_ * nGJ * (psr_radius / rr)**3 / ppc
          !do fi = 1, factor
            !dx = random(dseed)
            !dy = random(dseed)
            !dz = random(dseed)
            !call dumpParticlesHere(xi, yi, zi, dx, dy, dz, weight)
          !end do
        !end if ! E.B limiter
      end if ! closed zone
    end do ! npart
  end do ! ri

  ! max radius for spherical absorption and absorbing layer size in cartesian absorption
  rmax_sph = 0.5 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz) - ds_abs + 4.0
  smax_car = ds_abs / 4.0
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
            if (&
              & (r_g .lt. (psr_radius - rmin_dr)) .or.&
              & ((boundary_x .eq. 2) .and. (r_g .gt. rmax_sph)) .or.&
              & ((boundary_x .eq. 0) .and. ((x_g .lt. smax_car) .or.&
                                          & (x_g .gt. global_mesh%sx - smax_car) .or.&
                                          & (y_g .lt. smax_car) .or.&
                                          & (y_g .gt. global_mesh%sy - smax_car) .or.&
                                          & (z_g .lt. smax_car) .or.&
                                          & (z_g .gt. global_mesh%sz - smax_car))&
                & )&
              & ) then
              species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
            end if
          end do
        end do
      end do
    end do
  end do
end subroutine userParticleBoundaryConditions

  !subroutine userParticleBoundaryConditions(step)
    !implicit none
    !integer, optional, intent(in) :: step
    !integer                       :: s, ti, tj, tk, p
    !real                          :: x_g, y_g, z_g, r_g
    !integer                       :: n_part, n, sign
    !real                          :: x_loc, y_loc, z_loc, dx, dy, dz
    !integer(kind=2)               :: xi, yi, zi, xi_eb, yi_eb, zi_eb
    !real                          :: x_glob, y_glob, z_glob, weight, ppc, dens, sig, rmax, rmin
    !real                          :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
    !real                          :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma, rmax_sph, smax_car
    !real                          :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
    !real                          :: bz_at_surface, mu_x, mu_y, mu_z, mu_dot_n, phase, dummy_
    !logical                       :: dummy_flag
    !#ifdef GCA
      !real                          :: vE_x, vE_y, vE_z
      !real                          :: e0_SQR, b0_SQR
    !#endif

    !phase = step * psr_omega
    !mu_x = sin(psr_angle) * cos(phase)
    !mu_y = sin(psr_angle) * sin(phase)
    !mu_z = cos(psr_angle)

    !#ifdef RADIATION
      !! disable cooling until certain point
      !if (step .lt. cooling_on * psr_period) then
        !species(1)%cool_sp = .false.
        !species(2)%cool_sp = .false.
      !else
        !species(1)%cool_sp = .true.
        !species(2)%cool_sp = .true.
      !end if
    !#endif

    !nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))

    !ppc = 0.5 * ppc0
    !if (sigGJ_limiter .ne. 0) then
      !! compute number density and write to `lg_arr`
      !call computeDensity(1, reset=.true., ds=0, charge=.false.)
      !call computeDensity(2, reset=.false., ds=0, charge=.false.)
    !end if
    !n_part = INT(2.0 * (4.0 * M_PI / 3.0) * ((psr_radius + shell_width + inj_dr)**3 - (psr_radius + inj_dr)**3) * ppc)
    !do n = 1, n_part
      !call randomPointInSphericalShell(psr_radius + inj_dr, psr_radius + inj_dr + shell_width, x_glob, y_glob, z_glob)
      !rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
      !nx = x_glob / rr
      !ny = y_glob / rr
      !nz = z_glob / rr

      !mu_dot_n = mu_x * nx + mu_y * ny + mu_z * nz

      !x_glob = x_glob + xc_g
      !y_glob = y_glob + yc_g
      !z_glob = z_glob + zc_g
      !call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
      !! if particle is within the current MPI meshblock
      !if (dummy_flag) then
        !call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
        !! limiter on sigma
        !dummy_flag = .true.
        !if (sigGJ_limiter .ne. 0) then
          !density = lg_arr(xi, yi, zi)
          !b_sqr = bx(xi, yi, zi)**2 + by(xi, yi, zi)**2 + bz(xi, yi, zi)**2
          !if (density .gt. 0) then
            !sig = b_sqr * sigma * ppc0 / density
          !else
            !sig = sigma
          !end if
          !dummy_flag = (sig .gt. sigGJ_limiter * (psr_radius / rr)**3)
        !end if ! sigma limiter enabled
        !if (dummy_flag) then
          !! initialize at rest
          !call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
          !call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
          !b_sqr = sqrt(bx0**2 + by0**2 + bz0**2)
          !sign = 1
          !if (bx0 * nx + by0 * ny + bz0 * nz .lt. 0) then
            !sign = -1
            !bx0 = -bx0; by0 = -by0; bz0 = -bz0
          !end if
          !nx = bx0 / b_sqr
          !ny = by0 / b_sqr
          !nz = bz0 / b_sqr
          !dummy_ = abs(ex0 * bx0 + ey0 * by0 + ez0 * bz0) / (bx0**2 + by0**2 + bz0**2)
          !if (dummy_ .gt. edotb_thr_closed) then
            !weight = inj_mult * dummy_ * nGJ / ppc
            !if (prtl_kick .eq. 0) then
              !call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
              !call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
            !else
              !#ifndef GCA
                !u_ = nx * prtl_kick
                !v_ = ny * prtl_kick
                !w_ = nz * prtl_kick
                !call createParticle(1, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
                !call createParticle(2, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
              !#else
                !! kick along ExB:
                !ex0 = ex(xi, yi, zi); ey0 = ey(xi, yi, zi); ez0 = ez(xi, yi, zi)
                !bx0 = bx(xi, yi, zi); by0 = by(xi, yi, zi); bz0 = bz(xi, yi, zi)
                !b0_SQR = bx0**2 + by0**2 + bz0**2
                !e0_SQR = ex0**2 + ey0**2 + ez0**2

                !dummy_ = 1.0 / (b0_SQR + TINYFLD)

                !vE_x = (bz0 * ey0 - by0 * ez0) * dummy_
                !vE_y = (-bz0 * ex0 + bx0 * ez0) * dummy_
                !vE_z = (by0 * ex0 - bx0 * ey0) * dummy_

                !dummy_ = 1.0 / sqrt(abs(1.0 - vE_x**2 - vE_y**2 - vE_z**2) + TINYFLD)

                !u_ = vE_x * dummy_
                !v_ = vE_y * dummy_
                !w_ = vE_z * dummy_

                !dummy_ = sqrt(prtl_kick**2 - 1.0)

                !u_ = u_ + nx * dummy_
                !v_ = v_ + ny * dummy_
                !w_ = w_ + nz * dummy_

                !call createParticleFromAttributes(1, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
                                                  !& xi_past=xi, yi_past=yi, zi_past=zi,&
                                                  !& dx_past=dx, dy_past=dy, dz_past=dz,&
                                                  !& u=u_, v=v_, w=w_,&
                                                  !& u_eff=u_, v_eff=v_, w_eff=w_, u_par=sign*dummy_, u_perp=0.0,&
                                                  !& ind=species(1)%cntr_sp, proc=mpi_rank + 2 * mpi_size, weight=weight)
                !species(1)%cntr_sp = species(1)%cntr_sp + 1
                !call createParticleFromAttributes(2, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
                                                  !& xi_past=xi, yi_past=yi, zi_past=zi,&
                                                  !& dx_past=dx, dy_past=dy, dz_past=dz,&
                                                  !& u=u_, v=v_, w=w_,&
                                                  !& u_eff=u_, v_eff=v_, w_eff=w_, u_par=sign*dummy_, u_perp=0.0,&
                                                  !& ind=species(2)%cntr_sp, proc=mpi_rank + 2 * mpi_size, weight=weight)
                !species(2)%cntr_sp = species(2)%cntr_sp + 1
              !#endif
            !end if
          !end if ! E.B limiter
        !end if ! SIGMA limiter
      !end if ! current MPI block
    !end do
  
    !! max radius for spherical absorption and absorbing layer size in cartesian absorption
    !rmax_sph = 0.5 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz) - ds_abs + 4.0
    !smax_car = ds_abs / 4.0
    !! remove particles falling into the star
    !do s = 1, nspec
      !do ti = 1, species(s)%tile_nx
        !do tj = 1, species(s)%tile_ny
          !do tk = 1, species(s)%tile_nz
            !do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              !x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dx(p)
              !y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dy(p)
              !z_g = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dz(p)
              !r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
              !if (&
                !& (r_g .lt. (psr_radius - rmin_dr)) .or.&
                !& ((boundary_x .eq. 2) .and. (r_g .gt. rmax_sph)) .or.&
                !& ((boundary_x .eq. 0) .and. ((x_g .lt. smax_car) .or.&
                                            !& (x_g .gt. global_mesh%sx - smax_car) .or.&
                                            !& (y_g .lt. smax_car) .or.&
                                            !& (y_g .gt. global_mesh%sy - smax_car) .or.&
                                            !& (z_g .lt. smax_car) .or.&
                                            !& (z_g .gt. global_mesh%sz - smax_car))&
                  !& )&
                !& ) then
                !species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
              !end if
            !end do
          !end do
        !end do
      !end do
    !end do
  !end subroutine userParticleBoundaryConditions

  !subroutine userParticleBoundaryConditions(step)
    !implicit none
    !integer, optional, intent(in) :: step
    !integer                       :: s, ti, tj, tk, p
    !real                          :: x_g, y_g, z_g, r_g
    !integer                       :: n_part, n, sign
    !real                          :: x_loc, y_loc, z_loc, dx, dy, dz
    !integer(kind=2)               :: xi, yi, zi, xi_eb, yi_eb, zi_eb
    !real                          :: x_glob, y_glob, z_glob, weight, ppc, dens, sig, rmax, rmin
    !real                          :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
    !real                          :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma, rmax_sph, smax_car
    !real                          :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
    !real                          :: bz_at_surface, mu_x, mu_y, mu_z, mu_dot_n, phase, dummy_
    !logical                       :: dummy_flag
    !#ifdef GCA
      !real                          :: vE_x, vE_y, vE_z
      !real                          :: e0_SQR, b0_SQR
    !#endif

    !phase = step * psr_omega
    !mu_x = sin(psr_angle) * cos(phase)
    !mu_y = sin(psr_angle) * sin(phase)
    !mu_z = cos(psr_angle)

    !#ifdef RADIATION
      !! disable cooling until certain point
      !if (step .lt. cooling_on * psr_period) then
        !species(1)%cool_sp = .false.
        !species(2)%cool_sp = .false.
      !else
        !species(1)%cool_sp = .true.
        !species(2)%cool_sp = .true.
      !end if
    !#endif

    !nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))

    !! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    !ppc = 0.5 * ppc0
    !if (sigGJ_limiter .ne. 0) then
      !! compute number density and write to `lg_arr`
      !call computeDensity(1, reset=.true., ds=0, charge=.false.)
      !call computeDensity(2, reset=.false., ds=0, charge=.false.)
    !end if
    !n_part = INT(2.0 * (4.0 * M_PI / 3.0) * ((psr_radius + shell_width + inj_dr)**3 - (psr_radius + inj_dr)**3) * ppc)
    !do n = 1, n_part
      !call randomPointInSphericalShell(psr_radius + inj_dr, psr_radius + inj_dr + shell_width, x_glob, y_glob, z_glob)
      !rr = sqrt(x_glob**2 + y_glob**2 + z_glob**2)
      !nx = x_glob / rr
      !ny = y_glob / rr
      !nz = z_glob / rr

      !mu_dot_n = mu_x * nx + mu_y * ny + mu_z * nz

      !bz_at_surface = (3.0 * nz * mu_dot_n - mu_z)

      !dummy_flag = (((random(dseed) * 2.0) .lt. abs(bz_at_surface)) .and. (bz_at_surface .gt. 0.0))
      !! injection in the open zone
      !if (dummy_flag) then
        !x_glob = x_glob + xc_g
        !y_glob = y_glob + yc_g
        !z_glob = z_glob + zc_g
        !call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
        !! if particle is within the current MPI meshblock
        !if (dummy_flag) then
          !call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
          !dummy_flag = .true.

          !! limiter on sigma
          !if (sigGJ_limiter .ne. 0) then
            !density = lg_arr(xi, yi, zi)
            !b_sqr = bx(xi, yi, zi)**2 + by(xi, yi, zi)**2 + bz(xi, yi, zi)**2
            !if (density .gt. 0) then
              !sig = b_sqr * sigma * ppc0 / density
            !else
              !sig = sigma
            !end if
            !dummy_flag = (sig .gt. sigGJ_limiter * (psr_radius / rr)**3)
          !end if

          !if (dummy_flag) then
            !! kick along local b-field:
            !bx0 = bx(xi, yi, zi); by0 = by(xi, yi, zi); bz0 = bz(xi, yi, zi)
            !b_sqr = sqrt(bx0**2 + by0**2 + bz0**2)
            !sign = 1
            !if (bx0 * nx + by0 * ny + bz0 * nz .lt. 0) then
              !sign = -1
              !bx0 = -bx0; by0 = -by0; bz0 = -bz0
            !end if
            !nx = bx0 / b_sqr
            !ny = by0 / b_sqr
            !nz = bz0 / b_sqr

            !weight = inj_mult * nGJ / ppc

            !#ifndef GCA
              !u_ = nx * prtl_kick
              !v_ = ny * prtl_kick
              !w_ = nz * prtl_kick
              !call createParticle(1, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
              !call createParticle(2, xi, yi, zi, dx, dy, dz, u_, v_, w_, weight=weight)
            !#else
              !! kick along ExB:
              !ex0 = ex(xi, yi, zi); ey0 = ey(xi, yi, zi); ez0 = ez(xi, yi, zi)
              !bx0 = bx(xi, yi, zi); by0 = by(xi, yi, zi); bz0 = bz(xi, yi, zi)
              !b0_SQR = bx0**2 + by0**2 + bz0**2
              !e0_SQR = ex0**2 + ey0**2 + ez0**2

              !dummy_ = 1.0 / (b0_SQR + TINYFLD)

              !vE_x = (bz0 * ey0 - by0 * ez0) * dummy_
              !vE_y = (-bz0 * ex0 + bx0 * ez0) * dummy_
              !vE_z = (by0 * ex0 - bx0 * ey0) * dummy_

              !dummy_ = 1.0 / sqrt(abs(1.0 - vE_x**2 - vE_y**2 - vE_z**2) + TINYFLD)

              !u_ = vE_x * dummy_
              !v_ = vE_y * dummy_
              !w_ = vE_z * dummy_

              !dummy_ = sqrt(prtl_kick**2 - 1.0)

              !u_ = u_ + nx * dummy_
              !v_ = v_ + ny * dummy_
              !w_ = w_ + nz * dummy_

              !call createParticleFromAttributes(1, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
                                                !& xi_past=xi, yi_past=yi, zi_past=zi,&
                                                !& dx_past=dx, dy_past=dy, dz_past=dz,&
                                                !& u=u_, v=v_, w=w_,&
                                                !& u_eff=u_, v_eff=v_, w_eff=w_, u_par=sign*dummy_, u_perp=0.0,&
                                                !& ind=species(1)%cntr_sp, proc=mpi_rank + 2 * mpi_size, weight=weight)
              !species(1)%cntr_sp = species(1)%cntr_sp + 1
              !call createParticleFromAttributes(2, xi=xi, yi=yi, zi=zi, dx=dx, dy=dy, dz=dz,&
                                                !& xi_past=xi, yi_past=yi, zi_past=zi,&
                                                !& dx_past=dx, dy_past=dy, dz_past=dz,&
                                                !& u=u_, v=v_, w=w_,&
                                                !& u_eff=u_, v_eff=v_, w_eff=w_, u_par=sign*dummy_, u_perp=0.0,&
                                                !& ind=species(2)%cntr_sp, proc=mpi_rank + 2 * mpi_size, weight=weight)
              !species(2)%cntr_sp = species(2)%cntr_sp + 1
            !#endif
          !end if
        !end if
      !end if ! open zone injection

      !! injection in the closed zone
      !if (bz_at_surface .le. 0.0) then
        !x_glob = x_glob + xc_g
        !y_glob = y_glob + yc_g
        !z_glob = z_glob + zc_g
        !call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
        !! if particle is within the current MPI meshblock
        !if (dummy_flag) then
          !call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
          !! initialize at rest
          !call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
          !call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
          !dummy_ = abs(ex0 * bx0 + ey0 * by0 + ez0 * bz0) / (bx0**2 + by0**2 + bz0**2)
          !if (dummy_ .gt. edotb_thr_closed) then
            !weight = dummy_ * nGJ / ppc
            !call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
            !call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
          !end if ! E.B limiter
        !end if ! current MPI block
      !end if ! closed zone
    !end do
  
    !! max radius for spherical absorption and absorbing layer size in cartesian absorption
    !rmax_sph = 0.5 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz) - ds_abs + 4.0
    !smax_car = ds_abs / 4.0
    !! remove particles falling into the star
    !do s = 1, nspec
      !do ti = 1, species(s)%tile_nx
        !do tj = 1, species(s)%tile_ny
          !do tk = 1, species(s)%tile_nz
            !do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              !x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dx(p)
              !y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dy(p)
              !z_g = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dz(p)
              !r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
              !if (&
                !& (r_g .lt. (psr_radius - rmin_dr)) .or.&
                !& ((boundary_x .eq. 2) .and. (r_g .gt. rmax_sph)) .or.&
                !& ((boundary_x .eq. 0) .and. ((x_g .lt. smax_car) .or.&
                                            !& (x_g .gt. global_mesh%sx - smax_car) .or.&
                                            !& (y_g .lt. smax_car) .or.&
                                            !& (y_g .gt. global_mesh%sy - smax_car) .or.&
                                            !& (z_g .lt. smax_car) .or.&
                                            !& (z_g .gt. global_mesh%sz - smax_car))&
                  !& )&
                !& ) then
                !species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
              !end if
            !end do
          !end do
        !end do
      !end do
    !end do
  !end subroutine userParticleBoundaryConditions
  
  !subroutine userParticleBoundaryConditions(step)
    !implicit none
    !integer, optional, intent(in) :: step
    !integer                       :: s, ti, tj, tk, p
    !real                          :: x_g, y_g, z_g, r_g
    !integer                       :: n_part, n, sign
    !real                          :: x_loc, y_loc, z_loc, dx, dy, dz
    !integer(kind=2)               :: xi, yi, zi, xi_eb, yi_eb, zi_eb
    !real                          :: x_glob, y_glob, z_glob, weight, ppc, dens, sig, rmax, rmin
    !real                          :: e_dot_b, b_sqr, delta_er, bx0, by0, bz0, ex0, ey0, ez0
    !real                          :: u_, v_, w_, nx, ny, nz, rr, vx, vy, vz, gamma, rmax_sph, smax_car
    !real                          :: dens_GJ, e_b_scale, j_dot_b, density, jx0, jy0, jz0
    !real                          :: bz_at_surface, mu_x, mu_y, mu_z, mu_dot_n, phase
    !logical                       :: dummy_flag
    !real                          :: dummy_

    !phase = step * psr_omega
    !mu_x = sin(psr_angle) * cos(phase)
    !mu_y = sin(psr_angle) * sin(phase)
    !mu_z = cos(psr_angle)
    !nGJ = 2 * psr_omega * B_norm * psr_bstar / (CC * abs(unit_ch))

    !! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    !ppc = 0.5 * ppc0
    !n_part = INT(2.0 * (4.0 * M_PI / 3.0) * ((psr_radius + shell_width + inj_dr)**3 - (psr_radius + inj_dr)**3) * ppc)
    !do n = 1, n_part
      !call randomPointInSphericalShell(psr_radius + inj_dr, psr_radius + inj_dr + shell_width, x_glob, y_glob, z_glob)
      !x_glob = x_glob + xc_g
      !y_glob = y_glob + yc_g
      !z_glob = z_glob + zc_g
      !call globalToLocalCoords(x_glob, y_glob, z_glob, x_loc, y_loc, z_loc, containedQ=dummy_flag)
      !! if particle is within the current MPI meshblock
      !if (dummy_flag) then
        !call localToCellBasedCoords(x_loc, y_loc, z_loc, xi, yi, zi, dx, dy, dz)
        !! initialize at rest
        !call interpFromEdges(dx, dy, dz, xi, yi, zi, ex, ey, ez, ex0, ey0, ez0)
        !call interpFromFaces(dx, dy, dz, xi, yi, zi, bx, by, bz, bx0, by0, bz0)
        !dummy_ = abs(ex0 * bx0 + ey0 * by0 + ez0 * bz0) / (bx0**2 + by0**2 + bz0**2)
        !if (dummy_ .gt. edotb_thr_closed) then
          !weight = dummy_ * nGJ / ppc
          !call createParticle(1, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
          !call createParticle(2, xi, yi, zi, dx, dy, dz, 0.0, 0.0, 0.0, weight=weight)
        !end if ! E.B limiter
      !end if ! current MPI block
    !end do
  
    !! max radius for spherical absorption and absorbing layer size in cartesian absorption
    !rmax_sph = 0.5 * MIN(global_mesh%sx, global_mesh%sy, global_mesh%sz) - ds_abs + 4.0
    !smax_car = ds_abs / 4.0
    !! remove particles falling into the star
    !do s = 1, nspec
      !do ti = 1, species(s)%tile_nx
        !do tj = 1, species(s)%tile_ny
          !do tk = 1, species(s)%tile_nz
            !do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              !x_g = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dx(p)
              !y_g = REAL(species(s)%prtl_tile(ti, tj, tk)%yi(p) + this_meshblock%ptr%y0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dy(p)
              !z_g = REAL(species(s)%prtl_tile(ti, tj, tk)%zi(p) + this_meshblock%ptr%z0)&
                  !& + species(s)%prtl_tile(ti, tj, tk)%dz(p)
              !r_g = sqrt((x_g - xc_g)**2 + (y_g - yc_g)**2 + (z_g - zc_g)**2)
              !if (&
                !& (r_g .lt. (psr_radius - rmin_dr)) .or.&
                !& ((boundary_x .eq. 2) .and. (r_g .gt. rmax_sph)) .or.&
                !& ((boundary_x .eq. 0) .and. ((x_g .lt. smax_car) .or.&
                                            !& (x_g .gt. global_mesh%sx - smax_car) .or.&
                                            !& (y_g .lt. smax_car) .or.&
                                            !& (y_g .gt. global_mesh%sy - smax_car) .or.&
                                            !& (z_g .lt. smax_car) .or.&
                                            !& (z_g .gt. global_mesh%sz - smax_car))&
                  !& )&
                !& ) then
                !species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
              !end if
            !end do
          !end do
        !end do
      !end do
    !end do
  !end subroutine userParticleBoundaryConditions

