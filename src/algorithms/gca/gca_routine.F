#include "../../defs.F90"

! save interpolation of fields at `t = n`
ex0_n = ex0;  ey0_n = ey0;  ez0_n = ez0
bx0_n = bx0;  by0_n = by0;  bz0_n = bz0
bb0_n = sqrt(bx0_n**2 + by0_n**2 + bz0_n**2)

! compute position at `t = n-1/2`
! and save position at `t = n`
#if defined(oneD) || defined(twoD) || defined(threeD)
  x_ = (REAL(pt_xi(p) + pt_xi_past(p)) + (pt_dx(p) + pt_dx_past(p))) * 0.5
  xi_ = FLOOR(x_);   dx_ = x_ - REAL(xi_);
  x_n = REAL(pt_xi(p)) + pt_dx(p)
  pt_xi_past(p) = pt_xi(p);     pt_dx_past(p) = pt_dx(p)
#endif

#if defined(twoD) || defined(threeD)
  y_ = (REAL(pt_yi(p) + pt_yi_past(p)) + (pt_dy(p) + pt_dy_past(p))) * 0.5
  yi_ = FLOOR(y_);   dy_ = y_ - REAL(yi_)
  y_n = REAL(pt_yi(p)) + pt_dy(p)
  pt_yi_past(p) = pt_yi(p);     pt_dy_past(p) = pt_dy(p)
#endif

#if defined(threeD)
  z_ = (REAL(pt_zi(p) + pt_zi_past(p)) + (pt_dz(p) + pt_dz_past(p))) * 0.5
  zi_ = FLOOR(z_);   dz_ = z_ - REAL(zi_)
  z_n = REAL(pt_zi(p)) + pt_dz(p)
  pt_zi_past(p) = pt_zi(p);     pt_dz_past(p) = pt_dz(p)
#endif

#ifdef oneD
  yi_ = 0_2; dy_ = 0.5
  zi_ = 0_2; dz_ = 0.5
#elif twoD
  zi_ = 0_2; dz_ = 0.5
#endif

! save velocities at `t = n-1/2`
u0 = pt_u(p); v0 = pt_v(p); w0 = pt_w(p)
Gamma_ = sqrt(1.0 + u0**2 + v0**2 + w0**2)

! compute fields at `t = n-1/2`
#ifdef oneD
  lind = xi_
#elif twoD
  lind = xi_ + (NGHOST + yi_) * iy
#elif threeD
  lind = xi_ + (NGHOST + yi_) * iy + (NGHOST + zi_) * iz
#endif
dx = dx_; dy = dy_; dz = dz_
include "interp_efield.F"
include "interp_bfield.F"

e0_SQR = ex0**2 + ey0**2 + ez0**2
b0_SQR = bx0**2 + by0**2 + bz0**2
#include "gca_ecrossb_drift.F"

bb0 = sqrt(b0_SQR)
dummy_ = 1.0 / sqrt(b0_SQR + TINYFLD)

! at time `t = n-1/2`
bpar_x = bx0 * dummy_
bpar_y = by0 * dummy_
bpar_z = bz0 * dummy_

if (pt_proc(p) .ge. mpi_size) then
  ! particle was doing GCA last timestep
  u_par = pt_u_par(p)
  uperp_g = 0.0 !pt_u_perp(p)

  ! to recover the direction later
  uperp_g_x = u0 - u_par * bpar_x - vE_x * Gamma_
  uperp_g_y = v0 - u_par * bpar_y - vE_y * Gamma_
  uperp_g_z = w0 - u_par * bpar_z - vE_z * Gamma_
  dummy_ = sqrt(uperp_g_x**2 + uperp_g_y**2 + uperp_g_z**2)
  uperp_g_x = uperp_g_x * uperp_g / dummy_
  uperp_g_y = uperp_g_y * uperp_g / dummy_
  uperp_g_z = uperp_g_z * uperp_g / dummy_
else
  ! particle was doing Boris last timestep
  u_par = bpar_x * u0 + bpar_y * v0 + bpar_z * w0

  uperp_g_x = u0 - u_par * bpar_x - vE_x * Gamma_
  uperp_g_y = v0 - u_par * bpar_y - vE_y * Gamma_
  uperp_g_z = w0 - u_par * bpar_z - vE_z * Gamma_
  uperp_g = sqrt(uperp_g_x**2 + uperp_g_y**2 + uperp_g_z**2)
end if

doBorisQ = ((.not. species(s)%gca_sp) .or.&                                                   ! if species cannot do GCA
          & (bb0 * bb0_n .eq. 0.0) .or.&                                                      ! if B = 0
          & ((CC**2 * Gamma_ / (abs(q_over_m) * bb0 * B_norm)) .ge. gca_rhomin) .or.&         ! if Larmor too large
          & (uperp_g / Gamma_ .gt. gca_vperpmax) .or. &                                       ! if uperp_g too large
          & (sqrt(e0_SQR) .ge. gca_eoverbmin * bb0))                                          ! of E > B * fraction

! FOR PARTICLE DIAGNOSTICS >>>>
if (gca_follow) then
  if (((pt_proc(p) .eq. gca_follow_proc) .or. (pt_proc(p) .eq. gca_follow_proc + 2 * mpi_size)) .and.&
    & (s .eq. gca_follow_sp) .and.&
    & (pt_ind(p) .eq. gca_follow_ind)) then
    dummy_ = sqrt(1.0 + u0**2 + v0**2 + w0**2)
    print *, 'particle...', (pt_proc(p) .gt. mpi_size)
    print *, 'gamma:', dummy_
    print *, 'E, B:', sqrt(e0_SQR), bb0, sqrt(e0_SQR) / bb0
    print *, 'radius:', sqrt(real((xi_ + dx_ + this_meshblock%ptr%x0 - 280.0)**2 +&
                                & (yi_ + dy_ + this_meshblock%ptr%y0 - 280.0)**2 +&
                                & (zi_ + dz_ + this_meshblock%ptr%z0 - 280.0)**2))

    print *, 'upar, uperp, uperp/gamma, gammaE * vE:',&
                  & u_par, uperp_g, uperp_g / Gamma_, gammaE * sqrt(vE_x**2 + vE_y**2 + vE_z**2)
    print *, 'mu:', uperp_g**2 / (bb0 * gammaE)
    print *, 'rhoL:', CC**2 * dummy_ / (abs(q_over_m) * bb0 * B_norm)
    print *, 'rhoL_real:', CC**2 * uperp_g / (abs(q_over_m) * bb0 * B_norm)
    print *, 'omegaC_real:', abs(q_over_m) * bb0 * B_norm / (CC * uperp_g)

    if (doBorisQ) then
      print *, 'will do Boris now'
    else
      print *, 'will do GCA now'
    end if

    print *, '-----------------------------------------------------'
  end if
end if

if (sqrt(1.0 + u0**2 + v0**2 + w0**2) .gt. 1e8) then
  print *, 'bad particle!!!!', pt_proc(p), s, pt_ind(p), sqrt(1.0 + u0**2 + v0**2 + w0**2)
  stop
end if
! <<<< FOR PARTICLE DIAGNOSTICS

! GCA/Boris condition . . . . . . . . . . . . . . . . . . . . . . . . . . . .
if (doBorisQ) then
  if (pt_proc(p) .ge. mpi_size) pt_proc(p) = pt_proc(p) - 2 * mpi_size
  ! use Boris here
  ex0 = ex0_n; ey0 = ey0_n; ez0 = ez0_n
  bx0 = bx0_n; by0 = by0_n; bz0 = bz0_n
  include "boris_push.F"
  pt_u(p) = u0;     pt_v(p) = v0;     pt_w(p) = w0
  over_e_temp = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
  include "boris_update.F"
  ! compute "effective" velocities for the deposit
  pt_u_eff(p) = u0 * CC * over_e_temp
  pt_v_eff(p) = v0 * CC * over_e_temp
  pt_w_eff(p) = w0 * CC * over_e_temp
else
  if (pt_proc(p) .lt. mpi_size) pt_proc(p) = pt_proc(p) + 2 * mpi_size
  ! use GCA here
  ! at time `t = n-1/2`
  dummy_ = 1.0 / (uperp_g + TINYFLD)
  bperp_x = uperp_g_x * dummy_
  bperp_y = uperp_g_y * dummy_
  bperp_z = uperp_g_z * dummy_

  ! compute E_|| at `t = n`
  epar0_n = (ex0_n * bx0_n + ey0_n * by0_n + ez0_n * bz0_n) / bb0_n
  ! update u_||
  u_par = u_par + q_over_m * B_norm * epar0_n * CCINV

  ! `mu` at time `t = n-1/2` (in units of mc^2/2)
  mu_ = uperp_g**2 / (bb0 * gammaE)

  ! . . . . Solve nonlinear equation on coordinates at `t = n+1` . . . . . . . . . . . . .
  ! compute `vE` and `gammaE` at time `t = n` - - >
  ex0 = ex0_n; ey0 = ey0_n; ez0 = ez0_n
  bx0 = bx0_n; by0 = by0_n; bz0 = bz0_n
  #include "gca_ecrossb_drift.F"

  vE_x_n = vE_x; vE_y_n = vE_y; vE_z_n = vE_z; gammaE_n = gammaE
  ! < - - compute `vE` and `gammaE` at time `t = n`

  ! `Gamma` at `t = n`
  Gamma_n = gammaE_n * sqrt(1.0 + u_par**2 + mu_ * bb0_n * gammaE_n)

  ! initial guess
  #if defined(oneD) || defined(twoD) || defined(threeD)
    x_n1 = x_n
  #endif

  #if defined(twoD) || defined(threeD)
    y_n1 = y_n
  #endif

  #if defined(threeD)
    z_n1 = z_n
  #endif

  do iter = 1, 4
    ! convert coordinates at `t = n+1` to `xi` and `dx`
    #if defined(oneD) || defined(twoD) || defined(threeD)
      xi_ = FLOOR(x_n1); dx_ = x_n1 - REAL(xi_)
    #endif

    #if defined(twoD) || defined(threeD)
      yi_ = FLOOR(y_n1); dy_ = y_n1 - REAL(yi_)
    #endif

    #if defined(threeD)
      zi_ = FLOOR(z_n1); dz_ = z_n1 - REAL(zi_)
    #endif

    ! compute fields at `t = n+1`
    #ifdef oneD
      lind = xi_
    #elif twoD
      lind = xi_ + (NGHOST + yi_) * iy
    #elif threeD
      lind = xi_ + (NGHOST + yi_) * iy + (NGHOST + zi_) * iz
    #endif
    dx = dx_; dy = dy_; dz = dz_
    include "interp_efield.F"
    include "interp_bfield.F"

    e0_SQR = ex0**2 + ey0**2 + ez0**2
    b0_SQR = bx0**2 + by0**2 + bz0**2
    #include "gca_ecrossb_drift.F"

    vE_x_n1 = vE_x;   vE_y_n1 = vE_y;   vE_z_n1 = vE_z;   gammaE_n1 = gammaE
    bx0_n1 = bx0;     by0_n1 = by0;     bz0_n1 = bz0
    bb0_n1 = sqrt(b0_SQR)

    ! `Gamma` at `t = n`
    Gamma_n1 = gammaE_n1 * sqrt(1.0 + u_par**2 + mu_ * bb0_n1 * gammaE_n1)

    ! update coordinates
    #if defined(oneD) || defined(twoD) || defined(threeD)
      x_n1 = x_n + 0.5 * ((bx0_n / bb0_n) / Gamma_n + (bx0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                 & 0.5 * (vE_x_n + vE_x_n1) * CC
    #endif

    #if defined(twoD) || defined(threeD)
      y_n1 = y_n + 0.5 * ((by0_n / bb0_n) / Gamma_n + (by0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                 & 0.5 * (vE_y_n + vE_y_n1) * CC
    #endif

    #if defined(threeD)
      z_n1 = z_n + 0.5 * ((bz0_n / bb0_n) / Gamma_n + (bz0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                 & 0.5 * (vE_z_n + vE_z_n1) * CC
    #endif
  end do
  ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  ! compute "effective" velocities for the deposit
  #ifdef oneD
    pt_u_eff(p) = x_n1 - x_n
    pt_v_eff(p) = 0.5 * ((by0_n / bb0_n) / Gamma_n + (by0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                & 0.5 * (vE_y_n + vE_y_n1) * CC
    pt_w_eff(p) = 0.5 * ((bz0_n / bb0_n) / Gamma_n + (bz0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                & 0.5 * (vE_z_n + vE_z_n1) * CC
  #elif twoD
    pt_u_eff(p) = x_n1 - x_n
    pt_v_eff(p) = y_n1 - y_n
    pt_w_eff(p) = 0.5 * ((bz0_n / bb0_n) / Gamma_n + (bz0_n1 / bb0_n1) / Gamma_n1) * u_par * CC +&
                & 0.5 * (vE_z_n + vE_z_n1) * CC
  #elif threeD
    pt_u_eff(p) = x_n1 - x_n
    pt_v_eff(p) = y_n1 - y_n
    pt_w_eff(p) = z_n1 - z_n
  #endif

  ! convert coordinates at `t = n+1` to `xi` and `dx` ...
  ! ... compute position at `t = n+1/2`
  ! ... and update particle coordinates
  #if defined(oneD) || defined(twoD) || defined(threeD)
    xi_ = FLOOR(x_n1);   dx_ = x_n1 - REAL(xi_)
    pt_xi(p) = xi_;       pt_dx(p) = dx_
    x_ = (REAL(pt_xi(p) + pt_xi_past(p)) + (pt_dx(p) + pt_dx_past(p))) * 0.5
    xi_ = FLOOR(x_);     dx_ = x_ - REAL(xi_)
  #endif

  #if defined(twoD) || defined(threeD)
    yi_ = FLOOR(y_n1);   dy_ = y_n1 - REAL(yi_)
    pt_yi(p) = yi_;       pt_dy(p) = dy_
    y_ = (REAL(pt_yi(p) + pt_yi_past(p)) + (pt_dy(p) + pt_dy_past(p))) * 0.5
    yi_ = FLOOR(y_);     dy_ = y_ - REAL(yi_)
  #endif

  #if defined(threeD)
    zi_ = FLOOR(z_n1);   dz_ = z_n1 - REAL(zi_)
    pt_zi(p) = zi_;       pt_dz(p) = dz_
    z_ = (REAL(pt_zi(p) + pt_zi_past(p)) + (pt_dz(p) + pt_dz_past(p))) * 0.5
    zi_ = FLOOR(z_);     dz_ = z_ - REAL(zi_)
  #endif

  ! interpolate fields at `t = n+1/2`
  #ifdef oneD
    lind = xi_
  #elif twoD
    lind = xi_ + (NGHOST + yi_) * iy
  #elif threeD
    lind = xi_ + (NGHOST + yi_) * iy + (NGHOST + zi_) * iz
  #endif
  dx = dx_; dy = dy_; dz = dz_
  include "interp_efield.F"
  include "interp_bfield.F"

  e0_SQR = ex0**2 + ey0**2 + ez0**2
  b0_SQR = bx0**2 + by0**2 + bz0**2
  #include "gca_ecrossb_drift.F"

  ! FOR PARTICLE DIAGNOSTICS >>>>
  if (gca_follow) then
    if (((pt_proc(p) .eq. gca_follow_proc) .or. (pt_proc(p) .eq. gca_follow_proc + 2 * mpi_size)) .and.&
      & (s .eq. gca_follow_sp) .and.&
      & (pt_ind(p) .eq. gca_follow_ind)) then
      print *, 'mu [prev]:', mu_
      print *, 'uperp_g [prev]:', uperp_g
      print *, 'uperp_g [new]:', sqrt(mu_ * bb0 * gammaE)
      print *, '-----------------------------------------------------'
    end if
  end if
  ! <<<< FOR PARTICLE DIAGNOSTICS

  ! at time `t = n+1/2`
  bb0 = sqrt(b0_SQR)
  uperp_g = sqrt(mu_ * bb0 * gammaE)

  dummy_ = 1.0 / sqrt(b0_SQR + TINYFLD)

  bpar_x = bx0 * dummy_
  bpar_y = by0 * dummy_
  bpar_z = bz0 * dummy_

  Gamma_ = gammaE * sqrt(1.0 + u_par**2 + mu_ * bb0 * gammaE)

  ! recover cartesian velocities
  pt_u(p) = bpar_x * u_par + vE_x * Gamma_ + uperp_g * bperp_x
  pt_v(p) = bpar_y * u_par + vE_y * Gamma_ + uperp_g * bperp_y
  pt_w(p) = bpar_z * u_par + vE_z * Gamma_ + uperp_g * bperp_z

  ! save par/perp velocities
  pt_u_par(p) = u_par
  pt_u_perp(p) = uperp_g
end if ! boris/gca criterion
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
