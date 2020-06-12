! save interpolation of fields at `t = n`
ex0_n = ex0;  ey0_n = ey0;  ez0_n = ez0
bx0_n = bx0;  by0_n = by0;  bz0_n = bz0

! compute position at `t = n-1/2`
! and save position at `t = n`
#if defined(oneD) || defined(twoD) || defined(threeD)
  x_ = (REAL(pt_xi(p) + pt_xi_past(p)) + (pt_dx(p) + pt_dx_past(p))) * 0.5
  xi_ = INT(x_, 2);   dx_ = x_ - REAL(xi_);
  x_n = REAL(pt_xi(p)) + pt_dx(p)
#endif

#if defined(twoD) || defined(threeD)
  y_ = (REAL(pt_yi(p) + pt_yi_past(p)) + (pt_dy(p) + pt_dy_past(p))) * 0.5
  yi_ = INT(y_, 2);   dy_ = y_ - REAL(yi_)
  y_n = REAL(pt_yi(p)) + pt_dy(p)
#endif

#if defined(threeD)
  z_ = (REAL(pt_zi(p) + pt_zi_past(p)) + (pt_dz(p) + pt_dz_past(p))) * 0.5
  zi_ = INT(z_, 2);   dz_ = z_ - REAL(zi_)
  z_n = REAL(pt_zi(p)) + pt_dz(p)
#endif

#ifdef oneD
  yi_ = 0_2; dy_ = 0.5
  zi_ = 0_2; dz_ = 0.5
#elif twoD
  zi_ = 0_2; dz_ = 0.5
#endif

! save velocities at `t = n-1/2`
u0 = pt_u(p); v0 = pt_v(p); w0 = pt_w(p)

! compute fields at `t = n-1/2`
#ifdef oneD
  lind = xi_
#elif twoD
  lind = xi_ + (NGHOST + yi_) * iy
#elif threeD
  lind = xi_ + (NGHOST + yi_) * iy + (NGHOST + zi_) * iz
#endif

dx = dx_; dy = dy_; dz = dz_

! these "functions" take
! ... coordinates and linear index: `dx`, `dy`, `dz` and `lind` ...
! ... and "return" `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0`
include "interp_efield.F"
include "interp_bfield.F"

e0sqr = ex0**2 + ey0**2 + ez0**2
b0sqr = bx0**2 + by0**2 + bz0**2

! compute `vE` and `gammaE` at time `t = n-1/2`
! this "function" takes
! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
! ... and returns the `vE` components and the `gammaE`
include "gca_ecrossb_drift.F"

b0sqr = sqrt(b0sqr)

! at time `t = n-1/2`
bpar_x = bx0 / b0sqr
bpar_y = by0 / b0sqr
bpar_z = bz0 / b0sqr

u_par = bpar_x * u0 + bpar_y * v0 + bpar_z * w0

Gamma_ = sqrt(1.0 + u0**2 + v0**2 + w0**2)

uperp_g_x = u0 - u_par * bpar_x - vE_x * Gamma_
uperp_g_y = v0 - u_par * bpar_y - vE_y * Gamma_
uperp_g_z = w0 - u_par * bpar_z - vE_z * Gamma_
uperp_g = sqrt(uperp_g_x**2 + uperp_g_y**2 + uperp_g_z**2)

! GCA/Boris condition
if ((.not. species(s)%gca_sp) .or.&
  & ((CC**2 * uperp_g / (b0sqr * B_norm)) .ge. gca_rhomin) .or.&
  & (sqrt(e0sqr) .ge. gca_eoverbmin * b0sqr)) then
  ! use Boris here
  pt_wei(p) = 1.0

  ex0 = ex0_n; ey0 = ey0_n; ez0 = ez0_n
  bx0 = bx0_n; by0 = by0_n; bz0 = bz0_n
  ! this "function" takes
  ! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
  ! ... and the velocities: `u0`, `v0`, `w0` ...
  ! ... and returns the updated velocities `u0`, `v0`, `w0`
  include "boris_push.F"
  pt_u(p) = u0; pt_v(p) = v0; pt_w(p) = w0
  over_e_temp = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
  ! save old position
  pt_xi_past(p) = pt_xi(p)
  pt_yi_past(p) = pt_yi(p)
  pt_zi_past(p) = pt_zi(p)

  pt_dx_past(p) = pt_dx(p)
  pt_dy_past(p) = pt_dy(p)
  pt_dz_past(p) = pt_dz(p)
  ! this "function" takes
  ! ... inverse energy: `over_e_temp` ...
  ! ... reads the velocities from: `pt_*(p)` ...
  ! ... and updates the particle position `pt_*(p)`
  include "boris_update.F"
else
  ! use GCA here
  pt_wei(p) = -1.0

  ! at time `t = n-1/2`
  bperp_x = uperp_g_x / uperp_g
  bperp_y = uperp_g_y / uperp_g
  bperp_z = uperp_g_z / uperp_g

  ! compute E_|| at `t = n`
  b0sqr_n = sqrt(bx0_n**2 + by0_n**2 + bz0_n**2)
  epar0_n = (ex0_n * bx0_n + ey0_n * by0_n + ez0_n * bz0_n) / b0sqr_n
  ! update u_||
  u_par = u_par + q_over_m * B_norm * epar0_n * CCINV

  ! `mu` at time `t = n-1/2` (in units of mc^2/2)
  mu_ = uperp_g**2 / (b0sqr * gammaE)
  ! solve non-linear equation on `x^n+1`

  ! this "function" takes
  ! ... `mu` at `t = n-1/2` ...
  ! ... coordinates at `t = n` ...
  ! ... u_|| at `t = n+1/2` ...
  ! ... fields at `t = n`: `ex0_n`, ... `bx0_n`
  ! ... and "returns" the new coordinates `x_n1`, `y_n1`, `z_n1`
  include "gca_nonlin.F"

  ! save old position
  pt_xi_past(p) = pt_xi(p)
  pt_yi_past(p) = pt_yi(p)
  pt_zi_past(p) = pt_zi(p)

  pt_dx_past(p) = pt_dx(p)
  pt_dy_past(p) = pt_dy(p)
  pt_dz_past(p) = pt_dz(p)

  ! convert coordinates at `t = n+1` to `xi` and `dx` ...
  ! ... compute position at `t = n+1/2`
  ! ... and update particle coordinates
  #if defined(oneD) || defined(twoD) || defined(threeD)
    xi_ = INT(x_n1, 2);   dx_ = x_n1 - REAL(xi_)
    pt_xi(p) = xi_;       pt_dx(p) = dx_
    x_ = (REAL(pt_xi(p) + pt_xi_past(p)) + (pt_dx(p) + pt_dx_past(p))) * 0.5
    xi_ = INT(x_, 2);     dx_ = x_ - REAL(xi_)
  #endif

  #if defined(twoD) || defined(threeD)
    yi_ = INT(y_n1, 2);   dy_ = y_n1 - REAL(yi_)
    pt_yi(p) = yi_;       pt_dy(p) = dy_
    y_ = (REAL(pt_yi(p) + pt_yi_past(p)) + (pt_dy(p) + pt_dy_past(p))) * 0.5
    yi_ = INT(y_, 2);     dy_ = y_ - REAL(yi_)
  #endif

  #if defined(threeD)
    zi_ = INT(z_n1, 2);   dz_ = z_n1 - REAL(zi_)
    pt_zi(p) = zi_;       pt_dz(p) = dz_
    z_ = (REAL(pt_zi(p) + pt_zi_past(p)) + (pt_dz(p) + pt_dz_past(p))) * 0.5
    zi_ = INT(z_, 2);     dz_ = z_ - REAL(zi_)
  #endif

  ! interpolate fields at `t = n+1/2`
  dx = dx_; dy = dy_; dz = dz_

  ! compute fields at `t = n-1/2`
  ! these "functions" take
  ! ... coordinates and linear index: `dx`, `dy`, `dz` and `lind` ...
  ! ... and "return" `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0`
  include "interp_efield.F"
  include "interp_bfield.F"

  e0sqr = ex0**2 + ey0**2 + ez0**2
  b0sqr = bx0**2 + by0**2 + bz0**2

  ! compute `vE` and `gammaE` at time `t = n-1/2`
  ! this "function" takes
  ! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
  ! ... their norms `e0sqr` and `b0sqr` ...
  ! ... and returns the `vE` components and the `gammaE`
  include "gca_ecrossb_drift.F"

  ! at time `t = n+1/2`
  uperp_g = sqrt(mu_ * b0sqr * gammaE)

  bpar_x = bx0 / sqrt(b0sqr)
  bpar_y = by0 / sqrt(b0sqr)
  bpar_z = bz0 / sqrt(b0sqr)

  Gamma_ = gammaE * sqrt(1.0 + u_par**2 + mu_ * b0sqr * gammaE)

  ! recover cartesian velocities
  pt_u(p) = bpar_x * u_par + vE_x * Gamma_ + uperp_g * bperp_x
  pt_v(p) = bpar_y * u_par + vE_y * Gamma_ + uperp_g * bperp_y
  pt_w(p) = bpar_z * u_par + vE_z * Gamma_ + uperp_g * bperp_z
end if
