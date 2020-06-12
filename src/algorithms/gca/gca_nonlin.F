! this "function" takes
! ... `mu` at `t = n-1/2` ...
! ... coordinates at `t = n` ...
! ... u_|| at `t = n+1/2` ...
! ... fields at `t = n`: `ex0_n`, ... `bx0_n`
! ... and "returns" the new coordinates `x_n1`, `y_n1`, `z_n1`

! compute `vE` and `gammaE` at time `t = n` - - >
ex0 = ex0_n; ey0 = ey0_n; ez0 = ez0_n
bx0 = bx0_n; by0 = by0_n; bz0 = bz0_n
! this "function" takes
! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
! ... and returns the `vE` components and the `gammaE`
include "gca_ecrossb_drift.F"
vE_x_n = vE_x; vE_y_n = vE_y; vE_z_n = vE_z; gammaE_n = gammaE
! < - - compute `vE` and `gammaE` at time `t = n`

! `Gamma` at `t = n`
Gamma_n = gammaE_n * sqrt(1.0 + u_par**2 + mu_ * b0sqr_n * gammaE_n)

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
    xi_ = INT(x_n1, 2); dx_ = x_n1 - REAL(xi_)
  #endif

  #if defined(twoD) || defined(threeD)
    yi_ = INT(y_n1, 2); dy_ = y_n1 - REAL(yi_)
  #endif

  #if defined(threeD)
    zi_ = INT(z_n1, 2); dz_ = z_n1 - REAL(zi_)
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
  ! these "functions" take
  ! ... coordinates and linear index: `dx`, `dy`, `dz` and `lind` ...
  ! ... and "return" `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0`
  include "interp_efield.F"
  include "interp_bfield.F"

  e0sqr = ex0**2 + ey0**2 + ez0**2
  b0sqr = bx0**2 + by0**2 + bz0**2

  ! `ex0`, `ey0` and `ez0` now contain `E` at `t = n+1`
  ! compute `vE` and `gammaE` at `t = n+1`
  ! this "function" takes
  ! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
  ! ... their norms `e0sqr` and `b0sqr` ...
  ! ... and returns the `vE` components and the `gammaE`
  include "gca_ecrossb_drift.F"

  vE_x_n1 = vE_x;   vE_y_n1 = vE_y;   vE_z_n1 = vE_z;   gammaE_n1 = gammaE
  bx0_n1 = bx0;     by0_n1 = by0;     bz0_n1 = bz0
  b0sqr_n1 = sqrt(b0sqr)

  ! `Gamma` at `t = n`
  Gamma_n1 = gammaE_n1 * sqrt(1.0 + u_par**2 + mu_ * b0sqr_n1 * gammaE_n1)

  ! update coordinates
  #if defined(oneD) || defined(twoD) || defined(threeD)
    x_n1 = x_n + 0.5 * ((bx0_n / b0sqr_n) / Gamma_n + (bx0_n1 / b0sqr_n1) / Gamma_n1) * u_par +&
               & 0.5 * (vE_x_n + vE_x_n1)
  #endif

  #if defined(twoD) || defined(threeD)
    y_n1 = y_n + 0.5 * ((by0_n / b0sqr_n) / Gamma_n + (by0_n1 / b0sqr_n1) / Gamma_n1) * u_par +&
               & 0.5 * (vE_y_n + vE_y_n1)
  #endif

  #if defined(threeD)
    z_n1 = z_n + 0.5 * ((bz0_n / b0sqr_n) / Gamma_n + (bz0_n1 / b0sqr_n1) / Gamma_n1) * u_par +&
               & 0.5 * (vE_z_n + vE_z_n1)
  #endif
end do
