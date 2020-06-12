! this "function" takes
! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
! ... their norms `e0sqr` and `b0sqr` ...
! ... and returns the `vE` components and the `gammaE`

wE_x = (bz0 * ey0 - by0 * ez0) / (e0sqr + b0sqr)
wE_y = (-bz0 * ex0 + bx0 * ez0) / (e0sqr + b0sqr)
wE_z = (by0 * ex0 - bx0 * ey0) / (e0sqr + b0sqr)
wEsqr = wE_x**2 + wE_y**2 + wE_z**2

dummy_ = sign(1.0, wEsqr - 0.01)

! if `wEsqr < 0.01` -- using taylor expansion
dummy_ = 0.5 * (1.0 - dummy_) * (1.0 + wEsqr + 2.0 * wEsqr**2) + 0.25 * (1.0 + dummy_) * (1.0 - sqrt(1.0 - 4.0 * wEsqr)) / max(wEsqr, 0.001)

vE_x = wE_x * dummy_
vE_y = wE_y * dummy_
vE_z = wE_z * dummy_

gammaE = 1.0 / sqrt(1.0 - vE_x**2 - vE_y**2 - vE_z**2)
