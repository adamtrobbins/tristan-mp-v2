! this "function" takes
! ... the field quantities: `bx0`, `by0`, `bz0`, `ex0`, `ey0`, `ez0` ...
! ... their norms `e0_SQR` and `b0_SQR` ...
! ... and returns the `vE` components and the `gammaE`

! dummy_ = 1.0 / (e0_SQR + b0_SQR + TINYFLD)
!
! wE_x = (bz0 * ey0 - by0 * ez0) * dummy_
! wE_y = (-bz0 * ex0 + bx0 * ez0) * dummy_
! wE_z = (by0 * ex0 - bx0 * ey0) * dummy_
! wE_SQR = wE_x**2 + wE_y**2 + wE_z**2
!
! dummy_ = sign(1.0, wE_SQR - 0.01)
!
! ! if `wE_SQR < 0.01` -- using taylor expansion
! dummy_ = 0.5 * (1.0 - dummy_) * (1.0 + wE_SQR + 2.0 * wE_SQR**2) + 0.25 * (1.0 + dummy_) * (1.0 - sqrt(max(1.0 - 4.0 * wE_SQR, 0.0))) / max(wE_SQR, 0.001)
!
! vE_x = wE_x * dummy_
! vE_y = wE_y * dummy_
! vE_z = wE_z * dummy_

dummy_ = 1.0 / (e0_SQR + b0_SQR + TINYFLD)

vE_x = (bz0 * ey0 - by0 * ez0) * dummy_
vE_y = (-bz0 * ex0 + bx0 * ez0) * dummy_
vE_z = (by0 * ex0 - bx0 * ey0) * dummy_

gammaE = 1.0 / sqrt(abs(1.0 - vE_x**2 - vE_y**2 - vE_z**2) + TINYFLD)
