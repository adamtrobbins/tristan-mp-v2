!b_x
c000 = 0.5 * (bx(lind, 1, 1) + bx(lind - iy, 1, 1))
c100 = 0.5 * (bx(lind + 1, 1, 1) + bx(lind + 1 - iy, 1, 1))
c010 = 0.5 * (bx(lind, 1, 1) + bx(lind + iy, 1, 1))
c110 = 0.5 * (bx(lind + 1, 1, 1) + bx(lind + iy + 1, 1, 1))
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
bx0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

!b_y
c000 = 0.5 * (by(lind - 1, 1, 1) + by(lind, 1, 1))
c100 = 0.5 * (by(lind, 1, 1) + by(lind + 1, 1, 1))
c010 = 0.5 * (by(lind + iy - 1, 1, 1) + by(lind + iy, 1, 1))
c110 = 0.5 * (by(lind + iy, 1, 1) + by(lind + iy + 1, 1, 1))
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
by0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

!b_z
c000 = 0.25 * (bz(lind - 1 - iy, 1, 1) + bz(lind - 1, 1, 1) + &
              &bz(lind - iy, 1, 1) + bz(lind, 1, 1))
c100 = 0.25 * (bz(lind - iy, 1, 1) + bz(lind, 1, 1) + &
              &bz(lind + 1 - iy, 1, 1) + bz(lind + 1, 1, 1))
c010 = 0.25 * (bz(lind - 1, 1, 1) + bz(lind + iy - 1, 1, 1) + &
              &bz(lind, 1, 1) + bz(lind + iy, 1, 1))
c110 = 0.25 * (bz(lind, 1, 1) + bz(lind + iy, 1, 1) + &
              &bz(lind + 1, 1, 1) + bz(lind + iy + 1, 1, 1))
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
bz0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)
