!e_x
c000 = 0.5 * (ex(lind, 1, 1) + ex(lind - 1, 1, 1))
c100 = 0.5 * (ex(lind, 1, 1) + ex(lind + 1, 1, 1))
c010 = 0.5 * (ex(lind + iy, 1, 1) + ex(lind + iy - 1, 1, 1))
c110 = 0.5 * (ex(lind + iy, 1, 1) + ex(lind + iy + 1, 1, 1))
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
ex0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

!e_y
c000 = 0.5 * (ey(lind, 1, 1) + ey(lind - iy, 1, 1))
c100 = 0.5 * (ey(lind + 1, 1, 1) + ey(lind + 1 - iy, 1, 1))
c010 = 0.5 * (ey(lind, 1, 1) + ey(lind + iy, 1, 1))
c110 = 0.5 * (ey(lind + 1, 1, 1) + ey(lind + iy + 1, 1, 1))
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
ey0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

!e_z
c000 = ez(lind, 1, 1)
c100 = ez(lind + 1, 1, 1)
c010 = ez(lind + iy, 1, 1)
c110 = ez(lind + iy + 1, 1, 1)
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
ez0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)
