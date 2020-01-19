! !e_x
! c000 = 0.5 * (ex(pt_xi(p),pt_yi(p),pt_zi(p)) + ex(pt_xi(p)-1,pt_yi(p),pt_zi(p)))
! c100 = 0.5 * (ex(pt_xi(p),pt_yi(p),pt_zi(p)) + ex(pt_xi(p) + 1,pt_yi(p),pt_zi(p)))
! c010 = 0.5 * (ex(pt_xi(p),pt_yi(p) + 1,pt_zi(p)) + ex(pt_xi(p)-1,pt_yi(p) + 1,pt_zi(p)))
! c110 = 0.5 * (ex(pt_xi(p),pt_yi(p) + 1,pt_zi(p)) + ex(pt_xi(p) + 1,pt_yi(p) + 1,pt_zi(p)))
! c00 = c000 * (1-pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1-pt_dx(p)) + c110 * pt_dx(p)
! ex0 = c00 * (1-pt_dy(p)) + c10 * pt_dy(p)
!
! !e_y
! c000 = 0.5 * (ey(pt_xi(p),pt_yi(p),pt_zi(p)) + ey(pt_xi(p),pt_yi(p)-1,pt_zi(p)))
! c100 = 0.5 * (ey(pt_xi(p) + 1,pt_yi(p),pt_zi(p)) + ey(pt_xi(p) + 1,pt_yi(p)-1,pt_zi(p)))
! c010 = 0.5 * (ey(pt_xi(p),pt_yi(p),pt_zi(p)) + ey(pt_xi(p),pt_yi(p) + 1,pt_zi(p)))
! c110 = 0.5 * (ey(pt_xi(p) + 1,pt_yi(p),pt_zi(p)) + ey(pt_xi(p) + 1,pt_yi(p) + 1,pt_zi(p)))
! c00 = c000 * (1-pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1-pt_dx(p)) + c110 * pt_dx(p)
! ey0 = c00 * (1-pt_dy(p)) + c10 * pt_dy(p)
!
! !e_z
! c000 = ez(pt_xi(p),pt_yi(p),pt_zi(p))
! c100 = ez(pt_xi(p) + 1,pt_yi(p),pt_zi(p))
! c010 = ez(pt_xi(p),pt_yi(p) + 1,pt_zi(p))
! c110 = ez(pt_xi(p) + 1,pt_yi(p) + 1,pt_zi(p))
! c00 = c000 * (1-pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1-pt_dx(p)) + c110 * pt_dx(p)
! ez0 = c00 * (1-pt_dy(p)) + c10 * pt_dy(p)

! !e_x
! c000 = 0.5 * (ex(lind, 1, 1) + ex(lind - 1, 1, 1))
! c100 = 0.5 * (ex(lind, 1, 1) + ex(lind + 1, 1, 1))
! c010 = 0.5 * (ex(lind + mx, 1, 1) + ex(lind + mx - 1, 1, 1))
! c110 = 0.5 * (ex(lind + mx, 1, 1) + ex(lind + mx + 1, 1, 1))
! c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
! ex0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)
!
! !e_y
! c000 = 0.5 * (ey(lind, 1, 1) + ey(lind - mx, 1, 1))
! c100 = 0.5 * (ey(lind + 1, 1, 1) + ey(lind + 1 - mx, 1, 1))
! c010 = 0.5 * (ey(lind, 1, 1) + ey(lind + mx, 1, 1))
! c110 = 0.5 * (ey(lind + 1, 1, 1) + ey(lind + mx + 1, 1, 1))
! c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
! ey0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)
!
! !e_z
! c000 = ez(lind, 1, 1)
! c100 = ez(lind + 1, 1, 1)
! c010 = ez(lind + mx, 1, 1)
! c110 = ez(lind + mx + 1, 1, 1)
! c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
! ez0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

!e_x
c000 = 0.5 * (pt_ex(lind, 1, 1) + pt_ex(lind - 1, 1, 1))
c100 = 0.5 * (pt_ex(lind, 1, 1) + pt_ex(lind + 1, 1, 1))
c010 = 0.5 * (pt_ex(lind + mx, 1, 1) + pt_ex(lind + mx - 1, 1, 1))
c110 = 0.5 * (pt_ex(lind + mx, 1, 1) + pt_ex(lind + mx + 1, 1, 1))
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
ex0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

!e_y
c000 = 0.5 * (pt_ey(lind, 1, 1) + pt_ey(lind - mx, 1, 1))
c100 = 0.5 * (pt_ey(lind + 1, 1, 1) + pt_ey(lind + 1 - mx, 1, 1))
c010 = 0.5 * (pt_ey(lind, 1, 1) + pt_ey(lind + mx, 1, 1))
c110 = 0.5 * (pt_ey(lind + 1, 1, 1) + pt_ey(lind + mx + 1, 1, 1))
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
ey0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

!e_z
c000 = pt_ez(lind, 1, 1)
c100 = pt_ez(lind + 1, 1, 1)
c010 = pt_ez(lind + mx, 1, 1)
c110 = pt_ez(lind + mx + 1, 1, 1)
c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
ez0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)
