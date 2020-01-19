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

ex0 = 0.5 * (&
        & (1.0 - pt_dy(p)) *&
          & ((1.0 - pt_dx(p)) * ex(lind-1,1,1) + ex(lind,1,1) + pt_dx(p) * ex(1 + lind,1,1)) +&
        & pt_dy(p) *&
          & ((1.0 - pt_dx(p)) * ex(lind-1+mx,1,1) + ex(lind+mx,1,1) + pt_dx(p)*ex(1 + lind + mx,1,1))&
      & )

ey0 = 0.5 * (&
        & (1.0 - pt_dy(p)) *&
          & ((1.0 - pt_dx(p)) * (ey(lind,1,1) + ey(lind-mx,1,1)) +&
            & pt_dx(p) * (ey(lind+1,1,1) + ey(lind+1-mx,1,1))) +&
        & pt_dy(p) * (&
          & (1.0 - pt_dx(p)) * (ey(lind,1,1) + ey(lind + mx,1,1)) +&
            & pt_dx(p) * (ey(lind+1,1,1) + ey(lind+1+mx,1,1)))&
      & )

ez0 =  (1.0 - pt_dx(p)) * (1.0 - pt_dy(p)) * ez(lind,1,1) +&
        & (pt_dx(p) - pt_dx(p) * pt_dy(p)) * ez(lind+1,1,1) +&
        & (1.0 - pt_dx(p)) * pt_dy(p) * ez(lind+mx,1,1) +&
        & pt_dx(p) * pt_dy(p) * ez(lind + 1 + mx,1,1)
