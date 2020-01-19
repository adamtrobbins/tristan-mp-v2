! !b_x
! c000 = 0.5 * (bx(lind, 1, 1) + bx(lind - mx, 1, 1))
! c100 = 0.5 * (bx(lind + 1, 1, 1) + bx(lind + 1 - mx, 1, 1))
! c010 = 0.5 * (bx(lind, 1, 1) + bx(lind + mx, 1, 1))
! c110 = 0.5 * (bx(lind + 1, 1, 1) + bx(lind + mx + 1, 1, 1))
! c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
! bx0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)
!
! !b_y
! c000 = 0.5 * (by(lind - 1, 1, 1) + by(lind, 1, 1))
! c100 = 0.5 * (by(lind, 1, 1) + by(lind + 1, 1, 1))
! c010 = 0.5 * (by(lind + mx - 1, 1, 1) + by(lind + mx, 1, 1))
! c110 = 0.5 * (by(lind + mx, 1, 1) + by(lind + mx + 1, 1, 1))
! c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
! by0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)
!
! !b_z
! c000 = 0.25 * (bz(lind - 1 - mx, 1, 1) + bz(lind - 1, 1, 1) + &
!               &bz(lind - mx, 1, 1) + bz(lind, 1, 1))
! c100 = 0.25 * (bz(lind - mx, 1, 1) + bz(lind, 1, 1) + &
!               &bz(lind + 1 - mx, 1, 1) + bz(lind + 1, 1, 1))
! c010 = 0.25 * (bz(lind - 1, 1, 1) + bz(lind + mx - 1, 1, 1) + &
!               &bz(lind, 1, 1) + bz(lind + mx, 1, 1))
! c110 = 0.25 * (bz(lind, 1, 1) + bz(lind + mx, 1, 1) + &
!               &bz(lind + 1, 1, 1) + bz(lind + mx + 1, 1, 1))
! c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
! c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
! bz0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

bx0 = 0.5 * ((1.0 - pt_dy(p)) * ((1.0 - pt_dx(p)) * (bx(lind,1,1) + bx(lind-mx,1,1)) +&
        & pt_dx(p) * (bx(lind+1,1,1) + bx(lind+1-mx,1,1))) +&
        & pt_dy(p) * ((1.0 - pt_dx(p)) * (bx(lind,1,1) + bx(lind + mx,1,1)) +&
        & pt_dx(p) * (bx(lind+1,1,1) + bx(lind+1+mx,1,1))))

by0 = 0.5 * ((1.0 - pt_dy(p)) * ((1.0 - pt_dx(p)) * by(,1,1) - by(lind,1,1) - pt_dx(p)*by(lind+1,1,1)) +&
        & pt_dy(p) * (((1.0 - pt_dx(p)) * by(lind-1+mx,1,1)) + by(lind+mx,1,1) +&
        & pt_dx(p) * by(lind+1+mx,1,1)))

bz0 = 0.25 * ((1.0 - pt_dy(p)) * (((1.0 - pt_dx(p)) * (bz(lind-1,1,1) + bz(lind,1,1) +&
        & bz(lind-1-mx,1,1) + bz(lind-mx,1,1))) +&
        & pt_dx(p) * (bz(lind,1,1) + bz(lind+1,1,1) + bz(lind-mx,1,1) +&
          & bz(lind+1-mx,1,1))) +&
        & pt_dy(p) * (((1.0 - pt_dx(p)) * bz(lind-1,1,1)) + bz(lind,1,1) + pt_dx(p) * bz(lind+1,1,1) +&
          & bz(lind-1+mx,1,1) - pt_dx(p) * bz(lind-1+mx,1,1) + bz(lind+mx,1,1) +&
          & pt_dx(p) * bz(lind+1+mx,1,1)))
