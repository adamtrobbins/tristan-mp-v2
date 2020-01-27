! linear indices used isntead of 3d

#ifndef threeD
  ! c000=0.5*(bx(pt_xi(p),pt_yi(p),pt_zi(p))+bx(pt_xi(p),pt_yi(p)-1,pt_zi(p)))
  ! c100=0.5*(bx(pt_xi(p)+1,pt_yi(p),pt_zi(p))+bx(pt_xi(p)+1,pt_yi(p)-1,pt_zi(p)))
  ! c010=0.5*(bx(pt_xi(p),pt_yi(p),pt_zi(p))+bx(pt_xi(p),pt_yi(p)+1,pt_zi(p)))
  ! c110=0.5*(bx(pt_xi(p)+1,pt_yi(p),pt_zi(p))+bx(pt_xi(p)+1,pt_yi(p)+1,pt_zi(p)))
  ! c00=c000*(1.0-pt_dx(p))+c100*pt_dx(p)
  ! c10=c010*(1.0-pt_dx(p))+c110*pt_dx(p)
  ! bx0=c00*(1.0-pt_dy(p))+c10*pt_dy(p)
  !
  ! c000=0.5*(by(pt_xi(p)-1,pt_yi(p),pt_zi(p))+by(pt_xi(p),pt_yi(p),pt_zi(p)))
  ! c100=0.5*(by(pt_xi(p),pt_yi(p),pt_zi(p))+by(pt_xi(p)+1,pt_yi(p),pt_zi(p)))
  ! c010=0.5*(by(pt_xi(p)-1,pt_yi(p)+1,pt_zi(p))+by(pt_xi(p),pt_yi(p)+1,pt_zi(p)))
  ! c110=0.5*(by(pt_xi(p),pt_yi(p)+1,pt_zi(p))+by(pt_xi(p)+1,pt_yi(p)+1,pt_zi(p)))
  ! c00=c000*(1.0-pt_dx(p))+c100*pt_dx(p)
  ! c10=c010*(1.0-pt_dx(p))+c110*pt_dx(p)
  ! by0=c00*(1.0-pt_dy(p))+c10*pt_dy(p)
  !
  ! c000=0.25*(bz(pt_xi(p)-1,pt_yi(p)-1,pt_zi(p))+bz(pt_xi(p)-1,pt_yi(p),pt_zi(p))+&
  ! &bz(pt_xi(p),pt_yi(p)-1,pt_zi(p))+bz(pt_xi(p),pt_yi(p),pt_zi(p)))
  ! c100=0.25*(bz(pt_xi(p),pt_yi(p)-1,pt_zi(p))+bz(pt_xi(p),pt_yi(p),pt_zi(p))+&
  ! &bz(pt_xi(p)+1,pt_yi(p)-1,pt_zi(p))+bz(pt_xi(p)+1,pt_yi(p),pt_zi(p)))
  ! c010=0.25*(bz(pt_xi(p)-1,pt_yi(p),pt_zi(p))+bz(pt_xi(p)-1,pt_yi(p)+1,pt_zi(p))+&
  ! &bz(pt_xi(p),pt_yi(p),pt_zi(p))+bz(pt_xi(p),pt_yi(p)+1,pt_zi(p)))
  ! c110=0.25*(bz(pt_xi(p),pt_yi(p),pt_zi(p))+bz(pt_xi(p),pt_yi(p)+1,pt_zi(p))+&
  ! &bz(pt_xi(p)+1,pt_yi(p),pt_zi(p))+bz(pt_xi(p)+1,pt_yi(p)+1,pt_zi(p)))
  ! c00=c000*(1.0-pt_dx(p))+c100*pt_dx(p)
  ! c10=c010*(1.0-pt_dx(p))+c110*pt_dx(p)
  ! bz0=c00*(1.0-pt_dy(p))+c10*pt_dy(p)
  
  !b_x
  c000 = 0.5 * (bx(lind, -NGHOST, 0) + bx(lind - iy, -NGHOST, 0))
  c100 = 0.5 * (bx(lind + 1, -NGHOST, 0) + bx(lind + 1 - iy, -NGHOST, 0))
  c010 = 0.5 * (bx(lind, -NGHOST, 0) + bx(lind + iy, -NGHOST, 0))
  c110 = 0.5 * (bx(lind + 1, -NGHOST, 0) + bx(lind + iy + 1, -NGHOST, 0))
  c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
  c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
  bx0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

  !b_y
  c000 = 0.5 * (by(lind - 1, -NGHOST, 0) + by(lind, -NGHOST, 0))
  c100 = 0.5 * (by(lind, -NGHOST, 0) + by(lind + 1, -NGHOST, 0))
  c010 = 0.5 * (by(lind + iy - 1, -NGHOST, 0) + by(lind + iy, -NGHOST, 0))
  c110 = 0.5 * (by(lind + iy, -NGHOST, 0) + by(lind + iy + 1, -NGHOST, 0))
  c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
  c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
  by0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)

  !b_z
  c000 = 0.25 * (bz(lind - 1 - iy, -NGHOST, 0) + bz(lind - 1, -NGHOST, 0) + &
                &bz(lind - iy, -NGHOST, 0) + bz(lind, -NGHOST, 0))
  c100 = 0.25 * (bz(lind - iy, -NGHOST, 0) + bz(lind, -NGHOST, 0) + &
                &bz(lind + 1 - iy, -NGHOST, 0) + bz(lind + 1, -NGHOST, 0))
  c010 = 0.25 * (bz(lind - 1, -NGHOST, 0) + bz(lind + iy - 1, -NGHOST, 0) + &
                &bz(lind, -NGHOST, 0) + bz(lind + iy, -NGHOST, 0))
  c110 = 0.25 * (bz(lind, -NGHOST, 0) + bz(lind + iy, -NGHOST, 0) + &
                &bz(lind + 1, -NGHOST, 0) + bz(lind + iy + 1, -NGHOST, 0))
  c00 = c000 * (1 - pt_dx(p)) + c100 * pt_dx(p)
  c10 = c010 * (1 - pt_dx(p)) + c110 * pt_dx(p)
  bz0 = c00 * (1 - pt_dy(p)) + c10 * pt_dy(p)
#else
  !b_x
  c000 = 0.25 * (bx(lind, -NGHOST, -NGHOST) + bx(lind - iy, -NGHOST, -NGHOST) +&
               & bx(lind - iz, -NGHOST, -NGHOST) + bx(lind - iy - iz, -NGHOST, -NGHOST))
  c100 = 0.25 * (bx(lind + 1, -NGHOST, -NGHOST) + bx(lind + 1 - iy, -NGHOST, -NGHOST) +&
               & bx(lind + 1 - iz, -NGHOST, -NGHOST) + bx(lind + 1 - iy - iz, -NGHOST, -NGHOST))
  c001 = 0.25 * (bx(lind, -NGHOST, -NGHOST) + bx(lind + iz, -NGHOST, -NGHOST) +&
               & bx(lind - iy, -NGHOST, -NGHOST) + bx(lind - iy + iz, -NGHOST, -NGHOST))
  c101 = 0.25 * (bx(lind + 1, -NGHOST, -NGHOST) + bx(lind + 1 + iz, -NGHOST, -NGHOST) +&
               & bx(lind + 1 - iy, -NGHOST, -NGHOST) + bx(lind + 1 - iy + iz, -NGHOST, -NGHOST))
  c010 = 0.25 * (bx(lind, -NGHOST, -NGHOST) + bx(lind + iy, -NGHOST, -NGHOST) +&
               & bx(lind - iz, -NGHOST, -NGHOST) + bx(lind + iy - iz, -NGHOST, -NGHOST))
  c110 = 0.25 * (bx(lind + 1, -NGHOST, -NGHOST) + bx(lind + 1 - iz, -NGHOST, -NGHOST) +&
               & bx(lind + 1 + iy - iz, -NGHOST, -NGHOST) + bx(lind + 1 + iy, -NGHOST, -NGHOST))
  c011 = 0.25 * (bx(lind, -NGHOST, -NGHOST) + bx(lind + iy, -NGHOST, -NGHOST) +&
               & bx(lind + iy + iz, -NGHOST, -NGHOST) + bx(lind + iz, -NGHOST, -NGHOST))
  c111 = 0.25 * (bx(lind + 1, -NGHOST, -NGHOST) + bx(lind + 1 + iy, -NGHOST, -NGHOST) +&
               & bx(lind + 1 + iy + iz, -NGHOST, -NGHOST) + bx(lind + 1 + iz, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - pt_dx(p)) + c100 * pt_dx(p)
  c01 = c001 * (1.0 - pt_dx(p)) + c101 * pt_dx(p)
  c10 = c010 * (1.0 - pt_dx(p)) + c110 * pt_dx(p)
  c11 = c011 * (1.0 - pt_dx(p)) + c111 * pt_dx(p)
  c0 = c00 * (1.0 - pt_dy(p)) + c10 * pt_dy(p)
  c1 = c01 * (1.0 - pt_dy(p)) + c11 * pt_dy(p)
  bx0 = c0 * (1.0 - pt_dz(p)) + c1 * pt_dz(p)

  !b_y
  c000 = 0.25 * (by(lind - 1 - iz, -NGHOST, -NGHOST) + by(lind - 1, -NGHOST, -NGHOST) +&
               & by(lind - iz, -NGHOST, -NGHOST) + by(lind, -NGHOST, -NGHOST))
  c100 = 0.25 * (by(lind - iz, -NGHOST, -NGHOST) + by(lind, -NGHOST, -NGHOST) +&
               & by(lind + 1 - iz, -NGHOST, -NGHOST) + by(lind + 1, -NGHOST, -NGHOST))
  c001 = 0.25 * (by(lind - 1, -NGHOST, -NGHOST) + by(lind - 1 + iz, -NGHOST, -NGHOST) +&
               & by(lind, -NGHOST, -NGHOST) + by(lind + iz, -NGHOST, -NGHOST))
  c101 = 0.25 * (by(lind, -NGHOST, -NGHOST) + by(lind + iz, -NGHOST, -NGHOST) +&
               & by(lind + 1, -NGHOST, -NGHOST) + by(lind + 1 + iz, -NGHOST, -NGHOST))
  c010 = 0.25 * (by(lind - 1 + iy - iz, -NGHOST, -NGHOST) + by(lind - 1 + iy, -NGHOST, -NGHOST) +&
               & by(lind + iy - iz, -NGHOST, -NGHOST) + by(lind + iy, -NGHOST, -NGHOST))
  c110 = 0.25 * (by(lind + iy - iz, -NGHOST, -NGHOST) + by(lind + iy, -NGHOST, -NGHOST) +&
               & by(lind + 1 + iy - iz, -NGHOST, -NGHOST) + by(lind + 1 + iy, -NGHOST, -NGHOST))
  c011 = 0.25 * (by(lind - 1 + iy, -NGHOST, -NGHOST) + by(lind - 1 + iy + iz, -NGHOST, -NGHOST) +&
               & by(lind + iy, -NGHOST, -NGHOST) + by(lind + iy + iz, -NGHOST, -NGHOST))
  c111 = 0.25 * (by(lind + iy, -NGHOST, -NGHOST) + by(lind + iy + iz, -NGHOST, -NGHOST) +&
               & by(lind + 1 + iy, -NGHOST, -NGHOST) + by(lind + 1 + iy + iz, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - pt_dx(p)) + c100 * pt_dx(p)
  c01 = c001 * (1.0 - pt_dx(p)) + c101 * pt_dx(p)
  c10 = c010 * (1.0 - pt_dx(p)) + c110 * pt_dx(p)
  c11 = c011 * (1.0 - pt_dx(p)) + c111 * pt_dx(p)
  c0 = c00 * (1.0 - pt_dy(p)) + c10 * pt_dy(p)
  c1 = c01 * (1.0 - pt_dy(p)) + c11 * pt_dy(p)
  by0 = c0 * (1.0 - pt_dz(p)) + c1 * pt_dz(p)

  !b_z
  c000 = 0.25 * (bz(lind - 1 - iy, -NGHOST, -NGHOST) + bz(lind - 1, -NGHOST, -NGHOST) +&
             & bz(lind - iy, -NGHOST, -NGHOST) + bz(lind, -NGHOST, -NGHOST))
  c100 = 0.25 * (bz(lind - iy, -NGHOST, -NGHOST) + bz(lind, -NGHOST, -NGHOST) +&
             & bz(lind + 1 - iy, -NGHOST, -NGHOST) + bz(lind + 1, -NGHOST, -NGHOST))
  c001 = 0.25 * (bz(lind - 1 - iy + iz, -NGHOST, -NGHOST) + bz(lind - 1 + iz, -NGHOST, -NGHOST) +&
             & bz(lind - iy + iz, -NGHOST, -NGHOST) + bz(lind + iz, -NGHOST, -NGHOST))
  c101 = 0.25 * (bz(lind - iy + iz, -NGHOST, -NGHOST) + bz(lind + iz, -NGHOST, -NGHOST) +&
             & bz(lind + 1 - iy + iz, -NGHOST, -NGHOST) + bz(lind + 1 + iz, -NGHOST, -NGHOST))
  c010 = 0.25 * (bz(lind - 1, -NGHOST, -NGHOST) + bz(lind - 1 + iy, -NGHOST, -NGHOST) +&
             & bz(lind, -NGHOST, -NGHOST) + bz(lind + iy, -NGHOST, -NGHOST))
  c110 = 0.25 * (bz(lind, -NGHOST, -NGHOST) + bz(lind + iy, -NGHOST, -NGHOST) +&
             & bz(lind + 1, -NGHOST, -NGHOST) + bz(lind + 1 + iy, -NGHOST, -NGHOST))
  c011 = 0.25 * (bz(lind - 1 + iz, -NGHOST, -NGHOST) + bz(lind - 1 + iy + iz, -NGHOST, -NGHOST) +&
             & bz(lind + iz, -NGHOST, -NGHOST) + bz(lind + iy + iz, -NGHOST, -NGHOST))
  c111 = 0.25 * (bz(lind + iz, -NGHOST, -NGHOST) + bz(lind + iy + iz, -NGHOST, -NGHOST) +&
             & bz(lind + 1 + iz, -NGHOST, -NGHOST) + bz(lind + 1 + iy + iz, -NGHOST, -NGHOST))
  c00 = c000 * (1.0 - pt_dx(p)) + c100 * pt_dx(p)
  c01 = c001 * (1.0 - pt_dx(p)) + c101 * pt_dx(p)
  c10 = c010 * (1.0 - pt_dx(p)) + c110 * pt_dx(p)
  c11 = c011 * (1.0 - pt_dx(p)) + c111 * pt_dx(p)
  c0 = c00 * (1.0 - pt_dy(p)) + c10 * pt_dy(p)
  c1 = c01 * (1.0 - pt_dy(p)) + c11 * pt_dy(p)
  bz0 = c0 * (1.0 - pt_dz(p)) + c1 * pt_dz(p)
#endif
