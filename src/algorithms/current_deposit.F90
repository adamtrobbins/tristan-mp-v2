#include "../defs.F90"

module m_currentdeposit
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_fields
  use m_particles
  implicit none
contains
  subroutine depositCurrents()
    implicit none
    ! integer :: s, p
    ! integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    ! real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:),&
    !                                        & pt_u(:), pt_v(:), pt_w(:)
    ! real                                  :: xr, yr, zr, x1, y1, z1, x2, y2, z2
    ! real                                  :: gamma_inv
    ! integer(kind=2)                       :: i1, i2, j1, j2, k1, k2
    ! do s = 1, nspec
    !   if (spp_(s)%ch_sp .eq. 0) cycle
    !   pt_xi => sp_(s)%xi; pt_yi => sp_(s)%yi; pt_zi => sp_(s)%zi
    !   pt_dx => sp_(s)%dx; pt_dy => sp_(s)%dy; pt_dz => sp_(s)%dz
    !   pt_u => sp_(s)%u; pt_v => sp_(s)%v; pt_w => sp_(s)%w
    !   ! FIX1: vectorize & assume aligned the deposit
    !   do p = 1, spp_(s)%npart_sp
    !     ! push the particle back
    !     gamma_inv = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
    !     x2 = pt_xi(p) + pt_dx(p); y2 = pt_yi(p) + pt_dy(p); z2 = pt_zi(p) + pt_dz(p)
    !     x1 = x2 - pt_u(p) * CC * gamma_inv; y1 = y2 - pt_u(p) * CC * gamma_inv; z1 = z2 - pt_u(p) * CC * gamma_inv
		!
    !     i1 = INT(x1, 2); i2 = pt_xi(p)
    !     j1 = INT(y1, 2); j2 = pt_yi(p)
    !     k1 = INT(z1, 2); k2 = pt_zi(p)
		!
    !     xr = min(REAL(min(i1, i2) + 1), max(REAL(max(i1, i2)), 0.5 * (x1 + x2)))
    !   	yr = min(REAL(min(j1, j2) + 1), max(REAL(max(j1, j2)), 0.5 * (y1 + y2)))
    !   	zr = min(REAL(min(k1, k2) + 1), max(REAL(max(k1, k2)), 0.5 * (z1 + z2)))
    !   end do
    !   pt_xi => null(); pt_yi => null(); pt_zi => null()
    !   pt_dx => null(); pt_dy => null(); pt_dz => null()
    !   pt_u => null(); pt_v => null(); pt_w => null()
    ! end do
  end subroutine depositCurrents
end module m_currentdeposit
