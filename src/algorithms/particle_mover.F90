#include "../defs.F90"

module m_mover
  use m_globalnamespace
  use m_communications
  use m_domain
  use m_particles
  implicit none
contains
  subroutine moveParticles()
    implicit none
    integer                            :: s, p
    real(mprec)                        :: g_temp, e_temp
    real(mprec), pointer, contiguous   :: pt_x(:), pt_y(:), pt_z(:), pt_u(:), pt_v(:), pt_w(:)

    do s = 1, nspec
      pt_x => sp_(s)%x; pt_y => sp_(s)%y; pt_z => sp_(s)%z
      pt_u => sp_(s)%u; pt_v => sp_(s)%v; pt_w => sp_(s)%w
      if (spp_(s)%m_sp .eq. 0) then
        ! routine for massless particles
        !$omp simd
        !dir$ vector aligned
        do p = 1, spp_(s)%npart_sp
          e_temp = sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
          pt_x(p) = pt_x(p) + CC * pt_u(p) / e_temp
          pt_y(p) = pt_y(p) + CC * pt_v(p) / e_temp
          pt_z(p) = pt_z(p) + CC * pt_w(p) / e_temp
        end do
      else
        ! routine for massive particles
        !$omp simd
        !dir$ vector aligned
        do p = 1, spp_(s)%npart_sp
          g_temp = sqrt(1. + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)
          pt_x(p) = pt_x(p) + CC * pt_u(p) / g_temp
          pt_y(p) = pt_y(p) + CC * pt_v(p) / g_temp
          pt_z(p) = pt_z(p) + CC * pt_w(p) / g_temp
        end do
      end if
      pt_x => null(); pt_y => null(); pt_z => null()
      pt_u => null(); pt_v => null(); pt_w => null()
    end do
    
  end subroutine moveParticles
end module m_mover
