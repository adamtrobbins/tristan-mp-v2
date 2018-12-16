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
    real(mprec), pointer, contiguous   :: pt_x(:), pt_y(:), pt_z(:), pt_u(:), pt_v(:), pt_w(:)
    do s = 1, nspec
      pt_x => sp_(s)%x
      pt_y => sp_(s)%y
      pt_z => sp_(s)%z
      pt_u => sp_(s)%u
      pt_v => sp_(s)%v
      pt_w => sp_(s)%w
      if (spp_(s)%m_sp .eq. 0) then
        ! routine for massless particles
        !dir$ assume_aligned pt_x:64, pt_y:64, pt_z:64, pt_u:64, pt_v:64, pt_w:64
        pt_x = pt_x + CC * pt_u
        pt_y = pt_y + CC * pt_v
        pt_z = pt_z + CC * pt_w
      else
        ! routine for massive particles
      end if
    end do
  end subroutine moveParticles
end module m_mover
