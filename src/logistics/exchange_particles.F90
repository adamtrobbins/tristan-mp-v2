#include "../defs.F90"

module m_exchangeparts
  use m_globalnamespace
  use m_communications
  use m_domain
  use m_particles
  implicit none
contains


  subroutine exchangeParticles()
    implicit none
    real(mprec), pointer, contiguous   :: pt_x(:), pt_y(:), pt_z(:)
    integer :: s, p

    ! send -x | receive +x
    do s = 1, nspec
      pt_x => sp_(s)%x; pt_y => sp_(s)%y; pt_z => sp_(s)%z
      !$omp simd
      !dir$ vector aligned
      do p = 1, spp_(s)%npart_sp
        if (sp_)
      end do
    end do

  end subroutine exchangeParticles
end module m_exchangeparts
