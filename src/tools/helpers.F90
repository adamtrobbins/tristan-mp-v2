#include "../defs.F90"

module m_helpers
  use m_globalnamespace
  use m_domain
  use m_particles
  use m_fields
  use m_communications
  implicit none
contains
  function rnkToInd(rnk)
    implicit none
    integer, intent(in)   :: rnk
    integer, dimension(3) :: rnkToInd
    if ((rnk .lt. 0) .or. (rnk .ge. mpi_size)) then
      rnkToInd = (/-1, -1, -1/)
    else
      rnkToInd(3) = rnk / (sizex * sizey)
      rnkToInd(2) = (rnk - sizex * sizey * rnkToInd(3)) / sizex
      rnkToInd(1) = rnk - sizex * sizey * rnkToInd(3) - sizex * rnkToInd(2)
    end if
  end function rnkToInd

  function indToRnk(ind)
    implicit none
    integer, intent(in)                 :: ind(3)
    integer                             :: ind_(3), indToRnk
    ind_ = ind
    if (boundary_x .eq. 1) ind_(1) = modulo(ind_(1), sizex)
    if (boundary_y .eq. 1) ind_(2) = modulo(ind_(2), sizey)
    if (boundary_z .eq. 1) ind_(3) = modulo(ind_(3), sizez)
    if ((ind_(1) .lt. 0) .or. (ind_(1) .ge. sizex) .or.&
      & (ind_(2) .lt. 0) .or. (ind_(2) .ge. sizey) .or.&
      & (ind_(3) .lt. 0) .or. (ind_(3) .ge. sizez)) then
      indToRnk = -1
    else
      indToRnk = ind_(3) * sizex * sizey + ind_(2) * sizex + ind_(1)
    end if
  end function indToRnk

  subroutine interpFlds(dx, dy, dz, i, j, k, ex0, ey0, ez0, bx0, by0, bz0)
    implicit none
    integer, intent(in) :: i, j, k
    real, intent(in) :: dx, dy, dz
    real, intent(out) :: ex0, ey0, ez0, bx0, by0, bz0
    
  end subroutine
end module m_helpers
