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
    integer                             :: ind_t(3), indToRnk
    ind_t = ind
    if (boundary_x .eq. 1) then
      if (ind_t(1) .lt. 0) ind_t(1) = ind_t(1) + sizex
      ind_t(1) = modulo(ind_t(1), sizex)
    end if
    if (boundary_y .eq. 1) then
      if (ind_t(2) .lt. 0) ind_t(2) = ind_t(2) + sizey
      ind_t(2) = modulo(ind_t(2), sizey)
    end if
    if (boundary_z .eq. 1) then
      if (ind_t(3) .lt. 0) ind_t(3) = ind_t(3) + sizez
      ind_t(3) = modulo(ind_t(3), sizez)
    end if
    if ((ind_t(1) .lt. 0) .or. (ind_t(1) .ge. sizex) .or.&
      & (ind_t(2) .lt. 0) .or. (ind_t(2) .ge. sizey) .or.&
      & (ind_t(3) .lt. 0) .or. (ind_t(3) .ge. sizez)) then
      indToRnk = -1
    else
      indToRnk = ind_t(3) * sizex * sizey + ind_t(2) * sizex + ind_t(1)
    end if
  end function indToRnk

  subroutine interpFlds(dx, dy, dz, i, j, k, ex0, ey0, ez0, bx0, by0, bz0)
    implicit none
    integer(kind=2), intent(in) :: i, j, k
    real, intent(in)            :: dx, dy, dz
    real, intent(out)           :: ex0, ey0, ez0, bx0, by0, bz0
    real                        :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                 & c00, c01, c10, c11, c0, c1

    ! e_x
    c000 = 0.5 * (ex(     i,      j,      k) + ex(  i - 1,      j,      k))
    c100 = 0.5 * (ex(     i,      j,      k) + ex(  i + 1,      j,      k))
    c010 = 0.5 * (ex(     i,  j + 1,      k) + ex(  i - 1,  j + 1,      k))
    c110 = 0.5 * (ex(     i,  j + 1,      k) + ex(  i + 1,  j + 1,      k))
    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    #ifndef threeD
      ex0 = c0
    #else
      c001 = 0.5 * (ex(     i,      j,  k + 1) + ex(  i - 1,      j,  k + 1))
      c101 = 0.5 * (ex(     i,      j,  k + 1) + ex(  i + 1,      j,  k + 1))
      c011 = 0.5 * (ex(     i,  j + 1,  k + 1) + ex(  i - 1,  j + 1,  k + 1))
      c111 = 0.5 * (ex(     i,  j + 1,  k + 1) + ex(  i + 1,  j + 1,  k + 1))
      c01 = c001 * (1 - dx) + c101 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c1 = c01 * (1 - dy) + c11 * dy
      ex0 = c0 * (1 - dz) + c1 * dz
    #endif

    ! e_y
    c000 = 0.5 * (ey(     i,      j,      k) + ey(      i,  j - 1,      k))
    c100 = 0.5 * (ey( i + 1,      j,      k) + ey(  i + 1,  j - 1,      k))
    c010 = 0.5 * (ey(     i,      j,      k) + ey(      i,  j + 1,      k))
    c110 = 0.5 * (ey( i + 1,      j,      k) + ey(  i + 1,  j + 1,      k))
    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    #ifndef threeD
      ey0 = c0
    #else
      c001 = 0.5 * (ey(     i,      j,  k + 1) + ey(      i,  j - 1,  k + 1))
      c101 = 0.5 * (ey( i + 1,      j,  k + 1) + ey(  i + 1,  j - 1,  k + 1))
      c011 = 0.5 * (ey(     i,      j,  k + 1) + ey(      i,  j + 1,  k + 1))
      c111 = 0.5 * (ey( i + 1,      j,  k + 1) + ey(  i + 1,  j + 1,  k + 1))
      c01 = c001 * (1 - dx) + c101 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c1 = c01 * (1 - dy) + c11 * dy
      ey0 = c0 * (1 - dz) + c1 * dz
    #endif

    ! e_z
    #ifndef threeD
      c000 = ez(     i,      j,      k)
      c100 = ez( i + 1,      j,      k)
      c010 = ez(     i,  j + 1,      k)
      c110 = ez( i + 1,  j + 1,      k)
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      ez0 = c00 * (1 - dy) + c10 * dy
    #else
      c000 = 0.5 * (ez(     i,      j,      k) + ez(      i,      j,  k - 1))
      c100 = 0.5 * (ez( i + 1,      j,      k) + ez(  i + 1,      j,  k - 1))
      c010 = 0.5 * (ez(     i,  j + 1,      k) + ez(      i,  j + 1,  k - 1))
      c110 = 0.5 * (ez( i + 1,  j + 1,      k) + ez(  i + 1,  j + 1,  k - 1))
      c001 = 0.5 * (ez(     i,      j,      k) + ez(      i,      j,  k + 1))
      c101 = 0.5 * (ez( i + 1,      j,      k) + ez(  i + 1,      j,  k + 1))
      c011 = 0.5 * (ez(     i,  j + 1,      k) + ez(      i,  j + 1,  k + 1))
      c111 = 0.5 * (ez( i + 1,  j + 1,      k) + ez(  i + 1,  j + 1,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      ez0 = c0 * (1 - dz) + c1 * dz
    #endif

    ! b_x
    #ifndef threeD
      c000 = 0.5 * (bx(			i,			j,			k) + bx(			i,	j - 1,			k))
      c100 = 0.5 * (bx(	i + 1,			j,			k) + bx(	i + 1,	j - 1,			k))
      c010 = 0.5 * (bx(			i,			j,			k) + bx(			i,	j + 1,			k))
      c110 = 0.5 * (bx(	i + 1,			j,			k) + bx(	i + 1,	j + 1,			k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      bx0 = c00 * (1 - dy) + c10 * dy
    #else
      c000 = 0.25 * (bx(			i,			j,			k) + bx(			i,	j - 1,			k) + bx(			i,			j,	k - 1) + bx(			i,	j - 1,	k - 1))
      c100 = 0.25 * (bx(	i + 1,			j,			k) + bx(	i + 1,	j - 1,			k) + bx(	i + 1,			j,	k - 1) + bx(	i + 1,	j - 1,	k - 1))
      c001 = 0.25 * (bx(			i,			j,			k) + bx(			i,			j,	k + 1) + bx(			i,	j - 1,			k) + bx(			i,	j - 1,	k + 1))
      c101 = 0.25 * (bx(	i + 1,			j,			k) + bx(	i + 1,			j,	k + 1) + bx(	i + 1,	j - 1,			k) + bx(	i + 1,	j - 1,	k + 1))
      c010 = 0.25 * (bx(			i,			j,			k) + bx(			i,	j + 1,			k) + bx(			i,			j,	k - 1) + bx(			i,	j + 1,	k - 1))
      c110 = 0.25 * (bx(	i + 1,			j,			k) + bx(	i + 1,			j,	k - 1) + bx(	i + 1,	j + 1,	k - 1) + bx(	i + 1,	j + 1,			k))
      c011 = 0.25 * (bx(			i,			j,			k) + bx(			i,	j + 1,			k) + bx(			i,	j + 1,	k + 1) + bx(			i,			j,	k + 1))
      c111 = 0.25 * (bx(	i + 1,			j,			k) + bx(	i + 1,	j + 1,			k) + bx(	i + 1,	j + 1,	k + 1) + bx(	i + 1,			j,	k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      bx0 = c0 * (1 - dz) + c1 * dz
    #endif

    ! b_y
    #ifndef threeD
      c000 = 0.5 * (by(	i - 1,			j,			k) + by(			i,			j,			k))
      c100 = 0.5 * (by(			i,			j,			k) + by(	i + 1,			j,			k))
      c010 = 0.5 * (by(	i - 1,	j + 1,			k) + by(			i,	j + 1,			k))
      c110 = 0.5 * (by(			i,	j + 1,			k) + by(	i + 1,	j + 1,			k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      by0 = c00 * (1 - dy) + c10 * dy
    #else
      c000 = 0.25 * (by(	i - 1,			j,	k - 1) + by(	i - 1,			j,			k) + by(			i,			j,	k - 1) + by(			i,			j,			k))
      c100 = 0.25 * (by(			i,			j,	k - 1) + by(			i,			j,			k) + by(	i + 1,			j,	k - 1) + by(	i + 1,			j,			k))
      c001 = 0.25 * (by(	i - 1,			j,			k) + by(	i - 1,			j,	k + 1) + by(			i,			j,			k) + by(			i,			j,	k + 1))
      c101 = 0.25 * (by(			i,			j,			k) + by(			i,			j,	k + 1) + by(	i + 1,			j,			k) + by(	i + 1,			j,	k + 1))
      c010 = 0.25 * (by(	i - 1,	j + 1,	k - 1) + by(	i - 1,	j + 1,			k) + by(			i,	j + 1,	k - 1) + by(			i,	j + 1,			k))
      c110 = 0.25 * (by(			i,	j + 1,	k - 1) + by(			i,	j + 1,			k) + by(	i + 1,	j + 1,	k - 1) + by(	i + 1,	j + 1,			k))
      c011 = 0.25 * (by(	i - 1,	j + 1,			k) + by(	i - 1,	j + 1,	k + 1) + by(			i,	j + 1,			k) + by(			i,	j + 1,	k + 1))
      c111 = 0.25 * (by(			i,	j + 1,			k) + by(			i,	j + 1,	k + 1) + by(	i + 1,	j + 1,			k) + by(	i + 1,	j + 1,	k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      by0 = c0 * (1 - dz) + c1 * dz
    #endif

    ! b_z
    #ifndef threeD
      c000 = 0.25 * (bz(	i - 1,	j - 1,			k) + bz(	i - 1,			j,			k) + bz(			i,	j - 1,			k) + bz(			i,			j,			k))
      c100 = 0.25 * (bz(			i,	j - 1,			k) + bz(			i,			j,			k) + bz(	i + 1,	j - 1,			k) + bz(	i + 1,			j,			k))
      c010 = 0.25 * (bz(	i - 1,			j,			k) + bz(	i - 1,	j + 1,			k) + bz(			i,			j,			k) + bz(			i,	j + 1,			k))
      c110 = 0.25 * (bz(			i,			j,			k) + bz(			i,	j + 1,			k) + bz(	i + 1,			j,			k) + bz(	i + 1,	j + 1,			k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      bz0 = c00 * (1 - dy) + c10 * dy
    #else
      c000 = 0.25 * (bz(	i - 1,	j - 1,			k) + bz(	i - 1,			j,			k) + bz(			i,	j - 1,			k) + bz(			i,			j,			k))
      c100 = 0.25 * (bz(			i,	j - 1,			k) + bz(			i,			j,			k) + bz(	i + 1,	j - 1,			k) + bz(	i + 1,			j,			k))
      c001 = 0.25 * (bz(	i - 1,	j - 1,	k + 1) + bz(	i - 1,			j,	k + 1) + bz(			i,	j - 1,	k + 1) + bz(			i,			j,	k + 1))
      c101 = 0.25 * (bz(			i,	j - 1,	k + 1) + bz(			i,			j,	k + 1) + bz(	i + 1,	j - 1,	k + 1) + bz(	i + 1,			j,	k + 1))
      c010 = 0.25 * (bz(	i - 1,			j,			k) + bz(	i - 1,	j + 1,			k) + bz(			i,			j,			k) + bz(			i,	j + 1,			k))
      c110 = 0.25 * (bz(			i,			j,			k) + bz(			i,	j + 1,			k) + bz(	i + 1,			j,			k) + bz(	i + 1,	j + 1,			k))
      c011 = 0.25 * (bz(	i - 1,			j,	k + 1) + bz(	i - 1,	j + 1,	k + 1) + bz(			i,			j,	k + 1) + bz(			i,	j + 1,	k + 1))
      c111 = 0.25 * (bz(			i,			j,	k + 1) + bz(			i,	j + 1,	k + 1) + bz(	i + 1,			j,	k + 1) + bz(	i + 1,	j + 1,	k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      bz0 = c0 * (1 - dz) + c1 * dz
    #endif
  end subroutine
end module m_helpers
