#include "../defs.F90"

module m_helpers
  use m_globalnamespace
  use m_domain
  use m_particles
  use m_fields
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

  subroutine computeDensity(s)
    implicit none
    integer, intent(in)                   :: s
    integer                               :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    scalar_int_array(:,:,:) = 0
    do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
        do tk = 1, species(s)%tile_nz
          pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
          pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
          pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
          ! FIX1 vectorize/align
          do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
            scalar_int_array(pt_xi(p), pt_yi(p), pt_zi(p)) = scalar_int_array(pt_xi(p), pt_yi(p), pt_zi(p)) + 1
          end do
          pt_xi => null(); pt_yi => null(); pt_zi => null()
        end do
      end do
    end do
  end subroutine computeDensity

  subroutine interpFromEdges(dx, dy, dz, i, j, k, &
                           & fx, fy, fz, &
                           & intfx, intfy, intfz)
    implicit none
    integer(kind=2), intent(in)   :: i, j, k
    real, intent(in)              :: dx, dy, dz
    real, intent(in)              :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(in)              :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(in)              :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(out)             :: intfx, intfy, intfz
    real                          :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                   & c00, c01, c10, c11, c0, c1
    ! f_x
    c000 = 0.5 * (fx(     i,      j,      k) + fx(  i - 1,      j,      k))
    c100 = 0.5 * (fx(     i,      j,      k) + fx(  i + 1,      j,      k))
    c010 = 0.5 * (fx(     i,  j + 1,      k) + fx(  i - 1,  j + 1,      k))
    c110 = 0.5 * (fx(     i,  j + 1,      k) + fx(  i + 1,  j + 1,      k))
    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    #ifndef threeD
      intfx = c0
    #else
      c001 = 0.5 * (fx(     i,      j,  k + 1) + fx(  i - 1,      j,  k + 1))
      c101 = 0.5 * (fx(     i,      j,  k + 1) + fx(  i + 1,      j,  k + 1))
      c011 = 0.5 * (fx(     i,  j + 1,  k + 1) + fx(  i - 1,  j + 1,  k + 1))
      c111 = 0.5 * (fx(     i,  j + 1,  k + 1) + fx(  i + 1,  j + 1,  k + 1))
      c01 = c001 * (1 - dx) + c101 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c1 = c01 * (1 - dy) + c11 * dy
      intfx = c0 * (1 - dz) + c1 * dz
    #endif

    ! f_y
    c000 = 0.5 * (fy(     i,      j,      k) + fy(      i,  j - 1,      k))
    c100 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,  j - 1,      k))
    c010 = 0.5 * (fy(     i,      j,      k) + fy(      i,  j + 1,      k))
    c110 = 0.5 * (fy( i + 1,      j,      k) + fy(  i + 1,  j + 1,      k))
    c00 = c000 * (1 - dx) + c100 * dx
    c10 = c010 * (1 - dx) + c110 * dx
    c0 = c00 * (1 - dy) + c10 * dy
    #ifndef threeD
      intfy = c0
    #else
      c001 = 0.5 * (fy(     i,      j,  k + 1) + fy(      i,  j - 1,  k + 1))
      c101 = 0.5 * (fy( i + 1,      j,  k + 1) + fy(  i + 1,  j - 1,  k + 1))
      c011 = 0.5 * (fy(     i,      j,  k + 1) + fy(      i,  j + 1,  k + 1))
      c111 = 0.5 * (fy( i + 1,      j,  k + 1) + fy(  i + 1,  j + 1,  k + 1))
      c01 = c001 * (1 - dx) + c101 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c1 = c01 * (1 - dy) + c11 * dy
      intfy = c0 * (1 - dz) + c1 * dz
    #endif

    ! f_z
    #ifndef threeD
      c000 = fz(     i,      j,      k)
      c100 = fz( i + 1,      j,      k)
      c010 = fz(     i,  j + 1,      k)
      c110 = fz( i + 1,  j + 1,      k)
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfz = c00 * (1 - dy) + c10 * dy
    #else
      c000 = 0.5 * (fz(     i,      j,      k) + fz(      i,      j,  k - 1))
      c100 = 0.5 * (fz( i + 1,      j,      k) + fz(  i + 1,      j,  k - 1))
      c010 = 0.5 * (fz(     i,  j + 1,      k) + fz(      i,  j + 1,  k - 1))
      c110 = 0.5 * (fz( i + 1,  j + 1,      k) + fz(  i + 1,  j + 1,  k - 1))
      c001 = 0.5 * (fz(     i,      j,      k) + fz(      i,      j,  k + 1))
      c101 = 0.5 * (fz( i + 1,      j,      k) + fz(  i + 1,      j,  k + 1))
      c011 = 0.5 * (fz(     i,  j + 1,      k) + fz(      i,  j + 1,  k + 1))
      c111 = 0.5 * (fz( i + 1,  j + 1,      k) + fz(  i + 1,  j + 1,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfz = c0 * (1 - dz) + c1 * dz
    #endif
  end subroutine interpFromEdges

  subroutine interpFromFaces(dx, dy, dz, i, j, k, &
                           & fx, fy, fz, &
                           & intfx, intfy, intfz)
    implicit none
    integer(kind=2), intent(in)   :: i, j, k
    real, intent(in)              :: dx, dy, dz
    real, intent(in)              :: fx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(in)              :: fy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(in)              :: fz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                                      & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                                      & fldBoundZ)
    real, intent(out)             :: intfx, intfy, intfz
    real                          :: c000, c100, c001, c101, c010, c110, c011, c111,&
                                   & c00, c01, c10, c11, c0, c1
    ! f_x
    #ifndef threeD
      c000 = 0.5 * (fx(      i,      j,      k) + fx(      i,  j - 1,      k))
      c100 = 0.5 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j - 1,      k))
      c010 = 0.5 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k))
      c110 = 0.5 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfx = c00 * (1 - dy) + c10 * dy
    #else
      c000 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j - 1,      k) +&
                   & fx(      i,      j,  k - 1) + fx(      i,  j - 1,  k - 1))
      c100 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j - 1,      k) +&
                   & fx(  i + 1,      j,  k - 1) + fx(  i + 1,  j - 1,  k - 1))
      c001 = 0.25 * (fx(      i,      j,      k) + fx(      i,      j,  k + 1) +&
                   & fx(      i,  j - 1,      k) + fx(      i,  j - 1,  k + 1))
      c101 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,      j,  k + 1) +&
                   & fx(  i + 1,  j - 1,      k) + fx(  i + 1,  j - 1,  k + 1))
      c010 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k) +&
                   & fx(      i,      j,  k - 1) + fx(      i,  j + 1,  k - 1))
      c110 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,      j,  k - 1) +&
                   & fx(  i + 1,  j + 1,  k - 1) + fx(  i + 1,  j + 1,      k))
      c011 = 0.25 * (fx(      i,      j,      k) + fx(      i,  j + 1,      k) +&
                   & fx(      i,  j + 1,  k + 1) + fx(      i,      j,  k + 1))
      c111 = 0.25 * (fx(  i + 1,      j,      k) + fx(  i + 1,  j + 1,      k) +&
                   & fx(  i + 1,  j + 1,  k + 1) + fx(  i + 1,      j,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfx = c0 * (1 - dz) + c1 * dz
    #endif

    ! b_y
    #ifndef threeD
      c000 = 0.5 * (fy(  i - 1,      j,      k) + fy(      i,      j,      k))
      c100 = 0.5 * (fy(      i,      j,      k) + fy(  i + 1,      j,      k))
      c010 = 0.5 * (fy(  i - 1,  j + 1,      k) + fy(      i,  j + 1,      k))
      c110 = 0.5 * (fy(      i,  j + 1,      k) + fy(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfy = c00 * (1 - dy) + c10 * dy
    #else
      c000 = 0.25 * (fy(  i - 1,      j,  k - 1) + fy(  i - 1,      j,      k) +&
                   & fy(      i,      j,  k - 1) + fy(      i,      j,      k))
      c100 = 0.25 * (fy(      i,      j,  k - 1) + fy(      i,      j,      k) +&
                   & fy(  i + 1,      j,  k - 1) + fy(  i + 1,      j,      k))
      c001 = 0.25 * (fy(  i - 1,      j,      k) + fy(  i - 1,      j,  k + 1) +&
                   & fy(      i,      j,      k) + fy(      i,      j,  k + 1))
      c101 = 0.25 * (fy(      i,      j,      k) + fy(      i,      j,  k + 1) +&
                   & fy(  i + 1,      j,      k) + fy(  i + 1,      j,  k + 1))
      c010 = 0.25 * (fy(  i - 1,  j + 1,  k - 1) + fy(  i - 1,  j + 1,      k) +&
                   & fy(      i,  j + 1,  k - 1) + fy(      i,  j + 1,      k))
      c110 = 0.25 * (fy(      i,  j + 1,  k - 1) + fy(      i,  j + 1,      k) +&
                   & fy(  i + 1,  j + 1,  k - 1) + fy(  i + 1,  j + 1,      k))
      c011 = 0.25 * (fy(  i - 1,  j + 1,      k) + fy(  i - 1,  j + 1,  k + 1) +&
                   & fy(      i,  j + 1,      k) + fy(      i,  j + 1,  k + 1))
      c111 = 0.25 * (fy(      i,  j + 1,      k) + fy(      i,  j + 1,  k + 1) +&
                   & fy(  i + 1,  j + 1,      k) + fy(  i + 1,  j + 1,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfy = c0 * (1 - dz) + c1 * dz
    #endif

    ! b_z
    #ifndef threeD
      c000 = 0.25 * (fz(  i - 1,  j - 1,      k) + fz(  i - 1,      j,      k) +&
                   & fz(      i,  j - 1,      k) + fz(      i,      j,      k))
      c100 = 0.25 * (fz(      i,  j - 1,      k) + fz(      i,      j,      k) +&
                   & fz(  i + 1,  j - 1,      k) + fz(  i + 1,      j,      k))
      c010 = 0.25 * (fz(  i - 1,      j,      k) + fz(  i - 1,  j + 1,      k) +&
                   & fz(      i,      j,      k) + fz(      i,  j + 1,      k))
      c110 = 0.25 * (fz(      i,      j,      k) + fz(      i,  j + 1,      k) +&
                   & fz(  i + 1,      j,      k) + fz(  i + 1,  j + 1,      k))
      c00 = c000 * (1 - dx) + c100 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      intfz = c00 * (1 - dy) + c10 * dy
    #else
      c000 = 0.25 * (fz(  i - 1,  j - 1,      k) + fz(  i - 1,      j,      k) +&
                   & fz(      i,  j - 1,      k) + fz(      i,      j,      k))
      c100 = 0.25 * (fz(      i,  j - 1,      k) + fz(      i,      j,      k) +&
                   & fz(  i + 1,  j - 1,      k) + fz(  i + 1,      j,      k))
      c001 = 0.25 * (fz(  i - 1,  j - 1,  k + 1) + fz(  i - 1,      j,  k + 1) +&
                   & fz(      i,  j - 1,  k + 1) + fz(      i,      j,  k + 1))
      c101 = 0.25 * (fz(      i,  j - 1,  k + 1) + fz(      i,      j,  k + 1) +&
                   & fz(  i + 1,  j - 1,  k + 1) + fz(  i + 1,      j,  k + 1))
      c010 = 0.25 * (fz(  i - 1,      j,      k) + fz(  i - 1,  j + 1,      k) +&
                   & fz(      i,      j,      k) + fz(      i,  j + 1,      k))
      c110 = 0.25 * (fz(      i,      j,      k) + fz(      i,  j + 1,      k) +&
                   & fz(  i + 1,      j,      k) + fz(  i + 1,  j + 1,      k))
      c011 = 0.25 * (fz(  i - 1,      j,  k + 1) + fz(  i - 1,  j + 1,  k + 1) +&
                   & fz(      i,      j,  k + 1) + fz(      i,  j + 1,  k + 1))
      c111 = 0.25 * (fz(      i,      j,  k + 1) + fz(      i,  j + 1,  k + 1) +&
                   & fz(  i + 1,      j,  k + 1) + fz(  i + 1,  j + 1,  k + 1))
      c00 = c000 * (1 - dx) + c100 * dx
      c01 = c001 * (1 - dx) + c101 * dx
      c10 = c010 * (1 - dx) + c110 * dx
      c11 = c011 * (1 - dx) + c111 * dx
      c0 = c00 * (1 - dy) + c10 * dy
      c1 = c01 * (1 - dy) + c11 * dy
      intfz = c0 * (1 - dz) + c1 * dz
    #endif
  end subroutine interpFromFaces
end module m_helpers
