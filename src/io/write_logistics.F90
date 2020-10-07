#include "../defs.F90"

module m_writelogistics
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers, only: computeDensity, computeMomentum, interpFromFaces, interpFromEdges
  #ifdef GCA
    use m_helpers, only: computeDensityGCA
  #endif
  use m_exchangearray, only: exchangeArray
  implicit none

  ! input parameters
  integer                   :: output_flds_istep, output_dens_smooth
  logical                   :: write_derivatives

  character(len=STR_MAX)    :: fld_vars(100)
  integer                   :: n_fld_vars
contains
  subroutine defineFieldVarsToOutput()
    implicit none
    integer     :: s
    ! initialize field variables
    !   total number of fields (excluding particle densities)
    n_fld_vars = 0
    do s = 1, nspec
      fld_vars(0 * nspec + s) = 'dens' // STR(s)
      fld_vars(1 * nspec + s) = 'enrg' // STR(s)
      fld_vars(2 * nspec + s) = 'momX' // STR(s)
      fld_vars(3 * nspec + s) = 'momY' // STR(s)
      fld_vars(4 * nspec + s) = 'momZ' // STR(s)
      n_fld_vars = n_fld_vars + 5
      #ifdef GCA
        fld_vars(5 * nspec + s) = 'dgca' // STR(s)
        n_fld_vars = n_fld_vars + 1
      #endif
    end do

    fld_vars(n_fld_vars + 1 : n_fld_vars + 1 + 12) =&
                                 & (/'ex   ', 'ey   ', 'ez   ',&
                                   & 'bx   ', 'by   ', 'bz   ',&
                                   & 'jx   ', 'jy   ', 'jz   ',&
                                   & 'xx   ', 'yy   ', 'zz   '/)
    n_fld_vars = n_fld_vars + 12
    if (write_derivatives) then
      fld_vars(n_fld_vars + 1 : n_fld_vars + 1 + 4) = (/'curlBx', 'curlBy', 'curlBz', 'divE'/)
      n_fld_vars = n_fld_vars + 4
    end if
  end subroutine defineFieldVarsToOutput

  ! writes a field specified by `fld_var` from gridcell `i,j,k` ...
  ! ... to `sm_arr(i1, j1, k1)` with proper interpolation etc for the output
  subroutine selectFieldForOutput(fld_var, i1, j1, k1, i, j, k, writing_lgarrQ)
    implicit none
    character(len=STR_MAX), intent(in)  :: fld_var
    integer(kind=2), intent(in)         :: i1, j1, k1, i, j, k
    logical, intent(in)                 :: writing_lgarrQ
    real                                :: ex0, ey0, ez0, bx0, by0, bz0, jx0, jy0, jz0
    real                                :: dx1, dx2, dy1, dy2, dz1, dz2, divE
    select case (trim(fld_var))
    case('ex')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ex0 = ex(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ex0 * B_norm
    case('ey')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ey0 = ey(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ey0 * B_norm
    case('ez')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
      #else
        ez0 = ez(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = ez0 * B_norm
    case('bx')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        bx0 = bx(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = bx0 * B_norm
    case('by')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        by0 = by(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = by0 * B_norm
    case('bz')
      #ifndef DEBUG
        call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
      #else
        bz0 = bz(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = bz0 * B_norm
    case('jx')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jx0 = jx(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jx0 * B_norm
    case('jy')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jy0 = jy(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jy0 * B_norm
    case('jz')
      #ifndef DEBUG
        call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
      #else
        jz0 = jz(i, j, k)
      #endif
      sm_arr(i1, j1, k1) = -jz0 * B_norm
    case('curlBx')
      #ifdef oneD
        dx1 = 0.0; dx2 = 0.0
      #elif twoD
        dx1 = (bz(    i,    j,k) - bz(    i,j - 1,    k))
        dx2 = (bz(i - 1,    j,k) - bz(i - 1,j - 1,    k))
      #elif threeD
        dx1 = (bz(    i,    j,    k) - bz(    i,j - 1,    k)) - (by(    i,    j,    k) - by(    i,    j,k - 1))
        dx2 = (bz(i - 1,    j,    k) - bz(i - 1,j - 1,    k)) - (by(i - 1,    j,    k) - by(i - 1,    j,k - 1))
      #endif
      sm_arr(i1, j1, k1) = B_norm * 0.5 * (dx1 + dx2)
    case('curlBy')
      #ifdef oneD
        dy1 = -(bz(i, j, k) - bz(i - 1, j, k))
        dy2 = dy1
      #elif twoD
        dy1 = -(bz(   i,    j,    k) - bz(i - 1,    j,    k))
        dy2 = -(bz(   i,j - 1,    k) - bz(i - 1,j - 1,    k))
      #elif threeD
        dy1 = (bx(    i,    j,    k) - bx(    i,    j,k - 1)) - (bz(    i,    j,    k) - bz(i - 1,    j,    k))
        dy2 = (bx(    i,j - 1,    k) - bx(    i,j - 1,k - 1)) - (bz(    i,j - 1,    k) - bz(i - 1,j - 1,    k))
      #endif
      sm_arr(i1, j1, k1) = B_norm * 0.5 * (dy1 + dy2)
    case('curlBz')
      #ifdef oneD
        dz1 = (by(i, j, k) - by(i - 1, j, k))
        dz2 = dz1
      #elif twoD
        dz1 = (by(i, j, k) - by(i - 1, j, k)) - (bx(i, j, k) - bx(i, j - 1, k))
        dz2 = dz1
      #elif threeD
        dz1 = (by(    i,    j,    k) - by(i - 1,    j,    k)) - (bx(    i,    j,    k) - bx(    i,j - 1,    k))
        dz2 = (by(    i,    j,k - 1) - by(i - 1,    j,k - 1)) - (bx(    i,    j,k - 1) - bx(    i,j - 1,k - 1))
      #endif
      sm_arr(i1, j1, k1) = B_norm * 0.5 * (dz1 + dz2)
    case('divE')
      divE = 0.0
      #if defined(oneD) || defined (twoD) || defined (threeD)
        divE = divE + (ex(i, j, k) - ex(i - 1, j, k))
      #endif
      #if defined(twoD) || defined (threeD)
        divE = divE + (ey(i, j, k) - ey(i, j - 1, k))
      #endif
      #if defined(threeD)
        divE = divE + (ez(i, j, k) - ez(i, j, k - 1))
      #endif
      sm_arr(i1, j1, k1) = B_norm * divE
    case('xx')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%x0 + i, 4)
    case('yy')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%y0 + j, 4)
    case('zz')
      sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%z0 + k, 4)
    case default
      if (((fld_var(1:4) .ne. 'dens') .and.&
         & (fld_var(1:4) .ne. 'enrg') .and.&
         & (fld_var(1:3) .ne. 'mom') .and.&
         & (fld_var(1:4) .ne. 'dgca')) .or.&
         & (.not. writing_lgarrQ)) then
        call throwError("ERROR: unrecognized `fldname`")
      else
        ! interpolating cell-centered values to nodes
        #ifdef oneD
          sm_arr(i1, j1, k1) = 0.5 * (lg_arr(i, j, k) + lg_arr(i - 1, j, k))
        #elif twoD
          sm_arr(i1, j1, k1) = 0.25 * (lg_arr(i, j, k) + lg_arr(i - 1, j, k) +&
                                     & lg_arr(i - 1, j - 1, k) + lg_arr(i, j - 1, k))
        #elif threeD
          sm_arr(i1, j1, k1) = 0.125 * (lg_arr(i, j, k) + lg_arr(i - 1, j - 1, k - 1) +&
                                      & lg_arr(i - 1, j, k) + lg_arr(i, j - 1, k) + lg_arr(i, j, k - 1) +&
                                      & lg_arr(i - 1, j - 1, k) + lg_arr(i, j - 1, k - 1) + lg_arr(i - 1, j, k - 1))
        #endif

      end if
    end select
  end subroutine selectFieldForOutput

  subroutine prepareFieldForOutput(fldname, writing_lgarrQ)
    implicit none
    character(len=STR_MAX), intent(in)    :: fldname
    logical, intent(out)                  :: writing_lgarrQ
    integer                               :: s

    if (fldname(1:4) .eq. 'dens') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      ! fill `lg_arr` with density of species `s`
      call computeDensity(s, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'enrg') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      ! fill `lg_arr` with energy density of species `s`
      call computeMomentum(s, 0, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'dgca') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      #ifndef GCA
        call throwError('ERROR: `dgca` not defined without GCA flag.')
      #else
        call computeDensityGCA(s, reset=.true., ds=output_dens_smooth)
        call exchangeArray()
      #endif
    else if (fldname(1:4) .eq. 'momX') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      call computeMomentum(s, 1, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'momY') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      call computeMomentum(s, 2, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else if (fldname(1:4) .eq. 'momZ') then
      writing_lgarrQ = .true.
      s = STRtoINT(fldname(5:5))
      call computeMomentum(s, 3, reset=.true., ds=output_dens_smooth)
      call exchangeArray()
    else
      writing_lgarrQ = .false.
    end if
  end subroutine prepareFieldForOutput

end module m_writelogistics
