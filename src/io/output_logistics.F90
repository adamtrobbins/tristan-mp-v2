#include "../defs.F90"

module m_outputlogistics
  use m_globalnamespace
  use m_outputnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_outputnamespace
  use m_readinput, only: getInput
  use m_helpers, only: computeDensity, computeMomentum, interpFromFaces, interpFromEdges
  #ifdef GCA
    use m_helpers, only: computeDensityGCA
  #endif
  use m_exchangearray, only: exchangeArray

  implicit none
contains

  subroutine initializeOutput()
    implicit none
    call getInput('output', 'enable', tot_output_enable, .true.)
    ! individual `.tot.`, `diag` and `spec` outputs
    call getInput('output', 'params_enable', params_enable, .true.)
    call getInput('output', 'prtl_enable', prtl_tot_enable, .true.)
    call getInput('output', 'flds_enable', flds_tot_enable, .true.)
    call getInput('output', 'spec_enable', spectra_enable, .true.)

    call getInput('output', 'diag_enable', diag_enable, .false.)

    call getInput('output', 'start', tot_output_start, 0)
    call getInput('output', 'interval', tot_output_interval, 10)
    call getInput('output', 'stride', tot_output_stride, 10)
    call getInput('output', 'istep', output_flds_istep, 1)
    call getInput('output', 'smooth_window', output_dens_smooth, 2)

    call getInput('output', 'spec_log_bins', spec_log_bins, .true.)
    call getInput('output', 'spec_min', spec_min, 1e-2)
    call getInput('output', 'spec_max', spec_max, 1e2)
    call getInput('output', 'spec_num', spec_num, 100)
    if (spec_log_bins) then
      spec_min = log(spec_min)
      spec_max = log(spec_max)
    endif

    call getInput('output', 'flds_at_prtl', flds_at_prtl_enable, .false.)
    call getInput('output', 'write_xdmf', xdmf_enable, .true.)
    call getInput('output', 'write_nablas', derivatives_enable, .true.)
    call getInput('output', 'write_momenta', momenta_enable, .true.)

    #if defined(HDF5) && defined(MPI08)
      h5comm = MPI_COMM_WORLD%MPI_VAL
      h5info = MPI_INFO_NULL%MPI_VAL
    #elif defined(HDF5) && defined(MPI)
      h5comm = MPI_COMM_WORLD
      h5info = MPI_INFO_NULL
    #endif
  end subroutine initializeOutput

  subroutine initializeSlice()
    implicit none
    integer                 :: i
    character(len=STR_MAX)  :: var_name
    call getInput('slice_output', 'enable', slice_output_enable, .false.)
    call getInput('slice_output', 'start', slice_output_start, 0)
    call getInput('slice_output', 'interval', slice_output_interval, 10)

    #ifndef threeD
      slice_output_enable = .false.
    #endif

    slice_axes(:) = -1
    slice_pos(:) = -1

    do i = 1, 100
      write (var_name, "(A7,I1)") "sliceX_", i
      call getInput('slice_output', var_name, slice_pos(nslices + 1), -1)
      if (slice_pos(nslices + 1) .ne. -1) then
        nslices = nslices + 1
        slice_axes(nslices) = 1
      else
        exit
      end if
    end do

    do i = 1, 100
      write (var_name, "(A7,I1)") "sliceY_", i
      call getInput('slice_output', var_name, slice_pos(nslices + 1), -1)
      if (slice_pos(nslices + 1) .ne. -1) then
        nslices = nslices + 1
        slice_axes(nslices) = 2
      else
        exit
      end if
    end do

    do i = 1, 100
      write (var_name, "(A7,I1)") "sliceZ_", i
      call getInput('slice_output', var_name, slice_pos(nslices + 1), -1)
      if (slice_pos(nslices + 1) .ne. -1) then
        nslices = nslices + 1
        slice_axes(nslices) = 3
      else
        exit
      end if
    end do

  end subroutine initializeSlice

  subroutine prepareOutput()
    ! DEP_PRT [particle-dependent]
    implicit none
    integer                   :: s
    integer                   :: ierr, ndown, pid
    ! initialize particle variables
    n_prtl_vars = 9
    prtl_vars(1:n_prtl_vars) = (/'x    ', 'y    ', 'z    ',&
                               & 'u    ', 'v    ', 'w    ',&
                               & 'wei  ', 'ind  ', 'proc '/)
    prtl_var_types(1:n_prtl_vars) = (/'real ', 'real ', 'real ',&
                                    & 'real ', 'real ', 'real ',&
                                    & 'real ', 'int  ', 'int  '/)
    if (flds_at_prtl_enable) then
      n_prtl_vars = n_prtl_vars + 6
      prtl_vars(10:n_prtl_vars) = (/'ex   ', 'ey   ', 'ez   ',&
                                  & 'bx   ', 'by   ', 'bz   '/)
      prtl_var_types(10:n_prtl_vars) = (/'real ', 'real ', 'real ',&
                                       & 'real ', 'real ', 'real '/)
      do s = 1, nspec
        prtl_vars(n_prtl_vars + s) = 'dens' // STR(s)
        prtl_var_types(n_prtl_vars + s) = 'real '
      end do
      n_prtl_vars = n_prtl_vars + nspec
    end if

    #ifdef PRTLPAYLOADS
      do pid = 1, 3
        prtl_vars(n_prtl_vars + pid) = 'pld' // STR(pid)
        prtl_var_types(n_prtl_vars + pid) = 'real '
      end do
      n_prtl_vars = n_prtl_vars + 3
    #endif

    call prepareSpectraForOutput()
    call defineFieldVarsToOutput()

    ! initialize domain output variables
    !   FIX1: maybe add # of particles per domain
    n_dom_vars = 6
    dom_vars(1 : 6) = (/'x0   ', 'y0   ', 'z0   ',&
                      & 'sx   ', 'sy   ', 'sz   '/)

  end subroutine prepareOutput

  subroutine prepareSpectraForOutput()
    implicit none
    real                      :: energy, u_, v_, w_
    integer                   :: s, i, ti, tj, tk, p, spec_index
    integer                   :: ierr
    real, allocatable, dimension(:,:)     :: spectra
    real, allocatable, dimension(:)       :: send_spec, recv_spec
    #ifdef GCA
      real, allocatable, dimension(:,:)     :: gca_spectra
    #endif

    ! compute spectra
    if (.not. allocated(glob_spectra)) then
      allocate(glob_spectra(nspec, spec_num))
    end if
    allocate(spectra(nspec, spec_num))
    allocate(send_spec(spec_num), recv_spec(spec_num))
    spectra(:,:) = 0

    #ifdef GCA
      if (.not. allocated(glob_gca_spectra)) then
        allocate(glob_gca_spectra(2 * nspec, spec_num))
      end if
      allocate(gca_spectra(2 * nspec, spec_num))
      gca_spectra(:,:) = 0
    #endif

    do s = 1, nspec
     do ti = 1, species(s)%tile_nx
       do tj = 1, species(s)%tile_ny
         do tk = 1, species(s)%tile_nz
           do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
            u_ = species(s)%prtl_tile(ti, tj, tk)%u(p)
            v_ = species(s)%prtl_tile(ti, tj, tk)%v(p)
            w_ = species(s)%prtl_tile(ti, tj, tk)%w(p)
            if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
              energy = sqrt(u_**2 + v_**2 + w_**2)
            else
              energy = sqrt(1.0 + u_**2 + v_**2 + w_**2) - 1.0
            end if
            if (spec_log_bins) energy = log(energy + 1e-8)
            if (energy .le. spec_min) then
              spec_index = 1
            else if (energy .ge. spec_max) then
              spec_index = spec_num
            else
              spec_index = INT(CEILING((energy - spec_min) * REAL(spec_num) / (spec_max - spec_min)))
              if (spec_index .lt. 1) spec_index = 1
              if (spec_index .gt. spec_num) spec_index = spec_num
            end if
            spectra(s, spec_index) = spectra(s, spec_index) + species(s)%prtl_tile(ti, tj, tk)%weight(p)

            #ifdef GCA
              if (species(s)%prtl_tile(ti, tj, tk)%proc(p) .ge. mpi_size) then
                ! particle doing GCA
                gca_spectra(nspec + s, spec_index) = gca_spectra(nspec + s, spec_index) +&
                                                   & species(s)%prtl_tile(ti, tj, tk)%weight(p)
              else
                ! particle doing BORIS
                gca_spectra(s, spec_index) = gca_spectra(s, spec_index) +&
                                           & species(s)%prtl_tile(ti, tj, tk)%weight(p)
              end if
            #endif
           end do
         end do
       end do
     end do
    end do

    ! send to root rank
    do s = 1, nspec
      send_spec(:) = spectra(s,:)
      call MPI_REDUCE(send_spec, recv_spec, spec_num, MPI_REAL,&
                    & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      glob_spectra(s,:) = recv_spec(:)

      #ifdef GCA
        send_spec(:) = gca_spectra(s,:)
        call MPI_REDUCE(send_spec, recv_spec, spec_num, MPI_REAL,&
                      & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        glob_gca_spectra(s,:) = recv_spec(:)

        send_spec(:) = gca_spectra(nspec + s,:)
        call MPI_REDUCE(send_spec, recv_spec, spec_num, MPI_REAL,&
                      & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        glob_gca_spectra(nspec + s,:) = recv_spec(:)
      #endif

      #ifdef RADIATION
        ! compute radiation spectra
        if (allocated(rad_spectra) .and. allocated(glob_rad_spectra)) then
          send_spec(:) = rad_spectra(s,:)
          call MPI_REDUCE(send_spec, recv_spec, spec_num, MPI_REAL,&
                        & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
          glob_rad_spectra(s,:) = recv_spec(:)
          rad_spectra(s,:) = 0.0
        end if
      #endif
    end do

    if (allocated(spectra)) deallocate(spectra)
    if (allocated(send_spec)) deallocate(send_spec)
    if (allocated(recv_spec)) deallocate(recv_spec)
    #ifdef GCA
      if (allocated(gca_spectra)) deallocate(gca_spectra)
    #endif
  end subroutine prepareSpectraForOutput

  subroutine defineFieldVarsToOutput()
    implicit none
    integer     :: s
    ! initialize field variables
    !   total number of fields (excluding particle densities)
    n_fld_vars = 0
    do s = 1, nspec
      fld_vars(0 * nspec + s) = 'dens' // STR(s)
      fld_vars(1 * nspec + s) = 'enrg' // STR(s)
      n_fld_vars = n_fld_vars + 2
      if (momenta_enable) then
        fld_vars(2 * nspec + s) = 'momX' // STR(s)
        fld_vars(3 * nspec + s) = 'momY' // STR(s)
        fld_vars(4 * nspec + s) = 'momZ' // STR(s)
        n_fld_vars = n_fld_vars + 3
      end if
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
    if (derivatives_enable) then
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

end module m_outputlogistics
