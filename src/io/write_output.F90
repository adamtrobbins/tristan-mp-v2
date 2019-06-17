#include "../defs.F90"

module m_writeoutput
  #ifdef HDF5
    use hdf5
  #endif
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  implicit none

  integer                 :: output_start, output_interval, output_stride, output_istep
  integer                 :: n_fld_vars, n_prtl_vars, n_dom_vars
  character(len=STR_MAX)  :: prtl_vars(100), prtl_var_types(100), fld_vars(100), dom_vars(100)
  integer, allocatable, dimension(:,:) :: glob_spectra


  !--- PRIVATE functions -----------------------------------------!
  #ifndef HDF5
    private :: writeParticles_Binary, writeFields_Binary,&
             & writeSpectra_Binary
  #else
    private :: writeParticles_hdf5, writeFields_hdf5,&
             & writeSpectra_hdf5, writeDomain_hdf5
  #endif
  private :: initializeOutput
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: n_fld_vars, n_prtl_vars
  private :: prtl_vars, prtl_var_types
  private :: fld_vars
  !...............................................................!
contains
  subroutine writeOutput(time)
    implicit none
    integer, intent(in)        :: time
    integer                    :: step, ierr

    call initializeOutput()

    step = output_index
    #ifndef HDF5
      call writeParticles_Binary(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeParticles_Binary()", .true.)
      call writeFields_Binary(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeFields_Binary()", .true.)
      call writeSpectra_Binary(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeSpectra_Binary()", .true.)
    #else
      call writeParticles_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeParticles_hdf5()", .true.)
      call writeFields_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeFields_hdf5()", .true.)
      call writeSpectra_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeSpectra_hdf5()", .true.)
      call writeDomain_hdf5(step, time)
        call printDiag((mpi_rank .eq. 0), "...writeDomain_hdf5()", .true.)
    #endif
    call printDiag((mpi_rank .eq. 0), "output()", .true.)
    output_index = output_index + 1
  end subroutine writeOutput

  subroutine initializeOutput()
    ! DEP_PRT [particle-dependent]
    implicit none
    real                      :: energy, u_, v_, w_
    integer                   :: s, i, ti, tj, tk, p, spec_index
    integer                   :: ierr
    integer, allocatable, dimension(:,:)  :: spectra
    integer, allocatable, dimension(:)    :: send_spec, recv_spec
    ! initialize particle variables
    n_prtl_vars = 8
    prtl_vars(1:n_prtl_vars) = (/'x    ', 'y    ', 'z    ', &
                               & 'u    ', 'v    ', 'w    ', &
                               & 'ind  ', 'proc '/)
    prtl_var_types(1:n_prtl_vars) = (/'real ', 'real ', 'real ', &
                                    & 'real ', 'real ', 'real ', &
                                    & 'int  ', 'int  '/)

    ! initialize field variables
    !   total number of fields (excluding particle densities)
    n_fld_vars = 12
    n_fld_vars = n_fld_vars + nspec
    do s = 1, nspec
      ! hopefully less than 10 species
      fld_vars(s) = 'dens' // STR(s) // ' '
    end do
    fld_vars(nspec + 1 : n_fld_vars) = (/'ex   ', 'ey   ', 'ez   ',&
                                       & 'bx   ', 'by   ', 'bz   ',&
                                       & 'jx   ', 'jy   ', 'jz   ',&
                                       & 'xx   ', 'yy   ', 'zz   '/)

    ! initialize domain output variables
    !   FIX1: maybe add # of particles per domain
    n_dom_vars = 6
    dom_vars(1 : 6) = (/'x0   ', 'y0   ', 'z0   ',&
                      & 'sx   ', 'sy   ', 'sz   '/)

    ! compute spectra
    if (.not. allocated(glob_spectra)) then
      allocate(glob_spectra(nspec, spec_num))
    end if
    allocate(spectra(nspec, spec_num))
    allocate(send_spec(spec_num), recv_spec(spec_num))

    spectra(:,:) = 0
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
             energy = log(energy)
             if (energy .le. spec_min) then
               spec_index = 1
             else if (energy .ge. spec_max) then
               spec_index = spec_num
             else
               spec_index = INT(CEILING((energy - spec_min) * REAL(spec_num) / (spec_max - spec_min)))
               if (spec_index .lt. 1) spec_index = 1
               if (spec_index .gt. spec_num) spec_index = spec_num
             end if
             spectra(s, spec_index) = spectra(s, spec_index) + 1
           end do
         end do
       end do
     end do
    end do

    ! send to root rank
    do s = 1, nspec
     send_spec(:) = spectra(s,:)
     call MPI_REDUCE(send_spec, recv_spec, spec_num, MPI_INTEGER,&
                   & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     glob_spectra(s,:) = recv_spec(:)
    end do

    if (allocated(spectra)) deallocate(spectra)
    if (allocated(send_spec)) deallocate(send_spec)
    if (allocated(recv_spec)) deallocate(recv_spec)
  end subroutine initializeOutput

  #ifndef HDF5
    !--- PRTL.TOT.***** structure ----------------------------------!
    ! HEADER:                                                      _
    !   timestep......................[4 bytes]                     |
    !   # of cpus.....................[4 bytes]                     |
    !   # of species..................[4 bytes]                     |
    !   # of variables................[4 bytes]                     |- disp_header
    !   variable names................[#var * 5 bytes]              |
    !   variable types................[#var * 5 bytes]              |
    !   # of particles per species....                              |
    !   ....summed over all ranks.....[#spec * 4 bytes]            _|
    ! BODY:
    !   species = 1...............[#of parts of species=1 * # of variables * 4 bytes]
    !     var = 1.................[#of parts per species * 4 bytes]
    !       rank = 1..............[#of parts on rank=1 per species * 4 bytes]
    !         XXX.................[4 bytes]
    !         XXX.................[4 bytes]
    !         ....................
    !         XXX.................[4 bytes]
    !       rank = 2..............[#of parts on rank=1 per species * 4 bytes]
    !       ......................
    !       rank = N..............[#of parts on rank=N per species * 4 bytes]
    !     var = 2.................[#of parts per species * 4 bytes]
    !     ........................
    !     var = V.................[#of parts per species * 4 bytes]
    !   species = 2...............[#of parts of species=2 * # of variables * 4 bytes]
    !   ..........................
    !   species = S...............[#of parts of species=S * # of variables * 4 bytes]
    !   ..........................
    !...............................................................!
    subroutine writeParticles_Binary(step, time)
      implicit none
      integer, intent(in)                 :: step, time
      character(len=STR_MAX)              :: stepchar, filename
      type(MPI_FILE)                      :: prtl_out_file
      integer                             :: ierr, i, p, s, rnk, j, temp, ti, tj, tk, temp_int
      integer(kind=MPI_OFFSET_KIND)       :: disp, disp_header
      integer                             :: npart_stride(nspec), npart_stride_global(nspec, mpi_size)
      integer                             :: npart_cum(nspec), npart_all(nspec)
      integer, allocatable, dimension(:)  :: temp_int_arr, stride_indices_arr, stride_ti_arr, stride_tj_arr, stride_tk_arr
      real, allocatable, dimension(:)     :: temp_real_arr
      real                                :: temp_real1, temp_real2

      disp_header = 4 * 4 + 5 * n_prtl_vars + 5 * n_prtl_vars + 4 * nspec

      ! preparation
      ! number of strided particles per each species
      do s = 1, nspec
        npart_stride(s) = 0
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                if (modulo(species(s)%prtl_tile(ti, tj, tk)%ind(p), output_stride) .eq. 0) then
                  npart_stride(s) = npart_stride(s) + 1
                end if
              end do ! p
            end do ! tk
          end do ! tj
        end do ! ti
      end do ! s

      call MPI_ALLGATHER(npart_stride, nspec, MPI_INTEGER,&
                       & npart_stride_global, nspec, MPI_INTEGER,&
                       & MPI_COMM_WORLD, ierr)

      ! computing offsets
      do s = 1, nspec
        npart_cum(s) = 0
        npart_all(s) = 0
        do rnk = 0, mpi_size - 1
          if (rnk .lt. mpi_rank) then
            npart_cum(s) = npart_cum(s) + npart_stride_global(s, rnk + 1)
          end if
          npart_all(s) = npart_all(s) + npart_stride_global(s, rnk + 1)
          #ifdef DEBUG
            if (npart_cum(s) .lt. 0) then
              call throwError('ERROR: # of particles to output exceeds int*4')
            end if
            if (npart_all(s) .lt. 0) then
              call throwError('ERROR: # of particles to output exceeds int*4')
            end if
          #endif
        end do
      end do

      ! create/open file
      write(stepchar, "(i5.5)") step
      filename = trim(output_dir_name) // '/prtl.tot.' // trim(stepchar)

      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,&
                       & MPI_MODE_WRONLY + MPI_MODE_CREATE,&
                       & MPI_INFO_NULL, prtl_out_file, ierr)

      disp = 0
      call MPI_FILE_SET_VIEW(prtl_out_file, disp, MPI_INTEGER,&
                           & MPI_INTEGER, "native",&
                           & MPI_INFO_NULL, ierr)

      ! create header
      if (mpi_rank .eq. 0) then
        call MPI_FILE_WRITE(prtl_out_file, time, 1, MPI_INTEGER,&
                          & MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_WRITE(prtl_out_file, mpi_size, 1, MPI_INTEGER,&
                          & MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_WRITE(prtl_out_file, nspec, 1, MPI_INTEGER,&
                          & MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_WRITE(prtl_out_file, n_prtl_vars, 1, MPI_INTEGER,&
                          & MPI_STATUS_IGNORE, ierr)
        ! variables
        do i = 1, n_prtl_vars
          call MPI_FILE_WRITE(prtl_out_file, prtl_vars(i), 5, MPI_CHARACTER,&
                            & MPI_STATUS_IGNORE, ierr)
        end do
        ! types of variables
        do i = 1, n_prtl_vars
          call MPI_FILE_WRITE(prtl_out_file, prtl_var_types(i), 5, MPI_CHARACTER,&
                            & MPI_STATUS_IGNORE, ierr)
        end do
        do s = 1, nspec
          call MPI_FILE_WRITE(prtl_out_file, npart_all(s), 1, MPI_INTEGER,&
                            & MPI_STATUS_IGNORE, ierr)
        end do
      end if

      ! writing particle data
      !   npart_stride(s)           : # of particles to output per species (for MPI block)
      !   stride_indices_arr(s)     : indices of particles to output for each species (for MPI block)
      !   stride_ti_arr(s)          : tile_i of particles to output for each species (for MPI block)
      !   stride_tj_arr(s)          : tile_j of particles to output for each species (for MPI block)
      !   stride_tk_arr(s)          : tile_k of particles to output for each species (for MPI block)
      !   npart_all(s)              : # of particles to output per species (for all MPI blocks)
      !   npart_cum(s)              : # of particles to output per species (for all MPI blocks with `rnk < mpi_rank`)
      !   temp_int_arr(j)           : contains integer data to write (for a local rank)
      !   temp_real_arr(j)          : contains real data to write (for a local rank)
      do s = 1, nspec
        allocate(stride_indices_arr(npart_stride(s)))
        allocate(stride_ti_arr(npart_stride(s)))
        allocate(stride_tj_arr(npart_stride(s)))
        allocate(stride_tk_arr(npart_stride(s)))
        j = 1
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                if (modulo(species(s)%prtl_tile(ti, tj, tk)%ind(p), output_stride) .eq. 0) then
                  stride_indices_arr(j) = p
                  stride_ti_arr(j) = ti
                  stride_tj_arr(j) = tj
                  stride_tk_arr(j) = tk
                  j = j + 1
                end if
              end do ! particles
            end do ! tk
          end do ! tj
        end do ! ti
        do i = 1, n_prtl_vars
          disp = disp_header +&
                  & (s - 1) * npart_all(s) * n_prtl_vars * 4 +&
                  & (i - 1) * npart_all(s) * 4 +&
                  & npart_cum(s) * 4
          call MPI_FILE_SET_VIEW(prtl_out_file, disp, MPI_INTEGER,&
                               & MPI_INTEGER, "native",&
                               & MPI_INFO_NULL, ierr)
          if (trim(prtl_var_types(i)) .eq. 'int') then
            allocate(temp_int_arr(npart_stride(s)))
            select case (trim(prtl_vars(i))) ! select integer variable
              case('ind')
                do j = 1, npart_stride(s)
                  temp = stride_indices_arr(j)
                  ti = stride_ti_arr(j)
                  tj = stride_tj_arr(j)
                  tk = stride_tk_arr(j)
                  temp_int = species(s)%prtl_tile(ti, tj, tk)%ind(temp)
                  temp_int_arr(j) = temp_int
                end do
              case('proc')
                do j = 1, npart_stride(s)
                  temp = stride_indices_arr(j)
                  ti = stride_ti_arr(j)
                  tj = stride_tj_arr(j)
                  tk = stride_tk_arr(j)
                  temp_int = species(s)%prtl_tile(ti, tj, tk)%proc(temp)
                  temp_int_arr(j) = temp_int
                end do
              case default
                call throwError('ERROR: unrecognized `prtl_vars`: `'//trim(prtl_vars(i))//'`')
            end select
            ! writing to a file
            call MPI_FILE_WRITE(prtl_out_file, temp_int_arr, npart_stride(s), MPI_INTEGER,&
                              & MPI_STATUS_IGNORE, ierr)
            deallocate(temp_int_arr)
          else if (trim(prtl_var_types(i)) .eq. 'real') then
            allocate(temp_real_arr(npart_stride(s)))
            select case (trim(prtl_vars(i))) ! select integer variable
              case('x')
                do j = 1, npart_stride(s)
                  temp = stride_indices_arr(j)
                  ti = stride_ti_arr(j)
                  tj = stride_tj_arr(j)
                  tk = stride_tk_arr(j)
                  temp_int = species(s)%prtl_tile(ti, tj, tk)%xi(temp)
                  temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dx(temp)
                  temp_real_arr(j) = REAL(this_meshblock%ptr%x0 + temp_int) + temp_real1
                end do
              case('y')
                do j = 1, npart_stride(s)
                  temp = stride_indices_arr(j)
                  ti = stride_ti_arr(j)
                  tj = stride_tj_arr(j)
                  tk = stride_tk_arr(j)
                  temp_int = species(s)%prtl_tile(ti, tj, tk)%yi(temp)
                  temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dy(temp)
                  temp_real_arr(j) = REAL(this_meshblock%ptr%y0 + temp_int) + temp_real1
                end do
              case('z')
                do j = 1, npart_stride(s)
                  temp = stride_indices_arr(j)
                  ti = stride_ti_arr(j)
                  tj = stride_tj_arr(j)
                  tk = stride_tk_arr(j)
                  temp_int = species(s)%prtl_tile(ti, tj, tk)%zi(temp)
                  temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dz(temp)
                  temp_real_arr(j) = REAL(this_meshblock%ptr%z0 + temp_int) + temp_real1
                end do
              case('u')
                do j = 1, npart_stride(s)
                  temp = stride_indices_arr(j)
                  ti = stride_ti_arr(j)
                  tj = stride_tj_arr(j)
                  tk = stride_tk_arr(j)
                  temp_real1 = species(s)%prtl_tile(ti, tj, tk)%u(temp)
                  temp_real_arr(j) = REAL(temp_real1, 4)
                end do
              case('v')
                do j = 1, npart_stride(s)
                  temp = stride_indices_arr(j)
                  ti = stride_ti_arr(j)
                  tj = stride_tj_arr(j)
                  tk = stride_tk_arr(j)
                  temp_real1 = species(s)%prtl_tile(ti, tj, tk)%v(temp)
                  temp_real_arr(j) = REAL(temp_real1, 4)
                end do
              case('w')
                do j = 1, npart_stride(s)
                  temp = stride_indices_arr(j)
                  ti = stride_ti_arr(j)
                  tj = stride_tj_arr(j)
                  tk = stride_tk_arr(j)
                  temp_real1 = species(s)%prtl_tile(ti, tj, tk)%w(temp)
                  temp_real_arr(j) = REAL(temp_real1, 4)
                end do
              case default
                call throwError('ERROR: unrecognized `prtl_vars`: `'//trim(prtl_vars(i))//'`')
            end select
            ! writing to a file
            call MPI_FILE_WRITE(prtl_out_file, temp_real_arr, npart_stride(s), MPI_REAL,&
                              & MPI_STATUS_IGNORE, ierr)
            deallocate(temp_real_arr)
          else
            call throwError('ERROR: unrecognized `prtl_var_types`: `'//trim(prtl_var_types(i))//'`')
          end if
        end do
        deallocate(stride_indices_arr)
        deallocate(stride_ti_arr)
        deallocate(stride_tj_arr)
        deallocate(stride_tk_arr)
      end do ! species

      call MPI_FILE_CLOSE(prtl_out_file, ierr)
    end subroutine writeParticles_Binary

    !--- FLDS.TOT.***** structure ----------------------------------!
    ! HEADER:                                                      _
    !   timestep......................[4 bytes]                     |
    !   # of cpus.....................[4 bytes]                     |
    !   # of fields...................[4 bytes]                     |
    !   dimensions..[fx,fy,fz]........[3 * 4 bytes]                 |- disp_header
    !   field names...................[#flds * 5 bytes]             |
    !   meshblock dimensions..........[# of cpus * 6 * 4 bytes]    _|
    ! BODY:
    !   field = 1.................[fx * fy * fz * 4 bytes]
    !     rank = 1................[fx * fy * fz (for rank = 1) * 4 bytes]
    !       XXX...................[4 bytes]
    !       XXX...................[4 bytes]
    !       ......................
    !       XXX...................[4 bytes]
    !     rank = 2................[fx * fy * fz (for rank = 2) * 4 bytes]
    !     ........................
    !     rank = N................[fx * fy * fz (for rank = N) * 4 bytes]
    !   field = 1.................[fx * fy * fz * 4 bytes]
    !   ..........................
    !   field = F.................[fx * fy * fz * 4 bytes]
    !   ..........................
    !...............................................................!
    subroutine writeFields_Binary(step, time)
      implicit none
      integer, intent(in)                 :: step, time
      character(len=STR_MAX)              :: stepchar, filename
      type(MPI_FILE)                      :: flds_out_file
      integer                             :: s, ierr, f_xyz, f, rnk, temp
      integer(kind=2)                     :: i, j, k
      integer(kind=MPI_OFFSET_KIND)       :: disp, disp_grid, disp_header
      integer                             :: nfld_cum, nfld_all
      real                                :: ex0, ey0, ez0, bx0, by0, bz0
      real                                :: jx0, jy0, jz0
      real, allocatable, dimension(:)     :: temp_real_arr

      ! FIX implement `output_istep` downsampling

      ! create/open file
      write(stepchar, "(i5.5)") step
      filename = trim(output_dir_name) // '/flds.tot.' // trim(stepchar)

      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,&
                      & MPI_MODE_WRONLY + MPI_MODE_CREATE,&
                      & MPI_INFO_NULL, flds_out_file, ierr)

      disp = 0
      call MPI_FILE_SET_VIEW(flds_out_file, disp, MPI_INTEGER,&
                          & MPI_INTEGER, "native",&
                          & MPI_INFO_NULL, ierr)

      ! create header
      if (mpi_rank .eq. 0) then
        call MPI_FILE_WRITE(flds_out_file, time, 1, MPI_INTEGER,&
                      & MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_WRITE(flds_out_file, mpi_size, 1, MPI_INTEGER,&
                      & MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_WRITE(flds_out_file, n_fld_vars, 1, MPI_INTEGER,&
                      & MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_WRITE(flds_out_file, global_mesh%sx, 1, MPI_INTEGER,&
                      & MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_WRITE(flds_out_file, global_mesh%sy, 1, MPI_INTEGER,&
                      & MPI_STATUS_IGNORE, ierr)
        call MPI_FILE_WRITE(flds_out_file, global_mesh%sz, 1, MPI_INTEGER,&
                      & MPI_STATUS_IGNORE, ierr)
        ! variables
        do f = 1, n_fld_vars
          call MPI_FILE_WRITE(flds_out_file, fld_vars(f), 5, MPI_CHARACTER,&
                        & MPI_STATUS_IGNORE, ierr)
        end do
        ! writing grid data
        do rnk = 0, mpi_size - 1
          call MPI_FILE_WRITE(flds_out_file, meshblocks(rnk + 1)%x0, 1, MPI_INTEGER,&
                            & MPI_STATUS_IGNORE, ierr)
          call MPI_FILE_WRITE(flds_out_file, meshblocks(rnk + 1)%y0, 1, MPI_INTEGER,&
                            & MPI_STATUS_IGNORE, ierr)
          call MPI_FILE_WRITE(flds_out_file, meshblocks(rnk + 1)%z0, 1, MPI_INTEGER,&
                            & MPI_STATUS_IGNORE, ierr)
          call MPI_FILE_WRITE(flds_out_file, meshblocks(rnk + 1)%sx, 1, MPI_INTEGER,&
                            & MPI_STATUS_IGNORE, ierr)
          call MPI_FILE_WRITE(flds_out_file, meshblocks(rnk + 1)%sy, 1, MPI_INTEGER,&
                            & MPI_STATUS_IGNORE, ierr)
          call MPI_FILE_WRITE(flds_out_file, meshblocks(rnk + 1)%sz, 1, MPI_INTEGER,&
                            & MPI_STATUS_IGNORE, ierr)
        end do
      end if

      disp_header = 4 * 3 + 4 * 3 + 5 * n_fld_vars + 6 * 4 * mpi_size

      ! writing field data
      !   nfld_cum              : cumulative # of grid points before `mpi_rank` (for all `rnk < mpi_rank`)
      !   nfld_all              : overall # of grid points for all ranks
      f_xyz = this_meshblock%ptr%sx * this_meshblock%ptr%sy * this_meshblock%ptr%sz
      allocate(temp_real_arr(f_xyz))
      ! computing displacements
      nfld_all = 0; nfld_cum = 0;
      do rnk = 0, mpi_size - 1
        nfld_all = nfld_all + meshblocks(rnk + 1)%sx * meshblocks(rnk + 1)%sy * meshblocks(rnk + 1)%sz
        if (rnk .lt. mpi_rank) then
          nfld_cum = nfld_cum + meshblocks(rnk + 1)%sx * meshblocks(rnk + 1)%sy * meshblocks(rnk + 1)%sz
        end if
      end do

      do f = 1, n_fld_vars
        if (fld_vars(f)(1:4) .eq. 'dens') then
          s = STRtoINT(fld_vars(f)(5:5))
          call computeDensity(s)
        end if
        disp = disp_header +&
                & (f - 1) * nfld_all * 4 +&
                & nfld_cum * 4
        call MPI_FILE_SET_VIEW(flds_out_file, disp, MPI_INTEGER,&
                                & MPI_INTEGER, "native",&
                                & MPI_INFO_NULL, ierr)
        temp = 1
        do i = 0, this_meshblock%ptr%sx - 1
          do j = 0, this_meshblock%ptr%sy - 1
            do k = 0, this_meshblock%ptr%sz - 1
              select case (trim(fld_vars(f)))
              case('ex')
                call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
                temp_real_arr(temp) = ex0 * B_norm
              case('ey')
                call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
                temp_real_arr(temp) = ey0 * B_norm
              case('ez')
                call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
                temp_real_arr(temp) = ez0 * B_norm
              case('bx')
                call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
                temp_real_arr(temp) = bx0 * B_norm
              case('by')
                call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
                temp_real_arr(temp) = by0 * B_norm
              case('bz')
                call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
                temp_real_arr(temp) = bz0 * B_norm
              case('jx')
                call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
                temp_real_arr(temp) = -jx0 * B_norm
              case('jy')
                call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
                temp_real_arr(temp) = -jy0 * B_norm
              case('jz')
                call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
                temp_real_arr(temp) = -jz0 * B_norm
              case default
                if (fld_vars(f)(1:4) .eq. 'dens') then
                  s = STRtoINT(fld_vars(f)(5:5))
                  temp_real_arr(temp) = REAL(scalar_int_array(i, j, k))
                end if
              end select
              temp = temp + 1
            end do
          end do
        end do

        call MPI_FILE_WRITE(flds_out_file, temp_real_arr, f_xyz, MPI_REAL,&
                          & MPI_STATUS_IGNORE, ierr)
      end do

      deallocate(temp_real_arr)
      call MPI_FILE_CLOSE(flds_out_file, ierr)
    end subroutine writeFields_Binary

    !--- SPEC.TOT.***** structure ----------------------------------!
    ! HEADER:
    !   timestep......................[4 bytes]
    !   # of species..................[4 bytes]
    !   # of bins.....................[4 bytes]
    !   MIN energy....................[4 bytes]
    !   MAX energy....................[4 bytes]
    ! BODY:
    !   species = 1...............[# of bins * 4 bytes]
    !     bin = 1.................[4 bytes]
    !     bin = 2.................[4 bytes]
    !     ........................
    !     rank = N................[4 bytes]
    !   species = 2...............[# of bins * 4 bytes]
    !   ..........................
    !   species = S...............[# of bins * 4 bytes]
    !   ..........................
    !...............................................................!
    subroutine writeSpectra_Binary(step, time)
      implicit none
      integer, intent(in)       :: step, time
      character(len=STR_MAX)    :: stepchar, filename
      integer                   :: ierr

      ! root rank writing to a file
      if (mpi_rank .eq. 0) then
        write(stepchar, "(i5.5)") step
        filename = trim(output_dir_name) // '/spec.tot.' // trim(stepchar)

        open(unit = UNIT_output, file = filename, form = "unformatted", access = "stream")
        write(unit = UNIT_output) time
        write(unit = UNIT_output) nspec
        write(unit = UNIT_output) spec_num
        write(unit = UNIT_output) spec_min
        write(unit = UNIT_output) spec_max
        do s = 1, nspec
          do i = 1, spec_num
            write(unit = UNIT_output) glob_spectra(s, i)
          end do
        end do

        close(unit = UNIT_output)
      end if

      if (allocated(glob_spectra)) deallocate(glob_spectra)
    end subroutine writeSpectra_Binary

  #else

  subroutine writeFields_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer(HID_T)                    :: file_id, dset_id(40), filespace(40), memspace, plist_id
    integer                           :: comm, info, error, f, s
    integer(kind=2)                   :: i, j, k
    logical                           :: writing_intQ, writing_densQ
    integer                           :: dataset_rank = 3
    integer(HSSIZE_T), dimension(3)   :: offsets
    integer(HSIZE_T), dimension(3)    :: global_dims, blocks
    real                              :: ex0, ey0, ez0, bx0, by0, bz0
    real                              :: jx0, jy0, jz0
    ! type(MPI_Info)                    :: FILE_INFO_TEMPLATE

    ! downsampling variables
    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    integer :: offset_i, offset_j, offset_k, i1, j1, k1
    integer :: n_i, n_j, n_k, glob_n_i, glob_n_j, glob_n_k

    ! for convenience
    this_x0 = this_meshblock%ptr%x0
    this_y0 = this_meshblock%ptr%y0
    this_z0 = this_meshblock%ptr%z0
    this_sx = this_meshblock%ptr%sx
    this_sy = this_meshblock%ptr%sy
    this_sz = this_meshblock%ptr%sz

    write(stepchar, "(i5.5)") step
    filename = trim(output_dir_name) // '/flds.tot.' // trim(stepchar)

    ! assuming `global_mesh%{x0,y0,z0} .eq. 0`
    if (output_istep .eq. 1) then
      offset_i = this_x0;   offset_j = this_y0;   offset_k = this_z0
      n_i = this_sx - 1;    n_j = this_sy - 1;    n_k = this_sz - 1
      glob_n_i = global_mesh%sx
      glob_n_j = global_mesh%sy
      glob_n_k = global_mesh%sz

      i_start = 0; j_start = 0; k_start = 0
    else
      offset_i = FLOOR(REAL(this_x0) / REAL(output_istep))
      offset_j = FLOOR(REAL(this_y0) / REAL(output_istep))

      i_start = CEILING(REAL(this_x0) / REAL(output_istep)) * output_istep - this_x0
      i_end = (CEILING(REAL(this_x0 + this_sx) / REAL(output_istep)) - 1) * output_istep - this_x0
      j_start = CEILING(REAL(this_y0) / REAL(output_istep)) * output_istep - this_y0
      j_end = (CEILING(REAL(this_y0 + this_sy) / REAL(output_istep)) - 1) * output_istep - this_y0

      n_i = (i_end - i_start) / output_istep
      n_j = (j_end - j_start) / output_istep

      glob_n_i = FLOOR(REAL(global_mesh%sx) / REAL(output_istep))
      glob_n_i = MAX(1, glob_n_i)

      glob_n_j = FLOOR(REAL(global_mesh%sy) / REAL(output_istep))
      glob_n_j = MAX(1, glob_n_j)

      #ifndef threeD
        k_start = 0; k_end = 0
        offset_k = 0; n_k = 0
        glob_n_k = 1
      #else
        offset_k = FLOOR(REAL(this_z0) / REAL(output_istep))
        k_start = CEILING(REAL(this_z0) / REAL(output_istep)) * output_istep - this_z0
        k_end = (CEILING(REAL(this_z0 + this_sz) / REAL(output_istep)) - 1) * output_istep - this_z0
        n_k = (k_end - k_start) / output_istep
        glob_n_k = FLOOR(REAL(global_mesh%sz) / REAL(output_istep))
        glob_n_k = MAX(1, glob_n_k)
      #endif
    end if

    offsets(1) = offset_i
    offsets(2) = offset_j
    offsets(3) = offset_k
    blocks(1) = n_i + 1
    blocks(2) = n_j + 1
    blocks(3) = n_k + 1
    global_dims(1) = glob_n_i
    global_dims(2) = glob_n_j
    global_dims(3) = glob_n_k

    ! mpi_f08 thing
    ! call MPI_INFO_CREATE(FILE_INFO_TEMPLATE, error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "access_style", "write_once", error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "collective_buffering", "true", error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "cb_block_size", "4194304", error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "cb_buffer_size", "16777216", error)
    ! call MPI_INFO_SET(FILE_INFO_TEMPLATE, "cb_nodes", "1", error)

    comm = MPI_COMM_WORLD%MPI_VAL
    info = MPI_INFO_NULL%MPI_VAL
    ! info = FILE_INFO_TEMPLATE%MPI_VAL

    call h5open_f(error)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    do f = 1, n_fld_vars
      call h5screate_simple_f(dataset_rank, global_dims, filespace(f), error)
    end do

    do f = 1, n_fld_vars
      if (fld_vars(f)(1:4) .eq. 'dens') then
        writing_densQ = .true.
        s = STRtoINT(fld_vars(f)(5:5))
        call computeDensity(s) ! filled `scalar_int_array` with density of species `s`
      else
        writing_densQ = .false.
      end if

      call h5dcreate_f(file_id, fld_vars(f), H5T_NATIVE_REAL, filespace(f), &
                     & dset_id(f), error)
      call h5sclose_f(filespace(f), error)
      call h5dget_space_f(dset_id(f), filespace(f), error)

      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      call h5screate_simple_f(dataset_rank, blocks, memspace, error)
      call h5sselect_hyperslab_f(filespace(f), H5S_SELECT_SET_F, offsets, blocks, error)

      ! Create dataset by interpolating fields
      do i1 = 0, n_i
        do j1 = 0, n_j
          do k1 = 0, n_k
            i = i_start + i1 * output_istep
            j = j_start + j1 * output_istep
            k = k_start + k1 * output_istep
            select case (trim(fld_vars(f)))
            case('ex')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
              scalar_real_array(i1, j1, k1) = ex0 * B_norm
            case('ey')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
              scalar_real_array(i1, j1, k1) = ey0 * B_norm
            case('ez')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
              scalar_real_array(i1, j1, k1) = ez0 * B_norm
            case('bx')
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
              scalar_real_array(i1, j1, k1) = bx0 * B_norm
            case('by')
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
              scalar_real_array(i1, j1, k1) = by0 * B_norm
            case('bz')
              call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
              scalar_real_array(i1, j1, k1) = bz0 * B_norm
            case('jx')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
              scalar_real_array(i1, j1, k1) = -jx0 * B_norm
            case('jy')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
              scalar_real_array(i1, j1, k1) = -jy0 * B_norm
            case('jz')
              call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
              scalar_real_array(i1, j1, k1) = -jz0 * B_norm
            case('xx')
              scalar_real_array(i1, j1, k1) = REAL(this_meshblock%ptr%x0 + i, 4)
            case('yy')
              scalar_real_array(i1, j1, k1) = REAL(this_meshblock%ptr%y0 + j, 4)
            case('zz')
              scalar_real_array(i1, j1, k1) = REAL(this_meshblock%ptr%z0 + k, 4)
            case default
              if ((fld_vars(f)(1:4) .ne. 'dens') .or. (.not. writing_densQ)) then
                call throwError("ERROR: unrecognized `fld_vars(f)`")
              else
                scalar_real_array(i1, j1, k1) = REAL(scalar_int_array(i, j, k), 4)
              end if
            end select
          end do
        end do
      end do

      ! Write the dataset collectively
      call h5dwrite_f(dset_id(f), H5T_NATIVE_REAL, scalar_real_array(0 : n_i, 0 : n_j, 0 : n_k), global_dims, error, &
                    & file_space_id = filespace(f), mem_space_id = memspace, xfer_prp = h5p_default_f)

      call h5dclose_f(dset_id(f), error)
      call h5sclose_f(filespace(f), error)
    end do

    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine writeFields_hdf5

  subroutine writeParticles_hdf5(step, time)
    ! DEP_PRT [particle-dependent]
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    character(len=7)                  :: dsetname
    integer(HID_T)                    :: file_id, dset_id(100), filespace(100), memspace, plist_id
    integer                           :: comm, info, error, ierr
    integer                           :: rnk, s, p, j, ln_, ti, tj, tk, temp, temp_int
    integer                           :: dataset_rank = 1
    integer(HID_T)                    :: h5type
    integer(HSSIZE_T), dimension(1)   :: offsets
    integer(HSIZE_T), dimension(1)    :: global_dims, blocks
    integer                           :: npart_stride(nspec), npart_stride_global(nspec, mpi_size)
    integer, allocatable, dimension(:):: temp_int_arr, stride_indices_arr, stride_ti_arr, stride_tj_arr, stride_tk_arr
    real, allocatable, dimension(:)   :: temp_real_arr
    real                              :: temp_real1, temp_real2
    logical                           :: writing_intQ

    ! preparation
    ! number of strided particles per each species
    do s = 1, nspec
      npart_stride(s) = 0
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (modulo(species(s)%prtl_tile(ti, tj, tk)%ind(p), output_stride) .eq. 0) then
                npart_stride(s) = npart_stride(s) + 1
              end if
            end do ! p
          end do ! tk
        end do ! tj
      end do ! ti
    end do ! s

    call MPI_ALLGATHER(npart_stride, nspec, MPI_INTEGER,&
                    & npart_stride_global, nspec, MPI_INTEGER,&
                    & MPI_COMM_WORLD, ierr)

    write(stepchar, "(i5.5)") step
    filename = trim(output_dir_name) // '/prtl.tot.' // trim(stepchar)

    ! mpi_f08 thing
    comm = MPI_COMM_WORLD%MPI_VAL
    info = MPI_INFO_NULL%MPI_VAL

    call h5open_f(error)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    do s = 1, nspec
      offsets(1) = 0
      global_dims(1) = 0
      do rnk = 0, mpi_size - 1
        if (rnk .lt. mpi_rank) then
          offsets(1) = offsets(1) + npart_stride_global(s, rnk + 1)
        end if
        global_dims(1) = global_dims(1) + npart_stride_global(s, rnk + 1)
      end do
      blocks(1) = npart_stride(s)

      ! saving the particle (and tile) indices to output for a given species
      allocate(stride_indices_arr(npart_stride(s)))
      allocate(stride_ti_arr(npart_stride(s)))
      allocate(stride_tj_arr(npart_stride(s)))
      allocate(stride_tk_arr(npart_stride(s)))
      j = 1
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
              if (modulo(species(s)%prtl_tile(ti, tj, tk)%ind(p), output_stride) .eq. 0) then
                stride_indices_arr(j) = p
                stride_ti_arr(j) = ti
                stride_tj_arr(j) = tj
                stride_tk_arr(j) = tk
                j = j + 1
              end if
            end do ! particles
          end do ! tk
        end do ! tj
      end do ! ti

      do p = 1, n_prtl_vars
        call h5screate_simple_f(dataset_rank, global_dims, filespace(p), error)
      end do

      do p = 1, n_prtl_vars
        ! dataset name `var_name` + `species #`
        ln_ = len(trim(prtl_vars(p)))
        dsetname(1 : ln_ + 1) = trim(prtl_vars(p)) // trim(STR(s))
        dsetname(ln_ + 2 : 7) = ' '
        ! creating dataset for a given type
        if (trim(prtl_var_types(p)) .eq. 'int') then
          writing_intQ = .true.
          h5type = H5T_NATIVE_INTEGER
          ! Create dataset
          allocate(temp_int_arr(npart_stride(s)))
          select case (trim(prtl_vars(p))) ! select integer variable
            case('ind')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                ti = stride_ti_arr(j)
                tj = stride_tj_arr(j)
                tk = stride_tk_arr(j)
                temp_int = species(s)%prtl_tile(ti, tj, tk)%ind(temp)
                temp_int_arr(j) = temp_int
              end do
            case('proc')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                ti = stride_ti_arr(j)
                tj = stride_tj_arr(j)
                tk = stride_tk_arr(j)
                temp_int = species(s)%prtl_tile(ti, tj, tk)%proc(temp)
                temp_int_arr(j) = temp_int
              end do
            case default
              call throwError('ERROR: unrecognized `prtl_vars`: `'//trim(prtl_vars(p))//'`')
          end select
        else if (trim(prtl_var_types(p)) .eq. 'real') then
          writing_intQ = .false.
          h5type = H5T_NATIVE_REAL
          ! Create dataset
          allocate(temp_real_arr(npart_stride(s)))
          select case (trim(prtl_vars(p))) ! select real variable
            case('x')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                ti = stride_ti_arr(j)
                tj = stride_tj_arr(j)
                tk = stride_tk_arr(j)
                temp_int = species(s)%prtl_tile(ti, tj, tk)%xi(temp)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dx(temp)
                temp_real_arr(j) = REAL(this_meshblock%ptr%x0 + temp_int) + temp_real1
              end do
            case('y')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                ti = stride_ti_arr(j)
                tj = stride_tj_arr(j)
                tk = stride_tk_arr(j)
                temp_int = species(s)%prtl_tile(ti, tj, tk)%yi(temp)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dy(temp)
                temp_real_arr(j) = REAL(this_meshblock%ptr%y0 + temp_int) + temp_real1
              end do
            case('z')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                ti = stride_ti_arr(j)
                tj = stride_tj_arr(j)
                tk = stride_tk_arr(j)
                temp_int = species(s)%prtl_tile(ti, tj, tk)%zi(temp)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%dz(temp)
                temp_real_arr(j) = REAL(this_meshblock%ptr%z0 + temp_int) + temp_real1
              end do
            case('u')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                ti = stride_ti_arr(j)
                tj = stride_tj_arr(j)
                tk = stride_tk_arr(j)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%u(temp)
                temp_real_arr(j) = REAL(temp_real1, 4)
              end do
            case('v')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                ti = stride_ti_arr(j)
                tj = stride_tj_arr(j)
                tk = stride_tk_arr(j)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%v(temp)
                temp_real_arr(j) = REAL(temp_real1, 4)
              end do
            case('w')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                ti = stride_ti_arr(j)
                tj = stride_tj_arr(j)
                tk = stride_tk_arr(j)
                temp_real1 = species(s)%prtl_tile(ti, tj, tk)%w(temp)
                temp_real_arr(j) = REAL(temp_real1, 4)
              end do
            case default
              call throwError('ERROR: unrecognized `prtl_vars`: `'//trim(prtl_vars(p))//'`')
          end select
        else
          call throwError('ERROR: unrecognized `prtl_var_types`: `'//trim(prtl_var_types(p))//'`')
        end if

        call h5dcreate_f(file_id, dsetname, h5type, filespace(p),&
                       & dset_id(p), error)
        call h5sclose_f(filespace(p), error)
        call h5dget_space_f(dset_id(p), filespace(p), error)

        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        call h5screate_simple_f(dataset_rank, blocks, memspace, error)
        call h5sselect_hyperslab_f(filespace(p), H5S_SELECT_SET_F, offsets, blocks, error)

        ! Write the dataset collectively
        if (writing_intQ) then
          call h5dwrite_f(dset_id(p), h5type, temp_int_arr, global_dims, error,&
                        & file_space_id = filespace(p), mem_space_id = memspace,&
                        & xfer_prp = h5p_default_f)
          deallocate(temp_int_arr)
        else
          call h5dwrite_f(dset_id(p), h5type, temp_real_arr, global_dims, error,&
                        & file_space_id = filespace(p), mem_space_id = memspace,&
                        & xfer_prp = h5p_default_f)
          deallocate(temp_real_arr)
        end if

        call h5sclose_f(filespace(p), error)
        call h5dclose_f(dset_id(p), error)
      end do
      deallocate(stride_indices_arr)
      deallocate(stride_ti_arr)
      deallocate(stride_tj_arr)
      deallocate(stride_tk_arr)
    end do

    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine writeParticles_hdf5

  subroutine writeSpectra_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: error, s, i, datarank
    integer(HID_T)                    :: file_id, dset_id, dspace_id
    integer(HSIZE_T), dimension(1)    :: data_dims
    character(len=2)                  :: dsetname
    real, allocatable, dimension(:)   :: bin_data

    datarank = 1
    data_dims(1) = spec_num

    ! only root rank writes the spectra file
    if (mpi_rank .eq. 0) then
      ! saving the energy bins
      allocate(bin_data(spec_num))
      do i = 1, spec_num
        bin_data(i) = spec_min + (REAL(i - 1, 4) / REAL(spec_num, 4)) * (spec_max - spec_min)
      end do

      write(stepchar, "(i5.5)") step
      filename = trim(output_dir_name) // '/spec.tot.' // trim(stepchar)

      ! Initialize FORTRAN interface
      call h5open_f(error)
      ! Create a new file using default properties
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      do s = 1, nspec
        ! writing bins:
        dsetname = 'e' // trim(STR(s))
        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                       & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL, bin_data, data_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        ! writing spectra:
        dsetname = 'n' // trim(STR(s))
        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, &
                       & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, glob_spectra(s,:), data_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
      end do

      ! Close the file
      call h5fclose_f(file_id, error)
      ! Close FORTRAN interface
      call h5close_f(error)

      if (allocated(bin_data)) deallocate(bin_data)
    end if

    if (allocated(glob_spectra)) deallocate(glob_spectra)
  end subroutine writeSpectra_hdf5

  subroutine writeDomain_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: error, s, i, datarank, d, rnk
    integer(HID_T)                    :: file_id, dset_id, dspace_id
    integer(HSIZE_T), dimension(1)    :: data_dims
    integer, allocatable, dimension(:):: domain_data

    datarank = 1
    data_dims(1) = mpi_size

    ! only root rank writes the domain file
    if (mpi_rank .eq. 0) then
      write(stepchar, "(i5.5)") step
      filename = trim(output_dir_name) // '/domain.' // trim(stepchar)

      ! Initialize FORTRAN interface
      call h5open_f(error)
      ! Create a new file using default properties
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      allocate(domain_data(mpi_size))

      do d = 1, n_dom_vars
        select case (trim(dom_vars(d)))
          case('x0')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%x0
            end do
          case('y0')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%y0
            end do
          case('z0')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%z0
            end do
          case('sx')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%sx
            end do
          case('sy')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%sy
            end do
          case('sz')
            do rnk = 0, mpi_size - 1
              domain_data(rnk + 1) = meshblocks(rnk + 1)%sz
            end do
          case default
            call throwError('ERROR: unrecognized `dom_vars`: `'//trim(dom_vars(d))//'`')
        end select

        call h5screate_simple_f(datarank, data_dims, dspace_id, error)
        call h5dcreate_f(file_id, dom_vars(d), H5T_NATIVE_INTEGER, dspace_id, &
                       & dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, domain_data, data_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)
      end do

      ! Close the file
      call h5fclose_f(file_id, error)
      ! Close FORTRAN interface
      call h5close_f(error)

      if (allocated(domain_data)) deallocate(domain_data)
    end if
  end subroutine writeDomain_hdf5

  #endif

end module m_writeoutput
