#include "../defs.F90"

module m_writeoutput
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  implicit none

  integer               :: output_stride, output_interval, output_istep

  !--- PRIVATE functions -----------------------------------------!
  private :: writeParticles, writeFields,&
           & computeDensity
  !...............................................................!
contains
  subroutine writeOutput(step, time)
    implicit none
    integer, intent(in)        :: step, time
    integer                    :: ierr
    call writeParticles(step, time)
    call writeFields(step, time)
    call printDiag((mpi_rank .eq. 0), TAB // "output()" // TAB // TAB // TAB // "[OK]")
  end subroutine writeOutput

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
  subroutine writeParticles(step, time)
    implicit none
    integer, intent(in)                 :: step, time
    character(len=STR_MAX)              :: stepchar, filename
    type(MPI_FILE)                      :: prtl_out_file
    integer                             :: ierr, i, p, s, rnk, nvars, j, temp, ti, tj, tk, temp_int
    integer(kind=MPI_OFFSET_KIND)       :: disp, disp_header
    character(len=STR_MAX)              :: vars(100), var_types(100)
    integer                             :: npart_stride(nspec), npart_stride_global(nspec, mpi_size)
    integer                             :: npart_cum(nspec), npart_all(nspec)

    integer, allocatable, dimension(:)  :: temp_int_arr, stride_indices_arr, stride_ti_arr, stride_tj_arr, stride_tk_arr
    real, allocatable, dimension(:)     :: temp_real_arr
		real                                :: temp_real1, temp_real2

    ! body
    nvars = 8
    vars(1:nvars) = (/'x    ', 'y    ', 'z    ', &
                    & 'u    ', 'v    ', 'w    ', &
                    & 'ind  ', 'proc '/)
    var_types(1:nvars) = (/'real ', 'real ', 'real ', &
                         & 'real ', 'real ', 'real ', &
                         & 'int  ', 'int  '/)

    disp_header = 4 * 4 + 5 * nvars + 5 * nvars + 4 * nspec

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
      do rnk = 1, mpi_rank
        npart_cum(s) = npart_cum(s) + npart_stride_global(s, rnk)
        #ifdef DEBUG
          if (npart_cum(s) .lt. 0) then
            call throwError('ERROR: # of particles to output exceeds int*4')
          end if
        #endif
      end do
      npart_all(s) = npart_cum(s)
      do rnk = mpi_rank + 1, mpi_size
        npart_all(s) = npart_all(s) + npart_stride_global(s, rnk)
        #ifdef DEBUG
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
      call MPI_FILE_WRITE(prtl_out_file, nvars, 1, MPI_INTEGER,&
                    & MPI_STATUS_IGNORE, ierr)
      ! variables
      do i = 1, nvars
        call MPI_FILE_WRITE(prtl_out_file, vars(i), 5, MPI_CHARACTER,&
                      & MPI_STATUS_IGNORE, ierr)
      end do
      ! types of variables
      do i = 1, nvars
        call MPI_FILE_WRITE(prtl_out_file, var_types(i), 5, MPI_CHARACTER,&
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
      do i = 1, nvars
        disp = disp_header +&
                & (s - 1) * npart_all(s) * nvars * 4 +&
                & (i - 1) * npart_all(s) * 4 +&
                & npart_cum(s) * 4
        call MPI_FILE_SET_VIEW(prtl_out_file, disp, MPI_INTEGER,&
                                & MPI_INTEGER, "native",&
                                & MPI_INFO_NULL, ierr)
        if (trim(var_types(i)) .eq. 'int') then
          allocate(temp_int_arr(npart_stride(s)))
          select case (trim(vars(i))) ! select integer variable
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
              call throwError('ERROR: unrecognized vars: `'//trim(vars(i))//'`')
          end select
          ! writing to a file
          call MPI_FILE_WRITE(prtl_out_file, temp_int_arr, npart_stride(s), MPI_INTEGER,&
                            & MPI_STATUS_IGNORE, ierr)
          deallocate(temp_int_arr)
        else if (trim(var_types(i)) .eq. 'real') then
          allocate(temp_real_arr(npart_stride(s)))
          select case (trim(vars(i))) ! select integer variable
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
              call throwError('ERROR: unrecognized vars: `'//trim(vars(i))//'`')
          end select
          ! writing to a file
          call MPI_FILE_WRITE(prtl_out_file, temp_real_arr, npart_stride(s), MPI_REAL,&
                            & MPI_STATUS_IGNORE, ierr)
          deallocate(temp_real_arr)
        else
          call throwError('ERROR: unrecognized var_types: `'//trim(var_types(i))//'`')
        end if
      end do
      deallocate(stride_indices_arr)
			deallocate(stride_ti_arr)
			deallocate(stride_tj_arr)
			deallocate(stride_tk_arr)
    end do ! species

    call MPI_FILE_CLOSE(prtl_out_file, ierr)
  end subroutine writeParticles

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
  subroutine writeFields(step, time)
    implicit none
    integer, intent(in)                 :: step, time
    character(len=STR_MAX)              :: stepchar, filename
    type(MPI_FILE)                      :: flds_out_file
    integer                             :: s, ierr, f_xyz, f, rnk, nflds, temp
    integer(kind=2)                     :: i, j, k
    integer(kind=MPI_OFFSET_KIND)       :: disp, disp_grid, disp_header
    character(len=STR_MAX)              :: flds(100)
    integer                             :: nfld_cum, nfld_all
    real                                :: ex0, ey0, ez0, bx0, by0, bz0
    real, allocatable, dimension(:)     :: temp_real_arr

    ! FIX implement `output_istep` downsampling

    ! body
    nflds = 6 + nspec
    do s = 1, nspec
      ! hopefully less than 10 species
      flds(s) = 'dens' // STR(s) // ' '
    end do
    flds(nspec + 1 : nspec + 6) = (/'ex   ', 'ey   ', 'ez   ',&
                                  & 'bx   ', 'by   ', 'bz   '/)

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
      call MPI_FILE_WRITE(flds_out_file, nflds, 1, MPI_INTEGER,&
                    & MPI_STATUS_IGNORE, ierr)
      call MPI_FILE_WRITE(flds_out_file, global_mesh%sx, 1, MPI_INTEGER,&
                    & MPI_STATUS_IGNORE, ierr)
      call MPI_FILE_WRITE(flds_out_file, global_mesh%sy, 1, MPI_INTEGER,&
                    & MPI_STATUS_IGNORE, ierr)
      call MPI_FILE_WRITE(flds_out_file, global_mesh%sz, 1, MPI_INTEGER,&
                    & MPI_STATUS_IGNORE, ierr)
      ! variables
      do f = 1, nflds
        call MPI_FILE_WRITE(flds_out_file, flds(f), 5, MPI_CHARACTER,&
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

    disp_header = 4 * 3 + 4 * 3 + 5 * nflds + 6 * 4 * mpi_size

    ! writing field data
    !   nfld_cum              : cumulative # of grid points before `mpi_rank` (for all `rnk < mpi_rank`)
    !   nfld_all              : overall # of grid points for all ranks
    f_xyz = this_meshblock%ptr%sx * this_meshblock%ptr%sy * this_meshblock%ptr%sz
    allocate(temp_real_arr(f_xyz))
    ! computing displacements
    nfld_all = 0; nfld_cum = 0;
    do rnk = 0, mpi_size - 1
      nfld_all = nfld_all + meshblocks(rnk + 1)%sx * meshblocks(rnk + 1)%sy * meshblocks(rnk + 1)%sz
      if (rnk < mpi_rank) then
        nfld_cum = nfld_cum + meshblocks(rnk + 1)%sx * meshblocks(rnk + 1)%sy * meshblocks(rnk + 1)%sz
      end if
    end do

    do f = 1, nflds
      disp = disp_header +&
              & (f - 1) * nfld_all * 4 +&
              & nfld_cum * 4
      call MPI_FILE_SET_VIEW(flds_out_file, disp, MPI_INTEGER,&
                              & MPI_INTEGER, "native",&
                              & MPI_INFO_NULL, ierr)
      temp = 1
      if (flds(f)(1:4) .eq. 'dens') then
        s = STRtoINT(flds(f)(5:5))
        call computeDensity(s)
      end if
      do i = 0, this_meshblock%ptr%sx - 1
        do j = 0, this_meshblock%ptr%sy - 1
          do k = 0, this_meshblock%ptr%sz - 1
            call interpFlds(0.0, 0.0, 0.0, i, j, k, ex0, ey0, ez0, bx0, by0, bz0)
            select case (trim(flds(f)))
            case('ex')
              temp_real_arr(temp) = ex0 * B_norm
            case('ey')
              temp_real_arr(temp) = ey0 * B_norm
            case('ez')
              temp_real_arr(temp) = ez0 * B_norm
            case('bx')
              temp_real_arr(temp) = bx0 * B_norm
            case('by')
              temp_real_arr(temp) = by0 * B_norm
            case('bz')
              temp_real_arr(temp) = bz0 * B_norm
            case default
              if (flds(f)(1:4) .eq. 'dens') then
                temp_real_arr(temp) = scalar_array(i, j, k)
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
  end subroutine writeFields

  subroutine computeDensity(s)
    implicit none
    integer, intent(in)                   :: s
    integer                               :: p, ti, tj, tk
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
		scalar_array(:,:,:) = 0
		do ti = 1, species(s)%tile_nx
			do tj = 1, species(s)%tile_ny
				do tk = 1, species(s)%tile_nz
					pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
					pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
					pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi
					! FIX1 vectorize/align
					!$omp simd
					!dir$ vector aligned
					do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
						scalar_array(pt_xi(p), pt_yi(p), pt_zi(p)) = scalar_array(pt_xi(p), pt_yi(p), pt_zi(p)) + 1
					end do
					pt_xi => null(); pt_yi => null(); pt_zi => null()
				end do
			end do
		end do
  end subroutine computeDensity

end module m_writeoutput
