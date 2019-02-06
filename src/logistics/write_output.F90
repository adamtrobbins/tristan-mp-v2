#include "../defs.F90"

module m_writeoutput
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_communications
  use m_domain
  use m_particles
  #ifdef HDF5
    use HDF5
  #endif
  implicit none

  integer :: output_stride, output_interval

  !--- PRIVATE functions -----------------------------------------!
  private :: writeParticles
  !...............................................................!
contains
  subroutine writeOutput(step)
    implicit none
    integer, intent(in)        :: step
    integer                    :: ierr
    call writeParticles(step)
  end subroutine writeOutput

  !--- PRTL.TOT.***** structure ----------------------------------!
  ! HEADER:                                                      _
  !   # of cpus.....................[4 bytes]                     |
  !   # of species..................[4 bytes]                     |
  !   # of variables................[4 bytes]                     |
  !   variable names................[#var * 5 bytes]              |- disp_header
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
  subroutine writeParticles(step)
    implicit none
    integer, intent(in)                 :: step
    character(len=STR_MAX)              :: stepchar, filename
    integer                             :: prtl_out_file, ierr, i, p, s, rnk, nvars, j, temp
    integer(kind=MPI_OFFSET_KIND)       :: disp, disp_header
    character(len=STR_MAX)              :: vars(100), var_types(100)
    integer                             :: npart_stride(nspec), npart_stride_global(nspec, mpi_size)
    integer                             :: npart_cum(nspec), npart_all(nspec)

    integer, allocatable, dimension(:)  :: temp_int_arr, stride_indices_arr
    real, allocatable, dimension(:)     :: temp_real_arr

    ! body
    nvars = 8
    vars(1:nvars) = (/'x    ', 'y    ', 'z    ', &
                    & 'u    ', 'v    ', 'w    ', &
                    & 'ind  ', 'proc '/)
    var_types(1:nvars) = (/'real ', 'real ', 'real ', &
                         & 'real ', 'real ', 'real ', &
                         & 'int  ', 'int  '/)

    disp_header = 4 * 3 + 5 * nvars + 5 * nvars + 4 * nspec

    ! preparation
    ! number of strided particles per each species
    do s = 1, nspec
      npart_stride(s) = 0
      do p = 1, spp_(s)%npart_sp
        if (modulo(sp_(s)%ind(p), output_stride) .eq. 0) then
          npart_stride(s) = npart_stride(s) + 1
        end if
      end do
    end do

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
    !   npart_stride(s)           : # of particles to output per species (for a local rank)
    !   stride_indices_arr(s)     : indices of particles to output for each species (for local rank)
    !   npart_all(s)              : # of particles to output per species (for all ranks)
    !   npart_cum(s)              : # of particles to output per species (for all ranks with `rnk < mpi_rank`)
    !   temp_int_arr(j)           : contains integer data to write (for a local rank)
    !   temp_real_arr(j)          : contains real data to write (for a local rank)
    do s = 1, nspec
      if (npart_stride(s) .eq. 0) continue
      allocate(stride_indices_arr(npart_stride(s)))
      j = 1
      do p = 1, spp_(s)%npart_sp
        if (modulo(sp_(s)%ind(p), output_stride) .eq. 0) then
          stride_indices_arr(j) = p
          j = j + 1
        end if
      end do
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
                temp_int_arr(j) = sp_(s)%ind(stride_indices_arr(j))
              end do
            case('proc')
              do j = 1, npart_stride(s)
                temp_int_arr(j) = sp_(s)%proc(stride_indices_arr(j))
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
                temp_real_arr(j) = REAL(this_meshblock%ptr%x0 - 1 + sp_(s)%xi(temp)) + sp_(s)%dx(temp)
              end do
            case('y')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                temp_real_arr(j) = REAL(this_meshblock%ptr%y0 - 1 + sp_(s)%yi(temp)) + sp_(s)%dy(temp)
              end do
            case('z')
              do j = 1, npart_stride(s)
                temp = stride_indices_arr(j)
                temp_real_arr(j) = REAL(this_meshblock%ptr%z0 - 1 + sp_(s)%zi(temp)) + sp_(s)%dz(temp)
              end do
            case('u')
              do j = 1, npart_stride(s)
                temp_real_arr(j) = REAL(sp_(s)%u(stride_indices_arr(j)), 4)
              end do
            case('v')
              do j = 1, npart_stride(s)
                temp_real_arr(j) = REAL(sp_(s)%v(stride_indices_arr(j)), 4)
              end do
            case('w')
              do j = 1, npart_stride(s)
                temp_real_arr(j) = REAL(sp_(s)%w(stride_indices_arr(j)), 4)
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
    end do

    call MPI_FILE_CLOSE(prtl_out_file, ierr)
  end subroutine writeParticles

end module m_writeoutput
