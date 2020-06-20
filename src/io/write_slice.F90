#include "../defs.F90"

module m_writeslice
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
  use m_exchangearray

  implicit none

  logical                           :: slice_enable
  integer                           :: slice_start, slice_interval
  integer, private                  :: n_fld_vars, n_prtl_vars, n_dom_vars
  character(len=STR_MAX), private   :: fld_vars(100), dom_vars(100)
  logical                           :: slice_xdmf

  !--- PRIVATE functions -----------------------------------------!
  #ifdef HDF5
    private :: writeSlices_hdf5, writeXDMF_hdf5
  #endif
  private :: initializeSliceOutput
  !...............................................................!

contains
  subroutine writeSlices(time)
    implicit none
    integer, intent(in)        :: time
    integer                    :: step, ierr

    call initializeSliceOutput()

    step = slice_index
    #ifdef HDF5
      call writeSlices_hdf5(step, time)
        call printReport((mpi_rank .eq. 0), "...writeSlices_hdf5()", .true.)
    #endif
    call printDiag((mpi_rank .eq. 0), "slices()", .true.)
    slice_index = slice_index + 1
  end subroutine writeSlices

  subroutine initializeSliceOutput()
    ! DEP_PRT [particle-dependent]
    implicit none
    integer                   :: s
    integer                   :: ierr
    ! initialize field variables
    !   total number of fields (excluding particle densities)
    n_fld_vars = 12
    n_fld_vars = n_fld_vars + 2 * nspec
    do s = 1, nspec
      ! hopefully less than 10 species
      fld_vars(s) = 'dens' // STR(s)
    end do
    do s = 1, nspec
      ! hopefully less than 10 species
      fld_vars(nspec + s) = 'enrg' // STR(s)
    end do
    fld_vars(2 * nspec + 1 : n_fld_vars) = (/'ex   ', 'ey   ', 'ez   ',&
                                           & 'bx   ', 'by   ', 'bz   ',&
                                           & 'jx   ', 'jy   ', 'jz   ',&
                                           & 'xx   ', 'yy   ', 'zz   '/)

    ! initialize domain output variables
    !   FIX1: maybe add # of particles per domain
    n_dom_vars = 6
    dom_vars(1 : 6) = (/'x0   ', 'y0   ', 'z0   ',&
                      & 'sx   ', 'sy   ', 'sz   '/)
  end subroutine initializeSliceOutput

  #ifdef HDF5
  subroutine writeXDMF_hdf5(step, time, ni, nj, nk)
    implicit none
    integer, intent(in)               :: step, time, ni, nj, nk
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: var

    ! write(stepchar, "(i5.5)") step
    ! filename = trim(output_dir_name) // '/flds.tot.' // trim(stepchar) // '.xdmf'
    !
    ! open (UNIT_xdmf, file=filename, status="replace", access="stream", form="formatted")
    ! write (UNIT_xdmf, "(A)")&
    !                        & '<?xml version="1.0" ?>'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '<Xdmf Version="2.0">'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '  <Domain>'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '    <Grid Name="domain" GridType="Uniform">'
    ! write (UNIT_xdmf, "(A,A,A,I10,I10,I10,A)")&
    !                        & '      <Topology TopologyType=',&
    !                           & '"3DCoRectMesh"', ' Dimensions="', &
    !                           & nk, nj, ni, '"/>'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
    ! write (UNIT_xdmf, "(A,A)")&
    !                        & '        <DataItem Format="XML" Dimensions="3"',&
    !                           & ' NumberType="Float" Precision="4">'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '          0.0 0.0 0.0'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '        </DataItem>'
    ! write (UNIT_xdmf, "(A,A)")&
    !                        & '        <DataItem Format="XML" Dimensions="3"',&
    !                           & ' NumberType="Float" Precision="4">'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '          1.0 1.0 1.0'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '        </DataItem>'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '      </Geometry>'
    !
    ! do var = 1, n_fld_vars
    !   write (UNIT_xdmf, "(A)")&
    !                        & '      <Attribute Name="' // trim(fld_vars(var)) // '" Center="Node">'
    !   write (UNIT_xdmf, "(A,I10,I10,I10,A)")&
    !                        & '        <DataItem Format="HDF" Dimensions="',&
    !                        & nk, nj, ni,&
    !                        & '" NumberType="Float" Precision="4">'
    !   write (UNIT_xdmf, "(A,A)")&
    !                        & '          flds.tot.' // trim(stepchar) // ':/',&
    !                        & trim(fld_vars(var))
    !   write (UNIT_xdmf, "(A)")&
    !                        & '        </DataItem>'
    !   write (UNIT_xdmf, "(A)")&
    !                        & '      </Attribute>'
    ! end do
    !
    ! write (UNIT_xdmf, "(A)")&
    !                        & '    </Grid>'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '  </Domain>'
    ! write (UNIT_xdmf, "(A)")&
    !                        & '</Xdmf>'
    ! close (UNIT_xdmf)
  end subroutine writeXDMF_hdf5

  subroutine writeSlices_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: error
    integer(HID_T)                    :: plist_id, file_id, dset_id
    integer(HID_T)                    :: filespace, memspace
    integer                           :: dset_rank = 3
    integer(HSIZE_T), dimension(3)    :: global_dims, chunk_dims
    integer(HSIZE_T), dimension(3)    :: count_h5
    integer(HSSIZE_T), dimension(3)   :: offset_h5
    integer(HSIZE_T), dimension(3)    :: stride_h5
    integer(HSIZE_T), dimension(3)    :: block_h5

    real, allocatable                 :: dummy_array(:,:,:)

    integer :: i, j, k
    integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz, ymid
    ! integer(HID_T)                    :: file_id, dset_id(40), filespace(40), memspace, plist_id
    ! integer                           :: error, f, s
    ! integer(kind=2)                   :: i, j, k
    ! logical                           :: writing_intQ, writing_lgarrQ

    ! integer(HSSIZE_T), dimension(3)   :: offsets
    ! integer(HSIZE_T), dimension(3)    :: global_dims, blocks
    ! real                              :: ex0, ey0, ez0, bx0, by0, bz0
    ! real                              :: jx0, jy0, jz0

    ! integer                           :: ymid
    !

    ! integer :: i_start, i_end, j_start, j_end, k_start, k_end
    ! integer :: offset_i, offset_j, offset_k, i1, j1, k1
    ! integer :: n_i, n_j, n_k, glob_n_i, glob_n_j, glob_n_k

    this_x0 = this_meshblock%ptr%x0
    this_y0 = this_meshblock%ptr%y0
    this_z0 = this_meshblock%ptr%z0
    this_sx = this_meshblock%ptr%sx
    this_sy = this_meshblock%ptr%sy
    this_sz = this_meshblock%ptr%sz

    ymid = INT(global_mesh%sy / 2)

    chunk_dims(1) = this_sx
    chunk_dims(2) = 1
    chunk_dims(3) = this_sz

    global_dims(1) = global_mesh%sx
    global_dims(2) = 1
    global_dims(3) = global_mesh%sz

    write(stepchar, "(i5.5)") step
    filename = trim(slice_dir_name) // '/slices.' // trim(stepchar)

    ! Initialize HDF5 library and Fortran interfaces
    call H5open_f(error)
    ! Setup file access property list with parallel I/O access
    call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call H5Pset_fapl_mpio_f(plist_id, h5comm, h5info, error)

    ! Create the file collectively
    call H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call H5Pclose_f(plist_id, error)

    ! Create the data space for the dataset
    call H5Screate_simple_f(dset_rank, global_dims, filespace, error)
    call H5Screate_simple_f(dset_rank, chunk_dims, memspace, error)

    ! Create chunked dataset
    call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    call H5Pset_chunk_f(plist_id, dset_rank, chunk_dims, error)
    call H5Dcreate_f(file_id, fld_vars(1), H5T_NATIVE_REAL, filespace, dset_id, error, plist_id)
    call H5Sclose_f(filespace, error)

    stride_h5(:) = 1
    if ((ymid - this_y0 .ge. 0) .and. (ymid - this_y0 .le. this_sy - 1)) then
      count_h5(:) = 1
    else
      count_h5(:) = 0
    end if
    block_h5(:) = chunk_dims(:)
    offset_h5(1) = this_x0
    offset_h5(2) = 0
    offset_h5(3) = this_z0

    ! Select hyperslab in the file
    call H5Dget_space_f(dset_id, filespace, error)
    call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_h5, count_h5, error, stride_h5, block_h5)

    allocate(dummy_array(chunk_dims(1), chunk_dims(2), chunk_dims(3)))

    dummy_array(:,:,:) = -10.0

    if ((ymid - this_y0 .ge. 0) .and. (ymid - this_y0 .le. this_sy - 1)) then
      j = ymid - this_y0
      do i = 0, this_sx - 1
        do k = 0, this_sz - 1
          dummy_array(i + 1, 1, k + 1) = bz(i, j, k)
        end do
      end do
    end if

    ! Create property list for collective dataset write
    call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    ! Write the dataset collectively
    call H5Dwrite_f(dset_id, H5T_NATIVE_REAL, dummy_array, global_dims, error,&
                  & file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    deallocate(dummy_array)

    ! Close dataspaces
    call H5Sclose_f(filespace, error)
    call H5Sclose_f(memspace, error)

    ! Close the dataset
    call H5Dclose_f(dset_id, error)

    ! Close the property list
    call H5Pclose_f(plist_id, error)

    ! Close the file
    call H5Fclose_f(file_id, error)

    ! Close FORTRAN interfaces and HDF5 library
    call H5close_f(error)
  end subroutine writeSlices_hdf5

  ! subroutine writeSlices_hdf5(step, time)
  !   implicit none
  !   integer, intent(in)               :: step, time
  !   character(len=STR_MAX)            :: stepchar, filename
  !   integer(HID_T)                    :: file_id, dset_id(40), filespace(40), memspace, plist_id
  !   integer                           :: error, f, s
  !   integer(kind=2)                   :: i, j, k
  !   logical                           :: writing_intQ, writing_lgarrQ
  !   integer                           :: dataset_rank = 3
  !   integer(HSSIZE_T), dimension(3)   :: offsets
  !   integer(HSIZE_T), dimension(3)    :: global_dims, blocks
  !   real                              :: ex0, ey0, ez0, bx0, by0, bz0
  !   real                              :: jx0, jy0, jz0
  !
  !   integer                           :: ymid
  !
  !   ! downsampling variables
  !   integer :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
  !   integer :: i_start, i_end, j_start, j_end, k_start, k_end
  !   integer :: offset_i, offset_j, offset_k, i1, j1, k1
  !   integer :: n_i, n_j, n_k, glob_n_i, glob_n_j, glob_n_k
  !
  !   ! for convenience
  !   this_x0 = this_meshblock%ptr%x0
  !   this_y0 = this_meshblock%ptr%y0
  !   this_z0 = this_meshblock%ptr%z0
  !   this_sx = this_meshblock%ptr%sx
  !   this_sy = this_meshblock%ptr%sy
  !   this_sz = this_meshblock%ptr%sz
  !
  !   ymid = INT(global_mesh%sy / 2)
  !
  !   write(stepchar, "(i5.5)") step
  !   filename = trim(slice_dir_name) // '/slices.' // trim(stepchar)
  !
  !   ! RIGHT NOW THIS ONLY WORKS for 3D case
  !   ! ... with a cut in XZ plane Y = middle
  !   ! ... with an istep = 1
  !
  !   ! assuming `global_mesh%{x0,y0,z0} .eq. 0`
  !   offset_i = this_x0;   offset_j = this_y0;   offset_k = this_z0
  !   glob_n_i = global_mesh%sx
  !   glob_n_j = global_mesh%sy
  !   glob_n_k = global_mesh%sz
  !
  !   i_start = 0; j_start = 0; k_start = 0
  !
  !   ! if ((mpi_rank .eq. 0) .and. slice_xdmf) then
  !   !   call writeXDMF_hdf5(step, time, glob_n_i, glob_n_j, glob_n_k)
  !   ! end if
  !   n_i = this_sx - 1;    n_j = this_sy - 1;    n_k = this_sz - 1
  !
  !   offsets(1) = offset_i
  !   offsets(2) = offset_j
  !   offsets(3) = offset_k
  !   if ((ymid .ge. this_y0) .or. (ymid .lt. this_y0 + this_sx - 1)) then
  !     blocks(1) = n_i + 1
  !     blocks(2) = 1
  !     blocks(3) = n_k + 1
  !   else
  !     blocks(1) = 0
  !     blocks(2) = 0
  !     blocks(3) = 0
  !   end if
  !   global_dims(1) = glob_n_i
  !   global_dims(2) = glob_n_j
  !   global_dims(3) = glob_n_k
  !
  !   call h5open_f(error)
  !   call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  !   call h5pset_fapl_mpio_f(plist_id, h5comm, h5info, error)
  !   call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  !   call h5pclose_f(plist_id, error)
  !
  !   do f = 1, n_fld_vars
  !     call h5screate_simple_f(dataset_rank, global_dims, filespace(f), error)
  !   end do
  !
  !   do f = 1, n_fld_vars
  !     if (fld_vars(f)(1:4) .eq. 'dens') then
  !       writing_lgarrQ = .true.
  !       s = STRtoINT(fld_vars(f)(5:5))
  !       call computeDensity(s, reset=.true.) ! filled `lg_arr` with density of species `s`
  !       call exchangeArray()
  !     else if (fld_vars(f)(1:4) .eq. 'enrg') then
  !       writing_lgarrQ = .true.
  !       s = STRtoINT(fld_vars(f)(5:5))
  !       call computeEnergy(s, reset=.true.) ! filled `lg_arr` with energies of species `s`
  !       call exchangeArray()
  !     else
  !       writing_lgarrQ = .false.
  !     end if
  !
  !     call h5dcreate_f(file_id, fld_vars(f), H5T_NATIVE_REAL, filespace(f), &
  !                    & dset_id(f), error)
  !     call h5sclose_f(filespace(f), error)
  !     call h5dget_space_f(dset_id(f), filespace(f), error)
  !
  !     call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  !     call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  !
  !     call h5screate_simple_f(dataset_rank, blocks, memspace, error)
  !     call h5sselect_hyperslab_f(filespace(f), H5S_SELECT_SET_F, offsets, blocks, error)
  !
  !     ! Create dataset by interpolating fields
  !     j1 = 0
  !     do i1 = 0, n_i
  !       do k1 = 0, n_k
  !         i = i_start + i1
  !         j = ymid - this_y0
  !         k = k_start + k1
  !         select case (trim(fld_vars(f)))
  !         case('ex')
  !           call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
  !           sm_arr(i1, j1, k1) = ex0 * B_norm
  !         case('ey')
  !           call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
  !           sm_arr(i1, j1, k1) = ey0 * B_norm
  !         case('ez')
  !           call interpFromEdges(0.0, 0.0, 0.0, i, j, k, ex, ey, ez, ex0, ey0, ez0)
  !           sm_arr(i1, j1, k1) = ez0 * B_norm
  !         case('bx')
  !           call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
  !           sm_arr(i1, j1, k1) = bx0 * B_norm
  !         case('by')
  !           call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
  !           sm_arr(i1, j1, k1) = by0 * B_norm
  !         case('bz')
  !           call interpFromFaces(0.0, 0.0, 0.0, i, j, k, bx, by, bz, bx0, by0, bz0)
  !           sm_arr(i1, j1, k1) = bz0 * B_norm
  !         case('jx')
  !           call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
  !           sm_arr(i1, j1, k1) = -jx0 * B_norm
  !         case('jy')
  !           call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
  !           sm_arr(i1, j1, k1) = -jy0 * B_norm
  !         case('jz')
  !           call interpFromEdges(0.0, 0.0, 0.0, i, j, k, jx, jy, jz, jx0, jy0, jz0)
  !           sm_arr(i1, j1, k1) = -jz0 * B_norm
  !         case('xx')
  !           sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%x0 + i, 4)
  !         case('yy')
  !           sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%y0 + j, 4)
  !         case('zz')
  !           sm_arr(i1, j1, k1) = REAL(this_meshblock%ptr%z0 + k, 4)
  !         case default
  !           if (((fld_vars(f)(1:4) .ne. 'dens') .and. (fld_vars(f)(1:4) .ne. 'enrg')) .or.&
  !              & (.not. writing_lgarrQ)) then
  !             call throwError("ERROR: unrecognized `fld_vars(f)`")
  !           else
  !             sm_arr(i1, j1, k1) = lg_arr(i, j, k)
  !           end if
  !         end select
  !       end do
  !     end do
  !
  !     ! Write the dataset collectively
  !     call h5dwrite_f(dset_id(f), H5T_NATIVE_REAL, sm_arr(0 : n_i, j1 : j1, 0 : n_k), global_dims, error, &
  !                   & file_space_id = filespace(f), mem_space_id = memspace, xfer_prp = plist_id)
  !
  !     call h5dclose_f(dset_id(f), error)
  !     call h5sclose_f(filespace(f), error)
  !   end do
  !
  !   call h5sclose_f(memspace, error)
  !   call h5pclose_f(plist_id, error)
  !   call h5fclose_f(file_id, error)
  !   call h5close_f(error)
  !
  ! end subroutine writeSlices_hdf5
  #endif

end module m_writeslice
