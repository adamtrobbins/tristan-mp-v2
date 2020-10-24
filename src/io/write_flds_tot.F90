#include "../defs.F90"

module m_writefldstot
  #ifdef HDF5
    use hdf5
  #endif
  use m_globalnamespace, only: output_dir_name, mpi_rank, h5comm, h5info
  use m_aux
  use m_errors
  use m_domain
  use m_fields
  use m_outputlogistics, only: prepareFieldForOutput, selectFieldForOutput,&
                            & fld_vars, n_fld_vars, output_flds_istep
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  #ifdef HDF5
    private :: writeFields_hdf5, writeXDMF_hdf5
  #endif
  !...............................................................!
contains

  #ifdef HDF5

  subroutine writeXDMF_hdf5(step, time, ni, nj, nk)
    implicit none
    integer, intent(in)               :: step, time, ni, nj, nk
    character(len=STR_MAX)            :: stepchar, filename
    integer                           :: var

    write(stepchar, "(i5.5)") step
    filename = trim(output_dir_name) // '/flds.tot.' // trim(stepchar) // '.xdmf'

    open (UNIT_xdmf, file=filename, status="replace", access="stream", form="formatted")
    write (UNIT_xdmf, "(A)")&
                           & '<?xml version="1.0" ?>'
    write (UNIT_xdmf, "(A)")&
                           & '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
    write (UNIT_xdmf, "(A)")&
                           & '<Xdmf Version="2.0">'
    write (UNIT_xdmf, "(A)")&
                           & '  <Domain>'
    write (UNIT_xdmf, "(A)")&
                           & '    <Grid Name="domain" GridType="Uniform">'
    write (UNIT_xdmf, "(A,A,A,I10,I10,I10,A)")&
                           & '      <Topology TopologyType=',&
                              & '"3DCoRectMesh"', ' Dimensions="', &
                              & nk, nj, ni, '"/>'
    write (UNIT_xdmf, "(A)")&
                           & '      <Geometry GeometryType="ORIGIN_DXDYDZ">'
    write (UNIT_xdmf, "(A,A)")&
                           & '        <DataItem Format="XML" Dimensions="3"',&
                              & ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)")&
                           & '          0.0 0.0 0.0'
    write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
    write (UNIT_xdmf, "(A,A)")&
                           & '        <DataItem Format="XML" Dimensions="3"',&
                              & ' NumberType="Float" Precision="4">'
    write (UNIT_xdmf, "(A)")&
                           & '          1.0 1.0 1.0'
    write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
    write (UNIT_xdmf, "(A)")&
                           & '      </Geometry>'

    do var = 1, n_fld_vars
      write (UNIT_xdmf, "(A)")&
                           & '      <Attribute Name="' // trim(fld_vars(var)) // '" Center="Node">'
      write (UNIT_xdmf, "(A,I10,I10,I10,A)")&
                           & '        <DataItem Format="HDF" Dimensions="',&
                           & nk, nj, ni,&
                           & '" NumberType="Float" Precision="4">'
      write (UNIT_xdmf, "(A,A)")&
                           & '          flds.tot.' // trim(stepchar) // ':/',&
                           & trim(fld_vars(var))
      write (UNIT_xdmf, "(A)")&
                           & '        </DataItem>'
      write (UNIT_xdmf, "(A)")&
                           & '      </Attribute>'
    end do

    write (UNIT_xdmf, "(A)")&
                           & '    </Grid>'
    write (UNIT_xdmf, "(A)")&
                           & '  </Domain>'
    write (UNIT_xdmf, "(A)")&
                           & '</Xdmf>'
    close (UNIT_xdmf)
  end subroutine writeXDMF_hdf5

  subroutine writeFields_hdf5(step, time)
    implicit none
    integer, intent(in)               :: step, time
    character(len=STR_MAX)            :: stepchar, filename
    integer(HID_T)                    :: file_id, dset_id(40), filespace(40), memspace, plist_id
    integer                           :: error, f, s
    integer(kind=2)                   :: i, j, k
    logical                           :: writing_intQ, writing_lgarrQ
    integer                           :: dataset_rank = 3
    integer(HSSIZE_T), dimension(3)   :: offsets
    integer(HSIZE_T), dimension(3)    :: global_dims, blocks
    real                              :: ex0, ey0, ez0, bx0, by0, bz0
    real                              :: jx0, jy0, jz0

    ! downsampling variables
    integer           :: this_x0, this_y0, this_z0, this_sx, this_sy, this_sz
    integer           :: i_start, i_end, j_start, j_end, k_start, k_end
    integer           :: offset_i, offset_j, offset_k
    integer(kind=2)   :: i1, j1, k1
    integer           :: n_i, n_j, n_k, glob_n_i, glob_n_j, glob_n_k

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
    if (output_flds_istep .eq. 1) then
      offset_i = this_x0;   offset_j = this_y0;   offset_k = this_z0
      n_i = this_sx - 1;    n_j = this_sy - 1;    n_k = this_sz - 1
      glob_n_i = global_mesh%sx
      glob_n_j = global_mesh%sy
      glob_n_k = global_mesh%sz

      i_start = 0; j_start = 0; k_start = 0
    else
      i_start = 0; i_end = 0
      offset_i = 0; n_i = 0
      glob_n_i = 1

      j_start = 0; j_end = 0
      offset_j = 0; n_j = 0
      glob_n_j = 1

      k_start = 0; k_end = 0
      offset_k = 0; n_k = 0
      glob_n_k = 1
      #if defined(oneD) || defined (twoD) || defined (threeD)
        offset_i = CEILING(REAL(this_x0) / REAL(output_flds_istep))
        i_start = CEILING(REAL(this_x0) / REAL(output_flds_istep)) * output_flds_istep - this_x0
        i_end = (CEILING(REAL(this_x0 + this_sx) / REAL(output_flds_istep)) - 1) * output_flds_istep - this_x0
        n_i = (i_end - i_start) / output_flds_istep
        glob_n_i = CEILING(REAL(global_mesh%sx) / REAL(output_flds_istep))
        glob_n_i = MAX(1, glob_n_i)
      #endif
      #if defined(twoD) || defined (threeD)
        offset_j = CEILING(REAL(this_y0) / REAL(output_flds_istep))
        j_start = CEILING(REAL(this_y0) / REAL(output_flds_istep)) * output_flds_istep - this_y0
        j_end = (CEILING(REAL(this_y0 + this_sy) / REAL(output_flds_istep)) - 1) * output_flds_istep - this_y0
        n_j = (j_end - j_start) / output_flds_istep
        glob_n_j = CEILING(REAL(global_mesh%sy) / REAL(output_flds_istep))
        glob_n_j = MAX(1, glob_n_j)
      #endif
      #if defined(threeD)
        offset_k = CEILING(REAL(this_z0) / REAL(output_flds_istep))
        k_start = CEILING(REAL(this_z0) / REAL(output_flds_istep)) * output_flds_istep - this_z0
        k_end = (CEILING(REAL(this_z0 + this_sz) / REAL(output_flds_istep)) - 1) * output_flds_istep - this_z0
        n_k = (k_end - k_start) / output_flds_istep
        glob_n_k = CEILING(REAL(global_mesh%sz) / REAL(output_flds_istep))
        glob_n_k = MAX(1, glob_n_k)
      #endif
    end if

    if ((mpi_rank .eq. 0) .and. write_xdmf) then
      call writeXDMF_hdf5(step, time, glob_n_i, glob_n_j, glob_n_k)
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

    call h5open_f(error)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, h5comm, h5info, error)
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    do f = 1, n_fld_vars
      call h5screate_simple_f(dataset_rank, global_dims, filespace(f), error)
    end do

    do f = 1, n_fld_vars
      call prepareFieldForOutput(fld_vars(f), writing_lgarrQ)

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
            i = i_start + i1 * output_flds_istep
            j = j_start + j1 * output_flds_istep
            k = k_start + k1 * output_flds_istep
            call selectFieldForOutput(fld_vars(f), i1, j1, k1, i, j, k, writing_lgarrQ)
          end do
        end do
      end do

      ! Write the dataset collectively
      call h5dwrite_f(dset_id(f), H5T_NATIVE_REAL, sm_arr(0 : n_i, 0 : n_j, 0 : n_k), global_dims, error, &
                    & file_space_id = filespace(f), mem_space_id = memspace, xfer_prp = plist_id)

      call h5dclose_f(dset_id(f), error)
      call h5sclose_f(filespace(f), error)
    end do

    call h5sclose_f(memspace, error)
    call h5pclose_f(plist_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)
  end subroutine writeFields_hdf5

  #endif

end module m_writefldstot
