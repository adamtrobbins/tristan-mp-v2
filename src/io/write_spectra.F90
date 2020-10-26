#include "../defs.F90"

module m_writespectra
  #ifdef HDF5
    use hdf5
    use m_globalnamespace, only: h5comm, h5info
  #endif
  use m_globalnamespace, only: output_dir_name, mpi_rank
  use m_outputnamespace, only: glob_spectra,&
                             & spec_num, spec_min, spec_max, spec_log_bins
  use m_aux
  use m_errors, only: throwError
  use m_domain
  use m_fields
  use m_particles
  use m_exchangearray, only: exchangeArray

  #ifdef RADIATION
    use m_outputnamespace, only: glob_rad_spectra
  #endif

  #ifdef GCA
    use m_outputnamespace, only: glob_gca_spectra
  #endif

  implicit none

  !--- PRIVATE functions -----------------------------------------!
  #ifdef HDF5
    private :: writeSpectra_hdf5
  #endif
  !...............................................................!

contains
  subroutine writeSpectra(step, time)
    implicit none
    integer, intent(in)               :: step, time
    #ifdef HDF5
      call writeSpectra_hdf5(step, time)
    #endif
    call printDiag((mpi_rank .eq. 0), "...writeSpectra()", .true.)
  end subroutine writeSpectra

  #ifdef HDF5
    subroutine writeSpectra_hdf5(step, time)
      implicit none
      integer, intent(in)               :: step, time
      character(len=STR_MAX)            :: stepchar, filename
      integer                           :: error, s, i, datarank
      integer(HID_T)                    :: file_id, dset_id, dspace_id
      integer(HSIZE_T), dimension(1)    :: data_dims
      character(len=6)                  :: dsetname
      real, allocatable, dimension(:)   :: bin_data

      datarank = 1
      data_dims(1) = spec_num

      ! only root rank writes the spectra file
      if (mpi_rank .eq. 0) then
        ! saving the energy bins
        allocate(bin_data(spec_num))
        do i = 1, spec_num
          bin_data(i) = spec_min + (REAL(i - 0.5) / REAL(spec_num)) * (spec_max - spec_min)
          if (spec_log_bins) then
            bin_data(i) = exp(bin_data(i))
          endif
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
          call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                         & dset_id, error)
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL, glob_spectra(s,:), data_dims, error)
          call h5dclose_f(dset_id, error)
          call h5sclose_f(dspace_id, error)

          #ifdef GCA
            ! writing spectra:
            dsetname = 'nbor' // trim(STR(s))
            call h5screate_simple_f(datarank, data_dims, dspace_id, error)
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                           & dset_id, error)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, glob_gca_spectra(s,:), data_dims, error)
            call h5dclose_f(dset_id, error)
            call h5sclose_f(dspace_id, error)

            ! writing spectra:
            dsetname = 'ngca' // trim(STR(s))
            call h5screate_simple_f(datarank, data_dims, dspace_id, error)
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                           & dset_id, error)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, glob_gca_spectra(nspec + s,:), data_dims, error)
            call h5dclose_f(dset_id, error)
            call h5sclose_f(dspace_id, error)
          #endif

          #ifdef RADIATION
            if (allocated(glob_rad_spectra)) then
              ! writing bins:
              dsetname = 'er' // trim(STR(s))
              call h5screate_simple_f(datarank, data_dims, dspace_id, error)
              call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                             & dset_id, error)
              call h5dwrite_f(dset_id, H5T_NATIVE_REAL, bin_data, data_dims, error)
              call h5dclose_f(dset_id, error)
              call h5sclose_f(dspace_id, error)

              ! writing spectra:
              dsetname = 'nr' // trim(STR(s))
              call h5screate_simple_f(datarank, data_dims, dspace_id, error)
              call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
                             & dset_id, error)
              call h5dwrite_f(dset_id, H5T_NATIVE_REAL, glob_rad_spectra(s,:), data_dims, error)
              call h5dclose_f(dset_id, error)
              call h5sclose_f(dspace_id, error)
              glob_rad_spectra(s,:) = 0.0
            end if
          #endif
        end do

        ! Close the file
        call h5fclose_f(file_id, error)
        ! Close FORTRAN interface
        call h5close_f(error)

        if (allocated(bin_data)) deallocate(bin_data)
      end if

      if (allocated(glob_spectra)) deallocate(glob_spectra)
      #ifdef GCA
        if (allocated(glob_gca_spectra)) deallocate(glob_gca_spectra)
      #endif
    end subroutine writeSpectra_hdf5
  #endif

end module m_writespectra
