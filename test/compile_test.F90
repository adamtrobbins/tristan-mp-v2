program CompileTest
  use mpi_f08
  use hdf5
  implicit none

  integer :: mpi_rank, mpi_size, ierr

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)

  write (*, '(A,I0,A,I0)') 'process ', mpi_rank, ' of ', mpi_size

  call h5open_f(ierr)
  call h5close_f(ierr)

  call MPI_FINALIZE(ierr)
end program CompileTest