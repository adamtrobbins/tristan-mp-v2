#include "defs.F90"

module m_mainloop
  use m_globalnamespace
  use m_aux
  use m_communications
  use m_writeoutput
  use m_mover
  use m_exchangeparts
  use m_errors
  implicit none

  integer       :: timestep
contains
  subroutine mainloop()
    implicit none
    integer :: ierr

    ! ADD needs to be changed for restart
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call print_diag((mpi_rank .eq. 0), "Starting mainloop()")

    do timestep = 1, final_timestep
      call print_diag((mpi_rank .eq. 0), TAB // "time: " // STR(timestep))
      call print_diag((mpi_rank .eq. 0), TAB // "------------------------------------")
      call print_diag(.true., STR(mpi_rank) // " : " // STR(spp_(1)%npart_sp))

      call moveParticles()

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call exchangeParticles()
      call clearGhostParticles()
      if (mod(timestep, output_interval) .eq. 0) then
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call writeOutput(INT(timestep / output_interval))
      end if
      call print_diag((mpi_rank .eq. 0), "")
    end do

  end subroutine mainloop
end module m_mainloop
