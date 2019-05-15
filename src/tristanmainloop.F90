#include "defs.F90"

module m_mainloop
  use m_globalnamespace
  use m_aux
  use m_writeoutput
  use m_fldsolver
  use m_mover
  use m_currentdeposit
  use m_exchangeparts
  use m_exchangefields
  use m_exchangecurrents
  use m_particlelogistics
  use m_filtering
  use m_userfile
  use m_errors
  implicit none

  integer       :: timestep

  real(kind=8)  :: t_fullstep_1, t_fullstep_2, &
                 & t_movestep_1, t_movestep_2, &
                 & t_depositstep_1, t_depositstep_2, &
                 & t_filterstep_1, t_filterstep_2, &
                 & t_outputstep_1, t_outputstep_2


  !--- PRIVATE functions -----------------------------------------!
  private :: makeReport
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: t_fullstep_1, t_fullstep_2, &
           & t_movestep_1, t_movestep_2, &
           & t_depositstep_1, t_depositstep_2, &
           & t_filterstep_1, t_filterstep_2, &
           & t_outputstep_1, t_outputstep_2
  !...............................................................!
contains
  subroutine mainloop()
    implicit none
    integer       :: ierr, i
    integer       :: s, ti, tj, tk, p

    ! ADD needs to be changed for restart
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call printReport((mpi_rank .eq. 0), "Starting mainloop()")

    do timestep = 0, final_timestep
        t_fullstep_1 = MPI_WTIME()

      ! MAINLOOP >
      call exchangeFields(.true., .true.)
      call advanceBHalfstep()
      call exchangeFields(.false., .true.)

        t_movestep_1 = MPI_WTIME()
      call moveParticles()
        t_movestep_2 = MPI_WTIME()

      call advanceBHalfstep()
      call exchangeFields(.false., .true.)

      call advanceEFullstep()
      call exchangeFields(.true., .false.)

        t_depositstep_1 = MPI_WTIME()
      call depositCurrents()
        t_depositstep_2 = MPI_WTIME()

      call exchangeCurrents()

        t_filterstep_1 = MPI_WTIME()
      call filterCurrents()
        t_filterstep_2 = MPI_WTIME()

      call addCurrents()
      call exchangeFields(.true., .false.)

      call exchangeParticles()
      call clearGhostParticles()

      call userDriveParticles()

      t_outputstep_1 = MPI_WTIME(); t_outputstep_2 = MPI_WTIME()
      if (mod(timestep, output_interval) .eq. 0) then
        t_outputstep_1 = MPI_WTIME()
        call writeOutput(INT(timestep / output_interval), timestep)
        t_outputstep_2 = MPI_WTIME()
      end if
      ! </ MAINLOOP

        t_fullstep_2 = MPI_WTIME()

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (ierr .eq. MPI_SUCCESS) then
        call makeReport(timestep)
      end if
    end do
  end subroutine mainloop

  subroutine makeReport(tstep)
    implicit none
    integer, intent(in)           :: tstep
    integer                       :: ierr
    real                          :: fullstep
    real(kind=8), allocatable     :: dt_fullstep(:), dt_movestep(:), &
                                   & dt_depositstep(:), dt_filterstep(:), &
                                   & dt_outputstep(:)

    allocate(dt_fullstep(mpi_size), dt_movestep(mpi_size))
    allocate(dt_depositstep(mpi_size), dt_filterstep(mpi_size))
    allocate(dt_outputstep(mpi_size))

    call MPI_GATHER(t_fullstep_2 - t_fullstep_1, 1, MPI_REAL8,&
                  & dt_fullstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_movestep_2 - t_movestep_1, 1, MPI_REAL8,&
                  & dt_movestep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_depositstep_2 - t_depositstep_1, 1, MPI_REAL8,&
                  & dt_depositstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_filterstep_2 - t_filterstep_1, 1, MPI_REAL8,&
                  & dt_filterstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER(t_outputstep_2 - t_outputstep_1, 1, MPI_REAL8,&
                  & dt_outputstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. 0) then
      fullstep = SUM(dt_fullstep) * 1000 / mpi_size
      call printReport(.true., "timestep: " // STR(tstep))
      call printTime(dt_fullstep, "Full_step: ")
      call printTime(dt_movestep, "  move_step: ", fullstep)
      call printTime(dt_depositstep, "  deposit_step: ", fullstep)
      call printTime(dt_filterstep, "  filter_step: ", fullstep)
      call printTime(dt_outputstep, "  output_step: ", fullstep)
      print *, ""
    end if

    ! full # of particles ...

    deallocate(dt_fullstep, dt_movestep)
    deallocate(dt_depositstep, dt_filterstep)
    deallocate(dt_outputstep)
  end subroutine makeReport

end module m_mainloop
