#include "defs.F90"

module m_mainloop
  use m_globalnamespace
  use m_aux
  use m_communications
  use m_writeoutput
  use m_fldsolver
  use m_mover
  use m_exchangeparts
  use m_exchangefields
  use m_errors
  implicit none

  integer       :: timestep

  real(kind=8)  :: t_fullstep_1, t_fullstep_2

  !--- PRIVATE functions -----------------------------------------!
  private :: makeReport
  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: t_fullstep_1, t_fullstep_2
  !...............................................................!
contains
  subroutine mainloop()
    implicit none
    integer       :: ierr, i

    ! ADD needs to be changed for restart
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call printReport((mpi_rank .eq. 0), "Starting mainloop()")

    do timestep = 0, final_timestep
      t_fullstep_1 = MPI_WTIME()
      ! MAINLOOP >

      call exchangeFields(.true., .true.)
      call advanceBHalfstep()
      call exchangeFields(.false., .true.)

      call moveParticles()

      call advanceBHalfstep()
      call exchangeFields(.false., .true.)

      call advanceEFullstep()
      call exchangeFields(.true., .false.)
      
      ! deposit

      call exchangeParticles()
      call clearGhostParticles()

      if (mod(timestep, output_interval) .eq. 0) then
        call writeOutput(INT(timestep / output_interval), timestep)
      end if

      ! </ MAINLOOP

      t_fullstep_2 = MPI_WTIME()

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (ierr .eq. MPI_SUCCESS) then
        call makeReport(timestep)
      end if
    end do
  end subroutine mainloop

  #ifdef DEBUG
    subroutine showParticles()
      implicit none
      integer :: s, p
      print *, "PRINTING PARTICLES FOR RNK:", mpi_rank
      do s = 1, nspec
        print *, trim(TAB) // "PRINTING SPECIES:", s
        do p = 1, spp_(s)%npart_sp
          print *, p, sp_(s)%xi(p) + sp_(s)%dx(p), sp_(s)%yi(p) + sp_(s)%dy(p),&
                 & ISIGN(1, sp_(s)%proc(p)) * (ABS(sp_(s)%proc(p)) * 100 + sp_(s)%ind(p))
        end do
      end do
    end subroutine showParticles

    subroutine showField()
      implicit none
      integer :: j
      print *, "PRINTING BX FIELD:", mpi_rank
      do j = this_meshblock%ptr%sy - 1 + NGHOST, -NGHOST, -1
        print *, bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, j, 0)
      end do
    end subroutine showField
  #endif

  subroutine makeReport(tstep)
    implicit none
    integer, intent(in)           :: tstep
    integer :: ierr
    real(kind=8), allocatable     :: dt_fullstep(:)

    allocate(dt_fullstep(mpi_size))

    call MPI_GATHER(t_fullstep_2 - t_fullstep_1, 1, MPI_REAL8,&
                  & dt_fullstep, 1, MPI_REAL8,&
                  & 0, MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. 0) then
      call printReport(.true., "timestep: " // STR(tstep) // TAB // TAB // TAB // TAB // "[OK]")
      call printReport(.true., "-------------")
      call printTime(dt_fullstep, "Full_step: ")
      call printReport(.true., "")
    end if

    ! full # of particles ...
    ! call MPI_REDUCE(spp_(1)%npart_sp, allparts, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    deallocate(dt_fullstep)
  end subroutine makeReport

  subroutine printTime(dt_arr, msg)
    implicit none
    character(len=*), intent(in)          :: msg
    real(kind=8), intent(in)              :: dt_arr(:)
    real                                  :: dt_mean, dt_max, dt_min
    integer                               :: pcent_max, pcent_min
    character(len=STR_MAX)                :: tab_
    dt_mean = SUM(dt_arr) * 1000 / mpi_size
    dt_max = MAXVAL(dt_arr) * 1000
    dt_min = MINVAL(dt_arr) * 1000
    pcent_max = (dt_max - dt_mean) * 200 / (dt_max + dt_mean)
    pcent_min = (dt_mean - dt_min) * 200 / (dt_mean + dt_min)
    if (len(STR(dt_mean)) .lt. 4) then
      tab_ = TAB // TAB
    else
      tab_ = TAB
    end if
    call printReport(.true., msg // STR(dt_mean) // trim(tab_) // " +" //&
                                            & STR(pcent_max) // "% | -" // STR(pcent_min) // "%")
  end subroutine printTime

end module m_mainloop
