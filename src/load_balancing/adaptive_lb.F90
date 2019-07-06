#include "../defs.F90"

module m_adaptivelb
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_particlelogistics
  use m_fields
  use m_helpers
  implicit none

  logical :: alb_x, alb_y, alb_z
  integer :: alb_sxmin, alb_symin, alb_szmin
  integer :: alb_int_x, alb_int_y, alb_int_z
  integer :: alb_start_x, alb_start_y, alb_start_z

  !--- PRIVATE functions -----------------------------------------!

  !...............................................................!

  !--- PRIVATE variables -----------------------------------------!
  private :: metaRedistInX, metaRedistInY, metaRedistInZ,&
           & computeLoad, balanceLoad
  !...............................................................!
contains
  subroutine redistributeDomain(step)
    implicit none
    integer, intent(in) :: step
    #ifdef ALB
    ! do meta redistribution (determine new dimensions)
    print *, step, alb_x, alb_start_x, alb_int_x, modulo(step, alb_int_x)
    if ((step .ge. alb_start_x) .and. (modulo(step, alb_int_x) .eq. 0)) then
      if (alb_x) then
        call metaRedistInX()
      end if
    end if

    ! exchange all particles
    ! exchange all fields
    ! reallocate all quantities
    #endif
  end subroutine redistributeDomain

  ! the following subroutine computes the load
  !   (e.g. # of particles) for any given domain
  ! this can be arbitrary function for which the algorithm
  !   balances the domain distribution
  subroutine computeLoad(load)
    implicit none
    integer, intent(out) :: load
    integer              :: s, ti, tj, tk
    load = 0
    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            load = load + species(s)%prtl_tile(ti, tj, tk)%npart_sp
          end do
        end do
      end do
    end do
  end subroutine computeLoad

  ! this routine defines the algorithm for load balancing
  !   it can be used in any direction for already predefined domain slabs
  subroutine balanceLoad(load0, s0_old, load1, s1_old,&
                       & s0_new, s1_new, smin)
    implicit none
    integer, intent(in)   :: smin
    integer, intent(in)   :: load0, s0_old, load1, s1_old
    integer, intent(out)  :: s0_new, s1_new

    if (load0 .lt. load1) then
      s0_new = INT(s0_old + s1_old * REAL(load1 - load0) / REAL(2.0 * load1))
      s1_new = s0_old + s1_old - s0_new
      if (s1_new .lt. smin) then
        s1_new = smin
        s0_new = s0_old + s1_old - s1_new
      end if
    else
      s1_new = INT(s1_old + s0_old * REAL(load0 - load1) / REAL(2.0 * load0))
      s0_new = s1_old + s0_old - s1_new
      if (s0_new .lt. smin) then
        s0_new = smin
        s1_new = s1_old + s0_old - s0_new
      end if
    end if

    #ifdef DEBUG
      if ((s0_new .lt. smin) .or. (s1_new .lt. smin)) then
        call throwError("ERROR: `s0/s1 < smin` in `balanceLoad()`")
      end if
      if ((s0_new + s1_new .ne. s0_old + s1_old)) then
        call throwError("ERROR: `s0 + s1` wrong in `balanceLoad()`")
      end if
    #endif
  end subroutine balanceLoad

  subroutine metaRedistInX()
    implicit none
    integer, allocatable  :: group_x0(:), group_x1(:)
    integer               :: delta_i, i, j, k, cnt
    integer               :: sx0_old, sx1_old, sx0_new, sx1_new
    integer               :: ierr, load, load_x0, load_x1
    integer, allocatable  :: load_glob(:)
    ! type(mesh), allocatable :: new_meshblocks(:)

    allocate(load_glob(mpi_size))
    allocate(group_x0(sizey * sizez), group_x1(sizey * sizez))

    call computeLoad(load)
    call MPI_ALLGATHER(load, 1, MPI_INTEGER,&
                     & load_glob, 1, MPI_INTEGER,&
                     & MPI_COMM_WORLD, ierr)

    do delta_i = 0, sizex - 2
      ! select the left and right domain slabs
      cnt = 1
      do k = 0, sizez - 1
        do j = 0, sizey - 1
          i = delta_i
          group_x0(cnt) = indToRnk([i, j, k]) + 1
          i = delta_i + 1
          group_x1(cnt) = indToRnk([i, j, k]) + 1
          cnt = cnt + 1
        end do
      end do
      ! now all the actions are between these two slabs

      #ifdef DEBUG
        do cnt = 1, sizey * sizez
          if ((new_meshblocks(group_x0(cnt))%sx .ne. new_meshblocks(group_x0(1))%sx) .or.&
            & (new_meshblocks(group_x0(cnt))%x0 .ne. new_meshblocks(group_x0(1))%x0) .or.&
            & (new_meshblocks(group_x1(cnt))%sx .ne. new_meshblocks(group_x1(1))%sx) .or.&
            & (new_meshblocks(group_x1(cnt))%x0 .ne. new_meshblocks(group_x1(1))%x0) .or.&
            & (new_meshblocks(group_x0(cnt))%x0 + new_meshblocks(group_x0(cnt))%sx .ne.&
              & new_meshblocks(group_x1(cnt))%x0)) then
            call throwError('ERROR: wrong dimensions of slabs in `metaRedistInX()`')
          end if
        end do
      #endif

      ! cumulative loads for each slab
      load_x0 = 0; load_x1 = 0
      do cnt = 1, sizey * sizez
        load_x0 = load_x0 + load_glob(group_x0(cnt))
        load_x1 = load_x1 + load_glob(group_x1(cnt))
      end do
      ! old dimensions
      sx0_old = new_meshblocks(group_x0(1))%sx
      sx1_old = new_meshblocks(group_x1(1))%sx

      ! computing new dimensions
      call balanceLoad(load_x0, sx0_old, load_x1, sx1_old, sx0_new, sx1_new, alb_sxmin)

      do cnt = 1, sizey * sizez
        new_meshblocks(group_x0(cnt))%sx = sx0_new
        new_meshblocks(group_x1(cnt))%sx = sx1_new
        new_meshblocks(group_x1(cnt))%x0 = new_meshblocks(group_x0(cnt))%x0 + new_meshblocks(group_x0(cnt))%sx
      end do
    end do

    ! meshblocks(:) = new_meshblocks(:)

    deallocate(group_x0)
    deallocate(group_x1)
    deallocate(load_glob)
  end subroutine metaRedistInX

  subroutine metaRedistInY()
    implicit none
  end subroutine metaRedistInY

  subroutine metaRedistInZ()
    implicit none
  end subroutine metaRedistInZ

end module m_adaptivelb
