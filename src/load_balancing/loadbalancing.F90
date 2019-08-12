#include "../defs.F90"

module m_loadbalancing
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  use m_staticlb
  use m_adaptivelb
  implicit none

  !--- PRIVATE variables/functions -------------------------------!
  private :: metaRedistInX, metaRedistInY, metaRedistInZ,&
           & accumulateLoads, balanceLoad
  !...............................................................!

contains
  subroutine redistributeMeshblocksSLB(spat_load_ptr)
    implicit none
    ! pointer to a user defined function ... 
    ! ... that computes the spatial loads
    procedure (spatialDistribution), pointer, intent(in) :: spat_load_ptr
    integer :: ntimes
    integer :: nit
    ntimes = 10
    
    do nit = 1, ntimes
      if (slb_x) then
        call computeLoadSLB(lb_load, spat_load_ptr)
        call accumulateLoads()
        new_meshblocks(:) = meshblocks(:)
        call metaRedistInX(slb_sxmin)
        meshblocks(:) = new_meshblocks(:)
      end if
      if (slb_y) then
        ! ...
      end if
      if (slb_z) then
        ! ...
      end if
    end do
    
  end subroutine redistributeMeshblocksSLB

  ! accumulate loads from all the sources
  subroutine accumulateLoads()
    implicit none
    integer :: ierr
    if (.not. allocated(lb_load_glob)) allocate(lb_load_glob(mpi_size))
    call MPI_ALLGATHER(lb_load, 1, MPI_INTEGER,&
                     & lb_load_glob, 1, MPI_INTEGER,&
                     & MPI_COMM_WORLD, ierr)
  end subroutine accumulateLoads
  
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
  
  subroutine metaRedistInX(sxmin)
    implicit none 
    integer, intent(in)   :: sxmin
    integer               :: delta_i, i, j, k, cnt
    integer               :: sx0_old, sx1_old, sx0_new, sx1_new
    integer               :: load_x0, load_x1
    
    if (.not. allocated(lb_group_x0)) allocate(lb_group_x0(sizey * sizez))
    if (.not. allocated(lb_group_x1)) allocate(lb_group_x1(sizey * sizez))

    do delta_i = 0, sizex - 2
      ! select the left and right domain slabs
      cnt = 1
      do k = 0, sizez - 1
        do j = 0, sizey - 1
          i = delta_i
          lb_group_x0(cnt) = indToRnk([i, j, k]) + 1
          i = delta_i + 1
          lb_group_x1(cnt) = indToRnk([i, j, k]) + 1
          cnt = cnt + 1
        end do
      end do
      ! now all the actions are between these two slabs

      #ifdef DEBUG
        do cnt = 1, sizey * sizez
          if ((new_meshblocks(lb_group_x0(cnt))%sx .ne. new_meshblocks(lb_group_x0(1))%sx) .or.&
            & (new_meshblocks(lb_group_x0(cnt))%x0 .ne. new_meshblocks(lb_group_x0(1))%x0) .or.&
            & (new_meshblocks(lb_group_x1(cnt))%sx .ne. new_meshblocks(lb_group_x1(1))%sx) .or.&
            & (new_meshblocks(lb_group_x1(cnt))%x0 .ne. new_meshblocks(lb_group_x1(1))%x0) .or.&
            & (new_meshblocks(lb_group_x0(cnt))%x0 + new_meshblocks(lb_group_x0(cnt))%sx .ne.&
              & new_meshblocks(lb_group_x1(cnt))%x0)) then
            call throwError('ERROR: wrong dimensions of slabs in `metaRedistInX()`')
          end if
        end do
      #endif

      ! cumulative loads for each slab
      load_x0 = 0; load_x1 = 0
      do cnt = 1, sizey * sizez
        load_x0 = load_x0 + lb_load_glob(lb_group_x0(cnt))
        load_x1 = load_x1 + lb_load_glob(lb_group_x1(cnt))
      end do
      ! old dimensions
      sx0_old = new_meshblocks(lb_group_x0(1))%sx
      sx1_old = new_meshblocks(lb_group_x1(1))%sx

      ! computing new dimensions
      call balanceLoad(load_x0, sx0_old, load_x1, sx1_old, sx0_new, sx1_new, sxmin)

      do cnt = 1, sizey * sizez
        new_meshblocks(lb_group_x0(cnt))%sx = sx0_new
        new_meshblocks(lb_group_x1(cnt))%sx = sx1_new
        new_meshblocks(lb_group_x1(cnt))%x0 =&
          & new_meshblocks(lb_group_x0(cnt))%x0 + new_meshblocks(lb_group_x0(cnt))%sx
      end do
    end do

    ! meshblocks(:) = new_meshblocks(:)

    deallocate(lb_group_x0)
    deallocate(lb_group_x1)
  end subroutine metaRedistInX

  subroutine metaRedistInY(symin)
    implicit none
    integer, intent(in)   :: symin
  end subroutine metaRedistInY

  subroutine metaRedistInZ(szmin)
    implicit none
    integer, intent(in)   :: szmin
  end subroutine metaRedistInZ

end module m_loadbalancing
