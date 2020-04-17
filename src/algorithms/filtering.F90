#include "../defs.F90"

module m_filtering
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_fields
  use m_exchangefields
  use m_exchangecurrents

  !--- PRIVATE variables/functions -------------------------------!
  private :: filterInX, filterInY, filterInAll
  #ifdef threeD
    private :: filterInZ
  #endif
  !...............................................................!
contains
  subroutine filterCurrents()
    implicit none
    integer :: n_pass, iter
    n_pass = NGHOST
    iter = 0
    do while (.true.)
      if (n_pass .ge. nfilter) then
        if (nfilter .gt. 0) then
          ! filter `(nfilter - iter * NGHOST)` times
          call filterInAll(nfilter - iter * NGHOST)
        end if
        exit
      else
        ! filter `(NGHOST)` times
        call filterInAll(NGHOST)
        n_pass = n_pass + NGHOST
        iter = iter + 1
      end if
    end do
    call printDiag((mpi_rank .eq. 0), "filterCurrents()", .true.)
  end subroutine filterCurrents

  ! subroutine filterEfield(arr)
  ! ...
  ! end subroutine filterEfield

  subroutine filterInAll(do_n_times)
    implicit none
    integer, intent(in) :: do_n_times
    call filterInX(jx, do_n_times)
    call filterInX(jy, do_n_times)
    call filterInX(jz, do_n_times)
    call exchangeCurrents(.true.)

    call filterInY(jx, do_n_times)
    call filterInY(jy, do_n_times)
    call filterInY(jz, do_n_times)
    call exchangeCurrents(.true.)

    #ifdef threeD
      call filterInZ(jx, do_n_times)
      call filterInZ(jy, do_n_times)
      call filterInZ(jz, do_n_times)
      call exchangeCurrents(.true.)
    #endif
  end subroutine

  subroutine filterInX(arr, do_n_times)
    implicit none
    real, intent(inout) :: arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                             & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                             & fldBoundZ)
    integer, intent(in) :: do_n_times
    real                :: tmp2, tmp1
    integer             :: i, j, k, n_pass
    integer             :: imin, imax, jmin, jmax, kmin, kmax
    if (do_n_times .gt. NGHOST) then
      call throwError('ERROR: `filterInX()` called with `do_n_times` > NGHOST.')
    end if
    do n_pass = 1, do_n_times
      imin = -NGHOST + n_pass
      imax = this_meshblock%ptr%sx - 1 + NGHOST - n_pass
      jmin = -NGHOST + n_pass
      jmax = this_meshblock%ptr%sy - 1 + NGHOST - n_pass
      
      #ifndef threeD
        kmin = 0 
        kmax = 0
      #else
        kmin = -NGHOST + n_pass
        kmax = this_meshblock%ptr%sz - 1 + NGHOST - n_pass
      #endif
      do k = kmin, kmax
        do j = jmin, jmax
          tmp2 = arr(imin - 1, j, k)
          do i = imin, imax, 2
            tmp1 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k)
            arr(i - 1, j, k) = tmp2
            tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i + 1, j, k) + 0.25 * arr(i + 2, j, k)
            arr(i, j, k) = tmp1
          end do
          arr(imax + 1, j, k) = tmp2
        end do
      end do
    end do
  end subroutine

  subroutine filterInY(arr, do_n_times)
    implicit none
    real, intent(inout) :: arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                             & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                             & fldBoundZ)
    integer, intent(in) :: do_n_times
    real                :: tmp2, tmp1
    integer             :: i, j, k, n_pass
    integer             :: imin, imax, jmin, jmax, kmin, kmax
    if (do_n_times .gt. NGHOST) then
      call throwError('ERROR: `filterInY()` called with `do_n_times` > NGHOST.')
    end if
    do n_pass = 1, do_n_times
      imin = -NGHOST + n_pass
      imax = this_meshblock%ptr%sx - 1 + NGHOST - n_pass
      jmin = -NGHOST + n_pass
      jmax = this_meshblock%ptr%sy - 1 + NGHOST - n_pass
      #ifndef threeD
        kmin = 0 
        kmax = 0
      #else
        kmin = -NGHOST + n_pass
        kmax = this_meshblock%ptr%sz - 1 + NGHOST - n_pass
      #endif
      do k = kmin, kmax
        do i = imin, imax
          tmp2 = arr(i, jmin - 1, k)
          do j = jmin, jmax, 2
            tmp1 = 0.25 * arr(i, j - 1, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j + 1, k)
            arr(i, j - 1, k) = tmp2
            tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i, j + 1, k) + 0.25 * arr(i, j + 2, k)
            arr(i, j, k) = tmp1
          end do
          arr(i, jmax + 1, k) = tmp2
        end do
      end do
    end do
  end subroutine

  #ifdef threeD
    subroutine filterInZ(arr, do_n_times)
      implicit none
      real, intent(inout) :: arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                               & fldBoundZ)
      integer, intent(in) :: do_n_times
      real                :: tmp2, tmp1
      integer             :: i, j, k, n_pass
      integer             :: imin, imax, jmin, jmax, kmin, kmax
      if (do_n_times .gt. NGHOST) then
        call throwError('ERROR: `filterInZ()` called with `do_n_times` > NGHOST.')
      end if
      do n_pass = 1, do_n_times
        imin = -NGHOST + n_pass
        imax = this_meshblock%ptr%sx - 1 + NGHOST - n_pass
        jmin = -NGHOST + n_pass
        jmax = this_meshblock%ptr%sy - 1 + NGHOST - n_pass
        kmin = -NGHOST + n_pass
        kmax = this_meshblock%ptr%sz - 1 + NGHOST - n_pass
        do j = jmin, jmax
          do i = imin, imax
            tmp2 = arr(i, j, kmin - 1)
            do k = kmin, kmax, 2
              tmp1 = 0.25 * arr(i, j, k - 1) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j, k + 1)
              arr(i, j, k - 1) = tmp2
              tmp2 = 0.25 * arr(i, j, k) + 0.5 * arr(i, j, k + 1) + 0.25 * arr(i, j, k + 2)
              arr(i, j, k) = tmp1
            end do
            arr(i, j, kmax + 1) = tmp2
          end do
        end do
      end do
    end subroutine
  #endif
end module m_filtering
