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
  #ifdef FASTFILTERING
    integer               :: ntimes_main, ntimes_sec, filter_main_w, filter_sec_w
    real, allocatable     :: window_main(:), window_sec(:)

    private :: nfilter_main, nfitler_sec
  #endif

  private :: filterInX, filterInY, filterInAll
  #ifdef threeD
    private :: filterInZ
  #endif
  !...............................................................!
contains

  #ifdef FASTFILTERING
    subroutine computeWindowOfSizeN(window, N)
      implicit none
      ! this function basically computes ...
      ! ... the nonzero elements of the following sparse matrix (raised to the N-th power):
      !
      ! | a b 0 0 0 ... 0 0 0 0 0 |**N
      ! | b a b 0 0 ... 0 0 0 0 0 |
      ! | 0 b a b 0 ... 0 0 0 0 0 |
      ! | 0 0 b a b ... 0 0 0 0 0 |
      ! | 0 0 0 b a ... 0 0 0 0 0 |
      ! | ...       ...       ... |
      ! | 0 0 0 0 0 ... a b 0 0 0 |
      ! | 0 0 0 0 0 ... b a b 0 0 |
      ! | 0 0 0 0 0 ... 0 b a b 0 |
      ! | 0 0 0 0 0 ... 0 0 b a b |
      ! | 0 0 0 0 0 ... 0 0 0 b a |
      !
      ! where `a = 1/2` and `b = 1/4`
      ! In simplest case when `N = 1` -> `window = (1/4, 1/2, 1/4)`
      integer, intent(in) :: N
      real, intent(inout) :: window(-N : N)
      real                :: w_, coeff
      integer             :: k, i, j
      real                :: a, b
      a = 0.5; b = 0.25

      do k = -N, N
        w_ = 0.0
        do i = 0, N
          do j = 0, N - i
            if (N - i - 2 * j .eq. k) then
              coeff = factorial(N) / (factorial(i) * factorial(j) * factorial(N - i - j))
              w_ = w_ + coeff * a**i * b**(n - i)
            end if
          end do
        end do
        window(k) = w_
      end do
    end subroutine computeWindowOfSizeN

    subroutine initializeFilters()
      implicit none
      ! find the window sizes
      if (nfilter .le. NGHOST) then
        if (nfilter .gt. 0) then
          ntimes_main = 1
        else
          ntimes_main = 0
        end if
        ntimes_sec = 0
        filter_main_w = nfilter
        filter_sec_w = 0
      else
        ntimes_main = INT(nfilter / NGHOST)
        if (MOD(nfilter, NGHOST) .eq. 0) then
          ntimes_sec = 0
          filter_sec_w = 0
        else
          ntimes_sec = 1
          filter_sec_w = nfilter - ntimes_main * NGHOST
        end if
        filter_main_w = NGHOST
      end if
      ! allocate arrays to store the weights ...
      ! ... & compute the window weights
      if (filter_main_w .gt. 0) then
        allocate(window_main(-filter_main_w : filter_main_w))
        call computeWindowOfSizeN(window_main, filter_main_w)
      end if
      if (filter_sec_w .gt. 0) then
        allocate(window_sec(-filter_sec_w : filter_sec_w))
        call computeWindowOfSizeN(window_sec, filter_sec_w)
      end if
    end subroutine initializeFilters

    subroutine filterCurrents()
      implicit none
      integer :: iter
      do iter = 1, ntimes_main
        call filterInAll(window_main, filter_main_w)
      end do
      if (ntimes_sec .eq. 1) then
        call filterInAll(window_sec, filter_sec_w)
      end if
      call printDiag((mpi_rank .eq. 0), "filterCurrents()", .true.)
    end subroutine filterCurrents

    subroutine filterInAll(window, window_size)
      implicit none
      integer, intent(in) :: window_size
      real, intent(in)    :: window(-window_size : window_size)
      call filterInX(jx, window, window_size)
      call filterInX(jy, window, window_size)
      call filterInX(jz, window, window_size)
      call exchangeCurrents(.true.)

      ! call filterInY(jx, window, window_size)
      ! call filterInY(jy, window, window_size)
      ! call filterInY(jz, window, window_size)
      ! call exchangeCurrents(.true.)

      ! #ifdef threeD
      !   call filterInZ(jx, window, window_size)
      !   call filterInZ(jy, window, window_size)
      !   call filterInZ(jz, window, window_size)
      !   call exchangeCurrents(.true.)
      ! #endif
    end subroutine

    subroutine filterInX(arr, window, window_size)
      implicit none
      real, intent(inout) :: arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                               & fldBoundZ)
      integer, intent(in) :: window_size
      real, intent(in)    :: window(-window_size : window_size)
      real                :: filter_tmp(0 : window_size), dot_product_
      integer             :: i, j, k, t_, i_, k_
      integer             :: imin_, imax_

      k = 0
      do j = -NGHOST, this_meshblock%ptr%sy - 1 + NGHOST
        filter_tmp(0 : window_size - 1) = arr(-NGHOST : window_size - 1 - NGHOST, j, k)
        t_ = -1
        do i = window_size - NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
          if (t_ .gt. window_size - 1) then
            t_ = -1
          end if
          dot_product_ = 0.0
          k_ = -window_size
          imin_ = MAX(i - window_size, -NGHOST)
          imax_ = MIN(i + window_size, this_meshblock%ptr%sx - 1 + NGHOST)
          do i_ = imin_, imax_
            dot_product_ = dot_product_ + arr(i_, j, k) * window(k_)
            k_ = k_ + 1
          end do
          if (t_ .eq. -1) then
            filter_tmp(window_size) = dot_product_
          else
            filter_tmp(t_) = dot_product_
          end if
          arr(i - window_size, j, k) = filter_tmp(t_ + 1)
          t_ = t_ + 1
        end do
      end do
    end subroutine

  #else
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

      ! FIXTHISBACK
      ! call filterInY(jx, do_n_times)
      ! call filterInY(jy, do_n_times)
      ! call filterInY(jz, do_n_times)
      ! call exchangeCurrents(.true.)
      !
      ! #ifdef threeD
      !   call filterInZ(jx, do_n_times)
      !   call filterInZ(jy, do_n_times)
      !   call filterInZ(jz, do_n_times)
      !   call exchangeCurrents(.true.)
      ! #endif
    end subroutine

    subroutine filterInX(arr, do_n_times)
      implicit none
      real, intent(inout) :: arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                               & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                               & fldBoundZ)
      integer, intent(in) :: do_n_times
      real                :: tmp2, tmp1
      integer             :: i, j, k, n_pass
      if (do_n_times .gt. NGHOST) then
        call throwError('ERROR: `filterInX()` called with `do_n_times` > NGHOST.')
      end if
      do n_pass = 1, do_n_times
        #ifndef threeD
          k = 0
          do j = -NGHOST + n_pass, this_meshblock%ptr%sy - 1 + NGHOST - n_pass
            tmp2 = arr(-NGHOST + n_pass - 1, j, k)
            i = -NGHOST + n_pass
            do while (i .le. this_meshblock%ptr%sx - 1 + NGHOST - n_pass)
              tmp1 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k)
              arr(i - 1, j, k) = tmp2
              i = i + 1
              tmp2 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k)
              arr(i - 1, j, k) = tmp1
              i = i + 1
            end do
            arr(this_meshblock%ptr%sx - 1 + NGHOST - n_pass, j, k) = tmp2
          end do
        #else
          do k = -NGHOST + n_pass, this_meshblock%ptr%sz - 1 + NGHOST - n_pass
            do j = -NGHOST + n_pass, this_meshblock%ptr%sy - 1 + NGHOST - n_pass
              tmp2 = arr(-NGHOST + n_pass - 1, j, k)
              i = -NGHOST + n_pass
              do while (i .le. this_meshblock%ptr%sx - 1 + NGHOST - n_pass)
                tmp1 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k)
                arr(i - 1, j, k) = tmp2
                i = i + 1
                tmp2 = 0.25 * arr(i - 1, j, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i + 1, j, k)
                arr(i - 1, j, k) = tmp1
                i = i + 1
              end do
              arr(this_meshblock%ptr%sx - 1 + NGHOST - n_pass, j, k) = tmp2
            end do
          end do
        #endif
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
      if (do_n_times .gt. NGHOST) then
        call throwError('ERROR: `filterInY()` called with `do_n_times` > NGHOST.')
      end if
      do n_pass = 1, do_n_times
        #ifndef threeD
          k = 0
          do i = -NGHOST + n_pass, this_meshblock%ptr%sx - 1 + NGHOST - n_pass
            tmp2 = arr(i, -NGHOST + n_pass - 1, k)
            j = -NGHOST + n_pass
            do while (j .lt.  this_meshblock%ptr%sy - 1 + NGHOST - n_pass)
              tmp1 = 0.25 * arr(i, j - 1, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j + 1, k)
              arr(i, j - 1, k) = tmp2
              j = j + 1
              tmp2 = 0.25 * arr(i, j - 1, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j + 1, k)
              arr(i, j - 1, k) = tmp1
              j = j + 1
            end do
            arr(i, this_meshblock%ptr%sy - 1 + NGHOST - n_pass, k) = tmp2
          end do
        #else
          do k = -NGHOST + n_pass, this_meshblock%ptr%sz - 1 + NGHOST - n_pass
            do i = -NGHOST + n_pass, this_meshblock%ptr%sx - 1 + NGHOST - n_pass
              tmp2 = arr(i, -NGHOST + n_pass - 1, k)
              j = -NGHOST + n_pass
              do while (j .lt.  this_meshblock%ptr%sy - 1 + NGHOST - n_pass)
                tmp1 = 0.25 * arr(i, j - 1, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j + 1, k)
                arr(i, j - 1, k) = tmp2
                j = j + 1
                tmp2 = 0.25 * arr(i, j - 1, k) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j + 1, k)
                arr(i, j - 1, k) = tmp1
                j = j + 1
              end do
              arr(i, this_meshblock%ptr%sy - 1 + NGHOST - n_pass, k) = tmp2
            end do
          end do
        #endif
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
        if (do_n_times .gt. NGHOST) then
          call throwError('ERROR: `filterInZ()` called with `do_n_times` > NGHOST.')
        end if
        do n_pass = 1, do_n_times
          do j = -NGHOST + n_pass, this_meshblock%ptr%sy - 1 + NGHOST - n_pass
            do i = -NGHOST + n_pass, this_meshblock%ptr%sx - 1 + NGHOST - n_pass
              tmp2 = arr(i, j, -NGHOST + n_pass - 1)
              k = -NGHOST + n_pass
              do while (k .lt. this_meshblock%ptr%sz - 1 + NGHOST - n_pass)
                tmp1 = 0.25 * arr(i, j, k - 1) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j, k + 1)
                arr(i, j, k - 1) = tmp2
                k = k + 1
                tmp2 = 0.25 * arr(i, j, k - 1) + 0.5 * arr(i, j, k) + 0.25 * arr(i, j, k + 1)
                arr(i, j, k - 1) = tmp1
                k = k + 1
              end do
              arr(i, j, this_meshblock%ptr%sz - 1 + NGHOST - n_pass) = tmp2
            end do
          end do
        end do
      end subroutine
    #endif
  #endif

end module m_filtering
