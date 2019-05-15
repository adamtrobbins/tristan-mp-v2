#include "../defs.F90"

module m_aux
  use m_globalnamespace
  implicit none
  integer                  :: dseed
  ! integer, dimension(0:15)        :: state
  ! integer                         :: rand_ind

  interface STR
    module procedure intToStr
    module procedure realToStr
  end interface STR

  !--- PRIVATE functions -----------------------------------------!
  private :: intToStr, realToStr
  !...............................................................!
contains
  subroutine printDiag(bool, msg, prepend)
    implicit none
    character(len=*), intent(in)  :: msg
    logical, intent(in)           :: bool
    logical, optional, intent(in) :: prepend
    character(len=STR_MAX)        :: dummy
    integer                       :: sz, i, ierr
    #ifdef DEBUG
      if (bool) then
        sz = len(trim(msg))
        if (present(prepend)) then
          if (prepend) then
            dummy = '...'
            sz = sz + 3
          end if
        else
          dummy = ''
        end if
        dummy = trim(dummy) // trim(msg)
        do i = 1, (38 - sz)
          dummy = trim(dummy) // '.'
        end do
        dummy = trim(dummy) // '[OK]'
        print *, trim(dummy)
      end if
    #endif
  end subroutine printDiag

  subroutine printReport(bool, msg, prepend)
    implicit none
    character(len=*), intent(in)  :: msg
    logical, intent(in)           :: bool
    logical, optional, intent(in) :: prepend
    character(len=STR_MAX)        :: dummy
    integer                       :: sz, i, ierr
    if (bool) then
      sz = len(trim(msg))
      if (present(prepend)) then
        if (prepend) then
          dummy = '...'
          sz = sz + 3
        end if
      else
        dummy = ''
      end if
      dummy = trim(dummy) // trim(msg)
      do i = 1, (38 - sz)
        dummy = trim(dummy) // '.'
      end do
      dummy = trim(dummy) // '[OK]'
      print *, trim(dummy)
    end if
  end subroutine printReport

  subroutine printTime(dt_arr, msg, fullstep)
    implicit none
    character(len=*), intent(in)          :: msg
    character(len=STR_MAX)                :: dummy, dummy1
    real(kind=8), intent(in)              :: dt_arr(:)
    real, optional, intent(in)            :: fullstep
    real                                  :: dt_mean, dt_max, dt_min
    integer                               :: pcent_max, pcent_min, sz, sz1, i
    dt_mean = SUM(dt_arr) * 1000 / mpi_size
    dt_max = MAXVAL(dt_arr) * 1000
    dt_min = MINVAL(dt_arr) * 1000
    pcent_max = (dt_max - dt_mean) * 200 / (dt_max + dt_mean)
    pcent_min = (dt_mean - dt_min) * 200 / (dt_mean + dt_min)
    if (present(fullstep)) then
      if (dt_mean / fullstep .lt. 1e-3) then
        dt_mean = 0; pcent_max = 0; pcent_min = 0
      end if
    end if

    sz = len(msg)
    dummy(1 : sz) = msg
    do i = sz + 1, 19
      dummy(i : i) = ' '
    end do

    dummy1 = trim(STR(dt_mean))
    sz = len(trim(dummy1))
    dummy(20 : 20 + sz - 1) = trim(dummy1)
    do i = 20 + sz, 29
      dummy(i : i) = ' '
    end do
    dummy(30 : 36) = '[ms] (+'

    dummy1 = trim(STR(pcent_max))
    sz = len(trim(dummy1))
    dummy(37 : 37 + sz - 1) = trim(dummy1)
    dummy(37 + sz : 37 + sz + 1) = ' %'
    do i = 37 + sz + 2, 40
      dummy(i : i) = ' '
    end do

    dummy(41 : 44) = ' | -'
    dummy1 = trim(STR(pcent_min))
    sz = len(trim(dummy1))
    dummy(45 : 45 + sz - 1) = trim(dummy1)
    dummy(45 + sz : 45 + sz + 2) = ' %)'
    if (present(fullstep)) then
      do i = 45 + sz + 3, 49
        dummy(i : i) = ' '
      end do
      dummy(54 : 54) = '['
      dummy1 = trim(STR(dt_mean * 100 / fullstep))
      sz1 = len(trim(dummy1))
      dummy(55 : 55 + sz1 - 1) = trim(dummy1)
      dummy(55 + sz1 : 55 + sz1 + 2) = ' %]'
      do i = 55 + sz1 + 3, 80
        dummy(i : i) = ' '
      end do
    else
      do i = 45 + sz + 3, 80
        dummy(i : i) = ' '
      end do
    end if
    print *, dummy(1:80)
  end subroutine printTime

  function intToStr(my_int) result(string)
    implicit none
    integer, intent(in)       :: my_int
    character(:), allocatable :: string
    character(len=STR_MAX)    :: temp
    write(temp, '(i0)') my_int
    string = trim(temp)
  end function intToStr

  function realToStr(my_real) result(string)
    implicit none
    real, intent(in)          :: my_real
    character(:), allocatable :: string
    character(len=STR_MAX)    :: temp
    write(temp, '(G0.2)') my_real
    string = trim(temp)
  end function realToStr

  function STRtoINT(my_str) result(my_int)
    implicit none
    character(len=*), intent(in)  :: my_str
    integer                       :: my_int
    read (my_str, *) my_int
  end function STRtoINT

  !***********************************************************************
  ! init_random_seed() subroutine enables to avoid the repeating series of
  ! number given by random_number.
  ! Reference: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
  !***********************************************************************

  real function random(dseed)
    implicit none
    integer, intent(in) :: dseed
    call random_number(random)
  end function random

  subroutine initializeRandomSeed(rank)
    implicit none
    integer, intent(in) :: rank
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(COUNT = clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(PUT = seed)
    deallocate(seed)
  end subroutine initializeRandomSeed

  !******************************

  ! real function random(dseed)
  ! 	integer :: a, b, c, d, e
  ! 	integer, intent(in) :: dseed
  !   random = rand(dseed)
  !   print *, random
  !   return
  ! 	! a = state(rand_ind)
  ! 	! c = state(IAND((rand_ind + 13), 15))
  ! 	! b = IEOR(IEOR(a,c),IEOR(ISHFT(a, 16), ISHFT(c,15)))
  ! 	! c = state(IAND((rand_ind + 9), 15))
  ! 	! c = IEOR(c, ISHFT(c, -11))
  ! 	! state(rand_ind) = IEOR(b, c)
  ! 	! a = state(rand_ind)
  ! 	! d = IEOR(a, IAND(ISHFT(a, 5), 3661901092))
  ! 	! rand_ind = IAND(rand_ind + 15, 15)
  ! 	! a = state(rand_ind)
  ! 	! state(rand_ind) = IEOR(IEOR(IEOR(a, b), IEOR(d, ISHFT(a, 2))), IEOR(ISHFT(b, 18), ISHFT(c, 28)))
  ! 	! e = state(rand_ind)
  ! 	! e = abs(e)
  ! 	! random= real(e) / 2147483648.
  ! 	! return
  ! 	! real(dprec), intent(inout)  :: dseed
  !   ! random = 0.1
  ! 	! integer      :: I
  ! 	! real(dprec)  :: S2P31, S2P31M, seed
  ! 	! DATA            S2P31M/2147483647.D0/, S2P31/2147483648.D0/
  !   !
  ! 	! seed = dseed
  !   !
  ! 	! seed = DMOD(16807.D0 * seed, S2P31M)
  ! 	! random = seed / S2P31
  ! 	! dseed = seed
  ! 	! return
  ! end function

  ! subroutine initializeRandomSeed(rank)
  ! 	implicit none
  !   integer, intent(in) :: rank
  !   dseed = 123457
  !   call srand(dseed + rank)
  ! 	! dseed = dseed + rank
  ! end subroutine initializeRandomSeed

end module m_aux
