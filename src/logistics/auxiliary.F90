#include "../defs.F90"

module m_aux
  use m_globalnamespace
  implicit none
  real(dprec)                     :: dseed
  integer, dimension(0:15)        :: state
  integer                         :: rand_ind

contains
  real function random(dseed)
    ! FIX look at precision here
  	implicit none
  	integer                     :: a, b, c, d, e
  	real(dprec), intent(in)     :: dseed
  	a = state(rand_ind)
  	c = state(IAND((rand_ind + 13), 15))
  	b = IEOR(IEOR(a,c), IEOR(ISHFT(a, 16), ISHFT(c, 15)))
  	c = state(IAND((rand_ind + 9), 15))
  	c = IEOR(c, ISHFT(c, -11))
  	state(rand_ind) = IEOR(b, c)
  	a = state(rand_ind)
  	d = IEOR(a, IAND(ISHFT(a, 5), 3661901092_8))
  	rand_ind = IAND(rand_ind + 15, 15)
  	a = state(rand_ind)
  	state(rand_ind) = IEOR(IEOR(IEOR(a, b), IEOR(d, ISHFT(a, 2))), IEOR(ISHFT(b, 18), ISHFT(c, 28)))
  	e = state(rand_ind)
  	e = abs(e)
  	random = REAL(e) / 2147483648.
  end function

  subroutine initializeRandomSeed(rank)
  	implicit none
    integer, intent(in)  :: rank
  	integer              :: n
  	real                 :: temp
  	dseed = 123457.D0
  	dseed = dseed + rank
  	state(0) = 970391000 + rank * 100
  	state(1) = 1062140000 + rank * 100
  	state(2) = 2010910000 + rank * 100
  	state(3) = 2002910000 + rank * 100
  	state(4) = 534192249 + rank * 100
  	state(5) = 16109577 + rank * 100
  	state(6) = 96071135 + rank * 100
  	state(7) = 866448509 + rank * 100
  	state(8) = 724170151 + rank * 100
  	state(9) = 896443319 + rank * 100
  	state(10) = 181263722 + rank * 100
  	state(11) = 661530307 + rank * 100
  	state(12) = 42831488 + rank * 100
  	state(13) = 361611 + rank * 100
  	state(14) = 58206989 + rank * 100
  	state(15) = 461110074 + rank * 100
  	rand_ind = 0
  	do n = 1, 10000000
  		temp = random(dseed)
    enddo
  end subroutine initializeRandomSeed
end module m_aux
