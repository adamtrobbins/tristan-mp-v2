#include "defs.F90"

program tristan
  use m_initialize
  use m_mainloop
  use m_finalize

    use m_bincoupling

  implicit none
  !----- main code --------------------------
  integer :: num_couples, i
  type(couple), allocatable :: coupled_pairs(:)

  call initializeAll()

  call coupleParticlesOnTile(0, 0, 0,&
                           & (/1/), 1,&
                           & (/2/), 1,&
                           & coupled_pairs, num_couples)

  if (this_meshblock%ptr%rnk .eq. 0) then
    print *, "num pairs:", num_couples
    do i = 1, num_couples
      print *, i, coupled_pairs(i)%part_1%spec, coupled_pairs(i)%part_1%index, coupled_pairs(i)%part_2%spec, coupled_pairs(i)%part_2%index
    end do
  end if
  ! call mainloop()
  call finalizeAll()

  !..... main code ..........................
end program tristan
