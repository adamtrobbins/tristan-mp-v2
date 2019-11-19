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

  call coupleParticlesOnTile(1, 1, 1,&
                           & (/1/), 1,&
                           & (/2/), 1,&
                           & coupled_pairs, num_couples)

  if (this_meshblock%ptr%rnk .eq. 0) then
    print *, "num pairs:", num_couples
    print *, "num spec1:", species(1)%prtl_tile(1, 1, 1)%npart_sp
    print *, "num spec2:", species(2)%prtl_tile(1, 1, 1)%npart_sp
    do i = 1, num_couples
      print *, i, coupled_pairs(i)%part_1%spec, coupled_pairs(i)%part_1%index, coupled_pairs(i)%part_2%spec, coupled_pairs(i)%part_2%index
    end do
  end if
  ! call mainloop()
  call finalizeAll()

  !..... main code ..........................
end program tristan
