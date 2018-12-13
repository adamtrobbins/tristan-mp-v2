#include "../defs.F90"

module m_communications
  use m_globalnamespace
  use m_domain

  implicit none

  integer           :: rank, size0, sizex, sizey, sizez
  integer, private  :: ierr, statsize
contains
  ! subroutine initializeCommunications()
  !   implicit none
  !   sizex = getIntInput('node_configuration', 'sizex', 1)
  !   sizey = getIntInput('node_configuration', 'sizey', 1)
  !   #ifdef threeD
  !     sizez = getIntInput('node_configuration', 'sizez', 1)
  !   #else
  !     sizez = 1
  !   #endif
  !   size0 = sizex * sizey * sizez
  !   allocate(all_local_grids(size0))
  !   global_grid%x0 =
  ! end subroutine initializeCommunications
end module m_communications
