#include "../defs.F90"

module m_domain
  use m_globalnamespace

  implicit none

  ! type :: box
  !   integer*8       :: x0, y0, z0
  !   integer*8       :: sx, sy, sz
  !   integer         :: nghost
  ! end type box
  !
  ! type(box)               :: local_grid
  ! type(box)               :: global_grid
  ! type(box), allocatable  :: all_local_grids(:)
contains
  subroutine initializeDomain()
    implicit none
  end subroutine initializeDomain
end module m_domain
