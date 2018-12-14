#include "../defs.F90"

module m_domain
  use m_globalnamespace
  implicit none

  type :: box
    integer    :: x0, y0, z0
    integer    :: sx, sy, sz
    integer    :: nghost
  end type box

  type(box), pointer                :: this_meshblock
  type(box)                         :: global_mesh
  type(box), allocatable, target    :: meshblocks(:)
end module m_domain
