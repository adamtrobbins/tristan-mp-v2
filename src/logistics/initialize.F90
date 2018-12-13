#include "../defs.F90"

module m_initialize
  use m_globalnamespace
  use m_readinput
  use m_communications
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: firstRankInitialize
  !...............................................................!
contains
  ! initialize all the necessary functions
  subroutine initializeAll()
    implicit none
    call readCommandlineArgs()
    call initializeCommunications()
    call firstRankInitialize()
  end subroutine initializeAll

  subroutine firstRankInitialize()
    ! create output/restart directories
    !   if does not already exist
    call system('mkdir -p '//trim(output_dir_name))
    call system('mkdir -p '//trim(restart_dir_name))
  end subroutine firstRankInitialize
end module m_initialize
