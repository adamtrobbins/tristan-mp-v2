#include "defs.F90"

module m_qednamespace
  use m_globalnamespace
  implicit none

  #ifdef BWPAIRPRODUCTION
    ! BW pair production parameters
    real            :: BW_tau
    integer         :: BW_interval, BW_electron_sp, BW_positron_sp
    integer         :: BW_algorithm
  #endif

  #ifdef COMPTONSCATTERING
    ! Compton scattering parameters
    real            :: Compton_tau
    integer         :: Compton_interval, Compton_algorithm
    logical         :: Compton_el_recoil
    real(kind=8)    :: Thomson_lim
  #endif

  #ifdef PAIRANNIHILATION
    ! Pair annihilation parameters
    real    :: Annihilation_tau
    integer :: Annihilation_interval, Annihilation_photon_sp
    ! integer :: Annihilation_algorithm
  #endif

end module m_qednamespace
