#include "defs.F90"

module m_qednamespace
  use m_globalnamespace
  implicit none

  real :: QED_tau0

  #ifdef BWPAIRPRODUCTION
    ! BW pair production parameters
    integer         :: BW_interval, BW_electron_sp, BW_positron_sp
    integer         :: BW_algorithm
  #endif

  #ifdef COMPTONSCATTERING
    ! Compton scattering parameters
    integer         :: Compton_interval, Compton_algorithm
    logical         :: Compton_el_recoil
    real(kind=8)    :: Thomson_lim
  #endif

  #ifdef PAIRANNIHILATION
    ! Pair annihilation parameters
    integer :: Annihilation_interval, Annihilation_photon_sp
    ! integer :: Annihilation_algorithm
  #endif

end module m_qednamespace
