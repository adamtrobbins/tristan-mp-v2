#include "../defs.F90"

module m_momentumbinning
#ifdef DOWNSAMPLING

  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  implicit none

  ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  ! Spherical momenta binning...
  ! ... auxiliary types
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! binning is logarithmic in energy,..
  ! ... uniform in `theta`,..
  ! ... and approximately uniform in `phi`
  ! each `energy` bin contains `theta`-bins...
  ! ... which further contains `phi`-bins
  type :: phiBin
    ! left and right bounds of the sub-bin (for debugging)
    real                 :: phi_min, phi_max, phi_mid
    ! number of particles in the sub-bin
    integer              :: npart
    ! particle indices in this sub-bin
    integer, allocatable :: indices(:)
  end type phiBin

  type :: thetaBin
    ! upper and lower boundes of the bin (for debugging)
    real                       :: theta_min, theta_max, theta_mid
    ! number of sub-bins in phi
    integer                    :: n_phi_bins
    ! sub-bins in phi
    type(phiBin), allocatable  :: phi_bins(:)
  end type thetaBin

  type :: momentumBin
    ! upper and lower bounds for the energy bin
    real                          :: e_min, e_max
    ! polar bin opening angle (half of it)
    real                          :: th0_bin
    ! number of sub-bins in theta
    integer                       :: n_theta_bins
    ! sub-bins in theta
    type(thetaBin), allocatable   :: theta_bins(:)
  end type momentumBin
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! number of bins
  integer                         :: n_angular_bins, n_energy_bins
  ! min & max energy for merging particles
  real                            :: dwn_energy_min, dwn_energy_max
  ! log or linearly distributed energy bins
  logical                         :: log_e_bins

  !--- PRIVATE variables/functions -------------------------------!
  private :: findEnergyBin, findThetaBin, findPhiBin
  private :: initializeThetaBins, initializePhiBins
  !...............................................................!
contains

  subroutine binParticlesOnTile(momentum_bins, tile, ax1, ax2, ang)
    implicit none
    ! FIX: this is for photons only
    type(momentumBin), allocatable, intent(inout) :: momentum_bins(:)
    type(particle_tile), intent(in)               :: tile
    real, intent(in)                              :: ax1, ax2, ang
    integer                                       :: s
    logical                                       :: masslessQ

    integer :: p
    integer :: energy_ind, theta_ind, phi_ind
    real    :: prtl_ux, prtl_uy, prtl_uz, prtl_energy, prtl_gamma
    real    :: u_theta, u_phi
    integer :: dummy_int

    s = tile%spec
    if ((species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0)) then
      masslessQ = .true.
    else
      masslessQ = .false.
    end if

    do p = 1, tile%npart_sp
      prtl_ux = tile%u(p); prtl_uy = tile%v(p); prtl_uz = tile%w(p)
      ! rotation (not really random, because axis and angles are passed)
      call rotateRandomlyIn3D(prtl_ux, prtl_uy, prtl_uz, ax1, ax2, ang)
      prtl_energy = sqrt(prtl_ux**2 + prtl_uy**2 + prtl_uz**2)
      if (.not. masslessQ) then
        prtl_gamma = sqrt(1.0 + prtl_ux**2 + prtl_uy**2 + prtl_uz**2)
      end if
      prtl_ux = prtl_ux / prtl_energy
      prtl_uy = prtl_uy / prtl_energy
      prtl_uz = prtl_uz / prtl_energy
      if (.not. masslessQ) then
        ! binning by `lorentz-factor - 1` (kinetic energy)
        prtl_energy = prtl_gamma - 1.0
      end if
      if ((prtl_energy .ge. momentum_bins(0)%e_min) .and.&
        & (prtl_energy .lt. momentum_bins(n_energy_bins - 1)%e_max)) then
        u_theta = asin(prtl_uz)
        u_phi = atan2(prtl_uy, prtl_ux)
        if (u_phi .lt. 0) u_phi = u_phi + 2 * M_PI
        call findEnergyBin(momentum_bins,&
                        & prtl_energy, energy_ind)
        call findThetaBin(momentum_bins(energy_ind),&
                        & u_theta, theta_ind)
        dummy_int = momentum_bins(energy_ind)%theta_bins(theta_ind)%n_phi_bins
        call findPhiBin(momentum_bins(energy_ind)%theta_bins(theta_ind),&
                      & u_phi, phi_ind)
        dummy_int =&
            & momentum_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%npart
        momentum_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%npart = dummy_int + 1
        momentum_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%&
                                  &indices(dummy_int + 1) = p
      end if
    end do
  end subroutine binParticlesOnTile

  ! - - - finding bins - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine findEnergyBin(momentum_bins, energy, en_ind)
    implicit none
    type(momentumBin), intent(in), allocatable  :: momentum_bins(:)
    real, intent(in)                            :: energy
    integer, intent(out)                        :: en_ind
    integer                                     :: e_b
    en_ind = -1
    do e_b = 0, n_energy_bins - 1
      if ((energy .lt. momentum_bins(e_b)%e_max) .and.&
         &(energy .ge. momentum_bins(e_b)%e_min))then
        en_ind = e_b
        exit
      end if
    end do
    #ifdef DEBUG
      if ((en_ind .lt. 0) .or. (en_ind .ge. n_energy_bins)) then
        call throwError('Something is wrong in `findEnergyBin()`')
      end if
    #endif
  end subroutine findEnergyBin

  subroutine findThetaBin(momentum_bin, u_theta, th_ind)
    implicit none
    type(momentumBin), intent(in) :: momentum_bin
    real, intent(in)              :: u_theta
    integer, intent(out)          :: th_ind
    real                          :: d_theta

    d_theta = (M_PI - 2 * momentum_bin%th0_bin) / momentum_bin%n_theta_bins
    if (u_theta .le. -0.5 * M_PI + momentum_bin%th0_bin) then
      ! if on the southern pole bin
      th_ind = 0
    else if (u_theta .gt. 0.5 * M_PI - momentum_bin%th0_bin) then
      ! if on the northern pole bin
      th_ind = momentum_bin%n_theta_bins + 1
    else
      th_ind = INT((u_theta + 0.5 * M_PI - momentum_bin%th0_bin) / d_theta) + 1
      th_ind = MIN(th_ind, momentum_bin%n_theta_bins)
      #ifdef DEBUG
        if ((th_ind .gt. momentum_bin%n_theta_bins) .or. (th_ind .le. 0)) then
          call throwError('Something is wrong in `findThetaBin()`')
        end if
      #endif
    end if
  end subroutine findThetaBin

  subroutine findPhiBin(theta_bin, u_phi, ph_bin)
    implicit none
    type(thetaBin), intent(in)  :: theta_bin
    real, intent(in)            :: u_phi
    integer, intent(out)        :: ph_bin
    ph_bin = INT(u_phi * theta_bin%n_phi_bins / (2 * M_PI))
    ph_bin = MAX(0, MIN(ph_bin, theta_bin%n_phi_bins - 1))
    #ifdef DEBUG
      if ((ph_bin .lt. 0) .or. (ph_bin .ge. theta_bin%n_phi_bins)) then
        print *, u_phi, ph_bin, theta_bin%theta_min, theta_bin%theta_max
        call throwError('Something is wrong in `findPhiBin()`')
      end if
    #endif
  end subroutine findPhiBin

  ! - - - initializing bins - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine initializeMomentumBins(momentum_bins, nparts_in_tile)
    implicit none
    integer, intent(in)                         :: nparts_in_tile
    type(momentumBin), intent(out), allocatable :: momentum_bins(:)
    integer                                     :: e_b
    real, allocatable                           :: normbins(:)
    real                                        :: e_min, e_max

    if (log_e_bins) then
      ! initializing randomized energy bins ...
      ! ... intervals of which have lognormal distribution
      call log_normal(n_energy_bins + 1, normbins)
      ! randomize maximum energy bin from `E` to `10 * E`...
      ! ... for more randomness
      e_min = dwn_energy_min
      e_max = 10**log10(dwn_energy_max) + random(dseed)
    else
      ! initializing randomized energy bins ...
      ! ... intervals of which have linear distribution
      call lin_normal(n_energy_bins + 1, normbins)
      ! randomize maximum energy bin from `E` to `2 * E`...
      ! ... for more randomness
      e_min = dwn_energy_min
      e_max = dwn_energy_max * (1.0 + random(dseed))
    end if

    allocate(momentum_bins(0 : n_energy_bins - 1))
    do e_b = 0, n_energy_bins - 1
      momentum_bins(e_b)%e_min = (e_min + (e_max - e_min) * normbins(e_b + 1))
      momentum_bins(e_b)%e_max = (e_min + (e_max - e_min) * normbins(e_b + 2))
      momentum_bins(e_b)%th0_bin = 0.5 * M_PI / n_angular_bins
      momentum_bins(e_b)%n_theta_bins = n_angular_bins
      call initializeThetaBins(momentum_bins(e_b), nparts_in_tile)
    end do
  end subroutine initializeMomentumBins

  subroutine initializeThetaBins(momentum_bin, nparts_in_tile)
    implicit none
    type(momentumBin), intent(inout)            :: momentum_bin
    integer, intent(in)                         :: nparts_in_tile
    real                                        :: d_theta
    integer                                     :: th_b, ph_b, dummy4
    real                                        :: dummy1, dummy2, dummy3

    ! bins go from `0 -> n_theta_bins + 1`...
    ! ... with `n_theta_bins + 2` bins overall
    allocate(momentum_bin%theta_bins(0 : momentum_bin%n_theta_bins + 1))

    d_theta = (M_PI - 2 * momentum_bin%th0_bin) / momentum_bin%n_theta_bins

    do th_b = 0, momentum_bin%n_theta_bins + 1
      if (th_b .eq. 0) then
        ! south polar bin
        dummy1 = -0.5 * M_PI
        dummy2 = -0.5 * M_PI + momentum_bin%th0_bin
        dummy3 = -0.5 * M_PI
        dummy4 = 1
      else if (th_b .eq. momentum_bin%n_theta_bins + 1) then
        ! north polar bin
        dummy1 = 0.5 * M_PI - momentum_bin%th0_bin
        dummy2 = 0.5 * M_PI
        dummy3 = 0.5 * M_PI
        dummy4 = 1
      else
        dummy1 = -0.5 * M_PI + momentum_bin%th0_bin + d_theta * (th_b - 1)
        dummy2 = -0.5 * M_PI + momentum_bin%th0_bin + d_theta * th_b
        dummy3 = 0.5 * (dummy1 + dummy2)
        dummy4 = INT(2 * momentum_bin%n_theta_bins * cos(dummy3))
      end if
      momentum_bin%theta_bins(th_b)%theta_min = dummy1
      momentum_bin%theta_bins(th_b)%theta_max = dummy2
      momentum_bin%theta_bins(th_b)%theta_mid = dummy3
      momentum_bin%theta_bins(th_b)%n_phi_bins = dummy4
      call initializePhiBins(momentum_bin, momentum_bin%theta_bins(th_b), nparts_in_tile)
    end do
  end subroutine initializeThetaBins

  subroutine initializePhiBins(momentum_bin, theta_bin, nparts_in_tile)
    implicit none
    type(momentumBin), intent(inout)        :: momentum_bin
    type(thetaBin), intent(inout)           :: theta_bin
    integer, intent(in)                     :: nparts_in_tile
    integer                                 :: ph_b
    real                                    :: d_phi

    ! `phi` bins go from `0 -> n_ph_bins - 1`...
    ! ... with `n_ph_bins` bins overall
    allocate(theta_bin%phi_bins(0 : theta_bin%n_phi_bins - 1))
    if (abs(theta_bin%theta_mid) .ne. M_PI * 0.5) then
      d_phi = 2 * M_PI / REAL(theta_bin%n_phi_bins)
    else
      d_phi = 2 * M_PI
    end if
    do ph_b = 0, theta_bin%n_phi_bins - 1
      theta_bin%phi_bins(ph_b)%phi_min = d_phi * ph_b
      theta_bin%phi_bins(ph_b)%phi_max = d_phi * (ph_b + 1)
      theta_bin%phi_bins(ph_b)%phi_mid = d_phi * (ph_b + 0.5)
      theta_bin%phi_bins(ph_b)%npart = 0
      allocate(theta_bin%phi_bins(ph_b)%indices(nparts_in_tile))
    end do
  end subroutine initializePhiBins

#endif
end module m_momentumbinning
