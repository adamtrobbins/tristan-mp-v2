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
  type :: phi_bin
    ! left and right bounds of the sub-bin (for debugging)
    real                 :: phi_min, phi_max
    ! number of particles in the sub-bin
    integer              :: npart
    ! particle indices in this sub-bin
    integer, allocatable :: indices(:)
  end type phi_bin

  type :: theta_bin
    ! upper and lower boundes of the bin (for debugging)
    real                       :: theta_min, theta_max
    ! number of sub-bins in phi
    integer                    :: n_phi_bins
    ! sub-bins in phi
    type(phi_bin), allocatable :: phi_bins(:)
  end type theta_bin

  type :: energy_bin
    ! upper and lower bounds for the energy bin
    real                          :: e_min, e_max
    ! number of sub-bins in theta
    integer                       :: n_theta_bins
    ! sub-bins in theta
    type(theta_bin), allocatable  :: theta_bins(:)
  end type energy_bin
  ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

  ! polar bin opening angle (half of it)
  real                            :: th0_bin
  ! number of bins
  integer                         :: n_angular_bins, n_energy_bins
  ! min & max energy for merging particles
  real                            :: dwn_energy_min, dwn_energy_max
  type(energy_bin), allocatable   :: energy_bins(:)

  !--- PRIVATE variables/functions -------------------------------!
  private :: th0_bin
  private :: findEnergyBin, findThetaBin, findPhiBin
  private :: initializeEnergyBins, initializeThetaBins, initializePhiBins
  !...............................................................!
contains

  subroutine binParticlesOnTile(tile)
    implicit none
    integer :: p
    integer :: energy_ind, theta_ind, phi_ind
    real    :: UUx, UUy, UUz, EE
    real    :: log10_e_max, log10_e_min
    real    :: u_theta, u_phi
    integer :: dummy_int

    ! FIX: this is for photons only

    type(particle_tile), intent(in) :: tile
    call initializeEnergyBins(tile%npart_sp)
    do p = 1, tile%npart_sp
      UUx = tile%u(p); UUy = tile%v(p); UUz = tile%w(p)
      EE = sqrt(UUx**2 + UUy**2 + UUz**2)
      UUx = UUx / EE; UUy = UUy / EE; UUz = UUz / EE
      if ((EE .ge. energy_bins(0)%e_min) .and.&
        & (EE .lt. energy_bins(n_energy_bins - 1)%e_max)) then
        u_theta = asin(UUz)
        u_phi = atan2(UUy, UUx)
        if (u_phi .lt. 0) u_phi = u_phi + 2 * M_PI
        call findEnergyBin(EE, energy_ind)
        call findThetaBin(u_theta,&
                        & energy_bins(energy_ind)%n_theta_bins,&
                        & theta_ind)
        dummy_int = energy_bins(energy_ind)%theta_bins(theta_ind)%n_phi_bins
        call findPhiBin(u_phi,&
                      & dummy_int,&
                      & phi_ind)

        dummy_int =&
            & energy_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%npart
        energy_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%&
                                  &indices(dummy_int) = p
        energy_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%npart = dummy_int + 1
      end if
    end do
  end subroutine binParticlesOnTile

  subroutine findEnergyBin(energy, en_ind)
    implicit none
    integer, intent(out)  :: en_ind
    integer               :: e_b
    do e_b = 0, n_energy_bins - 1
      if ((energy .lt. energy_bins(e_b)%e_max) .and.&
         &(energy .ge. energy_bins(e_b)%e_min))then
        en_ind = e_b - 1
        exit
      end if
    end do
  end subroutine findEnergyBin

  subroutine findThetaBin(u_theta, n_th_bins, th_ind)
    implicit none
    real, intent(in)      :: u_theta
    integer, intent(in)   :: n_th_bins
    integer, intent(out)  :: th_ind
    real                  :: d_theta

    d_theta = (M_PI - 2 * th0_bin) / n_th_bins
    if (u_theta .le. -0.5 * M_PI + th0_bin) then
      ! if on the southern pole bin
      th_ind = 0
    else if (u_theta .gt. 0.5 * M_PI - th0_bin) then
      ! if on the northern pole bin
      th_ind = n_th_bins + 1
    else
      th_ind = INT((u_theta + 0.5 * M_PI - th0_bin) / d_theta) + 1
      #ifdef DEBUG
        if ((th_ind .ge. n_th_bins) .or. (th_ind .le. 0)) then
          call throwError('Something is wrong in `findThetaBin()`')
        end if
      #endif
    end if
  end subroutine findThetaBin

  subroutine findPhiBin(u_phi, n_ph_bins, ph_bin)
    implicit none
    real, intent(in)      :: u_phi
    integer, intent(in)   :: n_ph_bins
    integer, intent(out)  :: ph_bin
    ph_bin = INT(u_phi * n_ph_bins / (2 * M_PI))
    #ifdef DEBUG
      if ((ph_bin .lt. 0) .or. (ph_bin .ge. n_ph_bins)) then
        call throwError('Something is wrong in `findPhiBin()`')
      end if
    #endif
  end subroutine findPhiBin

  subroutine initializeEnergyBins(nparts_in_tile)
    implicit none
    integer, intent(in) :: nparts_in_tile
    integer             :: e_b
    real, allocatable   :: lognorm(:)
    real                :: e_min, e_max

    th0_bin = 0.5 * M_PI / n_angular_bins

    ! initializing randomized energy bins ...
    ! ... intervals of which have lognormal distribution
    call log_normal(n_energy_bins + 1, lognorm)
    ! randomize maximum energy bin from `E` to `10 * E`...
    ! ... for more randomness

    e_min = dwn_energy_min
    e_max = 10**log10(dwn_energy_max) + random(dseed)

    allocate(energy_bins(0 : n_energy_bins - 1))
    do e_b = 0, n_energy_bins - 1
      energy_bins(e_b)%e_min = (e_min + (e_max - e_min) * lognorm(e_b))
      energy_bins(e_b)%e_max = (e_min + (e_max - e_min) * lognorm(e_b + 1))
      energy_bins(e_b)%n_theta_bins = n_angular_bins
      call initializeThetaBins(energy_bins(e_b)%theta_bins,&
                             & energy_bins(e_b)%n_theta_bins,&
                             & nparts_in_tile)
    end do
  end subroutine initializeEnergyBins

  subroutine initializeThetaBins(th_bins, n_th_bins, nparts_in_tile)
    implicit none
    type(theta_bin), allocatable, intent(out)   :: th_bins(:)
    integer, intent(in)                         :: n_th_bins, nparts_in_tile
    real                                        :: d_theta
    integer                                     :: th_b, ph_b
    real                                        :: theta_mid

    ! bins go from `0 -> n_th_bins + 1`...
    ! ... with `n_th_bins + 2` bins overall
    allocate(th_bins(0 : n_th_bins + 1))

    d_theta = (M_PI - 2 * th0_bin) / n_th_bins

    ! south polar bin
    th_bins(0)%theta_min = -0.5 * M_PI
    th_bins(0)%theta_max = -0.5 * M_PI + th0_bin
    th_bins(0)%n_phi_bins = 1
    call initializePhiBins(th_bins(0)%phi_bins,&
                         & th_bins(0)%n_phi_bins, -0.5 * M_PI,&
                         & nparts_in_tile)
    ! intermediate `theta` bins
    do th_b = 1, n_th_bins
      th_bins(th_b)%theta_min = -0.5 * M_PI + th0_bin + d_theta * (th_b - 1)
      th_bins(th_b)%theta_max = -0.5 * M_PI + th0_bin + d_theta * th_b
      theta_mid = 0.5 * (th_bins(th_b)%theta_min + th_bins(th_b)%theta_max)
      th_bins(th_b)%n_phi_bins = 2 * n_th_bins * cos(theta_mid)
      call initializePhiBins(th_bins(th_b)%phi_bins,&
                           & th_bins(th_b)%n_phi_bins, theta_mid,&
                           & nparts_in_tile)
    end do
    ! north polar bin
    th_bins(n_th_bins + 1)%theta_min = -0.5 * M_PI - th0_bin
    th_bins(n_th_bins + 1)%theta_max = -0.5 * M_PI
    th_bins(n_th_bins + 1)%n_phi_bins = 1
    call initializePhiBins(th_bins(n_th_bins + 1)%phi_bins,&
                         & th_bins(n_th_bins + 1)%n_phi_bins, 0.5 * M_PI,&
                         & nparts_in_tile)
  end subroutine initializeThetaBins

  subroutine initializePhiBins(ph_bins, n_ph_bins, th_mid, nparts_in_tile)
    implicit none
    type(phi_bin), allocatable, intent(out) :: ph_bins(:)
    integer, intent(in)                     :: n_ph_bins, nparts_in_tile
    real, intent(in)                        :: th_mid
    integer                                 :: ph_b
    real                                    :: d_phi, d_phi_0

    ! number of `phi` bins in the equatorial plane
    d_phi_0 = M_PI / n_th_bins

    ! `phi` bins go from `0 -> n_ph_bins - 1`...
    ! ... with `n_ph_bins` bins overall
    allocate(ph_bins(0 : n_ph_bins - 1))
    if (abs(th_mid) .ne. M_PI * 0.5) then
      d_phi = d_phi_0 / cos(th_mid)
    else
      d_phi = 2 * M_PI
    end if
    do ph_b = 0, n_ph_bins - 1
      ph_bins(ph_b)%phi_min = d_phi * ph_b
      ph_bins(ph_b)%phi_max = d_phi * (ph_b + 1)
      ph_bins(ph_b)%npart = 0
      allocate(ph_bins(ph_b)%indices(nparts_in_tile))
    end do
  end subroutine initializePhiBins

  subroutine deinitializeEnergyBins()
    implicit none
    integer :: e_b, th_b, ph_b
    do e_b = 0, n_energy_bins - 1
      do th_b = 0, energy_bins(e_b)%n_theta_bins
        do ph_b = 0, energy_bins(e_b)%theta_bins(th_b)%n_phi_bins
          deallocate(energy_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%indices)
        end do
        deallocate(energy_bins(e_b)%theta_bins(th_b)%phi_bins)
      end do
      deallocate(energy_bins(e_b)%theta_bins)
    end do
    deallocate(energy_bins(e_b))
  end subroutine deinitializeEnergyBins

#endif
end module m_momentumbinning
