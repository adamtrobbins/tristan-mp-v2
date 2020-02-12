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
    real                 :: phi_min, phi_max
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

  !--- PRIVATE variables/functions -------------------------------!
  ! private :: findEnergyBin, findThetaBin, findPhiBin
  private :: initializeThetaBins, initializePhiBins
  !...............................................................!
contains

  ! subroutine binParticlesOnTile(tile)
  !   implicit none
  !   integer :: p
  !   integer :: energy_ind, theta_ind, phi_ind
  !   real    :: UUx, UUy, UUz, EE
  !   real    :: log10_e_max, log10_e_min
  !   real    :: u_theta, u_phi
  !   integer :: dummy_int
  !
  !   ! FIX: this is for photons only
  !
  !   type(particle_tile), intent(in) :: tile
  !
  !   do p = 1, tile%npart_sp
  !     UUx = tile%u(p); UUy = tile%v(p); UUz = tile%w(p)
  !     EE = sqrt(UUx**2 + UUy**2 + UUz**2)
  !     UUx = UUx / EE; UUy = UUy / EE; UUz = UUz / EE
  !     if ((EE .ge. energy_bins(0)%e_min) .and.&
  !       & (EE .lt. energy_bins(n_energy_bins - 1)%e_max)) then
  !       u_theta = asin(UUz)
  !       u_phi = atan2(UUy, UUx)
  !       if (u_phi .lt. 0) u_phi = u_phi + 2 * M_PI
  !       call findEnergyBin(EE, energy_ind)
  !       call findThetaBin(u_theta,&
  !                       & energy_bins(energy_ind)%n_theta_bins,&
  !                       & theta_ind)
  !       dummy_int = energy_bins(energy_ind)%theta_bins(theta_ind)%n_phi_bins
  !       call findPhiBin(u_phi,&
  !                     & dummy_int,&
  !                     & phi_ind)
  !
  !       dummy_int =&
  !           & energy_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%npart
  !       energy_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%&
  !                                 &indices(dummy_int) = p
  !       energy_bins(energy_ind)%theta_bins(theta_ind)%phi_bins(phi_ind)%npart = dummy_int + 1
  !     end if
  !   end do
  ! end subroutine binParticlesOnTile
  !
  ! ! - - - finding bins - - - - - - - - - - - - - - - - - - - - - - - -
  ! subroutine findEnergyBin(energy, en_ind)
  !   implicit none
  !   real, intent(in)      :: energy
  !   integer, intent(out)  :: en_ind
  !   integer               :: e_b
  !   en_ind = -1
  !   do e_b = 0, n_energy_bins - 1
  !     if ((energy .lt. energy_bins(e_b)%e_max) .and.&
  !        &(energy .ge. energy_bins(e_b)%e_min))then
  !       en_ind = e_b
  !       exit
  !     end if
  !   end do
  !   #ifdef DEBUG
  !     if ((en_ind .lt. 0) .or. (en_ind .ge. n_energy_bins)) then
  !       call throwError('Something is wrong in `findEnergyBin()`')
  !     end if
  !   #endif
  ! end subroutine findEnergyBin
  !
  ! subroutine findThetaBin(u_theta, n_th_bins, th_ind)
  !   implicit none
  !   real, intent(in)      :: u_theta
  !   integer, intent(in)   :: n_th_bins
  !   integer, intent(out)  :: th_ind
  !   real                  :: d_theta
  !
  !   d_theta = (M_PI - 2 * th0_bin) / n_th_bins
  !   if (u_theta .le. -0.5 * M_PI + th0_bin) then
  !     ! if on the southern pole bin
  !     th_ind = 0
  !   else if (u_theta .gt. 0.5 * M_PI - th0_bin) then
  !     ! if on the northern pole bin
  !     th_ind = n_th_bins + 1
  !   else
  !     th_ind = INT((u_theta + 0.5 * M_PI - th0_bin) / d_theta) + 1
  !     #ifdef DEBUG
  !       if ((th_ind .gt. n_th_bins) .or. (th_ind .le. 0)) then
  !         call throwError('Something is wrong in `findThetaBin()`')
  !       end if
  !     #endif
  !   end if
  ! end subroutine findThetaBin
  !
  ! subroutine findPhiBin(u_phi, n_ph_bins, ph_bin)
  !   implicit none
  !   real, intent(in)      :: u_phi
  !   integer, intent(in)   :: n_ph_bins
  !   integer, intent(out)  :: ph_bin
  !   ph_bin = INT(u_phi * n_ph_bins / (2 * M_PI))
  !   #ifdef DEBUG
  !     if ((ph_bin .lt. 0) .or. (ph_bin .ge. n_ph_bins)) then
  !       call throwError('Something is wrong in `findPhiBin()`')
  !     end if
  !   #endif
  ! end subroutine findPhiBin

  ! - - - initializing bins - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine initializeMomentumBins(momentum_bins, nparts_in_tile)
    implicit none
    integer, intent(in)                         :: nparts_in_tile
    type(momentumBin), intent(out), allocatable :: momentum_bins(:)
    integer                                     :: e_b
    real, allocatable                           :: lognorm(:)
    real                                        :: e_min, e_max

    ! initializing randomized energy bins ...
    ! ... intervals of which have lognormal distribution
    call log_normal(n_energy_bins + 1, lognorm)
    ! randomize maximum energy bin from `E` to `10 * E`...
    ! ... for more randomness
    e_min = dwn_energy_min
    e_max = 10**log10(dwn_energy_max) + random(dseed)

    allocate(momentum_bins(0 : n_energy_bins - 1))
    do e_b = 0, n_energy_bins - 1
      momentum_bins(e_b)%e_min = (e_min + (e_max - e_min) * lognorm(e_b))
      print *, 'daaaa', momentum_bins(e_b)%e_min
      momentum_bins(e_b)%e_max = (e_min + (e_max - e_min) * lognorm(e_b + 1))
      print *, 'daaaa', momentum_bins(e_b)%e_max
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
      ! call initializePhiBins(momentum_bin, momentum_bin%theta_bins(th_b), nparts_in_tile)
    end do
  end subroutine initializeThetaBins

  subroutine initializePhiBins(momentum_bin, theta_bin, nparts_in_tile)
    implicit none
    type(momentumBin), intent(inout)        :: momentum_bin
    type(thetaBin), intent(out)             :: theta_bin
    integer, intent(in)                     :: nparts_in_tile
    integer                                 :: ph_b
    real                                    :: d_phi, d_phi_0

    ! size of `phi` bins in the equatorial plane
    d_phi_0 = M_PI / momentum_bin%n_theta_bins

    ! `phi` bins go from `0 -> n_ph_bins - 1`...
    ! ... with `n_ph_bins` bins overall
    allocate(theta_bin%phi_bins(0 : theta_bin%n_phi_bins - 1))
    if (abs(theta_bin%theta_mid) .ne. M_PI * 0.5) then
      d_phi = d_phi_0 / cos(theta_bin%theta_mid)
    else
      d_phi = 2 * M_PI
    end if
    do ph_b = 0, theta_bin%n_phi_bins - 1
      theta_bin%phi_bins(ph_b)%phi_min = d_phi * ph_b
      theta_bin%phi_bins(ph_b)%phi_max = d_phi * (ph_b + 1)
      theta_bin%phi_bins(ph_b)%npart = 0
      allocate(theta_bin%phi_bins(ph_b)%indices(nparts_in_tile))
    end do
  end subroutine initializePhiBins

  ! subroutine deinitializeEnergyBins()
  !   implicit none
  !   integer :: e_b, th_b, ph_b
  !   do e_b = 0, n_energy_bins - 1
  !     do th_b = 0, energy_bins(e_b)%n_theta_bins
  !       do ph_b = 0, energy_bins(e_b)%theta_bins(th_b)%n_phi_bins
  !         deallocate(energy_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%indices)
  !       end do
  !       deallocate(energy_bins(e_b)%theta_bins(th_b)%phi_bins)
  !     end do
  !     deallocate(energy_bins(e_b)%theta_bins)
  !   end do
  !   deallocate(energy_bins)
  ! end subroutine deinitializeEnergyBins

#endif
end module m_momentumbinning
