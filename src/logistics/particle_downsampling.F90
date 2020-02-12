#include "../defs.F90"

module m_particledownsampling
#ifdef DOWNSAMPLING

  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_errors
  use m_domain
  use m_particles
  use m_particlelogistics
  use m_momentumbinning
  implicit none

  integer           :: dwn_interval

  !--- PRIVATE variables/functions -------------------------------!
  private :: downsampleParticles, downsampleOnTile,&
           & downsampleBinnedParticles
  !...............................................................!
contains
  subroutine downsamplingStep(timestep)
    implicit none
    integer, intent(in) :: timestep
    if (modulo(timestep, dwn_interval) .eq. 0) then
      call downsampleParticles()
    end if
    call printDiag((mpi_rank .eq. 0), "downsamplingStep()", .true.)
  end subroutine downsamplingStep

  subroutine downsampleParticles()
    implicit none
    integer :: s, ti, tj, tk
    integer :: bin_limit

    ! if # of particles on a tile is less than this limit...
    ! ... the algorithm won't be downsampling
    bin_limit = INT(sqrt(2.0 * n_energy_bins * n_angular_bins**2))

    do s = 1, nspec
      if (species(s)%dwn_sp) then
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              if (species(s)%prtl_tile(ti, tj, tk)%npart_sp .gt. bin_limit) then
                call downsampleOnTile(species(s)%prtl_tile(ti, tj, tk))
              end if
            end do
          end do
        end do
      end if
    end do

  end subroutine downsampleParticles

  subroutine downsampleOnTile(tile)
    implicit none
    type(particle_tile), intent(inout)  :: tile
    integer                             :: energy_ind, theta_ind, phi_ind
    type(momentumBin), allocatable      :: momentum_bins(:)

    call initializeMomentumBins(momentum_bins, tile%npart_sp)
    call binParticlesOnTile(momentum_bins, tile)
    call downsampleBinnedParticles(momentum_bins, tile)
  end subroutine downsampleOnTile

  subroutine downsampleBinnedParticles(momentum_bins, tile)
    implicit none
    type(particle_tile), intent(inout)          :: tile
    type(momentumBin), allocatable, intent(in)  :: momentum_bins(:)
    integer                                     :: e_b, th_b, ph_b, p_ind, p, npart
    real                                        :: en, u, v, w, theta, phi

    do e_b = 0, n_energy_bins - 1
      do th_b = 0, momentum_bins(e_b)%n_theta_bins + 1
        do ph_b = 0, momentum_bins(e_b)%theta_bins(th_b)%n_phi_bins - 1
          npart = momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%npart
          do p_ind = 1, npart
            p = momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%indices(p_ind)
            ! for testing purposes
            #ifdef DEBUG
              if ((p .le. 0) .or. (p .gt. tile%npart_sp)) then
                print *, p, tile%npart_sp, npart
                call printBins(momentum_bins, e_b, th_b, ph_b)
                call throwError('Something went terribly wrong during the binning.')
              end if
            #endif
            tile%ind(p) = p_ind + 100 * ph_b + 100**2 * e_b
            #ifdef DEBUG
              u = tile%u(p)
              v = tile%v(p)
              w = tile%w(p)
              en = sqrt(u**2 + v**2 + w**2)
              u = u / en; v = v / en; w = w / en
              theta = asin(w)
              phi = atan2(v, u)
              if (phi .lt. 0) phi = phi + 2 * M_PI
              if ((en .ge. momentum_bins(e_b)%e_max) .or.&
                & (en .lt. momentum_bins(e_b)%e_min) .or.&
                & (theta .ge. momentum_bins(e_b)%theta_bins(th_b)%theta_max) .or.&
                & (theta .lt. momentum_bins(e_b)%theta_bins(th_b)%theta_min) .or.&
                & (phi .ge. momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%phi_max) .or.&
                & (phi .lt. momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%phi_min)) then
                print *, p, npart, en, theta, phi, u, v, w
                print *, momentum_bins(e_b)%e_min, momentum_bins(e_b)%e_max
                print *, momentum_bins(e_b)%theta_bins(th_b)%theta_min,&
                       & momentum_bins(e_b)%theta_bins(th_b)%theta_max
                print *, momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%phi_min,&
                       & momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%phi_max

                print *, ">>>>>>"
                call printBins(momentum_bins, e_b, th_b, ph_b)
                call throwError('Something went wrong during particle binning.')
              end if
            #endif
          end do
        end do
      end do
    end do
  end subroutine downsampleBinnedParticles

  subroutine printBins(momentum_bins, e_b, th_b, ph_b)
    implicit none
    type(momentumBin), allocatable, intent(in)  :: momentum_bins(:)
    integer, intent(in)                         :: e_b, th_b, ph_b

    print *, 'INDICES', e_b, th_b, ph_b
    print *, 'E', momentum_bins(e_b)%e_min,&
                & momentum_bins(e_b)%e_max,&
                & momentum_bins(e_b)%th0_bin,&
                & momentum_bins(e_b)%n_theta_bins
    print *, ""
    print *, 'TH', momentum_bins(e_b)%theta_bins(th_b)%theta_min,&
                 & momentum_bins(e_b)%theta_bins(th_b)%theta_max,&
                 & momentum_bins(e_b)%theta_bins(th_b)%theta_mid,&
                 & momentum_bins(e_b)%theta_bins(th_b)%n_phi_bins
    print *, ""
    print *, 'PH', momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%phi_min,&
                 & momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%phi_max,&
                 & momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%npart

  end subroutine printBins

#endif
end module m_particledownsampling
