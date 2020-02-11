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

  real              :: dwn_interval
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
    bin_limit = int(sqrt(2 * n_energy_bins * n_angular_bins**2))

    do s = 1, nspec
      if (species(s)%dwn_sp) then
        do ti = 1, species(s)%tile_nx
          do tj = 1, species(s)%tile_ny
            do tk = 1, species(s)%tile_nz
              if (species(s)%prtl_tile(ti, tj, tk)%npart_sp .gt. bin_limit) then
                call downsampleOnTile(s, ti, tj, tk)
              end if
            end do
          end do
        end do
      end if
    end do

  end subroutine downsampleParticles

  subroutine downsampleOnTile(s, ti, tj, tk)
    implicit none
    integer, intent(in) :: s, ti, tj, tk
    ! particles are now binned into `energy_bins`
    call binParticlesOnTile(species(s)%prtl_tile(ti, tj, tk))
    call downsampleBinnedParticles()
    call deinitializeEnergyBins()
  end subroutine downsampleOnTile

  subroutine downsampleBinnedParticles(s, ti, tj, tk)
    implicit none
    integer, intent(in) :: s, ti, tj, tk
    integer             :: e_b, th_b, ph_b, p_ind, p, npart
    do e_b = 0, n_energy_bins - 1
      do th_b = 0, energy_bins(e_b)%n_theta_bins
        do ph_b = 0, energy_bins(e_b)%theta_bins(th_b)%n_phi_bins
          npart = energy_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%npart
          do p_ind = 1, npart
            p = energy_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%indices(p_ind)
            ! for testing purposes
            species(s)%prtl_tile(ti, tj, tk)%ind = p_ind + 100 * ph_b + 100**2 * e_b
          end do
        end do
      end do
    end do
  end subroutine downsampleBinnedParticles

#endif
end module m_particledownsampling
