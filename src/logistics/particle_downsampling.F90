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

  ! Auxiliary type for a group of merging particles
  type :: particleDwnGroup
    ! total momentum components and energy of the group
    real                  :: tot_px, tot_py, tot_pz, tot_en
    ! total weight of the group
    integer(kind=2)       :: tot_wei
    ! indices of particles on a tile contained in the group
    integer, allocatable  :: indices(:)
    ! size of the group
    integer               :: size
    ! average momentum unit vector of the bin
    real                  :: bin_px, bin_py, bin_pz
  end type particleDwnGroup

  integer           :: dwn_interval, dwn_maxweight

  !--- PRIVATE variables/functions -------------------------------!
  private :: downsampleParticles, downsampleOnTile,&
           & downsampleAllBins, downsampleBin,&
           & mergeParticlesInGroup
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
    call downsampleAllBins(momentum_bins, tile)
  end subroutine downsampleOnTile

  subroutine downsampleAllBins(momentum_bins, tile)
    implicit none
    type(particle_tile), intent(inout)            :: tile
    type(momentumBin), allocatable, intent(inout) :: momentum_bins(:)
    integer                                       :: e_b, th_b, ph_b, p_ind, p, npart
    #ifdef DEBUG
      real                                        :: en, u, v, w, theta, phi
    #endif

    ! loop through all the bins...
    ! ... in all 3 values (energy, theta, phi)
    do e_b = 0, n_energy_bins - 1
      do th_b = 0, momentum_bins(e_b)%n_theta_bins + 1
        do ph_b = 0, momentum_bins(e_b)%theta_bins(th_b)%n_phi_bins - 1
          npart = momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%npart
          #ifdef DEBUG
            do p_ind = 1, npart
              p = momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%indices(p_ind)
              if ((p .le. 0) .or. (p .gt. tile%npart_sp)) then
                call throwError('Wrong index in `downsampleAllBins()`.')
              end if
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
                call throwError('Wrong binning in `downsampleAllBins()`.')
              end if
              tile%ind(p) = ph_b + 100 * th_b + 100**2 * e_b
            end do
          #endif
          call downsampleBin(tile,&
                      & momentum_bins(e_b)%theta_bins(th_b)%theta_mid,&
                      & momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%phi_mid,&
                      & momentum_bins(e_b)%theta_bins(th_b)%phi_bins(ph_b)%indices, npart)
        end do
      end do
    end do
  end subroutine downsampleAllBins

  ! on each bin we are forming groups of particles...
  ! ... with cumulative weights less than `dwn_maxweight`...
  ! ... and sending them to merge into separate routine
  subroutine downsampleBin(tile, theta_mid, phi_mid, indices, npart)
    implicit none
    type(particle_tile), intent(inout)  :: tile
    integer, allocatable, intent(inout) :: indices(:)
    integer, intent(inout)              :: npart
    real, intent(in)                    :: theta_mid, phi_mid
    type(particleDwnGroup)              :: group
    integer       :: p_ind, p
    real          :: en

    allocate(group%indices(npart))
    group%indices(:) = -1
    group%size = 0
    group%tot_px = 0.0; group%tot_py = 0.0; group%tot_pz = 0.0
    group%tot_en = 0.0; group%tot_wei = 0_2

    group%bin_px = cos(theta_mid) * cos(phi_mid)
    group%bin_px = cos(theta_mid) * sin(phi_mid)
    group%bin_pz = sin(theta_mid)

    ! FIX: this is for photons only

    p_ind = 1
    do while (p_ind .le. npart)
      p = indices(p_ind)
      if (REAL(tile%weight(p)) .gt. sqrt(REAL(dwn_maxweight))) then
        ! particle to heavy to merge
        indices(p_ind) = indices(npart)
        npart = npart - 1
        cycle
      else
        group%tot_wei = group%tot_wei + tile%weight(p)
        group%tot_px = group%tot_px + tile%weight(p) * tile%u(p)
        group%tot_py = group%tot_py + tile%weight(p) * tile%v(p)
        group%tot_pz = group%tot_pz + tile%weight(p) * tile%w(p)
        en = sqrt(tile%u(p)**2 + tile%v(p)**2 + tile%w(p)**2)
        group%tot_en = group%tot_en + tile%weight(p) * en

        group%indices(group%size + 1) = p
        group%size = group%size + 1

        if ((group%tot_wei .ge. dwn_maxweight) .or. (p_ind .eq. npart)) then
          ! once there are enough particles in the group...
          ! ... send a group of these particles to merge...
          ! ... then reset the quantities
          if (group%size .gt. 2) then
            call mergeParticlesInGroup(group, tile)
          end if
          group%indices(:) = -1
          group%size = 0
          group%tot_px = 0.0; group%tot_py = 0.0; group%tot_pz = 0.0
          group%tot_en = 0.0; group%tot_wei = 0_2
        end if

        p_ind = p_ind + 1
      end if
    end do
  end subroutine downsampleBin

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !   The algorithm is adapted from Vranic et al. 2014 [1411.2248v1]
  subroutine mergeParticlesInGroup(group, tile)
    implicit none
    type(particleDwnGroup), intent(in)  :: group
    type(particle_tile), intent(inout)  :: tile
    integer(kind=2) :: wA, wB
    integer         :: p, p_ind, s
    real            :: pxA, pxB, pyA, pyB, pzA, pzB, rnd
    real            :: dxA, dyA, dzA, dxB, dyB, dzB
    integer(kind=2) :: xAi, yAi, zAi, xBi, yBi, zBi
    real            :: pA, pB, enA, enB, tot_p, cos_th, sin_th
    real            :: temp1_x, temp1_y, temp1_z, temp1
    real            :: temp2_x, temp2_y, temp2_z, temp2

    ! FIX: this is for photons only

    ! There are two directions in this problem...
    ! ... parallel to the total momentum...
    ! ... and perpendicular (which is arbitrary)

    ! We choose the perpendicular direction,..
    ! ... so that particles are created in the plane...
    ! ... containing both the total momentum...
    ! ... and the average momentum of the bin

    temp1_x = group%bin_pz * group%tot_py - group%bin_py * group%tot_pz
    temp1_y = -group%bin_pz * group%tot_px + group%bin_px * group%tot_pz
    temp1_z = group%bin_py * group%tot_px - group%bin_px * group%tot_py

    temp2_x = group%tot_pz * temp1_y - group%tot_py * temp1_z
    temp2_y = -group%tot_pz * temp1_x + group%tot_px * temp1_z
    temp2_z = group%tot_py * temp1_x - group%tot_px * temp1_y
    temp2 = sqrt(temp2_x**2 + temp2_y**2 + temp2_z**2)

    ! Now vector `temp2` is perpendicular to `tot_p`...
    ! ... but also lies in the plane containing the...
    ! ... average momentum of the bin

    tot_p = sqrt(group%tot_px**2 + group%tot_py**2 + group%tot_pz**2)

    ! new weights
    wA = INT(FLOOR(group%tot_wei / 2.0), 2)
    wB = INT(CEILING(group%tot_wei / 2.0), 2)

    ! new energies & momenta (magnitudes)
    enA = group%tot_en / (2.0 * wA)
    pA = enA
    enB = group%tot_en / (2.0 * wB)
    pB = enB

    ! direction between new particle momenta...
    ! ... and the total group momentum
    cos_th = tot_p / group%tot_en
    sin_th = sqrt(1.0 - cos_th**2)

    ! new momenta
    pxA = (group%tot_px / tot_p) * cos_th * pA +& ! parallel component
        & (temp2_x / temp2) * sin_th * pA   ! perp component
    pyA = (group%tot_py / tot_p) * cos_th * pA +&
        & (temp2_y / temp2) * sin_th * pA
    pzA = (group%tot_pz / tot_p) * cos_th * pA +&
        & (temp2_z / temp2) * sin_th * pA

    pxB = (group%tot_px / tot_p) * cos_th * pB -& ! parallel component
        & (temp2_x / temp2) * sin_th * pB   ! perp component
    pyB = (group%tot_py / tot_p) * cos_th * pB -&
        & (temp2_y / temp2) * sin_th * pB
    pzB = (group%tot_pz / tot_p) * cos_th * pB -&
        & (temp2_z / temp2) * sin_th * pB

    #ifdef DEBUG
      ! check weight conservation
      if (wA + wB .ne. group%tot_wei) then
        call throwError('Weight is not conserved in `mergeParticlesInGroup()`')
      end if
      if (pxA + pxB .ne. group%tot_px) then
        print *, pxA, pxB, pxA + pxB, group%tot_px, cos_th
        call throwError('Px is not conserved in `mergeParticlesInGroup()`')
      end if
      if (pyA + pyB .ne. group%tot_py) then
        call throwError('Px is not conserved in `mergeParticlesInGroup()`')
      end if
      if (pzA + pzB .ne. group%tot_pz) then
        call throwError('Px is not conserved in `mergeParticlesInGroup()`')
      end if
      if (enA + enB .ne. group%tot_en) then
        call throwError('Energy is not conserved in `mergeParticlesInGroup()`')
      end if
      if ((pxA**2 + pyA**2 + pzA**2 .ne. pA**2) .or.&
        & (pxB**2 + pyB**2 + pzB**2 .ne. pB**2)) then
        call throwError('Wrong momenta projections in `mergeParticlesInGroup()`')
      end if
    #endif

    ! take two random particles to position the new ones
    p_ind = INT((random(dseed) * group%size + 1))
    p = group%indices(p_ind)
    xAi = tile%xi(p); dxA = tile%dx(p)
    yAi = tile%yi(p); dyA = tile%dy(p)
    zAi = tile%zi(p); dzA = tile%dz(p)

    p = p_ind
    do while (p .eq. p_ind)
      p = INT((random(dseed) * group%size + 1))
    end do
    p = group%indices(p)
    xBi = tile%xi(p); dxB = tile%dx(p)
    yBi = tile%yi(p); dyB = tile%dy(p)
    zBi = tile%zi(p); dzB = tile%dz(p)

    ! "nullify" merged particles
    do p_ind = 1, group%size
      p = group%indices(p_ind)
      tile%proc(p) = -1
    end do

    ! inject new particles
    s = tile%spec
    call createParticle(s, xAi, yAi, zAi, dxA, dyA, dzA, pxA, pyA, pzA, weight=wA)
    call createParticle(s, xBi, yBi, zBi, dxB, dyB, dzB, pxB, pyB, pzB, weight=wB)
  end subroutine mergeParticlesInGroup
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
end module m_particledownsampling
