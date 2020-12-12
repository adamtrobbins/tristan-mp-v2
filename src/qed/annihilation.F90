#include "../defs.F90"

module m_annihilation
#ifdef PAIRANNIHILATION

  use m_globalnamespace
  use m_qednamespace
  use m_aux
  use m_errors
  use m_particlelogistics
  implicit none

  type, private :: particleID
    ! particle is uniquely identified by its:
    ! ... `s`, `ti, tj, tk` and `p`
    integer :: s, ti, tj, tk, p
    real    :: wei = -1
  end type particleID

  type, private :: particleGroup
    integer                       :: npart
    type(particleID), allocatable :: prtls(:)
  end type particleGroup

  type, private :: particlePair
    type(particleID) :: prtl1, prtl2
  end type particlePair

  !--- PRIVATE variables/functions -------------------------------!
  private :: pairAnnihilationWithGroups, pairAnnihilationWithGroups_mc,&
           & breakDownParticles, shuffleGroup, computeAnnihilationCrossSection,&
           & annihilatePairs

  !...............................................................!
contains
  subroutine pairAnnihilation()
    implicit none
    type(particleGroup), allocatable    :: lec_cell_bins(:,:,:), pos_cell_bins(:,:,:)
    type(particle_tile)                 :: electrons_tile, positrons_tile
    integer   :: ann_electrons(20), ann_positrons(20)
    integer   :: s, s_lec, s_pos, npart_lec, npart_pos, s0, si, nn
    integer   :: ti, tj, tk, nx_bin, ny_bin, nz_bin, pi, pj, pk
    integer   :: p, p_ind, s_ind, x1_tile, y1_tile, z1_tile

    ! find all the species that participate in the annihilation process
    s_lec = 0; s_pos = 0
    do s = 1, nspec
      if (species(s)%annihilation_sp) then
        if (species(s)%ch_sp .lt. 0) then
          ann_electrons(s_lec + 1) = s
          s_lec = s_lec + 1
        else if (species(s)%ch_sp .gt. 0) then
          ann_positrons(s_pos + 1) = s
          s_pos = s_pos + 1
        else
          call throwError("Only charged particles can participate in annihilation.")
        end if
      end if
    end do

    if ((s_lec * s_pos .eq. 0) .and. (s_lec + s_pos .ne. 0)) then
      call throwError("Pair annihilation requires at least 2 species with opposite charges to participate.")
    end if

    ! this assumes that all the species have the same tiles
    ! ... so picking just a random `s0`
    s0 = ann_electrons(1)
    ! loop over all tiles
    do ti = 1, species(s0)%tile_nx
      do tj = 1, species(s0)%tile_ny
        do tk = 1, species(s0)%tile_nz
          if ((Annihilation_sporadic) .and. (random(dseed) * Annihilation_interval .gt. 1.0)) then
            cycle
          end if

          ! number of cells on a tile
          nx_bin = species(s0)%prtl_tile(ti, tj, tk)%x2 - species(s0)%prtl_tile(ti, tj, tk)%x1
          ny_bin = species(s0)%prtl_tile(ti, tj, tk)%y2 - species(s0)%prtl_tile(ti, tj, tk)%y1
          nz_bin = species(s0)%prtl_tile(ti, tj, tk)%z2 - species(s0)%prtl_tile(ti, tj, tk)%z1

          ! count total # of electrons in the tile
          npart_lec = 0
          do si = 1, s_lec
            npart_lec = npart_lec + species(ann_electrons(si))%prtl_tile(ti, tj, tk)%npart_sp
          end do
          if (allocated(lec_cell_bins)) deallocate(lec_cell_bins)
          allocate(lec_cell_bins(nx_bin, ny_bin, nz_bin))

          if (npart_lec .eq. 0) then
            cycle
          end if

          ! count total # of positrons in the tile
          npart_pos = 0
          do si = 1, s_pos
            npart_pos = npart_pos + species(ann_positrons(si))%prtl_tile(ti, tj, tk)%npart_sp
          end do
          if (allocated(pos_cell_bins)) deallocate(pos_cell_bins)
          allocate(pos_cell_bins(nx_bin, ny_bin, nz_bin))

          if (npart_pos .eq. 0) then
            cycle
          end if

          ! initialize cell-based bins
          do pi = 1, nx_bin
            do pj = 1, ny_bin
              do pk = 1, nz_bin
                lec_cell_bins(pi, pj, pk)%npart = 0
                allocate(lec_cell_bins(pi, pj, pk)%prtls(npart_lec))
                pos_cell_bins(pi, pj, pk)%npart = 0
                allocate(pos_cell_bins(pi, pj, pk)%prtls(npart_pos))
              end do
            end do
          end do

          ! put particles into correct groups according to their cells
          x1_tile = species(s0)%prtl_tile(ti, tj, tk)%x1
          y1_tile = species(s0)%prtl_tile(ti, tj, tk)%y1
          z1_tile = species(s0)%prtl_tile(ti, tj, tk)%z1
          do si = 1, s_lec
            do p = 1, species(ann_electrons(si))%prtl_tile(ti, tj, tk)%npart_sp
              pi = species(ann_electrons(si))%prtl_tile(ti, tj, tk)%xi(p) - x1_tile + 1
              pj = species(ann_electrons(si))%prtl_tile(ti, tj, tk)%yi(p) - y1_tile + 1
              pk = species(ann_electrons(si))%prtl_tile(ti, tj, tk)%zi(p) - z1_tile + 1

              lec_cell_bins(pi, pj, pk)%npart = lec_cell_bins(pi, pj, pk)%npart + 1
              nn = lec_cell_bins(pi, pj, pk)%npart
              lec_cell_bins(pi, pj, pk)%prtls(nn)%s = ann_electrons(si)
              lec_cell_bins(pi, pj, pk)%prtls(nn)%ti = ti
              lec_cell_bins(pi, pj, pk)%prtls(nn)%tj = tj
              lec_cell_bins(pi, pj, pk)%prtls(nn)%tk = tk
              lec_cell_bins(pi, pj, pk)%prtls(nn)%p = p
            end do
          end do

          do si = 1, s_pos
            do p = 1, species(ann_positrons(si))%prtl_tile(ti, tj, tk)%npart_sp
              pi = species(ann_positrons(si))%prtl_tile(ti, tj, tk)%xi(p) - x1_tile + 1
              pj = species(ann_positrons(si))%prtl_tile(ti, tj, tk)%yi(p) - y1_tile + 1
              pk = species(ann_positrons(si))%prtl_tile(ti, tj, tk)%zi(p) - z1_tile + 1

              pos_cell_bins(pi, pj, pk)%npart = pos_cell_bins(pi, pj, pk)%npart + 1
              nn = pos_cell_bins(pi, pj, pk)%npart
              pos_cell_bins(pi, pj, pk)%prtls(nn)%s = ann_positrons(si)
              pos_cell_bins(pi, pj, pk)%prtls(nn)%ti = ti
              pos_cell_bins(pi, pj, pk)%prtls(nn)%tj = tj
              pos_cell_bins(pi, pj, pk)%prtls(nn)%tk = tk
              pos_cell_bins(pi, pj, pk)%prtls(nn)%p = p
            end do
          end do

          ! at this point positrons and electrons on a tile are distributed ...
          ! ... into groups based on their cells
          ! loop over all cells on a tile
          do pi = 1, nx_bin
            do pj = 1, ny_bin
              do pk = 1, nz_bin
                if ((lec_cell_bins(pi, pj, pk)%npart .ge. 1) .and. (pos_cell_bins(pi, pj, pk)%npart .ge. 1)) then
                  call pairAnnihilationWithGroups(lec_cell_bins(pi, pj, pk), pos_cell_bins(pi, pj, pk))
                end if
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine pairAnnihilation

  subroutine pairAnnihilationWithGroups(electron_group, positron_group)
    implicit none
    type(particleGroup), intent(in)     :: electron_group, positron_group
    if (Annihilation_algorithm .eq. 1) then
      call throwError('Annihilation algorithm currently only supports MC pairing.')
    else if (Annihilation_algorithm .eq. 2) then
      call pairAnnihilationWithGroups_mc(electron_group, positron_group)
    else
      call throwError('Unrecognized annihilation algorithm.')
    end if
  end subroutine pairAnnihilationWithGroups

  subroutine pairAnnihilationWithGroups_mc(electron_group, positron_group)
    implicit none
    type(particleGroup), intent(in)     :: electron_group, positron_group
    type(particleGroup)                 :: splitted_lec_group, splitted_pos_group
    real                                :: lec_weight, pos_weight
    type(particlePair), allocatable     :: ep_pairs(:)
    integer                             :: s1, s2, p1, p2
    integer                             :: num_couples, n, ti1, tj1, tk1, ti2, tj2, tk2
    real                                :: P_12, wei1, wei2, P_corr, wei_split_tot, rnd

    call breakDownParticles(electron_group, splitted_lec_group, lec_weight)
    call breakDownParticles(positron_group, splitted_pos_group, pos_weight)
    call shuffleGroup(splitted_lec_group)
    call shuffleGroup(splitted_pos_group)

    num_couples = min(splitted_lec_group%npart, splitted_pos_group%npart)
    allocate(ep_pairs(num_couples))

    ! couple particles
    wei_split_tot = 0.0
    do n = 1, num_couples
      ep_pairs(n)%prtl1 = splitted_lec_group%prtls(n)
      ep_pairs(n)%prtl2 = splitted_pos_group%prtls(n)
      wei_split_tot = wei_split_tot + min(ep_pairs(n)%prtl1%wei, ep_pairs(n)%prtl2%wei)
    end do

    ! make it independent of ppc0 & qed step:
    P_corr = (3.0 / 8.0) * QED_tau0 * REAL(Annihilation_interval) * CC / ppc0
    ! match with binary pairing:
    P_corr = P_corr * max(lec_weight, pos_weight)
    ! `min(wei_1, wei_2)` is what could be scattered in an ideal pairing world, ...
    ! ... `wei_split_tot` is what is actually available to scatter due to ...
    ! ... non-ideal pairing:
    P_corr = P_corr * min(lec_weight, pos_weight) / wei_split_tot

    do n = 1, num_couples
      ! compute P_12 for each pair of e+e-
      call computeAnnihilationCrossSection(ep_pairs(n), P_12)
      ! to match the optical depth with the binary pairing case:
      P_12 = P_12 * P_corr
      #ifdef DEBUG
        if ((P_12 .lt. 0.0) .or. (P_12 .gt. 1.0)) then
          print *, 'P_12 = ', P_12
          call throwError('Annihilation cross section P_12 out of bounds!')
        end if
      #else
        if ((P_12 .gt. 1.0)) then
          print '(1X,A,ES10.3,A)', 'Warning: Annihilation cross section P_12 = ', P_12, ' > 1 !!'
        endif
      #endif
      rnd = random(dseed)
      if (rnd .le. P_12) then
        ! pair produce
        call annihilatePairs(ep_pairs(n))
        ! schedule particles for deletion
        s1 = ep_pairs(n)%prtl1%s
        p1 = ep_pairs(n)%prtl1%p
        ti1 = ep_pairs(n)%prtl1%ti
        tj1 = ep_pairs(n)%prtl1%tj
        tk1 = ep_pairs(n)%prtl1%tk
        wei1 = ep_pairs(n)%prtl1%wei
        s2 = ep_pairs(n)%prtl2%s
        p2 = ep_pairs(n)%prtl2%p
        ti2 = ep_pairs(n)%prtl2%ti
        tj2 = ep_pairs(n)%prtl2%tj
        tk2 = ep_pairs(n)%prtl2%tk
        wei2 = ep_pairs(n)%prtl2%wei
        #ifdef DEBUG
          if ((wei1 .ne. 1) .or. (wei2 .ne. 1)) then
            call throwError("Only particles in groups with weight 1 can annihilate.")
          end if
        #endif
        species(s1)%prtl_tile(ti1, tj1, tk1)%weight(p1) = species(s1)%prtl_tile(ti1, tj1, tk1)%weight(p1) - wei1
        species(s2)%prtl_tile(ti2, tj2, tk2)%weight(p2) = species(s2)%prtl_tile(ti2, tj2, tk2)%weight(p2) - wei2
        if (species(s1)%prtl_tile(ti1, tj1, tk1)%weight(p1) .lt. 1e-6) then
          species(s1)%prtl_tile(ti1, tj1, tk1)%proc(p1) = -1
        end if
        if (species(s2)%prtl_tile(ti2, tj2, tk2)%weight(p2) .lt. 1e-6) then
          species(s2)%prtl_tile(ti2, tj2, tk2)%proc(p2) = -1
        end if
      end if
    end do
  end subroutine pairAnnihilationWithGroups_mc

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! . . . . Physics functions . . . .
  subroutine computeAnnihilationCrossSection(ep_pair, P_12)
    implicit none
    type(particlePair), intent(in)  :: ep_pair
    real, intent(out)               :: P_12
    integer                         :: s1, s2, p1, p2
    integer                         :: ti1, tj1, tk1, ti2, tj2, tk2
    real(kind=8)                    :: sigma_ann
    real(kind=8)                    :: lec_u, lec_v, lec_w, pos_u, pos_v, pos_w, lec_gamma, pos_gamma
    real(kind=8)                    :: COM_beta_u, COM_beta_v, COM_beta_w
    real(kind=8)                    :: COM_u, COM_v, COM_w, COM_gamma
    real(kind=8)                    :: gamma_2, gamma_2_sqr, u_2
    real                            :: wei1, wei2

    s1 = ep_pair%prtl1%s
    p1 = ep_pair%prtl1%p
    ti1 = ep_pair%prtl1%ti
    tj1 = ep_pair%prtl1%tj
    tk1 = ep_pair%prtl1%tk
    wei1 = ep_pair%prtl1%wei
    s2 = ep_pair%prtl2%s
    p2 = ep_pair%prtl2%p
    ti2 = ep_pair%prtl2%ti
    tj2 = ep_pair%prtl2%tj
    tk2 = ep_pair%prtl2%tk
    wei2 = ep_pair%prtl2%wei

    #ifdef DEBUG
      if ((wei1 .ne. wei2) .or. (wei1 .ne. 1.0)) then
        call throwError('Unequal weights in `computeAnnihilationCrossSection`: '//STR(wei1)//':'//STR(wei2))
      end if
    #endif

    lec_u = REAL(species(s1)%prtl_tile(ti1, tj1, tk1)%u(p1), 8)
    lec_v = REAL(species(s1)%prtl_tile(ti1, tj1, tk1)%v(p1), 8)
    lec_w = REAL(species(s1)%prtl_tile(ti1, tj1, tk1)%w(p1), 8)
    lec_gamma = sqrt(1d0 + lec_u**2 + lec_v**2 + lec_w**2)
    pos_u = REAL(species(s2)%prtl_tile(ti2, tj2, tk2)%u(p2), 8)
    pos_v = REAL(species(s2)%prtl_tile(ti2, tj2, tk2)%v(p2), 8)
    pos_w = REAL(species(s2)%prtl_tile(ti2, tj2, tk2)%w(p2), 8)
    pos_gamma = sqrt(1d0 + pos_u**2 + pos_v**2 + pos_w**2)

    gamma_2 = lec_gamma * pos_gamma - (lec_u * pos_u + lec_v * pos_v + lec_w * pos_w)
    gamma_2_sqr = gamma_2 * gamma_2

    ! using assymptotic relations for `gamma_2 >> 1` and `gamma_2 ~ 1` ...
    ! ... with an error of <0.01%
    if (gamma_2 .lt. 1.01) then
      sigma_ann = 1d0 / sqrt(1d0 - 1d0 / gamma_2_sqr)
    else if (gamma_2 .gt. 100) then
      sigma_ann = ((log(2d0 * gamma_2) - 1d0) / gamma_2) + ((3d0 * log(2d0 * gamma_2) - 2d0) / gamma_2_sqr)
    else
      u_2 = sqrt(gamma_2_sqr - 1d0)
      sigma_ann = ((gamma_2_sqr + 4d0 * gamma_2 + 1d0) * log(gamma_2 + u_2) / (gamma_2_sqr - 1d0) -&
                   (gamma_2 + 3d0) / u_2) / (gamma_2 + 1d0)
    end if

    ! take into account the relative velocity
    P_12 = sigma_ann * REAL(sqrt(gamma_2**2 - 1d0) / (lec_gamma * pos_gamma))
  end subroutine computeAnnihilationCrossSection

  subroutine annihilatePairs(ep_pair)
    implicit none
    type(particlePair), intent(in)    :: ep_pair
    integer                         :: s1, s2, p1, p2
    integer                         :: ti1, tj1, tk1, ti2, tj2, tk2
    real(kind=8)                    :: lec_u, lec_v, lec_w, pos_u, pos_v, pos_w, lec_gamma, pos_gamma
    real                            :: wei1, wei2

    s1 = ep_pair%prtl1%s
    p1 = ep_pair%prtl1%p
    ti1 = ep_pair%prtl1%ti
    tj1 = ep_pair%prtl1%tj
    tk1 = ep_pair%prtl1%tk
    wei1 = ep_pair%prtl1%wei
    s2 = ep_pair%prtl2%s
    p2 = ep_pair%prtl2%p
    ti2 = ep_pair%prtl2%ti
    tj2 = ep_pair%prtl2%tj
    tk2 = ep_pair%prtl2%tk
    wei2 = ep_pair%prtl2%wei

    #ifdef DEBUG
      if ((wei1 .ne. wei2) .or. (wei1 .ne. 1.0)) then
        call throwError('Unequal weights in `computeAnnihilationCrossSection`: '//STR(wei1)//':'//STR(wei2))
      end if
    #endif

    lec_u = REAL(species(s1)%prtl_tile(ti1, tj1, tk1)%u(p1), 8)
    lec_v = REAL(species(s1)%prtl_tile(ti1, tj1, tk1)%v(p1), 8)
    lec_w = REAL(species(s1)%prtl_tile(ti1, tj1, tk1)%w(p1), 8)
    lec_gamma = sqrt(1d0 + lec_u**2 + lec_v**2 + lec_w**2)
    pos_u = REAL(species(s2)%prtl_tile(ti2, tj2, tk2)%u(p2), 8)
    pos_v = REAL(species(s2)%prtl_tile(ti2, tj2, tk2)%v(p2), 8)
    pos_w = REAL(species(s2)%prtl_tile(ti2, tj2, tk2)%w(p2), 8)
    pos_gamma = sqrt(1d0 + pos_u**2 + pos_v**2 + pos_w**2)


  end subroutine annihilatePairs

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! . . . . Technical functions . . . .
  subroutine LorentzBoost(beta_frame_x, beta_frame_y, beta_frame_z,&
                        & beta_frame_sq, gamma_frame,&
                        & old_k_0,&
                        & old_k_x, old_k_y, old_k_z,&
                        & new_k_x, new_k_y, new_k_z)
    implicit none
    real(kind=8), intent(in)  :: beta_frame_x, beta_frame_y, beta_frame_z
    real(kind=8), intent(in)  :: beta_frame_sq, gamma_frame
    real(kind=8), intent(in)  :: old_k_x, old_k_y, old_k_z, old_k_0
    real(kind=8), intent(out) :: new_k_x, new_k_y, new_k_z
    real(kind=8)              :: gamma_frame_m1

    if (beta_frame_sq .gt. 0.0d0) then
      gamma_frame_m1 = gamma_frame - 1.0d0
      new_k_x = -beta_frame_x * gamma_frame * old_k_0 +&
              & old_k_x * (1.0d0 + (beta_frame_x**2 / beta_frame_sq) * gamma_frame_m1) +&
              & old_k_y * (beta_frame_x * beta_frame_y / beta_frame_sq) * gamma_frame_m1 +&
              & old_k_z * (beta_frame_x * beta_frame_z / beta_frame_sq) * gamma_frame_m1
      new_k_y = -beta_frame_y * gamma_frame * old_k_0 +&
              & old_k_x * (beta_frame_y * beta_frame_x / beta_frame_sq) * gamma_frame_m1 +&
              & old_k_y * (1.0d0 + (beta_frame_y**2 / beta_frame_sq) * gamma_frame_m1) +&
              & old_k_z * (beta_frame_y * beta_frame_z / beta_frame_sq) * gamma_frame_m1
      new_k_z = -beta_frame_z * gamma_frame * old_k_0 +&
              & old_k_x * (beta_frame_z * beta_frame_x / beta_frame_sq) * gamma_frame_m1 +&
              & old_k_y * (beta_frame_z * beta_frame_y / beta_frame_sq) * gamma_frame_m1 +&
              & old_k_z * (1.0d0 + (beta_frame_z**2 / beta_frame_sq) * gamma_frame_m1)
    else
      new_k_x = old_k_x; new_k_y = old_k_y; new_k_z = old_k_z
    end if
  end subroutine LorentzBoost

  subroutine breakDownParticles(group, set, set_weight)
    implicit none
    type(particleGroup)                           :: group
    real, intent(out)                             :: set_weight
    type(particleGroup), intent(out)              :: set
    integer                                       :: set_size_, set_size
    type(particleGroup)                           :: set_
    integer                                       :: s, si, i, p, q, pi, ti, tj, tk
    real                                          :: wei

    ! computing number of particles in the set
    set_size_ = 0
    do pi = 1, group%npart
      s = group%prtls(pi)%s
      ti = group%prtls(pi)%ti
      tj = group%prtls(pi)%tj
      tk = group%prtls(pi)%tk
      p = group%prtls(pi)%p
      set_size_ = set_size_ + CEILING(species(s)%prtl_tile(ti, tj, tk)%weight(p))
    end do
    set_%npart = 0
    allocate(set_%prtls(set_size_))
    ! assigning particles in the set
    i = 1
    set_weight = 0.0
    do pi = 1, group%npart
      s = group%prtls(pi)%s
      ti = group%prtls(pi)%ti
      tj = group%prtls(pi)%tj
      tk = group%prtls(pi)%tk
      p = group%prtls(pi)%p
      wei = species(s)%prtl_tile(ti, tj, tk)%weight(p)
      ! distribute the weight so that all weights are 1
      set_weight = set_weight + FLOOR(wei)
      do q = 1, FLOOR(wei)
        set_%prtls(i)%s = s
        set_%prtls(i)%ti = group%prtls(pi)%ti
        set_%prtls(i)%tj = group%prtls(pi)%tj
        set_%prtls(i)%tk = group%prtls(pi)%tk
        set_%prtls(i)%p = p
        set_%prtls(i)%wei = 1.0
        set_%npart = set_%npart + 1
        i = i + 1
      end do
    end do
    set_size = i - 1
    allocate(set%prtls(set_size))
    set%npart = set_size
    set%prtls(1 : set_size) = set_%prtls(1 : set_size)
  end subroutine breakDownParticles

  ! Knuth's algorithm to randomly shuffle a group
  subroutine shuffleGroup(group)
    implicit none
    type(particleGroup), intent(inout)  :: group
    type(particleID)                    :: temp
    integer                             :: i, j
    do i = 1, group%npart - 1
      j = randomInt(dseed, i, group%npart + 1)
      temp = group%prtls(i)
      group%prtls(i) = group%prtls(j)
      group%prtls(j) = temp
    end do
  end subroutine shuffleGroup
#endif
end module m_annihilation
