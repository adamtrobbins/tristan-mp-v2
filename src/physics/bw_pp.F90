#include "../defs.F90"

module m_bwpairproduction
#ifdef BWPAIRPRODUCTION

  use m_globalnamespace
  use m_aux
  use m_errors
  use m_bincoupling
  use m_particlelogistics
  implicit none

  ! BW parameters
  real    :: BW_tau
  integer :: BW_interval, BW_electron_sp, BW_positron_sp
  integer :: BW_algorithm

  !--- PRIVATE variables/functions -------------------------------!
  private :: bwOnTile_bin, bwOnTile_mc, PPfromTwoPhotons
  private :: generateRandomThetaBW, dSdO_BW
  !...............................................................!
contains
  subroutine bwPairProduction()
    implicit none
    integer   :: bw_species_1(20), bw_species_2(20)
    integer   :: s, si_1, si_2
    integer   :: ti, tj, tk

    ! find all the species that participate in BW process
    si_1 = 0; si_2 = 0
    do s = 1, nspec
      if (species(s)%bw_sp .eq. 1) then
        bw_species_1(si_1 + 1) = s
        si_1 = si_1 + 1
      else if (species(s)%bw_sp .eq. 2) then
        bw_species_2(si_2 + 1) = s
        si_2 = si_2 + 1
      end if
    end do

    if ((si_1 .ne. 0) .and. (si_2 .ne. 0)) then
      ! if there are two groups of BW photons

      ! assuming all tiles are equal
      s = bw_species_1(1)
      ! loop through all the tiles
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            if (BW_algorithm .eq. 1) then
              call bwOnTile_bin(ti, tj, tk,&
                              & bw_species_1(1 : si_1), si_1,&
                              & bw_species_2(1 : si_2), si_2)
            else if (BW_algorithm .eq. 2) then
              call bwOnTile_mc(ti, tj, tk,&
                             & bw_species_1(1 : si_1), si_1,&
                             & bw_species_2(1 : si_2), si_2)
            end if
          end do
        end do
      end do
    else if ((si_1 .ne. 0) .or. (si_2 .ne. 0)) then
      ! if only one group of BW photons

      if ((si_1 .eq. 0) .and. (si_2 .ne. 0)) then
        ! if only group #2
        bw_species_1(1 : si_2) = bw_species_2(1 : si_2)
        si_1 = si_2
      end if

      ! assuming all tiles are equal
      s = bw_species_1(1)
      ! loop through all the tiles
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            if (BW_algorithm .eq. 1) then
              call bwOnTile_bin(ti, tj, tk,&
                              & bw_species_1(1 : si_1), si_1,&
                              & n_sp_2 = 0)
            else if (BW_algorithm .eq. 2) then
              call bwOnTile_mc(ti, tj, tk,&
                             & bw_species_1(1 : si_1), si_1,&
                             & n_sp_2 = 0)
            end if
          end do
        end do
      end do
    end if
  end subroutine bwPairProduction

  subroutine bwOnTile_bin(ti, tj, tk,&
                        & sp_arr_1, n_sp_1,&
                        & sp_arr_2, n_sp_2)
    implicit none
    integer, intent(in)                           :: ti, tj, tk
    integer, intent(in)                           :: n_sp_1 ! # of species in set
    integer, intent(in)                           :: sp_arr_1(n_sp_1)
    integer, intent(in)                           :: n_sp_2 ! # of species in set
    integer, optional, intent(in)                 :: sp_arr_2(n_sp_2)
    type(spec_ind_pair), allocatable              :: set(:), set_1(:), set_2(:)
    integer                                       :: set_p_1, set_p_2, s1, s2, p1, p2
    integer                                       :: set_size, set_size_1, set_size_2
    type(couple)                                  :: pair_of_photons
    real                                          :: rnd, P_12, P_1
    logical                                       :: thresholdQ

    if (n_sp_2 .ne. 0) then
      ! two separate groups of photons interacting with each other
      if (random(dseed) .gt. 0.5) then ! shuffle sets for more randomness
        call prtlToSet(ti, tj, tk, sp_arr_1, n_sp_1, set_1, set_size_1)
        call prtlToSet(ti, tj, tk, sp_arr_2, n_sp_2, set_2, set_size_2)
      else
        call prtlToSet(ti, tj, tk, sp_arr_2, n_sp_2, set_1, set_size_1)
        call prtlToSet(ti, tj, tk, sp_arr_1, n_sp_1, set_2, set_size_2)
      end if

      call shuffleSet(set_1, set_size_1)
      call shuffleSet(set_2, set_size_2)

      do set_p_1 = 1, set_size_1
        s1 = set_1(set_p_1)%spec
        p1 = set_1(set_p_1)%index
        ! check if particle is scheduled for deletion
        if (species(s1)%prtl_tile(ti, tj, tk)%proc(p1) .lt. 0) cycle
        rnd = random(dseed)
        P_1 = 0.0
        do set_p_2 = 1, set_size_2
          s2 = set_2(set_p_2)%spec
          p2 = set_2(set_p_2)%index
          ! check if particle is scheduled for deletion
          if (species(s2)%prtl_tile(ti, tj, tk)%proc(p2) .lt. 0) cycle
          pair_of_photons%part_1 = set_1(set_p_1)
          pair_of_photons%part_2 = set_2(set_p_2)
          ! compute `P_12`
          call computeBWCrossSection(ti, tj, tk, pair_of_photons,&
                                   & P_12, thresholdQ)
          P_1 = P_1 + P_12
          if ((P_1 .gt. rnd) .and. (thresholdQ)) then
            ! pair produce
            call PPfromTwoPhotons(ti, tj, tk, pair_of_photons)
            ! schedule particles for deletion
            species(s1)%prtl_tile(ti, tj, tk)%proc(p1) = -1
            species(s2)%prtl_tile(ti, tj, tk)%proc(p2) = -1
            exit
          end if
        end do
      end do
      if (allocated(set_1)) deallocate(set_1)
      if (allocated(set_2)) deallocate(set_2)
    else
      ! one group of photons interacting with each other
      call prtlToSet(ti, tj, tk, sp_arr_1, n_sp_1, set, set_size)
      call shuffleSet(set, set_size)

      do set_p_1 = 1, set_size
        s1 = set(set_p_1)%spec
        p1 = set(set_p_1)%index
        ! check if particle is scheduled for deletion
        if (species(s1)%prtl_tile(ti, tj, tk)%proc(p1) .lt. 0) cycle
        rnd = random(dseed)
        P_1 = 0.0
        do set_p_2 = set_p_1 + 1, set_size
          s2 = set(set_p_2)%spec
          p2 = set(set_p_2)%index
          ! check if particle is scheduled for deletion
          if (species(s2)%prtl_tile(ti, tj, tk)%proc(p2) .lt. 0) cycle
          pair_of_photons%part_1 = set(set_p_1)
          pair_of_photons%part_2 = set(set_p_2)
          ! compute `P_12`
          call computeBWCrossSection(ti, tj, tk, pair_of_photons,&
                                   & P_12, thresholdQ)
          P_1 = P_1 + P_12
          if ((P_1 .gt. rnd) .and. (thresholdQ)) then
            ! pair produce
            call PPfromTwoPhotons(ti, tj, tk, pair_of_photons)
            ! schedule particles for deletion
            species(s1)%prtl_tile(ti, tj, tk)%proc(p1) = -1
            species(s2)%prtl_tile(ti, tj, tk)%proc(p2) = -1
            exit
          end if
        end do
      end do
      if (allocated(set)) deallocate(set)
    end if
  end subroutine

  subroutine bwOnTile_mc(ti, tj, tk,&
                       & sp_arr_1, n_sp_1,&
                       & sp_arr_2, n_sp_2)
    implicit none
    integer, intent(in)           :: ti, tj, tk
    integer, intent(in)           :: n_sp_1 ! # of species in set
    integer, intent(in)           :: sp_arr_1(n_sp_1)
    integer, intent(in)           :: n_sp_2 ! # of species in set
    integer, optional, intent(in) :: sp_arr_2(n_sp_2)
    type(couple), allocatable     :: pairs_of_photons(:)
    integer                       :: num_pairs, ph, s1, s2, p1, p2
    real                          :: rnd, P_12
    logical                       :: thresholdQ
    integer                       :: num_1, num_2

    if (n_sp_2 .ne. 0) then
      ! two separate BW groups
      ! shuffle sets for more randomness
      if (random(dseed) .gt. 0.5) then
        call coupleParticlesOnTile(ti, tj, tk, sp_arr_1, n_sp_1, sp_arr_2, n_sp_2,&
                                 & pairs_of_photons, num_pairs, num_1, num_2)
      else
        call coupleParticlesOnTile(ti, tj, tk, sp_arr_2, n_sp_2, sp_arr_1, n_sp_1,&
                                 & pairs_of_photons, num_pairs, num_1, num_2)
      end if
    else
      ! one BW group
      call coupleParticlesOnTile(ti, tj, tk, sp_arr_1, n_sp_1, sp_arr_1, n_sp_1,&
                               & pairs_of_photons, num_pairs, num_1)
      num_2 = num_1
    end if

    do ph = 1, num_pairs
      ! compute P_12 for each pair of photons `pairs_of_photons(ph)`
      call computeBWCrossSection(ti, tj, tk, pairs_of_photons(ph),&
                               & P_12, thresholdQ)
      ! to match the optical depth with the binary pairing case:
      P_12 = P_12 * num_2
      rnd = random(dseed)
      if ((rnd .le. P_12) .and. (thresholdQ)) then
        ! pair produce
        call PPfromTwoPhotons(ti, tj, tk, pairs_of_photons(ph))
        ! schedule particles for deletion
        s1 = pairs_of_photons(ph)%part_1%spec
        p1 = pairs_of_photons(ph)%part_1%index
        s2 = pairs_of_photons(ph)%part_2%spec
        p2 = pairs_of_photons(ph)%part_2%index
        species(s1)%prtl_tile(ti, tj, tk)%proc(p1) = -1
        species(s2)%prtl_tile(ti, tj, tk)%proc(p2) = -1
      end if
    end do
  end subroutine

  subroutine computeBWCrossSection(ti, tj, tk, pair_of_photons,&
                                 & P_12, thresholdQ)
    implicit none
    integer, intent(in)      :: ti, tj, tk
    type(couple), intent(in) :: pair_of_photons
    logical, intent(out)     :: thresholdQ
    real, intent(out)        :: P_12
    real(kind=8)             :: k1_x, k1_y, k1_z
    real(kind=8)             :: k2_x, k2_y, k2_z
    real(kind=8)             :: cosphi, S, beta, beta2, fs
    real(kind=8)             :: E1, E2
    real(kind=8)             :: ph_u1, ph_v1, ph_w1, ph_u2, ph_v2, ph_w2
    integer                  :: s1, s2, p1, p2

    ! "extract" photons
    s1 = pair_of_photons%part_1%spec
    p1 = pair_of_photons%part_1%index
    s2 = pair_of_photons%part_2%spec
    p2 = pair_of_photons%part_2%index

    ph_u1 = REAL(species(s1)%prtl_tile(ti, tj, tk)%u(p1), 8)
    ph_v1 = REAL(species(s1)%prtl_tile(ti, tj, tk)%v(p1), 8)
    ph_w1 = REAL(species(s1)%prtl_tile(ti, tj, tk)%w(p1), 8)
    ph_u2 = REAL(species(s2)%prtl_tile(ti, tj, tk)%u(p2), 8)
    ph_v2 = REAL(species(s2)%prtl_tile(ti, tj, tk)%v(p2), 8)
    ph_w2 = REAL(species(s2)%prtl_tile(ti, tj, tk)%w(p2), 8)

    E1 = sqrt((ph_u1)**2 + (ph_v1)**2 + (ph_w1)**2)
    E2 = sqrt((ph_u2)**2 + (ph_v2)**2 + (ph_w2)**2)

    if (E1 * E2 .lt. 1.0) then
      thresholdQ = .false.
    else
      ! photon k-vectors
      k1_x = ph_u1 / E1; k1_y = ph_v1 / E1; k1_z = ph_w1 / E1
      k2_x = ph_u2 / E2; k2_y = ph_v2 / E2; k2_z = ph_w2 / E2
      cosphi = k1_x * k2_x + k1_y * k2_y + k1_z * k2_z
      S = E1 * E2 * (1.0 - cosphi) * 0.5
      thresholdQ = (S .gt. 1.0000001)
    end if
    if (thresholdQ) then
      beta2 = 1.0 - 1.0 / S
      beta = sqrt(beta2)
      fs = (1.0 - beta2) *&
         & (-2.0 * beta * (2.0 - beta2) + (3.0 - beta2**2) *&
         & log((1.0 + beta) / (1.0 - beta)))
      P_12 = REAL(BW_interval) * BW_tau * fs
    else
      P_12 = 0.0
    end if
  end subroutine computeBWCrossSection

  subroutine LorentzBoost(beta_frame_x, beta_frame_y, beta_frame_z,&
                        & beta_frame_sq, gamma_frame,&
                        & old_k_x, old_k_y, old_k_z,&
                        & new_k_x, new_k_y, new_k_z)
    implicit none
    real(kind=8), intent(in)  :: beta_frame_x, beta_frame_y, beta_frame_z
    real(kind=8), intent(in)  :: beta_frame_sq, gamma_frame
    real(kind=8), intent(in)  :: old_k_x, old_k_y, old_k_z
    real(kind=8), intent(out) :: new_k_x, new_k_y, new_k_z
    real(kind=8)              :: gamma_frame_m1

    if (beta_frame_sq .gt. 0.0) then
      gamma_frame_m1 = gamma_frame - 1.0
      new_k_x = -beta_frame_x * gamma_frame +&
              & old_k_x * (1.0 + (beta_frame_x**2 / beta_frame_sq) * gamma_frame_m1) +&
              & old_k_y * (beta_frame_x * beta_frame_y / beta_frame_sq) * gamma_frame_m1 +&
              & old_k_z * (beta_frame_x * beta_frame_z / beta_frame_sq) * gamma_frame_m1
      new_k_y = -beta_frame_y * gamma_frame +&
              & old_k_x * (beta_frame_y * beta_frame_x / beta_frame_sq) * gamma_frame_m1 +&
              & old_k_y * (1.0 + (beta_frame_y**2 / beta_frame_sq) * gamma_frame_m1) +&
              & old_k_z * (beta_frame_y * beta_frame_z / beta_frame_sq) * gamma_frame_m1
      new_k_z = -beta_frame_z * gamma_frame +&
              & old_k_x * (beta_frame_z * beta_frame_x / beta_frame_sq) * gamma_frame_m1 +&
              & old_k_y * (beta_frame_z * beta_frame_y / beta_frame_sq) * gamma_frame_m1 +&
              & old_k_z * (1.0 + (beta_frame_z**2 / beta_frame_sq) * gamma_frame_m1)
    else
      new_k_x = old_k_x; new_k_y = old_k_y; new_k_z = old_k_z
    end if

  end subroutine LorentzBoost

  subroutine PPfromTwoPhotons(ti, tj, tk, pair_of_photons)
    implicit none
    integer, intent(in)      :: ti, tj, tk
    type(couple), intent(in) :: pair_of_photons
    integer                  :: s1, s2, p1, p2
    real                     :: x_new, y_new, z_new, dx_new, dy_new, dz_new, rnd
    integer                  :: tile_x1, tile_x2, tile_y1, tile_y2, tile_z1, tile_z2
    integer(kind=2)          :: xi_new, yi_new, zi_new

    real(kind=8)             :: ph1_u, ph1_v, ph1_w, ph2_u, ph2_v, ph2_w, eps1, eps2
    real(kind=8)             :: k1_x, k1_y, k1_z, k2_x, k2_y, k2_z, cos_phi, SS
    real(kind=8)             :: gamma_prtl_CM, beta_prtl_CM
    real(kind=8)             :: prtl1_CM_u, prtl1_CM_v, prtl1_CM_w
    real(kind=8)             :: beta_CM_x, beta_CM_y, beta_CM_z, beta_CM_sq, gamma_CM
    real(kind=8)             :: k1_CM_x, k1_CM_y, k1_CM_z, k1_CM
    real(kind=8)             :: a_CM_x, a_CM_y, a_CM_z, a_CM
    real(kind=8)             :: b_CM_x, b_CM_y, b_CM_z, b_CM
    real(kind=8)             :: rand_theta_CM, cos_rand_theta_CM, sin_rand_theta_CM
    real(kind=8)             :: rand_phi_CM, cos_rand_phi_CM, sin_rand_phi_CM
    real(kind=8)             :: p1_CM_x, p1_CM_y, p1_CM_z
    real(kind=8)             :: prtl1_u, prtl1_v, prtl1_w
    real(kind=8)             :: prtl2_u, prtl2_v, prtl2_w

    ! preliminary stuff
    s1 = pair_of_photons%part_1%spec
    p1 = pair_of_photons%part_1%index
    s2 = pair_of_photons%part_2%spec
    p2 = pair_of_photons%part_2%index

    tile_x1 = species(s1)%prtl_tile(ti, tj, tk)%x1
    tile_x2 = species(s1)%prtl_tile(ti, tj, tk)%x2
    tile_y1 = species(s1)%prtl_tile(ti, tj, tk)%y1
    tile_y2 = species(s1)%prtl_tile(ti, tj, tk)%y2
    tile_z1 = species(s1)%prtl_tile(ti, tj, tk)%z1
    tile_z2 = species(s1)%prtl_tile(ti, tj, tk)%z2

    ph1_u = REAL(species(s1)%prtl_tile(ti, tj, tk)%u(p1), 8)
    ph1_v = REAL(species(s1)%prtl_tile(ti, tj, tk)%v(p1), 8)
    ph1_w = REAL(species(s1)%prtl_tile(ti, tj, tk)%w(p1), 8)
    ph2_u = REAL(species(s2)%prtl_tile(ti, tj, tk)%u(p2), 8)
    ph2_v = REAL(species(s2)%prtl_tile(ti, tj, tk)%v(p2), 8)
    ph2_w = REAL(species(s2)%prtl_tile(ti, tj, tk)%w(p2), 8)

    ! k-vectors in lab frame
    eps1 = sqrt(ph1_u**2 + ph1_v**2 + ph1_w**2)
    eps2 = sqrt(ph2_u**2 + ph2_v**2 + ph2_w**2)
    k1_x = ph1_u / eps1; k1_y = ph1_v / eps1; k1_z = ph1_w / eps1
    k2_x = ph2_u / eps2; k2_y = ph2_v / eps1; k2_z = ph2_w / eps2

    ! angle between photons in lab frame
    cos_phi = k1_x * k2_x + k1_y * k2_y + k1_z * k2_z
    ! `S` parameter
    SS = eps1 * eps2 * (1.0 - cos_phi) / 2.0
    ! Lorentz-factor of electron/positron in CoM frame
    gamma_prtl_CM = sqrt(SS)
    beta_prtl_CM = sqrt(1.0 - 1.0 / SS)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! 3-velocity of the CoM frame
    beta_CM_x = (eps1 * k1_x + eps2 * k2_x) / (eps1 + eps2)
    beta_CM_y = (eps1 * k1_y + eps2 * k2_y) / (eps1 + eps2)
    beta_CM_z = (eps1 * k1_z + eps2 * k2_z) / (eps1 + eps2)
    beta_CM_sq = beta_CM_x**2 + beta_CM_y**2 + beta_CM_z**2

    #ifdef DEBUG
      if (beta_CM_sq .ge. 1.0) then
        call throwError('`beta_CM_sq` >= 1 in `PPfromTwoPhotons()`.')
      end if
    #endif

    gamma_CM = 1.0 / sqrt(1.0 - beta_CM_sq)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Lorentz boost `k1` from lab to CoM frame
    call LorentzBoost(beta_CM_x, beta_CM_y, beta_CM_z,&
                    & beta_CM_sq, gamma_CM,&
                    & k1_x, k1_y, k1_z,&
                    & k1_CM_x, k1_CM_y, k1_CM_z)

    k1_CM = sqrt(k1_CM_x**2 + k1_CM_y**2 + k1_CM_z**2)
    k1_CM_x = k1_CM_x / k1_CM
    k1_CM_y = k1_CM_y / k1_CM
    k1_CM_z = k1_CM_z / k1_CM

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define a basis in in CoM frame: `k1_CM`, `a_CM` and `b_CM`
    if (k1_CM_x .ne. 0.0) then
      a_CM_x = -k1_CM_y / k1_CM_x; a_CM_y = 1.0; a_CM_z = 0.0
      a_CM = sqrt(a_CM_x**2 + a_CM_y**2)
      a_CM_x = a_CM_x / a_CM; a_CM_y = a_CM_y / a_CM
    else
      a_CM_x = 1.0; a_CM_y = 0.0; a_CM_z = 0.0
    end if
    b_CM_x = a_CM_z * k1_CM_y - a_CM_y * k1_CM_z
    b_CM_y = -a_CM_z * k1_CM_x + a_CM_x * k1_CM_z
    b_CM_z = a_CM_y * k1_CM_x - a_CM_x * k1_CM_y
    b_CM = sqrt(b_CM_x**2 + b_CM_y**2 + b_CM_z**2)
    b_CM_x = b_CM_x / b_CM; b_CM_y = b_CM_y / b_CM; b_CM_z = b_CM_z / b_CM

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Generate random vector in the CoM frame ...
    ! ... respecting the differential cross section ...
    ! ... at angle `theta` w.r.t. `k1_CM`
    call generateRandomThetaBW(rand_theta_CM, SS)
    rand_phi_CM = 2.0 * M_PI * random(dseed)
    cos_rand_theta_CM = cos(rand_theta_CM)
    sin_rand_theta_CM = sin(rand_theta_CM)
    cos_rand_phi_CM = cos(rand_phi_CM)
    sin_rand_phi_CM = sin(rand_phi_CM)

    p1_CM_x = k1_CM_x * cos_rand_theta_CM +&
            & a_CM_x * sin_rand_theta_CM * cos_rand_phi_CM +&
            & b_CM_x * sin_rand_theta_CM * sin_rand_phi_CM
    p1_CM_y = k1_CM_y * cos_rand_theta_CM +&
            & a_CM_y * sin_rand_theta_CM * cos_rand_phi_CM +&
            & b_CM_y * sin_rand_theta_CM * sin_rand_phi_CM
    p1_CM_z = k1_CM_z * cos_rand_theta_CM +&
            & a_CM_z * sin_rand_theta_CM * cos_rand_phi_CM +&
            & b_CM_z * sin_rand_theta_CM * sin_rand_phi_CM

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Lorentz boost of particle 4-velocity from CoM to lab frame
    beta_CM_x = -beta_CM_x
    beta_CM_y = -beta_CM_y
    beta_CM_z = -beta_CM_z

    ! electron/positron 4-velocity
    prtl1_CM_u = gamma_prtl_CM * beta_prtl_CM * p1_CM_x
    prtl1_CM_v = gamma_prtl_CM * beta_prtl_CM * p1_CM_y
    prtl1_CM_w = gamma_prtl_CM * beta_prtl_CM * p1_CM_z

    call LorentzBoost(beta_CM_x, beta_CM_y, beta_CM_z,&
                    & beta_CM_sq, gamma_CM,&
                    & prtl1_CM_u, prtl1_CM_v, prtl1_CM_w,&
                    & prtl1_u, prtl1_v, prtl1_w)
    ! same CoM 4-velocity with a "-" sign
    call LorentzBoost(beta_CM_x, beta_CM_y, beta_CM_z,&
                    & beta_CM_sq, gamma_CM,&
                    & -prtl1_CM_u, -prtl1_CM_v, -prtl1_CM_w,&
                    & prtl2_u, prtl2_v, prtl2_w)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Create electron and positron with the generated momenta
    ! choose random coordinate on a tile
    call generateCoordInRegion(REAL(tile_x1), REAL(tile_x2),&
                             & REAL(tile_y1), REAL(tile_y2),&
                             & REAL(tile_z1), REAL(tile_z2),&
                             & x_new, y_new, z_new,&
                             & xi_new, yi_new, zi_new,&
                             & dx_new, dy_new, dz_new)

    ! create an electron/positron pair in the same location
    call createParticle(BW_electron_sp, xi_new, yi_new, zi_new,&
                                      & dx_new, dy_new, dz_new,&
                                      & REAL(prtl1_u, 4), REAL(prtl1_v, 4), REAL(prtl1_w, 4))
    call createParticle(BW_positron_sp, xi_new, yi_new, zi_new,&
                                      & dx_new, dy_new, dz_new,&
                                      & REAL(prtl2_u, 4), REAL(prtl2_v, 4), REAL(prtl2_w, 4))
  end subroutine PPfromTwoPhotons

  subroutine generateRandomThetaBW(SS, theta_final)
    implicit none
    real(kind=8), intent(in)  :: SS
    real(kind=8), intent(out) :: theta_final
    real(kind=8)              :: rand_prob, rand_theta
    integer                   :: iter
    iter = 0
    do while (.true.)
      rand_prob = random(dseed)
      rand_theta = M_PI * random(dseed)
      if (rand_prob .le. dSdO_BW(SS, rand_theta)) then
        exit
      end if
      iter = iter + 1
      #ifdef DEBUG
        if (iter .gt. 10000) then
          call throwError('Too many iterations in `generateRandomThetaBW()`.')
        end if
      #endif
    end do
    theta_final = rand_theta
  end subroutine generateRandomThetaBW

  real function dSdO_BW(s, theta)
    implicit none
    real(kind=8), intent(in)  :: s, theta
    real(kind=8)              :: beta, beta2, beta4
    beta2 = 1.0 - 1.0 / s
    beta4 = beta2**2
    beta = sqrt(beta2)
    dSdO_BW = (beta * (-8.0 - 8.0 * beta2 + 11.0 * beta4 -&
            & 4.0 * beta2 * (-2.0 + beta2) * cos(2.0 * theta) +&
            & beta4 * cos(4.0 * theta)) * sin(theta)) /&
            & (4.0 * (-(beta * (-2.0 + beta2)) + (-3.0 + beta4) * asinh(beta / sqrt(1.0 - beta2))) *&
            & (-2.0 + beta2 + beta2 * cos(2.0 * theta))**2)
  end function dSdO_BW

#endif
end module m_bwpairproduction
