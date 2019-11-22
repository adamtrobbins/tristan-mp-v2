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
      call prtlToSet(ti, tj, tk, sp_arr_1, n_sp_1, set_1, set_size_1)
      call shuffleSet(set_1, set_size_1)

      call prtlToSet(ti, tj, tk, sp_arr_2, n_sp_2, set_2, set_size_2)
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

    if (n_sp_2 .ne. 0) then
      call coupleParticlesOnTile(ti, tj, tk, sp_arr_1, n_sp_1, sp_arr_2, n_sp_2,&
                               & pairs_of_photons, num_pairs)
    else
      call coupleParticlesOnTile(ti, tj, tk, sp_arr_1, n_sp_1, sp_arr_1, n_sp_1,&
                               & pairs_of_photons, num_pairs)
    end if

    do ph = 1, num_pairs
      ! compute P_12 for each pair of photons `pairs_of_photons(ph)`
      call computeBWCrossSection(ti, tj, tk, pairs_of_photons(ph),&
                               & P_12, thresholdQ)
      ! to match the optical depth with the binary pairing case:
      P_12 = P_12 * num_pairs * 2
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

  subroutine PPfromTwoPhotons(ti, tj, tk, pair_of_photons)
    implicit none
    integer, intent(in)      :: ti, tj, tk
    type(couple), intent(in) :: pair_of_photons
    integer                  :: s1, s2, p1, p2
    real                     :: x_new, y_new, z_new, dx_new, dy_new, dz_new, rnd
    integer                  :: tile_x1, tile_x2, tile_y1, tile_y2, tile_z1, tile_z2
    integer(kind=2)          :: xi_new, yi_new, zi_new

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

    ! choose random coordinate on a tile
    call generateCoordInRegion(REAL(tile_x1), REAL(tile_x2),&
                             & REAL(tile_y1), REAL(tile_y2),&
                             & REAL(tile_z1), REAL(tile_z2),&
                             & x_new, y_new, z_new,&
                             & xi_new, yi_new, zi_new,&
                             & dx_new, dy_new, dz_new)

    ! create an electron/positron pair in the same location
    call createParticle(BW_electron_sp, xi_new, yi_new, zi_new,&
                                      & dx_new, dy_new, dz_new, 0.0, 0.0, 0.0)
    call createParticle(BW_positron_sp, xi_new, yi_new, zi_new,&
                                      & dx_new, dy_new, dz_new, 0.0, 0.0, 0.0)
  end subroutine PPfromTwoPhotons

#endif
end module m_bwpairproduction
