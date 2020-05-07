#include "../defs.F90"

module m_compton
#ifdef COMPTONSCATTERING

  use m_globalnamespace
  use m_aux
  use m_errors
  use m_bincoupling
  use m_particlelogistics
  implicit none

  ! Compton scattering parameters
  real    :: Compton_tau
  integer :: Compton_interval
  integer :: Compton_algorithm

  !--- PRIVATE variables/functions -------------------------------!
  private :: comptonOnTile_bin, comptonOnTile_mc, scatterPhoton
  private :: generateRandomCosTheta, Fu_KN
  !...............................................................!
contains
  subroutine comptonScattering()
    implicit none
    integer   :: compton_species_1(20), compton_species_2(20)
    integer   :: s, si_1, si_2
    integer   :: ti, tj, tk

    ! find all the species that participate in Compton scattering
    si_1 = 0; si_2 = 0
    do s = 1, nspec
      ! go through the list and assign species either to 
      ! group1 (electrons/positrons) or to group2 (photons).
      if (species(s)%compton_sp) then
        ! electrons or positrons:
        if ( (species(s)%m_sp .eq. 1.0) .and. (abs(species(s)%ch_sp) .eq. 1.0) ) then
          compton_species_1(si_1 + 1) = s
          si_1 = si_1 + 1
        ! the photons:
        else if ( (species(s)%m_sp .eq. 0) .and. (species(s)%ch_sp .eq. 0) ) then
          compton_species_2(si_2 + 1) = s
          si_2 = si_2 + 1
        end if
      endif
    end do

    ! proceed if scattering is turned on at least for one electron/positron
    ! species and one photon species:
    if ((si_1 .ne. 0) .and. (si_2 .ne. 0)) then

      ! assuming all tiles are equal
      s = compton_species_1(1)
      ! loop through all the tiles
      do ti = 1, species(s)%tile_nx
        do tj = 1, species(s)%tile_ny
          do tk = 1, species(s)%tile_nz
            if (Compton_algorithm .eq. 1) then
              ! TODO: not yet implemented...
              !call comptonOnTile_bin(ti, tj, tk,&
              !                & compton_species_1(1 : si_1), si_1,&
              !                & compton_species_2(1 : si_2), si_2)
            else if (Compton_algorithm .eq. 2) then
              call comptonOnTile_mc(ti, tj, tk,&
                             & compton_species_1(1 : si_1), si_1,&
                             & compton_species_2(1 : si_2), si_2)
            end if
          end do
        end do
      end do
    endif
  end subroutine comptonScattering

  subroutine comptonOnTile_bin(ti, tj, tk,&
                        & sp_arr_1, n_sp_1,&
                        & sp_arr_2, n_sp_2)
    implicit none
    integer, intent(in)       :: ti, tj, tk
    integer, intent(in)       :: n_sp_1 ! # of species in set
    integer, intent(in)       :: sp_arr_1(n_sp_1)  ! el/positrons
    integer, intent(in)       :: n_sp_2 ! # of species in set
    integer, intent(in)       :: sp_arr_2(n_sp_2) ! photons

    !TODO
  end subroutine comptonOnTile_bin

  subroutine comptonOnTile_mc(ti, tj, tk,&
                       & sp_arr_1, n_sp_1,&
                       & sp_arr_2, n_sp_2)
    implicit none
    integer, intent(in)       :: ti, tj, tk
    integer, intent(in)       :: n_sp_1 ! # of species in set
    integer, intent(in)       :: sp_arr_1(n_sp_1)  ! el/positrons
    integer, intent(in)       :: n_sp_2 ! # of species in set
    integer, intent(in)       :: sp_arr_2(n_sp_2) ! photons
    type(couple), allocatable :: el_photon_pairs(:)
    integer                   :: num_pairs, el_ph, s1, s2, p1, p2
    real                      :: rnd, P_12
    logical                   :: if_KleinNishina
    integer                   :: num_1, num_2
    real(kind=8)              :: el_gamma, el_beta, kel_x, kel_y, kel_z
    real(kind=8)              :: eph, kph_x, kph_y, kph_z, norm
    real(kind=8)              :: eph_RF, kph_RF_x, kph_RF_y, kph_RF_z
    real(kind=8)              :: eph_new, kph_new_x, kph_new_y, kph_new_z
    real, pointer             :: u_el, v_el, w_el, u_ph, v_ph, w_ph

    ! shuffle sets for more randomness
    if (random(dseed) .gt. 0.5) then
      call coupleParticlesOnTile(ti, tj, tk, sp_arr_1, n_sp_1, sp_arr_2, n_sp_2,&
                               & el_photon_pairs, num_pairs, num_1, num_2)
    else
      call coupleParticlesOnTile(ti, tj, tk, sp_arr_2, n_sp_2, sp_arr_1, n_sp_1,&
                               & el_photon_pairs, num_pairs, num_1, num_2)
    end if

    do el_ph = 1, num_pairs
      ! "extract" the el-photon pair:
      s1 = el_photon_pairs(el_ph)%part_1%spec
      p1 = el_photon_pairs(el_ph)%part_1%index
      s2 = el_photon_pairs(el_ph)%part_2%spec
      p2 = el_photon_pairs(el_ph)%part_2%index

      u_el => species(s1)%prtl_tile(ti, tj, tk)%u(p1)
      v_el => species(s1)%prtl_tile(ti, tj, tk)%v(p1)
      w_el => species(s1)%prtl_tile(ti, tj, tk)%w(p1)
      u_ph => species(s2)%prtl_tile(ti, tj, tk)%u(p2)
      v_ph => species(s2)%prtl_tile(ti, tj, tk)%v(p2)
      w_ph => species(s2)%prtl_tile(ti, tj, tk)%w(p2)

      kel_x = REAL(u_el, 8); kel_y = REAL(v_el, 8); kel_z = REAL(w_el, 8)
      kph_x = REAL(u_ph, 8); kph_y = REAL(v_ph, 8); kph_z = REAL(w_ph, 8)

      el_gamma = 1.0 + kel_x**2 + kel_y**2 + kel_z**2 ! here gamma^2
      el_beta = sqrt(1.0 - 1.0 / el_gamma)
      norm = 1.0 / sqrt(el_gamma - 1.0)
      el_gamma = sqrt(el_gamma)
      kel_x = kel_x * norm; kel_y = kel_y * norm; kel_z = kel_z * norm

      eph = sqrt(kph_x**2 + kph_y**2 + kph_z**2)
      norm = 1.0 / eph
      kph_x = kph_x * norm; kph_y = kph_y * norm; kph_z = kph_z * norm

      ! boost photon momentum into electron frame:
      call boostPhoton(el_gamma, el_beta, kel_x, kel_y, kel_z, &
                             & eph, kph_x, kph_y, kph_z, &
                             & eph_RF, kph_RF_x, kph_RF_y, kph_RF_z)

      ! compute cross section:
      call computeComptonCrossSection(eph, eph_RF, el_gamma, P_12, if_KleinNishina)

      ! to match the optical depth with the binary pairing case:
      P_12 = P_12 * num_2
      rnd = random(dseed)
      if (rnd .le. P_12) then
        ! TODO: Split particles if el and photon weight are not equal!

        ! scatter the photon in the electron rest frame:
        call scatterPhoton(if_KleinNishina, eph_RF, kph_RF_x, kph_RF_y, kph_RF_z)

        ! boost back into lab frame:
        kel_x = -kel_x; kel_y = -kel_y; kel_z = -kel_z
        call boostPhoton(el_gamma, el_beta, kel_x, kel_y, kel_z, &
                            & eph_RF, kph_RF_x, kph_RF_y, kph_RF_z, &
                            & eph_new, kph_new_x, kph_new_y, kph_new_z)
        ! obtain the recoil on the electron via momentum conservation:
        kph_new_x = kph_new_x * eph_new
        kph_new_y = kph_new_y * eph_new
        kph_new_z = kph_new_z * eph_new
        u_el = u_el - (kph_new_x - u_ph)
        v_el = v_el - (kph_new_y - v_ph)
        w_el = w_el - (kph_new_z - w_ph)
        ! store the new photon momentum:
        u_ph = kph_new_x
        v_ph = kph_new_y
        w_ph = kph_new_z
      end if
      u_el => null(); v_el => null(); w_el => null()
      u_ph => null(); v_ph => null(); w_ph => null()
    end do
  end subroutine comptonOnTile_mc

  subroutine computeComptonCrossSection(eph, eph_RF, el_gamma, P_12, if_KleinNishina)
    implicit none
    real(kind=8), intent(in)  :: eph, eph_RF, el_gamma
    real, intent(out)         :: P_12
    logical, intent(out)      :: if_KleinNishina
    real(kind=8)              :: over_eph, f_KN

    if ( eph_RF .lt. 1e-2 ) then
      if_KleinNishina = .false.  ! use classical Thomson cross-section
      P_12 = REAL(Compton_interval) * Compton_tau
    else
      if_KleinNishina = .true.  ! Klein-Nishina
      over_eph = 1.0 / eph_RF
      f_KN = ((1.0 - 2.0 * over_eph - 2.0 * over_eph**2) * log(1.0 + 2.0 * eph_RF) + &
           & 0.5 + 4.0 * over_eph - 0.5 / (1.0 + 2.0 * eph_RF)**2) * over_eph * 3.0 / 8.0
      P_12 = REAL(Compton_interval) * Compton_tau * f_KN
    end if
    ! Now transform cross section back into lab frame:
    P_12 = P_12 * eph_RF / (el_gamma * eph)
  end subroutine computeComptonCrossSection

  subroutine boostPhoton(gam, beta, n_x, n_y, n_z, &
                       & eph, k_x, k_y, k_z, &
                       & eph1, k1_x, k1_y, k1_z)
    implicit none
    real(kind=8), intent(in) :: gam, beta, n_x, n_y, n_z 
    ! the input momentum (energy & direction of motion):
    real(kind=8), intent(in) :: eph, k_x, k_y, k_z
    ! the transformed momentum:
    real(kind=8), intent(out) :: eph1, k1_x, k1_y, k1_z
    real(kind=8)              :: norm, costheta

    costheta = n_x * k_x + k_y * n_y + k_z * n_z
    ! The transformed photon energy and direction of motion:
    eph1 = eph * gam * (1.0 - beta * costheta)
    norm = eph / eph1
    k1_x = norm * (k_x + ((gam - 1.0) * costheta - gam * beta) * n_x)
    k1_y = norm * (k_y + ((gam - 1.0) * costheta - gam * beta) * n_y)
    k1_z = norm * (k_z + ((gam - 1.0) * costheta - gam * beta) * n_z)
    #ifdef DEBUG
      ! consistency check:
      if (sqrt(k1_x**2 + k1_y**2 + k1_z**2) .ne. 1.0) then
        print *, 'Error in photon boost! Transformed K vector does not have correct norm.'
        stop
      end if
    #endif
  end subroutine boostPhoton

  subroutine scatterPhoton(if_KleinNishina, eph_RF, kph_RF_x, kph_RF_y, kph_RF_z)
    implicit none
    logical, intent(in)          :: if_KleinNishina
    real(kind=8), intent(inout)  :: eph_RF, kph_RF_x, kph_RF_y, kph_RF_z
    real(kind=8)                 :: a_RF_x, a_RF_y, a_RF_z, norm
    real(kind=8)                 :: b_RF_x, b_RF_y, b_RF_z
    real(kind=8)                 :: rand_costheta_RF, rand_sintheta_RF
    real(kind=8)                 :: rand_phi_RF, rand_cosphi_RF, rand_sinphi_RF

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define a basis in the electron frame:
    if (kph_RF_x .ne. 0.0) then
      a_RF_x = -kph_RF_y / kph_RF_x; a_RF_y = 1.0; a_RF_z = 0.0
      norm = 1.0 / sqrt(a_RF_x**2 + a_RF_y**2)
      a_RF_x = a_RF_x * norm; a_RF_y = a_RF_y * norm
    else
      a_RF_x = 1.0; a_RF_y = 0.0; a_RF_z = 0.0
    end if
    b_RF_x = a_RF_z * kph_RF_y - a_RF_y * kph_RF_z
    b_RF_y = a_RF_x * kph_RF_z - a_RF_z * kph_RF_x
    b_RF_z = a_RF_y * kph_RF_x - a_RF_x * kph_RF_y
    norm = 1.0 / sqrt(b_RF_x**2 + b_RF_y**2 + b_RF_z**2)
    b_RF_x = b_RF_x * norm; b_RF_y = b_RF_y * norm; b_RF_z = b_RF_z * norm
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Generate random vector in the electron frame ...
    ! ... respecting the differential (Klein-Nishina) cross section ...
    call generateRandomCosTheta(eph_RF, rand_costheta_RF, if_KleinNishina)

    rand_phi_RF = 2.0 * M_PI * random(dseed)
    rand_sintheta_RF = sqrt( 1 - rand_costheta_RF**2 )
    rand_cosphi_RF = cos(rand_phi_RF)
    rand_sinphi_RF = sin(rand_phi_RF)

    ! update the momentum:
    kph_RF_x = rand_costheta_RF * kph_RF_x + &
             & rand_sintheta_RF * rand_cosphi_RF * a_RF_x + &
             & rand_sintheta_RF * rand_sinphi_RF * b_RF_x
    kph_RF_y = rand_costheta_RF * kph_RF_y + &
             & rand_sintheta_RF * rand_cosphi_RF * a_RF_y + &
             & rand_sintheta_RF * rand_sinphi_RF * b_RF_y
    kph_RF_z = rand_costheta_RF * kph_RF_z + &
             & rand_sintheta_RF * rand_cosphi_RF * a_RF_z + &
             & rand_sintheta_RF * rand_sinphi_RF * b_RF_z

    eph_RF = eph_RF / (1.0 + eph_RF * (1.0 - rand_costheta_RF))
  end subroutine scatterPhoton

  ! sample random u=cos(theta) in the electron frame for the Thomson or
  ! Klein-Nishina differential cross-section. The cross-sections are
  ! for un-polarized photons.
  subroutine generateRandomCosTheta(eph_RF, costheta_RF, if_KleinNishina)
    implicit none
    real(kind=8), intent(in)  :: eph_RF
    real(kind=8), intent(out) :: costheta_RF
    logical, intent(in)       :: if_KleinNishina
    real(kind=8)              :: rnd, u, Fu, umin, umax, u1
    real(kind=8)              :: c0, c1, c2
    integer                   :: iter
    iter = 0
    umin = -1.0
    umax = 1.0
   
    rnd = random(dseed)
    if ( .not. if_KleinNishina ) then
      ! `u` for Thomson can be sampled by transforming `rnd` with 
      ! an analytic formula, which is obtained by inverting the 
      ! cumulative distribution: 
      u1 = (4.0 * rnd - 2.0 + sqrt(5.0 + 16.0 * rnd * (rnd - 1.0)))**(1.0/3.0)
      u = (u1**2 - 1.0) / u1
    else 
      ! generate random costheta for Klein-Nishina
      ! by solving iteratively (via bisection) for F(u=cos(theta)) = rnd \in [0,1]
      c0 = 1.0 / (4.0 + 2.0 * (1.0 + eph_RF) / (1.0 / eph_RF + 2.0)**2 + &
         & (eph_RF - 2.0 / eph_RF - 2.0) * log(1.0 + 2.0 * eph_RF))
      c1 = 1.0 - 1.0 / eph_RF - 0.5 * eph_RF / (1.0 + 2.0 * eph_RF)**2
      c2 = eph_RF - 2.0 / eph_RF - 2.0
      do while (( umax - umin) > 1e-14)
        u = 0.5 * ( umin + umax )
        Fu = Fu_KN(eph_RF, u, c0, c1, c2)
        if ( (Fu - rnd) .lt. 0 ) then
          umin = u
        else
          umax = u
        endif 
        iter = iter + 1
        #ifdef DEBUG
          if (iter .gt. 10000) then
            call throwError('Too many iterations in `generateRandomCosTheta()` in `compton` module.')
          end if
        #endif
      end do
    endif
    costheta_RF = u
  end subroutine generateRandomCosTheta

  ! cumulative distribution F \in [0,1] of u = cos(theta) for Klein-Nishina
  ! with the photon energy in electron rest frame, eph_RF, as parameter.
  real(8) function Fu_KN(eph_RF, u, c0, c1, c2)
    implicit none
    real(kind=8), intent(in)  :: eph_RF, u, c0, c1, c2 ! u = cos(theta)
    real(kind=8)              :: g

    g = 1.0 / (1.0 + eph_RF - eph_RF*u)
    
    Fu_KN = c0 * (c1 + u + 0.5 * eph_RF * g**2 + (1.0 / eph_RF + 2.0) * g + &
          & c2 * log((1.0 + 2.0 * eph_RF) * g))
  end function Fu_KN

#endif
end module m_compton
