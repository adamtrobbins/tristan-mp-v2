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
  logical :: Compton_el_recoil
  real    :: Thomson_lim

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
      ! group1 (electrons/positrons) or to group2 (photons):
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
    real(kind=8)              :: el_gamma, el_beta, pel_x, pel_y, pel_z
    real(kind=8)              :: eph, kph_x, kph_y, kph_z
    real(kind=8)              :: eph_RF, kph_RF_x, kph_RF_y, kph_RF_z
    real(kind=8)              :: eph_new, kph_new_x, kph_new_y, kph_new_z
    real, pointer             :: u_el, v_el, w_el, u_ph, v_ph, w_ph

    ! couple the electrons/positrons (group1) and photons (group2): 
    call coupleParticlesOnTile(ti, tj, tk, sp_arr_1, n_sp_1, sp_arr_2, n_sp_2,&
                             & el_photon_pairs, num_pairs, num_1, num_2)

    do el_ph = 1, num_pairs
      ! "extract" the el-photon pair:
      s1 = el_photon_pairs(el_ph)%part_1%spec
      p1 = el_photon_pairs(el_ph)%part_1%index
      s2 = el_photon_pairs(el_ph)%part_2%spec
      p2 = el_photon_pairs(el_ph)%part_2%index

      #ifdef DEBUG
        if ((species(s1)%m_sp .ne. 1) .or. (abs(species(s1)%ch_sp) .ne. 1)) then
            call throwError('Wrong particle assigned to electron/positron group in Compton particle coupling!')
        endif
        if ((species(s2)%m_sp .ne. 0) .or. (species(s2)%ch_sp .ne. 0)) then
          call throwError('Wrong particle assigned to photon group in Compton particle coupling!')
        end if
      #endif

      u_el => species(s1)%prtl_tile(ti, tj, tk)%u(p1)
      v_el => species(s1)%prtl_tile(ti, tj, tk)%v(p1)
      w_el => species(s1)%prtl_tile(ti, tj, tk)%w(p1)
      u_ph => species(s2)%prtl_tile(ti, tj, tk)%u(p2)
      v_ph => species(s2)%prtl_tile(ti, tj, tk)%v(p2)
      w_ph => species(s2)%prtl_tile(ti, tj, tk)%w(p2)

      pel_x = REAL(u_el, 8); pel_y = REAL(v_el, 8); pel_z = REAL(w_el, 8)
      kph_x = REAL(u_ph, 8); kph_y = REAL(v_ph, 8); kph_z = REAL(w_ph, 8)

      el_gamma = 1.0 + pel_x**2 + pel_y**2 + pel_z**2 ! here gamma^2
      el_beta = sqrt(1.0 - 1.0 / el_gamma)
      el_gamma = sqrt(el_gamma)
      eph = sqrt(kph_x**2 + kph_y**2 + kph_z**2)

      ! boost photon momentum into electron frame:
      call boostPhoton(el_gamma, el_beta, pel_x, pel_y, pel_z, &
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
        pel_x = -pel_x; pel_y = -pel_y; pel_z = -pel_z
        call boostPhoton(el_gamma, el_beta, pel_x, pel_y, pel_z, &
                            & eph_RF, kph_RF_x, kph_RF_y, kph_RF_z, &
                            & eph_new, kph_new_x, kph_new_y, kph_new_z)
        ! obtain the recoil on the electron via momentum conservation:
        if (Compton_el_recoil) then
          u_el = u_el - (kph_new_x - u_ph)
          v_el = v_el - (kph_new_y - v_ph)
          w_el = w_el - (kph_new_z - w_ph)
        endif
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

    if ( eph_RF .lt. Thomson_lim ) then
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
    #ifdef DEBUG
      if ((P_12 .lt. 0.0) .or. (P_12 .gt. 1.0)) then
        print *, 'P_12 = ', P_12
        call throwError('Compton cross section P_12 out of bounds!')
      end if
    #endif
  end subroutine computeComptonCrossSection

  subroutine boostPhoton(gam, beta, p_x, p_y, p_z, &
                       & eph, k_x, k_y, k_z, &
                       & eph1, k1_x, k1_y, k1_z)
    implicit none
    real(kind=8), intent(in)  :: gam, beta, p_x, p_y, p_z 
    ! the input momentum:
    real(kind=8), intent(in)  :: eph, k_x, k_y, k_z
    ! the transformed momentum:
    real(kind=8), intent(out) :: eph1, k1_x, k1_y, k1_z
    real(kind=8)              :: p_dot_k
    real                      :: check

    p_dot_k = p_x * k_x + p_y * k_y + p_z * k_z
    ! The transformed photon momentum:
    eph1 = gam * eph - p_dot_k
    k1_x = k_x + (p_dot_k / (1.0 + gam) - eph) * p_x 
    k1_y = k_y + (p_dot_k / (1.0 + gam) - eph) * p_y 
    k1_z = k_z + (p_dot_k / (1.0 + gam) - eph) * p_z 
    #ifdef DEBUG
      check = eph1 / eph
      if ((check .lt. 0.0) .or. (check .gt. 2.0*gam)) then
        print *, 'Photon ratio energy1 / energy  = ', check
        call throwError('Error in boostPhoton().')
      end if
    #endif
  end subroutine boostPhoton

  subroutine scatterPhoton(if_KleinNishina, eph_RF, kph_RF_x, kph_RF_y, kph_RF_z)
    implicit none
    logical, intent(in)          :: if_KleinNishina
    real(kind=8), intent(inout)  :: eph_RF, kph_RF_x, kph_RF_y, kph_RF_z
    real(kind=8)                 :: a_RF_x, a_RF_y, a_RF_z
    real(kind=8)                 :: b_RF_x, b_RF_y, b_RF_z
    real(kind=8)                 :: c_RF_x, c_RF_y, c_RF_z
    real(kind=8)                 :: rand_costheta_RF, rand_sintheta_RF
    real(kind=8)                 :: rand_phi_RF, rand_cosphi_RF, rand_sinphi_RF, norm

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Define a basis in the electron frame:
    norm = 1.0 / eph_RF
    a_RF_x = kph_RF_x * norm; a_RF_y = kph_RF_y * norm; a_RF_z = kph_RF_z * norm
    if (a_RF_x .ne. 0.0) then
      b_RF_x = -a_RF_y / a_RF_x; b_RF_y = 1.0; b_RF_z = 0.0
      norm = 1.0 / sqrt(b_RF_x**2 + b_RF_y**2)
      b_RF_x = b_RF_x * norm; b_RF_y = b_RF_y * norm
    else
      b_RF_x = 1.0; b_RF_y = 0.0; b_RF_z = 0.0
    end if
    c_RF_x = b_RF_z * a_RF_y - b_RF_y * a_RF_z
    c_RF_y = b_RF_x * a_RF_z - b_RF_z * a_RF_x
    c_RF_z = b_RF_y * a_RF_x - b_RF_x * a_RF_y
    norm = 1.0 / sqrt(c_RF_x**2 + c_RF_y**2 + c_RF_z**2)
    c_RF_x = c_RF_x * norm; c_RF_y = c_RF_y * norm; c_RF_z = c_RF_z * norm
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Generate random vector in the electron frame ...
    ! ... respecting the differential (Klein-Nishina) cross section ...
    call generateRandomCosTheta(eph_RF, rand_costheta_RF, if_KleinNishina)
    ! just in case the returned costheta went out of bounds due to roundoffs:
    if (abs(rand_costheta_RF) .gt. 1.0d0) then
      rand_costheta_RF = sign(1.0d0, rand_costheta_RF)
    endif
    rand_phi_RF = 2.0 * M_PI * random(dseed)
    rand_sintheta_RF = sqrt(1.0d0 - rand_costheta_RF**2)
    rand_cosphi_RF = cos(rand_phi_RF)
    rand_sinphi_RF = sin(rand_phi_RF)

    ! update the momentum:
    eph_RF = eph_RF / (1.0d0 + eph_RF * (1.0d0 - rand_costheta_RF))
    kph_RF_x = eph_RF * (rand_costheta_RF * a_RF_x + &
                       & rand_sintheta_RF * rand_cosphi_RF * b_RF_x + &
                       & rand_sintheta_RF * rand_sinphi_RF * c_RF_x)
    kph_RF_y = eph_RF * (rand_costheta_RF * a_RF_y + &
                       & rand_sintheta_RF * rand_cosphi_RF * b_RF_y + &
                       & rand_sintheta_RF * rand_sinphi_RF * c_RF_y)
    kph_RF_z = eph_RF * (rand_costheta_RF * a_RF_z + &
                       & rand_sintheta_RF * rand_cosphi_RF * b_RF_z + &
                       & rand_sintheta_RF * rand_sinphi_RF * c_RF_z)
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
    real(kind=8)              :: c0, c1, c2, thresh
    integer                   :: iter
    iter = 0
    umin = -1.0
    umax = 1.0
    thresh = 1e-10
   
    rnd = random(dseed)
    if ( .not. if_KleinNishina ) then
      ! `u` for Thomson can be sampled by transforming `rnd` with 
      ! an analytic formula, which is obtained by inverting the 
      ! cumulative distribution: 
      u1 = (4.0 * rnd - 2.0 + sqrt(5.0 + 16.0 * rnd * (rnd - 1.0)))**(1.0/3.0)
      u = u1 - 1.0 / u1
    else 
      ! generate random costheta for Klein-Nishina
      ! by solving iteratively (via bisection) for F(u=cos(theta)) = rnd \in [0,1]
      c0 = 1.0 / (4.0 + 2.0 * (1.0 + eph_RF) / (1.0 / eph_RF + 2.0)**2 + &
         & (eph_RF - 2.0 / eph_RF - 2.0) * log(1.0 + 2.0 * eph_RF))
      c1 = 1.0 - 1.0 / eph_RF - 0.5 * eph_RF / (1.0 + 2.0 * eph_RF)**2
      c2 = eph_RF - 2.0 / eph_RF - 2.0
      do while (( umax - umin) .gt. thresh)
        u = 0.5 * ( umin + umax )
        Fu = Fu_KN(eph_RF, u, c0, c1, c2) - rnd
        if (abs(Fu) .lt. thresh) exit
        if ( Fu .lt. 0.0 ) then
          umin = u
        else
          umax = u
        endif 
        iter = iter + 1
        #ifdef DEBUG
          if (iter .gt. 200) then
            call throwError('Too many iterations in `generateRandomCosTheta()` in `compton` module.')
          end if
        #endif
      end do
    endif
    costheta_RF = u
    #ifdef DEBUG
      if ((REAL(costheta_RF) .lt. -1.0) .or. (REAL(costheta_RF) .gt. 1.0)) then
        print *, 'Cos(theta) = ', costheta_RF
        call throwError('Cos(theta) for Compton scattering went out of bounds.')
      end if
    #endif
  end subroutine generateRandomCosTheta

  ! cumulative distribution F \in [0,1] of u = cos(theta) for Klein-Nishina
  ! with the photon energy in electron rest frame, eph_RF, as parameter.
  real(8) function Fu_KN(eph_RF, u, c0, c1, c2)
    implicit none
    real(kind=8), intent(in)  :: eph_RF, u, c0, c1, c2
    real(kind=8)              :: g

    g = 1.0 / (1.0 + eph_RF - eph_RF*u)
    
    Fu_KN = c0 * (c1 + u + 0.5 * eph_RF * g**2 + (1.0 / eph_RF + 2.0) * g + &
          & c2 * log((1.0 + 2.0 * eph_RF) * g))
  end function Fu_KN

#endif
end module m_compton
