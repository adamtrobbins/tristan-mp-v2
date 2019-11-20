#ifdef PHOTONS
#if defined(PHOTONMERGE) || defined(PAIRS)
#ifdef twoD

module m_photon_interact
    use m_globaldata
	use m_system
	use m_aux
	use m_communications
	use m_fields
	use m_inputparser
	use m_fparser
	use m_particles
    use m_user

#else

module m_photon_interact_3d
    use m_globaldata_3d
	use m_system_3d
	use m_aux_3d
	use m_communications_3d
	use m_fields_3d
	use m_inputparser_3d
	use m_fparser_3d
    use m_particles_3d
    use m_user_3d

#endif

!-------------------------------------------------------------------------------
! 						subroutine interact_photons()
!
!   Look for Breit-Wheeler pair formation events and create pairs
!
!-------------------------------------------------------------------------------
    implicit none

    public :: interact_photons

    contains

    subroutine interact_photons()
        implicit none
        integer, allocatable :: c_num_photons(:)
        integer :: ind_ph
        logical, parameter :: tmstop = .true.

        if (modulo(lap, ph_step) .eq. 0) then
            allocate(c_num_photons(mx * my * mz))
            c_num_photons(:) = 0
            call sort_photons(c_num_photons)
            call timer(36)
#ifdef PHOTONMERGE
            if (merge_method .eq. 2) then
                call merge_photons(c_num_photons)
            else if (merge_method .eq. 1) then
                call merge_photons_conserv(c_num_photons)
            end if
#endif
            call timer(36, tmstop = tmstop)
#ifdef PAIRS
            call check_interactions(c_num_photons)
            if(debug) print *, rank, ": photon_interact done", "step=", lap, photons
#endif
            call clean_photons_up()
            if(allocated(c_num_photons)) deallocate(c_num_photons)
        else
            call timer(36)
            call timer(36, tmstop = tmstop)
        end if

    end subroutine interact_photons

    subroutine sort_photons(c_num_photons)
        implicit none
        integer, allocatable, intent(inout) :: c_num_photons(:)
        integer :: ind_ph, cell_num
        integer :: ind_temp, ind_temp2, temp_var1, temp_var2

        !-------------------------------------------------------------------------------
        !   Sorting photons here
        !-------------------------------------------------------------------------------
        ! counting the total number of photons in each cell
        do ind_ph = 1, photons
            ! only works in 2D
            cell_num = int(ph(ind_ph)%x) + int(ph(ind_ph)%y) * iy + 1
            c_num_photons(cell_num) = c_num_photons(cell_num) + 1
        end do

        ! counting the cumulative number of photons for each cell (without that cell itself)
        temp_var1 = c_num_photons(1)
        c_num_photons(1) = 0
        do ind_temp = 2, mx * my * mz
            temp_var2 = c_num_photons(ind_temp)
            c_num_photons(ind_temp) = temp_var1 + c_num_photons(ind_temp - 1)
            temp_var1 = temp_var2
        end do

        do ind_ph = 1, photons
            ! only works in 2D
            cell_num = int(ph(ind_ph)%x) + int(ph(ind_ph)%y) * iy + 1
            c_num_photons(cell_num) = c_num_photons(cell_num) + 1
            if (c_num_photons(cell_num) .ge. maxhlf - 1) then
                print *, 'ERROR in sort_photons'
                print *, rank, c_num_photons(cell_num), cell_num, int(ph(ind_ph)%x), int(ph(ind_ph)%y), iy, maxhlf
                stop
            end if
            call copyprt(ph(ind_ph), tempp(c_num_photons(cell_num)))
            ! tempp(c_num_photons(cell_num)) = ph(ind_ph)
        end do

        do ind_ph = 1, photons
            call copyprt(tempp(ind_ph), ph(ind_ph))
            ! ph(ind_ph) = tempp(ind_ph)
        end do

#ifdef PHOTDEBUG
        ! FOR DEBUG PURPOSES
        if (c_num_photons(mx * my * mz) .ne. photons) then
            print *, c_num_photons(mx * my * mz), photons
            print *, "The cumulative number of photons was computed incorrectly! Problem in photon_interact()."
            stop
        end if
#endif
        ! c_num_photons now contains the cumulative number of photons INCLUDING the one in that particular cell
        !-------------------------------------------------------------------------------
        !   / Sorting photons here
        !-------------------------------------------------------------------------------
    end subroutine sort_photons

#ifdef PHOTONMERGE

    subroutine merge_photons(c_num_photons)
        implicit none
        integer :: photon_bin_limit
        ! real(dprec), parameter :: M_PI = 3.1415926535897931d0
        real, parameter :: M_PI = 3.1415927
        real, allocatable :: bins_phi(:), bins_costheta(:), bins_energy_sq(:)
        integer, allocatable, intent(in) :: c_num_photons(:)
        real, allocatable :: lognorm(:)

        integer :: index_phi_1, index_costheta_1, index_energy_1
        integer :: index_phi_2, index_costheta_2, index_energy_2

        integer :: ind_1, ind_2, sum_ch
        integer :: i, cell_num, temp_var1, temp_var2
        real :: j
        real :: phi_kick, new_energy
        real :: E1, E2
        real :: ebin_max, log10_ebin_max, log10_ebin_min

        photon_bin_limit = int(ph_bins_num**(3. / 2.))

        allocate(bins_phi(ph_bins_num))
        allocate(bins_costheta(ph_bins_num))
        allocate(bins_energy_sq(ph_bins_num))

        ! initializing randomized momentum bins
        !      intervals have lognormal distribution
        call log_normal(ph_bins_num, lognorm)
        ! randomize maximum energy bin from epsph_max to 10 * epsph_max
        log10_ebin_max = log10(epsph_max) + random(dseed)
        ! ebin_max = 10.0**log10_ebin_max
        log10_ebin_min = log10(epsph_min)
        do i = 1, ph_bins_num
            bins_phi(i) = 2. * M_PI * lognorm(i) / lognorm(ph_bins_num)
            bins_costheta(i) = -1. + 2. * lognorm(i) / lognorm(ph_bins_num)
            ! energy squared
            bins_energy_sq(i) = 10.0**(2. * log10_ebin_min + &
                            &(1. * lognorm(i)) * 2. * &
                            &(log10_ebin_max - log10_ebin_min))
        end do
        bins_phi(ph_bins_num) = 2. * M_PI + 1e-5
        bins_costheta(ph_bins_num) = 1. + 1e-5
        ebin_max = sqrt(bins_energy_sq(ph_bins_num))
        bins_energy_sq(ph_bins_num) = bins_energy_sq(ph_bins_num) + 1.
        ! / initializing randomized momentum bins

        do cell_num = 1, mx * my * mz
            if (cell_num == 1) then
                if (c_num_photons(1) .lt. photon_bin_limit) then
                    cycle
                end if
                temp_var1 = 1
                temp_var2 = c_num_photons(1)
            else
                if (c_num_photons(cell_num) - c_num_photons(cell_num - 1) .lt. photon_bin_limit) then
                    cycle
                end if
                temp_var1 = c_num_photons(cell_num - 1) + 1
                temp_var2 = c_num_photons(cell_num)
            end if
            ! at this point temp_var1 contains the index of first photon in cell, and temp_var2 - the index of the last one

            ! random rotation in x-y plane (defined at every cell)
            phi_kick = random(dseed) * 2. * M_PI
            do ind_1 = temp_var1, temp_var2 - 1
                if ((ph(ind_1)%splitlev .eq. 1) .or. (ph(ind_1)%ch .le. 0)) cycle ! photon is scheduled for deletion
                E1 = sqrt((ph(ind_1)%u)**2 + (ph(ind_1)%v)**2 + (ph(ind_1)%w)**2)
                if (E1 .ge. ebin_max) cycle ! photon energy too high

                call find_bins(ind_1, ph_bins_num, bins_phi, bins_costheta, bins_energy_sq, &
                                &index_phi_1, index_costheta_1, index_energy_1, phi_kick)

                do ind_2 = ind_1 + 1, temp_var2
                    if ((ph(ind_2)%splitlev .eq. 1) .or. (ph(ind_2)%ch .le. 0)) cycle ! photon is scheduled for deletion
                    E2 = sqrt((ph(ind_2)%u)**2 + (ph(ind_2)%v)**2 + (ph(ind_2)%w)**2)
                    if (E2 .ge. ebin_max) cycle ! photon energy too high

                    if (ph(ind_1)%ch + ph(ind_2)%ch > max_photon_weight) cycle ! cumulative "charge" would be too high

                    call find_bins(ind_2, ph_bins_num, bins_phi, bins_costheta, bins_energy_sq, &
                                    &index_phi_2, index_costheta_2, index_energy_2, phi_kick)

                    if ((index_phi_1 .eq. index_phi_2) .and.&
                        & (index_costheta_1 .eq. index_costheta_2) .and.&
                        & (index_energy_1 .eq. index_energy_2)) then ! merge photons
                        sum_ch = ph(ind_1)%ch + ph(ind_2)%ch
                        ph(ind_1)%u = (ph(ind_1)%u * ph(ind_1)%ch + ph(ind_2)%u * ph(ind_2)%ch) / sum_ch
                        ph(ind_1)%v = (ph(ind_1)%v * ph(ind_1)%ch + ph(ind_2)%v * ph(ind_2)%ch) / sum_ch
                        ph(ind_1)%w = (ph(ind_1)%w * ph(ind_1)%ch + ph(ind_2)%w * ph(ind_2)%ch) / sum_ch

                        ph(ind_1)%ch = sum_ch
                        new_energy = sqrt(ph(ind_1)%u**2 + ph(ind_1)%v**2 + ph(ind_1)%w**2)

                        ! schedule second photon for deletion
                        ph(ind_2)%splitlev = 1
                        ph(ind_2)%ch = 0
                        if (new_energy .lt. epsph_min) then
                            ph(ind_1)%splitlev = 1
                            exit
                        end if
                    end if
                end do
            end do
        end do
        if(allocated(bins_phi)) deallocate(bins_phi)
        if(allocated(bins_costheta)) deallocate(bins_costheta)
        if(allocated(bins_energy_sq)) deallocate(bins_energy_sq)
    end subroutine merge_photons

    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
    !
    !   The algorithm is adapted from Vranic et al. 2014 [1411.2248v1]
    !       momentum bins are:
    !                uniform in direction: phi in x-y and cos(theta) in z
    !                linear in absolute value (squared) and not fixed in upper bound
    !                   this ensures the distribution function is not condensed down
    !
    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
    subroutine merge_photons_conserv(c_num_photons)
        integer, allocatable, intent(in) :: c_num_photons(:)
        real, parameter :: M_PI = 3.1415927
        integer :: cell_num, ind_ph
        integer :: i
        integer :: temp_var1, temp_var2, photon_bin_limit
        integer :: index_phi, index_theta, index_energy
        integer :: index_1, index_2, index_temp
        integer, allocatable :: mom_bins(:,:,:,:)

        real :: tot_px, tot_py, tot_pz, tot_e, tot_p

        integer :: tot_ch, wa_temp, wb_temp
        real :: Ea_temp, Eb_temp, Pax_temp, Pay_temp, Paz_temp, Pbx_temp, Pby_temp, Pbz_temp
        real :: costh_temp, sinth_temp
        real :: dx_temp, dy_temp, dz_temp
        real :: phi_mid, costheta_mid
        real :: e1x, e1y, e1z, e2x, e2y, e2z, e2_dot_d, e1sqr
        real :: Energy_ph

        real, allocatable :: bins_phi(:), bins_costheta(:), bins_energy_sq(:)
        real, allocatable :: lognorm(:)
        real :: ebin_max, log10_ebin_max, log10_ebin_min
        integer :: maxnum_temp

        photon_bin_limit = int(ph_bins_num**(3. / 2.))

        ! initializing randomized momentum bins
        !      intervals have lognormal distribution
        call log_normal(ph_bins_num, lognorm)
        log10_ebin_max = log10(epsph_max)
        log10_ebin_min = log10(epsph_min)

        allocate(bins_phi(ph_bins_num))
        allocate(bins_costheta(ph_bins_num))
        allocate(bins_energy_sq(ph_bins_num))

        do i = 1, ph_bins_num
            bins_phi(i) = 2. * M_PI * lognorm(i) / lognorm(ph_bins_num)
            bins_costheta(i) = -1. + 2. * lognorm(i) / lognorm(ph_bins_num)
            ! energy squared
            ! log-bins:
            ! bins_energy_sq(i) = 10.0**(2. * log10_ebin_min + &
            !                 &(1. * lognorm(i)) * 2. * &
            !                 &(log10_ebin_max - log10_ebin_min))
            ! lin-bins with unfixed max energy:
            bins_energy_sq(i) = (epsph_min + (epsph_max - epsph_min) * lognorm(i))**2
        end do
        bins_phi(ph_bins_num) = 2. * M_PI + 1e-5
        bins_costheta(ph_bins_num) = 1. + 1e-5
        ebin_max = sqrt(bins_energy_sq(ph_bins_num))
        bins_energy_sq(ph_bins_num) = bins_energy_sq(ph_bins_num) + 1e-5
        ! / initializing randomized momentum bins

        ! just to check/
        ! if (rank .eq. 0) then
        !     print *, 'Printing bins'
        !     print *, 'phi | costheta | E^2'
        !     do i = 1, ph_bins_num
        !         print *, bins_phi(i), '|', bins_costheta(i), '|', bins_energy_sq(i)
        !     end do
        ! end if
        ! /just to check

        do cell_num = 1, mx * my * mz
            ! making sure there's enough photons in a cell
            if (cell_num == 1) then
                if (c_num_photons(1) .lt. photon_bin_limit) then
                    cycle
                end if
                temp_var1 = 1
                temp_var2 = c_num_photons(1)
            else
                if (c_num_photons(cell_num) - c_num_photons(cell_num - 1) .lt. photon_bin_limit) then
                    cycle
                end if
                temp_var1 = c_num_photons(cell_num - 1) + 1
                temp_var2 = c_num_photons(cell_num)
            end if
            ! at this point temp_var1 contains the index of first photon in cell, and temp_var2 - the index of the last one

            if (.not. allocated(mom_bins)) allocate(mom_bins(ph_bins_num, ph_bins_num, ph_bins_num, temp_var2 - temp_var1 + 5))
            if (allocated(mom_bins)) then
              mom_bins(:,:,:,:) = 0 ! initialize with 0's
            else
              print *, 'ERROR: mom_bins not allocated'
              print *, ph_bins_num, temp_var2 - temp_var1 + 5
              stop
            end if

            ! surf through all photons of a given cell and write all indices
            !   mom_bins(i, j, k, :) contains array
            !       mom_bins(i, j, k, 1) -> number of all photons in the given i,j,k momentum bin
            !       mom_bins(i, j, k, >1) -> indices of all photons in the given i,j,k momentum bin

            do ind_ph = temp_var1, temp_var2 ! all the photons in a given cell
                if ((ph(ind_ph)%splitlev .eq. 1) .or. (ph(ind_ph)%ch .le. 0)) cycle ! photon is scheduled for deletion
                Energy_ph = sqrt((ph(ind_ph)%u)**2 + (ph(ind_ph)%v)**2 + (ph(ind_ph)%w)**2)
                if (Energy_ph .ge. ebin_max) cycle ! photon energy too high
                if (ph(ind_ph)%ch .gt. max_photon_weight) cycle ! photon weight too high
                call find_bins(ind_ph, ph_bins_num, bins_phi, bins_costheta, bins_energy_sq, &
                                &index_phi, index_theta, index_energy, 0.0)
#ifdef PHOTDEBUG
                ! for debug purposes
                if ((index_phi .gt. ph_bins_num) .or. &
                    &(index_theta .gt. ph_bins_num) .or. &
                    &(index_energy .gt. ph_bins_num)) then
                    print *, 'ERROR: index out of bounds in merge, probably find_bins() error'
                    print *, ph(ind_ph)%u, ph(ind_ph)%v, ph(ind_ph)%w, Energy_ph, ebin_max
                    print *, index_phi, index_theta, index_energy, ph_bins_num
                    do i = 1, ph_bins_num
                        print *, sqrt(bins_energy_sq(i)), bins_costheta(i), bins_phi(i)
                    end do
                    stop
                end if
#endif
                mom_bins(index_phi, index_theta, index_energy, 1) = mom_bins(index_phi, index_theta, index_energy, 1) + 1
                i = 1 + mom_bins(index_phi, index_theta, index_energy, 1)
#ifdef PHOTDEBUG
                if (i .ge. temp_var2 - temp_var1 + 5) then
                    print *, "ERROR: wrong 'i' in merge_photons, mom_bins", rank, i, temp_var2 - temp_var1
                end if
#endif
                mom_bins(index_phi, index_theta, index_energy, i) = ind_ph
            end do

            ! <editor-fold
            ! now mom_bins contains for every momentum bin indices of particles in it
            do index_phi = 1, ph_bins_num
                do index_theta = 1, ph_bins_num
                    do index_energy = 1, ph_bins_num
                        if (mom_bins(index_phi, index_theta, index_energy, 1) .lt. 4) then ! if less than 4 photons in the momentum bin
                            cycle
                        end if
#ifdef PHOTDEBUG
                        ! for debug purposes
                        if (mom_bins(index_phi, index_theta, index_energy, 1) .gt. temp_var2 - temp_var1 + 5) then
                            print *, 'ERROR: wrong mom_bins allocation in merge'
                            print *, index_phi, index_theta, index_energy, ph_bins_num, temp_var2 - temp_var1
                            stop
                        end if
#endif
                        ! computing total momentum + energy + weight in the given momentum bin
                        !   for all the photons in a given bin
                        tot_px = 0
                        tot_py = 0
                        tot_pz = 0
                        tot_e = 0
                        tot_ch = 0
                        do index_temp = 2, mom_bins(index_phi, index_theta, index_energy, 1) + 1
                            ind_ph = mom_bins(index_phi, index_theta, index_energy, index_temp)
#ifdef PHOTDEBUG
                            ! for debug purposes
                            if (ind_ph .eq. 0) then
                                print *, 'ERROR: ind_ph .eq. 0 in merge'
                                stop
                            end if
#endif
                            ! delete all the photons in that bin (later we will keep 2 of them)
                            ph(ind_ph)%splitlev = 1

                            tot_px = tot_px + ph(ind_ph)%u * ph(ind_ph)%ch
                            tot_py = tot_py + ph(ind_ph)%v * ph(ind_ph)%ch
                            tot_pz = tot_pz + ph(ind_ph)%w * ph(ind_ph)%ch
                            tot_e = tot_e + sqrt(ph(ind_ph)%u**2 + ph(ind_ph)%v**2 + ph(ind_ph)%w**2) * ph(ind_ph)%ch
                            tot_ch = tot_ch + ph(ind_ph)%ch
                        end do
#ifdef PHOTDEBUG
                        ! for debug purposes
                        if (tot_ch .eq. 0) then
                            print *, 'ERROR: tot_ch .eq. 0 in merge'
                            stop
                        end if
#endif
                        tot_p = sqrt(tot_px**2 + tot_py**2 + tot_pz**2)
                        ! determining weights of two new photons
                        if (int(tot_ch / 2) * 2 .ne. tot_ch) then
                            ! odd total weight
                            wa_temp = int(tot_ch / 2)
                            wb_temp = int(tot_ch / 2) + 1
                        else
                            ! even total weight
                            wa_temp = int(tot_ch / 2)
                            wb_temp = int(tot_ch / 2)
                        end if
                        ! energies of two new photons
                        Ea_temp = tot_e / (2. * wa_temp)
                        Eb_temp = tot_e / (2. * wb_temp)
                        costh_temp = tot_p / tot_e
#ifdef PHOTDEBUG
                        ! for debug purposes
                        if (costh_temp .gt. 1.001) then
                            print *, 'ERROR: costh_temp .gt. 1 in merge'
                            print *, costh_temp, tot_p, tot_e
                            stop
                        end if
#endif
                        ! hack
                        if (costh_temp .gt. 1) then
                            costh_temp = 1.
                            sinth_temp = 0.
                        else
                            sinth_temp = sqrt(abs(1. - costh_temp**2))
                        end if

                        ! choosing d-vector
                        if (index_phi .eq. 1) then
                            phi_mid = bins_phi(index_phi) * 0.5
                        else
                            phi_mid = (bins_phi(index_phi) + bins_phi(index_phi - 1)) * 0.5
                        end if
                        if (index_theta .eq. 1) then
                            costheta_mid = bins_costheta(index_theta) * 0.5
                        else
                            costheta_mid = (bins_costheta(index_theta) + bins_costheta(index_theta - 1)) * 0.5
                        end if
#ifdef PHOTDEBUG
                        ! for debug purposes
                        if (costheta_mid .gt. 1.001) then
                            print *, 'ERROR: costheta_mid .gt. 1 in merge'
                            print *, costheta_mid, bins_costheta(index_theta)
                            stop
                        end if
                        if (phi_mid .gt. 2. * M_PI + 0.001) then
                            print *, 'ERROR: phi_mid .gt. 2pi in merge'
                            print *, phi_mid, bins_phi(index_phi)
                            stop
                        end if
#endif
                        ! hack
                        if (costheta_mid .gt. 1.) then
                            costheta_mid = 1.
                        end if
                        if (phi_mid .gt. 2. * M_PI) then
                            phi_mid = 2. * M_PI - 1e-4
                        end if

                        dx_temp = cos(phi_mid) * sqrt(abs(1. - costheta_mid**2))
                        dy_temp = sin(phi_mid) * sqrt(abs(1. - costheta_mid**2))
                        dz_temp = costheta_mid

                        e2x = tot_px / tot_p
                        e2y = tot_py / tot_p
                        e2z = tot_pz / tot_p
                        e2_dot_d = dx_temp * e2x + dy_temp * e2y + dz_temp * e2z
                        e1x = dx_temp - e2x * e2_dot_d
                        e1y = dy_temp - e2y * e2_dot_d
                        e1z = dz_temp - e2z * e2_dot_d
                        e1sqr = sqrt(e1x**2 + e1y**2 + e1z**2)
                        e1x = e1x / e1sqr
                        e1y = e1y / e1sqr
                        e1z = e1z / e1sqr

                        Pax_temp = Ea_temp * (costh_temp * e2x + sinth_temp * e1x)
                        Pay_temp = Ea_temp * (costh_temp * e2y + sinth_temp * e1y)
                        Paz_temp = Ea_temp * (costh_temp * e2z + sinth_temp * e1z)

                        Pbx_temp = Eb_temp * (costh_temp * e2x - sinth_temp * e1x)
                        Pby_temp = Eb_temp * (costh_temp * e2y - sinth_temp * e1y)
                        Pbz_temp = Eb_temp * (costh_temp * e2z - sinth_temp * e1z)

#ifdef PHOTDEBUG
                        ! for debug purposes
                        !   checking that the energy, momentum and weights are conserved
                        if (wa_temp + wb_temp .ne. tot_ch) then
                            print *, 'wa_temp + wb_temp .ne. tot_ch', wa_temp, wb_temp, tot_ch
                            stop
                        end if
                        ! if (abs(Ea_temp * wa_temp + Eb_temp * wb_temp - tot_e) / abs(tot_e) .gt. 1.0e-3) then
                        !     print *, 'Energy not conserved in merge', &
                        !             &Ea_temp * wa_temp, Eb_temp * wb_temp, tot_e
                        !     ! stop
                        ! end if
                        ! if (abs(Pax_temp * wa_temp + Pbx_temp * wb_temp - tot_px) / abs(tot_px) .gt. 1.0e-3) then
                        !     print *, 'Px not conserved in merge', &
                        !             &Pax_temp * wa_temp, Pbx_temp * wb_temp, tot_px
                        !     ! stop
                        ! end if
                        ! if (abs(Pay_temp * wa_temp + Pby_temp * wb_temp - tot_py) / abs(tot_py) .gt. 1.0e-3) then
                        !     print *, 'Py not conserved in merge', &
                        !             &Pay_temp * wa_temp, Pby_temp * wb_temp, tot_py
                        !     ! stop
                        ! end if
                        ! if (abs(Paz_temp * wa_temp + Pbz_temp * wb_temp - tot_pz) / abs(tot_pz) .gt. 1.0e-3) then
                        !     print *, 'Pz not conserved in merge', &
                        !             &Paz_temp * wa_temp, Pbz_temp * wb_temp, tot_pz
                        !     ! stop
                        ! end if
#endif
                        ! now putting those two new particles
                        index_1 = mom_bins(index_phi, index_theta, index_energy, 2)
                        index_2 = mom_bins(index_phi, index_theta, index_energy, 3)
#ifdef PHOTDEBUG
                        if (index_2 .eq. index_1) then
                            print *, 'ERROR: ', index_1, 'eq', index_2, 'in merging'
                            stop
                        end if
#endif
                        if (Ea_temp .ge. epsph_min) then
                            ph(index_1)%splitlev = -1
                            ph(index_1)%u = Pax_temp
                            ph(index_1)%v = Pay_temp
                            ph(index_1)%w = Paz_temp
                            ph(index_1)%ch = wa_temp
                        end if
                        if (Eb_temp .ge. epsph_min) then
                            ph(index_2)%splitlev = -1
                            ph(index_2)%u = Pbx_temp
                            ph(index_2)%v = Pby_temp
                            ph(index_2)%w = Pbz_temp
                            ph(index_2)%ch = wb_temp
                        end if
#ifdef PHOTDEBUG
                        ! FOR DEBUG PURPOSES
                        if (isnan(ph(index_1)%u + ph(index_1)%v + ph(index_1)%w + ph(index_2)%u + ph(index_2)%v + ph(index_2)%w)) then
                            print *, "ERROR: NAN found in merging VEL"
                            print *, index_1, index_2
                            print *, ph(index_1)%u, ph(index_1)%v, ph(index_1)%w, ph(index_2)%u, ph(index_2)%v, ph(index_2)%w
                            print *, ph(index_1)%x, ph(index_1)%y, ph(index_1)%z, ph(index_2)%x, ph(index_2)%y, ph(index_2)%z
                            print *, index_phi, index_theta, index_energy
                            print *, Ea_temp, Eb_temp, costh_temp, sinth_temp
                            print *, e1x, e1y, e1z, e2x, e2y, e2z
                            stop
                        end if
                        if (isnan(ph(index_1)%x + ph(index_1)%y + ph(index_1)%z + ph(index_2)%x + ph(index_2)%y + ph(index_2)%z)) then
                            print *, "ERROR: NAN found in merging CRD"
                            print *, index_1, index_2
                            print *, ph(index_1)%u, ph(index_1)%v, ph(index_1)%w, ph(index_2)%u, ph(index_2)%v, ph(index_2)%w
                            print *, ph(index_1)%x, ph(index_1)%y, ph(index_1)%z, ph(index_2)%x, ph(index_2)%y, ph(index_2)%z
                            print *, index_phi, index_theta, index_energy
                            print *, Ea_temp, Eb_temp, costh_temp, sinth_temp
                            print *, e1x, e1y, e1z, e2x, e2y, e2z
                            stop
                        end if
#endif
                        ! the rest of the photons will be deleted
                    end do
                end do
            end do ! end surfing through all the bins of a given cell
            !</editor-fold>
            if(allocated(mom_bins)) deallocate(mom_bins)
        end do ! end surfing through cells

        if(allocated(bins_phi)) deallocate(bins_phi)
        if(allocated(bins_costheta)) deallocate(bins_costheta)
        if(allocated(bins_energy_sq)) deallocate(bins_energy_sq)
    end subroutine merge_photons_conserv


    subroutine log_normal(n_bins, lognorm)
        real, parameter :: M_PI = 3.1415927

        integer, intent(in) :: n_bins
        real, allocatable, intent(inout) :: lognorm(:)
        real :: x, y, z, sum
        integer :: i

        allocate(lognorm(n_bins))
        sum = 0.
        do i = 1, n_bins
            x = random(dseed)
            y = random(dseed)
            z = sqrt(-2. * log(x)) * cos(2. * M_PI * y) ! now z has standard normal distribution
            z = exp(0. + 1. * z) ! now z has lognormal distribution with certain sigma=1 and mu=0
            lognorm(i) = z
            sum = sum + z
        end do
        ! this ensures lognorm(max) = 1/
        ! lognorm(1) = lognorm(1) / sum
        ! do i = 2, n_bins
        !     lognorm(i) = lognorm(i) / sum + lognorm(i - 1)
        ! end do
        ! /this ensures lognorm(max) = 1

        ! this allows having lognorm(max) != 1/
        !   in this case bins are not fixed in upper limit
        lognorm(1) = lognorm(1) / (n_bins + 1.)
        do i = 2, n_bins
            lognorm(i) = lognorm(i) / (n_bins + 1.) + lognorm(i - 1)
        end do
        ! /this allows having lognorm(max) != 1
    end subroutine log_normal

    subroutine find_bins(ind_ph, n_bins, bins_phi, bins_costheta, bins_energy_sq, index_phi, index_costheta, index_energy, ph_kick)
        real, parameter :: M_PI = 3.1415927
        integer, intent(in) :: ind_ph
        integer, intent(in) :: n_bins
        integer, intent(inout) :: index_phi, index_costheta, index_energy
        real :: phi_ph, theta_ph, costheta_ph, energy_ph_sq, uv_squared
        real, allocatable, intent(in) :: bins_phi(:), bins_costheta(:), bins_energy_sq(:)
        integer :: i
        real, intent(in) :: ph_kick

        ! phi_ph = atan(ph(ind_ph)%v / ph(ind_ph)%u)
        ! if ((ph(ind_ph)%u > 0) .and. (ph(ind_ph)%v < 0)) then
        !     phi_ph = 2.0d0 * M_PI + phi_ph
        ! end if
        ! if (ph(ind_ph)%u < 0) then
        !     phi_ph = M_PI + phi_ph
        ! end if
        if ((ph(ind_ph)%u .eq. 0.) .and. (ph(ind_ph)%v .eq. 0.)) then
            phi_ph = 2. * M_PI * random(dseed)
        else
            phi_ph = atan2(ph(ind_ph)%v, ph(ind_ph)%u) + M_PI
            phi_ph = phi_ph + ph_kick
        end if
        if (phi_ph .gt. 2.0 * M_PI) phi_ph = phi_ph - 2.0 * M_PI

        energy_ph_sq = (ph(ind_ph)%u)**2 + (ph(ind_ph)%v)**2 + (ph(ind_ph)%w)**2
        costheta_ph = ph(ind_ph)%w / sqrt(energy_ph_sq)

        ! hack here
        if ((costheta_ph .gt. 1.) .and. (costheta_ph .lt. 1. + 1e-5)) then
            costheta_ph = 1.
        end if
        if ((costheta_ph .lt. -1.) .and. (costheta_ph .gt. -1. - 1e-5)) then
            costheta_ph = -1.
        end if

#ifdef PHOTDEBUG
        ! FOR DEBUG PURPOSES
        if ((costheta_ph .gt. 1.) .or. (costheta_ph .lt. -1.)) then
            print *, "ERROR: costheta_ph error in find_bins() subroutine."
            print *, costheta_ph
            print *, ph(ind_ph)%u, ph(ind_ph)%v, ph(ind_ph)%w, energy_ph_sq
            stop
        end if
        if ((phi_ph .gt. 2. * M_PI) .or. (phi_ph .lt. 0.)) then
            print *, "ERROR: phi_ph error in find_bins() subroutine."
            print *, phi_ph, ph_kick
            print *, ph(ind_ph)%u, ph(ind_ph)%v, ph(ind_ph)%w, energy_ph_sq
            stop
        end if
#endif

        index_phi = -1
        index_costheta = -1
        index_energy = -1

        do i = 1, n_bins
            if ((phi_ph .le. bins_phi(i)) .and. (index_phi .eq. -1)) then
                index_phi = i
            end if
            if ((costheta_ph .le. bins_costheta(i)) .and. (index_costheta .eq. -1)) then
                index_costheta = i
            end if
            if ((energy_ph_sq .le. bins_energy_sq(i)) .and. (index_energy .eq. -1)) then
                index_energy = i
            end if
            if ((index_phi .ne. -1) .and. (index_costheta .ne. -1) .and. (index_energy .ne. -1)) then
                exit
            end if
        end do
#ifdef PHOTDEBUG
        ! FOR DEBUG PURPOSES
        if ((index_phi .eq. -1) .or. (index_costheta .eq. -1) .or. (index_energy .eq. -1)) then
            print *, "ERROR: index_* error in find_bins() subroutine."
            print *, sqrt(energy_ph_sq), costheta_ph, phi_ph
            print *, index_energy, index_costheta, index_phi
            do i = 1, n_bins
                print *, sqrt(bins_energy_sq(i)), bins_costheta(i), bins_phi(i)
            end do
            stop
        end if
#endif
    end subroutine find_bins

#endif ! PHOTONMERGE

#ifdef PAIRS

    subroutine check_interactions(c_num_photons)
        implicit none
        integer, allocatable, intent(in) :: c_num_photons(:)
        integer :: cell_num, ind_1, ind_2
        integer :: ind_temp, ind_temp2, temp_var1, temp_var2
        ! logical, parameter :: tmstop = .true.
        logical :: threshold

        real :: rand_temp, sum_temp
        real :: tauIJ

        !-------------------------------------------------------------------------------
        !   Looking for interactions
        !-------------------------------------------------------------------------------

        do cell_num = 1, mx * my * mz
            if (cell_num == 1) then
                if (c_num_photons(1) .lt. 2) then
                    cycle
                end if
                temp_var1 = 1
                temp_var2 = c_num_photons(1)
            else
                if (c_num_photons(cell_num) - c_num_photons(cell_num - 1) .lt. 2) then
                    cycle
                end if
                temp_var1 = c_num_photons(cell_num - 1) + 1
                temp_var2 = c_num_photons(cell_num)
            end if
            ! at this point temp_var1 contains the index of first photon in cell, and temp_var2 - the index of the last one
#ifdef PLASMOIDCONTROL
            if (loc_density(int(ph(temp_var1)%x), int(ph(temp_var1)%y)) .gt. bw_dens_lim) then
                cycle
            end if
#endif
            ! the following algorithm uses binary probabilities for the BW photon-pair interaction
            do ind_1 = temp_var1, temp_var2 - 1 ! loop through all the particles in a given cell
                if ((ph(ind_1)%splitlev .eq. 1) .or. (ph(ind_1)%ch .le. 0)) then ! photon is scheduled for deletion
                    cycle
                end if
                do ind_2 = ind_1 + 1, temp_var2
                    if ((ph(ind_2)%splitlev .eq. 1) .or. (ph(ind_2)%ch .le. 0)) then ! photon is scheduled for deletion
                        cycle
                    end if

                    threshold = .false.
                    call setProbability(ind_1, ind_2, threshold, tauIJ)

                    if (threshold) then ! if the momenta satisfy BW limit condition
                        do while ((tauIJ > 0.0) .and. (ph(ind_1)%ch .gt. 0) .and. (ph(ind_2)%ch .gt. 0))
                            if (tauIJ < 1.0) then
                                rand_temp = random(dseed)
                                if (rand_temp <= tauIJ) then
                                    ph(ind_1)%ch = ph(ind_1)%ch - 1
                                    ph(ind_2)%ch = ph(ind_2)%ch - 1
                                    call create_pairs(ind_1, ind_2, threshold)
                                end if
                                exit
                            else
                                ph(ind_1)%ch = ph(ind_1)%ch - 1
                                ph(ind_2)%ch = ph(ind_2)%ch - 1
                                call create_pairs(ind_1, ind_2, threshold)
                                tauIJ = tauIJ - 1.
                            end if
                        end do
                    end if ! threshold
                    if (ph(ind_1)%ch .le. 0) then
                        ph(ind_1)%splitlev = 1
                    end if
                    if (ph(ind_2)%ch .le. 0) then
                        ph(ind_2)%splitlev = 1
                    end if
                end do
            end do
        end do
        !-------------------------------------------------------------------------------
        !   / Looking for interactions
        !-------------------------------------------------------------------------------
    end subroutine check_interactions

!-------------------------------------------------------------------------------
!
!   Auxiliary functions
!
!-------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------
    !   Calculating probability per unit timestep for two given photons
    !-------------------------------------------------------------------------------
    subroutine setProbability(ph_i, ph_j, threshold, tau_ij)
        integer, intent(in) :: ph_i, ph_j
        logical, intent(inout) :: threshold
        real, intent(inout) :: tau_ij
        real(dprec) :: k1_x, k1_y, k1_z
        real(dprec) :: k2_x, k2_y, k2_z

        real(dprec) :: cosphi, s, beta, beta2, fs
        real(dprec) :: E1, E2 ! photon energies
        real(dprec) :: ph_u1, ph_v1, ph_w1, ph_u2, ph_v2, ph_w2

        ph_u1 = real(ph(ph_i)%u, 8)
        ph_v1 = real(ph(ph_i)%v, 8)
        ph_w1 = real(ph(ph_i)%w, 8)
        ph_u2 = real(ph(ph_j)%u, 8)
        ph_v2 = real(ph(ph_j)%v, 8)
        ph_w2 = real(ph(ph_j)%w, 8)

        E1 = sqrt((ph_u1)**2 + (ph_v1)**2 + (ph_w1)**2)
        E2 = sqrt((ph_u2)**2 + (ph_v2)**2 + (ph_w2)**2)

        if (E1 * E2 .lt. 1.) then
            threshold = .false.
        else
            ! photon k-vectors
            k1_x = ph_u1 / E1
            k1_y = ph_v1 / E1
            k1_z = ph_w1 / E1
            k2_x = ph_u2 / E2
            k2_y = ph_v2 / E2
            k2_z = ph_w2 / E2
            cosphi = k1_x * k2_x + k1_y * k2_y + k1_z * k2_z
            s = E1 * E2 * (1. - cosphi) / 2.
            threshold = (s .gt. 1.0000001)
        end if
        if (threshold) then
            beta2 = 1. - 1. / s
            beta = sqrt(beta2)
            fs = (1. - beta2) * (-2. * beta * (2. - beta2) + (3. - beta2**2) * log((1. + beta) / (1. - beta)))
            tau_ij = ph_step * tau0 * fs * ph(ph_i)%ch * ph(ph_j)%ch ! to be set later to "real" value
#ifdef PAIRSDEBUG
            if (isnan(beta)) then
                print *, "ERROR: Something went wrong in setProbability()"
                print *, beta2, s, E1, E2, cosphi, threshold
                print *, ph(ph_i)%u, ph(ph_i)%v, ph(ph_i)%w, ph(ph_i)%ch
                print *, ph(ph_j)%u, ph(ph_j)%v, ph(ph_j)%w, ph(ph_j)%ch
                stop
            end if
            if (tau_ij .lt. 0.) then
                print *, "ERROR: tau_ij < 0 in setProbability"
                print *, beta2, s, E1, E2, cosphi, threshold
                print *, ph_step, tau0, fs
                print *, ph(ph_i)%u, ph(ph_i)%v, ph(ph_i)%w, ph(ph_i)%ch
                print *, ph(ph_j)%u, ph(ph_j)%v, ph(ph_j)%w, ph(ph_j)%ch
                stop
            end if
#endif
        else
            tau_ij = 0.
        end if
    end subroutine
#endif ! PAIRS
    !-------------------------------------------------------------------------------
    !   Deleting all the scheduled photons
    !-------------------------------------------------------------------------------
    subroutine clean_photons_up()
        implicit none
        integer :: ind_ph
        ind_ph = 1
        do while (ind_ph <= photons)
#ifdef PHOTDEBUG
            ! FOR DEBUG PURPOSES
            if ((ph(ind_ph)%splitlev .ne. 1) .and. (ph(ind_ph)%ch .le. 0)) then
                print *, "ERROR: error in clean_photons_up()."
                stop
            end if
#endif
            if ((ph(ind_ph)%splitlev .eq. 1) .or. (ph(ind_ph)%ch .le. 0)) then ! this means, the particle is scheduled for deletion
                call copyprt(ph(photons), ph(ind_ph))
                photons = photons - 1
            else
                ind_ph = ind_ph + 1
            end if
        end do
#ifdef PHOTDEBUG
        ! FOR DEBUG PURPOSES
        do ind_ph = 1, photons
            if (ph(ind_ph)%splitlev .eq. 1) then
                print *, "ERROR: photon with splitlev == 1 found in clean_photons_up"
                print *, ph(ind_ph)%x, ph(ind_ph)%y, ph(ind_ph)%z
                print *, ph(ind_ph)%ch, photons, ind_ph
                stop
            end if
            if (ph(ind_ph)%ch .le. 0) then
                print *, "ERROR: photon with ch <= 0 found in clean_photons_up"
                print *, ph(ind_ph)%x, ph(ind_ph)%y, ph(ind_ph)%z
                print *, ph(ind_ph)%ch, photons, ind_ph
                stop
            end if
        end do
#endif
    end subroutine clean_photons_up

#ifdef PAIRS
    !-------------------------------------------------------------------------------
    !   This function takes as arguments indices of photons, and create pairs from them
    !-------------------------------------------------------------------------------
    subroutine create_pairs(ph_i, ph_j, threshold)
        implicit none
        logical, intent(in) :: threshold
        integer, intent(in) :: ph_i, ph_j
        integer :: iter

        ! real, parameter :: M_PI = 3.1415927
        real(dprec), parameter :: M_PI = 3.1415926535897931d0

        real(dprec) :: rand_prob, rand_theta, rand_phi
        real(dprec) :: s, cosphi

        real(dprec) :: E1, E2 ! photon energy (m_e c^2)
        real(dprec) :: k1_x, k1_y, k1_z
        real(dprec) :: k2_x, k2_y, k2_z
        real(dprec) :: beta_x, beta_y, beta_z, beta_cm, gamma_cm

        ! auxilary vectors, k = center-of-momentum common k-vector
        ! a is perp to k in x-y plane
        ! b is perp to k and a
        real(dprec) :: k_vec, a_vec, b_vec
        real(dprec) :: k_vec_x, a_vec_x, b_vec_x
        real(dprec) :: k_vec_y, a_vec_y, b_vec_y
        real(dprec) :: k_vec_z, a_vec_z, b_vec_z
        real(dprec) :: p_el_x, p_el_y, p_el_z
        real(dprec) :: beta_el, gamma_el ! primarily in the CM frame
        real(dprec) :: u_el_x, u_el_y, u_el_z ! 3d 4-velocity in CM frame
        real(dprec) :: v_el_x, v_el_y, v_el_z ! 3d 4-velocity in normal frame

        real(dprec) :: cos_rand_theta, sin_rand_theta, cos_rand_phi, sin_rand_phi

        ! error handling variables
        real(dprec) :: er_k_vec_x, er_k_vec_y, er_k_vec_z, er_k_vec
        real(dprec) :: ph_u1, ph_v1, ph_w1, ph_u2, ph_v2, ph_w2

        ph_u1 = real(ph(ph_i)%u, 8)
        ph_v1 = real(ph(ph_i)%v, 8)
        ph_w1 = real(ph(ph_i)%w, 8)
        ph_u2 = real(ph(ph_j)%u, 8)
        ph_v2 = real(ph(ph_j)%v, 8)
        ph_w2 = real(ph(ph_j)%w, 8)

#ifdef PAIRSDEBUG
        if ((ph(ph_i)%ch .lt. 0) .or. (ph(ph_j)%ch .lt. 0)) then
            print *, "ERROR: Something went wrong for photons %ch."
            print *, ph_u1, ph_v1, ph_w1, ph(ph_i)%ch
            print *, ph_u2, ph_v2, ph_w2, ph(ph_j)%ch
            stop
        end if
#endif

        E1 = sqrt((ph_u1)**2 + (ph_v1)**2 + (ph_w1)**2)
        E2 = sqrt((ph_u2)**2 + (ph_v2)**2 + (ph_w2)**2)

        ! photon k-vectors
        k1_x = ph_u1 / E1
        k1_y = ph_v1 / E1
        k1_z = ph_w1 / E1
        k2_x = ph_u2 / E2
        k2_y = ph_v2 / E2
        k2_z = ph_w2 / E2
        cosphi = k1_x * k2_x + k1_y * k2_y + k1_z * k2_z
        ! square of energy for each of the photons in the CM frame (in m_e c^2, so it's the square of gamma-factor)
        s = E1 * E2 * (1. - cosphi) / 2.

#ifdef PAIRSDEBUG
        ! FOR DEBUG PURPOSES
        if (s <= 1.) then
            print *, "ERROR: Something went wrong in create_pairs() subroutine. s < 1."
            print *, threshold
            print *, s, cosphi, E1, E2
            print *, ph_u1, ph_v1, ph_w1, ph(ph_i)%ch
            print *, ph_u2, ph_v2, ph_w2, ph(ph_j)%ch
            stop
        end if
#endif

        beta_x = (E1 * k1_x + E2 * k2_x) / (E1 + E2)
        beta_y = (E1 * k1_y + E2 * k2_y) / (E1 + E2)
        beta_z = (E1 * k1_z + E2 * k2_z) / (E1 + E2)
        beta_cm = beta_x**2 + beta_y**2 + beta_z**2

#ifdef PAIRSDEBUG
        ! FOR DEBUG PURPOSES
        if (beta_cm > 1.) then
            print *, "ERROR: beta_cm > 1 in create_pairs() subroutine."
            print *, beta_cm, beta_x, beta_y, beta_z
            print *, s, cosphi, E1, E2
            print *, ph_u1, ph_v1, ph_w1, ph(ph_i)%ch
            print *, ph_u2, ph_v2, ph_w2, ph(ph_j)%ch
            ! stop
        end if
#endif
        ! hack here
        if (beta_cm .gt. 1.) then
            beta_cm = (2. - sqrt(beta_cm))**2
            print *, 'ERROR: Hack initiated in create_pairs(), not too happy with that. :('
        end if

        gamma_cm = 1. / sqrt(1. - beta_cm)
        beta_cm = sqrt(beta_cm)

#ifdef PAIRSDEBUG
        ! FOR DEBUG PURPOSES
        if (isnan(gamma_cm)) then
            print *, gamma_cm, beta_cm, 1. / sqrt(1. - beta_cm), 1. - beta_cm
            print *, beta_x, beta_y, beta_z
            stop
        endif
#endif
        !-------------------------------------------------------------------------------
        !   Forward Lorentz-transformation for photons
        !-------------------------------------------------------------------------------

        if (beta_cm .gt. 0.) then
            ! this is basically the k-vector of the first photon in CM frame
            k_vec_x = -beta_x * gamma_cm + &
                    &k1_x * (1. - (beta_x**2 / beta_cm**2) * (1. - gamma_cm)) + &
                    &k1_y * (gamma_cm - 1.) * beta_x * beta_y / beta_cm**2 + &
                    &k1_z * (gamma_cm - 1.) * beta_x * beta_z / beta_cm**2

            k_vec_y = -beta_y * gamma_cm +&
                    &k1_y * (1. - (beta_y**2 / beta_cm**2) * (1. - gamma_cm)) + &
                    &k1_x * (gamma_cm - 1.) * beta_y * beta_x / beta_cm**2 + &
                    &k1_z * (gamma_cm - 1.) * beta_y * beta_z / beta_cm**2

            k_vec_z = -beta_z * gamma_cm + &
                    &k1_z * (1. - (beta_z**2 / beta_cm**2) * (1. - gamma_cm)) + &
                    &k1_x * (gamma_cm - 1.) * beta_z * beta_x / beta_cm**2 + &
                    &k1_y * (gamma_cm - 1.) * beta_z * beta_y / beta_cm**2
        else
            k_vec_x = k1_x
            k_vec_y = k1_y
            k_vec_z = k1_z
        end if

        k_vec = sqrt(k_vec_x**2 + k_vec_y**2 + k_vec_z**2)
        k_vec_x = k_vec_x / k_vec
        k_vec_y = k_vec_y / k_vec
        k_vec_z = k_vec_z / k_vec

#ifdef PAIRSDEBUG
        ! FOR DEBUG PURPOSES
        ! error handling k-vector of the second photon in CM frame (has to be equal to the first one with opposite sign)
        er_k_vec_x = -beta_x * gamma_cm + k2_x * (1. - (beta_x**2 / beta_cm**2) * (1. - gamma_cm)) &
                &+ k2_y * (gamma_cm - 1.) * beta_x * beta_y / beta_cm**2 &
                &+ k2_z * (gamma_cm - 1.) * beta_x * beta_z / beta_cm**2
        er_k_vec_y = -beta_y * gamma_cm + k2_y * (1. - (beta_y**2 / beta_cm**2) * (1. - gamma_cm)) &
                &+ k2_x * (gamma_cm - 1.) * beta_y * beta_x / beta_cm**2 &
                &+ k2_z * (gamma_cm - 1.) * beta_y * beta_z / beta_cm**2
        er_k_vec_z = -beta_z * gamma_cm + k2_z * (1. - (beta_z**2 / beta_cm**2) * (1. - gamma_cm)) &
                &+ k2_x * (gamma_cm - 1.) * beta_z * beta_x / beta_cm**2 &
                &+ k2_y * (gamma_cm - 1.) * beta_z * beta_y / beta_cm**2

        er_k_vec = sqrt(er_k_vec_x**2 + er_k_vec_y**2 + er_k_vec_z**2)
        er_k_vec_x = er_k_vec_x / er_k_vec
        er_k_vec_y = er_k_vec_y / er_k_vec
        er_k_vec_z = er_k_vec_z / er_k_vec

        if (sqrt((k_vec_x + er_k_vec_x)**2 + (k_vec_y + er_k_vec_y)**2 + (k_vec_z + er_k_vec_z)**2) .gt. 1e-3) then
            print *, "k1:", k1_x, k1_y, k1_z
            print *, "k2:", k2_x, k2_y, k2_z
            print *, "k1':", k_vec_x, k_vec_y, k_vec_z
            print *, "k2':", er_k_vec_x, er_k_vec_y, er_k_vec_z
            print *, 'error is:', sqrt((k_vec_x + er_k_vec_x)**2 + (k_vec_y + er_k_vec_y)**2 + (k_vec_z + er_k_vec_z)**2)
            print *, "Problem in create_pairs() subroutine! Lorentz-transformation done wrong. Error > 1e-3."
            ! stop
        end if
#endif
        ! now since everything is fine, we can do stuff in CM frame

        if (abs(k_vec_x) .gt. 1e-6) then
            a_vec_z = 0.
            a_vec_y = 1.
            a_vec_x = -k_vec_y / k_vec_x
            a_vec = sqrt(a_vec_x**2 + a_vec_y**2)
            a_vec_x = a_vec_x / a_vec
            a_vec_y = a_vec_y / a_vec
        else
            a_vec_z = 0.
            a_vec_y = 0.
            a_vec_x = 1.
        end if
        b_vec_x = a_vec_z * k_vec_y - a_vec_y * k_vec_z
        b_vec_y = -a_vec_z * k_vec_x + a_vec_x * k_vec_z
        b_vec_z = a_vec_y * k_vec_x - a_vec_x * k_vec_y

        b_vec = sqrt(b_vec_x**2 + b_vec_y**2 + b_vec_z**2)
        b_vec_x = b_vec_x / b_vec
        b_vec_y = b_vec_y / b_vec
        b_vec_z = b_vec_z / b_vec

        !-------------------------------------------------------------------------------
        !   Random electron/positron vector in CM frame
        !-------------------------------------------------------------------------------

        iter = 0
        do while (.true.)
            rand_prob = random(dseed)
            rand_theta = M_PI * random(dseed)
#ifdef PAIRSDEBUG
            ! FOR DEBUG PURPOSES
            if (isnan(dsigma_breit_wheeler(s, rand_theta))) then
                print *, 'ERROR: Something wrong in create_pairs() subroutine.'
                print *, E1, E2, cosphi
                print *, iter, s, dsigma_breit_wheeler(s, rand_theta), rand_prob, rand_theta
                stop
            end if
#endif
            if(rand_prob <= dsigma_breit_wheeler(s, rand_theta)) then
                exit
            end if
            iter = iter + 1

            ! FOR DEBUG PURPOSES
            if (iter > 10000) then
                print *, "ERROR: Cannot exit loop in create_pairs() subroutine."
                print *, E1, E2, cosphi
                print *, iter, s, dsigma_breit_wheeler(s, rand_theta), rand_prob, rand_theta
                stop
            end if
        end do
        rand_phi = 2. * M_PI * random(dseed)
        cos_rand_theta = cos(rand_theta)
        cos_rand_phi = cos(rand_phi)
        sin_rand_theta = sin(rand_theta)
        sin_rand_phi = sin(rand_phi)

        ! at this point we have a rand_phi angle (0, 2pi) in the perpendicular plane
        ! and rand_theta from 0 to pi scattering angle
        ! note: the angles are in center-of-momentum frame

        p_el_x = k_vec_x * cos_rand_theta + a_vec_x * sin_rand_theta * cos_rand_phi + b_vec_x * sin_rand_theta * sin_rand_phi
        p_el_y = k_vec_y * cos_rand_theta + a_vec_y * sin_rand_theta * cos_rand_phi + b_vec_y * sin_rand_theta * sin_rand_phi
        p_el_z = k_vec_z * cos_rand_theta + a_vec_z * sin_rand_theta * cos_rand_phi + b_vec_z * sin_rand_theta * sin_rand_phi

        gamma_el = sqrt(s)
        beta_el = sqrt(1. - 1. / s)

        !-------------------------------------------------------------------------------
        !   Inverse Lorentz-transformation for e/p 4-velocities
        !-------------------------------------------------------------------------------
        beta_x = -beta_x
        beta_y = -beta_y
        beta_z = -beta_z

        ! creating positron
        u_el_x = gamma_el * beta_el * p_el_x
        u_el_y = gamma_el * beta_el * p_el_y
        u_el_z = gamma_el * beta_el * p_el_z

        if (beta_cm .gt. 0.) then
            v_el_x = -beta_x * gamma_cm * gamma_el + &
                &u_el_x * (1. + beta_x**2 * (gamma_cm - 1.) / beta_cm**2) + &
                &u_el_y * beta_x * beta_y * (gamma_cm - 1.) / beta_cm**2 + &
                &u_el_z * beta_x * beta_z * (gamma_cm - 1.) / beta_cm**2
            v_el_y = -beta_y * gamma_cm * gamma_el + &
                &u_el_y * (1. + beta_y**2 * (gamma_cm - 1.) / beta_cm**2) + &
                &u_el_x * beta_y * beta_x * (gamma_cm - 1.) / beta_cm**2 + &
                &u_el_z * beta_y * beta_z * (gamma_cm - 1.) / beta_cm**2
            v_el_z = -beta_z * gamma_cm * gamma_el + &
                &u_el_z * (1. + beta_z**2 * (gamma_cm - 1.) / beta_cm**2) + &
                &u_el_x * beta_z * beta_x * (gamma_cm - 1.) / beta_cm**2 + &
                &u_el_y * beta_z * beta_y * (gamma_cm - 1.) / beta_cm**2
        else
            v_el_x = u_el_x
            v_el_y = u_el_y
            v_el_z = u_el_z
        end if

        ions = ions + 1
        p(ions)%proc = rank
        totalpartnum = totalpartnum + 1
        p(ions)%ind = -2 * abs(totalpartnum)
        p(ions)%ch = 1.0 ! set to actual charge (FIX)
        p(ions)%x = (ph(ph_i)%x + ph(ph_j)%x) / 2.
        p(ions)%y = (ph(ph_i)%y + ph(ph_j)%y) / 2.
        p(ions)%z = (ph(ph_i)%z + ph(ph_j)%z) / 2.
        p(ions)%u = v_el_x
        p(ions)%v = v_el_y
        p(ions)%w = v_el_z

        ! creating electron
        u_el_x = -u_el_x
        u_el_y = -u_el_y
        u_el_z = -u_el_z

        if (beta_cm .gt. 0.) then
            v_el_x = -beta_x * gamma_cm * gamma_el + &
                &u_el_x * (1. + beta_x**2 * (gamma_cm - 1.) / beta_cm**2) + &
                &u_el_y * beta_x * beta_y * (gamma_cm - 1.) / beta_cm**2 + &
                &u_el_z * beta_x * beta_z * (gamma_cm - 1.) / beta_cm**2
            v_el_y = -beta_y * gamma_cm * gamma_el + &
                &u_el_y * (1. + beta_y**2 * (gamma_cm - 1.) / beta_cm**2) + &
                &u_el_x * beta_y * beta_x * (gamma_cm - 1.) / beta_cm**2 + &
                &u_el_z * beta_y * beta_z * (gamma_cm - 1.) / beta_cm**2
            v_el_z = -beta_z * gamma_cm * gamma_el + &
                &u_el_z * (1. + beta_z**2 * (gamma_cm - 1.) / beta_cm**2) + &
                &u_el_x * beta_z * beta_x * (gamma_cm - 1.) / beta_cm**2 + &
                &u_el_y * beta_z * beta_y * (gamma_cm - 1.) / beta_cm**2
        else
            v_el_x = u_el_x
            v_el_y = u_el_y
            v_el_z = u_el_z
        end if

        lecs = lecs + 1
        p(maxhlf + lecs)%proc = rank
        totalpartnum = totalpartnum + 1
        p(maxhlf + lecs)%ind = -2 * abs(totalpartnum)
        p(maxhlf + lecs)%ch = 1. ! set to actual charge (FIX)
        p(maxhlf + lecs)%x = (ph(ph_i)%x + ph(ph_j)%x) / 2.
        p(maxhlf + lecs)%y = (ph(ph_i)%y + ph(ph_j)%y) / 2.
        p(maxhlf + lecs)%z = (ph(ph_i)%z + ph(ph_j)%z) / 2.
        p(maxhlf + lecs)%u = v_el_x
        p(maxhlf + lecs)%v = v_el_y
        p(maxhlf + lecs)%w = v_el_z

#ifdef OUTSTATS
        dataE1E2cosphi(dataE1E2_cntr + 1, 1) = E1
        dataE1E2cosphi(dataE1E2_cntr + 1, 2) = E2
        dataE1E2cosphi(dataE1E2_cntr + 1, 3) = cosphi
        dataE1E2cosphi(dataE1E2_cntr + 1, 4) = 1. + v_el_x**2 + v_el_y**2 + v_el_z**2
        dataE1E2cosphi(dataE1E2_cntr + 1, 5) = ph(ph_i)%y + mycum
        dataE1E2_cntr = dataE1E2_cntr + 1
#endif
#ifdef FINDNAN
		if (isnan(p(maxhlf + lecs)%x + p(maxhlf + lecs)%y + p(maxhlf + lecs)%z +&
            & p(ions)%x + p(ions)%y + p(ions)%z)) then
			print *, 'ERROR: COORD problem in create_pairs'
		endif
		if (isnan(p(maxhlf + lecs)%u + p(maxhlf + lecs)%v + p(maxhlf + lecs)%w +&
            & p(ions)%u + p(ions)%v + p(ions)%w)) then
			print *, 'ERROR: VEL problem in create_pairs'
			print *, rank
            print *, s, cosphi, E1, E2, gamma_cm
            print *, gamma_cm, beta_cm, 1. / sqrt(1. - beta_cm), 1. - beta_cm
            print *, ph(ph_i)%u, ph(ph_i)%v, ph(ph_i)%w, ph(ph_i)%ch
            print *, ph(ph_j)%u, ph(ph_j)%v, ph(ph_j)%w, ph(ph_j)%ch
            print *, p(maxhlf + lecs)%u, p(maxhlf + lecs)%v, p(maxhlf + lecs)%w
            print *, p(ions)%u, p(ions)%v, p(ions)%w
            print *, rand_phi, rand_theta, iter
            print *, u_el_x, u_el_y, u_el_z, gamma_el
            print *, v_el_x, v_el_y, v_el_z
            stop
		endif
#endif
        ! if (mod(totalpartnum, 1000) == 0) then
        !     print *, "PRINTING BW PARAMETERS /"
        !     print *, "k1:", k1_x, k1_y, k1_z
        !     print *, "k2:", k2_x, k2_y, k2_z
        !     print *, "beta:", beta_cm, "gamma:", gamma_cm
        !     print *, "k:", k_vec_x, k_vec_y, k_vec_z
        !     print *, "s:", s, "theta:", rand_theta, "phi:", rand_phi
        !     print *, "io:", p(ions)%u, p(ions)%v, p(ions)%w
        !     print *, "el:", p(maxhlf + lecs)%u, p(maxhlf + lecs)%v, p(maxhlf + lecs)%w
        !     print *, "/ END PRINTING BW PARAMETERS"
        ! end if
    end subroutine create_pairs

    real function dsigma_breit_wheeler(s, theta)
        real(dprec), intent(in) :: s, theta
        real(dprec) :: beta, beta2, beta4
        beta2 = 1. - 1. / s
        beta4 = beta2**2
        beta = sqrt(beta2)
        dsigma_breit_wheeler = (beta * (-8. - 8. * beta2 + 11. * beta4 - &
            &4. * beta2 * (-2.0 + beta2) * cos(2. * theta) + &
            &beta4 * cos(4. * theta)) * sin(theta)) / &
            &(4. * (-(beta * (-2. + beta2)) + (-3. + beta4) * asinh(beta / sqrt(1. - beta2))) * &
            &(-2. + beta2 + beta2 * cos(2. * theta))**2)
    end function dsigma_breit_wheeler
#endif ! PAIRS
!-------------------------------------------------------------------------------
!
!   / Auxilary functions
!
!-------------------------------------------------------------------------------

#ifdef twoD
end module m_photon_interact
#else
end module m_photon_interact_3d
#endif ! twoD
#endif ! photonmerge || pairs
#endif ! PHOTONS
