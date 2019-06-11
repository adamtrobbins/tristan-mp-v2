#include "../defs.F90"

module m_thermalplasma
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_particles
  use m_fields
  use m_particlelogistics
  implicit none

  logical :: maxwell_generated = .false.

  type :: maxwellian
    ! tabulated maxwellian is only used with `T < 0.1 [me c^2]`
    !   so non-relativistic maxwellian is used with the parameter
    !     `beta_wave = beta / sqrt(T)` for simplicity
    real                            :: temperature, shift_gamma
    real, allocatable, dimension(:) :: DF_table
    real, allocatable, dimension(:) :: beta_table
    logical                         :: generated, shift_flag
    integer                         :: npoints, shift_dir
  end type maxwellian

  !--- PRIVATE functions -----------------------------------------!
  private :: tabulateMaxwellian, generateFromMaxwellian, deallocateMaxwellian
  !...............................................................!
contains
  ! See more details in Zenitani 2015
  !   arXiv:1504.03910v1

  ! this should only be used for `T << 1`
  subroutine tabulateMaxwellian(maxw, n)
    implicit none
    type(maxwellian), intent(inout) :: maxw
    integer, intent(in)     :: n
    integer                 :: iter
    real                    :: beta_wave1, beta_wave2, beta_wave_max, df
    ! `beta_wave` is `beta / sqrt(T)`

    if (maxw%generated) then
      call throwError('ERROR: maxwell table already generated.')
    else
      maxw%generated = .true.
    end if

    maxw%npoints = n
    allocate(maxw%DF_table(maxw%npoints))
    allocate(maxw%beta_table(maxw%npoints))

    beta_wave_max = 0.9 / sqrt(maxw%temperature)
    do iter = 1, maxw%npoints
      beta_wave1 = beta_wave_max * REAL(iter - 1) / REAL(maxw%npoints)
      beta_wave2 = beta_wave_max * REAL(iter) / REAL(maxw%npoints)
      df = (beta_wave2 - beta_wave1) * 0.5 *&
            & (exp(-beta_wave1**2 * 0.5) * beta_wave1**2 + exp(-beta_wave2**2 * 0.5) * beta_wave2**2)
      maxw%beta_table(iter) = beta_wave2
      if (iter .eq. 1) then
        maxw%DF_table(iter) = df
      else
        maxw%DF_table(iter) = maxw%DF_table(iter - 1) + df
      end if
    end do
    do iter = 1, maxw%npoints
      maxw%DF_table(iter) = maxw%DF_table(iter) / maxw%DF_table(maxw%npoints)
    end do
  end subroutine tabulateMaxwellian

  subroutine generateFromMaxwellian(maxw, u_, v_, w_)
    implicit none
    type(maxwellian), intent(in) :: maxw
    real, intent(out)            :: u_, v_, w_
    real                         :: U, ETA, X1, X2, X3, X4, X5, X6, X7, X8, dx1, dx2, BETA, gamma, gamma1
    logical                      :: flag1, flag2
    integer                      :: iter
    flag1 = .false.
    do while (.not. flag1)
      if (maxw%temperature .lt. 0.1) then
        ! using tabulated Maxwellian
        if (.not. maxw%generated) then
          call throwError('ERROR: maxwell table not generated yet.')
        end if
        X3 = random(dseed)
        do iter = 1, maxw%npoints
          if (maxw%DF_table(iter) .ge. X3) then
            if (iter .gt. 1) then
              dx1 = (maxw%DF_table(iter) - X3) / (maxw%DF_table(iter) - maxw%DF_table(iter - 1))
              dx2 = (X3 - maxw%DF_table(iter - 1)) / (maxw%DF_table(iter) - maxw%DF_table(iter - 1))
              U = maxw%beta_table(iter) * dx2 +&
                & maxw%beta_table(iter - 1) * dx1
            else
              dx2 = X3 / maxw%DF_table(iter)
              U = maxw%beta_table(iter) * dx2
            end if
            ! `U` now is in terms of `beta_wave = beta / sqrt(T)`
            !   step 1 (convert to `beta`):
            U = U * sqrt(maxw%temperature)
            !   step 2 (convert to 4-velocity):
            U = U / sqrt(1.0 - U**2)
            exit
          end if
        end do
      else
        ! using Sobol method
        flag2 = .false.
        do while (.not. flag2)
          X4 = random(dseed); X5 = random(dseed)
          X6 = random(dseed); X7 = random(dseed)
          if (X4 * X5 * X6 * X7 .eq. 0) cycle
          U = -maxw%temperature * log(X4 * X5 * X6)
          ETA = -maxw%temperature * log(X4 * X5 * X6 * X7)
          if (ETA**2 - U**2 .gt. 1) then
            flag2 = .true.
          end if
        end do
      end if

      ! generate projections
      X1 = random(dseed); X2 = random(dseed)
      u_ = U * (2.0 * X1 - 1.0)
      v_ = 2.0 * U * sqrt(X1 * (1.0 - X1)) * cos(2.0 * M_PI * X2)
      w_ = 2.0 * U * sqrt(X1 * (1.0 - X1)) * sin(2.0 * M_PI * X2)

      ! shift maxwellian
      if (maxw%shift_flag) then
        gamma = sqrt(1.0 + U**2)
        BETA = sqrt(1.0 - 1.0 / maxw%shift_gamma**2)
        if (maxw%shift_dir .eq. 1) then ! in +x
          u_ = maxw%shift_gamma * (u_ + BETA * gamma)
        else if (maxw%shift_dir .eq. -1) then ! in -x
          u_ = maxw%shift_gamma * (u_ - BETA * gamma)
        else if (maxw%shift_dir .eq. 2) then ! in +y
          v_ = maxw%shift_gamma * (v_ + BETA * gamma)
        else if (maxw%shift_dir .eq. -2) then ! in -y
          v_ = maxw%shift_gamma * (v_ - BETA * gamma)
        else if (maxw%shift_dir .eq. 3) then ! in +z
          w_ = maxw%shift_gamma * (w_ + BETA * gamma)
        else if (maxw%shift_dir .eq. -3) then ! in -z
          w_ = maxw%shift_gamma * (w_ - BETA * gamma)
        else
          call throwError('ERROR: unknown shift directions of maxwellian.')
        end if
        X8 = random(dseed)
        gamma1 = sqrt(1.0 + u_**2 + v_**2 + w_**2)
        if (0.5 * (gamma1 / gamma) / maxw%shift_gamma .gt. X8) then
          flag1 = .true.
        end if
      else
        flag1 = .true.
      end if
    end do
  end subroutine generateFromMaxwellian

  subroutine deallocateMaxwellian(maxw)
    implicit none
    type(maxwellian), intent(inout) :: maxw
    if (allocated(maxw%DF_table)) deallocate(maxw%DF_table)
    if (allocated(maxw%beta_table)) deallocate(maxw%beta_table)
  end subroutine deallocateMaxwellian

  subroutine fillRegionWithThermalPlasma(fillregion, fill_species, num_species, num_part,&
                                       & temperature, shift_gamma, shift_dir,&
                                       & spat_distr_ptr,&
                                       & dummy1, dummy2, dummy3)
    implicit none
    ! assuming that the charges of all species given in `fill_species` add up to `0`
    type(region), intent(in)         :: fillregion
    integer, intent(in)              :: num_part, num_species
    integer, intent(in)              :: fill_species(num_species)
    real, intent(in)                 :: temperature
    real, optional, intent(in)       :: shift_gamma
    integer, optional, intent(in)    :: shift_dir
    type(maxwellian)                 :: fill_maxwellian
    integer                          :: n, s, spec_
    integer(kind=2)                  :: xi_, yi_, zi_
    real                             :: u_, v_, w_, dx_, dy_, dz_
    real                             :: x_, y_, z_, rnd
    real                             :: x_glob

    procedure (spatialDistribution), pointer, intent(in), optional :: spat_distr_ptr
    real, intent(in), optional                                     :: dummy1, dummy2, dummy3

    fill_maxwellian%temperature = temperature
    fill_maxwellian%generated = .false.
    if (present(shift_gamma)) then
      fill_maxwellian%shift_gamma = shift_gamma
      fill_maxwellian%shift_flag = .true.
    else
      fill_maxwellian%shift_flag = .false.
    end if
    if (temperature .lt. 0.1) then
      call tabulateMaxwellian(fill_maxwellian, 2000)
    end if

    n = 0
    do while (n .lt. num_part)
      ! generate coords for all species
      rnd = random(dseed)
      x_ = fillregion%x_min + rnd * (fillregion%x_max - fillregion%x_min)
      xi_ = INT(x_, 2); dx_ = x_ - REAL(xi_)
      rnd = random(dseed)
      y_ = fillregion%y_min + rnd * (fillregion%y_max - fillregion%y_min)
      yi_ = INT(y_, 2); dx_ = y_ - REAL(yi_)
      #ifdef threeD
        rnd = random(dseed)
        z_ = fillregion%z_min + rnd * (fillregion%z_max - fillregion%z_min)
        zi_ = INT(z_, 2); dz_ = z_ - REAL(zi_)
      #else
        zi_ = 0; dz_ = 0.5
      #endif

      ! if spatial distribution function is present, compute it
      !   otherwise use uniform distribution
      if (present(spat_distr_ptr)) then
        x_glob = REAL(this_meshblock%ptr%x0) + x_
        rnd = spat_distr_ptr(x_glob = x_glob, dummy1 = dummy1, dummy2 = dummy2)
      else
        rnd = 1.0
      end if
      if (random(dseed) .lt. rnd) then
        do s = 1, num_species
          ! generate momenta for every species individually
          spec_ = fill_species(s)
          !   shift direction is opposite for opposite signed species
          if (present(shift_gamma)) then
            fill_maxwellian%shift_dir = INT(SIGN(1.0, species(spec_)%ch_sp)) * shift_dir
          end if
          call generateFromMaxwellian(fill_maxwellian, u_, v_, w_)
          call createParticle(spec_, xi_, yi_, zi_, dx_, dy_, dz_, u_, v_, w_)
        end do
      end if
      n = n + 1
    end do
    call deallocateMaxwellian(fill_maxwellian)
  end subroutine fillRegionWithThermalPlasma
end module m_thermalplasma
