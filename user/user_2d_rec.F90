#include "../src/defs.F90"

module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  #ifdef USROUTPUT
    use m_writeusroutput
  #endif
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real, private    :: nCS_over_nUP, current_width, upstream_T, cs_x, cs_x1, cs_x2
  real, private    :: boost_y_Gamma, boost_y_beta
  real, private    :: injector_x1, injector_x2, injector_sx, injector_betax, measure_x
  integer, private :: injector_reset_interval, open_boundaries
  integer, private :: cs_lecs, cs_ions, up_lecs, up_ions
  logical, private :: perturb, simple_bc
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    call getInput('problem', 'upstream_T', upstream_T)
    call getInput('problem', 'current_width', current_width)
    call getInput('problem', 'nCS_nUP', nCS_over_nUP, 0.0)
    if (boundary_x .ne. 1) then
      call getInput('problem', 'injector_sx', injector_sx)
      call getInput('problem', 'injector_betax', injector_betax)
    end if
    call getInput('problem', 'cs_lecs', cs_lecs, 1)
    call getInput('problem', 'cs_ions', cs_ions, 2)
    call getInput('problem', 'up_lecs', up_lecs, 1)
    call getInput('problem', 'up_ions', up_ions, 2)
    call getInput('problem', 'measure_x', measure_x, 0.2)
    call getInput('problem', 'open_boundaries', open_boundaries, -1)
    call getInput('problem', 'perturb', perturb, .false.)
    call getInput('problem', 'simple_bc', simple_bc, .true.)
    if (current_width .lt. 0.0) then
      call throwError("ERROR: `current_width` has to be > 0.")
    end if
    if (boundary_x .eq. 1) then
      ! double periodic
      cs_x1 = 0.25; cs_x2 = 0.75
    else
      ! single current sheet
      cs_x = 0.5
    end if
  end subroutine userReadInput

  function userSLBload(x_glob, y_glob, z_glob,&
                     & dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    return
  end function

  function userSpatialDistribution(x_glob, y_glob, z_glob,&
                                 & dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional  :: x_glob, y_glob, z_glob
    real, intent(in), optional  :: dummy1, dummy2, dummy3
    real                        :: rad2
    if (present(x_glob) .and. present(y_glob) .and.&
      & present(dummy1) .and. present(dummy2)) then
      if (present(dummy3) .and. (dummy3 .ne. 0.0)) then
        rad2 = (x_glob - dummy1)**2 + (y_glob - dummy3)**2
        userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2 *&
                                & (1.0 - exp(-rad2 / (5.0 * dummy2)**2))
      else
        userSpatialDistribution = 1.0 / (cosh((x_glob - dummy1) / dummy2))**2
      end if
    else
      call throwError("ERROR: variable not present in `userSpatialDistribution()`")
    end if
    return
  end function userSpatialDistribution

  subroutine userInitParticles()
    implicit none
    real                :: nUP, nCS
    type(region)        :: back_region
    real                :: sx_glob, sy_glob, shift_gamma, shift_beta, current_sheet_T
    integer             :: s, ti, tj, tk, p
    real                :: ux, uy, uz, gamma
    procedure (spatialDistribution), pointer :: spat_distr_ptr => null()
    spat_distr_ptr => userSpatialDistribution

    nUP = 0.5 * ppc0
    nCS = nUP * nCS_over_nUP

    sx_glob = REAL(global_mesh%sx)
    sy_glob = REAL(global_mesh%sy)

    back_region%x_min = 0.0
    back_region%y_min = 0.0
    back_region%x_max = sx_glob
    back_region%y_max = sy_glob
    call fillRegionWithThermalPlasma(back_region, (/up_lecs, up_ions/), 2, nUP, upstream_T)

    if (nCS_over_nUP .ne. 0) then
      shift_beta = sqrt(sigma) * c_omp / (current_width * nCS_over_nUP)
      if (shift_beta .ge. 1) then
        call throwError('ERROR: `shift_beta` >= 1 in `userInitParticles()`')
      end if
      shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)
      current_sheet_T = 0.5 * sigma / nCS_over_nUP

      if (boundary_x .eq. 1) then
        back_region%x_min = sx_glob * cs_x1 - 10 * current_width
        back_region%x_max = sx_glob * cs_x1 + 10 * current_width
        back_region%y_min = 0.0
        back_region%y_max = sy_glob
        call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS, current_sheet_T,&
                                       & shift_gamma = shift_gamma, shift_dir = 3,&
                                       & spat_distr_ptr = spat_distr_ptr,&
                                       & dummy1 = cs_x1 * sx_glob, dummy2 = current_width)
        back_region%x_min = sx_glob * cs_x2 - 10 * current_width
        back_region%x_max = sx_glob * cs_x2 + 10 * current_width
        call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS, current_sheet_T,&
                                       & shift_gamma = shift_gamma, shift_dir = -3,&
                                       & spat_distr_ptr = spat_distr_ptr,&
                                       & dummy1 = cs_x2 * sx_glob, dummy2 = current_width)
      else
        back_region%x_min = sx_glob * cs_x - 10 * current_width
        back_region%x_max = sx_glob * cs_x + 10 * current_width
        back_region%y_min = 0.0
        back_region%y_max = sy_glob
        if (perturb) then
          call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS, current_sheet_T,&
                                         & shift_gamma = shift_gamma, shift_dir = 3,&
                                         & spat_distr_ptr = spat_distr_ptr,&
                                         & dummy1 = cs_x * sx_glob, dummy2 = current_width, dummy3 = cs_x * sy_glob)
        else
          call fillRegionWithThermalPlasma(back_region, (/cs_lecs, cs_ions/), 2, nCS, current_sheet_T,&
                                         & shift_gamma = shift_gamma, shift_dir = 3,&
                                         & spat_distr_ptr = spat_distr_ptr,&
                                         & dummy1 = cs_x * sx_glob, dummy2 = current_width)
        end if
      end if
    end if
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob
    real    :: x_glob, sx_glob
    ex(:,:,:) = 0; ey(:,:,:) = 0; ez(:,:,:) = 0
    bx(:,:,:) = 0; by(:,:,:) = 0; bz(:,:,:) = 0
    jx(:,:,:) = 0; jy(:,:,:) = 0; jz(:,:,:) = 0

    if (boundary_x .ne. 1) then
      injector_x1 = injector_sx - 1.0e-5
      injector_x2 = REAL(global_mesh%sx) - injector_sx + 1.0e-5
      injector_reset_interval = INT(injector_sx / (injector_betax * CC))
    end if

    k = 0
    sx_glob = REAL(global_mesh%sx)
    do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
      i_glob = i + this_meshblock%ptr%x0
      x_glob = REAL(i_glob) + 0.5
      if (boundary_x .eq. 1) then
        by(i,:,:) = tanh((x_glob - cs_x1 * sx_glob) / current_width) -&
                  & tanh((x_glob - cs_x2 * sx_glob) / current_width) - 1.0
      else
        by(i,:,:) = tanh((x_glob - cs_x * sx_glob) / current_width)
      end if
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
    ! called after particles move and deposit ...
    ! ... and before the currents are added to the electric field
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
    ! ... dummy loop ...
    ! integer :: s, ti, tj, tk, p
    ! do s = 1, nspec
    !   do ti = 1, species(s)%tile_nx
    !     do tj = 1, species(s)%tile_ny
    !       do tk = 1, species(s)%tile_nz
    !         do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
    !           ...
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
  end subroutine userDriveParticles

  subroutine userExternalFields(xp, yp, zp,&
                              & ex_ext, ey_ext, ez_ext,&
                              & bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in)  :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ! some functions of xp, yp, zp
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    real                            :: nUP, old_x1, old_x2, x_glob
    real                            :: ux, uy, uz, gamma
    integer                         :: s, ti, tj, tk, p, nUP_tot
    integer                         :: injector_i1_glob, injector_i2_glob
    type(region)                    :: back_region
    integer, optional, intent(in)             :: step
    procedure (spatialDistribution), pointer  :: spat_distr_ptr => null()

    if (boundary_x .ne. 1) then

      ! reset the injector position every once in a while
      if ((modulo(step, injector_reset_interval) .eq. 0) .and. (step .gt. 0)) then
        injector_x1 = injector_x1 +&
                          & REAL(injector_reset_interval) * CC * injector_betax
        injector_x2 = injector_x2 -&
                          & REAL(injector_reset_interval) * CC * injector_betax
      end if

      ! move the injectors
      old_x1 = injector_x1; old_x2 = injector_x2
      injector_x1 = injector_x1 - injector_betax * CC
      injector_x2 = injector_x2 + injector_betax * CC

      injector_i1_glob = INT(injector_x1)
      injector_i2_glob = INT(injector_x2)

      if (modulo(step, injector_reset_interval) .eq. 0) then
        ! remove particles left and right from the injectors every once in a while
        do s = 1, nspec
          do ti = 1, species(s)%tile_nx
            do tj = 1, species(s)%tile_ny
              do tk = 1, species(s)%tile_nz
                do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
                  x_glob = REAL(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0)&
                         & + species(s)%prtl_tile(ti, tj, tk)%dx(p)
                  if ((x_glob .le. old_x1) .or. (x_glob .gt. old_x2)) then
                    species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
                  end if
                end do
              end do
            end do
          end do
        end do
      end if

      ! inject background particles at the injectors' positions
      nUP = 0.5 * ppc0

      ! left injector
      back_region%x_min = injector_x1
      back_region%x_max = old_x1
      back_region%y_min = 0.0
      back_region%y_max = REAL(global_mesh%sy)

      call fillRegionWithThermalPlasma(back_region, (/up_lecs, up_ions/), 2, nUP, upstream_T)

      ! right injector
      back_region%x_min = old_x2
      back_region%x_max = injector_x2
      back_region%y_min = 0.0
      back_region%y_max = REAL(global_mesh%sy)

      call fillRegionWithThermalPlasma(back_region, (/up_lecs, up_ions/), 2, nUP, upstream_T)
    end if
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    real                          :: sx_glob, x_glob, delta_x
    integer                       :: i, j, k
    integer                       :: i_glob, injector_i1_glob, injector_i2_glob
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
    logical                       :: updateE_, updateB_

    if (present(updateE)) then
      updateE_ = updateE
    else
      updateE_ = .true.
    end if

    if (present(updateB)) then
      updateB_ = updateB
    else
      updateB_ = .true.
    end if

    if ((step .ge. open_boundaries) .and. (open_boundaries .ge. 0)) then
      boundary_y = 0
    end if

    if (boundary_x .ne. 1) then

      injector_i1_glob = INT(injector_x1) + 1
      injector_i2_glob = INT(injector_x2)

      if ((injector_i1_glob .lt. this_meshblock%ptr%x0 + this_meshblock%ptr%sx) .or.&
        & (injector_i2_glob .ge. this_meshblock%ptr%x0)) then
        ! reset fields left and right from the injectors
        if (updateB_) then
          sx_glob = REAL(global_mesh%sx)
          do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
            i_glob = i + this_meshblock%ptr%x0
            x_glob = REAL(i_glob) + 0.5
            if (i_glob .le. injector_i1_glob) then
              if (simple_bc) then
                bx(i,:,:) = 0.0
                bz(i,:,:) = 0.0
                by(i,:,:) = tanh((x_glob - cs_x * sx_glob) / current_width)
              else
                delta_x = 4.0 * REAL(i_glob) / MAX(REAL(injector_i1_glob), 0.1)
                bx(i, :, :) = tanh(delta_x) * bx(injector_i1_glob - this_meshblock%ptr%x0, :, :)
                bz(i, :, :) = tanh(delta_x) * bz(injector_i1_glob - this_meshblock%ptr%x0, :, :)
                by(i, :, :) = (1.0 - tanh(delta_x)) * tanh((x_glob - cs_x * sx_glob) / current_width) +&
                        & tanh(delta_x) * by(injector_i1_glob - this_meshblock%ptr%x0, :, :)
              end if
            else if (i_glob .ge. injector_i2_glob) then
              if (simple_bc) then
                bx(i,:,:) = 0.0
                bz(i,:,:) = 0.0
                by(i,:,:) = tanh((x_glob - cs_x * sx_glob) / current_width)
              else
              delta_x = 4.0 * REAL(global_mesh%sx - 1 - i_glob) / MAX(REAL(global_mesh%sx - 1 - injector_i2_glob), 0.1)
              bx(i, :, :) = tanh(delta_x) * bx(injector_i2_glob - this_meshblock%ptr%x0, :, :)
              bz(i, :, :) = tanh(delta_x) * bz(injector_i2_glob - this_meshblock%ptr%x0, :, :)
              by(i, :, :) = (1.0 - tanh(delta_x)) * tanh((x_glob - cs_x * sx_glob) / current_width) +&
                      & tanh(delta_x) * by(injector_i2_glob - this_meshblock%ptr%x0, :, :)
              end if
            end if
          end do
        end if
        if (updateE_) then
          sx_glob = REAL(global_mesh%sx)
          do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
            i_glob = i + this_meshblock%ptr%x0
            if (i_glob .lt. injector_i1_glob) then
              if (simple_bc) then
                ex(i,:,:) = 0.0
                ey(i,:,:) = 0.0
                ez(i,:,:) = 0.0
              else
                delta_x = 4.0 * REAL(i_glob) / MAX(REAL(injector_i1_glob), 0.1)
                ex(i, :, :) = tanh(delta_x) * ex(injector_i1_glob - this_meshblock%ptr%x0, :, :)
                ey(i, :, :) = tanh(delta_x) * ey(injector_i1_glob - this_meshblock%ptr%x0, :, :)
                ez(i, :, :) = tanh(delta_x) * ez(injector_i1_glob - this_meshblock%ptr%x0, :, :)
              end if
            else if (i_glob .gt. injector_i2_glob) then
              if (simple_bc) then
                ex(i,:,:) = 0.0
                ey(i,:,:) = 0.0
                ez(i,:,:) = 0.0
              else
                delta_x = 4.0 * REAL(global_mesh%sx - 1 - i_glob) / MAX(REAL(global_mesh%sx - 1 - injector_i2_glob), 0.1)
                ex(i, :, :) = tanh(delta_x) * ex(injector_i2_glob - this_meshblock%ptr%x0, :, :)
                ey(i, :, :) = tanh(delta_x) * ey(injector_i2_glob - this_meshblock%ptr%x0, :, :)
                ez(i, :, :) = tanh(delta_x) * ez(injector_i2_glob - this_meshblock%ptr%x0, :, :)
              end if
            end if
          end do
        end if
      end if
    end if
  end subroutine userFieldBoundaryConditions
  !............................................................!

  !--- user-specific output -----------------------------------!
  #ifdef USROUTPUT
    subroutine userOutput(step)
      implicit none
      integer, optional, intent(in) :: step
      integer                       :: root_rank = 0
      real, allocatable             :: y_bins(:), ExB_arr(:), ExB_arr_global(:)
      ! real                          :: dr, x_glob, y_glob, z_glob, r_glob
      real                          :: dummy_x, dummy_y, dummy_z, dummy
      ! real, allocatable             :: sum_ExBr_f(:), sum_f(:), sum_ExBr_f_global(:), sum_f_global(:)
      ! integer                       :: ri, rnum = 50, i, j, k, ierr
      integer                       :: x_bin, yi, ynum, i, j, k, ierr

      ynum = INT(global_mesh%sy)

      ! allocate(y_bins(ynum))
      allocate(ExB_arr(ynum))
      allocate(ExB_arr_global(ynum))
      ExB_arr(:) = 0.0

      ! do yi = 0, ynum - 1
      !   y_bins(yi + 1) = 0.5 + REAL(yi)
      ! end do

      if (this_meshblock%ptr%x0 .eq. 0) then
        x_bin = INT(measure_x * global_mesh%sx)
        i = x_bin; k = 0
        do j = 0, this_meshblock%ptr%sy - 1
          dummy_x = -(ez(i,j,k) * by(i,j,k)) + ey(i,j,k) * bz(i,j,k)
          dummy = bx(i,j,k)**2 + by(i,j,k)**2 + bz(i,j,k)**2

          yi = j + this_meshblock%ptr%y0
          ExB_arr(yi + 1) = dummy_x / dummy
        end do
      end if

      call MPI_REDUCE(ExB_arr, ExB_arr_global, ynum, MPI_REAL, MPI_SUM, root_rank, MPI_COMM_WORLD, ierr)

      if (mpi_rank .eq. root_rank) then
        call writeUsrOutputTimestep(step)
        ! call writeUsrOutputArray('y', y_bins)
        call writeUsrOutputArray('ExB', ExB_arr_global)
        call writeUsrOutputEnd()
      end if
    end subroutine userOutput

    logical function userExcludeParticles(s, ti, tj, tk, p)
      implicit none
      integer, intent(in)       :: s, ti, tj, tk, p
      userExcludeParticles = .true.
    end function userExcludeParticles
  #endif
  !............................................................!
end module m_userfile
