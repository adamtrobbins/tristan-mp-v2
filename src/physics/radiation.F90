#include "../defs.F90"

module m_radiation
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_particlelogistics
  implicit none

  real              :: rad_gamma_c, rad_gamma_rad, rad_beta_rec
  real, allocatable :: rad_spectra(:), glob_rad_spectra(:)
  integer           :: rad_photon_ind

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine particleRadiate(u0, v0, w0, ui, vi, wi,&
                           & dx, dy, dz, xi, yi, zi,&
                           & bx, by, bz, ex, ey, ez)
    implicit none
    real, intent(inout)           :: u0, v0, w0
    real, intent(in)              :: ui, vi, wi
    real, intent(in)              :: bx, by, bz, ex, ey, ez
    real, intent(in)              :: dx, dy, dz
    integer(kind=2), intent(in)   :: xi, yi, zi
    real :: corr_

    real :: uci, vci, wci, kx, ky, kz, gci, betaci, over_gci
    real :: g_new, beta_new

    real :: e_bar_x, e_bar_y, e_bar_z, e_bar_sq, beta_dot_e
    real :: chiR, chiR_sq, kappaR_x, kappaR_y, kappaR_z
    real :: tau_rad, eph_rad, dummy_

    integer :: spec_index

    uci = 0.5 * (u0 + ui)
    vci = 0.5 * (v0 + vi)
    wci = 0.5 * (w0 + wi)
    
    gci = sqrt(1.0 + uci**2 + vci**2 + wci**2)
    over_gci = 1.0 / gci 
    betaci = sqrt(1.0 - over_gci**2)
    
    e_bar_x = ex + (vci * bz - wci * by) * over_gci
    e_bar_y = ey + (wci * bx - uci * bz) * over_gci
    e_bar_z = ez + (uci * by - vci * bx) * over_gci
    e_bar_sq = e_bar_x**2 + e_bar_y**2 + e_bar_z**2
    beta_dot_e = (ex * uci + ey * vci + ez * wci) * over_gci

    chiR_sq = abs(e_bar_sq - beta_dot_E**2)
    chiR = sqrt(chiR_sq)

    kappaR_x = (bz * e_bar_y - by * e_bar_z) + (ex * beta_dot_e)
    kappaR_y = (-bz * e_bar_x + bx * e_bar_z) + (ey * beta_dot_e)
    kappaR_z = (by * e_bar_x - bx * e_bar_y) + (ez * beta_dot_e)
    
    tau_rad = (rad_beta_rec * betaci) * (rad_gamma_c / rad_gamma_rad)**2 * (chiR * B_norm * CCINV)
    eph_rad = (gci / rad_gamma_c)**2 * chiR

    dummy_ = B_norm * rad_beta_rec / (rad_gamma_rad**2 * CC)
    #ifndef EMIT
      u0 = u0 + dummy_ * (kappaR_x - chiR_sq * gci * uci)
      v0 = v0 + dummy_ * (kappaR_y - chiR_sq * gci * vci)
      w0 = w0 + dummy_ * (kappaR_z - chiR_sq * gci * wci)
    #else
      u0 = u0 + dummy_ * kappaR_x
      v0 = v0 + dummy_ * kappaR_y
      w0 = w0 + dummy_ * kappaR_z

      g_new = sqrt(1.0 + u0**2 + v0**2 + w0**2)
      g_new = g_new - tau_rad * eph_rad
      beta_new = sqrt(1.0 - 1.0 / g_new**2)
      over_gci = 1.0 / sqrt(u0**2 + v0**2 + w0**2)
      kx = u0 * over_gci; ky = v0 * over_gci; kz = w0 * over_gci

      u0 = kx * g_new * beta_new
      v0 = ky * g_new * beta_new
      w0 = kz * g_new * beta_new
      if (random(dseed) .lt. tau_rad) then
        !if (g_new .gt. 1.5) then
        !  if (eph_rad .gt. g_new - 1.0) then
        !   eph_rad = g_new - 1.5
        !  end if
        call createParticle(rad_photon_ind, xi, yi, zi, dx, dy, dz,&
                          & kx * eph_rad, ky * eph_rad, kz * eph_rad) 
      end if
    #endif

    eph_rad = log(eph_rad)
    if (eph_rad .le. spec_min) then
      spec_index = 1
    else if (eph_rad .ge. spec_max) then
      spec_index = spec_num
    else
      spec_index = INT(CEILING((eph_rad - spec_min) * REAL(spec_num) / (spec_max - spec_min)))
      if (spec_index .lt. 1) spec_index = 1
      if (spec_index .gt. spec_num) spec_index = spec_num
    end if
    rad_spectra(spec_index) = rad_spectra(spec_index) + tau_rad
  end subroutine particleRadiate
end module m_radiation
