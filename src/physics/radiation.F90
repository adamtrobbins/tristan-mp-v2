#include "../defs.F90"

module m_radiation
  use m_globalnamespace
  use m_aux
  use m_errors
  implicit none

  real        :: rad_gamma_c, rad_gamma_rad, rad_beta_rec
  real        :: rad_beta_c, rad_beta_rad,&
              & rad_over_beta_c, rad_over_beta_rad
  real, allocatable :: rad_spectra(:), glob_rad_spectra(:)

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine particleRadiate(du, dv, dw,&
                           & u0, v0, w0,&
                           & ui, vi, wi,&
                           & bx, by, bz,&
                           & ex, ey, ez)
    implicit none
    real, intent(out) :: du, dv, dw
    real, intent(in)  :: u0, v0, w0, ui, vi, wi
    real, intent(in)  :: bx, by, bz, ex, ey, ez

    real :: uci, vci, wci, gci, over_gci

    real :: e_bar_x, e_bar_y, e_bar_z, e_bar_sq, beta_dot_e
    real :: chiR, chiR_sq, kappaR_x, kappaR_y, kappaR_z
    real :: tau_rad, eph_rad, dummy_

    integer :: spec_index

    uci = 0.5 * (u0 + ui)
    vci = 0.5 * (v0 + vi)
    wci = 0.5 * (w0 + wi)
    
    gci = sqrt(1.0 + uci**2 + vci**2 + wci**2)
    over_gci = 1.0 / gci 
    
    e_bar_x = ex + (vci * bz - wci * by) * over_gci
    e_bar_y = ey + (wci * bx - wci * bz) * over_gci
    e_bar_z = ez + (wci * by - wci * bx) * over_gci
    e_bar_sq = e_bar_x**2 + e_bar_y**2 + e_bar_z**2
    beta_dot_e = (ex * uci + ey * vci + ez * wci) * over_gci

    chiR_sq = e_bar_sq - beta_dot_E**2
    chiR = sqrt(chiR_sq)

    kappaR_x = (bz * e_bar_y - by * e_bar_z) + (ex * beta_dot_e)
    kappaR_y = (-bz * e_bar_x + bx * e_bar_z) + (ey * beta_dot_e)
    kappaR_z = (by * e_bar_x - bx * e_bar_y) + (ez * beta_dot_e)
    
    tau_rad = (rad_beta_rec * rad_beta_c * rad_over_beta_rad) * (rad_gamma_c / rad_gamma_rad)**2 * (chiR * B_norm * CCINV)
    eph_rad = (gci / rad_gamma_c)**2 * chiR * rad_over_beta_c

    dummy_ = rad_beta_rec * rad_over_beta_rad / (rad_gamma_rad**2 * CC)

    du = dummy_ * (kappaR_x - chiR_sq * gci * uci)
    dv = dummy_ * (kappaR_y - chiR_sq * gci * vci)
    dw = dummy_ * (kappaR_z - chiR_sq * gci * wci)
    
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
