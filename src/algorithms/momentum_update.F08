#ifndef VAY
#include "boris_push.F08"
#else
#include "vay_push.F08"
#endif

! RADIATION >
#ifdef RADIATION
#ifdef SYNCHROTRON
                if (species(s) % cool_sp) then
                  call particleRadiateSync(timestep, s, &
                                           pt_u(p), pt_v(p), pt_w(p), u_init, v_init, w_init, &
                                           dx_rad, dy_rad, dz_rad, xi_rad, yi_rad, zi_rad, pt_wei(p), &
                                           bx_rad, by_rad, bz_rad, ex_rad, ey_rad, ez_rad, &
                                           index=pt_ind(p), proc=pt_proc(p))
                end if
#endif
#ifdef INVERSECOMPTON
                if (species(s) % cool_sp) then
                  call particleRadiateIC(timestep, s, &
                                         pt_u(p), pt_v(p), pt_w(p), u_init, v_init, w_init, &
                                         dx_rad, dy_rad, dz_rad, xi_rad, yi_rad, zi_rad, pt_wei(p), &
                                         bx_rad, by_rad, bz_rad, ex_rad, ey_rad, ez_rad, &
                                         index=pt_ind(p), proc=pt_proc(p))
                end if
#endif
#endif
! </ RADIATION
