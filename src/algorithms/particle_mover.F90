#include "../defs.F90"

module m_mover
  use m_globalnamespace
  use m_aux
  use m_helpers
  use m_domain
  use m_particles
  implicit none
contains
  subroutine moveParticles()
    implicit none
    integer                               :: s, p, temp_i, ti, tj, tk
    real                                  :: g_temp, e_temp, temp_r
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:),&
                                           & pt_u(:), pt_v(:), pt_w(:)
    real                                  :: ex0, ey0, ez0, bx0, by0, bz0, q_over_m
    real                                  :: u0, v0, w0, u1, v1, w1, dummy_

    do s = 1, nspec
			do ti = 1, species(s)%tile_nx
				do tj = 1, species(s)%tile_ny
					do tk = 1, species(s)%tile_nz
			      pt_xi => species(s)%prtl_tile(ti, tj, tk)%xi
						pt_yi => species(s)%prtl_tile(ti, tj, tk)%yi
						pt_zi => species(s)%prtl_tile(ti, tj, tk)%zi

						pt_dx => species(s)%prtl_tile(ti, tj, tk)%dx
						pt_dy => species(s)%prtl_tile(ti, tj, tk)%dy
						pt_dz => species(s)%prtl_tile(ti, tj, tk)%dz

						pt_u => species(s)%prtl_tile(ti, tj, tk)%u
						pt_v => species(s)%prtl_tile(ti, tj, tk)%v
						pt_w => species(s)%prtl_tile(ti, tj, tk)%w

			      if (species(s)%m_sp .eq. 0) then
			        ! routine for massless particles
			        !$omp simd
			        !dir$ vector aligned
			        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
			          e_temp = sqrt(pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)

			          pt_dx(p) = pt_dx(p) + CC * pt_u(p) / e_temp
			          temp_i = INT(pt_dx(p))
			          temp_r = MAX(SIGN(1., pt_dx(p)) + temp_i, REAL(temp_i)) - 1
			          temp_i = INT(temp_r)
			          pt_xi(p) = pt_xi(p) + temp_i
			          pt_dx(p) = pt_dx(p) - temp_r

			          pt_dy(p) = pt_dy(p) + CC * pt_v(p) / e_temp
			          temp_i = INT(pt_dy(p))
			          temp_r = MAX(SIGN(1., pt_dy(p)) + temp_i, REAL(temp_i)) - 1
			          temp_i = INT(temp_r)
			          pt_yi(p) = pt_yi(p) + temp_i
			          pt_dy(p) = pt_dy(p) - temp_r

			          #ifdef threeD
			            pt_dz(p) = pt_dz(p) + CC * pt_w(p) / e_temp
			            temp_i = INT(pt_dz(p))
			            temp_r = MAX(SIGN(1., pt_dz(p)) + temp_i, REAL(temp_i)) - 1
			            temp_i = INT(temp_r)
			            pt_zi(p) = pt_zi(p) + temp_i
			            pt_dz(p) = pt_dz(p) - temp_r
			          #endif
			        end do
			      else
			        ! routine for massive particles
			        q_over_m = species(s)%ch_sp / species(s)%m_sp
			        !$omp simd
			        !dir$ vector aligned
			        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
			          ! advance the velocity >>
			          call interpFlds(pt_dx(p), pt_dy(p), pt_dz(p),&
			                        & pt_xi(p), pt_yi(p), pt_zi(p),&
			                        & ex0, ey0, ez0, bx0, by0, bz0)
			          dummy_ = 0.5 * CCINV * q_over_m * B_norm
			          ex0 = ex0 * dummy_; ey0 = ey0 * dummy_; ez0 = ez0 * dummy_
			          bx0 = bx0 * dummy_; by0 = by0 * dummy_; bz0 = bz0 * dummy_
			          ! BORIS PUSHER >
			          ! half acceleration:
			          u0 = CC * pt_u(p) + ex0
			          v0 = CC * pt_v(p) + ey0
			          w0 = CC * pt_w(p) + ez0
			          ! first half magnetic rotation:
			          g_temp = CC / sqrt(CC**2 + u0**2 + v0**2 + w0**2)
			          bx0 = g_temp * bx0
			      		by0 = g_temp * by0
			      		bz0 = g_temp * bz0
			          dummy_ = 2. / (1. + bx0 * bx0 + by0 * by0 + bz0 * bz0)
			      		u1 = (u0 + v0 * bz0 - w0 * by0) * dummy_
			      		v1 = (v0 + w0 * bx0 - u0 * bz0) * dummy_
			      		w1 = (w0 + u0 * by0 - v0 * bx0) * dummy_
			          ! second half magnetic rotation + half acceleration:
			          u0 = u0 + v1 * bz0 - w1 * by0 + ex0
			      		v0 = v0 + w1 * bx0 - u1 * bz0 + ey0
			      		w0 = w0 + u1 * by0 - v1 * bx0 + ez0
			          ! </ BORIS PUSHER
			          pt_u(p) = u0 * CCINV
			          pt_v(p) = v0 * CCINV
			          pt_w(p) = w0 * CCINV
			          ! <</ advance the velocity

			          g_temp = sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)

			          pt_dx(p) = pt_dx(p) + CC * pt_u(p) / g_temp
			          temp_i = INT(pt_dx(p))
			          temp_r = MAX(SIGN(1., pt_dx(p)) + temp_i, REAL(temp_i)) - 1
			          temp_i = INT(temp_r)
			          pt_xi(p) = pt_xi(p) + temp_i
			          pt_dx(p) = pt_dx(p) - temp_r

			          pt_dy(p) = pt_dy(p) + CC * pt_v(p) / g_temp
			          temp_i = INT(pt_dy(p))
			          temp_r = MAX(SIGN(1., pt_dy(p)) + temp_i, REAL(temp_i)) - 1
			          temp_i = INT(temp_r)
			          pt_yi(p) = pt_yi(p) + temp_i
			          pt_dy(p) = pt_dy(p) - temp_r

			          #ifdef threeD
			            pt_dz(p) = pt_dz(p) + CC * pt_w(p) / g_temp
			            temp_i = INT(pt_dz(p))
			            temp_r = MAX(SIGN(1., pt_dz(p)) + temp_i, REAL(temp_i)) - 1
			            temp_i = INT(temp_r)
			            pt_zi(p) = pt_zi(p) + temp_i
			            pt_dz(p) = pt_dz(p) - temp_r
			          #endif
			        end do
			      end if
			      pt_xi => null(); pt_yi => null(); pt_zi => null()
			      pt_dx => null(); pt_dy => null(); pt_dz => null()
			      pt_u => null(); pt_v => null(); pt_w => null()
					end do ! tk
				end do ! tj
			end do ! ti
    end do ! species
    call printDiag((mpi_rank .eq. 0), TAB // "moveParticles()" // TAB // TAB // TAB // "[OK]")
  end subroutine moveParticles
end module m_mover
