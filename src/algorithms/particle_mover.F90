#include "../defs.F90"

module m_mover
  use m_globalnamespace
  use m_aux
  use m_communications
  use m_domain
  use m_particles
  implicit none
contains
  subroutine moveParticles()
    implicit none
    integer                               :: s, p, temp_i
    real                                  :: g_temp, e_temp, temp_r
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:),&
                                           & pt_u(:), pt_v(:), pt_w(:)

    do s = 1, nspec
      pt_xi => sp_(s)%xi; pt_yi => sp_(s)%yi; pt_zi => sp_(s)%zi
      pt_dx => sp_(s)%dx; pt_dy => sp_(s)%dy; pt_dz => sp_(s)%dz
      pt_u => sp_(s)%u; pt_v => sp_(s)%v; pt_w => sp_(s)%w
      if (spp_(s)%m_sp .eq. 0) then
        ! routine for massless particles
        !$omp simd
        !dir$ vector aligned
        do p = 1, spp_(s)%npart_sp
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
        !$omp simd
        !dir$ vector aligned
        do p = 1, spp_(s)%npart_sp
          g_temp = sqrt(1. + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)

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
    end do
    call printDiag((mpi_rank .eq. 0), TAB // "moveParticles()" // TAB // TAB // TAB // "[OK]")
  end subroutine moveParticles
end module m_mover
