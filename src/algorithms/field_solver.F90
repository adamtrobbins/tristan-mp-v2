#include "../defs.F90"

module m_fldsolver
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_fields
  implicit none
contains
  subroutine advanceBHalfstep()
    implicit none
    integer :: i, j, k, ip1, jp1, kp1
    real :: const
    const = CORR * 0.5 * CC

    ! FIX0 needs to be changed for different BC-s

    #ifndef threeD
      k = 0
      do j = 0, this_meshblock%ptr%sy - 1
        jp1 = j + 1
        do i = 0, this_meshblock%ptr%sx - 1
          ip1 = i + 1
          bx(i, j, k) = bx(i, j, k) + const *&
                    & (-ez(i, jp1, k) + ez(i, j, k))
          by(i, j, k) = by(i, j, k) + const *&
                    & (ez(ip1, j, k) - ez(i, j, k))
          bz(i, j, k) = bz(i, j, k) + const *&
                    & (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
        enddo
      enddo
    #else
      do k = 0, this_meshblock%ptr%sz - 1
        kp1 = k + 1
        do j = 0, this_meshblock%ptr%sy - 1
          jp1 = j + 1
          do i = 0, this_meshblock%ptr%sx - 1
            ip1 = i + 1
            bx(i, j, k) = bx(i, j, k) + const *&
                      & (ey(i, j, kp1) - ey(i, j, k) - ez(i, jp1, k) + ez(i, j, k))
            by(i, j, k) = by(i, j, k) + const *&
                      & (ez(ip1, j, k) - ez(i, j, k) - ex(i, j, kp1) + ex(i, j, k))
            bz(i, j, k) = bz(i, j, k) + const *&
                      & (ex(i, jp1, k) - ex(i, j, k) - ey(ip1, j, k) + ey(i, j, k))
          enddo
        enddo
      enddo
    #endif
    call printDiag((mpi_rank .eq. 0), "advanceBHalfstep()", .true.)
  end subroutine advanceBHalfstep

  subroutine advanceEFullstep()
    implicit none
    integer :: i, j, k, im1, jm1, km1
    real :: const
    const = CORR * CC

    ! FIX0 needs to be changed for different BC-s

    #ifndef threeD
      k = 0
      do j = 0, this_meshblock%ptr%sy - 1
        jm1 = j - 1
        do i = 0, this_meshblock%ptr%sx - 1
          im1 = i - 1
          ex(i, j, k) = ex(i, j, k) + const *&
                    & (-bz(i, jm1, k) + bz(i, j, k))
          ey(i, j, k) = ey(i, j, k) + const *&
                    & (bz(im1, j, k) - bz(i, j, k))
          ez(i, j, k) = ez(i, j, k) + const *&
                    & (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
        enddo
      enddo
    #else
      do k = 0, this_meshblock%ptr%sz - 1
        km1 = k - 1
        do j = 0, this_meshblock%ptr%sy - 1
          jm1 = j - 1
          do i = 0, this_meshblock%ptr%sx - 1
            im1 = i - 1
            ex(i, j, k) = ex(i, j, k) + const *&
                      & (by(i, j, km1) - by(i, j, k) - bz(i, jm1, k) + bz(i, j, k))
            ey(i, j, k) = ey(i, j, k) + const *&
                      & (bz(im1, j, k) - bz(i, j, k) - bx(i, j, km1) + bx(i, j, k))
            ez(i, j, k) = ez(i, j, k) + const *&
                      & (bx(i, jm1, k) - bx(i, j, k) - by(im1, j, k) + by(i, j, k))
          enddo
        enddo
      enddo
    #endif
    call printDiag((mpi_rank .eq. 0), "advanceEFullstep()", .true.)
  end subroutine advanceEFullstep

  subroutine addCurrents()
    implicit none
    integer :: xmin, xmax, ymin, ymax, zmin, zmax
    xmin = 0;   xmax = this_meshblock%ptr%sx - 1
    ymin = 0;   ymax = this_meshblock%ptr%sy - 1
    zmin = 0;   zmax = this_meshblock%ptr%sz - 1
    ! "-" sign is taken care of in the deposit
    #ifndef threeD
      ex(xmin:xmax, ymin:ymax, 0) = &
          & ex(xmin:xmax, ymin:ymax, 0) + &
          & jx(xmin:xmax, ymin:ymax, 0)
      ey(xmin:xmax, ymin:ymax, 0) = &
          & ey(xmin:xmax, ymin:ymax, 0) + &
          & jy(xmin:xmax, ymin:ymax, 0)
      ez(xmin:xmax, ymin:ymax, 0) = &
          & ez(xmin:xmax, ymin:ymax, 0) + &
          & jz(xmin:xmax, ymin:ymax, 0)
    #else
      ex(xmin:xmax, ymin:ymax, zmin:zmax) = &
          & ex(xmin:xmax, ymin:ymax, zmin:zmax) + &
          & jx(xmin:xmax, ymin:ymax, zmin:zmax)
      ey(xmin:xmax, ymin:ymax, zmin:zmax) = &
          & ey(xmin:xmax, ymin:ymax, zmin:zmax) + &
          & jy(xmin:xmax, ymin:ymax, zmin:zmax)
      ez(xmin:xmax, ymin:ymax, zmin:zmax) = &
          & ez(xmin:xmax, ymin:ymax, zmin:zmax) + &
          & jz(xmin:xmax, ymin:ymax, zmin:zmax)
    #endif
    call printDiag((mpi_rank .eq. 0), "addCurrents()", .true.)
  end subroutine addCurrents
end module m_fldsolver
