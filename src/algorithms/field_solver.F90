#include "../defs.F90"

module m_fldsolver
  use m_globalnamespace
  use m_aux
  use m_communications
  use m_domain
  use m_fields
  implicit none
contains
  subroutine fillGhosts()
    implicit none
    integer :: i, j, k, imin, imax, jmin, jmax, kmin, kmax
    integer :: ind1, ind2, ind3
    integer :: send_cnt, status, ierr

    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifndef threeD
            if (ind3 .ne. 0) cycle
          #endif
          if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle

          ! highlight the region to send and save to `send_fld`
          if (ind1 .eq. 0) then
            imin = 0; imax = this_meshblock%ptr%sx - 1
          else if (ind1 .eq. -1) then
            imin = 0; imax = NGHOST - 1
          else if (ind1 .eq. 1) then
            imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx - 1
          end if
          if (ind2 .eq. 0) then
            jmin = 0; jmax = this_meshblock%ptr%sy - 1
          else if (ind2 .eq. -1) then
            jmin = 0; jmax = NGHOST - 1
          else if (ind2 .eq. 1) then
            jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy - 1
          end if
          if (ind3 .eq. 0) then
            kmin = 0; kmax = this_meshblock%ptr%sz - 1
          else if (ind3 .eq. -1) then
            kmin = 0; kmax = NGHOST - 1
          else if (ind3 .eq. 1) then
            kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz - 1
          end if
          #ifndef threeD
            kmin = 0; kmax = 0
          #endif

          ! write send/recv arrays in/from a given direction
          !     in 3D: 26 directions, in 2D: 8, in 1D: 2
          send_cnt = 0;
          do i = imin, imax
            do j = jmin, jmax
              do k = kmin, kmax
                send_fld(send_cnt + 1) = ex(i, j, k)
                send_cnt = send_cnt + 1
                send_fld(send_cnt + 1) = ey(i, j, k)
                send_cnt = send_cnt + 1
                send_fld(send_cnt + 1) = ez(i, j, k)
                send_cnt = send_cnt + 1
                send_fld(send_cnt + 1) = bx(i, j, k)
                send_cnt = send_cnt + 1
                send_fld(send_cnt + 1) = by(i, j, k)
                send_cnt = send_cnt + 1
                send_fld(send_cnt + 1) = bz(i, j, k)
                send_cnt = send_cnt + 1
              end do
            end do
          end do

          call MPI_SENDRECV(send_fld, send_cnt,&
                          & MPI_REAL, this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk, 200,&
                          & recv_fld, send_cnt,&
                          & MPI_REAL, this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk, 200,&
                          & MPI_COMM_WORLD, status, ierr)

          ! highlight the region to extract the `recv_fld`
          if (ind1 .eq. 0) then
            imin = 0; imax = this_meshblock%ptr%sx - 1
          else if (ind1 .eq. -1) then
            imin = -NGHOST; imax = -1
          else if (ind1 .eq. 1) then
            imin = this_meshblock%ptr%sx; imax = this_meshblock%ptr%sx + NGHOST - 1
          end if
          if (ind2 .eq. 0) then
            jmin = 0; jmax = this_meshblock%ptr%sy - 1
          else if (ind2 .eq. -1) then
            jmin = -NGHOST; jmax = -1
          else if (ind2 .eq. 1) then
            jmin = this_meshblock%ptr%sy; jmax = this_meshblock%ptr%sy + NGHOST - 1
          end if
          if (ind3 .eq. 0) then
            kmin = 0; kmax = this_meshblock%ptr%sz - 1
          else if (ind3 .eq. -1) then
            kmin = -NGHOST; kmax = -1
          else if (ind3 .eq. 1) then
            kmin = this_meshblock%ptr%sz; kmax = this_meshblock%ptr%sz + NGHOST - 1
          end if
          #ifndef threeD
            kmin = 0; kmax = 0
          #endif

          ! copy `recv_fld` to ghost cells
          send_cnt = 0
          do i = imin, imax
            do j = jmin, jmax
              do k = kmin, kmax
                ex(i, j, k) = recv_fld(send_cnt + 1)
                send_cnt = send_cnt + 1
                ey(i, j, k) = recv_fld(send_cnt + 1)
                send_cnt = send_cnt + 1
                ez(i, j, k) = recv_fld(send_cnt + 1)
                send_cnt = send_cnt + 1
                bx(i, j, k) = recv_fld(send_cnt + 1)
                send_cnt = send_cnt + 1
                by(i, j, k) = recv_fld(send_cnt + 1)
                send_cnt = send_cnt + 1
                bz(i, j, k) = recv_fld(send_cnt + 1)
                send_cnt = send_cnt + 1
              end do
            end do
          end do

        end do
      end do
    end do
  end subroutine fillGhosts
end module m_fldsolver
