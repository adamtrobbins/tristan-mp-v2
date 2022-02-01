#include "../defs.F90"

module m_exchangecurrents
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_fields
contains

  subroutine bufferSendArray( ind1, ind2, ind3,fill_ghosts, send_cnt)
    integer                       :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer, intent(in)           :: ind1, ind2, ind3
    logical, intent(in)           :: fill_ghosts
    integer, intent(out)          :: send_cnt

    ! highlight the region to send and save to `send_EB`
    if (.not. fill_ghosts) then
       !   sending ghost zones + normal zones
       if (ind1 .eq. 0) then
          imin = -NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
       else if (ind1 .eq. -1) then
          imin = -NGHOST; imax = NGHOST - 1
       else if (ind1 .eq. 1) then
          imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
       end if

       if (ind2 .eq. 0) then
          jmin = -NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
       else if (ind2 .eq. -1) then
          jmin = -NGHOST; jmax = NGHOST - 1
       else if (ind2 .eq. 1) then
          jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
       end if

       if (ind3 .eq. 0) then
          kmin = -NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
       else if (ind3 .eq. -1) then
          kmin = -NGHOST; kmax = NGHOST - 1
       else if (ind3 .eq. 1) then
          kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
       end if
#ifdef oneD
       jmin = 0; jmax = 0
       kmin = 0; kmax = 0
#elif twoD
       kmin = 0; kmax = 0
#endif
    else
       !   sending just the normal zones
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
#ifdef oneD
       jmin = 0; jmax = 0
       kmin = 0; kmax = 0
#elif twoD
       kmin = 0; kmax = 0
#endif
    end if

    send_cnt = 1
    do i = imin, imax
       do j = jmin, jmax
          do k = kmin, kmax
             send_EB(send_cnt + 0) = jx(i, j, k)
             send_EB(send_cnt + 1) = jy(i, j, k)
             send_EB(send_cnt + 2) = jz(i, j, k)
             send_cnt = send_cnt + 3
          end do
       end do
    end do
    send_cnt = send_cnt - 1
  end subroutine bufferSendArray

  subroutine extractRecvArray(ind1, ind2, ind3, fill_ghosts)
    implicit none
    integer                       :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer                       :: send_cnt
    integer, intent(in)           :: ind1, ind2, ind3
    logical, intent(in)           :: fill_ghosts    

    if (.not. fill_ghosts) then
       !   write to ghosts + normal zones
       if (ind1 .eq. 0) then
          imin = -NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
       else if (ind1 .eq. -1) then
          imin = -NGHOST; imax = NGHOST - 1
       else if (ind1 .eq. 1) then
          imin = this_meshblock%ptr%sx - NGHOST; imax = this_meshblock%ptr%sx + NGHOST - 1
       end if

       if (ind2 .eq. 0) then
          jmin = -NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
       else if (ind2 .eq. -1) then
          jmin = -NGHOST; jmax = NGHOST - 1
       else if (ind2 .eq. 1) then
          jmin = this_meshblock%ptr%sy - NGHOST; jmax = this_meshblock%ptr%sy + NGHOST - 1
       end if

       if (ind3 .eq. 0) then
          kmin = -NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
       else if (ind3 .eq. -1) then
          kmin = -NGHOST; kmax = NGHOST - 1
       else if (ind3 .eq. 1) then
          kmin = this_meshblock%ptr%sz - NGHOST; kmax = this_meshblock%ptr%sz + NGHOST - 1
       end if
#ifdef oneD
       jmin = 0; jmax = 0
       kmin = 0; kmax = 0
#elif twoD
       kmin = 0; kmax = 0
#endif
    else
       !   write to ghosts
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
#ifdef oneD
       jmin = 0; jmax = 0
       kmin = 0; kmax = 0
#elif twoD
       kmin = 0; kmax = 0
#endif
    end if
    

    send_cnt = 1
    do i = imin, imax
       do j = jmin, jmax
          do k = kmin, kmax
             if (.not. fill_ghosts) then
                ! add to existing values
                jx_buff(i, j, k) = jx_buff(i, j, k) + recv_fld(send_cnt + 0)
                jy_buff(i, j, k) = jy_buff(i, j, k) + recv_fld(send_cnt + 1)
                jz_buff(i, j, k) = jz_buff(i, j, k) + recv_fld(send_cnt + 2)
             else
                ! overwrite the existing values
                jx(i, j, k) = recv_fld(send_cnt + 0)
                jy(i, j, k) = recv_fld(send_cnt + 1)
                jz(i, j, k) = recv_fld(send_cnt + 2)
             end if
             send_cnt = send_cnt + 3
          end do
       end do
    end do
    
  end subroutine extractRecvArray

    
  subroutine exchangeCurrents(fill_ghosts_Q)
    implicit none
    integer           :: i, j, k, imin, imax, jmin, jmax, kmin, kmax
    integer           :: ind1, ind2, ind3, cntr, n_cntr
    integer           :: cnt, ierr
    integer           :: mpi_sendto, mpi_recvfrom, mpi_tag
    integer           :: mpi_offset
    logical           :: fill_ghosts
    logical           :: should_send, should_recv
    logical, optional, intent(in) :: fill_ghosts_Q

    #ifdef MPI08
      type(MPI_REQUEST), allocatable  :: mpi_req(:)
      type(MPI_STATUS)                :: istat
    #endif

    #ifdef MPI
      integer, allocatable            :: mpi_req(:)
      integer                         :: istat(MPI_STATUS_SIZE)
    #endif


    if (present(fill_ghosts_Q)) then
      fill_ghosts = fill_ghosts_Q
    else
      fill_ghosts = .false.
    end if

    if (.not. fill_ghosts) then
      jx_buff(:,:,:) = 0.0
      jy_buff(:,:,:) = 0.0
      jz_buff(:,:,:) = 0.0
    end if

    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
          #elif twoD
            if (ind3 .ne. 0) cycle
          #endif
          if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle

          mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
          mpi_recvfrom = this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk
          mpi_tag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1) + 50

          should_send = associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)
          should_recv = associated(this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr)
          if (should_send .and. should_recv) then

             call bufferSendArray(ind1, ind2, ind3, fill_ghosts, cnt)
             call MPI_SENDRECV(send_EB(1 : cnt), cnt, MPI_REAL, mpi_sendto, mpi_tag,&
                  & recv_fld(1 : cnt), cnt, MPI_REAL, mpi_recvfrom, mpi_tag,&
                  & MPI_COMM_WORLD, istat, ierr)
             call extractRecvArray(-ind1, -ind2, -ind3, fill_ghosts)
          else
             print *, "NOT SUPPOSED TO BE HERE"
             stop
          end if
        end do
      end do
    end do

    if (.not. fill_ghosts) then
      jx(:,:,:) = jx(:,:,:) + jx_buff(:,:,:)
      jy(:,:,:) = jy(:,:,:) + jy_buff(:,:,:)
      jz(:,:,:) = jz(:,:,:) + jz_buff(:,:,:)
    end if

    call printDiag("exchangeCurrents()", 2)
  end subroutine exchangeCurrents
end module m_exchangecurrents
