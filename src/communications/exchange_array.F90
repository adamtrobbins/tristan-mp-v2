#include "../defs.F90"

module m_exchangearray
  use m_globalnamespace
  use m_aux
  use m_errors
  use m_domain
  use m_fields
contains

  subroutine findCnt(ind1, ind2, ind3, send_cnt)
    integer                       :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer, intent(in)           :: ind1, ind2, ind3    
    integer, intent(out)          :: send_cnt 

    ! highlight the region to send and save to `send_EB`
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
            
    send_cnt = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
  end subroutine findCnt

  subroutine bufferSendArray( ind1, ind2, ind3, send_cnt)
    integer                       :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer, intent(in)           :: ind1, ind2, ind3
    integer, intent(out)          :: send_cnt

    ! highlight the region to send and save to `send_fld`
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

    send_cnt = 1
    do i = imin, imax
       do j = jmin, jmax
          do k = kmin, kmax
             send_EB(send_cnt + 0) = lg_arr(i, j, k)
             send_cnt = send_cnt + 1
          end do
       end do
    end do
    send_cnt = send_cnt - 1
  end subroutine bufferSendArray

  subroutine extractRecvArray(ind1, ind2, ind3)
    implicit none
    integer                       :: imin, imax, jmin, jmax, kmin, kmax, i, j, k
    integer                       :: send_cnt
    integer, intent(in)           :: ind1, ind2, ind3


    ! write received data to local memory
    ! highlight the region to extract the `recv_fld`
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
    
    ! copy `recv_fld` to ghost cells
    send_cnt = 1
    do i = imin, imax
       do j = jmin, jmax
          do k = kmin, kmax
             jx_buff(i, j, k) = jx_buff(i, j, k) + recv_fld(send_cnt)
             send_cnt = send_cnt + 1
          end do
       end do
    end do
    
  end subroutine extractRecvArray

    
  subroutine exchangeArray()
    implicit none
    integer           :: i, j, k, imin, imax, jmin, jmax, kmin, kmax
    integer           :: ind1, ind2, ind3, cntr, n_cntr
    integer           :: cnt, ierr
    integer           :: mpi_sendto, mpi_recvfrom, mpi_tag
    integer           :: mpi_offset
    logical           :: should_send, should_recv

    #ifdef MPI08
      type(MPI_REQUEST), allocatable  :: mpi_req(:)
      type(MPI_STATUS)                :: istat
    #endif

    #ifdef MPI
      integer, allocatable            :: mpi_req(:)
      integer                         :: istat(MPI_STATUS_SIZE)
    #endif

    jx_buff(:,:,:) = 0.0

    do ind1 = -1, 1
      do ind2 = -1, 1
        do ind3 = -1, 1
          if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
          #ifdef oneD
            if ((ind2 .ne. 0) .or. (ind3 .ne. 0)) cycle
          #elif twoD
            if (ind3 .ne. 0) cycle
          #endif

          mpi_tag = (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1) + 300

          should_send = associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)
          should_recv = associated(this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr)
          if (should_send .and. should_recv) then

             mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
             mpi_recvfrom = this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk

             call bufferSendArray(ind1, ind2, ind3, cnt)
             call MPI_SENDRECV(send_EB(1 : cnt), cnt, MPI_REAL, mpi_sendto, mpi_tag,&
                  & recv_fld(1 : cnt), cnt, MPI_REAL, mpi_recvfrom, mpi_tag,&
                  & MPI_COMM_WORLD, istat, ierr)
             call extractRecvArray(-ind1, -ind2, -ind3)
          else if ((.not. should_send) .and. should_recv) then
            mpi_recvfrom = this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk 
            call findCnt(ind1, ind2, ind3, cnt)
            call MPI_RECV(recv_fld(1 : cnt),cnt, MPI_REAL,mpi_recvfrom, mpi_tag, MPI_COMM_WORLD, istat, ierr)
            call extractRecvArray(-ind1, -ind2, -ind3)
         else if ((.not. should_recv) .and. should_send) then
            mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk 
            call bufferSendArray(ind1, ind2, ind3, cnt)
            call MPI_SEND(send_EB(1 : cnt), cnt, MPI_REAL,mpi_sendto, mpi_tag, MPI_COMM_WORLD, istat, ierr)
         end if
        end do
      end do
    end do

    lg_arr(:,:,:) = lg_arr(:,:,:) + jx_buff(:,:,:)

    call printDiag("exchangeArray()", 3)
  end subroutine exchangeArray
end module m_exchangearray

