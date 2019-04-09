#include "../defs.F90"

module m_exchangeparts
  use m_globalnamespace
  use m_aux
  use m_communications
  use m_domain
  use m_particles
  !--- PRIVATE functions -----------------------------------------!
  private :: copyToEnroute, copyFromEnroute,&
           & extractParticlesFromEnroute
  !...............................................................!
contains
  subroutine exchangeParticles()
    implicit none
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    integer, pointer, contiguous          :: pt_proc(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:)
    integer             :: s, p, send_x, send_y, send_z
    integer             :: mpi_sendto, mpi_recvfrom, mpi_sendtag, mpi_recvtag
    integer             :: ierr, ind1, ind2, ind3, cntr, temp_xyz
    type(MPI_STATUS)    :: istat
    integer             :: cnt_recv_enroute
    integer(kind=2)     :: new_xyz

    type(MPI_REQUEST), allocatable    :: mpi_req(:)
    type(MPI_STATUS), allocatable     :: mpi_stat(:)
    logical, allocatable              :: mpi_sendflags(:), mpi_recvflags(:)
    logical                           :: quit_loop

    allocate(mpi_req(sendrecv_neighbors))
    allocate(mpi_sendflags(sendrecv_neighbors))
    allocate(mpi_recvflags(sendrecv_neighbors))

    do s = 1, nspec ! loop over species
      enroute_bot%get(:,:,:)%cnt_send = 0

      pt_xi => sp_(s)%xi; pt_yi => sp_(s)%yi; pt_zi => sp_(s)%zi
      pt_dx => sp_(s)%dx; pt_dy => sp_(s)%dy; pt_dz => sp_(s)%dz
      pt_proc => sp_(s)%proc
      ! FIX1 make sure this is vectorized
      !$omp simd
      !dir$ vector aligned
      do p = 1, spp_(s)%npart_sp
        ! send_* = -1 / 0 / +1
        send_z = 0
        send_x = (ISIGN(1, pt_xi(p) - this_meshblock%ptr%sx) + 1) / 2 - (ISIGN(1, -pt_xi(p) - 1) + 1) / 2
        send_y = (ISIGN(1, pt_yi(p) - this_meshblock%ptr%sy) + 1) / 2 - (ISIGN(1, -pt_yi(p) - 1) + 1) / 2
        #ifdef threeD
          send_z = (ISIGN(1, pt_zi(p) - this_meshblock%ptr%sz) + 1) / 2 - (ISIGN(1, -pt_zi(p) - 1) + 1) / 2
        #endif
        ! FIX1 check for null() boundaries
        if ((send_x .ne. 0) .or. (send_y .ne. 0) .or. (send_z .ne. 0)) then
          ! this particle should be sent somewhere
          if (.not. associated(this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr)) then
            ! make ghost particle
            pt_proc(p) = -pt_proc(p) - 1
            cycle
          end if
          ! copy this particle to temporary `enroute_bot` array
          enroute_bot%get(send_x, send_y, send_z)%cnt_send = enroute_bot%get(send_x, send_y, send_z)%cnt_send + 1
          cntr = enroute_bot%get(send_x, send_y, send_z)%cnt_send
          call copyToEnroute(s, p, enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr))
          ! make ghost particle
          pt_proc(p) = -pt_proc(p) - 1

          ! shift coordinates to fit the new grid
          new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%xi
          temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx
          new_xyz = -(send_x - 1) * (2 + send_x) * (new_xyz * (send_x + 1) - (temp_xyz - 1) * send_x) / 2
          enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%xi = new_xyz

          new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%yi
          temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy
          new_xyz = -(send_y - 1) * (2 + send_y) * (new_xyz * (send_y + 1) - (temp_xyz - 1) * send_y) / 2
          enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%yi = new_xyz

          #ifdef threeD
            new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%zi
            temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz
            new_xyz = -(send_z - 1) * (2 + send_z) * (new_xyz * (send_z + 1) - (temp_xyz - 1) * send_z) / 2
            enroute_bot%get(send_x, send_y, send_z)%send_enroute(cntr)%zi = new_xyz
          #endif
        end if
      end do
      pt_xi => null(); pt_yi => null(); pt_zi => null()
      pt_dx => null(); pt_dy => null(); pt_dz => null()
      pt_proc => null()

      cntr = 0
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
            #ifndef threeD
              if (ind3 .ne. 0) cycle
            #endif
            if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
            cntr = cntr + 1
            mpi_sendtag = 100 * (mpi_rank + 1) + (ind3 + 2) + 3 * (ind2 + 1) + 9 * (ind1 + 1)

            ! post non-blocking send requests
            mpi_sendto = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
            call MPI_ISEND(enroute_bot%get(ind1,ind2,ind3)%send_enroute(1:enroute_bot%get(ind1,ind2,ind3)%cnt_send),&
                         & enroute_bot%get(ind1,ind2,ind3)%cnt_send, myMPI_ENROUTE,&
                         & mpi_sendto, mpi_sendtag, MPI_COMM_WORLD, mpi_req(cntr), ierr)
          end do
        end do
      end do

      ! wait to send & receive all the MPI calls and write data to memory
      quit_loop = .false.
      mpi_sendflags(:) = .false.
      mpi_recvflags(:) = .false.
      do while (.not. quit_loop)
        ! try to receive while not done
        quit_loop = .true.
        cntr = 0
        do ind1 = -1, 1
          do ind2 = -1, 1
            do ind3 = -1, 1
              if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
              #ifndef threeD
                if (ind3 .ne. 0) cycle
              #endif
              if (.not. associated(this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr)) cycle
              cntr = cntr + 1

              if (.not. mpi_sendflags(cntr)) then
                ! check if the message has been sent
                quit_loop = .false.
                call MPI_TEST(mpi_req(cntr), mpi_sendflags(cntr), istat, ierr)
              end if

              mpi_recvfrom = this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk
              mpi_recvtag = 100 * (mpi_recvfrom + 1) + (-ind3 + 2) + 3 * (-ind2 + 1) + 9 * (-ind1 + 1)

              if (.not. mpi_recvflags(cntr)) then
                quit_loop = .false.
                call MPI_IPROBE(mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, mpi_recvflags(cntr), istat)
                if (mpi_recvflags(cntr)) then
                  call MPI_GET_COUNT(istat, myMPI_ENROUTE, cnt_recv_enroute, ierr)
                  call MPI_RECV(recv_enroute(1:cnt_recv_enroute), cnt_recv_enroute, myMPI_ENROUTE,&
                              & mpi_recvfrom, mpi_recvtag, MPI_COMM_WORLD, istat, ierr)

                  ! write received data to local memory
                  if (cnt_recv_enroute .gt. 0) then
                    call extractParticlesFromEnroute(cnt_recv_enroute, s)
                  end if ! if > 0 particles received
                end if ! if the message can be received -> get the size & receive it
              end if ! if the message has already been received

            end do ! ind3
          end do ! ind2
        end do ! ind1
      end do ! global loop

    end do ! loop over species
    call printDiag((mpi_rank .eq. 0), TAB // "exchangeParticles()" // TAB // TAB // "[OK]")
  end subroutine exchangeParticles

  subroutine clearGhostParticles()
    implicit none
    integer :: s, p
    do s = 1, nspec
      ! FIX1 make sure this is vectorized (function call etc)
      p = 1
      do while (p .le. spp_(s)%npart_sp)
        if (sp_(s)%proc(p) .lt. 0) then
          call removeParticle(s, p)
        else
          p = p + 1
        end if
      end do
    end do
    call printDiag((mpi_rank .eq. 0), TAB // "clearGhostParticles()" // TAB // TAB // "[OK]")
  end subroutine clearGhostParticles

  subroutine copyToEnroute(spec_id, prtl_id, enroute)
    implicit none
    integer, intent(in)               :: spec_id, prtl_id
    type(prtl_enroute), intent(inout) :: enroute
    enroute%xi = sp_(spec_id)%xi(prtl_id)
    enroute%yi = sp_(spec_id)%yi(prtl_id)
    enroute%zi = sp_(spec_id)%zi(prtl_id)
    enroute%dx = sp_(spec_id)%dx(prtl_id)
    enroute%dy = sp_(spec_id)%dy(prtl_id)
    enroute%dz = sp_(spec_id)%dz(prtl_id)
    enroute%u = sp_(spec_id)%u(prtl_id)
    enroute%v = sp_(spec_id)%v(prtl_id)
    enroute%w = sp_(spec_id)%w(prtl_id)
    enroute%ind = sp_(spec_id)%ind(prtl_id)
    enroute%proc = sp_(spec_id)%proc(prtl_id)
  end subroutine copyToEnroute

  subroutine copyFromEnroute(enroute, spec_id, prtl_id)
    implicit none
    type(prtl_enroute), intent(inout) :: enroute
    integer, intent(in)               :: spec_id, prtl_id
    sp_(spec_id)%xi(prtl_id) = enroute%xi
    sp_(spec_id)%yi(prtl_id) = enroute%yi
    sp_(spec_id)%zi(prtl_id) = enroute%zi
    sp_(spec_id)%dx(prtl_id) = enroute%dx
    sp_(spec_id)%dy(prtl_id) = enroute%dy
    sp_(spec_id)%dz(prtl_id) = enroute%dz
    sp_(spec_id)%u(prtl_id) = enroute%u
    sp_(spec_id)%v(prtl_id) = enroute%v
    sp_(spec_id)%w(prtl_id) = enroute%w
    sp_(spec_id)%ind(prtl_id) = enroute%ind
    #ifdef DEBUG
      sp_(spec_id)%proc(prtl_id) = mpi_rank
    #else
      sp_(spec_id)%proc(prtl_id) = enroute%proc
    #endif
  end subroutine copyFromEnroute

  subroutine extractParticlesFromEnroute(cnt, spec_id)
    implicit none
    integer, intent(in)            :: cnt, spec_id
    integer                        :: p
    do p = 1, cnt
      spp_(spec_id)%npart_sp = spp_(spec_id)%npart_sp + 1
      call copyFromEnroute(recv_enroute(p), spec_id, spp_(spec_id)%npart_sp)
    end do
  end subroutine extractParticlesFromEnroute

end module m_exchangeparts
