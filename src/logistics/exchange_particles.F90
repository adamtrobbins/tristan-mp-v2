#include "../defs.F90"

module m_exchangeparts
  use m_globalnamespace
  use m_aux
  use m_communications
  use m_domain
  use m_particles
contains

  subroutine exchangeParticles()
    implicit none
    integer(kind=2), pointer, contiguous  :: pt_xi(:), pt_yi(:), pt_zi(:)
    integer, pointer, contiguous          :: pt_proc(:)
    real, pointer, contiguous             :: pt_dx(:), pt_dy(:), pt_dz(:)
    integer                               :: s, p, send_x, send_y, send_z
    integer                               :: status(mpi_statsize), ierr, ind1, ind2, ind3, temp_cntr, temp_xyz
    integer                               :: cnt_recv_enroute
    integer(kind=2)                       :: new_xyz

    do s = 1, nspec
      ! loop over species
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            enroute_bot%get(ind1,ind2,ind3)%cnt_send = 0
          end do
        end do
      end do

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
          enroute_bot%get(send_x, send_y, send_z)%cnt_send = enroute_bot%get(send_x, send_y, send_z)%cnt_send + 1
          temp_cntr = enroute_bot%get(send_x, send_y, send_z)%cnt_send
          call copyToEnroute(s, p, enroute_bot%get(send_x, send_y, send_z)%send_enroute(temp_cntr))

          ! change coordinates
          new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(temp_cntr)%xi
          temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sx - 1
          new_xyz = -(send_x - 1) * (2 + send_x) * (new_xyz * (send_x + 1) - (temp_xyz - 1) * send_x) / 2
          enroute_bot%get(send_x, send_y, send_z)%send_enroute(temp_cntr)%xi = new_xyz

          new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(temp_cntr)%yi
          temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sy - 1
          new_xyz = -(send_y - 1) * (2 + send_y) * (new_xyz * (send_y + 1) - (temp_xyz - 1) * send_y) / 2
          enroute_bot%get(send_x, send_y, send_z)%send_enroute(temp_cntr)%yi = new_xyz

          #ifdef threeD
            new_xyz = enroute_bot%get(send_x, send_y, send_z)%send_enroute(temp_cntr)%zi
            temp_xyz = this_meshblock%ptr%neighbor(send_x, send_y, send_z)%ptr%sz - 1
            new_xyz = -(send_z - 1) * (2 + send_z) * (new_xyz * (send_z + 1) - (temp_xyz - 1) * send_z) / 2
            enroute_bot%get(send_x, send_y, send_z)%send_enroute(temp_cntr)%zi = new_xyz
          #endif
          ! make ghost particle
          pt_proc(p) = -pt_proc(p) - 1
        end if
      end do
      pt_xi => null(); pt_yi => null(); pt_zi => null()
      pt_dx => null(); pt_dy => null(); pt_dz => null()
      pt_proc => null()

      ! FIX1 check for null() boundaries
      do ind1 = -1, 1
        do ind2 = -1, 1
          do ind3 = -1, 1
            if ((ind1 .eq. 0) .and. (ind2 .eq. 0) .and. (ind3 .eq. 0)) cycle
            #ifndef threeD
              if (ind3 .ne. 0) cycle
            #endif
            ! send the # of particles going to one direction & get # incoming from opposite direction
            call MPI_SENDRECV(enroute_bot%get(ind1,ind2,ind3)%cnt_send, 1, MPI_INTEGER,&
                            & this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk, 100,&
                            & cnt_recv_enroute, 1, MPI_INTEGER,&
                            & this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk, 100,&
                            & MPI_COMM_WORLD, status, ierr)
            ! send the particles going to one direction & get those incoming from the opposite direction
            call MPI_SENDRECV(enroute_bot%get(ind1,ind2,ind3)%send_enroute, enroute_bot%get(ind1,ind2,ind3)%cnt_send,&
                            & myMPI_ENROUTE, this_meshblock%ptr%neighbor(ind1,ind2,ind3)%ptr%rnk, 100,&
                            & recv_enroute, cnt_recv_enroute,&
                            & myMPI_ENROUTE, this_meshblock%ptr%neighbor(-ind1,-ind2,-ind3)%ptr%rnk, 100,&
                            & MPI_COMM_WORLD, status, ierr)
            call extractParticlesFromEnroute(cnt_recv_enroute, s)
          end do
        end do
      end do
    end do
    call print_diag((mpi_rank .eq. 0), TAB // "exchangeParticles()" // TAB // TAB // "[OK]")
  end subroutine exchangeParticles

  subroutine clearGhostParticles()
    implicit none
    integer, pointer, contiguous       :: pt_proc(:)
    integer                            :: s, p
    do s = 1, nspec
      pt_proc => sp_(s)%proc
      ! FIX1 make sure this is vectorized (function call)
      !$omp simd
      !dir$ vector aligned
      do p = 1, spp_(s)%npart_sp
        if (pt_proc(p) .lt. 0) call removeParticle(s, p)
      end do
      pt_proc => null()
    end do
    call print_diag((mpi_rank .eq. 0), TAB // "clearGhostParticles()" // TAB // TAB // "[OK]")
  end subroutine clearGhostParticles

end module m_exchangeparts
