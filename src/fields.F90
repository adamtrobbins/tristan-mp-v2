#include "defs.F90"

!--- FIELDS ----------------------------------------------------!
! To store all the field related quantities
!   - field indices go
!       from `-NGHOST`
!       to `sx/sy/sz - 1 + NGHOST` inclusively
!...............................................................!

module m_fields
  use m_globalnamespace
  use m_domain
  implicit none

  real, allocatable     :: ex(:,:,:), ey(:,:,:), ez(:,:,:),&
                         & bx(:,:,:), by(:,:,:), bz(:,:,:)
  real, allocatable     :: jx(:,:,:), jy(:,:,:), jz(:,:,:)
  real, allocatable     :: jx_buff(:,:,:), jy_buff(:,:,:), jz_buff(:,:,:)
  real, allocatable     :: recv_fld(:), send_fld(:)
  integer               :: sendrecv_buffsz, sendrecv_offsetsz
  real, allocatable     :: lg_arr(:,:,:)
  real, allocatable     :: sm_arr(:,:,:)
  ! absorption layer thickness
  real                  :: ds_abs
contains

  subroutine initializeFields()
    implicit none
    if (allocated(ex)) deallocate(ex)
    if (allocated(ey)) deallocate(ey)
    if (allocated(ez)) deallocate(ez)
    if (allocated(bx)) deallocate(bx)
    if (allocated(by)) deallocate(by)
    if (allocated(bz)) deallocate(bz)
    if (allocated(jx)) deallocate(jx)
    if (allocated(jy)) deallocate(jy)
    if (allocated(jz)) deallocate(jz)
    if (allocated(jx_buff)) deallocate(jx_buff)
    if (allocated(jy_buff)) deallocate(jy_buff)
    if (allocated(jz_buff)) deallocate(jz_buff)
    #ifdef oneD
      allocate(ex(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(ey(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(ez(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(by(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(bz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jx_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jy_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(jz_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      allocate(lg_arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST, 0:0, 0:0))
      ! 20 = max # of fields sent in each direction
      sendrecv_offsetsz = NGHOST * 20
      ! 2 (~5) directions to send/recv in 1D
      sendrecv_buffsz = sendrecv_offsetsz * 5
    #elif twoD
      allocate(ex(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(ey(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(ez(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(by(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(bz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jx_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jy_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(jz_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))
      allocate(lg_arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST, 0:0))

     ! 20 = max # of fields sent in each direction
     sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz) * NGHOST * 20
     ! 8 (~10) directions to send/recv in 2D
     sendrecv_buffsz = sendrecv_offsetsz * 10
    #elif threeD
      allocate(ex(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(ey(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(ez(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(bx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(by(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(bz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jx(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jy(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jz(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jx_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jy_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(jz_buff(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                     & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      allocate(lg_arr(-NGHOST : this_meshblock%ptr%sx - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sy - 1 + NGHOST,&
                    & -NGHOST : this_meshblock%ptr%sz - 1 + NGHOST))
      ! 20 = max # of fields sent in each direction
      sendrecv_offsetsz = MAX0(this_meshblock%ptr%sx, this_meshblock%ptr%sy, this_meshblock%ptr%sz)**2 * NGHOST * 20
      ! 26 (~30) directions to send/recv in 3D
      sendrecv_buffsz = sendrecv_offsetsz * 30
    #endif

    allocate(sm_arr(0:this_meshblock%ptr%sx - 1, 0:this_meshblock%ptr%sy - 1, 0:this_meshblock%ptr%sz - 1))

    if (allocated(send_fld)) deallocate(send_fld)
    if (allocated(recv_fld)) deallocate(recv_fld)
    allocate(send_fld(sendrecv_buffsz))
    allocate(recv_fld(sendrecv_offsetsz))
  end subroutine initializeFields
end module m_fields
