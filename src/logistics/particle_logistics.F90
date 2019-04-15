#include "../defs.F90"

module m_particlelogistics
  use m_globalnamespace
  use m_aux
	use m_errors
  use m_domain
  use m_particles
contains
	subroutine copyParticleFromTo(s, p_from, p_to, ti, tj, tk)
		implicit none
		! within a single tile
    integer, intent(in)   :: s, p_from, p_to, ti, tj, tk
		species(s)%prtl_tile(ti, tj, tk)%xi(p_to) = species(s)%prtl_tile(ti, tj, tk)%xi(p_from)
		species(s)%prtl_tile(ti, tj, tk)%yi(p_to) = species(s)%prtl_tile(ti, tj, tk)%yi(p_from)
		species(s)%prtl_tile(ti, tj, tk)%zi(p_to) = species(s)%prtl_tile(ti, tj, tk)%zi(p_from)

		species(s)%prtl_tile(ti, tj, tk)%dx(p_to) = species(s)%prtl_tile(ti, tj, tk)%dx(p_from)
		species(s)%prtl_tile(ti, tj, tk)%dy(p_to) = species(s)%prtl_tile(ti, tj, tk)%dy(p_from)
		species(s)%prtl_tile(ti, tj, tk)%dz(p_to) = species(s)%prtl_tile(ti, tj, tk)%dz(p_from)

		species(s)%prtl_tile(ti, tj, tk)%u(p_to) = species(s)%prtl_tile(ti, tj, tk)%u(p_from)
		species(s)%prtl_tile(ti, tj, tk)%v(p_to) = species(s)%prtl_tile(ti, tj, tk)%v(p_from)
		species(s)%prtl_tile(ti, tj, tk)%w(p_to) = species(s)%prtl_tile(ti, tj, tk)%w(p_from)

		species(s)%prtl_tile(ti, tj, tk)%ind(p_to) = species(s)%prtl_tile(ti, tj, tk)%ind(p_from)
		species(s)%prtl_tile(ti, tj, tk)%proc(p_to) = species(s)%prtl_tile(ti, tj, tk)%proc(p_from)
  end subroutine copyParticleFromTo

  subroutine removeParticleFromTile(s, ti, tj, tk, p)
    implicit none
    integer, intent(in)   :: s, ti, tj, tk, p
    if (p .ne. species(s)%prtl_tile(ti, tj, tk)%npart_sp) then
			call copyParticleFromTo(s, species(s)%prtl_tile(ti, tj, tk)%npart_sp, p, ti, tj, tk)
		end if
    species(s)%prtl_tile(ti, tj, tk)%npart_sp = species(s)%prtl_tile(ti, tj, tk)%npart_sp - 1
  end subroutine removeParticleFromTile

  subroutine createParticle(s, xi, yi, zi, dx, dy, dz, u, v, w, &
		                      & ind, proc)
    implicit none
    integer, intent(in)           :: s
    integer(kind=2), intent(in)   :: xi, yi, zi
    real, intent(in)              :: dx, dy, dz, u, v, w
    integer                       :: p
		integer												:: ti, tj, tk
		integer, optional, intent(in) :: ind, proc
		ti = xi / species(s)%tile_sx + 1
		tj = yi / species(s)%tile_sy + 1
		tk = zi / species(s)%tile_sz + 1
		#ifdef DEBUG
			if ((ti .gt. species(s)%tile_nx) .or. &
				& (tj .gt. species(s)%tile_ny) .or. &
				& (tk .gt. species(s)%tile_nz)) then
				call throwError('ERROR: wrong ti, tj, tk in `createParticle`')
			end if
		#endif
		species(s)%prtl_tile(ti, tj, tk)%npart_sp = species(s)%prtl_tile(ti, tj, tk)%npart_sp + 1
		! FIX0 check overflow
    p = species(s)%prtl_tile(ti, tj, tk)%npart_sp

    species(s)%prtl_tile(ti, tj, tk)%xi(p) = xi
		species(s)%prtl_tile(ti, tj, tk)%dx(p) = dx

		species(s)%prtl_tile(ti, tj, tk)%yi(p) = yi
		species(s)%prtl_tile(ti, tj, tk)%dy(p) = dy

		species(s)%prtl_tile(ti, tj, tk)%zi(p) = zi
		species(s)%prtl_tile(ti, tj, tk)%dz(p) = dz

		species(s)%prtl_tile(ti, tj, tk)%u(p) = u
		species(s)%prtl_tile(ti, tj, tk)%v(p) = v
		species(s)%prtl_tile(ti, tj, tk)%w(p) = w

		if (present(ind) .and. present(proc)) then
			species(s)%prtl_tile(ti, tj, tk)%ind(p) = ind
			species(s)%prtl_tile(ti, tj, tk)%proc(p) = proc
		else
    	species(s)%prtl_tile(ti, tj, tk)%ind(p) = species(s)%cntr_sp
			species(s)%prtl_tile(ti, tj, tk)%proc(p) = mpi_rank
			species(s)%cntr_sp = species(s)%cntr_sp + 1
		end if
  end subroutine createParticle
end module m_particlelogistics
