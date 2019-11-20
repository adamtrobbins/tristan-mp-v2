#include "../defs.F90"

module m_bincoupling
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_particles
  use m_errors
  implicit none

  type :: spec_ind_pair
    ! object which identifies a particle ...
    !     ... by its species and index on a given tile
    integer :: spec, index
  end type spec_ind_pair

  type :: couple
    ! object that contains two particles ...
    !     ... given by their species and index ...
    !     ... on a single tile
    type(spec_ind_pair) :: part_1
    type(spec_ind_pair) :: part_2
  end type

  !--- PRIVATE functions -----------------------------------------!
  private :: shuffleSet
  !...............................................................!
contains
  ! Knuth's algorithm to randomly shuffle a set
  subroutine shuffleSet(set, set_size)
    implicit none
    integer, intent(in)                 :: set_size
    type(spec_ind_pair), intent(inout)  :: set(set_size)
    type(spec_ind_pair)                 :: temp
    integer                             :: i, j
    do i = 1, set_size - 1
      j = randomInt(dseed, i, set_size + 1)
      temp = set(i)
      set(i) = set(j)
      set(j) = temp
    end do
  end subroutine shuffleSet

  ! this routine pairs particles in two sets #1 and #2 randomly
  !   - two sets can either be equal, or have no intersection
  !   - sets are passed as an array of species (indices) at a given tile
  subroutine coupleParticlesOnTile(ti, tj, tk,&
                                 & sp_arr_1, n_sp_1,&
                                 & sp_arr_2, n_sp_2,&
                                 & coupled_pairs, num_couples)
    implicit none
    integer, intent(in)                     :: ti, tj, tk
    integer, intent(in)                     :: n_sp_1, n_sp_2 ! # of species in set #1 and #2
    integer, intent(in)                     :: sp_arr_1(n_sp_1), sp_arr_2(n_sp_2) ! array of species in set #1 and #2
    integer                                 :: num_1, num_2 ! total # of particles in sets #1 and #2
    integer                                 :: s, si, p, i, j
    type(spec_ind_pair), allocatable        :: set_1(:), set_2(:) ! set #1 and #2 saved as "tuples" of species and index
    integer, intent(out)                    :: num_couples
    type(couple), allocatable, intent(out)  :: coupled_pairs(:)

    ! auxiliary variables
    integer                                 :: common_species
    logical                                 :: pairing_correctQ, same_setsQ

    ! check the number of common elements
    common_species = 0
    do i = 1, n_sp_1
      do j = 1, n_sp_2
        if (sp_arr_1(i) .eq. sp_arr_2(j)) then
          common_species = common_species + 1
        end if
      end do
    end do

    if (common_species .eq. 0) then
      same_setsQ = .false.
    else
      if ((n_sp_1 .eq. n_sp_2) .and. (n_sp_1 .eq. common_species)) then
        same_setsQ = .true.
      else
        call throwError("Sets should be equal or have no intersaction in `coupleParticlesOnTile()`")
      end if
    end if

    if (same_setsQ) then ! if two sets are exactly the same (e.g. gamma+gamma)
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      print *, 'SAME SET PAIRING'
      ! computing number of particles in the set
      do si = 1, n_sp_1
        s = sp_arr_1(si)
        num_1 = num_1 + species(s)%prtl_tile(ti, tj, tk)%npart_sp
      end do
      ! assigning particles in the set
      allocate(set_1(num_1))
      i = 1
      do si = 1, n_sp_1
        s = sp_arr_1(si)
        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
          set_1(i)%spec = s
          set_1(i)%index = p
          i = i + 1
        end do
      end do
      ! shuffle the set
      call shuffleSet(set_1, num_1)
      num_couples = INT(num_1 / 2)
      allocate(coupled_pairs(num_couples))
      ! assign pairs
      do i = 1, num_couples
        coupled_pairs(i)%part_1 = set_1(i)
        coupled_pairs(i)%part_2 = set_1(num_1 - i + 1)
      end do
    else ! if two sets have no common elements (e.g. compton scattering)
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      print *, 'DIFFERENT SET PAIRING'
      ! computing number of particles in set #1
      do si = 1, n_sp_1
        s = sp_arr_1(si)
        num_1 = num_1 + species(s)%prtl_tile(ti, tj, tk)%npart_sp
      end do
      ! computing number of particles in set #2
      do si = 1, n_sp_2
        s = sp_arr_2(si)
        num_2 = num_2 + species(s)%prtl_tile(ti, tj, tk)%npart_sp
      end do
      ! assigning particles in set #1
      allocate(set_1(num_1))
      i = 1
      do si = 1, n_sp_1
        s = sp_arr_1(si)
        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
          set_1(i)%spec = s
          set_1(i)%index = p
          i = i + 1
        end do
      end do
      ! assigning particles in set #2
      allocate(set_2(num_2))
      i = 1
      do si = 1, n_sp_2
        s = sp_arr_2(si)
        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp
          set_2(i)%spec = s
          set_2(i)%index = p
          i = i + 1
        end do
      end do
      ! now we can simply work with `set_1` and `set_2`
      num_couples = min(num_1, num_2)
      if ((num_1 .eq. 1) .and. (num_2 .eq. 1)) then
        ! check most simple case
        if ((set_1(num_1)%spec .eq. set_2(num_2)%spec) .and.&
          & (set_1(num_1)%index .eq. set_2(num_2)%index)) then
          ! cannot pair a single particle to itself
          num_couples = 0
        else
          ! trivial pairing
          allocate(coupled_pairs(num_couples))
          coupled_pairs(1)%part_1 = set_1(num_1)
          coupled_pairs(1)%part_2 = set_2(num_2)
        end if
      else
        ! if pairing is non trivial
        allocate(coupled_pairs(num_couples))
        pairing_correctQ = .false.
        do while (.not. pairing_correctQ)
          ! on average this will be performed `e~3` times ...
          !     ... (if there are common elements in set #1 and #2)
          call shuffleSet(set_1, num_1)
          call shuffleSet(set_2, num_2)
          pairing_correctQ = .true.
          do i = 1, num_couples
            if ((set_1(i)%spec .eq. set_2(i)%spec) .and.&
              & (set_1(i)%index .eq. set_2(i)%index)) then
              ! check if the particle is coupled to itself
              pairing_correctQ = .false.
              exit
            end if
          end do
        end do
        do i = 1, num_couples
          coupled_pairs(i)%part_1 = set_1(i)
          coupled_pairs(i)%part_2 = set_2(i)
        end do
      end if
    end if

    if (allocated(set_1)) deallocate(set_1)
    if (allocated(set_2)) deallocate(set_2)
  end subroutine

end module m_bincoupling
