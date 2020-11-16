#include "../defs.F90"

module m_annihilation
#ifdef PAIRANNIHILATION

  use m_globalnamespace
  use m_qednamespace
  use m_aux
  use m_errors
  use m_bincoupling
  use m_particlelogistics

  use m_particlebinning, only: initializePositionBins

  use m_particlebinning, only: positionBin_XYZ
  implicit none

  !--- PRIVATE variables/functions -------------------------------!
  !...............................................................!
contains
  subroutine pairAnnihilation()
    implicit none
    integer   :: ann_electrons(20), ann_positrons(20)
    integer   :: s, s_lec, s_pos
    integer   :: ti, tj, tk

    ! find all the species that participate in the annihilation process
    s_lec = 0; s_pos = 0
    do s = 1, nspec
      if (species(s)%annihilation_sp) then
        if (species(s)%ch_sp .gt. 0) then
          ann_electrons(s_lec + 1) = s
          s_lec = s_lec + 1
        else
          ann_positrons(s_pos + 1) = s
          s_pos = s_pos + 1
        end if
      end if
    end do

    if ((s_lec * s_pos .eq. 0) .and. (s_lec + s_pos .ne. 0)) then
      call throwError("Pair annihilation requires at least 2 species with opposite charges to participate.")
    end if

    call initializePositionBins(species(s)%prtl_tile(ti, tj, tk), position_grid)
    call binParticlePositions(species(s)%prtl_tile(ti, tj, tk), position_grid)
    ...

    call pairAnnihilationOnTile(...)

  end subroutine pairAnnihilation

  subroutine pairAnnihilationOnTile_mc()
    implicit none
    initializePositionBins
  end subroutine pairAnnihilationOnTile_mc
#endif
end module m_annihilation
