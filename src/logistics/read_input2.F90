#include "../defs.F90"

module m_readinput
  use m_globalnamespace
  implicit none

  !--- PRIVATE functions -----------------------------------------!
  private :: str2int, str2real
  !...............................................................!
contains
  ! read input/output filename/directory
  subroutine readCommandlineArgs()
    implicit none
    integer                 :: i
    character(len=STR_MAX)  :: arg, arg1
    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      select case (arg)
        case ('-i', '--input')
          call get_command_argument(i + 1, arg1)
          input_file_name = trim(arg1)
        case ('-o', '--output')
          call get_command_argument(i + 1, arg1)
          output_dir_name = trim(arg1)
        case ('-r', '--restart')
          call get_command_argument(i + 1, arg1)
          restart_dir_name = trim(arg1)
      end select
    end do
  end subroutine readCommandlineArgs

  ! read float number
  function getIntInput(blockname, varname, def_val) result(val)
    implicit none
    character(len=*), intent(in)  :: blockname, varname
    integer, intent(in)           :: def_val
    integer                       :: val

    integer                       :: iostatus, k1 = 0, k2 = 0
    character(len=STR_MAX)        :: istream, value_str

    ! opening the input file
    open(unit = UNIT_input, file = trim(input_file_name),&
          & action = 'read', IOSTAT = iostatus)
    if (iostatus /= 0) stop 'ERROR: while reading input file'
    find_blockname: do while (.true.)
      read(UNIT_input, *, IOSTAT = iostatus) istream
      if (iostatus .lt. 0) exit
      if (trim(istream) .eq. '<'//trim(blockname)//'>') then
        ! found the right <blockname>
        find_varname: do while (.true.)
          read(UNIT_input, *, IOSTAT = iostatus) istream
          if (iostatus .lt. 0) exit find_blockname
          if (trim(istream) .eq. trim(varname)) then
            ! found the right variable
            backspace(UNIT_input)
            read(UNIT_input, '(a)') istream
            k1 = scan(istream, '=')
            k2 = scan(istream, '#')
            if (k2 .eq. 0) k2 = STR_MAX
            value_str = istream(k1+1:k2-1)
            value_str = adjustl(value_str)
            exit find_blockname
          end if
        end do find_varname
      end if
    end do find_blockname
    close(UNIT_input)

    if (k2 .eq. 0) then
      val = def_val
    else
      call str2int(value_str, val, iostatus)
      if (iostatus /= 0) then
        val = def_val
      end if
    end if
  end function getIntInput

  function getRealInput(blockname, varname, def_val) result(val)
    implicit none
    character(len=*), intent(in)  :: blockname, varname
    real(mprec), intent(in)       :: def_val
    real(mprec)                   :: val

    integer                       :: iostatus, k1 = 0, k2 = 0
    character(len=STR_MAX)        :: istream, value_str

    ! opening the input file
    open(unit = UNIT_input, file = trim(input_file_name),&
          & action = 'read', IOSTAT = iostatus)
    if (iostatus /= 0) stop 'ERROR: while reading input file'
    find_blockname: do while (.true.)
      read(UNIT_input, *, IOSTAT = iostatus) istream
      if (iostatus .lt. 0) exit
      if (trim(istream) .eq. '<'//trim(blockname)//'>') then
        ! found the right <blockname>
        find_varname: do while (.true.)
          read(UNIT_input, *, IOSTAT = iostatus) istream
          if (iostatus .lt. 0) exit find_blockname
          if (trim(istream) .eq. trim(varname)) then
            ! found the right variable
            backspace(UNIT_input)
            read(UNIT_input, '(a)') istream
            k1 = scan(istream, '=')
            k2 = scan(istream, '#')
            if (k2 .eq. 0) k2 = STR_MAX
            value_str = istream(k1+1:k2-1)
            value_str = adjustl(value_str)
            exit find_blockname
          end if
        end do find_varname
      end if
    end do find_blockname
    close(UNIT_input)
    val = def_val

    if (k2 .eq. 0) then
      val = def_val
    else
      call str2real(value_str, val, iostatus)
      if (iostatus /= 0) then
        val = def_val
      end if
    end if
  end function getRealInput

  subroutine str2int(str, int, stat)
    implicit none
    character(len=*), intent(in) :: str
    integer, intent(out)         :: int
    integer, intent(out)         :: stat
    read(str, *, iostat = stat) int
  end subroutine str2int
  subroutine str2real(str, spr, stat)
    implicit none
    character(len=*), intent(in) :: str
    real(mprec), intent(out)     :: spr
    integer, intent(out)         :: stat
    read(str, *, iostat = stat) spr
  end subroutine str2real
end module m_readinput
