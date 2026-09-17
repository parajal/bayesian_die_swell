! This is a fallback routine for compilers that do not support the fortran 2008
! intrinsic execute_command_line, but do have system, in particular ifort 14
! and older. In ifort 15 execute_command_line is supported.
! Note: optional arguments of execute_command_line are not supported!!
! Note: nagfor supports execute_command_line starting with 3.5.1, but
!       for older versions the use line needs to be uncommented for this
!       fallback routine to work.
! Note: NOT in a module, but ends up in the lib as .o file and is only linked
!       if the intrinsic execute_command_line does not exist.
subroutine execute_command_line ( cmd )
!  use F90_UNIX_PROC
  character(len=*), intent(in) :: cmd
  call system ( cmd )
end subroutine execute_command_line
