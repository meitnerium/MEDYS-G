Module get_command_line_module
! purpose: catch the command line input

use iso_varying_string
Implicit None
Private
Public:: get_input_name
Contains

subroutine get_input_name (input_file,error)
! warning: assumes that only the input file name is passed in command line.
  type(varying_string),intent(inout)::input_file
  Integer:: i
  CHARACTER(len=500) :: arg
  logical:: founded
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read command line arguments !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test        
i = 0
founded=.False.
DO
   CALL get_command_argument(i, arg)
   IF (LEN_TRIM(arg) == 0) EXIT
   If(i==1) then
     input_file=arg
     founded=.True.
   endif
   i = i+1
END DO
if(founded.eq..false.) Stop "The input file was not given as command line entry"
end subroutine

End Module
