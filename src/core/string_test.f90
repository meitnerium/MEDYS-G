Module string_module_test
!Purpose: test the string module

use basics
use string_module
Implicit none


Character(len=*),parameter:: char1="This is a sentence."
type(string)::strg1,strg2,strg3


contains

subroutine test_assign_string()
Implicit None
write(*,*) ""
write(*,*) "Begin of test_assign_string"
write(*,*) ""
strg1=char1
strg2=strg1
write(*,*) char1
write(*,*) strg1%value  
write(*,*) strg2%value  
write(*,*) ""
write(*,*) "End of test_assign_string"
write(*,*) ""

end subroutine

subroutine test_add_string()

write(*,*) ""
write(*,*) "Begin test_add_string"
write(*,*) ""

strg3=strg1+strg2
write(*,*) strg3%value

write(*,*) ""
write(*,*) "End test_add_string"
write(*,*) ""



end subroutine

subroutine test_len_string()

write(*,*) ""
write(*,*) "Begin of test_len_string"
write(*,*) ""

write(*,*) strg1%len()
write(*,*) ""
write(*,*) "End of test_len_string"
write(*,*) ""
end subroutine

End module

Program test_string
use string_module_test

call test_assign_string()

call test_len_string()
call test_add_string()
End program
