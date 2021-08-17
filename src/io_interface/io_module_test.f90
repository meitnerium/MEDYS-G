Module io_module_test

use io_module
use test_environment
Implicit None
Private
Public:: test_io_module

Contains

subroutine test_io_module()
  Implicit none
     type(string):: test_name
     type(string)::current_test
     type(io_tracker)::tracker1
     type(string):: file_name

    test_name="IO_TRACKER"
     call begin_test(test_name)

     current_test="open_new_file"
     call header_1(current_test)
     file_name="test.file"
     call tracker1%open_new_file(file_name)
     
     current_test="close_file"
     call header_1(current_test)
     call tracker1%close()
 
     current_test="open_old_file"
     call header_1(current_test)
     file_name="test.file"
     call tracker1%open_old_file(file_name)

     current_test="open while opened."
     call header_1(current_test)
     call tracker1%open_old_file(file_name)

     call end_test(test_name)

end subroutine

End module


Program io_module_test_program
use io_module_test

call test_io_module()
end program
