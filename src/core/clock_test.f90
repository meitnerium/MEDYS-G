Module clock_test_module

use test_environment
use clock_module
use string_module
Implicit None
Private
Public::test_clock


contains

subroutine test_clock()
     Implicit None
     type(string):: test_name
     type(clock)::clk
     type(string)::current_test    
     
     test_name="CLOCK"
     call begin_test(test_name)

     current_test="Reset a not started clock"
     call header_1(current_test)
     call clk%reset

     current_test="Start"
     call header_1(current_test)
     call clk%start

     current_test="Elapsed time"
     call header_1(current_test)
     write(*,*) "     ", clk%Elapsed_time()
     

     current_test="Reset a started clock"
     call header_1(current_test)
     call clk%reset
     
     call end_test(test_name)
end subroutine
end module

Program clock_test
use clock_test_module
Implicit None

call test_clock()
End program
