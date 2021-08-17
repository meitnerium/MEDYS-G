Module clock_module
! Purpose: Tools for timing
use precision_module

Implicit None
Private
Public::clock


Type clock
   Logical,private:: started_tag=.False.
   Real(kind=real_8),private:: ST! start time
 contains
Procedure, private:: is_new=>clock_is_new
Procedure, public:: start=>start_clock
Procedure,Public:: reset=>reset_clock
Procedure,Public:: elapsed_time

end type

Contains

Logical function clock_is_new(self)
     Class(clock),intent(in)::self
     If(self%started_tag) then
        clock_is_new=.True.
     else
        clock_is_new=.False.
     end if
end function

subroutine reset_clock(self)
     Class(clock),intent(inout)::self
     if(.not. self%is_new()) then
        self%ST=0.0_real_8
        self%started_tag=.False.
     End if

end subroutine

subroutine start_clock(self)
     Class(clock), intent(inout)::self
     if(.not. self%is_new()) call self%reset()
     call cpu_time(self%ST)
     self%started_tag=.True.
end subroutine

Real(kind=real_8) function elapsed_time(self)
     Class(clock),intent(in)::self
     Real(kind=real_8):: now
     call cpu_time(now)
     elapsed_time= now-self%ST 
end function


End module
