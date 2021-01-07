Module error_module

use precision_module,only:int_4
use string_module

Implicit None
Private
Public::err_tracker

type::err_tracker
  logical:: Q_tag=.False.
  integer(kind=int_4)::critical
  integer(kind=int_4)::warning
  type(string)::message
  contains
 procedure:: set=>set_error
! procedure:: get=>get_error
 procedure:: get_message
 procedure:: Q


end type

Contains

subroutine set_error(self,iostat,errmsg)
 Class(err_tracker),intent(inout)::self
integer(kind=int_4),intent(in)::iostat
 character(len=*),intent(in)::errmsg
if(iostat.ne.0) then
    self%warning=iostat
    self%message=errmsg
    self%Q_tag=.True.
end if

end subroutine

function get_message(self) result(res)
Class(err_tracker),intent(in)::self
Character(len=:),allocatable::res
allocate(res,source=self%message%value)
end function

logical function Q(self)
Class(err_tracker),intent(in)::self
Q=self%Q_tag
end function

end module
