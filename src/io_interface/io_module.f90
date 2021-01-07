Module io_module
! Purpose: interface to intel 
!   use precision_module
   use basics
   use string_module
   Implicit None

   Private
   Public:: io_tracker

   Integer(kind=int_4),parameter:: test_unit=53
   Integer(kind=int_4),parameter:: start_unit=101 
   Integer(kind=int_4),parameter:: max_n_unit=500
                                                
type io_tracker
     Logical,private:: opened_tag=.False.
     Integer(kind=int_4)::channel
     Integer(kind=int_4),public::error
!     Integer(kind=int_4),private::warning
  Contains
     Procedure,private,pass::set_unit
     Procedure::get_unit
     Procedure,public,pass::is_opened=>io_tracker_is_opened
     Procedure,public:: open_new_file
     Procedure,public:: open_old_file
     Procedure,public:: close=> close_file
     Procedure:: reset=> reset_io_tracker 
!     Procedure,public:: create_directory
End type

   Contains

   Subroutine set_unit (self)
      Implicit None
      class(io_tracker), intent(inout):: self
      Integer(kind=int_4):: counter
      logical:: i_open
      i_open=.True.
      counter=start_unit
      Do while (i_open)

      Inquire(unit=counter, opened=i_open)
        counter=counter+1
      End do
      self%channel=counter-1
      self%opened_tag=.True.
   End Subroutine set_unit

Integer(int_4) function get_unit(self)
   Implicit None 
   class(io_tracker),intent(in)::self
   get_unit=self%channel
end function
!
logical function io_tracker_is_opened(self) result(res)
Implicit None
!! False => alredy in use, stop or change
!! True => new one, can continue
class(io_tracker),intent(in)::self
If (self%opened_tag) then
    res= .True.
  else
    res=.False.
End if
End function

subroutine reset_io_tracker(self)
      Implicit None
      Class(io_tracker),intent(inout)::self
      self%opened_tag=.False.
      self%channel=0_int_4
end subroutine
!
  Subroutine open_new_file(self,file_name)
      Implicit None
      type(string), intent(in):: file_name
      Class(io_tracker),intent(inout):: self
      logical:: iexist,iopened
      integer(kind=int_8)::ierr
      If(self%is_opened()) Stop 'In open_new_file (io_tracker):: this io_tracker is already in use.'
      inquire(file=file_name%value,exist=iexist,opened=iopened)
      If(iexist) Stop "In open_new_file (io_tracker):: this file already exist."
      If(iopened) Stop "In open_new_file (io_tracker):: this file is already opened."
      call self%set_unit()
      Open(unit=self%channel, file=file_name%value, status='new', iostat=ierr)
      If (ierr .ne. 0) Stop "In open_new_file (io_tracker):: error while opening the file."
   end subroutine 
!
  Subroutine open_old_file(self,file_name)
      Implicit None
      type(string), intent(in):: file_name
      Class(io_tracker),intent(inout):: self
      logical:: iexist,iopened
      integer(kind=int_8)::ierr
      If(self%is_opened()) Stop 'In open_old_file (io_tracker):: this io_tracker is already in use.'
      inquire(file=file_name%value,exist=iexist,opened=iopened)
      If(.not. iexist) Stop "In open_old_file (io_tracker):: this file do not exist."
      If(iopened) Stop "In open_old_file (io_tracker):: this file is already opened."
      call self%set_unit()
      Open(unit=self%channel, file=file_name%value, status='old', iostat=ierr)
      If (ierr .ne. 0) Stop "In open_old_file (io_tracker):: error while opening the file."
   end subroutine 

subroutine close_file(self)
     implicit none
     class(io_tracker),intent(inout)::self
     If (.not. self%is_opened()) Stop "In close_file (io_tracker):: file have not been opened."
     rewind(self%get_unit())
     close(self%get_unit())
     call self%reset()
end subroutine
!
  Subroutine create_directory (self, directory_name)
!#if defined(__INTEL_COMPILER)
    use ifport  ! ifort module
!#endif    
    Implicit None
    Class(io_tracker),intent(inout)::self
    type(string)::directory_name
    logical:: i_exist
!#if defined(__INTEL_COMPILER)
      inquire(directory=directory_name%value,exist=i_exist)
!#elif defined(__gfortran__)
 !     inquire(file=temp_strg,exist=i_exist)
!#endif
      If (.not. i_exist) then
        i_exist = makedirqq(directory_name%value) 
       else
        Stop "In create_directory (io_tracker):: Directory already exist." 
      end if
  End subroutine create_directory 
!

End module 
