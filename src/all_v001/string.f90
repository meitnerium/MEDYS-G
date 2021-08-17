Module String_module
!Purpose: define a string object
Use basics
Implicit None
Private
Public:: string

! !!!!!!!!!!
Character(len=*),parameter,public:: string_module_version="00.00.001"
Character(len=*),parameter,public:: string_module_date="8 mai 2014"
! !!!!!!!!!!



type String
     logical:: allocated_tag=.False.
     Integer(kind=int_4)::err=0
     Integer(kind=int_4)::length
     Character(len=:), allocatable::value
        Contains
     procedure,private::is_allocated
     generic,public:: assignment(=)=> assign_string_from_string, assign_string_from_character
     procedure, private:: assign_string_from_string
     procedure, private:: assign_string_from_character 
     generic,public:: len=> len_string
     procedure,private,pass:: len_string
     generic,public:: operator(+)=> add_string_to_string!,add_string_to_character,add_character_to_string,add_character_to_character
     procedure,private,pass:: add_string_to_string     
!     procedure,private:: add_string_to_character     
!     procedure,private:: add_character_to_string  
!     procedure,private:: add_character_to_character     
!     procedure,public::reset=>reset_string
     generic,public:: operator(.eq.)=> compare_string_to_character,compare_string_to_string
     procedure,private:: compare_string_to_character
     procedure,private:: compare_string_to_string
End type


   Contains

logical function is_allocated(self)
   Implicit None
   Class(string),intent(in)::self
   If(self%allocated_tag) then
      is_allocated=.True.
   else
      is_allocated=.False.
   end if
end function

subroutine assign_string_from_character(string_out,string_in)
Implicit None
   Class(string),intent(out)::string_out
   Character(len=*),intent(in)::string_in
   allocate(string_out%value,source=trim(adjustl(string_in)))
   string_out%allocated_tag=.True.

End subroutine

subroutine assign_string_from_string(string_out,string_in)
Implicit None
   Class(string),intent(out)::string_out
   Type(string),intent(in)::string_in
   allocate(string_out%value,source=trim(adjustl(string_in%value)))
   string_out%allocated_tag=.True.
End subroutine

 type(string) function add_string_to_string(strg1,strg2) result(out_strg)
  Implicit None
   Class(string),intent(in)::strg1
   type(string),intent(in)::strg2
   Integer(kind=int_4)::n
   n=strg1%len()+strg2%len()
!   allocate(character(len=n)::out_strg%value)
   call allocate_string(out_strg,n)
   out_strg%value=strg1%value//strg2%value
end function 

Integer(kind=int_4) function len_string(strg) result(res)
   Implicit None
   Class(string),intent(in)::strg
   res=len(strg%value)
end function

subroutine allocate_string(strg,n)
   Implicit none
   Class(string),intent(inout)::strg
   Integer(kind=int_4)::n

   If(allocated(strg%value)) deallocate(strg%value)
   allocate(character(len=n)::strg%value)
   strg%allocated_tag=.False.   
end subroutine

logical function compare_string_to_character(self,char1) result(res)
     Implicit None
     class(string),intent(in)::self
     character(len=*),intent(in)::char1
     res= self%value.eq.Trim(adjustl(char1))
end function
logical function compare_string_to_string(self,strg1) result(res)
     Implicit None
     class(string),intent(in)::self,strg1
     res= self%value.eq.strg1%value
end function

logical function compare_character_to_string(char1,self) result(res)
     Implicit None
     class(string),intent(in)::self
     character(len=*),intent(in)::char1
     res= self%value.eq.Trim(adjustl(char1))
end function

End Module
