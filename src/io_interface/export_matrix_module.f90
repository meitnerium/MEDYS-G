Module Export_matrix_module


use basics
use io_module

Implicit None
Private
Public:: Export_Matrix

! pass any options with
type export_matrix_desc
end type


Interface Export_matrix
   module procedure export_integer_matrix
   module procedure export_real_matrix
   module procedure export_complex_matrix
end interface



contains 


  subroutine Export_complex_matrix(Matrix,io_desc,opt_desc)
     Complex(kind=comp_16),intent(in),dimension(:,:):: Matrix
     Type(io_tracker),intent(inout)::io_desc
     type(export_matrix_desc),intent(in),optional::opt_desc
     Integer(kind=int_4)::n,m,i,j
     n=size(matrix(:,1))
     m=size(matrix(1,:))
     ! format for human vizualization 
     Do i=1,n 
        Do j=1,m
           write(*,*) i,j, real(matrix(i,j)), aimag(matrix(i,j))
        End do
     end do
    
  end subroutine
  
  SUBROUTINE Export_real_matrix(Matrix,io_desc,opt_desc)
     Real(kind=real_8),intent(in),dimension(:,:):: Matrix
     Type(io_tracker),intent(inout)::io_desc
     type(export_matrix_desc),intent(in),optional::opt_desc
  integer::i,j
  DO I = LBOUND(matrix,1), UBOUND(matrix,1)
     WRITE(*,*) (matrix(I,J), J = LBOUND(matrix,2), UBOUND(matrix,2))
  END DO
  end subroutine

  SUBROUTINE Export_integer_matrix(Matrix,io_desc,opt_desc)
     Integer(kind=int_4),intent(in),dimension(:,:):: Matrix
     Type(io_tracker),intent(inout)::io_desc
     type(export_matrix_desc),intent(in),optional::opt_desc
  integer::i,j
  DO I = LBOUND(matrix,1), UBOUND(matrix,1)
     WRITE(*,*) (matrix(I,J), J = LBOUND(matrix,2), UBOUND(matrix,2))
  END DO
  end subroutine

End module
