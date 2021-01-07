module matrix_rearrangement_module
!Purpose: Module containing all procedures used to reorder matrices given by Columbus (only molecular orbitals). 
!3 possible cases: 1D(AO(or nk),MO), 2D(MO,MO) and 4D(MO,MO,MO,MO)

use basics
use string_module,only:string
use io_module
use error_module

implicit none
private
public::matrix_rearrangement

namelist/matrix_reordering_dim/reorder_size
  integer(kind=int_4)::reorder_size
namelist/matrix_reordering/order_modif_vector
  integer(kind=int_4),dimension(:),allocatable::order_modif_vector

type::matrix_rearrangement

  type(err_tracker)::error
    integer(kind=int_4)::reorder_size ! Size (or dimension) of the modification vector (number of MO wanted)    
    integer(kind=int_4),dimension(:),allocatable::order_modif_vector ! Vector giving the wanted MO order
    
  contains
  
    procedure,private::read_nml=>read_matrix_rearrangement_namelist ! Reading procedure for the reordering namelist
    procedure,public::matrix_1D_reordering=>matrix_MO_reordering_case_1D_real ! Reordering procedure

end type

contains

subroutine read_matrix_rearrangement_namelist(self,filename,order_modif_vector)
  implicit none

  class(matrix_rearrangement),intent(inout)::self
  type(string),intent(in)::filename
  type(io_tracker)::info

  integer(kind=int_4),dimension(:),allocatable,intent(out)::order_modif_vector
  character(len=16)::fname2="vector_array.nml"
  integer(kind=int_4)::n,nml_error_size,nml_error,i

  call info%open_old_file(filename)
  write(*,*)'importing size namelist'			!TEST
  n=info%get_unit()
  write(*,*)'n_unit = ',n
  read(n,nml=matrix_reordering_dim,iostat=nml_error_size)

  write(*,*)'Checkpoint = == == = ='			!TEST
  rewind(n)
  call info%close()
  call info%reset()

  if (nml_error_size .eq. 0) then
    self%reorder_size=reorder_size
  else
    write(*,*)'Iostat for size is :',nml_error_size
    call self%error%set(nml_error_size,"Error while reading reorder_size's namelist")
    !TODO changed for test purpose, what is this file !
    self%reorder_size=10
    write(*,*) 'ERROR MESSAGE for size reading'
  end if

write(*,*)''
write(*,*)'Reorder size = ',self%reorder_size
write(*,*)''

  allocate (order_modif_vector(self%reorder_size))

  open(113, file=fname2)
  read(113,fmt=*,iostat=nml_error) (order_modif_vector(i),i=1,self%reorder_size)

  close(113)
  write(*,*) order_modif_vector
  write(*,*)'iostat is : ',nml_error

!  call info%open_old_file(filename)
!  write(*,*)'importing reordering namelist'
!  n=info%get_unit()
!write(*,*)'n_unit = ',n
!  read(n,nml=matrix_reordering,iostat=nml_error)
!  rewind(n)
!  call info%close()
!  call info%reset()
!write(*,*)'The iostat value for order_modif_vector is ', nml_error

  if (nml_error.eq.0) then
    self%order_modif_vector=order_modif_vector
    write(*,*)'File read'
  else
    call self%error%set(nml_error,"Error while reading matrix_reordering's namelist")
    write(*,*) 'ERROR MESSAGE for order_modif_vector'
  end if

end subroutine read_matrix_rearrangement_namelist

subroutine matrix_MO_reordering_case_1D_real(self,matrix_MO1,filename)
  implicit none

  class(matrix_rearrangement)::self
  type(string),intent(in)::filename

  real(kind=real_8),dimension(:,:),intent(inout)::matrix_MO1
  real(kind=real_8),dimension(:,:),allocatable::matrix_MO2
  integer(kind=int_4)::dim_1,dim_2,i,j

  call self%read_nml(filename,self%order_modif_vector) 

  dim_1=size(matrix_MO1(:,1))
  dim_2=size(matrix_MO1(1,:))
  allocate(matrix_MO2(dim_1,dim_2))

  write(*,*)'dimension dim_1 et dim_2: ', dim_1,dim_2			!TEST

  open(12,file='ordermodifvect.dat',status='unknown')
  open(15,file='modified_matrix.dat',status='unknown')
write(*,*)'matrix_MO loop'						!TEST
  do i=1,dim_2
    matrix_MO2(:,i)=matrix_MO1(:,self%order_modif_vector(i))
    write(15,*) 'MO # :',i 
    write(15,*)''
    do j=1,dim_1
      write(15,*) matrix_MO2(j,i),j
    end do
  write(12,*) self%order_modif_vector(i)
  end do
  close(12)
  close(15)

  matrix_MO1=matrix_MO2
  deallocate(matrix_MO2)

end subroutine matrix_MO_reordering_case_1D_real

end module matrix_rearrangement_module

