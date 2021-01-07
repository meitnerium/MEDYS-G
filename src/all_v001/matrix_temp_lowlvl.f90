module matrix_temp_lowlvl
  use basics

  implicit none
  
  contains

subroutine reshape_matrix_3D_to_2D (matrix3D)
  implicit none
    complex(kind=comp_16),intent(in),dimension(:,:,:)::matrix3D 
    complex(kind=comp_16),dimension(:,:),allocatable::zmatrix
    integer(kind=int_4),dimension(3)::s3D
    integer(kind=int_4),dimension(2)::s2D
    integer(kind=int_4)::m,n,l,zm,zn !opened_unit
    
    s3D=shape(matrix3D) 
    m=s3D(1)
    n=s3D(2)
    l=s3D(3)

    write(*,*) 'Shape of matrix 3D before reshaping:',s3D

    s2D=shape(zmatrix)
    zm=s2D(1)
    zn=s2D(2)
!    zm=m+size(m)*(n-1)
!    m=mod(zm,size(m))
!    n=((zm)/(size(m))) + 1
    zn=l
    allocate(zmatrix(zm,zn))
    
    zmatrix=reshape(matrix3D,(/zm,zn/))
    write(*,*) 'Shape of new matrix 2D after reshaping:',s2D

end subroutine reshape_matrix_3D_to_2D

subroutine reshape_matrix_4D_to_2D (matrix4D,zmatrix)
  implicit none

    complex(kind=comp_16),dimension(:,:,:,:),intent(in)::matrix4D
    complex(kind=comp_16),dimension(:,:),allocatable,intent(out)::zmatrix
    integer(kind=int_4),dimension(4)::s4D
    integer(kind=int_4),dimension(2)::s2D
    integer(kind=int_4)::l,m,n,k,zm,zn
!    real(kind=real_8)::l,m,n,k,zm,zn
    
    s4D=shape(matrix4D) 
    m=s4D(1)
    n=s4D(2)
    l=s4D(3)
    k=s4D(4)

    write(*,*) 'Shape of matrix 4D before reshaping:',s4D

    s2D=shape(zmatrix)
    zm=s2D(1)
    zn=s2D(2)
    
!    zm=m+size(m)*(n-1)+size(m)*size(n)*(l-1)
!    m=mod(zm,size(m))
!    n=1+((mod(zm,size(n))-m)/size(m))
!    l=1+((zm-m-(size(m)*(n-1)))/size(m)*size(n))
    zn=k
    allocate(zmatrix(zm,zn))
    
    zmatrix=reshape(matrix4D,(/zm,zn/))
    write(*,*) 'Shape of new matrix 2D after reshaping:',s2D


end subroutine reshape_matrix_4D_to_2D

subroutine file_to_zmatrix (opened_unit,matrix)
  implicit none
  ! Purpose: reads the values contained into a matrix (in a file).
  ! warning: assume that the right file has been correctly opened with unit #opened_unit.

    complex(kind=comp_16),intent(inout),dimension(:,:),allocatable::matrix
    integer(kind=int_4),intent(in)::opened_unit
    integer(kind=int_4)::i,j,ierr,m,n
    real(kind=real_8)::re,im

    read(opened_unit,'(2i4)') m,n
    allocate(matrix(m,n))
    ierr=0

      do while (ierr /= -1)
         read(opened_unit,'(2i4,E30.20)',iostat=ierr) i,j,re,im
         matrix(i,j)=cmplx(re,im)
      end do

end subroutine file_to_zmatrix

subroutine zmatrix_to_file (opened_unit,matrix)
  implicit none
  ! purpose: writes the values contained in a matrix into a file.
  ! warning: assume that the right file has been correctly opened with unit #opened_unit.

    complex(kind=comp_16),intent(in),dimension(:,:)::matrix
    integer(kind=int_4),intent(in)::opened_unit
    integer(kind=int_4)::n,m,i,j
    integer(kind=int_4),dimension(2)::s

    s=shape(matrix)
    m=s(1)
    n=s(2)
    write(opened_unit,'(2i4)') m,n
     do i=1,n
       do j=1,m
          write(opened_unit,'(2i4,E30.20)') i, j, Real(matrix(i,j)),aimag(matrix(i,j))
       end do
     end do
  
end subroutine zmatrix_to_file

end module matrix_temp_lowlvl
