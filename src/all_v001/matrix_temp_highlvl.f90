module matrix_temp_highlvl
  use basics
  use matrix_temp_lowlvl
  use io_module
  use string_module
!****************************************************************************************************************************
!purpose: verifying if the matrices generated are the same as the ones in the previous version of MEDYS. 
!****************************************************************************************************************************

  implicit none

  contains

subroutine compare_qpint_to_v002(char_filename1,char_filename2)
!Purpose: Compares the results of matrices from version 1 and 2.
  implicit none

    type(io_tracker)::info
    type(string)::filename1
    type(string)::filename2
  
    character(len=*),intent(in)::char_filename1,char_filename2
    complex(kind=comp_16),dimension(:,:),allocatable::matrix_v001
    complex(kind=comp_16),dimension(:,:),allocatable::matrix_v002
    integer(kind=int_4)::n_unit,l,m,n,k
    integer(kind=int_4),dimension(2)::smat1,szmat2

    filename1=char_filename1
    filename2=char_filename2

    call info%open_old_file(filename1) 
    n_unit=info%get_unit()

    call file_to_zmatrix(n_unit,matrix_v001)
    call info%close()
    call info%reset()

    call info%open_old_file(filename2) 
    n_unit=info%get_unit()

    call file_to_zmatrix(n_unit,matrix_v002)
    call info%close()

    smat1=shape(matrix_v001)
    m=smat1(1)
    n=smat1(2)
    allocate(matrix_v001(m,n))

    szmat2=shape(matrix_v002)
    k=szmat2(1)
    l=szmat2(2)
    allocate(matrix_v002(k,l))

    write(*,*) 'Shape of version 1 - Shape of version 2 =',smat1-szmat2 ! Shows results in terminal

end subroutine compare_qpint_to_v002

subroutine export_supermatrix4D_to_matrix2D (matrix4D,char_filename) 
! Purpose: reshapes and exports a matrix from Etienne's code (MEDYS version v0.0.1)
  implicit none

    type(io_tracker)::info
    type(string)::filename
    
    character(len=*),intent(in)::char_filename
    complex(kind=comp_16),dimension(:,:,:,:),intent(in)::matrix4D
    complex(kind=comp_16),dimension(:,:),allocatable::matrix2D_v001
    integer(kind=int_4)::n_unit,m,n
    integer(kind=int_4),dimension(2)::s2Dmat1

    filename=char_filename

    call info%open_new_file(filename)
    n_unit=info%get_unit()

    call reshape_matrix_4D_to_2D(matrix4D,matrix2D_v001)
    s2Dmat1=shape(matrix2D_v001)
    m=s2Dmat1(1)
    n=s2Dmat1(2)
    allocate(matrix2D_v001(m,n))

    call zmatrix_to_file(n_unit,matrix2D_v001)
    call info%close()
    call info%reset()

end subroutine export_supermatrix4D_to_matrix2D

subroutine export_supermatrix3D_to_matrix2D (matrix3D,char_filename)
!Purpose: reshapes and exports a matrix from MEDYS version v0.0.2
  implicit none
   
    type(io_tracker)::info
    type(string)::filename

    character(len=*),intent(in)::char_filename
    complex(kind=comp_16),dimension(:,:,:),intent(in)::matrix3D
    complex(kind=comp_16),dimension(:,:),allocatable::matrix2D_v002
    integer(kind=int_4),dimension(2)::s2Dmat2
    integer(kind=int_4)::n_unit,k,l

    filename=char_filename

    call info%open_new_file(filename) 
    n_unit=info%get_unit()

    call reshape_matrix_3D_to_2D(matrix3D)
    s2Dmat2=shape(matrix2D_v002)
    k=s2Dmat2(1)
    l=s2Dmat2(2)
    allocate(matrix2D_v002(k,l))

    call zmatrix_to_file(n_unit,matrix2D_v002)
    call info%close() 
    call info%reset()

end subroutine export_supermatrix3D_to_matrix2D

end module matrix_temp_highlvl
