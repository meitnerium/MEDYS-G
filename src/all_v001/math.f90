Module math
! Contains common mathematical functions and operations used in molecular dynamics programs
use Basics
!use blas95
!use lapack95
!use f95_precision
use blas95
use lapack95
use f95_precision


contains

!****************************************************************************************
!****************************************************************************************
function chop(x)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
real(kind=real_8):: chop
real(kind=real_8),intent(in) :: x

if(dabs(x).lt.1.d-20)then
 chop=0.d0
else
 chop = x
endif

end function

!****************************************************************************************
!****************************************************************************************
subroutine MATexp(a,expa,lda,sfact)
! computes the exponential "expa" of a complex hermitian matrix "a" (times a scalar factor, "sfact").
! Note: the hermitian character of "a" is required to use the unitary property of its 
! eigenvectors (otherwise the computation would yield wrong results)
!****************************************************************************************
!****************************************************************************************

use basics

implicit none
integer(kind=int_4), intent(inout) :: lda
complex(kind=comp_16), intent(inout), dimension(lda,lda) :: a
complex(kind=comp_16), intent(out), dimension(lda,lda) :: expa
complex(kind=comp_16), intent(in) :: sfact

complex(kind=comp_16), allocatable, dimension(:,:) :: evec
complex(kind=comp_16), allocatable, dimension(:) :: eval
complex(kind=comp_16), allocatable, dimension(:,:) :: AIdiag,tmp
integer(kind=int_4) :: i,j

allocate(evec(lda,lda))
allocate(eval(lda))

 call Diagonalize(a,lda,evec,eval)

allocate(AIdiag(lda,lda))

!$OMP PARALLEL DO PRIVATE(j)
do i=1,lda
do j=1,lda
AIdiag(i,j)=dcmplx(0.d0,0.d0)
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do i=1,lda
AIdiag(i,i)=cdexp(sfact*eval(i))
enddo
!$OMP END PARALLEL DO

allocate(tmp(lda,lda))
!expa = MatmulP(evec,matmulP(AIdiag,dconjg(Transpose(evec)),lda,lda,lda),lda,lda,lda)

tmp = MatmulPAdj(AIdiag,evec,lda,lda,lda)
expa = MatmulP2(evec,tmp,lda,lda,lda)

deallocate(tmp)
deallocate(AIdiag)
 call cchop(expa,lda,lda)

end subroutine MATexp

!****************************************************************************************
! A parallel of matmul interfaced within the current program: matrix multiplication of Matrix M1 (n x m) and 
! M2 (m x l); the result being of dimension (n x l)
function MatMulP2(M1,M2,n,m,l)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
integer(kind=int_4), intent(in) :: n,m,l
complex(kind=comp_16), dimension(n,l) :: MatMulP2
complex(kind=comp_16), intent(inout), dimension(n,m) :: M1
complex(kind=comp_16), intent(inout), dimension(m,l) :: M2
complex(kind=comp_16) :: alpha,beta
integer(kind=int_4) :: nc,mc,lc
complex(kind=comp_16), allocatable, dimension(:,:) :: tmp
nc=n
mc=m
lc=l
alpha=dcmplx(1.d0,0.d0)
beta=dcmplx(0.d0,0.d0)
allocate(tmp(n,l))

tmp(:,:)=dcmplx(0.d0,0.d0)
 !call zgemm('N','N',nc,lc,mc,alpha,M1,nc,M2,mc,beta,tmp,nc)
 call gemm(M1,M2,tmp)

MatMulP2(:,:) = tmp(:,:)
end function

!****************************************************************************************
! A parallel of matmul interfaced within the current program: 
! matrix multiplication of Matrix M1 (n x m) and the adjoint of M2 (l x m); the result being of dimension (n x l)
function MatMulPAdj(M1,M2,n,m,l)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
integer(kind=int_4), intent(in) :: n,m,l
complex(kind=comp_16), dimension(n,l) :: MatMulPAdj
complex(kind=comp_16), intent(inout), dimension(n,m) :: M1
complex(kind=comp_16), intent(inout), dimension(l,m) :: M2
complex(kind=comp_16) :: alpha,beta
integer(kind=int_4) :: nc,mc,lc
complex(kind=comp_16), allocatable, dimension(:,:) :: tmp
nc=n
mc=m
lc=l
alpha=dcmplx(1.d0,0.d0)
beta=dcmplx(0.d0,0.d0)
allocate(tmp(n,l))

tmp(:,:)=dcmplx(0.d0,0.d0)

 !call zgemm('N','C',nc,lc,mc,alpha,M1,lc,M2,mc,beta,tmp,nc)
 call gemm(M1,M2,tmp)
MatMulPAdj(:,:) = tmp(:,:)
end function

!****************************************************************************************
!****************************************************************************************
Subroutine Diagonalize(a,lda,evec,eval)
! Computes the eigenvalues and eigenvectors of a hermitian complex matrix
!****************************************************************************************
!****************************************************************************************

use basics

IMPLICIT NONE
!variables for zheev routine in LAPACK
Integer(kind=int_4), intent(in) :: lda !dimensionality of the problem
complex(kind=comp_16), intent(inout), dimension(lda,lda) :: a !input matrix to diagonalize, gets replaced with matrix of eigenvectors
complex(kind=comp_16), allocatable, dimension(:) :: work !working array for zheev
real(kind=real_8), allocatable, dimension(:) :: rwork !working array for zheev
real(kind=real_8), allocatable, dimension(:) :: w !array containing the eigenvalues
Integer(kind=int_4) :: info !dimension of the allocatable working array and information output parameter
Integer :: lwork

!variables introduced to interface hgeev with the current program
complex(kind=comp_16), allocatable, dimension(:,:) :: aaa !the working matrix, this avoids the input matrix a to be overwritten
complex(kind=comp_16), intent(out), dimension(lda) :: eval !eigenvalues will be copied in this array
complex(kind=comp_16), intent(out), dimension(lda,lda) :: evec !column-eigenvectors will be copied in this array
! other variables
Integer(kind=int_4) :: i,j
!memory allocation
allocate(aaa(lda,lda))
allocate(w(lda))
allocate(rwork(lda*3-2))

!make a copy of the input matrix into an internal working matrix
w=0.d0
rwork=0.d0
info=0
!do i=1,lda
!do j=1,lda
!aaa(i,j)=a(i,j)
!enddo
!enddo
aaa=a

!first run with a negative value of lwork to determine the size of working arrays

allocate(work(1))
work=dcmplx(0.d0,0.d0)
lwork=-1
! call zgeev('V','N',lda,aaa,lda,w,vl,lda,vr,lda,work,lwork,rwork,info)
! call zheev('N','U',lda,aaa,lda,w,work,lwork,rwork,info)
!write(*,*) "TEST",lda,aaa,lda,w,work,lwork,rwork,info
 call zheev('V','U',lda,aaa,lda,w,work,lwork,rwork,info)

! call heev(aaa,w,'V','U',info)
!write(*,*)'lda',lda
!write(*,*)'w',w
!write(*,*)'work',work
!write(*,*)'lwork',lwork
!write(*,*)'rwork',rwork
!write(*,*)'info',info
!STOP'prezheev'

if(info.ne.0)then
write(*,*)'Error initializing length of work arrays for matrix diagonalization. info=',info
write(*,*)'Calculation aborted',lwork
goto 22
endif

aaa=a
lwork=idint(dreal(work(1)))
!write(*,*)'lwork',lwork
!write(*,*)'w',w
!write(*,*)'info',info
deallocate(work)
!second run for the true calculation, lwork being optimally determined
allocate(work(lwork))
work=dcmplx(0.d0,0.d0)
! call zgeev('V','N',lda,aaa,lda,w,vl,lda,vr,lda,work,lwork,rwork,info)
! call zheev('N','U',lda,aaa,lda,w,work,lwork,rwork,info)
info=0
!write(*,*)lda
!do i=1,lda
!  write(*,*)(aaa(i,j),j=1,lda)
!end do

 !write(*,*) 'test2',lda,aaa,lda,w,work,lwork,rwork,info
 call zheev('V','U',lda,aaa,lda,w,work,lwork,rwork,info)

! call heev(aaa,w,'V','U',info)
deallocate(work)

if(info.ne.0)then
write(*,*)'Error during matrix diagonalization. info=',info
write(*,*)'Calculation aborted'
goto 22
endif

 !call cchop(aaa,lda,lda)  !suppression des chop
! call cchop(w,lda,1)

do i=1,lda
!eval(i)=dcmplx(chop(w(i)),0.d0)
	eval(i)=dcmplx(w(i),0.d0)
	do j=1,lda
		evec(i,j)=aaa(i,j)
	enddo
enddo

goto 22
22 continue
end subroutine Diagonalize

!****************************************************************************************
!****************************************************************************************
subroutine cchop(mat,r,s)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
integer(kind=int_4), intent(in):: r,s
complex(kind=comp_16), intent(inout), dimension(r,s) :: mat
integer(kind=int_4):: i,j
real(kind=real_8), allocatable, dimension(:,:)::RE,IM
allocate(RE(r,s))
allocate(IM(r,s))

!$OMP PARALLEL DO PRIVATE(j)
do i=1,r
do j=1,s
RE(i,j)=dreal(mat(i,j))
IM(i,j)=dimag(mat(i,j))
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(j)
do i=1,r
do j=1,s
if(dabs(RE(i,j)).lt.1.d-15)then
RE(i,j)=0.d0
endif
if(dabs(IM(i,j)).lt.1.d-15)then
IM(i,j)=0.d0
endif
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(j)
do i=1,r
do j=1,s
mat(i,j)=dcmplx(RE(i,j),IM(i,j))
enddo
enddo
!$OMP END PARALLEL DO

end subroutine cchop

!****************************************************************************************
!****************************************************************************************
function identiMAT(lda)
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),intent(in)::lda
real(kind=real_8),dimension(lda,lda)::identiMAT
integer::i

identiMAT=0.d0
do i=1,lda
	identiMAT(i,i)=1.d0
end do
end function identiMAT

!****************************************************************************************
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
function ZidentiMAT(lda)
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),intent(in)::lda
complex(kind=comp_16),dimension(lda,lda)::ZidentiMAT
integer(kind=int_4)::i

ZidentiMAT=0.d0
do i=1,lda
	ZidentiMAT(i,i)=dcmplx(1.d0,0.d0)
end do
end function ZidentiMAT

!****************************************************************************************
!****************************************************************************************



subroutine MATinverse(a,ai,lda)
! computes the inverse "ai" of a complex hermitian matrix "a", if non singular.
! Note: the hermitian character of "a" is required to use the unitary property of its 
! eigenvectors (otherwise the computation would yield wrong results)
!****************************************************************************************
!****************************************************************************************

use basics

implicit none
integer(kind=int_4), intent(inout) :: lda
complex(kind=comp_16), intent(inout), dimension(lda,lda) :: a
complex(kind=comp_16), intent(out), dimension(lda,lda) :: ai

complex(kind=comp_16), allocatable, dimension(:,:) :: evec
complex(kind=comp_16), allocatable, dimension(:) :: eval
complex(kind=comp_16), allocatable, dimension(:,:) :: AIdiag,tmp
integer(kind=int_4) :: i,j
logical :: singul

allocate(evec(lda,lda))
allocate(eval(lda))

 call Diagonalize(a,lda,evec,eval)

!check if a is singular...
singul=.false.

do i=1,lda
if(cdabs(eval(i)).lt.1.d-15)then
singul=.true.
endif
enddo

if(singul.eq..true.)then
write(*,*)'Error during matrix inversion, singular matrix...'
write(*,*)'Calculation Aborted'
ai(:,:) = a(:,:)
goto 22
else

allocate(AIdiag(lda,lda))

do i=1,lda
do j=1,lda
AIdiag(i,j)=dcmplx(0.d0,0.d0)
enddo
enddo



do i=1,lda
AIdiag(i,i)=dcmplx(1.d0/dreal(eval(i)),0.d0)
enddo

allocate(tmp(lda,lda))
!ai=MatmulP(evec,MatmulP(AIdiag,dconjg(Transpose(evec)),lda,lda,lda),lda,lda,lda)

tmp = MatmulPAdj(AIdiag,evec,lda,lda,lda)
ai = MatmulP2(evec,tmp,lda,lda,lda)

deallocate(tmp)
deallocate(AIdiag)
endif
22 continue
 call cchop(ai,lda,lda)

end subroutine MATinverse

!****************************************************************************************
!****************************************************************************************
recursive function fact(n,x) result(res)
!Calcul la multifactorielle de n. Si x=1->n!, si x=2->n!!, ...
!warning: rapidly generate big numbers, be sure than int_4 and real_8 correctly set with 
!         expected result.
!warning: unless result is an integer, the output value is a real.
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4), intent(in) :: n
integer(kind=int_4), optional, intent(in) :: x
real(kind=real_8) :: res
integer(int_4)::x_=1
if (present(x))   x_=x

if (n<=0) then
	res=1
else
	res=n*fact(n-x_,x_)
endif
end function fact

!!!!!!!!!!!!!!!!!!!!
Function kronecker(ii,jj)
implicit none
integer::ii,jj
real(kind=real_8)::kronecker
if (ii==jj) then
	kronecker=1.d0
else
	kronecker=0.d0
end if
End function kronecker

end Module math
