Module routines
! Contains all general subroutines needed in the TDMCSCF program
use Basics
use MKL_DFTI
use OMP_LIB
use blas95
use f95_precision
use lapack95
use math
use LiH

contains

!****************************************************************************************
!****************************************************************************************
subroutine aointTOmointMP(hao,zao,hmo,zmo,lcao,nao,nmo)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
integer(kind=int_4), intent(in) :: nao,nmo
real(kind=real_8), intent(in), dimension(nao,nao) :: hao,zao
real(kind=real_8), intent(in), dimension(nao,nmo) :: lcao
real(kind=real_8), intent(out), dimension(nmo,nmo) :: hmo,zmo
real(kind=real_8) :: tmp
integer(kind=int_4) :: i,j,k,l,r,s,t,u

write(*,*)'Integral transformations: AOs to MOs'

!$OMP PARALLEL DO SHARED(zmo,lcao,hao,zao) PRIVATE(NTHREADS, i,j,r,s,tmp)
do i=1,nmo
 do j=1,nmo

 tmp=0.d0
  do r=1,nao
   do s=1,nao
    tmp=tmp+lcao(r,i)*lcao(s,j)*hao(r,s)
   enddo
  enddo

  hmo(i,j)=chop(tmp)

  tmp=0.d0
  do r=1,nao
   do s=1,nao
    tmp=tmp+lcao(r,i)*lcao(s,j)*zao(r,s)
   enddo
  enddo

  zmo(i,j)=chop(tmp)

 enddo
enddo
!$OMP END PARALLEL DO

end subroutine aointTOmointMP

!****************************************************************************************
!****************************************************************************************
subroutine H1e(H1,Ers,hmo,nmo,ncsf)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(in) :: nmo,ncsf
real(kind=real_8), intent(out), dimension(ncsf,ncsf) :: H1
real(kind=real_8), intent(in), dimension(nmo,nmo,ncsf,ncsf) :: Ers
real(kind=real_8), intent(in), dimension(nmo,nmo) :: hmo
integer(kind=int_4) :: r,s

H1(:,:)=0.d0

do r=1,nmo
do s=1,nmo

H1(:,:)=H1(:,:)+(hmo(r,s))*Ers(r,s,:,:)

enddo
enddo
 call rchop(H1,ncsf,ncsf)

end subroutine H1e

!****************************************************************************************
!****************************************************************************************
subroutine H1ett(H1,Ers,hmo,nmo,ncsf)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(in) :: nmo,ncsf
complex(kind=comp_16), intent(out), dimension(ncsf,ncsf) :: H1
real(kind=real_8), intent(in), dimension(nmo,nmo,ncsf,ncsf) :: Ers
complex(kind=comp_16), intent(in), dimension(nmo,nmo) :: hmo
integer(kind=int_4) :: r,s

H1(:,:)=0.d0

do r=1,nmo
do s=1,nmo

H1(:,:)=H1(:,:)+(hmo(r,s))*Ers(r,s,:,:)

enddo
enddo
 call cchop(H1,ncsf,ncsf)

end subroutine H1ett


!****************************************************************************************
!****************************************************************************************
subroutine H2e(H2,erstu,g,nmo,ncsf)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(in) :: nmo,ncsf
real(kind=real_8), intent(out), dimension(ncsf,ncsf) :: H2
real(kind=real_8), intent(in), dimension(nmo,nmo,nmo,nmo,ncsf,ncsf) :: erstu
real(kind=real_8), intent(in), dimension(nmo,nmo,nmo,nmo) :: g
integer(kind=int_4) :: r,s,t,u

H2(:,:)=0.d0

do r=1,nmo
do s=1,nmo
do t=1,nmo
do u=1,nmo
H2(:,:)=H2(:,:)+0.5d0*g(r,s,t,u)*erstu(r,s,t,u,:,:)
enddo
enddo
enddo
enddo
end subroutine H2e

!****************************************************************************************
!****************************************************************************************
subroutine H1eT(H1,hmo,zmo,Et,nmo)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(in) :: nmo
complex(kind=comp_16), intent(out), dimension(nmo,nmo) :: H1
complex(kind=comp_16), intent(in), dimension(nmo,nmo) :: hmo,zmo
real(kind=real_8), intent(in) :: Et
integer(kind=int_4) :: r,s

H1(:,:)=dcmplx(0.d0,0.d0)

do r=1,nmo
do s=1,nmo

H1(r,s)= hmo(r,s)+zmo(r,s)*dcmplx(Et,0.d0)

enddo
enddo
 call cchop(H1,nmo,nmo)

end subroutine H1eT

!****************************************************************************************
!****************************************************************************************
subroutine H2eT(H2,erstu,g,nmo,ncsf)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(in) :: nmo,ncsf
complex(kind=comp_16), intent(out), dimension(ncsf,ncsf) :: H2
real(kind=real_8), intent(in), dimension(nmo,nmo,nmo,nmo,ncsf,ncsf) :: erstu
complex(kind=comp_16), intent(in), dimension(nmo,nmo,nmo,nmo) :: g
integer(kind=int_4) :: r,s,t,u

H2(:,:)=dcmplx(0.d0,0.d0)

do r=1,nmo
do s=1,nmo
do t=1,nmo
do u=1,nmo
H2(:,:)=H2(:,:)+dcmplx(0.5d0,0.d0)*g(r,s,t,u)*dcmplx(erstu(r,s,t,u,:,:),0.d0)
enddo
enddo
enddo
enddo

end subroutine H2eT

!****************************************************************************************
!****************************************************************************************
subroutine H1eTQP(HQQ,HQP,HPQ,HPP,hmo,zmo,Et,nmo,nq,np,dt)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(in) :: nmo,nq,np
complex(kind=comp_16), intent(out), dimension(nq,nq) :: HQQ
complex(kind=comp_16), intent(out), dimension(nq,np) :: HQP
complex(kind=comp_16), intent(out), dimension(np,nq) :: HPQ
complex(kind=comp_16), intent(out), dimension(np,np) :: HPP
complex(kind=comp_16), intent(in), dimension(nmo,nmo) :: hmo,zmo
real(kind=real_8), intent(in) :: Et,dt
integer(kind=int_4) :: r,s

HQQ(:,:)=dcmplx(0.d0,0.d0)
HQP(:,:)=dcmplx(0.d0,0.d0)
HPQ(:,:)=dcmplx(0.d0,0.d0)
HPP(:,:)=dcmplx(0.d0,0.d0)

!$OMP PARALLEL DO
do r=1,nq

	do s=1,nq
	HQQ(r,s)= (hmo(r,s)*dcmplx(dt,0.d0)+zmo(r,s)*dcmplx(Et,0.d0))
	enddo

	do s=1,np
	HQP(r,s)= (hmo(r,s+nq)*dcmplx(dt,0.d0)+zmo(r,s+nq)*dcmplx(Et,0.d0))
	enddo

enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do r=1,np

	do s=1,nq
	HPQ(r,s)= (hmo(r+nq,s)*dcmplx(dt,0.d0)+zmo(r+nq,s)*dcmplx(Et,0.d0))
	enddo

	do s=1,np
	HPP(r,s)= (hmo(r+nq,s+nq)*dcmplx(dt,0.d0)+zmo(r+nq,s+nq)*dcmplx(Et,0.d0))
	enddo

enddo
!$OMP END PARALLEL DO

 call cchop(HQQ,nq,nq)
 call cchop(HQP,nq,np)
 call cchop(HPQ,np,nq)
 call cchop(HPP,np,np)

end subroutine H1eTQP

!****************************************************************************************
!****************************************************************************************
subroutine MatConst(MM1,MM2,MM3,MM4,HPQ,np,nq)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(np,nq) :: HPQ
complex(kind=comp_16), intent(out), dimension(nq,nq) :: MM1
complex(kind=comp_16), intent(out), dimension(nq,np) :: MM2
complex(kind=comp_16), intent(out), dimension(np,np) :: MM3
complex(kind=comp_16), intent(out), dimension(np,nq) :: MM4
complex(kind=comp_16), allocatable, dimension(:,:) :: tmp,XX
complex(kind=comp_16), allocatable, dimension(:,:) :: WW
complex(kind=comp_16), allocatable, dimension(:,:) :: one
integer(kind=int_4) :: i,j
complex(kind=comp_16) :: alpha, beta

allocate(WW(nq,nq))
alpha=dcmplx(0.25d0,0.d0)
beta=dcmplx(0.d0,0.d0)
 call zgemm('C','N',nq,nq,np,alpha,HPQ,np,HPQ,np,beta,WW,nq)

allocate(one(nq,nq))
one(:,:) = dcmplx(0.d0,0.d0)
do i=1,nq
one(i,i) = dcmplx(1.d0,0.d0)
enddo

allocate(tmp(nq,nq))
tmp = One + WW
deallocate(WW)

allocate(XX(nq,nq))
 call MATinverse(tmp,XX,nq)
deallocate(tmp)

MM1(:,:) = dcmplx(2.d0,0.d0)*XX(:,:) - one(:,:)
deallocate(one)

alpha=dcmplx(0.d0,-1.d0)
beta=dcmplx(0.d0,0.d0)
 call zgemm('N','C',nq,np,nq,alpha,XX,nq,HPQ,np,beta,MM2,nq)

MM3(:,:) = dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO
do i=1,np
MM3(i,i) = dcmplx(1.d0,0.d0)
enddo
!$OMP END PARALLEL DO
alpha=dcmplx(0.d0,-0.5d0)
beta=dcmplx(1.d0,0.d0)
 call zgemm('N','N',np,np,nq,alpha,HPQ,np,MM2,nq,beta,MM3,np)

alpha=dcmplx(0.d0,-1.d0)
beta=dcmplx(0.d0,0.d0)
 call zgemm('N','N',np,nq,nq,alpha,HPQ,np,XX,nq,beta,MM4,np)

deallocate(XX)

 call cchop(MM1,nq,nq)
 call cchop(MM2,np,nq)
 call cchop(MM3,np,np)
 call cchop(MM4,nq,np)

106   FORMAT(6E18.6E3)
end subroutine MatConst

!****************************************************************************************
!****************************************************************************************
subroutine readinput(basis_file,omega,E0,xmin,ymin,zmin,delta,theta,phy,pulsetime,nx,ny,nz,nt,nper,pulsed,xmax,ymax,zmax,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,tf,dt,nq,prop,ortho,opw)
!****************************************************************************************
!****************************************************************************************
!! Reads an imput file of the form
!&grille
!basis_file='LiH4e3OM'
!nx= 1
!xmin= 0.d0
!ny= 1
!ymin= 0.d0
!nz= 4
!zmin= -1.d0
!/
!&champ
!E0= 0.25d0
!omega= 0.3d0
!delta=      0.000E+000
!theta= 0.d0
!phy= 0.d0
!pulsed= .false.
!pulsetime=      0.6d0
!nper= 1
!nt= 5000
!/
!&job
!prop= .true.
!ortho= .true.
!opw= .true.
!/
! Then (re-)defines some parameters of the calculation
use basics
implicit none
! Parameters that are read
 character(len=32), intent(out) :: basis_file
real(kind=real_8), intent(out) :: omega,E0,xmin,ymin,zmin,delta,pulsetime,theta,phy
integer(kind=int_4), intent(out) :: nx,ny,nz,nt,nper,nq
logical,intent(out) :: pulsed,prop,ortho,opw
! Parameters that are defined
real(kind=real_8), intent(out) :: xmax,ymax,zmax,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,tf,dt
! Other Parameters
real(kind=real_8) :: nperf
integer(kind=int_4) :: nperi

Namelist/grille/basis_file,nx,xmin,ny,ymin,nz,zmin,nq
Namelist/champ/E0,omega,delta,theta,phy,pulsed,pulsetime,nper,nt
Namelist/job/prop,ortho,opw

open(1,file='input',status='old')
read(1,nml=grille)
read(1,nml=champ)
read(1,nml=job)
 close(1)

basis_file=trim(basis_file)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variables for the external field and time-grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(pulsed.eq..true.)then
!write(*,*)'Duration of the pulse: ',pulsetime,' fs'
! this converts the pulse duration expressed (in the input file) in femtoseconds  to atomic units (as used in the code)
pulsetime=pulsetime/2.41888432650516d-2
!redefintion of nper, if the field is pulsed, 
! that is because nper is an important quantity in the definition of the field (we do not want the user to calculate it each time)
nperf=pulsetime*omega/(2.d0*pi)
!write(*,*)"f=",nperf
nperi=idint(pulsetime*omega/(2.d0*pi))
!write(*,*)"i=",nperi
!write(*,*)"f-i",nint(nperf-dfloat(nperi))
!write(*,*)'nper have been changed such that the pulse duration contains an integer number of optical cycles',nper
nper=Max((2*nperi)+nint(nperf-dfloat(nperi)),2)
!write(*,*)nper
! a factor 2 is needed here because we define the pulse duration in terms of the width at half height,
! which for a squared sine envelope is half of the total interaction time
! we want nper times the optical period to give the final time, while keep writing inputs as much intuitive as possible
endif

tf=dfloat(nper*2)*pi/omega
dt=tf/dfloat(nt)
delta=delta*pi
theta=theta*pi
phy=phy*pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Space grid variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(nx.ge.2)then
!! Usual definition of the grid when at least 2 points are used
dx=2.d0*dabs(xmin)/dfloat(nx)
xmax=xmin+dfloat(nx-1)*dx
dkx=(2.d0*pi)/(dfloat(nx)*dx)
kxmin=-pi/dx
else
!! Definition with a single point, i.e. this dimension is removed in calculation and parameters are given trivial values 
!! (this allows to write the program as it was 3D but use it for 1D calculations without any change)
dx=1.d0
dkx=1.d0
xmin=0.d0
xmax=0.d0
kxmin=0.d0
endif

if(ny.ge.2)then
!! Usual definition of the grid when at least 2 points are used
dy=2.d0*dabs(ymin)/dfloat(ny)
ymax=ymin+dfloat(ny-1)*dy
dky=(2.d0*pi)/(dfloat(ny)*dy)
kymin=-pi/dy
else
!! Definition with a single point, i.e. this dimension is removed in calculation and parameters are given trivial values 
!! (this allows to write the program as it was 3D but use it for 1D calculations without any change)
dy=1.d0
dky=1.d0
ymin=0.d0
ymax=0.d0
kymin=0.d0
endif

if(nz.ge.2)then
!! Usual definition of the grid when at least 2 points are used
dz=2.d0*dabs(zmin)/dfloat(nz)
zmax=zmin+dfloat(nz-1)*dz
dkz=(2.d0*pi)/(dfloat(nz)*dz)
kzmin=-pi/dz
else
!! Definition with a single point, i.e. this dimension is removed in calculation and parameters are given trivial values 
!! (this allows to write the program as it was 3D but use it for 1D calculations without any change)
dz=1.d0
dkz=1.d0
zmin=0.d0
zmax=0.d0
kzmin=0.d0
endif

!! Display the spacial grid parameters on screen (or log file)
write(*,'(A19)') 'Spatial grid size'
write(*,'(A6,F7.2,A4,F7.2,A4,E10.3)')'axe X ',xmin,' to ',xmax,' dx=',dx
write(*,'(A6,F7.2,A4,F7.2,A4,E10.3)')'axe Y ',ymin,' to ',ymax,' dy=',dy
write(*,'(A6,F7.2,A4,F7.2,A4,E10.3)')'axe Z ',zmin,' to ',zmax,' dz=',dz
write(*,*)
write(*,'(A24)')'Grid in momentum space'
write(*,'(A6,F8.2,A4,F8.2,A5,E10.3)')'axe X ',kxmin,' to ',kxmin+dfloat(nx-1)*dkx,' dkx=',dkx
write(*,'(A6,F8.2,A4,F8.2,A5,E10.3)')'axe Y ',kymin,' to ',kymin+dfloat(ny-1)*dky,' dky=',dky
write(*,'(A6,F8.2,A4,F8.2,A5,E10.3)')'axe Z ',kzmin,' to ',kzmin+dfloat(nz-1)*dkz,' dkz=',dkz
write(*,*)
!! Display the field and time grid parameters on screen (or log file)
write(*,'(A22)')'Time-grid Parameters'
write(*,'(A3,F8.2,A4,F8.4)')'tf=',tf,' dt=',dt
write(*,*)
if(pulsed.eq..true.)then
write(*,'(A41)')'External field of sine squared waveform'
write(*,'(A36,F5.2)')'Approximate (input) duration in fs :',pulsetime*2.41888432650516d-2
write(*,'(A47,F5.2)')'Actual duration (after rounding if any) in fs :',0.5*tf*2.41888432650516d-2
else
write(*,'(A16)')'External field'
endif
write(*,'(A3,E10.3,A7,F6.2)')'E0=',E0,' omega=',omega
if(pulsed.eq..true.)then
write(*,'(A27, I2)')'Number of optical period : ',nper/2
else
write(*,'(A27, I2)')'Number of optical period : ',nper
endif
write(*,'(A53,F5.2)')'Phase with respect to envelope (in units of Pi rad) :',idint(delta/pi)
write(*,*)
write(*,*)'Orientation:'
write(*,'(A7,A7,A7)')' E_x ',' E_y ',' E_z '
write(*,'(A2,F5.2,A2,F5.2,A2,F5.2)')' ',EX(1.d0,theta,phy), ' ', EY(1.d0,theta,phy), ' ', EZ(1.d0,theta)
end subroutine readinput


!****************************************************************************************
!****************************************************************************************
Function EX(EE,theta,phy)
!****************************************************************************************
!****************************************************************************************

Use Basics
Implicit None
real(kind=real_8) :: EX
real(kind=real_8), intent(in):: EE,theta,phy
EX = EE*dsin(theta)*dcos(phy)

End Function
!****************************************************************************************
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
Function EY(EE,theta,phy)
!****************************************************************************************
!****************************************************************************************

Use Basics
Implicit None
real(kind=real_8) :: EY
real(kind=real_8), intent(in):: EE,theta,phy
EY = EE*dsin(theta)*dsin(phy)

End Function
!****************************************************************************************
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
Function EZ(EE,theta)
!****************************************************************************************
!****************************************************************************************

Use Basics
Implicit None
real(kind=real_8) :: EZ
real(kind=real_8), intent(in):: EE,theta
EZ = EE*dcos(theta)

End Function
!****************************************************************************************
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
Function Env(omega,t,pulsetime,pulsed,nper)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
real(kind=real_8) :: Env
Real(kind=real_8), intent(in) :: omega,t,pulsetime
integer(kind=int_4), intent(inout) :: nper
real(kind=real_8) :: sigma,tmax,FWHM
logical,intent(in) :: pulsed

FWHM=pulsetime
! tmax here denotes the time associated to maximum field amplitude
tmax=dfloat(nper)*Pi/omega
sigma=FWHM/(2.d0*dsqrt(dlog(2.d0)))
if(pulsed.eq..true.)then
	if(t.lt.2.d0*tmax)then

	Env=dsin(Pi*t/tmax/2.d0)**2.d0
	else
	Env=0.d0
	endif

else
Env=1.d0
endif


end function

!****************************************************************************************
!****************************************************************************************
Function EE(E0,omega,t,delta,pulsetime,pulsed,nper)
!****************************************************************************************
!****************************************************************************************

Use Basics

Implicit None
Real(kind=real_8) :: EE
Real(kind=real_8), intent(in) :: E0,t,omega,delta,pulsetime
integer(kind=int_4), intent(inout) :: nper
logical,intent(in) :: pulsed

if(pulsed.eq..true.)then
EE=Env(omega,t,pulsetime,pulsed,nper)*E0*sin(omega*t+delta)
else
if(t.le.dfloat(nper)*2.d0*Pi/omega)then
EE=Env(omega,t,pulsetime,pulsed,nper)*E0*sin(omega*t+delta)
else
EE=0.d0
endif
endif

End Function

!****************************************************************************************
!****************************************************************************************
Function AAA(E0,omega,tn,delta,pulsetime,pulsed,nper)
!****************************************************************************************
!****************************************************************************************
! Attention this function is only valid for an integer number of optical cycles
! This analytical integration also only stands for more than one cycle pulses,
! (that is because of the one over n-1 factor in the expressions)
! a generalization could be done later if one is interessed in such pulses
Use Basics

Implicit None
Real(kind=real_8) :: AAA
Real(kind=real_8), intent(in) :: E0,omega,delta,pulsetime,tn
integer(kind=int_4), intent(in) :: nper
logical,intent(in) :: pulsed
Integer(kind=int_4) :: n
Real(kind=real_8) :: tf

n=nper
tf = dfloat(2*nper)*Pi/omega


if(pulsed.eq..true.)then

AAA = E0*(-1.d0*dcos(delta)*dcos(omega*tn)/(2.d0*omega) &
& + dsin(delta)*dsin(omega*tn)/(2.d0*omega) &
& + (dfloat(n)*dcos(delta)/(4.d0*omega))*(dcos(dfloat(n+1)/dfloat(n)*omega*tn)/dfloat(n+1)&
& + dcos(dfloat(n-1)/dfloat(n)*omega*tn)/dfloat(n-1))&
& - (dfloat(n)*dsin(delta)/(4.d0*omega))*(dsin(dfloat(n+1)/dfloat(n)*omega*tn)/dfloat(n+1)&
& + dsin(dfloat(n-1)/dfloat(n)*omega*tn)/dfloat(n-1)))

else

AAA = -1.d0*E0*dcos(tn*omega)/omega

endif

if(tn.gt.tf)then
AAA = 0.d0
endif

End Function

!****************************************************************************************
!****************************************************************************************
subroutine buildMuPhi(MuPhi,MuChi,LCAO,nx,nao,nmo)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
Integer(kind=int_4), intent(inout) :: nx,nao,nmo
Complex(kind=comp_16), intent(in), dimension(nao,nx) :: MuChi
real(kind=real_8), intent(in), dimension(nao,nmo) :: LCAO
Complex(kind=comp_16), intent(out), dimension(nmo,nx) :: MuPhi
integer(kind=int_4) :: i,j,k,l
Complex(kind=comp_16) :: mutmp

do k=1,nx
do i=1,nmo
	mutmp = dcmplx(0.d0,0.d0)
	do j=1,nao
	mutmp = mutmp + dcmplx(LCAO(j,i),0.d0)*MuChi(j,k)
	enddo

MuPhi(i,k)=mutmp

enddo
enddo

end subroutine buildMuPhi

!****************************************************************************************
!****************************************************************************************
subroutine MatUSE(LL,MM1,MM2,MM3,MM4,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq) :: LL

complex(kind=comp_16), intent(inout), dimension(nq,nq) :: MM1
complex(kind=comp_16), intent(inout), dimension(nq,np) :: MM2
complex(kind=comp_16), intent(inout), dimension(np,np) :: MM3
complex(kind=comp_16), intent(inout), dimension(np,nq) :: MM4

complex(kind=comp_16), allocatable, dimension(:,:) :: tmp,tmp2,LLq,LLp
complex(kind=comp_16) :: alpha, beta
integer(kind=int_4) :: i,j
alpha=dcmplx(1.d0,0.d0)
beta=dcmplx(0.d0,0.d0)

! Copy all MM_I blocks in a work array
allocate(tmp(nq,nq))
allocate(tmp2(nq,nq))
allocate(LLq(nq,nq))
allocate(LLP(np,nq))

!$OMP PARALLEL DO
do i=1,nq
do j=1,nq
LLq(i,j)=LL(i,j)
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do i=1,np
do j=1,nq
LLp(i,j)=LL(i+nq,j)
enddo
enddo
!$OMP END PARALLEL DO

 call zgemm('N','N',nq,nq,nq,alpha,MM1,nq,LLq,nq,beta,tmp,nq)
 call zgemm('N','N',nq,nq,np,alpha,MM2,nq,LLp,np,beta,tmp2,nq)

!$OMP PARALLEL DO
do i=1,nq
do j=1,nq
LL(i,j)=tmp(i,j)+tmp2(i,j)
enddo
enddo
!$OMP END PARALLEL DO

deallocate(tmp)
deallocate(tmp2)
allocate(tmp(np,nq))
allocate(tmp2(np,nq))

 call zgemm('N','N',np,nq,np,alpha,MM3,np,LLp,np,beta,tmp2,np)
 call zgemm('N','N',np,nq,nq,alpha,MM4,np,LLq,nq,beta,tmp,np)

!$OMP PARALLEL DO
do i=1,np
do j=1,nq
LL(i+nq,j)=tmp(i,j)+tmp2(i,j)
enddo
enddo
!$OMP END PARALLEL DO


end subroutine MatUSE

!****************************************************************************************
!****************************************************************************************
subroutine U1QQl(LL,propag,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq) :: LL
complex(kind=comp_16), intent(inout), dimension(nq,nq) :: propag
complex(kind=comp_16), allocatable, dimension(:,:) :: tmp
integer(kind=int_4) :: i
complex(kind=comp_16) :: alpha, beta
alpha=dcmplx(1.d0,0.d0)
beta=dcmplx(0.d0,0.d0)

allocate(tmp(nq,nq))

tmp(1:nq,1:nq) = matmul(propag(1:nq,1:nq),LL(1:nq,1:nq))
!tmp(1:nq,1:nq) = matmul(LL(1:nq,1:nq),propag(1:nq,1:nq))
LL(1:nq,1:nq)=tmp(1:nq,1:nq)

deallocate(tmp)

end subroutine U1QQl

!****************************************************************************************
!****************************************************************************************
subroutine U1QQ(LL,propag,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq+np) :: LL

complex(kind=comp_16), intent(inout), dimension(nq,nq) :: propag

complex(kind=comp_16), allocatable, dimension(:,:) :: tmp,tmp2
integer(kind=int_4) :: i,j,k,n
complex(kind=comp_16) :: partsum
n=nq+np

allocate(tmp(nq+np,nq))
allocate(tmp2(nq+np,nq))

!$OMP PARALLEL DO
do i=1,nq+np
do j=1,nq
tmp(i,j) = LL(i,j)
enddo
enddo
!$OMP END PARALLEL DO

tmp2 = matmulP2(tmp,propag,n,nq,nq)

deallocate(tmp)

!$OMP PARALLEL DO
do i=1,nq+np
do j=1,nq
LL(i,j) = tmp2(i,j)
enddo
enddo
!$OMP END PARALLEL DO

deallocate(tmp2)

end subroutine U1QQ

!****************************************************************************************
!****************************************************************************************
subroutine U1PP(LL,propag,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq+np) :: LL

complex(kind=comp_16), intent(inout), dimension(np,np) :: propag

complex(kind=comp_16), allocatable, dimension(:,:) :: tmp,tmp2
integer(kind=int_4) :: i,j,k
complex(kind=comp_16) :: partsum

allocate(tmp(np+nq,np))
allocate(tmp2(np+nq,np))

!$OMP PARALLEL DO

do i=1,np+nq
do j=1,np

tmp(i,j) = LL(i,j+nq)

enddo
enddo
!$OMP END PARALLEL DO

tmp2 = matmulP2(tmp,propag,nq+np,np,np)

deallocate(tmp)

!$OMP PARALLEL DO

do i=1,np+nq
do j=1,np

LL(i,j+nq) = tmp2(i,j)

enddo
enddo
!$OMP END PARALLEL DO

deallocate(tmp2)

end subroutine U1PP

!****************************************************************************************
!****************************************************************************************
subroutine U1(LL,propag,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics
Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq+np) :: LL

complex(kind=comp_16), intent(inout), dimension(np+nq,np+nq) :: propag

complex(kind=comp_16), allocatable, dimension(:,:) :: tmp
integer(kind=int_4) :: i,j,k
complex(kind=comp_16) :: partsum

allocate(tmp(np+nq,np+nq))

!$OMP PARALLEL DO
do i=1,np+nq
do j=1,np+nq
tmp(i,j) = LL(i,j)
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO REDUCTION(+:partsum)
do i=1,np+nq
do j=1,np+nq

partsum = dcmplx(0.d0,0.d0)

	do k=1,np+nq
	partsum = partsum + tmp(i,k)*propag(k,j)
	enddo

LL(i,j) = partsum

enddo
enddo
!$OMP END PARALLEL DO

deallocate(tmp)

end subroutine U1

!****************************************************************************************
! Computes the result of LL^dagger . (h + z E(t)) . LL
! The term r.E is not used at this stage, thus Et is not needed here
subroutine transmo1(hmot,zmot,LL,nbase,tothmo,totzmo)
!****************************************************************************************
!****************************************************************************************

use basics
use math
implicit none
integer(kind=int_4), intent(inout) :: nbase
complex(kind=comp_16), intent(in), dimension(nbase,nbase) :: LL
complex(kind=comp_16), intent(out), dimension(nbase,nbase) :: hmot,zmot
complex(kind=comp_16), intent(in), dimension(nbase,nbase) :: tothmo,totzmo
complex(kind=comp_16) :: alpha,beta
complex(kind=comp_16), allocatable, dimension(:,:) :: tmp,tmp2
integer(kind=int_4) :: n
!allocate(tmp(nbase,nbase))
!tmp(:,:) = LL(:,:)


allocate(tmp2(nbase,nbase))
n = nbase
alpha=dcmplx(1.d0,0.d0)
beta=dcmplx(0.d0,0.d0)

 call zgemm('N','N',n,n,n,alpha,tothmo,n,LL,n,beta,tmp2,n)
 call zgemm('C','N',n,n,n,alpha,LL,n,tmp2,n,beta,hmot,n)
! call zgemm('N','C',n,n,n,alpha,tothmo,n,LL,n,beta,tmp2,n)
! call zgemm('N','N',n,n,n,alpha,LL,n,tmp2,n,beta,hmot,n)
 call zgemm('N','N',n,n,n,alpha,totzmo,n,LL,n,beta,tmp2,n)
 call zgemm('C','N',n,n,n,alpha,LL,n,tmp2,n,beta,zmot,n)
! call zgemm('N','C',n,n,n,alpha,totzmo,n,LL,n,beta,tmp2,n)
! call zgemm('N','N',n,n,n,alpha,LL,n,tmp2,n,beta,zmot,n)

!deallocate(tmp)
deallocate(tmp2)

end subroutine transmo1

!****************************************************************************************
!****************************************************************************************
subroutine getWP(LL,My_Desc_Handle3,nq,np,dkx,dky,dkz,nx,ny,nz,lin,WP)
!****************************************************************************************
!****************************************************************************************
use basics
use math
Implicit none
integer(kind=int_4), intent(in) :: lin
integer(kind=int_4), intent(inout) :: nq,np,nx,ny,nz
complex(kind=comp_16), intent(inout), dimension(nq+np,nq) :: LL
real(kind=real_8), intent(in) :: dkx,dky,dkz
type(DFTI_DESCRIPTOR), Pointer :: My_Desc_Handle3
complex(kind=comp_16), dimension(np) :: WP
real(kind=real_8) :: sfact
complex(kind=comp_16), allocatable, dimension(:) :: tmp
integer(kind=int_4) :: i,j,k,isgn

allocate(tmp(np))

!$OMP PARALLEL DO
do i=1,np
tmp(i) = LL(i+nq,lin)
enddo
!$OMP END PARALLEL DO

! Fourier Transform k -> r
isgn=1
sfact=(dkx*(2.d0*pi)**-0.5d0)*(dky*(2.d0*pi)**-0.5d0)*(dkz*(2.d0*pi)**-0.5d0)
 call FFT_calc(tmp,np,sfact,isgn,My_Desc_Handle3)
 call sort_GEN(tmp,nx,ny,nz)
WP(:) = tmp(:)

deallocate(tmp)

end subroutine getWP

!****************************************************************************************
!****************************************************************************************
subroutine getWP2(GIT,My_Desc_Handle3,np,dkx,dky,dkz,nx,ny,nz,WP)
!****************************************************************************************
!****************************************************************************************
use basics
use math
Implicit none
integer(kind=int_4), intent(inout) :: np,nx,ny,nz
complex(kind=comp_16), intent(inout), dimension(np) :: GIT
real(kind=real_8), intent(in) :: dkx,dky,dkz
type(DFTI_DESCRIPTOR), Pointer :: My_Desc_Handle3
complex(kind=comp_16), dimension(np) :: WP
real(kind=real_8) :: sfact
complex(kind=comp_16), allocatable, dimension(:) :: tmp
integer(kind=int_4) :: i,j,k,isgn

allocate(tmp(np))

!$OMP PARALLEL DO
do i=1,np
tmp(i) = GIT(i)
enddo
!$OMP END PARALLEL DO

! Fourier Transform k -> r
isgn=1
sfact=(dkx*(2.d0*pi)**-0.5d0)*(dky*(2.d0*pi)**-0.5d0)*(dkz*(2.d0*pi)**-0.5d0)
 call FFT_calc(tmp,np,sfact,isgn,My_Desc_Handle3)
 call sort_GEN(tmp,nx,ny,nz)
WP(:) = tmp(:)

deallocate(tmp)

end subroutine getWP2

!****************************************************************************************
!****************************************************************************************
subroutine U1QQ2(LL,propag,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq+np) :: LL

complex(kind=comp_16), intent(inout), dimension(nq,nq) :: propag

complex(kind=comp_16), allocatable, dimension(:,:) :: tmp,tmp2
integer(kind=int_4) :: i,j,k,n
complex(kind=comp_16) :: partsum
n=nq+np
allocate(tmp(nq,nq+np))
allocate(tmp2(nq,nq+np))

!$OMP PARALLEL DO
do i=1,nq
do j=1,nq+np
tmp(i,j) = LL(i,j)
enddo
enddo
!$OMP END PARALLEL DO

!tmp2 = matmulP(tmp,propag,nq+np,nq,nq)
tmp2 = matmulP2(propag,tmp,nq,nq,n)
deallocate(tmp)

!$OMP PARALLEL DO
do i=1,nq
do j=1,nq+np
LL(i,j) = tmp2(i,j)
enddo
enddo
!$OMP END PARALLEL DO

deallocate(tmp2)

end subroutine U1QQ2

!****************************************************************************************
!****************************************************************************************
subroutine U1PP2(LL,propag,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq+np) :: LL

complex(kind=comp_16), intent(inout), dimension(np,np) :: propag

complex(kind=comp_16), allocatable, dimension(:,:) :: tmp,tmp2
integer(kind=int_4) :: i,j,k
complex(kind=comp_16) :: partsum

allocate(tmp(np,np+nq))
allocate(tmp2(np,np+nq))

!$OMP PARALLEL DO
do i=1,np
do j=1,np+nq
tmp(i,j) = LL(i+nq,j)
enddo
enddo
!$OMP END PARALLEL DO

!tmp2 = matmulP2(propag,tmp,nq+np,np,np)
tmp2 = matmulP2(propag,tmp,np,np,np+nq)

deallocate(tmp)

!$OMP PARALLEL DO
do i=1,np
do j=1,np+nq
LL(i+nq,j) = tmp2(i,j)
enddo
enddo
!$OMP END PARALLEL DO

deallocate(tmp2)

end subroutine U1PP2

!****************************************************************************************
!****************************************************************************************
subroutine rDotE2(LL,At,My_Desc_Handle1,nq,np,dx,dkx,xmin)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq+np) :: LL
real(kind=real_8), intent(in) :: At,dx,dkx,xmin
type(DFTI_DESCRIPTOR), Pointer :: My_Desc_Handle1

real(kind=real_8) :: sfact
complex(kind=comp_16), allocatable, dimension(:) :: tmp,tmp2
integer(kind=int_4) :: i,j,k,lin,isgn
complex(kind=comp_16) :: partsum

do lin=1,nq+np

allocate(tmp(np))
allocate(tmp2(np))

!$OMP PARALLEL DO
do i=1,np
tmp(i) = LL(i+nq,lin)
enddo
!$OMP END PARALLEL DO

! Fourier Transform k -> r
isgn=-1
sfact=(dkx*(2.d0*pi)**-0.5)
 call FFT_calc(tmp,np,sfact,isgn,My_Desc_Handle1)

! Effect of Linear Electric Potential (currently under SFA)
!$OMP PARALLEL DO
do i=1,np/2
 tmp2(i)=tmp(i)*cdexp(dcmplx(0.d0,-1.d0*At*(xmin+dfloat(i-1+np/2)*dx)))
 tmp2(i+np/2)=tmp(i+np/2)*cdexp(dcmplx(0.d0,-1.d0*At*(xmin+dfloat(i-1)*dx)))
enddo
!$OMP END PARALLEL DO

deallocate(tmp)

! Fourier Transform r -> k
isgn=1
sfact=(dx*(2.d0*pi)**-0.5)

 call FFT_calc(tmp2,np,sfact,isgn,My_Desc_Handle1)

!$OMP PARALLEL DO
do i=1,np
LL(i+nq,lin) = tmp2(i)
enddo
!$OMP END PARALLEL DO

deallocate(tmp2)

enddo

end subroutine rDotE2

!****************************************************************************************
!****************************************************************************************
subroutine MatUSE2(LL,MM1,MM2,MM3,MM4,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(inout), dimension(nq+np,nq+np) :: LL

complex(kind=comp_16), intent(inout), dimension(nq,nq) :: MM1
complex(kind=comp_16), intent(inout), dimension(nq,np) :: MM2
complex(kind=comp_16), intent(inout), dimension(np,np) :: MM3
complex(kind=comp_16), intent(inout), dimension(np,nq) :: MM4

complex(kind=comp_16), allocatable, dimension(:,:) :: tmp,M12,M34
integer(kind=int_4) :: i,j,k,n
complex(kind=comp_16) :: partsum
n=nq+np
allocate(tmp(nq+np,nq+np))
allocate(M12(nq+np,nq+np))

M12(:,:)=dcmplx(0.d0,0.d0)

! 1) Store the old values in memory
!$OMP PARALLEL DO
do i = 1, nq+np
do j = 1, nq+np

tmp(i,j) = LL(i,j)

enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do i=1,nq
do j=1,nq
M12(i,j)=MM1(i,j)
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do i=1,nq
do j=1,np
M12(i,j+nq)=MM2(i,j)
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do i=1,np
do j=1,nq
M12(i+nq,j)=MM4(i,j)
enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
do i=1,np
do j=1,np
M12(i+nq,j+nq)=MM3(i,j)
enddo
enddo
!$OMP END PARALLEL DO

!LL = matmulP(tmp,M12,nq+np,nq+np,nq+np)
LL = MatmulP2(M12,tmp,n,n,n)

deallocate(tmp)
deallocate(M12)
 call cchop(LL,nq+np,nq+np)

end subroutine MatUSE2

!****************************************************************************************
!****************************************************************************************
subroutine checkcut(psi,nx,check)
!****************************************************************************************
!****************************************************************************************

Use Basics

Implicit None
integer(kind=int_4), intent(in) :: nx
complex(kind=comp_16), intent(in), dimension(nx) :: psi
logical, intent(out) :: check
real(kind=real_8) :: seuil
integer(kind=int_4) :: x

 check=.false.
seuil=5.d-7

x=nx/4
if(cdabs(psi(x)).gt.seuil)then
 check=.true.
endif

x=3*nx/4
if(cdabs(psi(x)).gt.seuil)then
 check=.true.
endif

end subroutine checkcut

!****************************************************************************************
!****************************************************************************************
Function Mask(ni,N)
!****************************************************************************************
!****************************************************************************************

Use Basics

Implicit None
complex(kind=comp_16) :: Mask
integer(kind=int_4), intent(in) :: ni,N
integer(kind=int_4) :: ct
real(kind=real_8):: tmp

 ct=N/4
tmp=0.d0
if(ni.le.ct)then
tmp=1.d0
endif
if((ni.gt.ct).and.(ni.le.2*ct))then
tmp=0.5d0*(dcos(dfloat(ni-ct-1)/dfloat(ct-1)*pi)+1)
endif
if((ni.gt.N-2*ct).and.(ni.le.N-ct))then
tmp=0.5d0*(-dcos(dfloat(ni-1-N+2*ct)/dfloat(ct-1)*pi)+1)
endif
if(ni.gt.N-ct)then
tmp=1.d0
endif
Mask=dcmplx(tmp,0.d0)

end Function

!****************************************************************************************
!****************************************************************************************
Subroutine Cut(psi,nx)
!****************************************************************************************
!****************************************************************************************

Use Basics
Implicit None
integer(kind=int_4), intent(in):: nx
complex(kind=comp_16), intent(inout), dimension(nx) :: psi
integer(kind=int_4) :: x

!OMP PARALLEL DO
do x=1,nx
psi(x)=Mask(x,nx)*psi(x)
enddo
!OMP END PARALLEL DO


End Subroutine Cut

!****************************************************************************************
!****************************************************************************************
Subroutine asympta(psi,psia,nx,dx,My_Desc_Handle1)
!****************************************************************************************
!****************************************************************************************

Use Basics

Implicit None
type(DFTI_DESCRIPTOR), Pointer :: My_Desc_Handle1
integer(kind=int_4), intent(in) :: nx
real(kind=real_8), intent(inout) :: dx
complex(kind=comp_16), intent(inout), dimension(nx) :: psi
complex(kind=comp_16), intent(inout), dimension(nx) :: psia
complex(kind=comp_16), allocatable, dimension(:) :: ctmp
integer(kind=int_4) :: x,isgn
real(kind=real_8) :: sfact

! Extraction of the Asymptotic component
allocate(ctmp(nx))
 
!$OMP PARALLEL DO
do x=1,nx
 ctmp(x)=(dcmplx(1.d0,0.d0)-Mask(x,nx))*psi(x)
enddo
!$OMP END PARALLEL DO

! Fourier Transform r -> k
isgn=1
sfact=(dx*(2.d0*pi)**-0.5)

 call FFT_calc(ctmp,nx,sfact,isgn,My_Desc_Handle1)


do x=1,nx
 psia(x)=cdabs(ctmp(x))**2.d0
enddo


deallocate(ctmp)

 call cut(psi,nx)

end subroutine asympta

!****************************************************************************************
!****************************************************************************************
subroutine rDotE(LL,At,My_Desc_Handle1,nq,np,dx,dkx,xmin,LLA,nanalyse)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nq,np,nanalyse
complex(kind=comp_16), intent(inout), dimension(nq+np,nq+np) :: LL
complex(kind=comp_16), intent(inout), dimension(nq,np) :: LLA
real(kind=real_8), intent(inout) :: At,dx,dkx,xmin
type(DFTI_DESCRIPTOR), Pointer :: My_Desc_Handle1
logical :: check
real(kind=real_8) :: sfact
complex(kind=comp_16), allocatable, dimension(:) :: tmp,tmp2,psia
integer(kind=int_4) :: i,j,k,lin,isgn
complex(kind=comp_16) :: partsum

do lin=1,nq+np

allocate(tmp(np))
allocate(tmp2(np))
allocate(psia(np))

!$OMP PARALLEL DO
do i=1,np
tmp(i) = LL(lin,i+nq)
enddo
!$OMP END PARALLEL DO

! Fourier Transform k -> r
isgn=-1
sfact=(dkx*(2.d0*pi)**-0.5)
 call FFT_calc(tmp,np,sfact,isgn,My_Desc_Handle1)

! call cut(tmp,np)

nanalyse=nanalyse+1
!!!!!!!!!$OMP PARALLEL DO
!!!!!!!!do i=1,np
!!!!!!!!psia(i) = LLA(lin,i)
!!!!!!!!enddo
!!!!!!!!!$OMP END PARALLEL DO
! call asympta(tmp,psia,np,dx,My_Desc_Handle1)

!!!!!!!!$OMP PARALLEL DO
!!!!!!!do i=1,np
!!!!!!!LLA(lin,i) = psia(i)
!!!!!!!enddo
!!!!!!!!$OMP END PARALLEL DO

deallocate(psia)

! Effect of Linear Electric Potential (currently under SFA)
!$OMP PARALLEL DO
do i=1,np/2
 tmp2(i)=tmp(i)*cdexp(dcmplx(0.d0,-1.d0*At*(xmin+dfloat(i-1+np/2)*dx)))
 tmp2(i+np/2)=tmp(i+np/2)*cdexp(dcmplx(0.d0,-1.d0*At*(xmin+dfloat(i-1)*dx)))
enddo
!$OMP END PARALLEL DO

deallocate(tmp)

! Fourier Transform r -> k
isgn=1
sfact=(dx*(2.d0*pi)**-0.5)

 call FFT_calc(tmp2,np,sfact,isgn,My_Desc_Handle1)

!$OMP PARALLEL DO
do i=1,np
LL(lin,i+nq) = tmp2(i)
enddo
!$OMP END PARALLEL DO

deallocate(tmp2)

enddo

end subroutine rDotE

!****************************************************************************************
!****************************************************************************************
subroutine storeAll(tothmo,totzmo,MuPhi,hmo,zmo,kxmin,dkx,nbase,nmo,nx)
!****************************************************************************************
!****************************************************************************************

use Basics
Implicit None
Integer(kind=int_4), intent(in) :: nx,nmo,nbase
complex(kind=comp_16), intent(out), dimension(nbase,nbase) :: tothmo,totzmo
real(kind=real_8), intent(in), dimension(nmo,nmo) :: hmo,zmo
complex(kind=comp_16), intent(in), dimension(nmo,nx) :: MuPhi
real(kind=real_8), intent(in) :: kxmin,dkx
Integer(kind=int_4) :: i,j,nq,np

nq=nmo
np=nx

tothmo(:,:) = dcmplx(0.d0,0.d0)
totzmo(:,:) = dcmplx(0.d0,0.d0)
!! QQ Block
do i=1,nmo
do j=1,nmo
tothmo(i,j)=dcmplx(hmo(i,j),0.d0)
totzmo(i,j)=dcmplx(zmo(i,j),0.d0)
enddo
enddo

do i=1,np
!! PP Block
tothmo(i+nq,i+nq)=dcmplx(0.5d0*(kxmin+dfloat(i-1)*dkx)**2.d0,0.d0)
enddo

do i=1,np
do j=1,nq

!! QP, PQ Blocks
!totzmo(i+nq,j)=MuPhi(j,i)
totzmo(i+nq,j)=dconjg(MuPhi(j,i))
totzmo(j,i+nq)=dconjg(totzmo(i+nq,j))

enddo
enddo

end subroutine storeAll

!****************************************************************************************
!****************************************************************************************
subroutine transchannels(LLA,LLAX,dkx,nmo,nx,My_Desc_Handle1)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(inout) :: nmo,nx
!complex(kind=comp_16), intent(inout), dimension(nmo,nx) :: LLA
!complex(kind=comp_16), intent(out), dimension(nmo,nx) :: LLAX
complex(kind=comp_16), intent(inout), dimension(nx,nmo) :: LLA
complex(kind=comp_16), intent(out), dimension(nx,nmo) :: LLAX

real(kind=real_8), intent(inout) :: dkx
type(DFTI_DESCRIPTOR), Pointer :: My_Desc_Handle1
real(kind=real_8) :: sfact
complex(kind=comp_16), allocatable, dimension(:) :: tmp
integer(kind=int_4) :: i,j,k,lin,isgn

allocate(tmp(nx))

do lin=1,nmo
!tmp(:) = LLA(lin,:)
tmp(:) = LLA(:,lin)

! Fourier Transform k -> r
isgn=1
sfact=(dkx*(2.d0*pi)**-0.5)
 call FFT_calc(tmp,nx,sfact,isgn,My_Desc_Handle1)
 call sort_1D(tmp,nx)

!LLAX(lin,:)=tmp(:)
LLAX(:,lin)=tmp(:)

enddo

end subroutine transchannels

!****************************************************************************************
!****************************************************************************************
subroutine Uvolkov(kxmin,dkx,dt,At,Uv,np)
!****************************************************************************************
!****************************************************************************************
use basics

Implicit none
integer(kind=int_4), intent(in) :: np
complex(kind=comp_16), intent(out), dimension(np,np) :: Uv
real(kind=real_8), intent(in) :: dkx,kxmin,dt,At
real(kind=real_8) :: kk
integer(kind=int_4) :: i

! Initialize the matrix
Uv(:,:) = dcmplx(0.d0,0.d0)

do i=1,np
kk=kxmin+dfloat(i-1)*dkx
Uv(i,i) = cdexp(dcmplx(0.d0,-0.5d0*dt*(kk**2.d0)))*cdexp(dcmplx(0.d0,-0.5d0*At*At*dt))*cdexp(dcmplx(0.d0,At*kk*dt))
enddo

end subroutine Uvolkov

!****************************************************************************************
!****************************************************************************************
subroutine Propag1Dx(LLn,rmin,dr,nr,idt,Er,NS,coord,znuc)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nr
complex(kind=comp_16), intent(inout), dimension(nr) :: LLn
real(kind=real_8), intent(in) :: rmin,dr,Er
complex(kind=comp_16), intent(in) :: idt

integer(kind=int_4),intent(in) :: ns
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc

integer(kind=int_4) :: i

!$OMP PARALLEL DO
do i=1,nr/2
LLn(i+nr/2)=LLn(i+nr/2)*cdexp(dcmplx((rmin+dfloat(i-1)*dr)*Er,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,rmin+dfloat(i-1)*dr,0.d0,0.d0),0.d0)*idt)
LLn(i)=LLn(i)*cdexp(dcmplx((rmin+dfloat(i-1+nr/2)*dr)*Er,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,rmin+dfloat(i-1+nr/2)*dr,0.d0,0.d0),0.d0)*idt)
enddo
!$OMP END PARALLEL DO

end subroutine propag1Dx

!****************************************************************************************
!****************************************************************************************
subroutine Propag1Dy(LLn,rmin,dr,nr,idt,Er,NS,coord,znuc)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nr
complex(kind=comp_16), intent(inout), dimension(nr) :: LLn
real(kind=real_8), intent(in) :: rmin,dr,Er
complex(kind=comp_16), intent(in) :: idt

integer(kind=int_4),intent(in) :: ns
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc

integer(kind=int_4) :: i

!$OMP PARALLEL DO
do i=1,nr/2
LLn(i+nr/2)=LLn(i+nr/2)*cdexp(dcmplx((rmin+dfloat(i-1)*dr)*Er,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,0.d0,rmin+dfloat(i-1)*dr,0.d0),0.d0)*idt)
LLn(i)=LLn(i)*cdexp(dcmplx((rmin+dfloat(i-1+nr/2)*dr)*Er,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,0.d0,rmin+dfloat(i-1+nr/2)*dr,0.d0),0.d0)*idt)
enddo
!$OMP END PARALLEL DO

end subroutine propag1Dy

!****************************************************************************************
!****************************************************************************************
subroutine Propag1Dz(LLn,rmin,dr,nr,idt,Er,NS,coord,znuc)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nr
complex(kind=comp_16), intent(inout), dimension(nr) :: LLn
real(kind=real_8), intent(in) :: rmin,dr,Er
complex(kind=comp_16), intent(in) :: idt

integer(kind=int_4),intent(in) :: ns
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc

integer(kind=int_4) :: i

!$OMP PARALLEL DO
do i=1,nr/2
LLn(i+nr/2)=LLn(i+nr/2)*cdexp(dcmplx((rmin+dfloat(i-1)*dr)*Er,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,0.d0,0.d0,rmin+dfloat(i-1)*dr),0.d0)*idt)
LLn(i)=LLn(i)*cdexp(dcmplx((rmin+dfloat(i-1+nr/2)*dr)*Er,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,0.d0,0.d0,rmin+dfloat(i-1+nr/2)*dr),0.d0)*idt)
enddo
!$OMP END PARALLEL DO

end subroutine propag1Dz

!****************************************************************************************
!****************************************************************************************
subroutine Propag2Dxy(LLn,xmin,dx,nx,ymin,dy,ny,idt,EEx,EEy,NS,coord,znuc)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nx,ny
complex(kind=comp_16), intent(inout), dimension(nx*ny) :: LLn
real(kind=real_8), intent(in) :: xmin,dx,EEx,ymin,dy,EEy
complex(kind=comp_16), intent(in) :: idt

integer(kind=int_4),intent(in) :: ns
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc

integer(kind=int_4) :: i,j

!$OMP PARALLEL DO
do i=1,nx/2

	do j=1,ny/2

LLn((j-1)*nx+i)=LLn((j-1)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i+nx/2-1)*dx,ymin+dfloat(j+ny/2-1)*dy,0.d0),0.d0)*idt)

LLn((j-1)*nx+i+nx/2)=LLn((j-1)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i-1)*dx,ymin+dfloat(j+ny/2-1)*dy,0.d0),0.d0)*idt)
LLn((j-1+ny/2)*nx+i)=LLn((j-1+ny/2)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i+nx/2-1)*dx,ymin+dfloat(j-1)*dy,0.d0),0.d0)*idt)

LLn((j-1+ny/2)*nx+i+nx/2)=LLn((j-1+ny/2)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i-1)*dx,ymin+dfloat(j-1)*dy,0.d0),0.d0)*idt)

	enddo

enddo
!$OMP END PARALLEL DO

end subroutine propag2Dxy

!****************************************************************************************
!****************************************************************************************
subroutine Propag2Dyz(LLn,xmin,dx,nx,ymin,dy,ny,idt,EEx,EEy,NS,coord,znuc)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nx,ny
complex(kind=comp_16), intent(inout), dimension(nx*ny) :: LLn
real(kind=real_8), intent(in) :: xmin,dx,EEx,ymin,dy,EEy
complex(kind=comp_16), intent(in) :: idt

integer(kind=int_4),intent(in) :: ns
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc

integer(kind=int_4) :: i,j

!$OMP PARALLEL DO
do i=1,nx/2

	do j=1,ny/2

LLn((j-1)*nx+i)=LLn((j-1)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,0.d0,xmin+dfloat(i+nx/2-1)*dx,ymin+dfloat(j+ny/2-1)*dy),0.d0)*idt)

LLn((j-1)*nx+i+nx/2)=LLn((j-1)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,0.d0,xmin+dfloat(i-1)*dx,ymin+dfloat(j+ny/2-1)*dy),0.d0)*idt)
LLn((j-1+ny/2)*nx+i)=LLn((j-1+ny/2)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,0.d0,xmin+dfloat(i+nx/2-1)*dx,ymin+dfloat(j-1)*dy),0.d0)*idt)

LLn((j-1+ny/2)*nx+i+nx/2)=LLn((j-1+ny/2)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,0.d0,xmin+dfloat(i-1)*dx,ymin+dfloat(j-1)*dy),0.d0)*idt)
	enddo

enddo
!$OMP END PARALLEL DO

end subroutine propag2Dyz

!****************************************************************************************
!****************************************************************************************
subroutine Propag2Dxz(LLn,xmin,dx,nx,ymin,dy,ny,idt,EEx,EEy,NS,coord,znuc)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nx,ny
complex(kind=comp_16), intent(inout), dimension(nx*ny) :: LLn
real(kind=real_8), intent(in) :: xmin,dx,EEx,ymin,dy,EEy
complex(kind=comp_16), intent(in) :: idt

integer(kind=int_4),intent(in) :: ns
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc

integer(kind=int_4) :: i,j

!$OMP PARALLEL DO
do i=1,nx/2

	do j=1,ny/2

LLn((j-1)*nx+i)=LLn((j-1)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i+nx/2-1)*dx,0.d0,ymin+dfloat(j+ny/2-1)*dy),0.d0)*idt)

LLn((j-1)*nx+i+nx/2)=LLn((j-1)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i-1)*dx,0.d0,ymin+dfloat(j+ny/2-1)*dy),0.d0)*idt)
LLn((j-1+ny/2)*nx+i)=LLn((j-1+ny/2)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i+nx/2-1)*dx,0.d0,ymin+dfloat(j-1)*dy),0.d0)*idt)

LLn((j-1+ny/2)*nx+i+nx/2)=LLn((j-1+ny/2)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i-1)*dx,0.d0,ymin+dfloat(j-1)*dy),0.d0)*idt)

	enddo

enddo
!$OMP END PARALLEL DO

end subroutine propag2Dxz

!****************************************************************************************
!****************************************************************************************
subroutine Propag3D(LLn,xmin,dx,nx,ymin,dy,ny,zmin,dz,nz,idt,EEx,EEy,EEz,NS,coord,znuc)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nx,ny,nz
complex(kind=comp_16), intent(inout), dimension(nx*ny*nz) :: LLn
real(kind=real_8), intent(in) :: xmin,dx,EEx,ymin,dy,EEy,zmin,dz,EEz
complex(kind=comp_16), intent(in) :: idt

integer(kind=int_4),intent(in) :: ns
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc

integer(kind=int_4) :: i,j,k

!$OMP PARALLEL DO
do i=1,nx/2

	do j=1,ny/2

		do k=1,nz/2

LLn((k-1)*nx*ny+(j-1)*nx+i)=LLn((k-1)*nx*ny+(j-1)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy + (zmin+dfloat(k+nz/2-1)*dz)*EEz ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i+nx/2-1)*dx,ymin+dfloat(j+ny/2-1)*dy,zmin+dfloat(k+nz/2-1)*dz),0.d0)*idt)

LLn((k-1)*nx*ny+(j-1)*nx+i+nx/2)=LLn((k-1)*nx*ny+(j-1)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy + (zmin+dfloat(k+nz/2-1)*dz)*EEz ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i-1)*dx,ymin+dfloat(j+ny/2-1)*dy,zmin+dfloat(k+nz/2-1)*dz),0.d0)*idt)
LLn((k-1)*nx*ny+(j-1+ny/2)*nx+i)=LLn((k-1)*nx*ny+(j-1+ny/2)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy + (zmin+dfloat(k+nz/2-1)*dz)*EEz ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i+nx/2-1)*dx,ymin+dfloat(j-1)*dy,zmin+dfloat(k+nz/2-1)*dz),0.d0)*idt)
LLn((k-1+nz/2)*nx*ny+(j-1)*nx+i)=LLn((k-1+nz/2)*nx*ny+(j-1)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy + (zmin+dfloat(k-1)*dz)*EEz ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i+nx/2-1)*dx,ymin+dfloat(j+ny/2-1)*dy,zmin+dfloat(k-1)*dz),0.d0)*idt)

LLn((k-1)*nx*ny+(j-1+ny/2)*nx+i+nx/2)=LLn((k-1)*nx*ny+(j-1+ny/2)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy + (zmin+dfloat(k+nz/2-1)*dz)*EEz ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i-1)*dx,ymin+dfloat(j-1)*dy,zmin+dfloat(k+nz/2-1)*dz),0.d0)*idt)
LLn((k-1+nz/2)*nx*ny+(j-1+ny/2)*nx+i)=LLn((k-1+nz/2)*nx*ny+(j-1+ny/2)*nx+i)*cdexp(dcmplx((xmin+dfloat(i+nx/2-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy + (zmin+dfloat(k-1)*dz)*EEz ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i+nx/2-1)*dx,ymin+dfloat(j-1)*dy,zmin+dfloat(k-1)*dz),0.d0)*idt)
LLn((k-1+nz/2)*nx*ny+(j-1)*nx+i+nx/2)=LLn((k-1+nz/2)*nx*ny+(j-1)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j+ny/2-1)*dy)*EEy + (zmin+dfloat(k-1)*dz)*EEz ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i-1)*dx,ymin+dfloat(j+ny/2-1)*dy,zmin+dfloat(k-1)*dz),0.d0)*idt)

LLn((k-1+nz/2)*nx*ny+(j-1+ny/2)*nx+i+nx/2)=LLn((k-1+nz/2)*nx*ny+(j-1+ny/2)*nx+i+nx/2)*cdexp(dcmplx((xmin+dfloat(i-1)*dx)*EEx + (ymin+dfloat(j-1)*dy)*EEy + (zmin+dfloat(k-1)*dz)*EEz ,0.d0)*idt)*cdexp(dcmplx(VNE(NS,coord,znuc,xmin+dfloat(i-1)*dx,ymin+dfloat(j-1)*dy,zmin+dfloat(k-1)*dz),0.d0)*idt)
		enddo

	enddo

enddo
!$OMP END PARALLEL DO

end subroutine propag3D

!****************************************************************************************
!****************************************************************************************
subroutine Upp(LL,nq,np,nx,ny,nz,My_Desc_Handle3,CC,CCI,xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,idt,EEx,EEy,EEz,NS,coord,znuc)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nq,np,nx,ny,nz,NS
complex(kind=comp_16), intent(inout), dimension(np+nq,nq) :: LL
complex(kind=comp_16), intent(in), dimension(np+nq,np+nq) :: CC,CCI
type(DFTI_DESCRIPTOR), Pointer :: My_Desc_Handle3
real(kind=real_8), intent(in) :: xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,EEx,EEy,EEz
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc
complex(kind=comp_16), intent(in) :: idt

complex(kind=comp_16), allocatable, dimension(:,:) :: LLc,LLp
complex(kind=comp_16), allocatable, dimension(:) :: LLn
complex(kind=comp_16) :: alpha, beta
Integer(kind=int_4) :: isgn,n,i,j,k,ndim,cas
real(kind=real_8) :: sfact,kx,ky,kz
alpha=dcmplx(1.d0,0.d0)
beta=dcmplx(0.d0,0.d0)

allocate(LLc(np+nq,nq))
allocate(LLp(np+nq,nq))
 call checkcase(nx,ny,nz,cas)

! Extract P components from LL
! call mkl_zomatcopy('C','N',np+nq,nq,alpha,LL,np+nq,LLc,np+nq)
!LLc(1:nq,:)=dcmplx(0.d0,0.d0)

LLc = dcmplx(0.d0,0.d0)
!$OMP PARALLEL DO
	do i = 1+nq, nq+np
	do j = 1, nq
LLc(i,j)=LL(i,j)
enddo
enddo
!$OMP END PARALLEL DO

! Tranformation : OPW -> PW
 call zgemm('N','N',np+nq,nq,np+nq,alpha,CC,np+nq,LLc,np+nq,beta,LLp,np+nq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Propagation under Upp : split-operator
allocate(LLn(np))

do n=1,nq

! Extract column n 
	 call zcopy(np, LLp(1+nq:nq+np,n),1,LLn,1)

! Apply Exp(-i T dt/2), diagonal in PW rep
!$OMP PARALLEL DO
	do i = 1, nx
	do j = 1, ny
	do k = 1, nz

	LLn((k-1)*nx*ny+(j-1)*nx+i) = LLn((k-1)*nx*ny+(j-1)*nx+i)*cdexp(idt*dcmplx(0.25d0*((kxmin+(i-1)*dkx)**2.d0) + 0.25d0*((kymin+(j-1)*dky)**2.d0) + 0.25d0*((kzmin+(k-1)*dkz)**2.d0),0.d0))

	enddo
	enddo
	enddo
!$OMP END PARALLEL DO

! Transformation : k -> r
	isgn=1
	sfact=1.d0
	 call FFT_calc(LLn,np,sfact,isgn,My_Desc_Handle3)

! Apply Exp(-i r.E dt)
	if(cas.eq.1)then
	 call Propag1Dx(LLn,xmin,dx,nx,idt,EEx,NS,coord,znuc)
	endif
	if(cas.eq.2)then
	 call Propag1Dy(LLn,ymin,dy,ny,idt,EEy,NS,coord,znuc)
	endif
	if(cas.eq.3)then
	 call Propag1Dz(LLn,zmin,dz,nz,idt,EEz,NS,coord,znuc)
	endif

	if(cas.eq.4)then
	 call Propag2Dxy(LLn,xmin,dx,nx,ymin,dy,ny,idt,EEx,EEy,NS,coord,znuc)
	endif
	if(cas.eq.5)then
	 call Propag2Dyz(LLn,ymin,dy,ny,zmin,dz,nz,idt,EEy,EEz,NS,coord,znuc)
	endif
	if(cas.eq.6)then
	 call Propag2Dxz(LLn,xmin,dx,nx,zmin,dz,nz,idt,EEx,EEz,NS,coord,znuc)
	endif

	if(cas.eq.7)then
	 call Propag3D(LLn,xmin,dx,nx,ymin,dy,ny,zmin,dz,nz,idt,EEx,EEy,EEz,NS,coord,znuc)
	endif

! Transformation : r -> k
	isgn=-1
	sfact=1.d0/dfloat(np)
	 call FFT_calc(LLn,np,sfact,isgn,My_Desc_Handle3)

! Apply Exp(-i T dt/2)
!$OMP PARALLEL DO
	do i = 1, nx
	do j = 1, ny
	do k = 1, nz

	LLn((k-1)*nx*ny+(j-1)*nx+i) = LLn((k-1)*nx*ny+(j-1)*nx+i)*cdexp(idt*dcmplx(0.25d0*((kxmin+(i-1)*dkx)**2.d0) + 0.25d0*((kymin+(j-1)*dky)**2.d0) + 0.25d0*((kzmin+(k-1)*dkz)**2.d0),0.d0))

	enddo
	enddo
	enddo
!$OMP END PARALLEL DO

! Update column n of LLp
	 call zcopy(np, LLn,1,LLp(1+nq:nq+np,n),1)

enddo

deallocate(LLn)

! Transformation : PW -> OPW
 call zgemm('N','N',np+nq,nq,np+nq,alpha,CCI,np+nq,LLp,np+nq,beta,LLc,np+nq)
deallocate(LLp)
! Update the P components of LL after propagation
! call mkl_zomatcopy('C','N',np,nq,alpha,LLc(1+nq:np+nq,:),np,LL(1+nq:np+nq,:),np)
!$OMP PARALLEL DO
	do i = 1+nq, nq+np
	do j = 1, nq
LL(i,j)=LLc(i,j)
enddo
enddo
!$OMP END PARALLEL DO

deallocate(LLc)
end subroutine Upp
!****************************************************************************************
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
subroutine checkcase(nx,ny,nz,cas)
!****************************************************************************************
!****************************************************************************************
use basics

implicit none
integer(kind=int_4),intent(in) :: nx,ny,nz
integer(kind=int_4),intent(out) :: cas
integer(kind=int_4) :: ndim
! Check the number of dimensions
ndim=0
if(nx.gt.1)then
ndim=ndim+1
endif
if(ny.gt.1)then
ndim=ndim+1
endif
if(nz.gt.1)then
ndim=ndim+1
endif

if(ndim.eq.1)then
	if(nx.ne.1)then
	cas=1
	endif
	if(ny.ne.1)then
	cas=2
	endif
	if(nz.ne.1)then
	cas=3
	endif

endif

if(ndim.eq.2)then
	if(nx.eq.1)then
	cas=5
	endif
	if(ny.eq.1)then
	cas=6
	endif
	if(nz.eq.1)then
	cas=4
	endif
endif

if(ndim.eq.3)then
cas=7
endif

end subroutine checkcase
!****************************************************************************************
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
subroutine writeorbitals(phi,phik,nx,ny,nz,xmin,ymin,zmin,dx,dy,dz,kxmin,kymin,kzmin,dkx,dky,dkz,nw,nq,np)
!****************************************************************************************
!****************************************************************************************
use basics
Integer(kind=int_4), intent(inout) :: nx,ny,nz,nw,nq,np
Real(kind=real_8), intent(in) :: xmin,dx,ymin,dy,zmin,dz,kxmin,dkx,kymin,dky,kzmin,dkz
complex(kind=comp_16), intent(in), dimension(np,nq) :: phi,phik
Real(kind=real_8) :: kx,ky,kz
Integer(kind=int_4) :: i,j,k,wn
complex(kind=comp_16), allocatable, dimension(:) :: ctmp
complex(kind=comp_16), allocatable, dimension(:,:) :: ctmpa
!variables needed to interface with the Fourier Transform Routines
type(DFTI_DESCRIPTOR), Pointer :: My_Desc_Handle1,My_Desc_Handle2,My_Desc_Handle3
Integer(kind=int_4) :: L1(1),L2(2),L3(3)
Integer(kind=int_4) :: isgn
real(kind=real_8), allocatable, dimension(:) :: rtmp
real(kind=real_8) :: sfact,test7,test8,test6

open(91,file='files/phix0',status='unknown')

!do i = 1, nx
do i = 1+nx/2, 1+nx/2
do j = 1+ny/2, 1+ny/2
!do k = 1+nz/2, 1+nz/2
do k = 1, nz

kx=xmin+(i-1)*dx
ky=ymin+(j-1)*dy
kz=zmin+(k-1)*dz

write(91,108)kx,ky,kz,(dreal(phi((k-1)*nx*ny+(j-1)*nx+i,wn)),wn=1,nq)

enddo
enddo
enddo

open(92,file='files/phik0',status='unknown')
open(93,file='files/phikt',status='unknown')

allocate(ctmp(np))
allocate(ctmpa(np,nq))
 call Init_fft_3D(nx,ny,nz,My_Desc_Handle3,L3)
sfact=(dx/dsqrt(2.d0*pi))*(dy/dsqrt(2.d0*pi))*(dz/dsqrt(2.d0*pi))
!sfact=(dkx/dsqrt(2.d0*pi))*(dky/dsqrt(2.d0*pi))*(dkz/dsqrt(2.d0*pi))
isgn=-1

do wn=1,nq
 ctmp(:)=phi(:,wn)


 call FFT_calc(ctmp,np,sfact,isgn,My_Desc_Handle3)
 call sort_GEN(ctmp,nx,ny,nz)

 ctmpa(:,wn)=ctmp(:)

enddo

allocate(rtmp(np))

do i =1, nq

rtmp(:)=cdabs(phi(:,i))**2.d0
 call integrate(rtmp, test8, dx*dy*dz, np)

rtmp(:)=cdabs(phik(:,i))**2.d0
 call integrate(rtmp, test7, dkx*dky*dkz, np)

rtmp(:)=cdabs(ctmpa(:,i))**2.d0
 call integrate(rtmp, test6, dkx*dky*dkz, np)
! call integrate(rtmp, test6, dx*dy*dz, np)

write(*,'(A16,I1,A2,F7.4,F7.4,F7.4)')'Norm of orbital ',i,':',test8, test7, test6
enddo


!do i=1,nx
do i = 1+nx/2, 1+nx/2
do j = 1+ny/2, 1+ny/2
!do k = 1+nz/2, 1+nz/2
do k=1,nz

kx=kxmin+(i-1)*dkx
ky=kymin+(j-1)*dky
kz=kzmin+(k-1)*dkz

write(92,108)kx,ky,kz,(cdabs(phik((k-1)*nx*ny+(j-1)*nx+i,wn))**2.d0, wn=1,nw)
write(93,108)kx,ky,kz,(cdabs(ctmpa((k-1)*nx*ny+(j-1)*nx+i,wn))**2.d0, wn=1,nw)
!write(93,108)xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,(cdabs(ctmpa((k-1)*nx*ny+(j-1)*nx+i,wn))**2.d0, wn=1,nw)
!write(92,108)kx,ky,kz,dreal(phik((k-1)*nx*ny+(j-1)*nx+i,nw)),dimag(phik((k-1)*nx*ny+(j-1)*nx+i,nw))
!write(93,108)kx,ky,kz,dreal(ctmp((k-1)*nx*ny+(j-1)*nx+i)),dimag(ctmp((k-1)*nx*ny+(j-1)*nx+i))
enddo
enddo
enddo

 close(92)
! close(93)

108   FORMAT(8E18.6E3)
end subroutine writeorbitals
!****************************************************************************************
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
subroutine transmuk(muk,mu,CC,norb,nq,np)
! Transforms muk using the GS Coeffs CC (input matrix is overwritten)
!****************************************************************************************
!****************************************************************************************
use basics
Integer(kind=int_4), intent(in) :: nq,np,norb
Complex(kind=comp_16), dimension(np,nq), intent(inout) :: muk
Complex(kind=comp_16), dimension(np+nq,np+nq), intent(in) :: CC
real(kind=real_8), dimension(norb,norb), intent(in) :: mu
Integer(kind=int_4) :: i,j,k,l
complex(kind=comp_16), allocatable, dimension(:) :: allmu,rowj

allocate(allmu(np+nq))
! Extraction of desired elements
do i=1,nq
! Q space
	call zcopy(nq, mu(1:nq,i),1,allmu(1:nq),1)
! P space
	call zcopy(np, muk(1:np,i),1,allmu(1+nq:nq+np),1)

! Multiplication by the (j+nq)th row-vector of CC^H
	do j=1,np
	allocate(rowj(j+nq))
	call zcopy(nq+j, CC(1:nq+j,j+nq),1,rowj,1)
	muk(j,i) = dotc(rowj, allmu(1:j+nq))
	deallocate(rowj)
	enddo

enddo
deallocate(allmu)

end subroutine transmuk
!****************************************************************************************
!****************************************************************************************

!****************************************************************************************
subroutine transmotest(matin,matout,LL,nq,np)
!****************************************************************************************
!****************************************************************************************

use basics
use math
implicit none
integer(kind=int_4), intent(inout) :: nq,np
complex(kind=comp_16), intent(in), dimension(nq+np,nq) :: LL
complex(kind=comp_16), intent(out), dimension(nq,nq) :: matout
real(kind=real_8), intent(in), dimension(nq,nq) :: matin
complex(kind=comp_16) :: alpha
integer(kind=int_4) :: i,j,r,s

do i=1,nq
do j=1,nq

alpha=dcmplx(0.d0,0.d0)


do r=1,nq
do s=1,nq
alpha = alpha + dconjg(LL(r,i))*LL(s,j)*dcmplx(matin(r,s),0.d0)
enddo
enddo


matout(i,j)=dcmplx(chop(dreal(alpha)),chop(dimag(alpha)))

enddo
enddo

end subroutine transmotest

!****************************************************************************************
!****************************************************************************************
Function VNE(NS,coord,znuc,xx,yy,zz)
!****************************************************************************************
!****************************************************************************************

Use Basics

Implicit None
real(kind=real_8) :: VNE
integer(kind=int_4),intent(in) :: NS
real(kind=real_8), intent(in), dimension(NS,3) :: coord
real(kind=real_8), intent(in), dimension(NS) :: znuc
real(kind=real_8), intent(in) :: xx,yy,zz

integer(kind=int_4) :: n
real(kind=real_8) :: Rn

VNE = 0.d0

do n=1,ns
Rn=dsqrt((xx-coord(n,1))**2.d0 + (yy-coord(n,2))**2.d0 + (zz-coord(n,3))**2.d0)
VNE = VNE + (-1.d0*znuc(n)/(Rn+0.25d0))
enddo

end function

end module routines
