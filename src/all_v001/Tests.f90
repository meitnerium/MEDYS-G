Module Tests
!Gere les tests
use Math
use Basics

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transmat_prim2AO(Mat,center,EXPO,lmn,NS,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use basics

implicit none

! variables I/O
Real(kind=real_8), Intent(out), Allocatable, Dimension(:,:):: Mat
Real(kind=real_8), Intent(out), Allocatable, Dimension(:,:):: center
Real(kind=real_8), Intent(out), Allocatable, Dimension(:)::EXPO
Integer(kind=int_4), Intent(out), Allocatable, Dimension(:,:)::lmn
Integer(kind=int_4), Intent(in)			::NS
Real(kind=real_8), Dimension(:,:),Intent(inout)	::ZET,ETA,coord
Integer(kind=int_4), Dimension(:),Intent(in)	::ICONU
Integer(kind=int_4), Dimension(:),Intent(in)	::NF
Integer(kind=int_4), Dimension(:),Intent(in)	::LMNP1
Integer(kind=int_4), Intent(in)			::NCONS


! variables internes 
integer(kind=int_4)				:: iorb,iprim,inR,counter1,counter2,counter3,IS,i,j,k,r,s,m,p,imu
integer(kind=int_4)				:: PrimSet,nuPrim,l
real(kind=real_8)				:: ai,ci,dkit



counter1=0 
counter2=0 !number of atomic orbitals
counter3=0 !number of primitives

Do IS=1,NS
Do j=1,NF(IS)
	counter1=counter1+1
	select case(LMNP1(counter1))
	!Cas orbital S
	case(1)
	counter2=counter2+1
		do i=1,iconu(counter1)
			counter3=counter3+1
		end do

	!Cas orbital P
	case(2)
	do l=1,3
		counter2=counter2+1
		do i=1,iconu(counter1)
			counter3=counter3+1
		end do
	end do

	!Cas orbital D
	case(3)
	do l=1,5
		counter2=counter2+1
		do i=1,iconu(counter1)
			counter3=counter3+1
		end do
	end do
	end select
end do
end do
allocate(Mat(counter2,counter3),lmn(counter3,3),center(counter3,3),EXPO(counter3))

counter1=0
counter2=0
counter3=0
Do IS=1,NS
Do j=1,NF(IS)
	counter1=counter1+1
	select case(LMNP1(counter1))
	case(1)
!		Cas d'une orbitale S
		counter2=counter2+1
		do i=1,iconu(counter1)
			counter3=counter3+1
			
			center(counter3,:)=coord(IS,:)
			
			lmn(counter3,:)=0

			Mat(counter2,counter3)=eta(i,counter1)

			EXPO(counter3)=zet(i,counter1)
		end do
	case(2)
!		Cas d'une orbitale P
		do l=1,3
		counter2=counter2+1
			do i=1,iconu(counter1)
				counter3=counter3+1
				
				center(counter3,:)=coord(IS,:)
			
				
				If(l==1)then
				lmn(counter3,:)=[1,0,0]

				Else If(l==2)then 
				lmn(counter3,:)=[0,1,0]

				Else If(l==3)then
				lmn(counter3,:)=[0,0,1]

				end if

				Mat(counter2,counter3)=eta(i,counter1)
				EXPO(counter3)=zet(i,counter1)

			end do
		end do

	case(3)
!		Cas d'une orbitale D
STOP 'LES ORBITALES D NE SONT PAS DEVELLOPEES, TRAVAILLEZ LA SUBROUTINE' 
		do l=1,5
		counter2=counter2+1
			do i=1,iconu(counter1)
				counter3=counter3+1

				center(counter3,:)=coord(IS,:)
			
				If(l==1)then
				lmn(counter3,:)=[1,0,0]

				Else If(l==2)then 
				lmn(counter3,:)=[0,1,0]

				Else If(l==3)then
				lmn(counter3,:)=[0,0,1]

				Else If(l==4)then 
				lmn(counter3,:)=[0,1,0]

				Else If(l==5)then
				lmn(counter3,:)=[0,0,1]

				end if

				Mat(counter2,counter3)=eta(i,counter1)
				EXPO(counter3)=zet(i,counter1)
				
			end do
		end do
	end select
end do
end do

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Cartesiangaussian(x,y,z,alpha,Rx,Ry,Rz,l,m,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use math
implicit none

Real(kind=real_8)		:: Cartesiangaussian
Real(kind=real_8),intent(in)	:: x,y,z,Rx,Ry,Rz,alpha
Integer(kind=int_4),intent(in)	:: l,m,n
Integer(kind=int_4)		:: l2,m2,n2

Cartesiangaussian =(x-Rx)**(l) *((y-Ry)**(m) )*((z-Rz)**(n) )* Exp(- alpha*((x-Rx)**(2.d0)+(y-Ry)**(2.d0)+(z-Rz)**(2.d0)) )
l2=2*l
m2=2*m
n2=2*n
!Normalization
Cartesiangaussian = Cartesiangaussian * (2.d0*alpha/pi)**(3.d0/4.d0)*((8.d0*alpha)**(l+m+n)*fact(l,1)*fact(m,1)*fact(n,1))/(fact(l2,1)*fact(m2,1)*fact(n2,1))
End function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Print_d_matrix(mat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Real(kind=real_8),dimension(:,:),intent(in)::mat
   Integer(kind=int_4)::nn1,nn2,i
   nn1=Size(Mat(:,1))
   nn2=Size(Mat(1,:))

 Do i=1,nn1 
     write(*,*) mat(i,:)
 End do
end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Print_i_matrix(mat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Integer(kind=int_4),dimension(:,:),intent(in)::mat
   Integer(kind=int_4)::nn1,nn2,i
   nn1=Size(Mat(:,1))
   nn2=Size(Mat(1,:))

 Do i=1,nn1 
     write(*,*) mat(i,:)
 End do

end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_transmat_prim2AO(Mat,prim_center,prim_EXPO,prim_lmn)
! Linear combinaison of primitives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real(kind=real_8),dimension(:,:),intent(in)::Mat,prim_center
Real(kind=real_8),dimension(:),intent(in)::prim_EXPO
Integer(kind=int_4),dimension(:,:),intent(in)::prim_lmn

 

write(*,*) "Mat: "
    call print_d_matrix(mat)
write(*,*) "prim_center: "
    call print_d_matrix(prim_center)
write(*,*) "prim_expo: "
    write(*,*) prim_expo
write(*,*) "prim_lmn: "
    call print_i_matrix(prim_lmn)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine export_orbitals(Mat,prim_center,prim_EXPO,prim_lmn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Real(kind=real_8),dimension(:,:),intent(in)::Mat,prim_center
Real(kind=real_8),dimension(:),intent(in)::prim_EXPO
Integer(kind=int_4),dimension(:,:),intent(in)::prim_lmn
Real(kind=real_8)::AOValue
Real(kind=real_8):: xmax,ymax,zmax,dx,dy,dz,rx,ry,rz
Integer(kind=int_4)::nx,ny,nz,ix,iy,iz,nao,iao,nprim,iprim

xmax=10
ymax=10
zmax=10
nx= 100
ny=100
nz=100

nao= size(Mat(:,1))
nprim= size(Mat(1,:))
dx=(2.d0*xmax)/(nx-1)
dy=(2.d0*ymax)/(ny-1)
dz=(2.d0*zmax)/(nz-1)
write(*,*)'AO nao',nao

open(unit=71,file="AO1_grid.dat",status="new")
open(unit=72,file="AO2_grid.dat",status="new")
open(unit=73,file="AO3_grid.dat",status="new")
open(unit=74,file="AO4_grid.dat",status="new")
open(unit=75,file="AO5_grid.dat",status="new")
open(unit=76,file="AO6_grid.dat",status="new")
open(unit=77,file="AO7_grid.dat",status="new")
open(unit=78,file="AO8_grid.dat",status="new")
open(unit=79,file="AO9_grid.dat",status="new")
open(unit=80,file="AO10_grid.dat",status="new")
open(unit=81,file="AO11_grid.dat",status="new")
open(unit=82,file="AO12_grid.dat",status="new")
open(unit=83,file="AO13_grid.dat",status="new")
open(unit=84,file="AO14_grid.dat",status="new")
open(unit=85,file="AO15_grid.dat",status="new")
open(unit=86,file="AO16_grid.dat",status="new")
open(unit=87,file="AO17_grid.dat",status="new")

Do iao=1,nao
Do ix=1,nx
Do iy=1,ny
Do iz=1,nz
AOValue=0.d0
	rx=-xmax+ix*dx
	ry=-ymax+iy*dy
	rz=-zmax+iz*dz
	Do iprim=1, nprim
		!write(*,*)  rx,ry ,rz, (/( mat(iao,iprim)*cartesianGaussian(rx,ry,rz,prim_expo(iprim),prim_center(iprim,1),prim_center(iprim,2),prim_center(iprim,3),prim_lm(iprim,1),prim_lmn(iprim,2),prim_lmn(iprim,3)),iao=1,nao)/)
		AOValue = AOValue + mat(iao,iprim)*cartesianGaussian(rx,ry,rz,prim_expo(iprim),prim_center(iprim,1),prim_center(iprim,2),prim_center(iprim,3),prim_lmn(iprim,1),prim_lmn(iprim,2),prim_lmn(iprim,3))
	end do

Select case(iao)
	case (1)
	write(71,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (2)
	write(72,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (3)
	write(73,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (4)
	write(74,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (5)
	write(75,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (6)
	write(76,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue	
	
	case (7)
	write(77,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (8)
	write(78,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (9)
	write(79,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (10)
	write(80,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (11)
	write(81,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (12)
	write(82,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (13)
	write(83,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (14)
	write(84,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (15)
	write(85,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (16)
	write(86,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

	case (17)
	write(87,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,AOValue

end select

end do
end do
end do
end do

close(71)
close(72)
close(73)
close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine export_MOs(Mat,prim_center,prim_EXPO,prim_lmn,lcao)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use math

Real(kind=real_8),dimension(:,:),intent(inout)::Mat,prim_center,lcao
Real(kind=real_8),dimension(:),intent(in)::prim_EXPO
Integer(kind=int_4),dimension(:,:),intent(in)::prim_lmn
Real(kind=real_8),Allocatable, Dimension(:,:):: MOmat,temp
Real(kind=real_8)::MOValue
Real(kind=real_8):: xmax,ymax,zmax,dx,dy,dz,rx,ry,rz
Integer(kind=int_4)::nx,ny,nz,ix,iy,iz,nao,iao,jao,nprim,iprim,i,j,k

xmax=10
ymax=10
zmax=10
nx= 100
ny=100
nz=100

nao= size(Mat(:,1))
nprim= size(Mat(1,:))
dx=(2.d0*xmax)/(nx-1)
dy=(2.d0*ymax)/(ny-1)
dz=(2.d0*zmax)/(nz-1)

!allocate(MOmat(nao,nao),temp(nao,nao))
!MOmat=0.d0
!temp=0.d0

!call gemm(mat,lcao,MOmat)


open(unit=71,file="MO1_grid.dat",status="new")
open(unit=72,file="MO2_grid.dat",status="new")
open(unit=73,file="MO3_grid.dat",status="new")
open(unit=74,file="MO4_grid.dat",status="new")
open(unit=75,file="MO5_grid.dat",status="new")
open(unit=76,file="MO6_grid.dat",status="new")
open(unit=77,file="MO7_grid.dat",status="new")
open(unit=78,file="MO8_grid.dat",status="new")
open(unit=79,file="MO9_grid.dat",status="new")
open(unit=80,file="MO10_grid.dat",status="new")
open(unit=81,file="MO11_grid.dat",status="new")
open(unit=82,file="MO12_grid.dat",status="new")
open(unit=83,file="MO13_grid.dat",status="new")
open(unit=84,file="MO14_grid.dat",status="new")

Do iao=1,nao
Do ix=1,nx
Do iy=1,ny
Do iz=1,nz
MOValue=0.d0
	rx=-xmax+ix*dx
	ry=-ymax+iy*dy
	rz=-zmax+iz*dz
	Do iprim=1, nprim
	Do jao=1, nao
		MOValue = MOValue + lcao(jao,iao)*mat(jao,iprim)*cartesianGaussian(rx,ry,rz,prim_expo(iprim),prim_center(iprim,1),prim_center(iprim,2),prim_center(iprim,3),prim_lmn(iprim,1),prim_lmn(iprim,2),prim_lmn(iprim,3))
	end do
	end do

Select case(iao)
	case (1)
	write(71,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (2)
	write(72,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (3)
	write(73,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (4)
	write(74,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (5)
	write(75,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (6)
	write(76,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue	
	
	case (7)
	write(77,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (8)
	write(78,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (9)
	write(79,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (10)
	write(80,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (11)
	write(81,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (12)
	write(82,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue
	
	case (13)
	write(83,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

	case (14)
	write(84,'(F9.5,F9.5,F9.5,F26.14)') rx,ry,rz,MOValue

end select

end do
end do
end do
end do

close(71)
close(72)
close(73)
close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)
close(83)
close(84)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine export_overlap(Mat,nki,dki,kmin,nao)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Complex(kind=comp_16),dimension(:,:,:,:),intent(in)::Mat
Integer(kind=int_4),dimension(3),intent(in)::nki
Real(kind=real_8),dimension(3),intent(in)::dki,kmin
Integer(kind=int_4),intent(in)::nao

Real(kind=real_8)::kx,ky,kz
Integer(kind=int_4)::ikx,iky,ikz,iao




open(unit=71,file="Overlap1_grid.dat",status="new")
open(unit=72,file="Overlap2_grid.dat",status="new")
open(unit=73,file="Overlap3_grid.dat",status="new")
open(unit=74,file="Overlap4_grid.dat",status="new")
open(unit=75,file="Overlap5_grid.dat",status="new")
open(unit=76,file="Overlap6_grid.dat",status="new")
open(unit=77,file="Overlap7_grid.dat",status="new")
open(unit=78,file="Overlap8_grid.dat",status="new")
open(unit=79,file="Overlap9_grid.dat",status="new")
open(unit=80,file="Overlap10_grid.dat",status="new")
open(unit=81,file="Overlap11_grid.dat",status="new")
open(unit=82,file="Overlap12_grid.dat",status="new")
open(unit=84,file="Overlap14_grid.dat",status="new")
open(unit=86,file="Overlap16_grid.dat",status="new")

Do iao=1,nao
Do ikx=1,nki(1)
Do iky=1,nki(2)
Do ikz=1,nki(3)

kx=kmin(1)+ikx*dki(1)
ky=kmin(2)+iky*dki(2)
kz=kmin(3)+ikz*dki(3)

Select case(iao)
	case (1)
	write(71,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (2)
	write(72,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (3)
	write(73,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (4)
	write(74,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (5)
	write(75,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (6)
	write(76,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (7)
	write(77,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (8)
	write(78,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (9)
	write(79,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (10)
	write(80,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (11)
	write(81,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (12)
	write(82,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (14)
	write(84,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (16)
	write(86,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

end select

end do
end do
end do
end do

close(71)
close(72)
close(73)
close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)
close(84)
close(86)
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine export_TransitionX(Mat,nki,dki,kmin,nao)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Complex(kind=comp_16),dimension(:,:,:,:),intent(in)::Mat
Integer(kind=int_4),dimension(3),intent(in)::nki
Real(kind=real_8),dimension(3),intent(in)::dki,kmin
Integer(kind=int_4),intent(in)::nao

Real(kind=real_8)::kx,ky,kz
Integer(kind=int_4)::ikx,iky,ikz,iao




!open(unit=71,file="TransitionX1_grid.dat",status="new")
!open(unit=72,file="TransitionX2_grid.dat",status="new")
!open(unit=73,file="TransitionX3_grid.dat",status="new")
!open(unit=74,file="TransitionX4_grid.dat",status="new")
open(unit=75,file="TransitionX5_grid.dat",status="new")
open(unit=76,file="TransitionX6_grid.dat",status="new")
open(unit=77,file="TransitionX7_grid.dat",status="new")
open(unit=78,file="TransitionX8_grid.dat",status="new")
open(unit=79,file="TransitionX9_grid.dat",status="new")
open(unit=80,file="TransitionX10_grid.dat",status="new")
open(unit=81,file="TransitionX11_grid.dat",status="new")
open(unit=82,file="TransitionX12_grid.dat",status="new")

Do iao=1,nao
Do ikx=1,nki(1)
Do iky=1,nki(2)
Do ikz=1,nki(3)

kx=kmin(1)+ikx*dki(1)
ky=kmin(2)+iky*dki(2)
kz=kmin(3)+ikz*dki(3)

Select case(iao)
	!case (1)
	!write(71,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	!case (2)
	!write(72,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	!case (3)
	!write(73,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

!	case (4)
!	write(74,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (5)
	write(75,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (6)
	write(76,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (7)
	write(77,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (8)
	write(78,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (9)
	write(79,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (10)
	write(80,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (11)
	write(81,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (12)
	write(82,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (14)
	write(84,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (16)
	write(86,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

end select
end do
end do
end do
end do

!close(71)
!close(72)
!close(73)
!close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine export_TransitionY(Mat,nki,dki,kmin,nao)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Complex(kind=comp_16),dimension(:,:,:,:),intent(in)::Mat
Integer(kind=int_4),dimension(3),intent(in)::nki
Real(kind=real_8),dimension(3),intent(in)::dki,kmin
Integer(kind=int_4),intent(in)::nao

Real(kind=real_8)::kx,ky,kz
Integer(kind=int_4)::ikx,iky,ikz,iao




open(unit=71,file="TransitionY1_grid.dat",status="new")
open(unit=72,file="TransitionY2_grid.dat",status="new")
open(unit=73,file="TransitionY3_grid.dat",status="new")
open(unit=74,file="TransitionY4_grid.dat",status="new")
open(unit=75,file="TransitionY5_grid.dat",status="new")
open(unit=76,file="TransitionY6_grid.dat",status="new")
open(unit=77,file="TransitionY7_grid.dat",status="new")
open(unit=78,file="TransitionY8_grid.dat",status="new")
open(unit=79,file="TransitionY9_grid.dat",status="new")
open(unit=80,file="TransitionY10_grid.dat",status="new")
open(unit=81,file="TransitionY11_grid.dat",status="new")
open(unit=82,file="TransitionY12_grid.dat",status="new")

Do iao=1,nao
Do ikx=1,nki(1)
Do iky=1,nki(2)
Do ikz=1,nki(3)

kx=kmin(1)+ikx*dki(1)
ky=kmin(2)+iky*dki(2)
kz=kmin(3)+ikz*dki(3)

Select case(iao)
	case (1)
	write(71,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (2)
	write(72,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (3)
	write(73,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (4)
	write(74,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (5)
	write(75,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (6)
	write(76,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (7)
	write(77,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (8)
	write(78,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (9)
	write(79,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (10)
	write(80,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (11)
	write(81,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (12)
	write(82,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

end select

end do
end do
end do
end do

close(71)
close(72)
close(73)
close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine export_TransitionZ(Mat,nki,dki,kmin,nao)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Complex(kind=comp_16),dimension(:,:,:,:),intent(in)::Mat
Integer(kind=int_4),dimension(3),intent(in)::nki
Real(kind=real_8),dimension(3),intent(in)::dki,kmin
Integer(kind=int_4),intent(in)::nao

Real(kind=real_8)::kx,ky,kz
Integer(kind=int_4)::ikx,iky,ikz,iao




!open(unit=71,file="TransitionZ1_grid.dat",status="new")
!open(unit=72,file="TransitionZ2_grid.dat",status="new")
!open(unit=73,file="TransitionZ3_grid.dat",status="new")
!open(unit=74,file="TransitionZ4_grid.dat",status="new")
open(unit=75,file="TransitionZ5_grid.dat",status="new")
open(unit=76,file="TransitionZ6_grid.dat",status="new")
open(unit=77,file="TransitionZ7_grid.dat",status="new")
open(unit=78,file="TransitionZ8_grid.dat",status="new")
open(unit=79,file="TransitionZ9_grid.dat",status="new")
open(unit=80,file="TransitionZ10_grid.dat",status="new")
open(unit=81,file="TransitionZ11_grid.dat",status="new")
open(unit=82,file="TransitionZ12_grid.dat",status="new")

Do iao=1,nao
Do ikx=1,nki(1)
Do iky=1,nki(2)
Do ikz=1,nki(3)

kx=kmin(1)+ikx*dki(1)
ky=kmin(2)+iky*dki(2)
kz=kmin(3)+ikz*dki(3)

Select case(iao)
!	case (1)
!	write(71,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

!	case (2)
!	write(72,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))
!
!	case (3)
!	write(73,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

!	case (4)
!	write(74,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (5)
	write(75,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (6)
	write(76,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (7)
	write(77,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (8)
	write(78,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (9)
	write(79,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (10)
	write(80,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (11)
	write(81,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

	case (12)
	write(82,'(F9.5,F9.5,F9.5,F21.14,F21.14)') kx,ky,kz,real(Mat(ikx,iky,ikz,iao)),aimag(Mat(ikx,iky,ikz,iao))

end select

end do
end do
end do
end do

!close(71)
!close(72)
!close(73)
!close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine export_MOTransitionX(Mat,nao)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Real(kind=real_8),dimension(:,:),intent(in)::Mat

Integer(kind=int_4),intent(in)::nao

Integer(kind=int_4)::iao,jao




open(unit=71,file="MOTransitionX1_grid.dat",status="new")
open(unit=72,file="MOTransitionX2_grid.dat",status="new")
open(unit=73,file="MOTransitionX3_grid.dat",status="new")
open(unit=74,file="MOTransitionX4_grid.dat",status="new")
open(unit=75,file="MOTransitionX5_grid.dat",status="new")
open(unit=76,file="MOTransitionX6_grid.dat",status="new")
open(unit=77,file="MOTransitionX7_grid.dat",status="new")
open(unit=78,file="MOTransitionX8_grid.dat",status="new")
open(unit=79,file="MOTransitionX9_grid.dat",status="new")
open(unit=80,file="MOTransitionX10_grid.dat",status="new")
open(unit=81,file="MOTransitionX11_grid.dat",status="new")
open(unit=82,file="MOTransitionX12_grid.dat",status="new")

Do iao=1,nao
Do jao=1,nao
Select case(iao)
	case (1)
	write(71,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (2)
	write(72,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (3)
	write(73,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (4)
	write(74,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (5)
	write(75,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (6)
	write(76,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (7)
	write(77,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (8)
	write(78,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (9)
	write(79,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (10)
	write(80,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (11)
	write(81,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (12)
	write(82,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (14)
	write(84,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

	case (16)
	write(86,'(F9.5,F9.5,F9.5,F21.14,F21.14)') Mat(iao,jao)

end select

end do
end do

close(71)
close(72)
close(73)
close(74)
close(75)
close(76)
close(77)
close(78)
close(79)
close(80)
close(81)
close(82)

end subroutine
end module Tests

