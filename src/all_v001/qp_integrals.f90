module qp_integrals

! construction des intégrales mixtes orbitales moléculaires / orbitales du continuum

use Basics
!use blas95
!use f95_precision
use blas95
use lapack95
use f95_precision
use math

contains

!****************************************************************************************
!****************************************************************************************
!subroutine SuperMatrix_OM_PW(ski,muxki,muyki,muzki,nR,nki,kimin,dki,norb_lect,lcao,NS,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord)
!TODO : deleted Nr since already in variables
subroutine SuperMatrix_OM_PW(ski,muxki,muyki,muzki,nki,kimin,dki,norb_lect,lcao,NS,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord,lp,ld)
!AJOUT JN: lp(3), ld(6)^
!!calcul des intégrales mixtes gaussiennes/ondes planes   (recouvrement et moment dipolaire)
!
!	input: 	-nR::		nombre de points de la base nucléaire
!		-NS::	nombre de centres nucléaires
!       	-nkx,nky,nkz::		nombre d'ondes planes 
!		-kxmin,kymin,kzmin::		
!		-dkx,dky,dkz::		
!		-norb::    	nombre d'orbitales atomiques 
!				(iorb(Norb),set(1),direction(1))
!		-LCAO::		coefficients LCAO
!		-NS/ETA/ZET/,NCONS,ICONU,NF,LMNP1/coord
!		
!	output:	-ski::	 	intégrale de recouvement <k|i>, i=indice de l'orbitale moléculaire
!		-muki::		moment de transition dipolaire <k|z|i>
!****************************************************************************************
!****************************************************************************************

use basics
use variables
!TODO : use variables a été ajouter

implicit none

! variables I/O
!Integer(kind=int_4), Intent(in)			::nR,NS
Integer(kind=int_4), Intent(in)			::NS
! TODO Nr is already on variables mod
Real(kind=real_8), Dimension(:,:),Intent(in)	::lcao
Real(kind=real_8), Dimension(:,:),Intent(inout)	::ZET,ETA,coord
Integer(kind=int_4), Dimension(:),Intent(in)	::ICONU
Integer(kind=int_4), Dimension(:),Intent(in)	::NF
Integer(kind=int_4), Dimension(:),Intent(in)	::LMNP1
Integer(kind=int_4), Intent(in)			::NCONS
Integer(kind=int_4), Dimension(3),Intent(in)	::nki,lp
Integer(kind=int_4), Dimension(6),Intent(in)	::ld
Real(kind=real_8), Dimension(3),Intent(in)	::kimin,dki
Integer(kind=int_4), Intent(in)			:: norb_lect

complex(kind=comp_16),dimension(:,:),intent(out):: ski,muxki,muyki,muzki

! variables internes 
integer(kind=int_4)				:: ak,bk,ck,iorb,iprim,inR,counter1,counter2,IS,i,j,k,r,s,m,p,imu,kpm
integer(kind=int_4)				:: PrimSet,nuPrim,l,lj
real(kind=real_8)				::center, kx,ky,kz,ai,ci,dkit
real(kind=real_8),dimension(3)			::Ri
complex(kind=comp_16),dimension(:,:),allocatable:: tempski,tempmuxki,tempmuyki,tempmuzki
complex(kind=comp_16),dimension(norb_lect,norb_lect)	::lcao_Con_complex

nk=nki(1)*nki(2)*nki(3)
allocate(tempski(nk,norb_lect),tempmuxki(nk,Norb_lect),tempmuyki(nk,norb_lect),tempmuzki(nk,norb_lect))
! initialisation

tempski=dcmplx(0.d0,0.d0)
tempmuxki=dcmplx(0.d0,0.d0)
tempmuyki=dcmplx(0.d0,0.d0)
tempmuzki=dcmplx(0.d0,0.d0)
ski = dcmplx(0.d0,0.d0)
muxki=dcmplx(0.d0,0.d0)
muyki=dcmplx(0.d0,0.d0)
muzki=dcmplx(0.d0,0.d0)


counter1=0
counter2=0
Do IS=1,NS
  Do j=1,NF(IS)
	counter1=counter1+1
	select case(LMNP1(counter1))
	case(1)
!		Cas d'une orbitale S
		counter2=counter2+1
		do k=1,nki(1)
		  kx=(k-1)*dki(1)+kimin(1)
		  do p=1,nki(2)
			ky=(p-1)*dki(2)+kimin(2)
		    do m=1,nki(3)
			  kz=(m-1)*dki(3)+kimin(3)
              kpm=k+nki(1)*(p-1)+nki(1)*nki(2)*(m-1) !TODO This is the correct version
		      do i=1,iconu(counter1)
			    tempski(kpm,counter2)=tempski(kpm,counter2)+generalSki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),0,0)
			    tempmuxki(kpm,counter2)=tempmuxki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),0,0,1)
			    tempmuyki(kpm,counter2)=tempmuyki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),0,0,2)
			    tempmuzki(kpm,counter2)=tempmuzki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),0,0,3)
		      end do
		    end do
		  end do
		end do
	case(2)
!		Cas d'une orbitale P
		do l=1,3
			counter2=counter2+1
			do k=1,nki(1)
				kx=(k-1)*dki(1)+kimin(1)

			do p=1,nki(2)
				ky=(p-1)*dki(2)+kimin(2)
	
			do m=1,nki(3)
				kz=(m-1)*dki(3)+kimin(3)
                kpm=k+nki(1)*(p-1)+nki(1)*nki(2)*(m-1) !TODO This is ok
			do i=1,iconu(counter1)

				tempski(kpm,counter2)=tempski(kpm,counter2)+generalSki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),1,lp(l))


				tempmuxki(kpm,counter2)=tempmuxki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),lp(l),1,1)

				tempmuyki(kpm,counter2)=tempmuyki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),lp(l),1,2)

				tempmuzki(kpm,counter2)=tempmuzki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),lp(l),1,3)
			end do
			end do
			end do
			end do
		end do

	case(3)
!		Cas d'une orbitale D
		do l=1,6
			counter2=counter2+1
			do k=1,nki(1)
				kx=(k-1)*dki(1)+kimin(1)

			do p=1,nki(2)
				ky=(p-1)*dki(2)+kimin(2)
	
			do m=1,nki(3)
				kz=(m-1)*dki(3)+kimin(3)
                kpm=k+nki(1)*(p-1)+nki(1)*nki(2)*(m-1)
			do i=1,iconu(counter1)
				tempski(kpm,counter2)=tempski(kpm,counter2)+generalSki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),2,ld(l))
				
				tempmuxki(kpm,counter2)=tempmuxki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),ld(l),2,1)

				!tempmuyki(k,p,m,counter2)=tempmuyki(k,p,m,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),l,2,2)
				tempmuyki(kpm,counter2)=tempmuyki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),ld(l),2,2)
         !TODO Typing error?
				tempmuzki(kpm,counter2)=tempmuzki(kpm,counter2)+generalMUki(kx,ky,kz,zet(i,counter1),eta(i,counter1),coord(IS,:),ld(l),2,3)

			end do
			end do
			end do
			end do
		end do
	end select
end do
end do

! transformation des intégrales <k|o|iao> -> <k|o|imo>
lcao_Con_complex=dcmplx(0.d0,0.d0)
!lcao_Con_complex=dcmplx(Transpose(lcao),0.d0)
lcao_Con_complex(:,:)=Transpose(lcao(:,:))

!TODO comlex vector vector multiplication?
!TODO comented for test purpose
! lcao_Con_complex is a matrix
! to be tested
do ak=1,nk
	!call gemv(lcao_Con_complex,tempski(ak,:),ski(ak,:))
   ski(ak,:)=MATMUL(lcao_Con_complex, tempski(ak,:))
   !write(*,*) "test on line 186 qp_integral after gemv "
	!call gemv(lcao_Con_complex,tempmuxki(ak,:),muxki(ak,:))
   muxki(ak,:)=MATMUL(lcao_Con_complex, tempmuxki(ak,:))
	!call gemv(lcao_Con_complex,tempmuyki(ak,:),muyki(ak,:))
   muyki(ak,:)=MATMUL(lcao_Con_complex, tempmuyki(ak,:))
	!call gemv(lcao_Con_complex,tempmuzki(ak,:),muzki(ak,:))
   muzki(ak,:)=MATMUL(lcao_Con_complex, tempmuzki(ak,:))
end do



deallocate(tempski,tempmuxki,tempmuyki,tempmuzki)
end subroutine SuperMatrix_OM_PW

!****************************************************************************************
!****************************************************************************************
! General Overlap function
function generalSki(kx,ky,kz,ai,ci,Ri,l,lj)
!****************************************************************************************
!****************************************************************************************
implicit none
integer(kind=int_4), intent(in):: lj,l
complex(kind=comp_16) :: generalSki
real(kind=real_8), intent(in) :: kx,ky,kz,ai,ci
real(kind=real_8),dimension(3),intent(in)::Ri

If(l==0) then
generalSki=primski(kx,ky,kz,ai,ci,Ri)
end if

If(l==1) then
generalSki=primpaki(kx,ky,kz,ai,ci,Ri,lj)
end if

If(l==2) then
generalSki=primdaki(kx,ky,kz,ai,ci,Ri,lj)
end if

end function


!****************************************************************************************
!****************************************************************************************
!Calcul de :
!	<k|g0,0,0>
!Pour un deplacement nul
!Maintenant normalise
function Gk(kx,ky,kz,ai,Ri)
!****************************************************************************************
!****************************************************************************************

use basics

implicit none
complex(kind=comp_16) :: Gk
real(kind=real_8), intent(in) :: kx,ky,kz,ai
real(kind=real_8),dimension(3),intent(in)::Ri

!Gk=(2.d0*pi*ai)**(3.d0/4.d0)*(cdexp(dcmplx(-ai*(kx**2.d0+ky**2.d0+kz**2.d0),-(kx*Ri(1)+ky*Ri(2)+kz*Ri(3)))))

Gk=(2.d0*pi*ai)**(-3.d0/4.d0)*(cdexp(dcmplx(-(kx**2.d0+ky**2.d0+kz**2.d0)/(4.d0*ai),-(kx*Ri(1)+ky*Ri(2)+kz*Ri(3)))))

end function
!****************************************************************************************
!****************************************************************************************
!Calcul de :
!	<k|g l,m,n>
!Pour un deplacement non nul
function Gklmn(kx,ky,kz,ai,Ri,l,m,n)
!****************************************************************************************
!****************************************************************************************

use basics

implicit none
complex(kind=comp_16) :: Gklmn
real(kind=real_8), intent(in) :: kx,ky,kz,ai
integer,intent(in)::l,m,n
real(kind=real_8),dimension(3),intent(in)::Ri

Gklmn= Hnk(kx,l,ai)*Hnk(ky,m,ai)*Hnk(kz,n,ai)*Gk(kx,ky,kz,ai,Ri) 

Gklmn=dcmplx(0.d0,dsqrt(4.d0*ai))**(l+m+n)*Gklmn
!Gklmn=1.d0/sqrt(fact(2*l-1,2)*fact(2*m-1,2)*fact(2*n-1,2))*Gklmn

end function
!****************************************************************************************
!****************************************************************************************
! Overlap between a general gaussian type s function centered on atom i (at R)
! and a plane wave, exp(i kx x +i ky y + i kz z)  <K|si>
function primski(kx,ky,kz,ai,ci,Ri)
!****************************************************************************************
!****************************************************************************************

use basics

implicit none
complex(kind=comp_16) :: primski
real(kind=real_8), intent(in) :: kx,ky,kz,ai,ci
real(kind=real_8),dimension(3),intent(in)::Ri
real(kind=real_8)::Riz

primski= dcmplx(ci,0.d0)*Gklmn(kx,ky,kz,ai,Ri,0,0,0)

end function

!****************************************************************************************
!****************************************************************************************
! Overlap between a general gaussian type pa, a=x,y,or z function centered on atom i (at Raz)
! and a plane wave, exp(i kx x +i ky y + i kz z) <K|pi>
function primpaki(kx,ky,kz,ai,ci,Ri,lj)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
integer(kind=int_4),intent(in)::lj !i=indice de la gaussienne,  lj=x,y ou z
complex(kind=comp_16) :: primpaki
real(kind=real_8), intent(in) :: kx,ky,kz,ai,ci
real(kind=real_8),dimension(3),intent(in)::Ri
integer::l,m,n
real(kind=real_8)::Riz
real(kind=real_8):: ka

l=0
m=0
n=0

If(lj==1)then
 l=1
end if

If(lj==2)then 
 m=1
end if

If(lj==3)then
 n=1
end if

primpaki = Gklmn(kx,ky,kz,ai,Ri,l,m,n)

primpaki = dcmplx(ci,0.d0)*primpaki

end function

!****************************************************************************************
!****************************************************************************************
! Overlap between a general gaussian type da function centered on atom i (at Raz)
! and a plane wave, exp(i kx x +i ky y + i kz z) <K|dai>, a= xy xz, yz x^2y^2,z^2 (voir ordre dans intégrales argos)
function primdaki(kx,ky,kz,ai,ci,Ri,lj)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
integer(kind=int_4),intent(in)::lj 
complex(kind=comp_16) :: primdaki
real(kind=real_8), intent(in) :: kx,ky,kz,ai,ci
real(kind=real_8),dimension(3),intent(in)::Ri
real(kind=real_8)::Riz
real(kind=real_8):: n0
integer::l,m,n

l=0
m=0
n=0
If(lj==1) then
	!dx2
	l=2
end if

If(lj==2) then
	!dy2
	m=2
end if

If(lj==3) then
	!dz2
	n=2
end if

If(lj==4) then
	!dxy
	l=1
	m=1	
end if

If(lj==5) then
	!dxz 
	l=1
	n=1
end if
If(lj==6) then
	!dyz
	m=1
	n=1
end if

primdaki= dcmplx(ci,0.d0)*Gklmn(kx,ky,kz,ai,Ri,l,m,n)

!TODO not sure !
!primpaki = dcmplx(ci,0.d0)*primpaki
!primpaki = dcmplx(ci,0.d0)*primdaki

end function

!****************************************************************************************
!****************************************************************************************
! General transition moment function
function generalMUki(kx,ky,kz,ai,ci,Ri,lj,l,imu)
!****************************************************************************************
!****************************************************************************************
implicit none
integer(kind=int_4), intent(in):: lj,l,imu
complex(kind=comp_16) :: generalMUki
real(kind=real_8), intent(in) :: kx,ky,kz,ai,ci
real(kind=real_8),dimension(3),intent(in)::Ri

If(l==0) then
generalMUki=primmuski(kx,ky,kz,ai,ci,Ri,imu)
end if

If(l==1) then
generalMUki=primmupaki(kx,ky,kz,ai,ci,Ri,imu,lj)
end if

If(l==2) then
generalMUki=primmudaki(kx,ky,kz,ai,ci,Ri,imu,lj)
end if

end function

!****************************************************************************************
!****************************************************************************************
! Transition dipole along x, y or z between a general gaussian type s function centered on atom A 
! and a plane wave, exp(i kx x +i ky y + i kz z)  <K|r|si>
function primmuski(kx,ky,kz,ai,ci,Ri,imu)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
complex(kind=comp_16) :: primmuski
real(kind=real_8), intent(in) :: kx,ky,kz,ai,ci
integer(kind=int_4),intent(in)::imu
real(kind=real_8),dimension(3),intent(in)::Ri
real(kind=real_8)::Riz,ki
integer::dl,dm,dn

dl=0
dm=0
dn=0
if (imu==1) then
	Riz=Ri(1)
	primmuski=Gklmn(kx,ky,kz,ai,Ri,dl+1,dm,dn)

else if (imu==2) then
	Riz=Ri(2)
	primmuski=Gklmn(kx,ky,kz,ai,Ri,dl,dm+1,dn)

else if (imu==3) then
	Riz=Ri(3)
	primmuski=Gklmn(kx,ky,kz,ai,Ri,dl,dm,dn+1)
end if

primmuski=Riz*Gklmn(kx,ky,kz,ai,Ri,dl,dm,dn) + (4.d0*ai)**(-1.d0/2.d0)*primmuski

primmuski= dcmplx(ci,0.d0)*primmuski

end function

!****************************************************************************************
!****************************************************************************************
! Transition dipole along z between a general gaussian type pa, lj=x,y,or z function centered on atom i (at Riz)
!(all nucleus on z axis) 
! and a plane wave, exp(i kx x +i ky y + i kz z) <K|z|pai>
function primmupaki(kx,ky,kz,ai,ci,Ri,imu,lj)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
integer(kind=int_4),intent(in)::lj
complex(kind=comp_16) :: primmupaki
real(kind=real_8), intent(in) :: kx,ky,kz,ai,ci
integer(kind=int_4),intent(in)::imu
real(kind=real_8),dimension(3),intent(in)::Ri
real(kind=real_8)::Riz,ki
integer::l,m,n,dl,dm,dn

dl=0
dm=0
dn=0
l=0
m=0
n=0

if (imu==1) then
	Riz=Ri(1)
	dl=1

else if (imu==2) then
	Riz=Ri(2)
	dm=1

else if (imu==3) then
	Riz=Ri(3)
	dn=1

end if


if(lj==1) then 
	!mupxk
	l=1
end if

if(lj==2) then 
	!mupyk
	m=1
end if

If(lj==3) then
	!mupzk
	n=1
end if


primmupaki=Gklmn(kx,ky,kz,ai,Ri,l+dl,m+dm,n+dn)

primmupaki=Riz*Gklmn(kx,ky,kz,ai,Ri,l,m,n) + (4.d0*ai)**(-1.d0/2.d0)*primmupaki

primmupaki=dcmplx(ci,0.d0)*primmupaki
end function

!****************************************************************************************
!****************************************************************************************
! Transition dipole along z between a general gaussian type da function centered on atom i (at z = Riz) 
!(all nucleus on z axis) 
! and a plane wave, exp(i kx x +i ky y + i kz z) <K|z|dai>
function primmudaki(kx,ky,kz,ai,ci,Ri,imu,lj)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
integer(kind=int_4),intent(in)::lj
complex(kind=comp_16) :: primmudaki
real(kind=real_8), intent(in) :: kx,ky,kz,ai,ci
integer(kind=int_4),intent(in)::imu
real(kind=real_8),dimension(3),intent(in)::Ri
real(kind=real_8)::Riz,ki
integer::dl,dm,dn,l,m,n
real(kind=real_8)::n0

l=0
m=0
n=0
dl=0
dm=0
dn=0

if (imu==1) then
	Riz=Ri(1)
	dl=1

else if (imu==2) then
	Riz=Ri(2)
	dm=1

else if (imu==3) then
	Riz=Ri(3)
	dn=1
end if

If(lj==1) then
	!dx2
	l=2
end if

If(lj==2) then
	!dy2
	m=2
end if

If(lj==3) then
	!dz2
	n=2
end if

If(lj==4) then
	!dxy
	l=1
	m=1	
end if

If(lj==5) then
	!dxz 
	l=1
	n=1
end if
If(lj==6) then
	!dyz
	m=1
	n=1
end if

	primmudaki= Gklmn(kx,ky,kz,ai,Ri,l+dl,m+dm,n+dn)
	
	primmudaki=Riz*Gklmn(kx,ky,kz,ai,Ri,l,m,n) + (4.d0*ai)**(-1.d0/2.d0)*primmudaki

	primmudaki=dcmplx(ci,0.d0)*primmudaki



end function


!call k2phiK(muEPS,muB_EPS,ski,muki_EPS,orb_Q,nk)
!****************************************************************************************
!****************************************************************************************
! Passage de <K|z|PhiI> à <PhiK|z|PhiI>
subroutine k2phiK(muA,muB,ski,muki,norb_lect)
!****************************************************************************************
!****************************************************************************************
use variables
!use blas95
!use f95_precision
integer(kind=int_4),intent(in)::norb_lect
complex(kind=comp_16),dimension(:,:),intent(in)::ski,muki
complex(kind=comp_16),dimension(:,:),intent(in)::muA
complex(kind=comp_16),dimension(:,:),intent(out)::muB
complex(kind=comp_16),allocatable,dimension(:,:)::tempMu
integer(kind=int_4)::i,j

allocate(tempMu(nk,norb_lect))

muB=dcmplx(0.d0,0.d0)
 call gemm(ski(:,:),muA,tempMu)
muB= muki(:,:)-tempMu              
 

end subroutine k2phiK

!****************************************************************************************
!****************************************************************************************
! Calcul Hqp,Hpq et Hqppq. Indépendant du temps
subroutine Hnd_CSF(muB,Ers,Hqp,Hpq,Hqppq,dki)
    !TODO nk est dans modules variables donc pas obligatoire de le passer ici  
!****************************************************************************************
!****************************************************************************************
use variables

complex(kind=comp_16),dimension(:,:),intent(out)::Hqp,Hpq
complex(kind=comp_16),dimension(:,:),intent(out)::Hqppq
complex(kind=comp_16),dimension(:,:),intent(in)::muB
real(kind=real_8),dimension(:,:,:,:),intent(in)::Ers
Real(kind=real_8),Dimension(3),Intent(in)	::dki

integer(kind=int_4)::l,ll,i,j,jj,k,debut,fin,r,s
integer(kind=int_4)::kkx,kky,kkz,kk,jkk
Hqp=dcmplx(0.d0,0.d0)
Hpq=dcmplx(0.d0,0.d0)
Hqppq=dcmplx(0.d0,0.d0)

do ll=1,dimQ	
	do j=1,dimP
	   do kk=1,nk 
		jkk=(j-1)*nk+kk
		do i=1,orb_Q
			Hqp(ll,jkk) = Hqp(ll,jkk) + (muB(kk,i)*Ers(ll,j+dimQ,i,Norb))
		end do
	   end do
	end do
end do
!		debut=1+ (j-1)*nk
!		fin=nk + (j-1)*nk		
!		do i=1,orb_Q 
!			Hqp(ll,debut:fin)=Hqp(ll,debut:fin) + (muB(:,i) * Ers(ll,j+dimQ,i,Norb))
!		end do


!TODO error typing?
!  end do
!end do

Hpq=transpose(conjg(Hqp))
call gemm(Hqp,Hpq,Hqppq)
!Hqppq=Hqppq * dVk(dki) ! dVk is not needed anylonger



end subroutine Hnd_CSF

!****************************************************************************************
!****************************************************************************************
! Calcul H0qq,H0pp,MUqq,MUpp. Indépendant du temps
subroutine Hd_CSF(hA,muA,Vee,Ers,ee,H0qq,H0pp,MUqq,MUpp)
!****************************************************************************************
!****************************************************************************************
use variables

real(kind=real_8),dimension(:,:),intent(out)::H0qq,H0pp,MUqq,Mupp
complex(kind=comp_16),dimension(:,:),intent(in)::muA,hA
real(kind=real_8),dimension(:,:,:,:),intent(in)::Ers
complex(kind=comp_16),dimension(dimQ,dimQ)::Mtemp
integer(kind=int_4)::r,s,u,v
real(kind=real_8),dimension(:,:,:,:),intent(in)::Vee
real(kind=real_8),dimension(:,:,:,:,:,:),intent(in)::ee
real(kind=real_8),dimension(dimQ,dimQ)::test

write(*,*)'orb_Q',orb_Q

H0qq=0.d0
MUqq=0.d0
H0pp=0.d0
MUpp=0.d0

do s=1,orb_Q
!write(*,*)'In Hd_CSF, s =',s
do r=1,orb_Q
	H0qq=H0qq + (hA(r,s)  * Ers(1:dimQ,1:dimQ,r,s))
	MUqq=MUqq + (muA(r,s)  * Ers(1:dimQ,1:dimQ,r,s))
	H0pp=H0pp + (hA(r,s)  * Ers(dimQ+1:dimCSF,dimQ+1:dimCSF,r,s))
	MUpp=MUpp + (muA(r,s)  * Ers(dimQ+1:dimCSF,dimQ+1:dimCSF,r,s))
	do v=1,orb_Q
	do u=1,orb_Q
		H0qq= H0qq + (1.d0/2.d0)*(Vee(r,s,u,v)*ee(1:dimQ,1:dimQ,r,s,u,v))
		H0pp= H0pp + (1.d0/2.d0)*(Vee(r,s,u,v)*ee(dimQ+1:dimCSF,dimQ+1:dimCSF,r,s,u,v))
	end do
	end do
end do
end do



end subroutine Hd_CSF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Hnk(kk,n,zeta)
! A function similar to the hermite polynomials
! Very useful to write Hamiltonian matrix elements in a generalized form
! Stops at n = 5 but n > 5 can be calculated using H_n+1 = -(kk/2zeta)H_n + d/dk (H_n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use basics
implicit none
real(kind=real_8) :: Hnk
real(kind=real_8), intent(in) :: zeta,kk
integer, intent(in) :: n

! default value
Hnk = 1.d0

if(n.eq.1)then
Hnk = -1.d0*kk/(2.d0*zeta)
endif

if(n.eq.2)then
Hnk = (kk*kk-2.d0*zeta)/((2.d0*zeta)**2.d0)
endif

if(n.eq.3)then
Hnk = (-1.d0*(kk**3.d0) + 6.d0*zeta*kk )/((2.d0*zeta)**3.d0)
endif

if(n.eq.4)then
Hnk = (12.d0*zeta-12.d0*zeta*kk*kk+kk**4.d0)/((2.d0*zeta)**4.d0)
endif

if(n.eq.5)then
Hnk = -1.d0*(60.d0*(zeta**2.d0)*kk -20.d0*zeta*(kk**3.d0) + kk**5.d0)/((2.d0*zeta)**5.d0)
endif
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function dVk(dki)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use basics
implicit none
real(kind=real_8) 				:: dVk
real(kind=real_8), dimension(3),intent(in) 	:: dki
real(kind=real_8)				:: dkx, dky, dkz

dkx=dki(1)
dky=dki(2)
dkz=dki(3)

if (dki(1) .eq. 0.d0) then 
	dkx=1.d0
end if

if (dki(2) .eq. 0.d0) then 
	dky=1.d0
end if

if (dki(3) .eq. 0.d0) then 
	dkz=1.d0
end if

dVk=dkx*dky*dkz

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module qp_integrals
