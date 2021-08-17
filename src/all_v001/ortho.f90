module ortho

! orthoganisation de Gram-Schmidt, changement de base entre PW, OPW, OOPW

use Basics
use variables
use math
use blas95
use lapack95
use f95_precision
use math

contains
!****************************************************************************************
!****************************************************************************************
subroutine OPW2PW_3D(gamm,ski,skiCon,nki,dki)  !A optimiser
!Passage de OPW vers PW. 
! /!\ N'othogonalise pas les paquets d'ondes entre eux
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(inout)::gamm
complex(kind=comp_16),dimension(:,:,:,:),intent(in)::ski,skiCon
Integer(kind=int_4),Dimension(3),Intent(in)		::nki
Real(kind=real_8),Dimension(3),Intent(in)		::dki

complex(kind=comp_16),dimension(dimP,nk)::gamm2,gammCorrection
complex(kind=comp_16)::temp1,temp2
integer::kp1,kp2,kp3,k1,k2,k3,u,j,i,kkp,kk

gammCorrection=dcmplx(0.d0,0.d0)

do k1=1,nki(1)
do k2=1,nki(2)
do k3=1,nki(3)
	kk=k1+nki(1)*(k2-1)+nki(1)*nki(2)*(k3-1)
	do i=1,dimP
		temp2=dcmplx(0.d0,0.d0)
		do kp1=1,nki(1)
		do kp2=1,nki(2)
		do kp3=1,nki(3)
			kkp=kp1+nki(1)*(kp2-1)+nki(1)*nki(2)*(kp3-1)
			temp1=dcmplx(0.d0,0.d0)
			do j=1,orb_Q
				temp1=temp1 + (ski(k1,k2,k3,j)*skiCon(j,kp3,kp2,kp1)*gamm(i,kkp))
			end do
			temp2= temp2 +temp1
		end do
		end do
		end do
		gammCorrection(i,kk)= temp2 * dki(1)*dki(2)*dki(3)
		gamm2(i,kk)= gamm(i,kk) - gammCorrection(i,kk)
	end do
end do
end do
end do
gamm=0.d0

gamm=gamm2



end subroutine OPW2PW_3D

!****************************************************************************************
!****************************************************************************************
subroutine OPW2PW(gamm,ski,skiCon)  !A optimiser
!Passage de OPW vers PW. 
! /!\ N'othogonalise pas les paquets d'ondes entre eux
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(inout)::gamm
complex(kind=comp_16),dimension(:,:),intent(in)::ski,skiCon

complex(kind=comp_16),dimension(dimP,nk)::gamm2,gammCorrection
complex(kind=comp_16)::temp1,temp2
integer::k,j,i,kk,dkx

gammCorrection=dcmplx(0.d0,0.d0)
do k=1,nk
	do i=1,dimP
		temp2=dcmplx(0.d0,0.d0)
		do kk=1,nk
			temp1=dcmplx(0.d0,0.d0)
			do j=1,Norb_cat
				temp1=temp1 + (ski(k,j)*skiCon(j,kk)*gamm(i,kk))
			end do
			temp2= temp2 +temp1
		end do
		gammCorrection(i,k)= temp2 * dkx
		gamm2(i,k)= gamm(i,k) - gammCorrection(i,k)
	end do
end do
gamm=0.d0
gamm=gamm2



end subroutine OPW2PW

!****************************************************************************************
!****************************************************************************************
subroutine OPWtoOPW2(coeffGS,muB,hB)  
!Orthogonalisation de muB et hB
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(inout)::muB,hB
complex(kind=comp_16),dimension(:,:),intent(in)::coeffGS

complex(kind=comp_16),dimension(Norb_cat,nk)::muB2,hB2
integer::i

do i=1,Norb_cat
	call gemv(coeffGS,muB(i,:),muB2(i,:))
end do

muB=muB2

end subroutine OPWtoOPW2

!****************************************************************************************
!****************************************************************************************
! Coefficients nécessaire à la ré-orthogalisation
subroutine grammSchmidt(psipsi,coeffC)
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(in)::psipsi
complex(kind=comp_16),dimension(:,:),intent(out)::coeffC

complex(kind=comp_16),allocatable,dimension(:)::Norm
complex(kind=comp_16),allocatable,dimension(:,:)::coeffD
logical,allocatable,dimension(:,:)::coeffD_L
logical,allocatable,dimension(:)::Norm_L
complex(kind=comp_16)::sommeC
integer::k,kk,kk2,l,k2,l2,taille

taille=int(sqrt(real(size(coeffC))))
allocate(Norm(taille),coeffD(taille,taille),coeffD_L(taille,taille),Norm_L(taille))

coeffC=dcmplx(0.d0,0.d0)
coeffD=dcmplx(0.d0,0.d0)
Norm=dcmplx(0.d0,0.d0)
coeffD_L=.FALSE.
Norm_L=.FALSE.

Norm(1)=dcmplx(1.d0,0.d0)/zsqrt(psipsi(1,1))
Norm_L(1)=.TRUE.
coeffC(1,1)=Norm(1)
do k=1,taille
	coeffD(1,k)=Norm(1)*psipsi(1,k)
	coeffD_L(1,k)=.TRUE.
end do
do k=2,taille
  do l=1,taille
	k2=k
	l2=l
	call Norm_SUB(k2,coeffD,coeffD_L,Norm,Norm_L,psipsi)
	call coeffD_SUB(k2,l2,coeffD,coeffD_L,Norm,Norm_L,psipsi)
  end do
end do
do l=1,taille
	do k=1,l
		sommeC=dcmplx(0.d0,0.d0)
		do kk=1,l-1
			sommeC= sommeC + coeffC(kk,k) * coeffD(kk,l)
		end do
	    coeffC(l,k)=kronecker(l,k) - sommeC * Norm(l)
	end do
end do
end subroutine grammSchmidt

!****************************************************************************************
!****************************************************************************************
recursive subroutine coeffD_SUB(k,l,coeffD,coeffD_L,Norm,Norm_L,psipsi)
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(in)::psipsi
complex(kind=comp_16),dimension(:,:),intent(inout)::coeffD
complex(kind=comp_16),dimension(:),intent(inout)::Norm
logical,dimension(:,:),intent(inout)::coeffD_L
logical,dimension(:),intent(inout)::Norm_L
integer,intent(inout)::k,l

complex(kind=comp_16)::sommeD
integer::kk,kk2

if (coeffD_L(k,l) .ne. .TRUE.) then
	sommeD=0.d0
	do kk=1,k-1
		kk2=kk
		call coeffD_SUB(kk2,l,coeffD,coeffD_L,Norm,Norm_L,psipsi)
		call coeffD_SUB(kk2,l,coeffD,coeffD_L,Norm,Norm_L,psipsi)
		sommeD= sommeD + conjg(coeffD(kk2,k))*coeffD(kk2,l)
	end do
	if (norm_L(k).eq. .FALSE.) then
		write(*,*)'Erreur dans coeffD_SUB',k,l
	end if
	coeffD(k,l)=Norm(k) * (psipsi(k,l)-sommeD)
	coeffD_L(k,l)=.TRUE.
end if


end subroutine

!****************************************************************************************
!****************************************************************************************
subroutine norm_SUB(l,coeffD,coeffD_L,Norm,Norm_L,psipsi)
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(in)::psipsi
complex(kind=comp_16),dimension(:,:),intent(in)::coeffD
complex(kind=comp_16),dimension(:),intent(inout)::Norm
logical,dimension(:,:),intent(inout)::coeffD_L
logical,dimension(:),intent(inout)::Norm_L
integer,intent(in)::l

integer::k
complex(kind=comp_16)::sommeN

if (norm_L(l) .ne. .TRUE.) then
	sommeN=0.d0
	do k=1,l-1
		if (coeffD_L(k,l) .eq. .FALSE.) then
			write(*,*)'Erreur dans norm_sub',l,k
		end if
		sommeN= sommeN + cdabs(coeffD(k,l))**2.d0
	end do
	Norm(l)=(psipsi(l,l) - sommeN)**(-1.d0/2.d0)
	Norm_L(l)=.TRUE.
end if


end subroutine

end module ortho
