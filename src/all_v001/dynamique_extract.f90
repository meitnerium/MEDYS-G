Module dynamique
!Gère la dynamique dans sa première version
use variables
use Basics
use math
use champ1
use ortho
use grille
use qp_integrals
use blas95
use lapack95
use f95_precision
use NewSubs
use observable

contains

!****************************************************************************************
!****************************************************************************************
subroutine new_dyna(fctQ,fctP,Upp,muEPS,muEPS_Sq,sauvChamp,matA,matAlpha,matPhi,eps,cEta.cZeta,&
lmn_vec,prim_center,lcCG,Ers)
!Actual electron dynamics calculations:  Many-electron wavepacket in TDCI form, 
! Full Gaussian version 
!INPUT: 
!       Hamiltonian matrices H0qq, H0pp: core H in Q (P) subspace
!                            Hqp,Hpq: Q-P coupling H (time-independent part), Hqppq=Hqp*Hpq
!       initial CI content (fctQ, fctP),  
!OUTPUT: Final CI vectors (fctQ, fctP, sauvFctQ(?), sauvFctP(?)), TD field (sauvChamp). 
!****************************************************************************************
!****************************************************************************************
!Real(kind=real_8),Intent(in)				::theta,phy,delta
Complex(kind=comp_16),Dimension(:),Intent(inout)	::fctQ
Complex(kind=comp_16),Dimension(:,:,:),Intent(inout)	::fctP,Upp
Complex(kind=comp_16), Dimension(:,:),Intent(in) 	:: muEPS
Complex(kind=comp_16),Dimension(:,:) ,Intent(inout)		:: muEPS_Sq

Real(kind=real_8), Dimension(:,:),Intent(inout)		:: matA
 Real(kind=real_8),dimension(:,:),Intent(inout):: matAlpha, matPhi

Real(kind=real_8),dimension(:,:),Intent(in)	:: cEta,   lcCG, prim_center
Real(kind=real_8),dimension(:),Intent(inout)::  cZeta,eps, sauvChamp
!Integer(kind=int_4),allocatable,dimension(:,:),Intent(in)::prim_lmn
Integer(kind=int_4), dimension(:,3),Intent(in):: lmn_vec
Real(kind=real_8), Dimension(:,:,:,:),Intent(in) 	:: Ers

! Complex(kind=comp_16),Dimension(:,:),Intent(in)		::eigenVectH0pp

Complex(kind=comp_16),Dimension(dimQ,dimQ)		::Mat33,Inv33

Integer							::i,j,t,r,s,k,u,v,a,counter,ik1,ik2,ik3,ntk,ik
 

Complex(kind=comp_16),Dimension(dimQ,dimQ) 		::Hqppq 
!Real(kind=real_8),Dimension(:,:),Intent(in)			::H0qq,MUqqne
!Real(kind=real_8),Dimension(:,:),Intent(in)			::H0pp,MUpp
Real(kind=real_8),Dimension(3)	,Intent(in)			::eps


Real(kind=real_8)					::PreStart,PreEnd,TotalPre
Real(kind=real_8)					::PropStart,PropEnd,TotalProp
Real(kind=real_8)					::TotalCI,TotalGamma,TotalVolk
Real(kind=real_8),Dimension(3)				::sauvkmin
Real(kind=real_8)	        			::t1,t2
write(*,*)'Debut dynamique'

sauvChamp=0.d0
! 
sauvChamp(1)=champNint(delta) 

write(*,*)nt,'number of timesteps'
open(51,file='champ.dat',status='unknown',form='formatted')
!////////////////////////////////////////////////////////
! start loop over time variable (tn), i.e. TIME-PROPAGATION
!////////////////////////////////////////////////////////

call Volkov_OM(1,1,lcCG,matA,matAlpha,matPhi,eps,cEta,cZeta,lmn_vec,prim_center,muEPS,muEPS_Sq)
call Hqppq_CSF(muEPS_Sq,Ers,Hqppq) 

call cpu_time ( t1 )

do tn=1,nt  ! begin time-loop
t=tn
call cpu_time ( t2 )
write ( *, * ) 'tn = ',tn,'Elapsed CPU time = ', t2 - t1 
!
! CalPre: Preliminary steps (calculation of Inv33=(1_QQ +1/4 HBARqppq )**(-1), Mat33=1_QQ -1/4 HBARqppq
!
	call CalPre(Mat33,Inv33,Hqppq,delta)
!
! fill arrays matA,matAlpha,matPhi
!
        call fill_MatA_Alpha_phi(tn,matA,matAlpha,matPhi)
!
 write(*,*) 'calling propagation, tn=',tn,'/',nt
	call new_propagation(fctQ,fctP,Upp,muEPS,muEPS_Sq,H0qq,H0pp,MUqq,MUpp,Mat33,Inv33, & 
matA,matAlpha,matPhi,eps,cEta.cZeta,lmn_vec,prim_center,lcCG,Ers)

	sauvChamp(tn+1)=champNint(delta)
! 
	write(51,*)tmin+(tn*pdt),Int0(delta),sauvChamp(tn+1)

end subroutine new_dyna















!****************************************************************************************
!****************************************************************************************
! Calcul de Hqppq. Indépendant du temps (nouveau pour version MEDYS*)
subroutine Hqppq_CSF(muEPS_Sq,Ers,Hqppq) 
    !TODO nk est dans modules variables donc pas obligatoire de le passer ici  
!****************************************************************************************
!****************************************************************************************
use variables

 
complex(kind=comp_16),dimension(:,:),intent(out)::Hqppq
complex(kind=comp_16),dimension(:,:),intent(in)::muEPS_Sq
real(kind=real_8),dimension(:,:,:,:),intent(in)::Ers


integer(kind=int_4)::l,ll,i,j,jj,k,debut,fin,r,s
integer(kind=int_4)::kkx,kky,kkz,kk,jkk

Hqppq=dcmplx(0.d0,0.d0)

do l=1,dimQ
  do ll=1,dimQ		
	do j=1,dimP 
		do r=1,orb_Q
                  do s=1,orb_Q
			Hqppq(l ,ll) = Hqppq(l ,ll) + muEPS_Sq(r,s)*Ers(l ,j+dimQ,r,Norb))*Ers(ll ,j+dimQ,s,Norb))
		end do
	   end do
	end do
   end do
end do

end subroutine Hqppq_CSF 


!***************************************************
!***************************************************
subroutine fill_MatA_Alpha_phi(tn,matA,matAlpha,matPhi)  !!!   sketch! to be finished
!***************************************************
!***************************************************
Real(kind=real_8),  Dimension(:,:),Intent(out)		:: matA
 Real(kind=real_8), dimension(:,:),Intent(out):: matAlpha, matPhi
Real(kind=real_8),   	:: temp,temp2


temp=0.d0
temp2=0.d0

matA(tn,tn)=Int0()  ! matA(i,i) est l'aire du champ de t(i) a t(i+1)
%
if(tn.ge.2)then
do k=1,tn-1
matA(tn-k,tn)=matA(tn-k,tn-k) +matA(tn-k+1,tn) ! matA(i,n) est l'aire du champ de t(i) a t(n+1)  NOTE: avec n_max=nt, on a nt+1 temps (nt intervalles) 
enddo
endif
do k=1,tn
do i=k,tn 
temp=temp+matA(k,i) ! integrale de t(k) a t(tn) de l'aire A
temp2=temp2+matA(k,i)*matA(k,i)  ! integrale de t(k) a t(tn) de A**2
enddo
matAlpha(k,tn)=temp*pdt
matPhi(k,tn)=temp2*pdt
enddo

end subroutine fill_MatA_Alpha_phi


!****************************************************************************************
!****************************************************************************************
subroutine CalPre(Mat33,Inv33,Hqppq,delta)
!Calculs préliminaires -> génère les matrices auxiliaires X et Y
!ENTREE:Produit Hqppq
!SORTIES:X(Inv33) ,Y(Mat33)
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(out)	::Mat33,Inv33
complex(kind=comp_16),dimension(:,:),intent(in)		::Hqppq
complex(kind=comp_16),allocatable,dimension(:,:)	::temp3,temp4
complex(kind=comp_16),allocatable,dimension(:)		::work
complex(kind=comp_16),allocatable,dimension(:,:)	::work2
real(kind=real_8),intent(in)::delta
integer(kind=int_4)					::l,ll,i,ii,j,k,info,r,s,lwork
integer,allocatable,dimension(:)			::ipiv
allocate(ipiv(dimQ))

allocate(work(dimQ),temp3(dimQ,dimQ),temp4(dimQ,dimQ))


temp3=ZidentiMAT(dimQ) + ((1.d0/4.d0)*Hqppq*(Int0(delta)**2))
info=1
Inv33=dcmplx(0.d0,0.d0)
 Inv33=temp3


 CALL ZGETRF(dimQ,dimQ,Inv33,dimQ,ipiv,INFO)

IF (INFO.EQ.0) THEN
 CALL ZGETRI(dimQ,Inv33,dimQ,ipiv,WORK,-1,INFO)

lwork=dint(dreal(work(1)))
deallocate(work)
allocate(work(lwork))
write(*,*)" ZGETRI(dimQ,Inv33,dimQ,ipiv,WORK,lwork,INFO)"
write(*,*)" dimQ",dimQ
 CALL ZGETRI(dimQ,Inv33,dimQ,ipiv,WORK,lwork,INFO)

if (info.ne.0) then
write(*,*)'ERREUR INVERSION DE MATRICE DANS CalPre'
end if
else 
write(*,*)'ERREUR'
end if

 
Mat33=dcmplx(0.d0,0.d0)
Mat33=ZidentiMAT(dimQ) - ((1.d0/4.d0)*Hqppq* (Int0(delta)**2))
 

deallocate(work)
deallocate(temp3)
deallocate(ipiv)

end subroutine CalPre









!****************************************************************************************
!****************************************************************************************
subroutine  new_propagation(fctQ,fctP,Upp,muEPS,muEPS_Sq,H0qq,H0pp,MUqq,MUpp,Mat33,Inv33, &
matA,matAlpha,matPhi,eps,cEta,cZeta,lmn_vec,prim_center,lcCG,Ers)
!
!
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(in)	   ::Mat33,Inv33
real(kind=real_8),dimension(:,:),intent(in)	       ::H0qq,H0pp,MUqq,Mupp

!real(kind=real_8),intent(in)			           ::theta,phy,delta
!Integer(kind=int_4),Intent(in)		               ::opt_GS
complex(kind=comp_16),dimension(:),intent(inout)   ::fctQ
complex(kind=comp_16),dimension(:,:,:),intent(inout) ::fctP,Upp
Complex(kind=comp_16), Dimension(:,:),Intent(in) 	:: muEPS,muEPS_Sq

Real(kind=real_8),  Dimension(:,:),Intent(inout)		:: matA
 Real(kind=real_8), dimension(:,:),Intent(inout):: matAlpha, matPhi

Real(kind=real_8), dimension(:,:),Intent(in)	:: cEta,   lcCG, prim_center
Real(kind=real_8), dimension(:),Intent(in)::  cZeta,eps
!Integer(kind=int_4),allocatable,dimension(:,:),Intent(in)::prim_lmn
Integer(kind=int_4), dimension(:,3),Intent(in):: lmn_vec
Real(kind=real_8),  Dimension(:,:,:,:),Intent(in) 	:: Ers

complex(kind=comp_16),allocatable,dimension(:)     ::CI1,CI2,CI3,CI4,CI5,CI
complex(kind=comp_16),allocatable,dimension(:,:)   ::muEPS_r_gam,Uqq  , gamm2!,gammCorrection
!complex(kind=comp_16),allocatable,dimension(:,:,:)   ::Upp
complex(kind=comp_16),allocatable,dimension(:,:,:,:)	   ::Volkov_OM_mat

integer						                   ::l,ll,i,j,jj,k,r,s,k1,k2,k3,kk,kkj
real(kind=real_8)				                   ::Et,normQ,normP,P_normfactor

real(kind=real_8),allocatable,dimension(:)	       ::fctPInt,fctPInt2

real(kind=real_8)                                  ::seuil !TODO to be deleted after test

    Character(len=5)                    ::nom
!    Character(len=23)                   ::canal
!    CHARACTER(LEN=20) :: FMT

allocate(Uqq(dimQ,dimQ) )!,Upp(nt,dimP,dimP))
Uqq=dcmplx(0.d0,0.d0)
Upp=dcmplx(0.d0,0.d0)
call Udiag(H0qq,MUqq,dimQ,Uqq,delta)
call Udiag(H0pp,MUpp,dimP,Upp(tn,:,:),delta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(Volkov_OM_mat(nt,nt,orb_Q,orb_Q), muEPS_r_gam(orb_Q,dimP))
allocate(CI1(dimQ),CI2(dimQ),CI3(dimQ),CI4(dimQ),CI5(dimQ),CI(dimQ))

!
!!!!!!Propagation within Q-space!!!  
!
muEPS_r_gam=dcmplx(0.d0,0.d0)
Volkov_OM_mat=dcmplx(0.d0,0.d0)


call Get_muEPS_r_gam(Upp,muEPS,muEPS_r_gam,lcCG,matA,matAlpha,matPhi,eps,cEta,cZeta,lmn_vec,prim_center,Volkov_OM_mat)

 CI1=dcmplx(0.d0,0.d0); CI2=dcmplx(0.d0,0.d0); CI3=dcmplx(0.d0,0.d0)
 CI4=dcmplx(0.d0,0.d0); CI5=dcmplx(0.d0,0.d0); CI=dcmplx(0.d0,0.d0)

!call gemv(HqpBar,BigGamma,CI1)   		!Etape1 !; CI1=CI1*dVk(dki) nolonger needed

do i=1,dimQ	
	do j=1,dimP 
		do r=1,orb_Q
			CI1(i) = CI1(i) + muEPS_r_gam(r,j)*Ers(i ,j+dimQ,r,Norb)) !on a donc besoin de muEPS_r_gam(orb_Q,dimP) (est-ce utile de mettre tn?)
		end do
	end do
end do
call gemv(Inv33,CI1,CI2)  				!Etape2 
call gemv(Inv33,fctQ,CI3) 				!Etape3 
call gemv(Mat33,CI3,CI4) 				!Etape4
 CI5=CI4 - (dcmplx(0.d0,1.d0) * CI2) 	!Etape5
call gemv(Uqq,CI5,CI)                   !Etape6
 
fctQ=CI					!propagation of Q-part of wp done

deallocate(CI4,CI5,CI)
!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!Propagation within P-space!!!
!
allocate(gamm2(dimP,orb_Q))
 CI1=dcmplx(0.d0,0.d0)
 CI1=dcmplx(0.d0,-1.d0/2.d0)*CI2 + CI3 		!nouvelle Etape10

do j=1,dimP 
		do r=1,orb_Q
                   do i=1,dimQ
			gamm2(j,r)=dcomplx(0.d0, -1.d0)*Conjg(Ers(i ,j+dimQ,r,Norb))*CI1(i) !nouvelles Etape11-12
                   enddo
		end do
	end do
call gemm(Upp(tn,:,:),gamm2,fctP(:,:,tn+1)) !nouvelle Etape13(a)

call Get_Observable_P_I(Upp,fctP,muEPS_Sq,Volkov_OM_mat,matA,P_J_Ion_amplitude,tn+1)

!    open(52,file='PSomme.dat',status='replace',form='formatted')
     do j=1,dimP
        write (nom,'(I5.5)') j
        open(90+j,file="Proba_ionisation_canal"// ADJUSTL(nom) //".dat",status='replace',form='formatted')
!        open(110+j,file="spectre2dx_canal"// ADJUSTL(nom) //".dat",status='replace',form='formatted')

    end do

        do j=1,dimP
          write(90+j,'(3X,I5.2,3X,F10.5,3X,ES24.14)')  tn+1, (tn+1)*delta, cdabs(P_J_Ion_amplitude(j))**2
        end do
end subroutine new_propagation

!****************************************************************************************
!****************************************************************************************
subroutine Udiag(H0,MU,lda,matrice,delta)
!Opérateurs diagonaux Uqq et Upp
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),intent(in)::lda
real(kind=real_8),dimension(:,:),intent(in)::H0,MU
real(kind=real_8),intent(in)::delta
complex(kind=comp_16),allocatable,dimension(:,:)::tempMAT,inversH,tempMAT2
complex(kind=comp_16),dimension(:,:),intent(out)::matrice
complex(kind=comp_16), allocatable,dimension(:) :: eval
complex(kind=comp_16), allocatable,dimension(:,:) :: evec
integer::i,l,info

allocate(tempMAT(lda,lda),tempMAT2(lda,lda),inversH(lda,lda),eval(lda),evec(lda,lda))
matrice=0.d0

tempMat=H0*pdt + Int0(delta)*MU
 call Diagonalize(tempMat,lda,evec,eval)
inversH=transpose(conjg(evec))
!write(*,*)'U0 Diag',lda
!do l=1,lda
!write(*,'(e12.5,3X,8e12.5)')real(eval(l)),(real(evec(i,l)),i=1,lda)
!end do

tempMat=dcmplx(0.d0,0.d0)
do l=1,lda
	tempMat(l,l)=exp(dcmplx(0.d0,-1.d0)*eval(l))
end do
call gemm(tempMat,inversH,tempMat2)
call gemm(evec,tempMat2,matrice)


end subroutine Udiag


!****************************************************************************************
!****************************************************************************************
subroutine calculEE(Ers,ee)
!Générateurs bi-électroniques
!****************************************************************************************
!****************************************************************************************
real(kind=real_8),dimension(:,:,:,:)::Ers
real(kind=real_8),dimension(:,:,:,:,:,:)::ee
integer::r,s,u,v

ee=0.d0
do v=1,norb
do u=1,norb
do s=1,norb
do r=1,norb
	!calculEE(:,:,r,s,u,v)=0.d0
	!call DGEMM('N','N',dimCSF,dimCSF,dimCSF,1.d0,Ers(:,:,r,s),dimCSF,Ers(:,:,u,v),dimCSF,0.d0,calculEE(:,:,r,s,u,v),dimCSF)
	call gemm(Ers(:,:,r,s),Ers(:,:,u,v),ee(:,:,r,s,u,v))
	if (s.eq.u) then
		ee(:,:,r,s,u,v)=ee(:,:,r,s,u,v) - Ers(:,:,r,v)
	end if
end do
end do
end do
end do

end subroutine calculEE

!****************************************************************************************
!****************************************************************************************
subroutine etatInitial(H0,lda,EI)
!****************************************************************************************
!****************************************************************************************
integer(kind=int_8),intent(in)			::lda
real(kind=real_8),dimension(:,:),intent(in)	::H0
complex(kind=comp_16),allocatable,dimension(:,:)::tempMAT,inversH,workk,tempMAT2
complex(kind=comp_16),dimension(:),intent(out)	::EI
complex(kind=comp_16), allocatable,dimension(:)	:: eval
complex(kind=comp_16), allocatable,dimension(:,:) :: evec
integer,allocatable,dimension(:)		::ipivv,isuppz
integer::l,info,i,j

allocate(tempMAT(lda,lda),tempMAT2(lda,lda),inversH(lda,lda),ipivv(lda),workk(lda,lda),isuppz(lda),eval(lda),evec(lda,lda))

tempMat=H0
EI=dcmplx(0.d0,0.d0)
evec=0.d0
eval=0.d0
!write(*,*)'H0',tempMat
 call Diagonalize(tempMat,lda,evec,eval)

Write(18,*)
Do I=1,lda

	Write (18,'(A,I3)')'Eigenvalue of the eigenvector #',I
	
	Write (18, '(F22.14)') Real(eval(I))
	
	if (dabs(aimag(eval(I))) .gt. 1.d-5) then
		Write(*,*)'eigenValue',I,'has a complex part'
		STOP 'error in etatInitial'
	endif


	Write (18,'(A,I3)') 'Eigenvector #',I
	Do J=1,lda
	        Write (18,'(2(F22.14))') (evec(j,i))
	End do

	Write (18,*)
End Do
EI=-evec(:,1)

end subroutine etatInitial

!****************************************************************************************
!****************************************************************************************
subroutine EigenVectors(H0,lda,EI)
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),intent(in)			::lda
real(kind=real_8),dimension(:,:),intent(in)	::H0
complex(kind=comp_16),allocatable,dimension(:,:)::tempMAT,inversH,workk,tempMAT2
complex(kind=comp_16),dimension(:,:),intent(out)::EI
complex(kind=comp_16), allocatable,dimension(:)	:: eval
complex(kind=comp_16), allocatable,dimension(:,:) :: evec
integer,allocatable,dimension(:)		::ipivv,isuppz
integer::l,info,i,j

allocate(tempMAT(lda,lda),tempMAT2(lda,lda),inversH(lda,lda),ipivv(lda),workk(lda,lda),isuppz(lda),eval(lda),evec(lda,lda))

tempMat=H0
EI=dcmplx(0.d0,0.d0)
evec=0.d0
eval=0.d0
!write(*,*)'H0',tempMat
 call Diagonalize(tempMat,lda,evec,eval)

Write(18,*)
Do I=1,lda

	Write (18,'(A,I3)')'Eigenvalue of the eigenvector #',I
	
	Write (18, '(F22.14)') Real(eval(I))
	
	if (dabs(aimag(eval(I))) .gt. 1.d-5) then
		Write(*,*)'eigenValue',I,'has a complex part'
		STOP 'error in etatInitial'
	endif


	Write (18,'(A,I3)') 'Eigenvector #',I
	Do J=1,lda
	        Write (18,'(2(F22.14))') (evec(j,i))
	End do

	Write (18,*)
End Do
EI=-evec

end subroutine EigenVectors

!!****************************************************************************************
!!****************************************************************************************
!subroutine CSF2States(FctP,eigenVectH0pp)
!!****************************************************************************************
!!****************************************************************************************
!Complex(kind=comp_16),Dimension(:,:),Intent(inout) 	::FctP
!Complex(kind=comp_16),Dimension(:,:),Intent(in) 	::eigenVectH0pp
!
!Complex(kind=comp_16),Allocatable,Dimension(:,:) 	::temp_gamma
!Integer(kind=int_4)					::sizeP,sizek,j
!
!sizeP=size(FctP(:,1)) 
!sizek=size(FctP(1,:)) 
!
!allocate(temp_gamma(sizeP,sizek))
!
!call gemm(eigenVectH0pp,FctP,temp_gamma)
!
!FctP = temp_gamma
!
!end subroutine CSF2States

!*******************************************
!*******************************************
subroutine Get_muEPS_r_gam(Upp,muEPS,muEPS_r_gam,lcCG,matA,matAlpha,matPhi,eps,cEta,cZeta,lmn_vec,prim_center,Volkov_OM_mat)
!*******************************************
!*******************************************
complex(kind=comp_16),dimension(:,:,:),intent(in) ::fctP,Upp
Complex(kind=comp_16), Dimension(:,:),Intent(in) 	:: muEPS

Real(kind=real_8),  Dimension(:,:),Intent(inout)		:: matA
 Real(kind=real_8), dimension(:,:),Intent(inout):: matAlpha, matPhi

Real(kind=real_8), dimension(:,:),Intent(in)	:: cEta,   lcCG, prim_center
Real(kind=real_8), dimension(:),Intent(in)::  cZeta, eps
!Integer(kind=int_4),allocatable,dimension(:,:),Intent(in)::prim_lmn
Integer(kind=int_4), dimension(:,3),Intent(in):: lmn_vec
Real(kind=real_8),  Dimension(:,:,:,:),Intent(in) 	:: Ers

 complex(kind=comp_16),allocatable,dimension(:,:)     ::, gamm2 !CI1,CI2,CI3,CI4,CI5,CI
complex(kind=comp_16), dimension(:,:),Intent(inout)   ::muEPS_r_gam !,gammCorrection
complex(kind=comp_16), dimension(:,:,:,:)	   ::Volkov_OM_mat
complex(kind=comp_16), dimension(dimP)	   ::temp
integer						                       :: nn,j,r,s,tnn

 allocate{gamm2(orb_Q,orb_Q)
if(tn.ge.2)then
do j=1,dimP 
do r=1,orb_Q
 do s=1,orb_Q
  do tnn=2,tn 
    call Volkov_OM(tn,tnn-1,lcCG,matA,matAlpha,matPhi,eps,cEta,cZeta,lmn_vec,prim_center,muEPS,gamm2) 
      Volkov_OM_mat(tn,tnn-1,:,:)=gamm2
    temp=fctP(:,s,tnn)  
      if(tnn.lt.tn)then
      do nn=tnn,tn-1
       call gemv(Upp(nn,:,:),temp,temp)
      end do
     endif
  muEPS_r_gam(r,j)=muEPS_r_gam(r,j)+temp(j)*matA(tnn,tnn)*Volkov_OM_mat(tn,tnn-1,r,s)
  enddo
 end do
end do
end do
muEPS_r_gam=matA(tn,tn)*muEPS_r_gam
endif

end subroutine Get_muEPS_r_gam

!***************************************************

!***************************************************
!***************************************************
subroutine Volkov_OM(tn,tnn,lcCG,matA,matAlpha,matPhi,eps,cEta,cZeta,lmn_vec,prim_center,muEPS,Volkov_OM_mat_res)
!***************************************************
!***************************************************
Real(kind=real_8),  Dimension(:,:),Intent(inout)		:: matA
 Real(kind=real_8), dimension(:,:),Intent(inout):: matAlpha, matPhi

Real(kind=real_8), dimension(:,:),Intent(in)	:: cEta,   lcCG, prim_center
Real(kind=real_8), dimension(:),Intent(in)::  cZeta,eps
!Integer(kind=int_4),allocatable,dimension(:,:),Intent(in)::prim_lmn
Integer(kind=int_4), dimension(:,3),Intent(in):: lmn_vec
complex(kind=comp_16), dimension(:,:),Intent(out)	   ::Volkov_OM_mat_res
complex(kind=comp_16),allocatable,dimension(:,:,:,:) 	   ::GMutG 
Complex(kind=comp_16), Dimension(:,:),Intent(in) 	:: muEPS
Complex(kind=comp_16), Dimension(: ) 	:: temp
Complex(kind=comp_16), Dimension(:,: ) 	:: temp,temp2!,GMutG 

Integer(kind=int_4),     Intent(in)::tn,tnn
integer						                       :: tn,j,r,s,tnn,ss

allocate(temp(orb_Q,totPrimCount), temp2(orb_Q,orb_Q),Volkov_OM_mat_res(orb_Q,orb_Q),GMutG(totPrimCount,totPrimCount,1,1))

!***************
call Get_GMutG(tn,tnn,matA,matAlpha,matPhi,eps,cEta,cZeta,lmn_vec,prim_center,GMutG)

temp=matmul(lcCG,GMutG(:,:,1,1)) !! for <G1|d (U_volkov) d |G2>
temp2=matmul(temp,transpose(lcCG)) 
temp=matmul(lcCG,GMutG(:,:,0,1))  !! for <G1|d^0 (U_volkov) d |G2>
temp3=matmul(temp,transpose(lcCG))
     temp=matmul(muEPS,temp3)
     Volkov_OM_mat_res=temp2-temp-conjg(transpose(temp))
temp=matmul(lcCG,GMutG(:,:,0,0)) !! for <G1|d^0 (U_volkov) d^0  |G2>
temp3=matmul(temp,transpose(lcCG)) 
temp2=matmul(muEPS,temp3)
     Volkov_OM_mat_res=Volkov_OM_mat_res+ matmul(temp2, muEPS)

end subroutine Volkov_OM
!***************************************************
!***************************************************
subroutine Get_GMutG(tn,tnn,matA,matAlpha,matPhi,eps,cEta.cZeta,lmn_vec,prim_center,GMutG)
!***************************************************
!***************************************************
Real(kind=real_8),  Dimension(:,:),Intent(inout)		:: matA
 Real(kind=real_8), dimension(:,:),Intent(inout):: matAlpha, matPhi
Real(kind=real_8)    		:: alpha,Anl

Real(kind=real_8), dimension(:,:),Intent(in)	:: cEta,   lcCG, prim_center
Real(kind=real_8), dimension(:),Intent(in)::  cZeta,eps
!Integer(kind=int_4),allocatable,dimension(:,:),Intent(in)::prim_lmn
Integer(kind=int_4), dimension(:,3),Intent(in):: lmn_vect
complex(kind=comp_16), dimension(:,:,:,:),Intent(out)	   ::Volkov_OM_mat
complex(kind=comp_16), dimension(:,:,:,:) 	   ::GMutG 
Complex(kind=comp_16), Dimension(:,:),Intent(in) 	:: muEPS
Complex(kind=comp_16), Dimension(: ) 	:: temp
Complex(kind=comp_16), Dimension(:,: ) 	:: temp,temp2!,GMutG 
Complex(kind=comp_16), Dimension(3 ) 	::gmu1dg00,gmu1dg01,gmu1dg10,gmu1dg11
Integer(kind=int_4), Intent(in)::tn,tnn,totPrimCount 
Integer(kind=int_4)           :: i,j,k,r,s,kp,kpp,delt

 
 delt=(tn-tnn) 
GMutG=dcmplx(0.d0,0.d0)
Anl=0.d0
alpha=0.d0
do r=1,totPrimCount
   do s=1,totPrimCount
 
      do i=1,3
         if(delt.ne.0)then
            Anl=matA(tnn,tn)*eps(i)
            alpha=matAlpha(tnn,tn)*eps(i)
         endif
         call Get_GMut1DG(cZeta(r),cZeta(s),prim_center(r,i),prim_center(s,i),lmn_vect(r,i),lmn_vect(s,i),0,0,delt,alpha,Anl,gmu1dg00(i))
         call Get_GMut1DG(cZeta(r),cZeta(s),prim_center(r,i),prim_center(s,i),lmn_vect(r,i),lmn_vect(s,i),0,1,delt,alpha,Anl,gmu1dg01(i))
         call Get_GMut1DG(cZeta(r),cZeta(s),prim_center(r,i),prim_center(s,i),lmn_vect(r,i),lmn_vect(s,i),1,0,delt,alpha,Anl,gmu1dg10(i))
         call Get_GMut1DG(cZeta(r),cZeta(s),prim_center(r,i),prim_center(s,i),lmn_vect(r,i),lmn_vect(s,i),1,1,delt,alpha,Anl,gmu1dg11(i))
      enddo
      do k=1,3
         select case(k)
           case(1)
           kp=2
           kpp=3
           case(2)
           kp=3
           kpp=1
           case(3)
           kp=1
           kpp=2
        end select case
!
      GMutG(r,s,1,1)=GMutG(r,s,1,1)+eps(kp)*eps(kpp)*( gmu1dg10(kp)*gmu1dg01(kpp)+gmu1dg10(kpp)*gmu1dg01(kp) )*gmu1dg00(k)
      GMutG(r,s,1,1)=GMutG(r,s,1,1)+eps(k)*eps(k)*gmu1dg11(k)*gmu1dg00(kp)*gmu1dg00(kpp)
!
      GMutG(r,s,0,1)=GMutG(r,s,0,1)+eps(k)*gmu1dg01(k)*gmu1dg00(kp)*gmu1dg00(kpp)
    enddo
      GMutG(r,s,0,0)= gmu1dg00(1)*gmu1dg00(2)*gmu1dg00(3)
enddo
GMutG=GMutG*cdexp(0.d0,-matPhi(tnn,tn) )

end subroutine Get_GMutG









!***************************************************
!***************************************************
!***************************************************
subroutine Get_GMut1DG(zet1,zet2,R1,R2,l1,l2,v1,v2,dlt,alpha,Anl,res)
!***************************************************
!***************************************************
Real(kind=real_8),Intent(in)    		:: alpha,Anl,zet1,zet2,R1,R2
Integer(kind=int_4) ::l1,l2,v1,v2,lv1,lv2,llvv,l12,dlt
Complex(kind=comp_16) 	:: res,rho 

Complex(kind=comp_16)  	:: temp,temp2,prefact,capB,zetC,Im_u
Real(kind=real_8)    		:: z1Mod,z1Phase,z2Mod,z2Phase,zCMod,zCPhase,a
integer						                       :: k,kp,ksum, k1,k2,k3,k4,p,pp

Im_u=dcmplx(0.d0, 1.d0)	!Def. found in Module basics : Real(kind=real_8),parameter :: pio2=pi/2.0_real_8

 

prefact=sqrt(sqrt(4*zet1*zet2/pi**2 ))*cdexp(0.d0,-(R2+alpha)*Anl)
prefact=sqrt(pi)*sqrt(((4*zet1)**l1)*((4*zet2)**l2))*prefact
select case(dlt)
case(:0)
a=0.d0
zetC=dcmplx(zet1+zet2,0.d0)
rho=dcmplx(1.d0,0.0)
pmax=0
ppmax=0
case default
a=1/(2.d0*dlt*pdt)
zetC=dcmplx(zet2,a)
z2Mod=abs(zetC)
z2Phase=-ATAN2(zetC)
rho=Im_u*a*zetC/z2Mod**2
zetC=  -Im_u*a*(1.d0+rho)+ zet1 ! - Im_u*a+ (a**2)zetC/z2Mod**2+zet1
pmax=l2
ppmax=v2
prefact=sqrt(a)*cdexp(0.d0, -pi/4)*prefact 
end select case


zCMod=abs(zetC)
zCPhase=ATAN2(zetC)
capB=-dcmplx(zet2*(R2-R1 +alpha) , -Anl/2.d0)*cdexp(0,d0,-zCPhase)/zCMod
prefact=cdexp(zetC*capB*capB)*prefact
prefact=cdexp(-zet1*(alpha + R2-R1)*(alpha + R2-R1), 0.0d0)*prefact


temp=dcmplx(1.d0,0.d0)
res=dcmplx(0.d0,0.0)


do k=0,pmax  

   do kp=0, ppmax
      p=k+kp
      if((dlt.ne.0).and.(Mod(p,2).eq.0)) then
         temp=comb(l2,k)*comb(v2,kp)*fact( p-1,2)/(sqrt(2)**p)*(sqrt(z2Mod)**(-(p+1))*cdexp(0.d0,-0.5d0*z2Phase*(p+1))
      endif
temp2=0.d0
      do k1= 0,l1
         do k2=0,v1
            do k3=0,pmax-k
               do k4=0,ppmax-kp
                 ksum=k1+k2+k3+k4
                 if(Mod(ksum,2).eq.0) then
                   temp2=temp2+comb(l1,k1)*comb(v1,k2)*comb(pmax-k,k3)*comb(ppmax-kp,k4)*(rho**(pmax-k+k4))*((capB+R2-R1+alpha)**(l1-k1))*((capB+R2+alpha)**(v1-k2))*(capB**(pmax-k-k3))*((R2-rho*capB)**(ppmax-kp-k4))*(fact(ksum-1,2)/(sqrt(2)**ksum) )*(sqrt(zCMod))**(-(ksum+1))*cdexp(0.d0,-0.5d0*zCPhase*(ksum+1))
                 endif
              enddo
            enddo
          enddo
       enddo
       res=res+temp*temp2
     enddo
  enddo
enddo
end subroutine Get_GMut1DG
!***************************************************
!***************************************************
function comb(n,k)
!***************************************************
!***************************************************
implicit none
integer, intent (in) :: n
integer, intent (in) :: k
integer :: comb
 
    comb = fact(n,1) / (fact(k,1) * fact(n - k,1))

end function


!***************************************************
!***************************************************
!***************************************************
subroutine Get_Observable_P_I(Upp,fctP,muEPS_Sq,Volkov_OM_mat,matA,P_J_Ion_amplitude,t_idx)
!***************************************************
!***************************************************
complex(kind=comp_16),dimension(:,:,:),intent(in) ::fctP,Upp
Real(kind=real_8),  Dimension(:,:),Intent(in)		:: matA
complex(kind=comp_16),allocatable,dimension(:,:,:,:),Intent(in)	   ::Volkov_OM_mat
Complex(kind=comp_16), Dimension(:,:),Intent(in) 	::  muEPS_Sq
Complex(kind=comp_16), Dimension(:) 	::  P_J_Ion_amplitude
Integer(kind=int_4), Intent(in)::t_idx
Integer(kind=int_4),  ::nn,r,s,tnn,tnp
Complex(kind=comp_16), Dimension(dimP) 	::  temp,temp2
Complex(kind=comp_16),  	::  ztemp 

P_J_Ion_amplitude=dcmplx(0.d0,0.d0)
if(t_idx.ge.2) then
  do j=1,dimP
     do r=1,orb_Q
        do s=1,orb_Q
           do tnn=2, t_idx-1
              do tnp=1,tnn-1
                if(tnn.lt.t_idx) then
                   temp=fctP(:,r,tnn)
                   temp2=fctP(:,s,tnp)
                   do nn=tnn,tn-1
                      call gemv(Upp(nn,:,:),temp,temp)
                   enddo
                   do nn=tnp,tn-1
                      call gemv(Upp(nn,:,:),temp2,temp2)
                   enddo
                 endif
               ztemp=matA(tnp,tnp)*matA(tnn,tnn)*conjg(temp(j))*temp2(j)*Volkov_OM_mat(tnn,tnp,r,s)  !  matA(tnp-1,tnp-1)*matA(tnn-1,tnn-1)  ????
               P_J_Ion_amplitude(j)= P_J_Ion_amplitude(j)+(ztemp+conjg(ztemp))
               enddo
               P_J_Ion_amplitude(j)= P_J_Ion_amplitude(j)+matA(tnn,tnn)*matA(tnn,tnn)*conjg(temp(j))*temp2(j)*muEPS_Sq(r,s)
            enddo
         enddo
      enddo
 enddo

 end subroutine Get_Observable_P_I        


!*******************************************
!*******************************************
subroutine Get_Observable_eMomentum_map(kvec,fctP,Upp,ek_Amp,lcCG,muEPS,matA,matAlpha,matPhi,eps,cEta,cZeta,lmn_vec,prim_center)
!*******************************************
!*******************************************
complex(kind=comp_16),dimension(:,:,:),intent(in) ::fctP,Upp
Complex(kind=comp_16), Dimension(:,:),Intent(in) 	:: muEPS

Real(kind=real_8),  Dimension(:,:),Intent(in)		:: matA
 Real(kind=real_8), dimension(:,:),Intent(in):: matAlpha, matPhi

Real(kind=real_8), dimension(:,:),Intent(in)	:: cEta,   lcCG, prim_center
Real(kind=real_8), dimension(:),Intent(in)::  cZeta, eps,kvec
!Integer(kind=int_4),allocatable,dimension(:,:),Intent(in)::prim_lmn
Integer(kind=int_4), dimension(:,3),Intent(in):: lmn_vec
!Real(kind=real_8),  Dimension(:,:,:,:),Intent(in) 	:: Ers

 complex(kind=comp_16),allocatable,dimension(:)     ::, gamm2 !CI1,CI2,CI3,CI4,CI5,CI
complex(kind=comp_16), dimension(dimP),Intent(out)   ::ek_Amp  
complex(kind=comp_16), dimension(:,:,:,:)	   ::Volkov_OM_mat
complex(kind=comp_16), dimension(dimP,orb_Q)	   ::temp
integer						                       :: nn,j,r,s,tnn !,totPrimCount

 allocate{gamm2(orb_Q)
  
ek_Amp=dcmplx(0.d0,0.d0)

if(tn.ge.2)then

  do tnn=2,tn 
    call emomentum_from_MO(tn-tnn,matA(tnn,tn),matAlpha(tnn,tn),matPhi(tnn,tn),eps,cZeta ,lmn_vec ,prim_center,lcCG,muEPS,kvec,gamm2) 
    temp=fctP(:,:,tnn)  
      if(tnn.lt.tn)then
      do nn=tnn,tn-1
       call gemv(Upp(nn,:,:),temp,temp)
     endif
    end do
  do j=1,dimP 
  ek_Amp(j)=ek_Amp(j)+matA(tnn,tnn)*DOT_PRODUCT(temp(j,:),gamm2)
  enddo
 
endif

end subroutine Get_Observable_eMomentum_map

!*******************************************
!*******************************************
subroutine emomentum_from_MO(it,Anl, Alpha_nl, Phi_nl,eps,cZeta ,lmn_vec ,prim_center,lcCG,muEPS,kvec,gamm_res) 

!*******************************************
!*******************************************
 
Complex(kind=comp_16), Dimension(:,:),Intent(in) 	:: muEPS
Real(kind=real_8),   Intent(in)		:: Anl, Alpha_nl,Phi_nl

Real(kind=real_8), dimension(:,:),Intent(in)	:: cEta,   lcCG, prim_center
Real(kind=real_8), dimension(:),Intent(in)::  cZeta, eps,kvec
 
Integer(kind=int_4), dimension(:,3),Intent(in):: lmn_vec
Integer(kind=int_4), dimension(3):: lmn_p
Integer(kind=int_4),  Intent(in):: it 
 complex(kind=comp_16), dimension(totPrimCount),     ::BigGamma1, BigGamma2

 complex(kind=comp_16), dimension(dimP), Intent(out)     ::  gamm_res 

complex(kind=comp_16), dimension(dimP)	   ::temp1,temp2, temp
integer						                       :: nn,i_prim,j,r,s,tnn 


BigGamma1=dcmplx(0.d0,0.d0)
BigGamma2=dcmplx(0.d0,0.d0)

do i_prim=1,totPrimCount

lmn_p=lmn_vec(i)
    do j=1,3
       select case(j)
       case(1)
        lmn_p(1)=lmn_p(1)+1
       case(2)
        lmn_p(2)=lmn_p(2)+1
       case(3)
        lmn_p(3)=lmn_p(3)+1
       end select case
       
       call uvolkov_on_CG(Anl,Alpha_nl, Phi_nl,cZeta(i_prim),lmn_p,prim_center(i_prim,:), eps,kvec,it, BigGamma1(i_prim))
       BigGamma1(i_prim)=eps(j)*BigGamma1(i_prim)/dsqrt(4.d0*cZeta(i_prim))
    enddo

    call uvolkov_on_CG(Anl,Alpha_nl, Phi_nl,cZeta(i_prim),lmn_vec(i_prim),prim_center(i_prim,:), eps,kvec,it, BigGamma2(i_prim))

    BigGamma1(i_prim)=BigGamma1(i_prim)-dot_product(eps,prim_center(i_prim,:) )*BigGamma2(i_prim)

enddo

call gemm(lcCG, BigGamma1,temp1)
call gemm(lcCG, BigGamma2,temp2)
call gemm(temp2,muEPS,temp)

gamm_res=temp1-temp

end subroutine emomentum_from_MO

!*******************************************
!*******************************************
subroutine uvolkov_on_CG(Anl,Alpha_nl, Phi_nl,zet,lmn,Rc, eps, veck,Dtime, res)

!*******************************************
!*******************************************
use basics
implicit none
!Complex(kind=comp_16), Dimension(:,:),Intent(in) 	:: muEPS
Real(kind=real_8),   Intent(in)		:: Anl, Alpha_nl,Phi_nl,zet

Real(kind=real_8), dimension(:,:),Intent(in)	:: cEta,   lcCG, prim_center
Real(kind=real_8), dimension(3),Intent(in)::   eps,veck,Rc
Real(kind=real_8),  ::   k_transl, R_transl,Atemp, Alpha_temp
Integer(kind=int_4), dimension(3),Intent(in):: lmn 
Integer(kind=int_4),  Intent(in):: Dtime 
 complex(kind=comp_16),  Intent(out)    ::res

integer						                       :: n,nn,i,j  
nn=0
res=dcmplx(1.d0,0.d0)
do i=1,3 
Atemp=Anl*eps(i)
Alpha_temp=Alpha_nl*eps(i)
k_transl=veck(i)-Atemp
R_transl=Rc(i)+Alpha_temp

n=lmn(i)

res= res*Hermite(k_transl/dsqrt(4.d0*zet),n )*cdexp( -k_transl*k_transl*( dcmplx(1.d0/(4.d0*zet), -Dtime*pdt/2.d0) )
res=res*cdexp(0.d0, R_transl*k_transl)
nn=nn+n
enddo

res=res*cdexp(0.d0,-pio2*nn-Phi_nl)/(dsqrt(dsqrt(twopi*zet)))**3  ! e^(-i pi/2) =-i

end subroutine uvolkov_on_CG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Hermite(kk,n)
! Hermite polynomials

! Stops at n = 10 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use basics
implicit none
real(kind=real_8) :: Hermite
real(kind=real_8), intent(in) :: kk
real(kind=real_8), dimension(10)   :: A
integer, intent(in) :: n
integer,   :: i

! default value
 
     Hermite=0.d0
     call Hermite_Coeff(n,A)

    if(n.le.10)then
     do i = 0, n  
        Hermite=Hermite+ A(i)*(kk**i)  
     end do 
    else
      print(*,*) `error: n too large' 
    endif
   stop
end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*******************************************
!*******************************************
Subroutine Hermite_Coeff(n,A) 
!*******************************************
!******************************************* 
  integer i,j,n
real(kind=real_8), dimension(10), intent(inout) :: A
real(kind=real_8), dimension(10,10)   :: B

  !Establish l0 and l1 coefficients
  B(0,0)=1.d0 ; B(1,0)=0.d0 ; B(1,1)=2.d0
  !Return if order is less than two
  if (n>1) then
    do i = 2, n
      B(i,0)=-2.d0*(i-1)*B(i-2,0)
      do j = 1, i
        !Basic recursion relation
        B(i,j)=2.d0*B(i-1,j-1)-2.d0*(i-1)*B(i-2,j)
      end do
    end do
    do i = 0, n
      A(i)=B(n,i)
    end do
  end if
 
end Subroutine Hermite_Coeff


END MODULE dynamique


