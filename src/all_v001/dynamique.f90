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

! TODO change call for Udiag,dyna,CalPre : add delta before theta
!TODO ajouter delta dans le call de dyna depuis le main

!****************************************************************************************
!****************************************************************************************
subroutine dyna(ski,skiCon,fctQ,fctP,sauvChamp,delta,theta,phy,nki,dki,kmin,&
Hqp,Hpq,Hqppq,H0qq,MUqq,H0pp,MUpp,eigenVectH0pp,coeffGS,opt_GS)
!Actual electron dynamics calculations:  Many-electron wavepacket in TDCI form  
!INPUT: 
!       Hamiltonian matrices H0qq, H0pp: core H in Q (P) subspace
!                            Hqp,Hpq: Q-P coupling H (time-independent part), Hqppq=Hqp*Hpq
!       initial CI content (fctQ, fctP),  
!OUTPUT: Final CI vectors (fctQ, fctP, sauvFctQ(?), sauvFctP(?)), TD field (sauvChamp). 
!****************************************************************************************
!****************************************************************************************

Real(kind=real_8),Intent(in)				::theta,phy,delta

Complex(kind=comp_16),Dimension(:,:),Intent(in)		::ski,skiCon,coeffGS
Integer(kind=int_4),Dimension(3),Intent(in)		::nki
Real(kind=real_8),Dimension(3),Intent(inout)		::dki,kmin
Complex(kind=comp_16),Dimension(:),Intent(inout)	::fctQ
Complex(kind=comp_16),Dimension(:,:),Intent(inout)	::fctP
Real(kind=real_8),Dimension(:),Intent(out)		::sauvChamp

Complex(kind=comp_16),Dimension(:,:),Intent(in)		::eigenVectH0pp
Complex(kind=comp_16),Dimension(dimP,nk)		::StateP
Complex(kind=comp_16),Dimension(dimQ,dimQ)		::Mat33,Inv33

Complex(kind=comp_16),Dimension(orb_Q,orb_Q)		::muEPS
Complex(kind=comp_16),Dimension(nk,orb_Q)		::muB_EPS
Integer							::i,j,t,r,s,k,u,v,a,counter,ik1,ik2,ik3,ntk,ik
Integer(kind=int_4),Intent(in)		                ::opt_GS

Complex(kind=comp_16),Dimension(:,:)			::Hqp
Complex(kind=comp_16),Dimension(:,:)			::Hpq
Complex(kind=comp_16),Dimension(:,:)			::Hqppq
Real(kind=real_8),Dimension(:,:)			::H0qq,MUqq
Real(kind=real_8),Dimension(:,:)			::H0pp,MUpp
Real(kind=real_8),Dimension(3)				::eps

complex(kind=comp_16),dimension(nk)			::obs

Real(kind=real_8)					::PreStart,PreEnd,TotalPre
Real(kind=real_8)					::PropStart,PropEnd,TotalProp
Real(kind=real_8)					::TotalCI,TotalGamma,TotalVolk
Real(kind=real_8),Dimension(3)				::sauvkmin
Real(kind=real_8)	        			::t1,t2
write(*,*)'Debut dynamique'



sauvkmin=kmin
sauvChamp=0.d0
! 
sauvChamp(1)=champNint()

!open(48,file='fctQ.dat',status='unknown',form='formatted')
!open(49,file='fctP.dat',status='unknown',form='formatted')
!open(50,file='champ.dat',status='unknown',form='formatted')
!open(51,file='norme.dat',status='unknown',form='formatted')

write(*,*)nt,'number of timesteps'
open(51,file='champ.dat',status='unknown',form='formatted')
!////////////////////////////////////////////////////////
! start loop over time variable (tn), i.e. TIME-PROPAGATION
!////////////////////////////////////////////////////////



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
! Actual time-propagation of fctP, fctQ from tn to tn+1
!
    write(*,*) 'calling propagation, tn=',tn,'/',nt
	call propagation(fctQ,fctP,Hqp,Hpq,H0qq,H0pp,MUqq,MUpp,Mat33,Inv33,ski,skiCon,dki,&
delta,theta,phy,nki,kmin,coeffGS,opt_GS)

	sauvChamp(tn+1)=champNint()
! 
	write(51,*)tmin+(tn*pdt),Int0(),sauvChamp(tn+1)
!
!Time-dependent observables
!

!TODO calculer population_State en transformant populationP
! pour sauver de la memoire
!	call CSF2State(FctP,eigenVectH0pp,StateP)
!	call population_State(tn,StateP,dki,nk,dimP)


	!call population_Q(FctQ,dimQ) !now calculated in propagation
!	call ionisation(FctP,dki,nk,dimP,FctQ,dimQ)
!	call r_WavePacket(FctP,dki,nki,kmin,dimP)
!	write(51,'(5e19.8e5)')normeQ+normeP,normeQ,normeP,(1.d0-normeQ-normeP)
!if (tn .eq. 250)then
!	 call photoelectron(FctP,dki,nki,kmin,dimP)
!	STOP'END OF TEST'
!endif

end do  ! end of time-loop

!TODO Why write this ?
!write(51,*)Int0(),tn
!kmin=sauvkmin
open(48,file='fctQ.dat',status='unknown',form='formatted')
open(49,file='fctP.dat',status='unknown',form='formatted')

do i=1,dimQ
	write(48,*)real(fctQ(i)),aimag(fctQ(i))
end do

 close(48)
 close(49)
 close(50)
 close(51)

 
 call photoelectron(FctP,dki,nki,kmin,dimP)
! call StateElectrons(FctP,eigenVectH0pp,dki,nki,kmin,dimP)
write(*,*)'after photoelectron'

write(*,*)'Fin de la dynamique' 


end subroutine dyna

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


temp3=ZidentiMAT(dimQ) + ((1.d0/4.d0)*Hqppq*(Int0()**2))
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
Mat33=ZidentiMAT(dimQ) - ((1.d0/4.d0)*Hqppq* (Int0()**2))
 

deallocate(work)
deallocate(temp3)
deallocate(ipiv)

end subroutine CalPre


!****************************************************************************************
!****************************************************************************************
subroutine propagation(fctQ,fctP,Hqp,Hpq,H0qq,H0pp,MUqq,MUpp,Mat33,Inv33,ski,skiCon,dki,delta,theta,phy,nki,kmin,coeffGS,opt_GS)
!
!
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(in)	   ::Mat33,Inv33
real(kind=real_8),dimension(:,:),intent(in)	       ::H0qq,H0pp,MUqq,Mupp
complex(kind=comp_16),dimension(:,:),intent(in)	   ::Hqp,Hpq
complex(kind=comp_16),dimension(:,:),intent(in)	   ::ski,skiCon,coeffGS
integer(kind=int_4),dimension(3),intent(in)	       ::nki
Real(kind=real_8),Dimension(3),Intent(inout)	   ::dki
Real(kind=real_8),Dimension(3),Intent(inout)	   ::kmin
real(kind=real_8),intent(in)			           ::theta,phy,delta
Integer(kind=int_4),Intent(in)		               ::opt_GS
complex(kind=comp_16),dimension(:),intent(inout)   ::fctQ
complex(kind=comp_16),dimension(:,:),intent(inout) ::fctP

complex(kind=comp_16),allocatable,dimension(:)     ::CI1,CI2,CI3,CI4,CI5,CI
complex(kind=comp_16),allocatable,dimension(:,:)   ::gamm,Uqq,gamm2,gammCorrection,Upp
complex(kind=comp_16),allocatable,dimension(:)	   ::BigGamma,BigGamma2,BigGamma3
complex(kind=comp_16),allocatable,dimension(:,:)   ::HqpBar
complex(kind=comp_16),allocatable,dimension(:,:)   ::HpqBar
complex(kind=comp_16),allocatable,dimension(:)	   ::ktemp
real(kind=real_8),allocatable,dimension(:)         :: normeQ_CSF
integer						                       ::l,ll,i,j,jj,k,r,s,k1,k2,k3,kk,kkj
real(kind=real_8)				                   ::Et,normQ,normP,P_normfactor

real(kind=real_8),allocatable,dimension(:)	       ::fctPInt,fctPInt2

real(kind=real_8)                                  ::seuil !TODO to be deleted after test
!
!!   Propagator for bound-state dynamics (in Q and/or P subspace) 
!    Uses H0qq(pp)=fireld-free Hamiltonian matrix in bound CSF basis
!         Muqq(pp)=Dipole matrix
!         H0*dt+Mu*Int0() = TD-Hamiltonian, Int0()=field(tn)*Delta t 
!

allocate(Uqq(dimQ,dimQ),Upp(dimP,dimP))
Uqq=dcmplx(0.d0,0.d0)
Upp=dcmplx(0.d0,0.d0)
call Udiag(H0qq,MUqq,dimQ,Uqq,delta)
call Udiag(H0pp,MUpp,dimP,Upp,delta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(BigGamma(dimP*nk),BigGamma2(dimP*nk))	
allocate(HqpBar(dimQ,nk*dimP),HpqBar(nk*dimP,dimQ))
allocate(normeQ_CSF(dimQ))
HqpBar=Hqp*Int0()
HpqBar=Hpq*Int0()



allocate(CI1(dimQ),CI2(dimQ),CI3(dimQ),CI4(dimQ),CI5(dimQ),CI(dimQ))
BigGamma=dcmplx(0.d0,0.d0)
do j=1,dimP
	BigGamma(1+ (j-1)*nk:j*nk)=FctP(j,:)   
end do
!
!!!!!!Propagation within Q-space!!!  
!
 CI1=dcmplx(0.d0,0.d0); CI2=dcmplx(0.d0,0.d0); CI3=dcmplx(0.d0,0.d0)
 CI4=dcmplx(0.d0,0.d0); CI5=dcmplx(0.d0,0.d0); CI=dcmplx(0.d0,0.d0)

call gemv(HqpBar,BigGamma,CI1)   		!Etape1 !; CI1=CI1*dVk(dki) nolonger needed
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
!steps 7,8,9 done (steps 1,2,3; relevant results kept in CI2, CI3)
!
 CI1=dcmplx(0.d0,0.d0); BigGamma2=dcmplx(0.d0,0.d0)
 CI1=(1.d0/2.d0)*CI2 + dcmplx(0.d0,1.d0)*CI3 		!Etape10
call gemv(HpqBar,CI1,BigGamma2) 			!Etape11
BigGamma2=BigGamma - BigGamma2 				!Etape12 
deallocate(CI1,CI2,CI3,BigGamma,HqpBar,HpqBar)
!
! Preparation for propagation among cationic bound-states (action of Upp) 
!
allocate(gamm2(dimP,nk))
do j=1,dimP
	gamm2(j,:)=BigGamma2(1+ (j-1)*nk:j*nk)
end do
deallocate(BigGamma2)
allocate(gamm(dimP,nk))

gamm=dcmplx(0.d0,0.d0)
do k=1,nk
	call gemv(Upp,gamm2(:,k),gamm(:,k))          	!Etape13(a)
end do
! 
select case(opt_GS)
  case(1)
	do j=1,dimP
  		call gemv(coeffGS,gamm(j,:),gamm2(j,:))
	end do
	deallocate(gamm)
	!
	!!!Propagation of ionized electron (Volkov is an option)
	!
	call PW3D_VolkovGS(gamm2,delta,theta,phy,nki,kmin,dki)	!Etape13(b)
	!  
	do j=1,dimP
  		call gemv(conjg(coeffGS),gamm2(j,:),fctP(j,:))
	end do
	!                         		!propagation of P-part of wp done
	deallocate(gamm2)
  case(0)
	!here CI1 is of dimension orb_Q for temporary use
	allocate(CI1(orb_Q),fctPInt(dimP),fctPInt2(dimP))
	CI1=0.d0; fctPInt=0.d0;   fctPInt2=0.d0
	normQ=0.d0; normP=0.d0; P_normfactor=1.d0
	do j=1,dimP
		do i=1,orb_Q
          CI1(i) = DOT_PRODUCT(ski(:,i), gamm(j,:))
		end do
		call gemv(ski,CI1,gamm2(j,:))
	end do
	gamm2=gamm-gamm2
	call PW3D_VolkovGS(gamm2,delta,theta,phy,nki,kmin,dki)	!Etape13(b)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Ajout par J-Nicolas 28/08/17
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	gamm=gamm2
	do j=1,dimP
		do i=1,orb_Q
          CI1(i) = DOT_PRODUCT(ski(:,i), gamm(j,:))
		end do
		call gemv(ski,CI1,gamm2(j,:))
	end do
	gamm2=gamm-gamm2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	fctP=gamm2
	deallocate(gamm2)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! OPW not orthonormalized. Need renormalise continuum wp
        ! Calculation of each CSF population done at the same time
	!!!!!!!!!!!!!!!!!!!!!!!!!!!
  normeQ_CSF=0.d0
	do j=1,dimQ
	    normeQ_CSF(j)=normeQ_CSF(j) + cdabs(fctQ(j))**2.d0
        normQ=normQ+normeQ_CSF(j)
	end do
	open(52,file='CSF_Q.dat',status='unknown',form='formatted',position='append')
	write(52,'(<dimQ>ES24.14)')(normeQ_CSF(i),i=1,dimQ)
	close(52)
  fctPInt=0.d0
  fctPint2=0.d0
  normP=0.d0
	do j=1,dimP
	  do k=1,nk
		fctPInt(j)=fctPInt(j)+cdabs(fctP(j,k))**2
	  end do

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!Mise en commentaire par J-Nicolas 28/08/17
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      do i=1,orb_Q
!        CI1(i) = DOT_PRODUCT(ski(:,i),fctP(j,:))
!        !write(*,*) CDABS(CI1(i))**2
!        fctPint2(j)=fctPint2(j)+CDABS(CI1(i))**2
!      end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      normP=normP+fctPInt(j)-fctPint2(j)
	end do
	open(50,file='ionisation.dat',status='unknown',form='formatted',position='append')
	write(50,*)normP,normQ+normP,1.d0-normQ-normP
    close(50)
    open(52,file='CSF_P.dat',status='unknown',form='formatted',position='append')
    write(52,'(<dimP>ES24.14)')(fctPInt(i)-fctPint2(i),i=1,dimP)
    close(52)
   !TODO faire une boucle pour trouver le seuil minimum avant d avoir une erreur
   !floating point invalid
   seuil=1.d-20
   !do while (seuil.gt.1.d-18)
   ! write(*,*) "seuil=",seuil
    if((normP.gt.seuil).and.(1.d0-normQ.gt.seuil)) then
!    if (normP.gt.0.d0) then
      P_normfactor=dsqrt((1.d0-normQ)/normP)
      fctP=P_normfactor*fctP
    end if
    !seuil=seuil*1.d-1
  !end do
end select




end subroutine propagation

!!****************************************************************************************
!!****************************************************************************************
!subroutine VolkovGS(gamm,coeffGS,MatKa,MatKb)
!!Opérateur Uvolkov avec orthogonalisatiom
!!****************************************************************************************
!!****************************************************************************************
!complex(kind=comp_16),dimension(:,:),intent(inout)::gamm
!complex(kind=comp_16),dimension(:,:),intent(in)::coeffGS
!complex(kind=comp_16),dimension(:,:),intent(in)::MatKa,MatKb
!complex(kind=comp_16),dimension(dimP,nk)::gamm2
!complex(kind=comp_16),dimension(nk)::Volkov
!real(kind=real_8)::kk
!integer::k,j



!do k=1,nk
!	kk=kxmin+ (k*dkx)
!	Volkov(k)=(cdexp(dcmplx(0.d0,-0.5d0*pdt*(kk**2.d0)))*cdexp(dcmplx(0.d0,-0.5d0*Int0()*Int0()*pdt))*cdexp(dcmplx(0.d0,Int0()*kk*pdt)))
!end do

!!gamm2=gamm
!do j=1,dimP
!!	call gemv(MatKb,gamm2(j,:),gamm(j,:))
!	gamm2(j,:)=Volkov(:)*gamm(j,:)
!!	call gemv(MatKa,gamm2(j,:),gamm(j,:))
!end do
!gamm=gamm2
!kxmin=kxmin+Int0()

!end subroutine VolkovGS

!****************************************************************************************
!****************************************************************************************
subroutine PW3D_VolkovGS_kmin(gamm,delta,theta,phy,nki,kmin,dki)
!Opérateur Uvolkov avec orthogonalisatiom
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(inout)::gamm
real(kind=real_8),intent(in)::theta,phy,delta
Integer(kind=int_4),Dimension(3),Intent(in)		::nki
Real(kind=real_8),Dimension(3),Intent(inout)		::kmin
Real(kind=real_8),Dimension(3),Intent(in)		::dki
complex(kind=comp_16),dimension(dimP,nk)::gamm2
complex(kind=comp_16)::Volkov
real(kind=real_8)::kkx,kky,kkz,kk,ksquare
integer::l,k,i,j,counter

Volkov=dcmplx(0.d0,0.d0)

do l=1,dimP
do k=1,nki(3)
	kkz=kmin(3)+ ((k-1)*dki(3))

do i=1,nki(2)
	kky=kmin(2)+ ((i-1)*dki(2))

do j=1,nki(1)
	kkx=kmin(1)+ ((j-1)*dki(1))
	
	counter=j+nki(1)*(i-1)+nki(1)*nki(2)*(k-1)

	kk=kkx*EX(Int0(delta),theta,phy) +kky*EY(Int0(delta),theta,phy) +kkz*EZ(Int0(delta),theta)
	ksquare=kkx**2.d0+kky**2.d0+kkz**2.d0

	Volkov=(cdexp(dcmplx(0.d0,-0.5d0*pdt*ksquare))*cdexp(dcmplx(0.d0,-0.5d0*Int0(delta)*Int0(delta)*pdt))*cdexp(dcmplx(0.d0,-kk*pdt)))
	gamm(l,counter)=Volkov*gamm(l,counter)
end do
end do
end do
end do

!gamm2=gamm
!do j=1,dimP
!!	call gemv(MatKb,gamm2(j,:),gamm(j,:))
!	gamm2(j,:)=Volkov(:)*gamm(j,:)
!!	call gemv(MatKa,gamm2(j,:),gamm(j,:))
!end do


!gamm=gamm2
kmin(1)=kmin(1)+Int0(delta)*EX(1.d0,theta,phy)
kmin(2)=kmin(2)+Int0(delta)*EY(1.d0,theta,phy)
kmin(3)=kmin(3)+Int0(delta)*EZ(1.d0,theta)

end subroutine


!****************************************************************************************
!****************************************************************************************
subroutine PW3D_VolkovGS(gamm2,delta,theta,phy,nki,kmin,dki)
!Opérateur Uvolkov avec orthogonalisatiom
!
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(dimP,nk),intent(inout)	::gamm2
real(kind=real_8),intent(in)				::theta,phy,delta
Integer(kind=int_4),Dimension(3),Intent(in)		::nki
Real(kind=real_8),Dimension(3),Intent(inout)		::kmin
Real(kind=real_8),Dimension(3),Intent(in)		::dki
complex(kind=comp_16),dimension(dimP,nk)		::gamm
complex(kind=comp_16)					::Volkov,lostgamm
real(kind=real_8)					::kkx,kky,kkz,kk,ksquare
integer							::l,k,i,j,counter,counter2
integer,dimension(3)					::dnk

dnk(1)=0.d0
dnk(2)=0.d0
dnk(3)=0.d0

if (dki(1).ne. 0)then
	dnk(1)=int(EX(Int0(delta),theta,phy)/dki(1))
end if
if (dki(2).ne. 0)then
	dnk(2)=int(EY(Int0(delta),theta,phy)/dki(2))
end if
if (dki(3).ne. 0)then
	dnk(3)=int(EZ(Int0(delta),theta)/dki(3))
end if

gamm=dcmplx(0.d0,0.d0)

if (Int0(delta) .gt. 0.0) then

do k=nki(3)-dnk(3),1,-1
kkz=kmin(3)+ ((k-1)*dki(3))

do i=nki(2)-dnk(2),1,-1
kky=kmin(2)+ ((i-1)*dki(2))

do j=nki(1)-dnk(1),1,-1
kkx=kmin(1)+ ((j-1)*dki(1))

	counter = j+nki(1)*(i-1)+nki(1)*nki(2)*(k-1)
	counter2 = j+dnk(1)+nki(1)*(i+dnk(2)-1)+nki(1)*nki(2)*(k+dnk(3)-1)

	kk = kkx*EX(Int0(delta),theta,phy) +kky*EY(Int0(delta),theta,phy) +kkz*EZ(Int0(delta),theta)
	ksquare = kkx**2.d0+kky**2.d0+kkz**2.d0

	Volkov = (cdexp(dcmplx(0.d0,-0.5d0*pdt*ksquare)))
	Volkov = Volkov*cdexp(dcmplx(0.d0,-0.5d0*Int0(delta)*Int0(delta)*pdt))*cdexp(dcmplx(0.d0,-kk*pdt))

	gamm2(:,counter2)=Volkov*gamm2(:,counter)
end do
end do
end do


else if (Int0(delta) .le. 0.d0) then

do k=-dnk(3)+1,nki(3)
kkz=kmin(3)+ ((k-1)*dki(3))

do i=-dnk(2)+1,nki(2)
kky=kmin(2)+ ((i-1)*dki(2))

do j=-dnk(1)+1,nki(1)
kkx=kmin(1)+ ((j-1)*dki(1))

	counter = j+nki(1)*(i-1)+nki(1)*nki(2)*(k-1)
	counter2 = j+dnk(1)+nki(1)*(i+dnk(2)-1)+nki(1)*nki(2)*(k+dnk(3)-1)

	kk = kkx*EX(Int0(delta),theta,phy) +kky*EY(Int0(delta),theta,phy) +kkz*EZ(Int0(delta),theta)
	ksquare = kkx**2.d0+kky**2.d0+kkz**2.d0

	Volkov = (cdexp(dcmplx(0.d0,-0.5d0*pdt*ksquare)))
	Volkov = Volkov*cdexp(dcmplx(0.d0,-0.5d0*Int0(delta)*Int0(delta)*pdt))*cdexp(dcmplx(0.d0,-kk*pdt))

	gamm2(:,counter2)=Volkov*gamm2(:,counter)
end do
end do
end do

end if

end subroutine PW3D_VolkovGS

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

END MODULE dynamique


