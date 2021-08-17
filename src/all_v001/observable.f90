Module observable
!G��re la dynamique dans sa premi��re version
use variables
use Basics
use math
use blas95
use lapack95
use f95_precision
!use MKL_DFTI
use ortho
use qp_integrals

contains

!****************************************************************************************
!****************************************************************************************
subroutine r_WavePacket(FctP,dki,nki,kmin,dimP)
!****************************************************************************************
!****************************************************************************************
Complex(kind=comp_16),dimension(:,:),intent(in)		::fctP
Real(kind=real_8),dimension(3),intent(in)		::dki,kmin
Integer(kind=int_4),dimension(3),intent(in)		::nki
Integer(kind=int_4),intent(in)				::dimP

Complex(kind=comp_16),dimension(:),allocatable		::rWave
Real(kind=real_8),dimension(:),allocatable		::rWavePacket_obs
Real(kind=real_8),dimension(3)				::rmin,dr
Real(kind=real_8)					::rx,ry,rz,kx,ky,kz
Integer(kind=int_4)					::i,j,k1,k2,k3,r1,r2,r3,kk,rr

rmin=0.d0
dr=0.d0

do i=1,3
	if (dki(i) .ne. 0.d0) then
		rmin(i)=-pi/(dki(i))
		dr(i)=2*dabs(rmin(i))/(nki(i)-1)
	endif
end do

open(30,file='rWavePacket.dat',status='unknown',form='formatted',position='append')
allocate(rWave(nk))
rWave=dcmplx(0.d0,0.d0)

do i=1,dimP
	do r3=1,nki(3)
		rz=rmin(3)+(r3-1)*dki(3)
	do r2=1,nki(2)
		ry=rmin(2)+(r2-1)*dki(2)
	do r1=1,nki(1)
		rx=rmin(1)+(r1-1)*dki(1)
		rr=r1+nki(1)*(r2-1)+nki(1)*nki(2)*(r3-1)
		
		do k3=1,nki(3)
			kz=kmin(3)+(k3-1)*dki(3)
		
		do k2=1,nki(2)
			ky=kmin(2)+(k2-1)*dki(2)
		
		do k1=1,nki(1)
			kx=kmin(1)+(k1-1)*dki(1)
			kk=k1+nki(1)*(k2-1)+nki(1)*nki(2)*(k3-1)
			rWave(rr)=rWave(rr)+fctP(i,kk)*(cdexp(dcmplx(0.d0,(sqrt((rx*kx)**2.d0+(ry*ky)**2.d0+(rz*kz)**2.d0)))))
			
		end do
		end do
		end do
	end do
	end do
	end do
end do

allocate(rWavePacket_obs(nk))
do i=1,nk
	rWavePacket_obs(i)=cdabs(rWave(i))**2.d0
end do

write(30,'(<nk>F34.25)')(rWavePacket_obs(i),i=1,nk)

deallocate(rWave,rWavePacket_obs)
close(30)

end subroutine
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
subroutine ionisation(FctP,dki,nk,dimP,FctQ,dimQ)
!****************************************************************************************
!****************************************************************************************
Complex(kind=comp_16),Dimension(:,:),Intent(in)		::fctP
Real(kind=real_8),Dimension(3),Intent(in)		::dki
Integer(kind=int_4),Intent(in)				::nk,dimP
Complex(kind=comp_16),Dimension(:),Intent(in)		::fctQ
Integer(kind=int_4),Intent(in)				::dimQ

Real(kind=real_8),Dimension(dimP)			::fctPInt
integer(kind=int_4)					::i,j,k
Real(kind=real_8),Dimension(dimQ)			::normeQ_CSF
Real(kind=real_8)					::normeQ,normeP

open(50,file='ionisationold.dat',status='unknown',form='formatted',position='append')

fctPInt=0.d0
normeP=0.d0
!TODO ramener ce calcul dans propagation et utiliser le nouveau calcul de normP
do j=1,dimP
	do k=1,nk
		fctPInt(j)=fctPInt(j) + cdabs(FctP(j,k))**2.d0
	end do
	fctPInt(j)=fctPInt(j)! dVk is not needed anylonger *(dVk(dki))
	normeP=fctPInt(j)+normeP
end do

normeQ_CSF=0.d0
normeQ=0.d0

do j=1,dimQ
	normeQ_CSF(j)=normeQ_CSF(j) + cdabs(fctQ(j))**2.d0
	normeQ=normeQ+normeQ_CSF(j)
end do

write(50,*)normeP,normeQ+normeP,1.d0-normeQ-normeP

close(50)
end subroutine

!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
subroutine population_State(tn,StateP,dki,nk,dimP)
!****************************************************************************************
!****************************************************************************************
! TODO state is not calculated correctly for H2+
Complex(kind=comp_16),Dimension(:,:),Intent(in)		::StateP
Real(kind=real_8),Dimension(3),Intent(in)		::dki
Integer(kind=int_4),Intent(in)				::nk,dimP

Real(kind=real_8),Dimension(dimP)			::fctStateInt
integer(kind=int_4)					::i,j,k,tn

open(52,file='Pop_State.dat',status='unknown',form='formatted',position='append')

fctStateInt=0.d0

do j=1,dimP
	do k=1,nk
		fctStateInt(j)=fctStateInt(j) + cdabs(StateP(j,k))**2.d0
	end do
	fctStateInt(j)=fctStateInt(j)!dVk is not needed anylonger *(dVk(dki))
end do

write(52,'(I6,(<dimP>ES24.14))')tn,(fctStateInt(i),i=1,dimP)

close(52)

end subroutine

!****************************************************************************************
!****************************************************************************************
subroutine population_P(FctP,dki,nk,dimP)
!****************************************************************************************
!****************************************************************************************
Complex(kind=comp_16),Dimension(:,:),Intent(in)		::fctP
Real(kind=real_8),Dimension(3),Intent(in)		::dki
Integer(kind=int_4),Intent(in)				::nk,dimP

Real(kind=real_8),Dimension(dimP)			::fctPInt
integer(kind=int_4)					::i,j,k

open(52,file='CSF_P.dat',status='unknown',form='formatted',position='append')

fctPInt=0.d0

do j=1,dimP
	do k=1,nk
		fctPInt(j)=fctPInt(j) + cdabs(FctP(j,k))**2.d0
		write(*,*) "fctPint(",j,")=",fctPInt(j)
	end do
	fctPInt(j)=fctPInt(j)! dVk is not needed anylonger *(dVk(dki))
end do

write(52,'(<dimP>ES24.14)')(fctPInt(i),i=1,dimP)

close(52)

end subroutine
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
subroutine population_Q(FctQ,dimQ)
!****************************************************************************************
!****************************************************************************************
Complex(kind=comp_16),Dimension(:),Intent(in)		:: fctQ
Integer(kind=int_4),Intent(in)			    	    :: dimQ

Real(kind=real_8),Dimension(dimQ)			        :: normeQ_CSF
integer(kind=int_4)					                :: j,i

open(52,file='CSF_Q.dat',status='unknown',form='formatted',position='append')

normeQ_CSF=0.d0

do j=1,dimQ
	normeQ_CSF(j)=normeQ_CSF(j) + cdabs(fctQ(j))**2.d0
end do

write(52,'(<dimQ>ES24.14)')(normeQ_CSF(i),i=1,dimQ)

close(52)

end subroutine
!****************************************************************************************

!****************************************************************************************
!****************************************************************************************
subroutine photoelectron(FctP,dki,nki,kmin,dimP)
!Ecriture fichiers -> MatLab/Mathematica
!****************************************************************************************
!****************************************************************************************
Complex(kind=comp_16),Dimension(:,:),Intent(inout)	::fctP
Real(kind=real_8),Dimension(3),Intent(inout)		::dki,kmin
Integer(kind=int_4),Dimension(3),Intent(in)		::nki
Integer(kind=int_4),Intent(in)				::dimP

Complex(kind=comp_16),allocatable,dimension(:,:)	::fctPSomme
Real(kind=real_8),allocatable,dimension(:,:)		::fctPSomme2
Complex(kind=comp_16),allocatable,dimension(:)        ::fctPSumcoerr

Real(kind=real_8),allocatable,dimension(:)		::Mx,My,Mz
Real(kind=real_8),allocatable,dimension(:,:)		::E
Integer							::i,j,ik,kx,ky,kz,ktot
Real(kind=real_8)					::kk,kkx,kky,kkz
Character(len=5)					::nom
Character(len=23)					::canal
CHARACTER(LEN=20) :: FMT,FMTX

!write(*,*)'start obs'
open(52,file='PSomme.dat',status='replace',form='formatted')
!do j=1,dimP
!    write (nom,'(I5.5)') j
!    open(90+j,file="spectre2d_canal"// ADJUSTL(nom) //".dat",status='replace',form='formatted')
!end do
open(70,file="spectre_tot_coer.dat",status='unknown',form='formatted')
open(80,file="ky.dat",status='unknown',form='formatted')
open(81,file="kzy.dat",status='unknown',form='formatted')

!open(61,file='spectre_canal_1.dat',status='unknown',form='formatted')
!open(62,file='spectre_canal_2.dat',status='unknown',form='formatted')
!open(63,file='spectre_canal_3.dat',status='unknown',form='formatted')
ktot=nki(1)*nki(2)*nki(3)

allocate(fctPSomme(dimP,ktot),fctPSomme2(dimP,ktot),E(ktot,4),Mx(nki(1)),My(nki(2)),Mz(nki(3)))
allocate(fctPSumcoerr(ktot))
fctPSomme=0.d0
fctPSomme2=0.d0
fctPSumcoerr=0.d0
E=0.d0


!write(*,*)kmin
!write(*,*)dki
!write(*,*)nki
if (nki(1).gt.1) then
  !open(91,file="kx.dat",status='replace',form='formatted')
  open(102,file="spectreXZ_tot_coer.dat",status='replace',form='formatted')
  open(61,file="kzx.dat",status='unknown',form='formatted')
  open(62,file="kx.dat",status='unknown',form='formatted')
  do j=1,dimP
    write (nom,'(I5.5)') j
    open(110+j,file="spectre2dXZ_canal"// ADJUSTL(nom) //".dat",status='replace',form='formatted')
  end do
end if
WRITE(FMT,*) nki(2)
WRITE(FMTX,*) nki(1)

do kz=1,nki(3)
  do ky=1,nki(2)
    do kx=1,nki(1)
      kkx=kmin(1)+real(kx)*dki(1)
      kky=kmin(2)+real(ky)*dki(2)
      kkz=kmin(3)+real(kz)*dki(3)
      kk=(kkx*kkx + kky*kky + kkz*kkz)**(1.d0/2.d0)
      ik=kx + nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1)
      fctPSomme(:,ik)=fctPSomme(:,ik) + FctP(:,ik)
      fctPSomme2(:,ik)=fctPSomme2(:,ik) + real(FctP(:,ik))**2.d0 + aimag(FctP(:,ik))**2.d0
      do i=1,dimP
        fctPSumcoerr(ik)=fctPSumcoerr(ik)+FctP(i,ik)
      end do
      E(ik,1)=kk*kk/(2.d0)
      E(ik,2)=kkx
      E(ik,3)=kky
      E(ik,4)=kkz
    end do
    if (nki(1).gt.1) then
      write(102,'('// ADJUSTL(FMTX) //'ES24.14)') (abs(fctPSumcoerr(i + nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1)))**(2.d0),i=1,nki(1))
      write(62,'('// ADJUSTL(FMTX) //'ES24.14)') (E(i+ nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1),2),i=1,nki(1))
      write(61,'('// ADJUSTL(FMTX) //'ES24.14)') (E(i+ nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1),4),i=1,nki(1))

      !write(62,'('// ADJUSTL(FMTX) //'ES24.14)') (E(1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1),4),i=1,nki(1))
      do j=1,dimP
        write(110+j,'('// ADJUSTL(FMTX) //'ES24.14)') (fctPSomme2(j,i+ nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1)),i=1,nki(1))
      end do
    end if

  end do
  write(70,'('// ADJUSTL(FMT) //'ES24.14)') (abs(fctPSumcoerr(1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1)))**(2.d0),i=1,nki(2))
  write(80,'('// ADJUSTL(FMT) //'ES24.14)') (E(1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1),3),i=1,nki(2))
  write(81,'('// ADJUSTL(FMT) //'ES24.14)') (E(1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1),4),i=1,nki(2))

  !do j=1,dimP
  !  write(90+j,'('// ADJUSTL(FMT) //'ES24.14)') (fctPSomme2(j,1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1)),i=1,nki(2))
  !end do
end do

!do j=1,dimP
!    close(90+j)
!end do
close(70)
close(80)
close(81)
do j=1,dimP


write (nom,'(I5.5)') j

! ajouter spectrek_canal ... avec E(ik,1)
canal='spectre_canal_'//trim(nom)//'.dat'

open(60,file=canal,status='unknown',form='formatted')
	if (nki(1) .eq. 1) then
	do i=1,ktot
		write(52,*) real(fctPSomme(j,i)),aimag(fctPSomme(j,i))
		write(60,*) E(i,3),E(i,4),fctPSomme2(j,i)
	end do

	else if (nki(2) .eq. 1) then
	do i=1,ktot
		write(52,*) real(fctPSomme(j,i)),aimag(fctPSomme(j,i))
		write(60,*) E(i,2),E(i,4),fctPSomme2(j,i)
	end do

	else if (nki(3) .eq. 1) then
	do i=1,ktot
		write(52,*) real(fctPSomme(j,i)),aimag(fctPSomme(j,i))
		write(60,*) E(i,2),E(i,3),fctPSomme2(j,i)
	end do
	endif

close(60)

end do

 close(52)
! close(61)
! close(62)
! close(63)
!TODO add other close
if (nki(1).gt.1) then
  do j=1,dimP
    close(110+j)
  end do
end if

end subroutine

!****************************************************************************************
!****************************************************************************************
subroutine StateElectrons(FctP,eigenVectH0pp,dki,nki,kmin,dimP)
!Ecriture fichiers -> MatLab/Mathematica
!****************************************************************************************
!****************************************************************************************
Complex(kind=comp_16),Dimension(:,:),Intent(in)		::fctP,eigenVectH0pp
Real(kind=real_8),Dimension(3),Intent(inout)		::dki,kmin
Integer(kind=int_4),Dimension(3),Intent(in)		::nki
Integer(kind=int_4),Intent(in)				::dimP

Complex(kind=comp_16),allocatable,dimension(:,:)	::StateElectron
Real(kind=real_8),allocatable,dimension(:,:)		::fctPSomme2
Real(kind=real_8),allocatable,dimension(:)		::Mx,My,Mz
Real(kind=real_8),allocatable,dimension(:,:)		::E
Integer							::i,j,ik,kx,ky,kz,ktot
Real(kind=real_8)					::kk,kkx,kky,kkz
 Character(len=5)					::nom 
 Character(len=23)					::canal

ktot=nki(1)*nki(2)*nki(3)

allocate(StateElectron(dimP,ktot),fctPSomme2(dimP,ktot),E(ktot,4),Mx(nki(1)),My(nki(2)),Mz(nki(3)))
StateElectron=0.d0
fctPSomme2=0.d0
E=0.d0

 call gemm(eigenVectH0pp,fctP,StateElectron)

do kz=1,nki(3)
do ky=1,nki(2)
do kx=1,nki(1)
	kkx=kmin(1)+real(kx)*dki(1)
	kky=kmin(2)+real(ky)*dki(2)
	kkz=kmin(3)+real(kz)*dki(3)
	kk=(kkx*kkx + kky*kky + kkz*kkz)**(1.d0/2.d0)
	ik=kx + nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1)
	fctPSomme2(:,ik)=fctPSomme2(:,ik) + real(StateElectron(:,ik))**2.d0 + aimag(StateElectron(:,ik))**2.d0
	E(ik,1)=kk*kk/(2.d0)
	E(ik,2)=kkx
	E(ik,3)=kky
	E(ik,4)=kkz
end do
end do
end do

do j=1,dimP

write (nom,'(I5.5)') j

 canal='spectre_state_'//trim(nom)//'.dat'


open(60,file=canal,status='unknown',form='formatted')
do i=1,ktot
	write(60,*) E(i,2),E(i,4),fctPSomme2(j,i)
end do
 close(60)

end do

! close(52)
! close(61)
! close(62)
! close(63)


end subroutine

!****************************************************************************************
!****************************************************************************************
subroutine norm(sauvChamp,sauvFctQ,sauvFctP,ski,skiCon,dki,nki,kmin)
!Ecriture fichiers -> MatLab/Mathematica
!****************************************************************************************
!****************************************************************************************
complex(kind=comp_16),dimension(:,:),intent(in)::sauvFctQ
complex(kind=comp_16),dimension(:,:,:),intent(inout)::sauvFctP
real(kind=real_8),dimension(:),intent(in)::sauvChamp
complex(kind=comp_16),dimension(:,:,:,:),intent(in):: ski,skiCon
Real(kind=real_8),Dimension(3),Intent(in)		::dki,kmin
Integer(kind=int_4),Dimension(3),Intent(in)		::nki

real(kind=real_8),allocatable,dimension(:)::normeQ,normeP
real(kind=real_8),allocatable,dimension(:,:)::normeQ_CSF,fctPInt,fctPSomme
 character(len=32)::filename
integer::t,j,i,k,kx,ky,kz,kkx,kky,kkz,kk,ik

allocate(normeQ(nt+1),normeP(nt+1),fctPInt(dimP,nt+1),normeQ_CSF(dimQ,nt+1),fctPSomme(nk,nt+1))


open(48,file='fctQ.dat',status='unknown',form='formatted')
open(49,file='fctP.dat',status='unknown',form='formatted')
open(50,file='champ2.dat',status='unknown',form='formatted')
open(51,file='norme.dat',status='unknown',form='formatted')
filename='obs_'//job
open(60,file=filename,status='unknown',form='formatted')
open(61,file='P1.dat',status='unknown',form='formatted')
open(62,file='P2.dat',status='unknown',form='formatted')
open(63,file='P3.dat',status='unknown',form='formatted')
open(64,file='P4.dat',status='unknown',form='formatted')
open(65,file='PSomme.dat',status='unknown',form='formatted')

nk=nki(1)*nki(2)*nki(3)

Write(*,*)nk
normeQ=0.d0
fctPInt=0.d0
normeP=0.d0
normeQ_CSF=0.d0
fctPSomme=0.d0

do t=1,nt,1
!	call OPW2PW_3D(sauvfctP(:,:,t),ski,skiCon,nki,dki)
	do j=1,dimQ
		normeQ_CSF(j,t)= cdabs(sauvFctQ(j,t))**2.d0
		normeQ(t)=normeQ(t) + normeQ_CSF(j,t)
	end do
	write(48,*)
!	write(48,*)'t',(i,i=1,dimQ)
	write(48,'(9e19.8e5)')(tmin+(t*pdt)),(normeQ_CSF(i,t),i=1,dimQ),normeQ(t),sauvChamp(t)
	write(48,*)

	do j=1,dimP
		do k=1,nk
			fctPInt(j,t)=fctPInt(j,t) + cdabs(sauvFctP(j,k,t))**2.d0
		end do
		fctPInt(j,t)=fctPInt(j,t)*(dki(1)*dki(2)*dki(3))
		normeP(t)=normeP(t) +fctPInt(j,t)
	end do
	write(49,*)
	write(49,'(11e19.8e5)')(tmin+(t*pdt)),(fctPInt(i,t),i=1,dimP),normeP(t),sauvChamp(t)
	write(49,*)
	write(50,*)(tmin+(t*pdt)),sauvChamp(t)
	write(51,'(5e19.8e5)')normeQ(t)+normeP(t),normeQ(t),normeP(t),(1.d0-normeQ(t)-normeP(t))
	write(60,'(6e19.8e5)')(tmin+(t*pdt)),normeQ_CSF(1,t),normeQ_CSF(2,t),fctPInt(1,t),fctPInt(2,t),sauvChamp(t)
	!call OPW2PW(sauvFctP(:,:,t),ski,skiCon)
!	call OPW2PW(sauvFctP(:,:,t),ski,skiCon,nki(1),dki(1))
!	call OPW2PW(sauvFctP(:,nk(1)+1:nk(2),t),skiy,skiCony,nk(2),dki(2))
!	call OPW2PW(sauvFctP(:,nk(2)+1:nk(3),t),skiz,skiConz,nk(3),dki(3))
	do j=1,dimP
	do kx=1,nki(1)
	do ky=1,nki(2)
	do kz=1,nki(3)
		kkx=kmin(1)+kx*dki(1)
		kky=kmin(2)+ky*dki(2)
		kkz=kmin(3)+kz*dki(3)
		kk=(kkx*kkx + kky*kky + kkz*kkz)**(1.d0/2.d0)
		ik=kx+nki(1)*(ky-1)+nki(1)*nki(2)*(kz-1)
		fctPSomme(ik,t)= fctPSomme(ik,t) + sauvFctP(j,ik,t)**2.d0 * kk
	end do
	end do
	end do
	end do

	write(61,'(1024e19.8e5)')(sauvFctP(1,k,t),k=1,nk)
	write(62,'(1024e19.8e5)')(sauvFctP(2,k,t),k=1,nk)
	write(63,'(1024e19.8e5)')(sauvFctP(3,k,t),k=1,nk)
	write(64,'(1024e19.8e5)')(sauvFctP(4,k,t),k=1,nk)
!	write(65,'(1024e19.8e5)')(FctPSomme(k,t),k=1,nk)
end do
write(*,*)kmin, dki
write(65,'(1024e19.8e5)')(FctPSomme(k,nt),k=1,nk)

 close(48)
 close(49)
 close(50)
 close(51)
 close(60)
 close(61)
 close(62)
 close(63)
 close(64)
 close(65)

end subroutine


!!****************************************************************************************
!!****************************************************************************************
!subroutine test()
!!Subroutine de test
!!****************************************************************************************
!!****************************************************************************************
!complex(kind=comp_16),dimension(1024)::vec1
!integer(kind=int_4)::n,Status
!type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
!real(kind=real_8)::somme1,somme2,somme3

!open(60,file='test.dat',status='unknown',form='formatted')
!somme1=0.d0
!do n=1,1024
!	vec1(n)=exp(-0.005*(n-512)**2)
!	somme1=somme1+abs(vec1(n))
!end do

!vec1=vec1/somme1
!somme1=0.d0
!do n=1,1024
!	write(60,*)n,real(vec1(n)),aimag(vec1(n))
!	somme1=somme1+abs(vec1(n))
!end do
!write(*,*)'somme1',somme1
! close(60)

!Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, 1024 )
!Status = DftiCommitDescriptor( My_Desc1_Handle )
!Status = DftiComputeForward( My_Desc1_Handle, vec1 )
!Status = DftiFreeDescriptor(My_Desc1_Handle)

!open(60,file='test2.dat',status='unknown',form='formatted')
!somme2=0.d0
!do n=1,1024
!	somme2=somme2+abs(vec1(n))
!end do

!do n=1,1024
!	write(60,*)n,real(vec1(n)),aimag(vec1(n))
!end do
!write(*,*)'somme2',somme2
! close(60)


!open(60,file='test3.dat',status='unknown',form='formatted')
!Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, 1024 )
!Status = DftiCommitDescriptor( My_Desc1_Handle )
!Status = DftiComputeForward( My_Desc1_Handle, vec1 )
!Status = DftiFreeDescriptor(My_Desc1_Handle)
!somme3=0.d0
!do n=1,1024
!	write(60,*)n,real(vec1(n)),aimag(vec1(n))
!	somme3=somme3+abs(vec1(n))
!end do
!write(*,*)'somme3',somme3
! close(60)
!end subroutine test

end module observable


