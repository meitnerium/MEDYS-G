Module nouveaudrt
!Génère la table DRT
use basics
use math

contains

!****************************************************************************************
!****************************************************************************************
subroutine dimension_DRT(S,Ne,Norb,lDRT,dimCSF)
!****************************************************************************************
!****************************************************************************************
Integer(kind=int_4),intent(in)::Ne,Norb
Integer(kind=int_4)::Norb_fact
Real(kind=real_8),intent(in)::S
Integer(kind=int_8),intent(out)::lDRT,dimCSF

Integer(kind=int_4)::atemp,btemp,ctemp,dtemp
Real(kind=real_8)::dimCSF_r

lDRT=0.d0
dimCSF=0.d0

btemp=2*S
atemp=(Ne - btemp)/2.0
 ctemp=Norb-(atemp+btemp)
dtemp=min(atemp,ctemp)
Norb_fact=Norb+1
lDRT=((atemp+1)*(ctemp+1)*((btemp)+(dtemp/2.0)+(1)))-((1.0/6.0)*(dtemp)*(dtemp+1)*(dtemp+2))
dimCSF_r=((btemp+1.d0)/(Norb+1.d0))*((fact(Norb_fact,1))/(fact(Norb_fact-atemp,1)*fact(atemp,1)))*(fact(Norb_fact,1)/(fact(Norb_fact-ctemp,1)*fact(ctemp,1)))
!dimCSF_r=((btemp+1.d0)/(Norb+1.d0))*((fact(Norb+1,1))/(fact(Norb+1-atemp,1)*fact(atemp,1)))*(fact(Norb+1,1)/(fact(Norb+1-ctemp,1)*fact(ctemp,1)))
write(*,*) "dimCSF_r in dimension_DRT = ",dimCSF_r
dimCSF=dimCSF_r
end subroutine

!****************************************************************************************
!****************************************************************************************
subroutine gen_drt(S,Ne,Norb_cat,Norb_sym,restriction,tabDRT,chemin,lDRT,dimCSF,FC,dimQ)
!subroutine gen_drt(S,Ne,Norb,res)
! génère une table DRT à partir de la valeur du spin, du nombre d'électrons et d'orbitales

!A DVLP augmenter le nombre de ionisation possible

!A DVLP la numérotation des "chemins"

!A DVLP supprimer les chemins entrainant l'ionisation des orbitales de coeur
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),intent(in)::Ne,Norb_cat,Norb_sym,FC
real(kind=real_8),intent(in)::S
logical,intent(in)::restriction
integer(kind=int_4),dimension(:,:),intent(out)::tabDRT
integer(kind=int_4),dimension(:,:),intent(out)::chemin
integer(kind=int_8),intent(inout)::lDRT,dimCSF,dimQ

integer(kind=int_4), allocatable, dimension(:,:)::abcAutoris,abc,vecD,vecJ,vecY
integer(kind=int_4), allocatable, dimension(:)::comptJ,abctemp,xx
integer(kind=int_4)::dimDRT,Norb
integer(kind=int_4)::atemp,btemp,ctemp,dtemp,pkp,km,km_temp,jm,id,indjm,jj,currentj,compteurj,cas,maxID,i,j,k,LimA
integer(kind=int_4),allocatable,dimension(:)::km_var,restantligne
integer(kind=int_4)::lDRT_reel

tabDRT=0.d0
 chemin=0.d0
Norb=Norb_cat+Norb_sym
Norb=Norb+1 !!!Pour les besoins du programme, non transmis au main
!write(*,*)'DEBUT DRT'
allocate(abcAutoris(3,4),comptJ(Norb),vecJ(Norb,10))
abcAutoris(:,1)=(/0,0,1/)
abcAutoris(:,2)=(/0,1,0/)
abcAutoris(:,3)=(/1,-1,1/)
abcAutoris(:,4)=(/1,0,0/)
compteurj=1
comptJ=0
comptJ(Norb)=1
vecJ=0
vecJ(Norb,:)=(/1/)
btemp=2*S
atemp=(Ne - btemp)/2.0
 ctemp=Norb-(atemp+btemp)-1
dtemp=min(atemp,ctemp)
!lDRT=((atemp+1)*(ctemp+1)*((btemp)+(dtemp/2.0)+(1)))-((1.0/6.0)*(dtemp)*(dtemp+1)*(dtemp+2))
!write(*,*)S,Ne,Norb-1,'lDRT_sub=',lDRT
allocate(vecD(lDRT,4),abctemp(3),abc(lDRT,3))
vecD=0
abc(1,:)=(/atemp,btemp,ctemp/)
!write(*,*)'abc(1)',abc(1,:)

!!!Table DRT :: vecteur ABC, vecteur J
do km=Norb-1,1,-1
	comptJ(km)=0
	do indjm=1,comptJ(km+1)
		jm=vecJ(km+1,indjm)
		!write(*,*)'jm=',jm
		if (km.eq.Norb-1 .and. restriction.eq..TRUE.) then
			maxID=3
		else
			maxID=4
		end if
		do id=1,maxID
			abctemp=abc(jm,:)-abcAutoris(:,id)
			cas=0
			if(km.gt.FC) then
				LimA=FC
			else
				LimA=0
			end if
			if(abctemp(1).ge.LimA .and. abctemp(2).ge.0 .and. abctemp(3).ge.0) then
				if (comptJ(km).gt.0) then
				do jj=1,comptJ(km)
					if(abctemp(1).ne.abc((vecJ(km+1,comptJ(km+1))+jj),1) .or. abctemp(2).ne.abc((vecJ(km+1,comptJ(km+1))+jj),2) .or. abctemp(3).ne.abc((vecJ(km+1,comptJ(km+1))+jj),3)) then
						cas=1
						currentj=compteurj+1
					else
						cas=2
						currentj=(vecJ(km+1,comptJ(km+1))+jj)
						exit
					end if
				end do
				else
					cas=1
					currentj=compteurj+1
				end if
			else
				cas=3
			end if
			select case (cas)
				case(1)
				  abc(currentj,1)=abctemp(1)
				  abc(currentj,2)=abctemp(2)
				  abc(currentj,3)=abctemp(3)
				  vecD(jm,id)=currentj
				  comptJ(km)=comptJ(km)+1
				  vecJ(km,comptJ(km))=currentj
				  compteurj=compteurj+1
				  !write(*,*)
				case(2)
				  vecD(jm,id)=currentj
				  !write(*,*)
				case(3)
				  !write(*,*)
			end select

		end do		
	end do
end do

pkp=0
do km=Norb,1,-1 !!!
	do jj=1,comptJ(km)
		pkp=pkp+1
	end do
end do

!!!!!!Table DRT: Vecteur Y et valeur X
allocate(xx(pkp),vecY(pkp,4))
xx=0
xx(pkp)=1
vecY=0
do i=pkp-1,1,-1
	xx(i)=0
	if (vecD(i,1).ne.0) xx(i)=xx(i)+xx(vecD(i,1))
	if (vecD(i,2).ne.0) xx(i)=xx(i)+xx(vecD(i,2))
	if (vecD(i,3).ne.0) xx(i)=xx(i)+xx(vecD(i,3))
	if (vecD(i,4).ne.0) xx(i)=xx(i)+xx(vecD(i,4))
	if (vecD(i,2).ne.0) then
		if (vecD(i,1).ne.0) then
			vecY(i,2)=vecY(i,1)
			if (vecD(i,1).ne.0) vecY(i,2)=vecY(i,2)+xx(vecD(i,1))
		end if
	end if
	if (vecD(i,3).ne.0) then
		if (vecD(i,2).eq.0) then
			if (vecD(i,1).ne.0) then
				vecY(i,3)=vecY(i,1)
				if (vecD(i,1).ne.0) vecY(i,3)=vecY(i,3)+xx(vecD(i,1))
			end if
		else
			vecY(i,3)=vecY(i,2)+xx(vecD(i,2))
		end if
	end if
	if (vecD(i,4).ne.0) then
		if (vecD(i,3).eq.0) then
			if (vecD(i,2).eq.0) then
				if (vecD(i,1).ne.0) then
					vecY(i,4)=vecY(i,1)
					if (vecD(i,1).ne.0) vecY(i,4)=vecY(i,4)+xx(vecD(i,1))
				end if
			else
				vecY(i,4)=vecY(i,2)
				if (vecD(i,2).ne.0) vecY(i,4)=vecY(i,4)+xx(vecD(i,2))
			end if
		else
			vecY(i,4)=vecY(i,3)
			if (vecD(i,3).ne.0) vecY(i,4)=vecY(i,4)+xx(vecD(i,3))
		end if
	end if
end do

pkp=0
lDRT_reel=0
do km=Norb,1,-1 !!!
	do jj=1,comptJ(km)
		i=vecJ(km,jj)
!		write(*,'(3i4,a2,4i4,a2,i4,a2,4i4)') (abc(i,k),k=1,3),'|',(vecD(i,j),j=1,4),'|',xx(i),'  ',(vecY(i,jm),jm=1,4)
		pkp=pkp+1
		lDRT_reel=lDRT_reel+1
	end do
	!write(*,*)
end do

dimDRT=pkp
!allocate(tabDRT(pkp,14))
allocate(km_var(lDRT_reel))
tabDRT(:,:)=0

do km=Norb,1,-1 !!!
	do jj=1,comptJ(km)
		i=vecJ(km,jj)
		km_temp=km
		tabDRT(i,:)=(/km_temp,i,abc(i,1),abc(i,2),abc(i,3),vecD(i,1),vecD(i,2),vecD(i,3),vecD(i,4),vecY(i,1),vecY(i,2),vecY(i,3),vecY(i,4),xx(i)/)
		km_var(i)=km
	end do
	!write(*,*)
end do
lDRT=i
dimCSF=xx(1)
dimQ=xx(2)
write(*,*) "dimDRT=",dimDRT
!do i=1,dimDRT
!write(*,'(2i4,3x,3i4,3x,4i4,3x,4i4,3x,i4)')(tabDRT(i,j),j=1,14)
!end do
!write(*,*)'FIN DRT'

!!!!Inventaire des CSF -> chemin (inventaire des CSF exprimés en noeuds)
allocate(restantligne(Norb))
 chemin=0
 restantligne=1
 call cascade(Norb,1,restantligne,km_var,xx,vecD,chemin)
!!!Ecriture
!write(*,*)
!do i=1,xx(1)
!	write(*,'(6i4)') (chemin(j,i),j=1,Norb)
!end do
!!!!!!!!!!!
!!!!!1test gen_gen
!allocate(Ers(xx(1),xx(1),norb-1,norb-1))  ! <1|E34|2>
! call gen_gen(tabDRT,chemin,norb-1,xx(1),Ers)
!do i=1,xx(1)
!do j=1,xx(1)
!if (Ers(i,j,1,3).ne.0.d0) then
!write(*,*)i,j,Ers(i,j,1,3)
!end if
!end do
!end do
!write(*,*)'km',tabDRT(:,1)
!write(*,*)'dimQ',dimQ
deallocate(abcAutoris,abc,vecD,vecJ,vecY,comptJ,abctemp,xx)

end subroutine gen_drt


!****************************************************************************************
!****************************************************************************************
recursive subroutine cascade(Norb,ligne,restantligne,km_var,xx,vecD,chemin)
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),dimension(:,:)::chemin
integer(kind=int_4), dimension(:,:)::vecD
integer(kind=int_4),dimension(:)::km_var,restantligne
integer(kind=int_4)::ligne,Norb
integer(kind=int_4),dimension(:)::xx
integer::i,j,debut,fin

!write(*,*)'ligne',ligne 
if (restantligne(Norb-km_var(ligne)+1).le.xx(1)) then  ! restantligne -> nombre de lignes restantes - nombre de ligne déja complétées
	debut=restantligne(Norb-km_var(ligne)+1)
	fin=debut + xx(ligne) -1
	chemin(Norb+1-km_var(ligne),debut:fin)=ligne
	restantligne(Norb-km_var(ligne)+1)=fin+1

	do i=1,4
		if (vecD(ligne,i).ne.0) then
			call cascade(Norb,vecD(ligne,i),restantligne,km_var,xx,vecD,chemin) 
		end if
	end do

end if
end subroutine

!****************************************************************************************
!****************************************************************************************
subroutine gen_gen(tabDRT,chemin,Norb,dimCSF,Ers)
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),dimension(:,:),intent(in)::tabDRT
integer(kind=int_4),dimension(:,:),intent(in)::chemin
integer(kind=int_4),intent(in)::Norb
integer(kind=int_8),intent(in)::dimCSF
real(kind=real_8),dimension(:,:,:,:),intent(out)::Ers !dimSCF,dimSCF,Norb,Norb
real(kind=real_8)::res
integer(kind=int_4),dimension(4)::b
integer(kind=int_4)::pointeur
integer::i,j,II,JJ,l

Ers=0.d0
!!!!!!Cas i,i
do i=1,Norb
	do II=1,dimCSF
		if (i.ne.1) then
			Ers(II,II,i,i)=(2.d0*tabDRT(chemin(Norb+1-i,II),3)+tabDRT(chemin(Norb+1-i,II),4) ) -(2.d0*tabDRT(chemin(Norb-i+2,II),3)+tabDRT(chemin(Norb-i+2,II),4))
		else
			Ers(II,II,i,i)=2.d0*tabDRT(chemin(Norb+1-i,II),3)+tabDRT(chemin(Norb+1-i,II),4)
		end if
	end do
end do

!!!!!!Cas i,i+1
do II=1,dimCSF-1
do JJ=II+1,dimCSF
	pointeur=lycos(Norb,chemin(:,II),chemin(:,JJ))
!	write(*,*)pointeur
	if (pointeur.ne.0) then
		b(1)=tabDRT(chemin(pointeur-1,II),4)
		b(2)=tabDRT(chemin(pointeur,II),4)
		b(3)=tabDRT(chemin(pointeur,JJ),4)
		b(4)=tabDRT(chemin(pointeur+1,JJ),4)
		call biblio(b,res)
		Ers(II,JJ,tabDRT(chemin(pointeur,II),1)-1,tabDRT(chemin(pointeur,II),1))=res
	end if
end do
end do

!!!!!!!Cas i,i+n avec n>1
do j=2,Norb
do i=1,norb-j
	l=i+j-1
	Ers(:,:,i,i+j)=matmul(Ers(:,:,i,l),Ers(:,:,l,i+j)) - matmul(Ers(:,:,l,i+j),Ers(:,:,i,l))
end do
end do

!!!!!!!Cas i,i+n avec n>1
!do j=2,Norb
!do i=1,norb-j
!	l=i+j-1
!	Ers(:,:,i,i+j)=matmul(Ers(:,:,i,l),Ers(:,:,l,i+j)) - matmul(Ers(:,:,l,i+j),Ers(:,:,i,l))
!end do
!end do
!!!!!!Reconstitution de l'ensemble i,j<->j,i
do i=1,Norb-1
do j=i+1,Norb
	Ers(:,:,j,i)=transpose(Ers(:,:,i,j))
end do
end do


end subroutine gen_gen

!****************************************************************************************
!****************************************************************************************
integer function lycos(Norb,vec1,vec2)
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),dimension(:),intent(in)::vec1,vec2
integer(kind=int_4),intent(in)::Norb
integer(kind=int_4)::compteur,pointeur
integer::i

compteur=0
pointeur=0
do i=1,Norb
	if (vec1(i).ne.vec2(i)) then
		compteur=compteur+1
		pointeur=i
	end if
end do

if (compteur.eq.1) then
	lycos=pointeur
else
	lycos=0
end if

end function 

!****************************************************************************************
!****************************************************************************************
subroutine biblio(b,res)
!****************************************************************************************
!****************************************************************************************
integer(kind=int_4),dimension(:),intent(in)::b
real(kind=real_8),intent(out)::res
integer(kind=int_4)::dBp,dB,dBd

dBp=b(2)-b(1)
 dB=b(3)-b(1)
dBd=b(4)-b(1)
res=0.d0

!Cas possibles dans les conditions supprimes?
if (b(1).ne.b(2).or.b(2).ne.b(3).or.b(3).ne.b(4)) then
select case(dBp)
case(0)
	select case(dBd)
	case(0)
		select case(dB)
			case(1)
				res=sqrt((b(1)+2.d0)/(b(1)+1.d0))
			case(-1)
				res=sqrt((b(1))/(b(1)+1.d0))
			case default
				write(*,*)'ERREUR dB de dBd de dBp'
				write(*,'(4i3)')b
		end select
	case(1)
		res=1.d0
	case(-1)
		res=1.d0
	case default
		write(*,*)'ERREUR dBd0 de dBp'
		write(*,'(4i3)')b
	end select
case(1)
	if (dB.eq.0) then
		select case(dBd)
		case(0)
			res=sqrt((b(1)+2.d0)/(b(1)+1.d0))
		case(1)
			res=1.d0
		case default
			write(*,*)'ERREUR dBd1 de dBp'
			write(*,'(4i3)')b
		end select
	end if
case(-1)
	if (dB.eq.0) then
		select case(dBd)
		case(0)
			res=sqrt((b(1))/(b(1)+1.d0))
		case(-1)
			res=1.d0
		case default
			write(*,*)'ERREUR dBd-1 de dBp'
			write(*,'(4i3)')b
		end select
	end if
case default
	write(*,*)'ERREUR dBp'
	write(*,'(4i3)')b
end select
end if

end subroutine

!*************************************************************************************************
subroutine gen_bielect(dimCSF,Norb,Ersmp,Ers,tabDRT,chemin,Ne,S)
!Routine tres peu efficace pour calculer les generateurs bielectroniques
!N'est pas trop lent pour le moment... Probleme potentiel plus tard
!*************************************************************************************************
integer(kind=int_4),intent(in)::Norb,Ne
integer(kind=int_8),intent(in)::dimCSF
 Real(kind=real_8), Intent(out)::S
integer(kind=int_4),dimension(:,:),intent(in)::tabDRT
integer(kind=int_4),dimension(:,:),intent(in)::chemin
real(kind=real_8),dimension(:,:,:,:),intent(in)::Ers
logical,dimension(dimCSF,dimCSF)::CheckCondon
real(kind=real_8),dimension(:,:,:,:,:,:),intent(out)::Ersmp
integer::r,q,m,p,II,JJ

!Verifier les regles de slater-condon
 CheckCondon=.TRUE.
 call ReglesSlaterCondon(dimCSF,Norb,tabDRT,chemin,Ne,S,CheckCondon)

!write(*,*)' Calcul des generateur bi-electroniques...'
do II=1,dimCSF
do JJ=1,II
	if (CheckCondon(II,JJ).eq..true.) then
		do r=1,Norb
		do q=1,Norb
		do m=1,Norb
		do p=1,Norb
			If (s.eq.m) then
				Ersmp(II,JJ,r,q,m,p)=Ers(II,JJ,r,q)*Ers(II,JJ,m,p)-Ers(II,JJ,r,p)
				Ersmp(JJ,II,r,q,m,p)=-Ersmp(II,JJ,r,q,m,p)
			Else
				Ersmp(II,JJ,r,q,m,p)=Ers(II,JJ,r,q)*Ers(II,JJ,m,p)
				Ersmp(JJ,II,r,q,m,p)=-Ersmp(II,JJ,r,q,m,p)
			End if
		end do
		end do
		end do
		end do
	end if
end do
end do

!open (13,FILE='BiElec.dat',STATUS='UNKNOWN')

!do II=1,dimCSF
!do JJ=1,dimCSF
!	do r=1,Norb
!	do s=1,Norb
!		do m=1,Norb
!		do p=1,Norb
!			if (Ersmp(II,JJ,r,s,m,p).ne.0) then
!!			write (13,('(I4, I4, I4, I4, I4, I4, F8.3)')) II,JJ,r,s,m,p,Ersmp(II,JJ,r,s,m,p)
!			write (13,('(I4, I4, F8.3)')) II,JJ,Ersmp(II,JJ,r,s,m,p)
!		end if
!		end do
!		end do
!	end do
!	end do
!end do
!end do

!write(*,*)' Done.'

! close(13)
!close(12)

end subroutine
!*************************************************************************************************

!*************************************************************************************************
subroutine ReglesSlaterCondon(dimCSF,Norb,tabDRT,chemin,Ne,S,CheckCondon)
!Retrouver les occupations par les coeficients a,b,c de la table DRT
!Verifier les differences
!Retourner le booleen associe
!Le but n'est pas de calculer les determinants mais de savoir si les generateurs seront nuls
!*************************************************************************************************
integer(kind=int_4),dimension(:,:),intent(in)::tabDRT
integer(kind=int_4),dimension(:,:),intent(in)::chemin
 Integer(kind=int_4),Intent(in)::Ne,Norb
 Integer(kind=int_8),Intent(in)::dimCSF
 Real(kind=real_8), Intent(in)::S
logical,dimension(:,:),intent(out)::CheckCondon
logical,dimension(Norb,2)::OccupationConfigII
logical,dimension(Norb,2)::OccupationConfigJJ
integer(kind=int_4),dimension(dimCSF,Norb)::aII,bII,cII,aJJ,bJJ,cJJ
integer(kind=int_4) :: II,JJ,NDiff,i,NelecII,NelecJJ,w
integer(kind=int_4) :: DaII,DaJJ,DbII,DbJJ,DcII,DcJJ
Real(kind=real_8)::Nup
!Me marche pas pour le dernier noeud, il faudrait automatiser...
!ATTENTION AUX DERNIER NIVEAU, TEST APRES LE LOOP DE REMPLISSAGE

!Les OccupationConfig disent s'il y a un electron dans quels spinorbitales
!OC(:,1) est up et OC(:,2) est down
 CheckCondon=.TRUE.
 DaII=0.d0
 DaJJ=0.d0
 DbII=0.d0
 DbJJ=0.d0
 DcII=0.d0
 DcJJ=0.d0
 NelecII=0.d0
 NelecJJ=0.d0

!Lire les coefficients a,b,c a partir de la table DRT
do II=1,dimCSF
do JJ=1,dimCSF
	do i=1,Norb
		aII(II,i)=tabDRT(chemin(i,II),3)
		bII(II,i)=tabDRT(chemin(i,II),4)
		cII(II,i)=tabDRT(chemin(i,II),5)
		aJJ(JJ,i)=tabDRT(chemin(i,JJ),3)
		bJJ(JJ,i)=tabDRT(chemin(i,JJ),4)
		cJJ(JJ,i)=tabDRT(chemin(i,JJ),5)
	end do
end do
end do

!Verifier les configurations par SC
do II=1,dimCSF
do JJ=1,II
do i=1,Norb-1
	!Initialisation a chaque loop
	Ndiff=0.d0
	OccupationConfigII=.FALSE.
 	OccupationConfigJJ=.FALSE.

	!Lire les coeficients a,b,c de la table DRT
	DaII=aII(II,i+1)-aII(II,i)
	DbII=bII(II,i+1)-bII(II,i)
	DcII=cII(II,i+1)-cII(II,i)
	DaJJ=aJJ(JJ,i+1)-aJJ(JJ,i)
	DbJJ=bJJ(JJ,i+1)-bJJ(JJ,i)
	DcJJ=cJJ(JJ,i+1)-cJJ(JJ,i)

	!Identifier quels spinorbitales sont remplies
	if(DbII.eq.0) then
		if (DcII.eq.0) then

			OccupationConfigII(i,:)=.true.
			NelecII=NelecII+2.d0
		else if (DcII.eq.1) then

			OccupationConfigII(i,:)=.false.
	end if
	else if(DbII.eq.1) then

		OccupationConfigII(i,1)=.true.
		NelecII=NelecII+1.d0	
	else if(DbII.eq.(-1)) then

		OccupationConfigII(i,1)=.false.	
		NelecII=NelecII+1.d0	
	end if
	if(DbJJ.eq.0) then
		if (DcJJ.eq.0) then

			OccupationConfigJJ(i,:)=.true.
			NelecJJ=NelecJJ+2.d0
		else if (DcII.eq.1) then

			OccupationConfigJJ(i,:)=.false.
		end if
	else if(DbJJ.eq.1) then
		OccupationConfigJJ(i,1)=.true.
		NelecJJ=NelecJJ+1.d0
	else if(DbJJ.eq.(-1)) then	
		OccupationConfigJJ(i,2)=.true.
		NelecJJ=NelecJJ+1.d0	
	end if
end do
	if (NelecII.eq.Ne-2) then

		!Il manque 2 electrons
		OccupationConfigII(Norb,:)=.true.
	else if (NelecII.eq.Ne-1) then

		!Il manque 1 electron
		!Determiner le spin de lelectron manquant
		Nup=0.d0
		do w=1,Norb-1
			if (OccupationConfigII(w,1).eq..true.) then

				Nup=Nup+1.d0
			end if
		end do
		if (Nup.eq.Ne/2+S) then

			OccupationConfigII(Norb,2)=.true.
		else if (Nup.eq.Ne/2+S-1) then
		
			OccupationConfigII(Norb,1)=.true.
		end if
	end if

	if (NelecII.eq.Ne-2) then

		!Il manque 2 electrons
		OccupationConfigII(Norb,:)=.true.
	else if (NelecII.eq.Ne-1) then

		!Il manque 1 electron
		!Determiner le spin de lelectron manquant
		Nup=0.d0
		do w=1,Norb-1
			if (OccupationConfigII(w,1).eq..true.) then

				Nup=Nup+1
			end if
		end do
		if (Nup.eq.Ne/2+S) then

			OccupationConfigII(Norb,2)=.true.
		else if (Nup.eq.Ne/2+S-1) then

			OccupationConfigII(Norb,1)=.true.
		end if
	end if
	!Le schema d'occupation est maintenant complet
	!Compter les differences
	call CompterDiff(OccupationConfigII,OccupationConfigJJ,Norb,Ndiff)

	!LE SEUL CHECK EST S'IL Y A MOINS DE 2 DIFFERENCES
	!LES DIFFERENCES SONT COMPTEES EN DOUBLE
	if (Ndiff.le.4) then

		!2 differences ou moins, il est possible que le generateur soit non-nul
		CheckCondon(II,JJ)=.true.
	else if (Ndiff.gt.4) then

		!Plus de 2 differences, le generateur est forcement nul
		CheckCondon(II,JJ)=.false.
	end if			
end do
end do

!Copier le triangle du bas de la matrice CheckCondon sur la partie supperieure


end subroutine
!*************************************************************************************************
!*************************************************************************************************

!*************************************************************************************************
subroutine CompterDiff(OccupationConfigII,OccupationConfigJJ,Norb,Ndiff)
!*************************************************************************************************
integer(kind=int_4),intent(in)::Norb
integer(kind=int_4),intent(out)::NDiff
logical,dimension(Norb,2),intent(in)::OccupationConfigII,OccupationConfigJJ
integer::i,j

Ndiff=0.d0

!Boucle simple pour comparer les deux matrices logiques
do i=1,Norb
	do j=1,2
		if (OccupationConfigII(i,j).ne.OccupationConfigJJ(i,j)) then
			Ndiff=Ndiff+1.d0			
		end if
	end do
end do


end subroutine
!*************************************************************************************************
!*************************************************************************************************

!*************************************************************************************************
subroutine CheckOccup(tabDRT,chemin)
!Calculates the occupation of all orbitals, in every CSF
!*************************************************************************************************
integer(kind=int_4),dimension(:,:),intent(in)::tabDRT
integer(kind=int_4),dimension(:,:),intent(in)::chemin
integer(kind=int_4) :: II,i,w,Ne,Norb,dimSCF,NeII
integer(kind=int_4) :: DaII,DbII,DcII
integer(kind=int_4),allocatable,dimension(:) :: Occupation

 DaII=0.d0
 DbII=0.d0
 DcII=0.d0

 Norb=tabDRT(1,1)-1
 dimSCF=tabDRT(1,14)
 Ne=2*tabDRT(1,3)+tabDRT(1,4)

allocate(Occupation(Norb))

open(45,file='OccupCSF.dat',status='unknown',form='formatted')
do II=1,dimSCF
NeII=0
Occupation=0
do i=1,Norb-1
	DaII=tabDRT(chemin(i,II),3)-tabDRT(chemin(i+1,II),3)
	DbII=tabDRT(chemin(i,II),4)-tabDRT(chemin(i+1,II),4)
	DcII=tabDRT(chemin(i,II),5)-tabDRT(chemin(i+1,II),5)
	if(DaII .eq. 1 .and. DbII .eq. 0 .and. DcII .eq. 0) then
		Occupation(Norb-i+1)=2
		NeII=NeII+2

	else if(DaII .eq. 0 .and. DbII .eq. 1 .and. DcII .eq. 0) then
		Occupation(Norb-i+1)=1
		NeII=NeII+1

	else if(DaII .eq. 1 .and. DbII .eq. -1 .and. DcII .eq. 1) then
		Occupation(Norb-i+1)=1
		NeII=NeII+1

	else if(DaII .eq. 0 .and. DbII .eq. 0 .and. DcII .eq. 1) then
		Occupation(Norb-i+1)=0
	end if
end do
Occupation(1)=(Ne-NeII)
write(45,*)II
write(45,*)Occupation
write(45,*)
end do

deallocate(Occupation)
 close(45)

end subroutine
!*************************************************************************************************
!*************************************************************************************************
!*************************************************************************************************
Subroutine get_CSF_occupation(EE)
!*************************************************************************************************

Real(kind=real_8),dimension(:,:,:,:),allocatable,intent(in)	::EE
Integer(kind=int_4)						::i,j,norb,ncsf

norb=size(EE(1,1,:,1)) 
ncsf=size(EE(:,1,1,1))

open(50,file='Occupation.dat',status='replace')

Do i=1,ncsf
	Write(50,*)'SCF #',i,' occupation of :'
	Write(50,'(<norb>F10.2)') (EE(i,i,j,j),j=1,norb)
	Write(50,*)
end do

close(50)

end subroutine

!*************************************************************************************************
!**********************************************************************************************


END MODULE nouveaudrt
