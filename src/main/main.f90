!Version * (fully Gaussian)
! Medys is a work done by 
! TT Nguyen-Dang Group at Chemistry Department of 
! Laval University
! Contributor are
! T.T. Nguyen-Dang
! F. Dion
! E. Blanchette
! J. Viau-Trudel
! A. Sainjon
!J.N. Vigneau
!Last edit: Jean-Nicolas Vigneau
!Date: 28/06/17
! Comments edited by TT Nguyen-Dang, (10/10/2015)
!
! Modifications  made starting 15/09/2020
! Restructured by TT Nguyen-Dang, to avoid using an explicit PW continuum orbital basis
! Fully Gaussian
!

program main

use variables
use basics
use string_module,only:string
use read_input
use grille
use qp_integrals
use read_columbus
use read_gamess
use nouveaudrt
use observable
!use Transformation
use Tests
use dynamique_new
!use matrix_temp_highlvl
use matrix_rearrangement_module

implicit none

 
!//////////////////////////////////////////////////////////////////////////////////////////////////
!Variable and array type declarations    
!//////////////////////////////////////////////////////////////////////////////////////////////////
! 
 Integer:: i,j,k,l,w,x,y,z
 Integer(kind=int_4) 					:: norb_lect,inR,nq,orb_min,orb_max,maxprim
 Character(len=32) 					:: sortie, nom
 Character(len=5) 					:: nom2
 Logical 						:: prop,ortho,opw

!for final storage of 1e, 2e integrals (over active+FC orbitals only)
 Real(kind=real_8),Allocatable,Dimension(:,:) 		:: SO,T1,V1,mu,lcao
 Real(kind=real_8),Allocatable,Dimension(:,:,:,:) 	:: vee,vee2

!temporary storage of integrals and LCAO MO's coeff (over full orbital basis)
 Real(kind=real_8),Allocatable,Dimension(:,:) 		:: SO_temp,T1_temp,V1_temp,mu_temp,mux_temp,muy_temp,muz_temp,lcao_temp
 Real(kind=real_8),Allocatable,Dimension(:,:,:,:) 	:: vee_temp,vee2_temp
! Real(kind=real_8),allocatable,dimension(:,:,:) 	:: gauss

!for DRT construction
 Real(kind=real_8) 					:: S
 Integer(kind=int_8) 					:: lDRT
 Integer(kind=int_4) 					:: FC
 Integer(kind=int_4), Allocatable,Dimension(:,:) 	:: tabDRT,tabDRT_tmp,chemin,chemin_tmp

!Hamiltonian matrix components in CSF basis
 Real(kind=real_8), Allocatable,Dimension(:,:,:,:) 	:: Ers
 Real(kind=real_8), Allocatable,Dimension(:,:,:,:,:,:) 	:: Ersmp
 Real(kind=real_8), Allocatable,Dimension(:,:) 		:: matrice,mux,muy,muz
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muEPS,hA,hA_temp
!!!! 
Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muEPS_Sq 
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: hbarB
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: hB,hB_temp
 Complex(kind=comp_16), Allocatable,Dimension(:) 	:: fctQ
 Complex(kind=comp_16), Allocatable,Dimension(:,:,:) 	:: fctP,Upp  ! FctP(dimP,orb_Q,nt) corresponds to M(J,r,t_n-t_l) of ttnd's notes, with n=l+1
!!!!
Complex(kind=comp_16),Allocatable,Dimension(:,:)	::Hqp ! may not be needed ??? (to be removed?)
Complex(kind=comp_16),Allocatable,Dimension(:,:)	::Hpqmain ! may not be needed ??? (to be removed?)
Complex(kind=comp_16),Allocatable,Dimension(:,:)	::eigenVectH0pp
Real(kind=real_8),Allocatable,Dimension(:,:)		::H0qq,MUqq
Real(kind=real_8),Allocatable,Dimension(:,:)		::H0pp,MUpp

! For bound to free (coupling) integrals
!
!! changes made on 28/10/2015: 3D k-grid -> 1D vector (will be of dimension nk)
! ****** all this section may be removed 
! Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: ski_temp
! Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: ski,skiCon,coeffGS,psipsi,SuperSki
! Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muxki,muyki,muzki
! Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muxki_temp,muyki_temp,muzki_temp
! Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muB_EPS,muki_EPS
! ********
! Integer(kind=int_4), Allocatable,Dimension(:,:) 	:: matPrim,AOrbDef ! check where those are used :  nowhere but allocated
! Real(kind=real_8), Allocatable,Dimension(:,:) 		:: CenterDef,MapCenter !   used  nowhere. CenterDef allocated
 Real(kind=real_8), Allocatable,Dimension(:,:)		:: ZET,ETA
 Real(kind=real_8), Allocatable,Dimension(:,:) 		:: coord
 Integer(kind=int_4), Allocatable,Dimension(:) 		:: ICONU,LMNP1
 Integer(kind=int_4), Allocatable,Dimension(:) 		:: NF
 Real(kind=real_8), Dimension(3)              		:: eps,kvec
 Real(kind=real_8), Allocatable,Dimension(:) 	  	:: znuc
 Integer(kind=int_4), Allocatable,Dimension(:) 		:: CMap
 Integer(kind=int_4) 					:: NS,NCONS,Iflg1,Iflg2
 Complex(kind=comp_16) 					:: idt

!plane wave basis parameters  ****** now used only for e_momentum spectra (post-dynamics) calculations
 Integer(kind=int_4), Dimension(3)			:: ni, nki
 Real(kind=real_8), Dimension(3)			:: imax, imin, dki, di, kmin
! ********

!For field and wf storage
 Real(kind=real_8) 					:: theta,phy!,delta
 Complex(kind=comp_16), Allocatable,Dimension(:)	:: sauvFctQ,eMomentum_gam  !to be checked
 Complex(kind=comp_16), Allocatable,Dimension(:,:)	:: sauvFctP,sauvHbarB !to be checked
 Real(kind=real_8), Allocatable,Dimension(:)		:: matA 
 !Real(kind=real_8),allocatable,dimension(:,:):: matAlpha, matPhi

! Primitive basis specs
!  
Real(kind=real_8),allocatable,dimension(:,:):: cEta, lcCG_temp, lcCG, prim_center
Real(kind=real_8),allocatable,dimension(:)::prim_EXPO, cZeta,sauvChamp
Integer(kind=int_4),allocatable,dimension(:,:)::prim_lmn
Integer(kind=int_4),allocatable,dimension(:,:):: lmn_vec
Integer(kind=int_4),dimension(3)::lp
Integer(kind=int_4),dimension(6)::ld
!Integer(kind=int_4) ::totPrimCount !declaration in variables
! timers
Real(kind=real_8):: timer1,timer2

Integer(kind=int_8) ::totmemallocated
!Integer(kind=int_4) ::k1,k2,k3,ki,QCsoft
Integer(kind=int_4) :: QCsoft

!  
integer(kind=int_4)::ii,jj,kk,ll,io_stat_vee2,io_stat_vee !i,j,k

! MOs reordering  
type(matrix_rearrangement)::reordering
type(string)::input_name
Character(len=32) 					:: argument
character(8)  :: date
character(10) :: time
character(5)  :: zone
integer,dimension(8) :: values
!integer(kind=int_4),dimension(:),allocatable::order_modif_vector
! TODO : add version as variable and add citation
DO i = 1, iargc()
  CALL getarg(i, argument)
  if ((argument == '-v') .or. (argument == '--version')) THEN
    WRITE (*,*) 'MEDYS V *  '
    stop
  end if
END DO
input_name='input'

!//////////////////////////////////////////////////////////////////////////////////////////////////
!Read in MODEL PARAMETERS, in particular:
!Ne=electron number, Norb_cat=number of active MO, Norb_sym=number of symbolic continuum MO
!FC=number of frozen MOs, S=spin q.number. 
!//////////////////////////////////////////////////////////////////////////////////////////////////

open(49,file='status.dat',status='replace',form='formatted')
write(49,*) "running"
close(49)
call param_debut(basis_file,nq,inR,Ne,Norb_cat,Norb_sym,FC,restriction,S,job,ni,imin,imax,di,dki,kmin,&
E0,omega,delta,theta,phy,pulsed,nper,pdt,tmin,nt,QCsoft,lp,ld,Iflg1,Iflg2)
!  subroutine param_debut(basis_file,nq,inR,Ne,Norb_cat,Norb_sym,FC,restriction,S,job,ni,imin,imax,di,E0&
!,omega,delta,theta,phy,pulsed,nper,pdt,tmin,nt,QCsoft,lp,ld)
 !TODO use correct value of pi ...
 delta=delta*pi
 

open(50,file='Parameters.dat',status='old',form='formatted',position='append')
Write(50,*)
!Write(50,'(A19)')'Time to treat input',timer2-timer1,'seconds'
!Write(50,*)
!
! norb_lect=number of MO read in from a QChem prior calculations (here with Columbus)
!QCsoft=0 (use Columbus data) is default value
!
write(*,*) "QCsoft=",QCsoft
select case(QCsoft)
case(0)
 call read_dimensions(norb_lect,basis_file)
case(1)
 call read_dimensions_gamess(norb_lect,basis_file,NCONS,maxprim,NS)
end select
 !!!!!!!!!!!!!!!HERE
orb_Q=FC+Norb_cat
Norb=orb_Q+Norb_sym

write(*,*) "Norb totales	", Norb 
write(*,*) "Norb liees 	", orb_Q
write(*,*) "Norb lues	", norb_lect
write(*,*)
!//////////////////////////////////////////////////////////////////////////////////////////////////
! read in informations on basis set used in QChem prior calculations (with Columbus-from argosls.sp)
!
! NCONS= number of contraction sets, NS=number of symmetry-distinct atoms (number of atoms in no-sym run)
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!select case(QCsoft)
!case(0)
 call scan_dimension_argos(basis_file,NCONS,maxprim,NS)
!case(1)
! call scan_dimension_gamess(basis_file,NCONS,maxprim,NS)
!end select
!write(*,*)"TTTTTTTEEEEEEEEESSSSSSSSSSSSSSTTTTTTTTTTTTTTTTTTTTT NCONS=",NCONS

! call lecture_NS(NS,basis_file,NCONS)- redundant call
allocate(ETA(maxprim,NCONS),ZET(maxprim,NCONS),ICONU(NCONS),NF(NS),LMNP1(NCONS),coord(NS,3),znuc(NS),CMap(NCONS))
!
!  maxprim=largest number of primitives among the NCONS contraction sets, ICONU(I)=number of primitives in Ith contraction set
!  NF(J)= number of AO sets for Jth type of symmetry-distinct atom, LMNP1(I)=value of (l+m+n+1) for Ith contraction set
!  coord(J,3)= position vector of the Jth atom of  symmetry-distinct type
!

! 
!memory usage calulation in bytes
totmemallocated=maxprim*NCONS*8*2 ! FOR ETA AND ZET
totmemallocated=totmemallocated+NCONS*4+NS*4+NCONS*4+3*NS*8+NS*8+4*NCONS
totmemallocated=totmemallocated+3*nPrim*4+2*6*nPrim*8+3*orb_Q*4+2*NS*8
write(*,*) 'total memory allocated on line 175 of main =',totmemallocated,' bytes'
!select case(QCsoft)
!case(0)
call lecture_gaus(norb_lect,basis_file,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord)
!case(1)
! write(*,*) "NCONS=",ncons
! call lecture_gaus_gamess(norb_lect,basis_file,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord)
!end select

Write(50,*)
Write(50,'(A)')"Read basis set and parameters from columbus: Done"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The code for reading information from argosls.sp is INCORRECT for COLUMBUS calculations with symmetry 
! scan_dimension_argos and lecture_gauss to be revised to this end.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Write(50,*)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
! k-grid size. nk=number of grid points of (1-3D) k-grid.
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 nk=1
do i=1,3
	nki(i)=ni(i)
	nk=nk*nki(i)
end do

! time step
idt=dcmplx(0.d0,-1.d0*pdt)


!///////////////////////////////////////////////////////////////
!
! INTEGRALS READING AND/OR CONSTRUCTION
!Read in 1e, 2e integrals (in full AO basis) and LCAO matrix  
!/////////////////////////////////////////////////////////////// 

 write(*,*) "test on line 189 of main.f90",norb_lect
allocate(SO_temp(norb_lect,norb_lect))
allocate(T1_temp(norb_lect,norb_lect))
allocate(V1_temp(norb_lect,norb_lect))
allocate(mux_temp(norb_lect,norb_lect))
allocate(muy_temp(norb_lect,norb_lect))
allocate(muz_temp(norb_lect,norb_lect))
allocate(lcao(norb_lect,norb_lect))
allocate(vee_temp(norb_lect,norb_lect,norb_lect,norb_lect))
allocate(vee2_temp(norb_lect,norb_lect,norb_lect,norb_lect))
totmemallocated=(totmemallocated+norb_lect*norb_lect*7+norb_lect**3*2)*8
write(*,*) 'total memory allocated on line 216 of main =',totmemallocated,' bytes'
!
SO_temp=0.d0
v1_temp=0.d0
T1_temp=0.d0
mux_temp=0.d0
muy_temp=0.d0
muz_temp=0.d0
vee_temp=0.d0
vee2_temp=0.d0
lcao=0.d0

Write(50,*)
Write(50,'(A)')"Reading integrals" 
Write(50,*)

!actual reading of integrals (over AO) and LCAO coefficients for all MO

!select case(QCsoft)
!case(0)
 call read_integrals(SO_temp,T1_temp,V1_temp,muz_temp,mux_temp,muy_temp,vee2_temp,norb_lect,basis_file)
 call read_coeff_lcao(lcao,norb_lect,basis_file)
!case(1)
! call read_integrals_gamess(SO_temp,T1_temp,V1_temp,muz_temp,mux_temp,muy_temp,vee2_temp,norb_lect,basis_file)
! call read_coeff_lcao_gamess(lcao,norb_lect,basis_file)
!stop
!end select
write(99,*) "[MO]"
do i=1,norb_lect
write(99,*) " Sym= A1"
write(99,*) " Ene= 0.000000"
write(99,*) " Spin= Alpha"
write(99,*) " Occup= 1"

 do j=1,norb_lect
   if (j.lt.10) then
     write(99,('(I1,2x,F13.10)')) j,lcao(j,i)
   else if (j.lt.100) then
     write(99,('(I2,2x,F13.10)')) j,lcao(j,i)
   else
     write(99,('(I3,2x,F13.10)')) j,lcao(j,i)
   end if
 end do
end do
Write(50,*)
Write(50,'(A)')'Read integrals: Done '

!do i=1,norb_lect
!    do j=1,ICONU(i)
!        WRITE(200+i,'(f16.10)')ETAtmp(j,i)*dexp(ZET(j,i)))
!    end do
!end do


! TODO possible to create molden file?
  write(50,*) 'norb_lect',norb_lect
do i=1,norb_lect
  write(50,*) ' Sym=A1'
  write(50,*) ' Ene=        -0.0000000000'
  write(50,*) ' Spin= Alpha'
  write(50,*) ' Occup=   1'

  do j=1,norb_lect
	  Write(50,'(I3,2x,F13.10)')j,lcao(j,i) 
  end do
end do

Write(50,*)
Write(50,'(A)')'Read LCAO Coefficients: Done' 

!TODO, this was commeted since I dont have .nml files
! uncommented, but how to construct this file??
! MO Reordering  
 call reordering%matrix_1D_reordering(lcao,input_name)

 write(100,*) "[MO]"
do i=1,norb_lect
write(100,*) " Sym= A1"
write(100,*) " Ene= 0.000000"
write(100,*) " Spin= Alpha"
write(100,*) " Occup= 1"

 do j=1,norb_lect
   if (j.lt.10) then
     write(100,('(I1,2x,F13.10)')) j,lcao(j,i)
   else if (j.lt.100) then
     write(100,('(I2,2x,F13.10)')) j,lcao(j,i)
   else
     write(100,('(I3,2x,F13.10)')) j,lcao(j,i)
   end if
 end do
end do

!
!  Integral (AO->MO) transformations. In full MO/AO basis.
    call aointTOmointn5(SO_temp,T1_temp,V1_temp,mux_temp,muy_temp,muz_temp,vee2_temp,lcao,norb_lect,vee_temp)

Write(50,*)
Write(50,'(A)')'Change of basis (from AO to MO) done'





!**************************************************************************************************************
!**************************************************************************************************************

100 do i=1,norb_lect
    do j=1,norb_lect
	if (dabs(mux_temp(i,j)) .lt. 1.d-10) then
		mux_temp(i,j)=0.d0
	end if

	if (dabs(muy_temp(i,j)) .lt. 1.d-10) then
		muy_temp(i,j)=0.d0
	end if

	if (dabs(muz_temp(i,j)) .lt. 1.d-10) then
		muz_temp(i,j)=0.d0
	end if

    end do	
    !write(30,'(<nk>F34.25)')(rWavePacket_obs(i),i=1,nk)
    write(150,'(<norb_lect>F34.25)') (mux_temp(i,j),j=1,norb_lect)
    write(151,'(<norb_lect>F34.25)') (muy_temp(i,j),j=1,norb_lect)
    write(152,'(<norb_lect>F34.25)') (muz_temp(i,j),j=1,norb_lect)
    end do

!
! call new sub JoinPrimInfo.  Need calculate total number of primitives first. These are presently for nosym QC results from COLUMBUS
! 
call PrimCount(NS,NF,ICONU, totPrimCount,LMNP1)
allocate(cEta(norb_lect,totPrimCount),cZeta(totPrimCount),lmn_vec(totPrimCount,3), prim_center(totPrimCount,3),lcCG_temp(norb_lect,totPrimCount))
call Join_Prim_Info(NS,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord,cEta,cZeta,lmn_vec,prim_center,lcao,lcCG_temp,lp,ld,norb_lect)

 

!//////////////////////////////////////////////////////////////////////////////////////////////////
!Restriction of integrals to active MO (+FC). New integral matrices dimensioned to orb_Q= Norb_cat + FC
!//////////////////////////////////////////////////////////////////////////////////////////////////
! 
!totmemallocated=totmemallocated+orb_Q*orb_Q*16*3+nk*orb_Q*16*3+orb_Q**4*8+orb_Q**2*8+nk*orb_Q*2*16
!write(*,*) 'total memory allocated on line 317 of main =',totmemallocated,' bytes'
allocate(mux(orb_Q,orb_Q))
allocate(muy(orb_Q,orb_Q))
allocate(muz(orb_Q,orb_Q))
allocate(lcCG(orb_Q,totPrimCount))
allocate(vee2(orb_Q,orb_Q,orb_Q,orb_Q),hA(orb_Q,orb_Q))

write(*,*) "test on line 318 after allocation"
!allocate(hA_temp(norb_lect,norb_lect))

!Fill in new itegral matrices
 call cpu_time(timer1)
!hA_temp=T1_temp+V1_temp

do i=1,orb_Q

	do j=1,orb_Q
		hA(i,j)=T1_temp(i,j)+V1_temp(i,j)
		write(*,*) 'i=',i,'j=',j,'hA(i,j)=',hA(i,j)
		mux(i,j)=mux_temp(i,j)
		muy(i,j)=muy_temp(i,j)
		muz(i,j)=muz_temp(i,j)
			do l=1,orb_Q
		do k=1,orb_Q
				vee2(i,j,k,l)=vee_temp(i,j,k,l)
				!write(*,*) 'i=',i,'j=',j,'k=',k,'l=',l,'vee=',vee2(i,j,k,l)
			end do
		end do
		!write(*,*) 'i=',i,'j=',j,'muz(i,j)=',muz(i,j)
	end do
  !      lcCG(i,:)=lcCG_temp(i,:)
end do
lcCG=lcCG_temp(1:orb_Q,:)
Write(50,*)
Write(50,'(A)')'Restrictions of integrals to active orbitals: Done '


! totmemallocated=totmemallocated+(nk*orb_Q*2+orb_Q**2)*16
! write(*,*) 'total memory allocated on line 350 of main =',totmemallocated,' bytes'

write(*,*) "test on line 352 after allocate"
! mu*eps,  
! eps=polarization vector, in polar form
 
eps(1)=EX(1.d0,theta,phy)
eps(2)=EY(1.d0,theta,phy)
eps(3)=EZ(1.d0,theta) 

muEPS =eps(1)*mux +eps(2)*muy +eps(3)*muz 

deallocate(SO_temp,v1_temp,T1_temp,mux_temp,muy_temp,muz_temp,vee_temp,vee2_temp)


!!//////////////////////////////////////////////////////////////////////////////////////////////////
!  CSF SPACE
!  DRT construction
!!//////////////////////////////////////////////////////////////////////////////////////////////////

 call dimension_DRT(S,Ne,Norb,lDRT,dimCSF)
 
allocate(tabDRT(lDRT,14),chemin(Norb+1,dimCSF))

totmemallocated=totmemallocated+(14*lDRT+(Norb+1)*dimCSF)*4
write(*,*)"allocating DRT need ",(14*lDRT+(Norb+1)*dimCSF)*4," bytes"
write(*,*)"total memory used = ",totmemallocated," bytes"
 call gen_drt(S,Ne,orb_Q,Norb_sym,restriction,tabDRT,chemin,lDRT,dimCSF,FC,dimQ)

allocate(tabDRT_tmp(lDRT,14),chemin_tmp(Norb+1,dimCSF))
 tabDRT_tmp=tabDRT
 chemin_tmp=chemin
deallocate(tabDRT,chemin)
allocate(tabDRT(lDRT,14),chemin(Norb+1,dimCSF))
tabDRT=tabDRT_tmp
 chemin=chemin_tmp
deallocate(chemin_tmp,tabDRT_tmp)

!! Exhibit constructed DRT
do i=1,lDRT
write(*,'(2i5,4x,3i5,4x,4i5,4x,4i6,4x,i6)')(tabDRT(i,j),j=1,14)
end do
!!//////////////////////////////////////////////////////////////////////////////////////////////////
!! Generators 
!!  Ers(I,J,r,s), Ersmp(I,J,r,s,m,p) (r,s=MO indices, I,J=CSF indices)
!! matrices in CSF basis
!!//////////////////////////////////////////////////////////////////////////////////////////////////
dimP=dimCSF-dimQ
allocate(Ers(dimCSF,dimCSF,Norb,Norb))
allocate(Ersmp(dimCSF,dimCSF,Norb,Norb,Norb,Norb))
Ers=0.d0
!********************************************************************************************************************
!actual calculations of the generators
write(*,*)
write(*,*)' Calcul des generateurs mono-electroniques...'
!1e 
 call gen_gen(tabDRT,chemin,Norb,dimCSF,Ers)
!Writes the occupation of the orbitals for each CSF
 call get_CSF_occupation(Ers)

write(*,*)' Calcul des generateurs bi-electroniques...'
!2e
 call gen_bielect(dimCSF,Norb,Ersmp,Ers,tabDRT,chemin,Ne,S)
 call calculEE(Ers,Ersmp)

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!AJOUT JN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

write(*,*)' Mono-électroniques:'
write(*,*)'dimCSF=',dimCSF
write(*,*)'dimQ=',dimQ
!do i=1,Norb
!  do j=1,Norb
!    write(*,*)'i,j=',i,j
!    do II=1,dimCSF
!      write(*,"(5E25.16)") (Ers(II,JJ,i,j),JJ=1,dimCSF)
!    end do
!  end do
!end do
!write(*,*)' Test de Ersmp'
!do w=1,3
!  do x=1,3
!    do y=1,3
!      do z=1,3
!        write(*,*) 'Ersmp(I,J,',w,',',x,',',y,',',z,')'
!        do I=1,dimCSF
!          write(*,"(5E25.16)") (Ersmp(I,J,w,x,y,z),J=1,dimCSF)
!        end do
!      end do
!    end do
!  end do
!end do

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!FIN AJOUT JN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



!//////////////////////////////////////////////////////////////////////////////////////////////////
! START DYNAMICS
!Call dyna=actual dynamical calculation (Propagation of TD wavefunction through TD CI coefficients)
!////////////////////////////////////////////////////////////////////////////////////////////////
write(*,*) "dimQ", dimQ
!JN allocate(sauvChamp(nt+1),fctQ(dimQ),fctP(dimP,orb_Q, nt), Upp(nt,dimP,dimP), matA(nt),eMomentum_gam(dimP)) ! ,matAlpha(nt,nt),matPhi(nt,nt)
allocate(sauvChamp(nt+1),fctQ(dimQ),fctP(dimP,orb_Q, nt+1), Upp(nt,dimP,dimP), matA(nt+1),eMomentum_gam(dimP)) !,matAlpha(nt,nt),matPhi(nt,nt) !JN
!allocate(Hqp(dimQ,dimP*nk),Hpq(dimP*nk,dimQ))

!todo: allocate muEPS_Sq in new_dyna and delete variable in main
allocate(H0qq(dimQ,dimQ),H0pp(dimP,dimP),MUqq(dimQ,dimQ),MUpp(dimP,dimP),eigenVectH0pp(dimP,dimP),muEPS_Sq(orb_Q,orb_Q)) !Hqppq(dimQ,dimQ),

matA=0.d0
!matAlpha=0.d0
! matPhi=0.d0

Write(*,*)'Post allocation'
!
! WP Initialisation
!
fctQ=dcmplx(0.d0,0.d0)
fctP=dcmplx(0.d0,0.d0)
fctQ(1)=dcmplx(1.d0,0.d0)
Write(*,*)'Hamiltonians'
!
! Hd_CSF calculates H0_QQ, H0_PP, mu_QQ, mu_PP  in Q(P)-space CSF basis 
!        from integrals hA, muEPS, Vee and matrices Ers, ee in CSF basis
!!!!!!
! Hnd_CSF calculate t-independent part of H_QP, H_PQ, and H_QP*H_PQ     !!!!  Replaced  
!         from integrals muB_EPS (1e only) and Ers matrices in CSF basis
!!!!!!!
!

call Hd_CSF(hA,muEPS,Vee2,Ers,Ersmp,H0qq,H0pp,MUqq,MUpp)
!
!call Hnd_CSF(muB_EPS,Ers,Hqp,Hpq,Hqppq,dki)  !!!!  Replaced by (new sub) Hqppq_CSF called and constructed in dynamiquef90 
!
!
!
open(18,file='eigenvectors.dat',status='unknown',form='formatted')
Write(18,'(A)')'the first column is the Real part of each coefficients, the second one is the Imaginary part'
Write(18,*)
Write(18,'(A)')'SUBSPACE P'
Write(18,'(I3,A)')dimP,' eigenvectors of the P subspace'
! 
 call EigenVectors(H0pp,dimP,eigenVectH0pp)
Write(18,*)
Write(18,'(A)')'SUBSPACE Q'
Write(18,'(I3,A)')dimQ,' eigenvectors of the Q subspace'
 call etatInitial(H0qq,dimQ,fctQ)

 close(18)
write(*,*)'Debut dynamique'

!
call new_dyna(fctQ,fctP,Upp,muEPS,muEPS_Sq,sauvChamp,matA,eps,cEta,cZeta,&
lmn_vec,prim_center,lcCG,Ers,H0qq,H0pp,MUqq,MUpp) !,matAlpha,matPhi
 
!    Write e-momentum spectra for each ionic channel if requested (iflg2=1). Others possible value of iflg2 are 0 (no spectrum generated), 2= spectra generated for selected J (to be done). Those are presently unavailable.
write(*,*) "Iflg2", Iflg2
select case(Iflg2)

case (1)

do j=1,dimP
         write (nom2,'(I5.5)') j
         open(110+j,file="spectre_k_canal"// ADJUSTL(nom2) //".dat",status='replace',form='formatted') ! attention: numéro d'unité dépassant 100 
!        open(110+j,file="spectre2dx_canal"// ADJUSTL(nom2) //".dat",status='replace',form='formatted')
         write(*,*)"J=,", j,"fctP(J,:,2)=" , fctP(J,:,2) 
end do
!write(*,*) "dimP = ",dimP
!write(*,*) "nki(1) = ",nki(1)
!write(*,*) "nki(2) = ",nki(2)
!write(*,*) "nki(3) = ",nki(3)

write(*,*)"test de lcCG, prim_center"

do i=1,totPrimCount
write(*,*)"MO 1, CG number", i,"lcCG(1,i)=",lcCG(1,i),"Rc(i)=",prim_center(i,:)
enddo

write(*,*)"eps=",eps( 1), eps( 2), eps( 3)  

do i=1,nki(1)
!        write(*,*) "k(1)",i,"/",nki(1)
	kvec(1)=(i-1)*dki(1)+kmin(1)
	do j=1,nki(2)
!           write(*,*) "k(2)",j,"/",nki(2)
	   kvec(2)=(j-1)*dki(2)+kmin(2)
	   do k=1,nki(3)
!             write(*,*) "k(3)",k,"/",nki(3)
             kvec(3)=(k-1)*dki(3)+kmin(3)
             call Get_Observable_eMomentum_map(kvec,fctP,Upp,eMomentum_gam,lcCG,muEPS,matA,eps,cEta,cZeta,lmn_vec,prim_center,tn) !,matAlpha,matPhi
!             write(*,*) "dimP = ",dimP
             do jj=1,dimP
!                write(*,*) jj,"/",dimP
                write(110+jj, '(3X,F10.5,3X,F10.5,3X,F10.5,3X,ES24.14)') kvec(1), kvec(2), kvec(3), cdabs(eMomentum_gam(jj))**2
!                write(*, '(3X,F10.5,3X,F10.5,3X,F10.5,3X,ES24.14)') kvec(1), kvec(2), kvec(3), cdabs(eMomentum_gam(jj))**2
           enddo
        enddo
         enddo
enddo
     
case(0)

end select

! A concevoir et faire: ecriture de fctP, matA, matAlpha, matPhi, cEta, cZeta, lmn_vec, prim_center,lcCG dans fichiers
! (a) pour calculs ultérieurs d'observables (b) pour restart
! On aurait besoin d'écrire UPP(nt) dans un fichier aussi 
!
 open(49,file='status.dat',status='replace',form='formatted')
write(49,*) "Medys terminated at "
! using keyword arguments
call date_and_time(date,time,zone,values)
call date_and_time(DATE=date,ZONE=zone)
call date_and_time(TIME=time)
call date_and_time(VALUES=values)
print '(a,2x,a,2x,a)', date, time, zone
print '(8i5)', values
close(49)
call cpu_time(timer2)
write(*,*) 'CPU TIME',timer2-timer1,' seconds'
end program

!    contains
!
!    ! Currently, this subroutine write Ecin for xz or yz plan
!    ! We must think a lot for using this subroutine for writting general 2d observable
!    subroutine write2d(Obs,dki,nki,kmin,dimP)
!    !Ecriture fichiers -> MatLab/Mathematica
!    !****************************************************************************************
!    !****************************************************************************************
!    Complex(kind=comp_16),Dimension(:),Intent(in)  ::Obs
!    Real(kind=real_8),Dimension(3),Intent(in)        ::dki,kmin
!    Integer(kind=int_4),Dimension(3),Intent(in)     ::nki
!    Integer(kind=int_4),Intent(in)              ::dimP
!
!    Complex(kind=comp_16),allocatable,dimension(:,:)    ::fctPSomme
!    Real(kind=real_8),allocatable,dimension(:,:)        ::fctPSomme2
!    Complex(kind=comp_16),allocatable,dimension(:)        ::fctPSumcoerr
!
!    Real(kind=real_8),allocatable,dimension(:)      ::Mx,My,Mz
!    Real(kind=real_8),allocatable,dimension(:,:)        ::E
!    Integer                         ::i,j,ik,kx,ky,kz,ktot
!    Real(kind=real_8)                   ::kk,kkx,kky,kkz
!    Character(len=5)                    ::nom
!    Character(len=23)                   ::canal
!    CHARACTER(LEN=20) :: FMT
!
!    !write(*,*)'start obs'
!    open(52,file='PSomme.dat',status='replace',form='formatted')
!    do j=1,dimP
!        write (nom,'(I5.5)') j
!        open(90+j,file="spectre2dy_canal"// ADJUSTL(nom) //".dat",status='replace',form='formatted')
!        open(110+j,file="spectre2dx_canal"// ADJUSTL(nom) //".dat",status='replace',form='formatted')
!
!    end do
!    open(70,file="spectre_tot_coer.dat",status='unknown',form='formatted')
!    open(83,file="kx.dat",status='unknown',form='formatted')
!    open(80,file="ky.dat",status='unknown',form='formatted')
!    open(81,file="kz.dat",status='unknown',form='formatted')
!    !open(61,file='spectre_canal_1.dat',status='unknown',form='formatted')
!    !open(62,file='spectre_canal_2.dat',status='unknown',form='formatted')
!    !open(63,file='spectre_canal_3.dat',status='unknown',form='formatted')
!    ktot=nki(1)*nki(2)*nki(3)
!
!
!    E=0.d0
!
!
!    !write(*,*)kmin
!    !write(*,*)dki
!    !write(*,*)nki
!    WRITE(FMT,*) nki(2)
!    do kz=1,nki(3)
!      do ky=1,nki(2)
!        do kx=1,nki(1)
!          kkx=kmin(1)+real(kx)*dki(1)
!          kky=kmin(2)+real(ky)*dki(2)
!          kkz=kmin(3)+real(kz)*dki(3)
!          kk=(kkx*kkx + kky*kky + kkz*kkz)**(1.d0/2.d0)
!          ik=kx + nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1)
!          E(ik,1)=kk*kk/(2.d0)
!          E(ik,2)=kkx
!          E(ik,3)=kky
!          E(ik,4)=kkz
!        end do
!        if (nki(1).gt.1) then
!          write(82,'('// ADJUSTL(FMT) //'ES24.14)') obsv(i+ nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1),i=1,nki(2))
!          write(83,'('// ADJUSTL(FMT) //'ES24.14)') (E(i+ nki(1)*(ky-1) + nki(1)*nki(2)*(kz-1),2),i=1,nki(2))
!          write(81,'('// ADJUSTL(FMT) //'ES24.14)') (E(1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1),4),i=1,nki(2))
!
!          do j=1,dimP
!            write(110+j,'('// ADJUSTL(FMT) //'ES24.14)') (fctPSomme2(j,1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1)),i=1,nki(2))
!          end do
!
!        end if
!      end do
!      ! kz.dat
!      !!writting ky.dat and observable for one kz value and all ky value as colomn
!      if (nki(2).gt.1) then
!        write(70,'('// ADJUSTL(FMT) //'ES24.14)') (abs(fctPSumcoerr(1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1)))**(2.d0),i=1,nki(2))
!        write(80,'('// ADJUSTL(FMT) //'ES24.14)') (E(1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1),3),i=1,nki(2))
!        write(81,'('// ADJUSTL(FMT) //'ES24.14)') (E(1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1),4),i=1,nki(2))
!
!
!        do j=1,dimP
!          write(90+j,'('// ADJUSTL(FMT) //'ES24.14)') (fctPSomme2(j,1+nki(1)*(i-1) + nki(1)*nki(2)*(kz-1)),i=1,nki(2))
!        end do
!      end if
!
!    end do
!
!    do j=1,dimP
!        close(90+j)
!        close(110+j)
!    end do
!    close(70)
!    close(80)
!    close(81)
!    do j=1,dimP
!
!
!    write (nom,'(I5.5)') j
!
!    ! ajouter spectrek_canal ... avec E(ik,1)
!    canal='spectre_canal_'//trim(nom)//'.dat'
!
!    open(60,file=canal,status='unknown',form='formatted')
!        if (nki(1) .eq. 1) then
!        do i=1,ktot
!            write(52,*) real(fctPSomme(j,i)),aimag(fctPSomme(j,i))
!            write(60,*) E(i,3),E(i,4),fctPSomme2(j,i)
!        end do
!
!        else if (nki(2) .eq. 1) then
!        do i=1,ktot
!            write(52,*) real(fctPSomme(j,i)),aimag(fctPSomme(j,i))
!            write(60,*) E(i,2),E(i,4),fctPSomme2(j,i)
!        end do
!
!        else if (nki(3) .eq. 1) then
!        do i=1,ktot
!            write(52,*) real(fctPSomme(j,i)),aimag(fctPSomme(j,i))
!            write(60,*) E(i,2),E(i,3),fctPSomme2(j,i)
!        end do
!        endif
!
!    close(60)
!
!    end do
!
!     close(52)
!    ! close(61)
!    ! close(62)
!    ! close(63)
!
!
!    end subroutine
