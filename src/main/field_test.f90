!Version 1.0.1
! Field is a work done by 
! F. Dion
!

program fieldtest

use variables
use basics
use string_module,only:string
use read_input
use champ1
use grille
use qp_integrals
use read_columbus
use read_gamess
use nouveaudrt
use observable
!use Transformation
use Tests
use dynamique
!use matrix_temp_highlvl
use matrix_rearrangement_module


use math
use champ1
use mkl95_blas
use mkl95_lapack
use mkl95_precision
use NewSubs



implicit none

 
!//////////////////////////////////////////////////////////////////////////////////////////////////
!Variable and array type declarations    
!//////////////////////////////////////////////////////////////////////////////////////////////////
! 
 Integer:: i,j,k,l,w,x,y,z
 Integer(kind=int_4) 					:: norb_lect,inR,nq,orb_min,orb_max,maxprim,opt_GS
 Character(len=32) 					:: sortie
 Logical 						:: prop,ortho,opw

!for final storage of 1e, 2e integrals (over active+FC orbitals only)
 Real(kind=real_8),Allocatable,Dimension(:,:) 		:: SO,T1,V1,mu,lcao,t
 Real(kind=real_8),Allocatable,Dimension(:,:,:,:) 	:: vee,vee2

!temporary storage of integrals and LCAO MO's coeff (over full orbital basis)
 Real(kind=real_8),Allocatable,Dimension(:,:) 		:: SO_temp,T1_temp,V1_temp,mu_temp,mux_temp,muy_temp,muz_temp,lcao_temp
 Real(kind=real_8),Allocatable,Dimension(:,:,:,:) 	:: vee_temp,vee2_temp
 Real(kind=real_8),allocatable,dimension(:,:,:) 	:: gauss

!for DRT construction
 Real(kind=real_8) 					:: S
 Integer(kind=int_4) 					:: lDRT
 Integer(kind=int_4) 					:: FC
 Integer(kind=int_4), Allocatable,Dimension(:,:) 	:: tabDRT,tabDRT_tmp,chemin,chemin_tmp

!Hamiltonian matrix components in CSF basis
 Real(kind=real_8), Allocatable,Dimension(:,:,:,:) 	:: Ers
 Real(kind=real_8), Allocatable,Dimension(:,:,:,:,:,:) 	:: Ersmp
 Real(kind=real_8), Allocatable,Dimension(:,:) 		:: matrice,mux,muy,muz
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muEPS,hA,hA_temp
!!!!  
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: hbarB
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: hB,hB_temp
 Complex(kind=comp_16), Allocatable,Dimension(:) 	:: fctQ
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: fctP
!!!!
Complex(kind=comp_16),Allocatable,Dimension(:,:)	::Hqp
Complex(kind=comp_16),Allocatable,Dimension(:,:)	::Hpq
Complex(kind=comp_16),Allocatable,Dimension(:,:)	::Hqppq,eigenVectH0pp
Real(kind=real_8),Allocatable,Dimension(:,:)		::H0qq,MUqq
Real(kind=real_8),Allocatable,Dimension(:,:)		::H0pp,MUpp

! For bound to free (coupling) integrals
!
!! changes made on 28/10/2015: 3D k-grid -> 1D vector (will be of dimension nk)
!
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: ski_temp
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: ski,skiCon,coeffGS,psipsi,SuperSki
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muxki,muyki,muzki
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muxki_temp,muyki_temp,muzki_temp
 Complex(kind=comp_16), Allocatable,Dimension(:,:) 	:: muB_EPS,muki_EPS
! 
 Integer(kind=int_4), Allocatable,Dimension(:,:) 	:: matPrim,AOrbDef
 Real(kind=real_8), Allocatable,Dimension(:,:) 		:: CenterDef,MapCenter
 Real(kind=real_8), Allocatable,Dimension(:,:)		:: ZET,ETA
 Real(kind=real_8), Allocatable,Dimension(:,:) 		:: coord
 Integer(kind=int_4), Allocatable,Dimension(:) 		:: ICONU,LMNP1
 Integer(kind=int_4), Allocatable,Dimension(:) 		:: NF
 Real(kind=real_8), Dimension(3)              		:: eps
 Real(kind=real_8), Allocatable,Dimension(:) 	  	:: znuc
 Integer(kind=int_4), Allocatable,Dimension(:) 		:: CMap
 Integer(kind=int_4) 					:: NS,NCONS
 Complex(kind=comp_16) 					:: idt

!plane wave basis parameters
 Integer(kind=int_4), Dimension(3)			:: ni, nki
 Real(kind=real_8), Dimension(3)			:: imax, imin, dki, di, kmin

!For field and wf storage
 Real(kind=real_8) 					:: theta,phy,delta
 Complex(kind=comp_16), Allocatable,Dimension(:)	:: sauvFctQ
 Complex(kind=comp_16), Allocatable,Dimension(:,:)	:: sauvFctP,sauvHbarB
 Real(kind=real_8), Allocatable,Dimension(:)		:: sauvChamp

! Primitive basis specs
Real(kind=real_8),allocatable,dimension(:,:)::Mat,prim_center
Real(kind=real_8),allocatable,dimension(:)::prim_EXPO
Integer(kind=int_4),allocatable,dimension(:,:)::prim_lmn
Integer(kind=int_4),dimension(3)::lp
Integer(kind=int_4),dimension(6)::ld

! timers
Real(kind=real_8):: timer1,timer2

Integer(kind=int_8) ::totmemallocated
Integer(kind=int_4) ::k1,k2,k3,ki,QCsoft

!  
integer(kind=int_4)::ii,jj,kk,ll,io_stat_vee2,io_stat_vee

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
    WRITE (*,*) 'MEDYS V0.99 '
    stop
  end if
END DO
input_name='input'

!//////////////////////////////////////////////////////////////////////////////////////////////////
!Read in MODEL PARAMETERS, in particular:
!Ne=electron number, Norb_cat=number of active MO, Norb_sym=number of symbolic continuum MO
!FC=number of frozen MOs, S=spin q.number. 
!//////////////////////////////////////////////////////////////////////////////////////////////////

 call param_debut(basis_file,nq,inR,Ne,Norb_cat,Norb_sym,FC,restriction,S,job,ni,imin,imax,di,dki,kmin,&
E0,omega,delta,theta,phy,pulsed,nper,pdt,tmin,nt,opt_GS,QCsoft,lp,ld)

 !TODO use correct value of pi ...
 delta=delta*pi
 

allocate(sauvChamp(nt+1),fctQ(dimQ),fctP(dimP,nk))


open(51,file='champ.dat',status='unknown',form='formatted')
!////////////////////////////////////////////////////////
! start loop over time variable (tn), i.e. TIME-PROPAGATION
!////////////////////////////////////////////////////////

do tn=1,nt  ! begin time-loop
t=tn
       sauvChamp(tn+1)=champNint(delta)
!
        write(51,*)tmin+(tn*pdt),Int0(delta),sauvChamp(tn+1)
end do

end program fieldtest
