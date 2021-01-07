Module read_gamess

! read informations from columbus files
use basics
use blas95
use lapack95
use f95_precision
contains
subroutine scan_dimension_gamess(basis_file,NCONS,maxprim,NS)
 character(len=32),intent(in)::basis_file
Integer(kind=int_4) :: NCONS,maxprim,NS
 character(len=64)::read_file,numeroR

  character(len=100)::output,output2,output3            !specific line to read
 character(len=5)::atomsymb,iorb
 character(len=100)::output_name,output_name2,output_name3,system
integer(kind=int_4)::NGEN,NAORDS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT,i,j
integer(kind=int_4)::NDPT,IS
integer(kind=int_4),dimension(NCONS)::SetCenter ! center of the contraction set
integer(kind=int_4)::counter1,ll
integer(kind=int_4), allocatable,dimension(:)::ICSU,NRCR,NC,MCONS,IGCS,MCRS,IGEN
real(kind=real_8), allocatable,dimension(:,:)::ZETtmp,ETAtmp
integer(kind=int_4)::nprimmax
real(kind=real_8),allocatable,dimension(:)::X,Y,Z,AN
 character(len=3),allocatable,dimension(:)::MTYPE
real,allocatable,dimension(:)::CHG
 character::NEANT
 character(len=21),parameter::FMT1 ='(t1,d16.8,t1,d16.8)'
 character(len=13)::pointeur
integer::nligne,reste
Integer(kind=int_4) ::ierror
 CHARACTER(len=8),dimension(:),allocatable :: ATOMNAME,ATOMNAME2

 CHARACTER(len=100)::line1,line2
integer(kind=int_4) :: SHELLNUM,PRIMNUM,REALSHELLNUM
CHARACTER(len=2) :: OTYPE
real :: EXPONENT,CONTRACTION


allocate(ATOMNAME(2),ATOMNAME2(2),X(2),Y(2),Z(2),AN(2))
allocate(MCONS(10),IGCS(10))
ZET=0.d0
ETA=0.d0
coord=0.d0
ICONU=0.d0
NF=0.d0
LMNP1=0.d0
i=1
REALSHELLNUM=0

output_name='BASIS/'//trim(basis_file)//'/gamess1e.out'
write(*,*) "output_name=",output_name
 OPEN (UNIT=10,FILE=output_name,STATUS='old',form='formatted',IOSTAT=ierror)
 write(*,*) 'open',output_name,'ierror=',ierror
 line1=' ATOM      ATOMIC                      COORDINATES (BOHR)'
 line2='     ATOMIC BASIS SET'
 openif: IF ( ierror == 0 ) THEN

    readloop: DO
        READ (10,'(a)',IOSTAT=ierror) output    !Get to next line


      IF ( output(1:50) .eq. line1(1:50) ) THEN   
          READ (10,'(a)',IOSTAT=ierror) output
          READ (10,9070) ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          write(*,*) "test",ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          !coord(i,:)=(/X(i),Y(i),Z(i)/)

          i=i+1
          READ (10,9070) ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          !coord(i,:)=(/X(i),Y(i),Z(i)/)


      end if
      IF ( output(1:50) .eq. line2(1:50) ) THEN   
          READ (10,'(a)',IOSTAT=ierror) output
          write(*,*) output
          READ (10,'(a)',IOSTAT=ierror) output
          write(*,*) output
          READ (10,'(a)',IOSTAT=ierror) output
          write(*,*) output
          READ (10,'(a)',IOSTAT=ierror) output
          write(*,*) output
          READ (10,'(a)',IOSTAT=ierror) output
          write(*,*) output
          READ (10,'(a)',IOSTAT=ierror) output
          write(*,*) output
          READ (10,'(a)',IOSTAT=ierror) output ! ATOM NAME
          write(*,*) output
          READ (10,'(a)',IOSTAT=ierror) output
          write(*,*) output
          do while (ierror.eq.0)  
            READ (10,9140,IOSTAT=ierror) SHELLNUM,OTYPE,PRIMNUM,EXPONENT,CONTRACTION
            write(*,*) 'ierror',ierror
            write(*,9140) SHELLNUM,OTYPE,PRIMNUM,EXPONENT,CONTRACTION
            write(*,*) SHELLNUM,OTYPE,PRIMNUM,EXPONENT,CONTRACTION
            if (SHELLNUM.gt.0) then
              write(*,*) 'All is ok '
              REALSHELLNUM=REALSHELLNUM+1
              write(*,*) 'REALSHELLNUM', REALSHELLNUM
              
            end if
          end do
      end if

 9140 FORMAT(4x,I3,3X,A2,4x,I4,F20.7,4F20.12)

!************************************************************************************************************
!   There is no more line to test. The IOSTAT is tested each line to determine when the file has been
!   entirely read
!************************************************************************************************************



        IF ( ierror /= 0 ) EXIT         !Exit if not valid
        !WRITE (*,'(93a)') output       !Echo to screen (93 = line lenght)


END DO readloop
 END IF openif
close(10)

9070 FORMAT(1X,A8,A2,3X,F5.1,3F20.10)
end subroutine scan_dimension_gamess
subroutine read_coeff_lcao_gamess(lcao,norb,basis_file)
!read coefficients
!****************************************************************************************
!****************************************************************************************


 character(len=6)::curseur
integer(kind=int_4),intent(in)::norb
integer(kind=int_4)::diviseur,reste,pointeur,i,j,k,inR,nbcol
real(kind=real_8)::coeff1,coeff2,coeff3,coeff4
real(kind=real_8),dimension(norb,norb),intent(out)::lcao
 character(len=32),intent(in)::basis_file
 character(len=64)::read_file,numeroR,output_name


write(*,*) "test"
end subroutine read_coeff_lcao_gamess
!****************************************************************************************
!****************************************************************************************


!****************************************************************************************
!****************************************************************************************
subroutine read_integrals_gamess(SO,T1,V1,mu,mux,muy,vee,norb,basis_file)
!read atomic integrals
!****************************************************************************************
!****************************************************************************************


Implicit None
! Variables du bloc de lecture
Integer(kind=int_4),intent(in) ::norb
Integer ::i,j,k,l,num,inR,f2
 character(len=7)::cable1,cable2,cable3,cable4,cable5
integer(kind=int_4)::table1,table2,table3,table4,table5,col2,col3,col4,col5
real(kind=real_8),dimension(:,:),intent(out)::SO,T1,V1,mu,mux,muy
real(kind=real_8),dimension(:,:,:,:),intent(out)::vee
real(kind=real_8)::col1
logical::Log_SO,Log_T1,Log_V1,Log_mu,Log_mux,Log_muy
 character(len=32),intent(in)::basis_file
 character(len=64)::read_file
 character(len=4)::numeroR
 

 INTEGER :: entier,entier2,i1,i2,IS1,JS1,KS1,LS1,IS2,JS2,KS2,LS2
 INTEGER ::nbasisfunctions,op,ip					!Number of basis functions
 Integer ::ierror,ierror2,ierror3					!I/O status
 INTEGER ::nlines,ibasisfunctions,iatom				!Number of lines read in
 INTEGER , allocatable ::indice(:)				!indice the represent an basis function
 REAL*8 :: temp1,temp2,NREC1,NREC2
 REAL*8 , allocatable :: energy(:),LCAO(:,:),integralOM(:,:),integralBNH(:,:),integralKE(:,:),integralZAM(:,:),integralX(:,:),integralY(:,:),integralZ(:,:),integral2E(:,:,:,:),LCAO2(:,:)		!energy
 CHARACTER(len=4),allocatable :: geo(:)				!
 character(len=100)::output,output2,output3			!specific line to read
 character(len=5)::atomsymb,iorb
 character(len=100)::output_name,output_name2,output_name3,system,H2,H2C2,H2O!Name of the file to open
 CHARACTER(len=100)::line1					!NUMBER OF CARTESIAN GAUSSIAN BASIS...
 CHARACTER(len=100)::line2,line3,line4,line5			!1 electron int..
 CHARACTER(len=100)::line6,line61				!eigenvectors
 CHARACTER(len=100)::line7,line8,line9				!electrostatic moment
 CHARACTER(len=100)::line10,line11				!2 electron integral
 CHARACTER(len=100)::line12,line13					!LCAO from .dat
 CHARACTER(len=100)::VEC,fEND

output_name='BASIS/'//trim(basis_file)//'/gamess1e.out'
output_name2='BASIS/'//trim(basis_file)//'/gamess2e.out'
output_name3='BASIS/'//trim(basis_file)//'/gamess.dat'
 !output_name='H2/H2_energy.out'	
 !output_name2='H2/acet_MCSCF_6311g_SP_NCORE_4_NACT_6_NELS_6_STSYM_AG_NSTATE_10.log'	
 !output_name3='H2/acet_MCSCF_6311g_SP_NCORE_4_NACT_6_NELS_6_STSYM_AG_NSTATE_10.dat'

line1=' TOTAL NUMBER OF BASIS SET SHELLS             ='		!Line to look for
line2=' OVERLAP MATRIX'						!
line3=' BARE NUCLEUS HAMILTONIAN INTEGRALS (H=T+V)'		!
line4=' KINETIC ENERGY INTEGRALS'				!
line5=' Z-ANGULAR MOMENTUM INTEGRALS'				!
line6='          EIGENVECTORS'					!
line61=' ...... END OF RHF CALCULATION ......'			!
line7='             X INTEGRALS'				!
line8='             Y INTEGRALS'				!
line9='             Z INTEGRALS'				!
line10='          2 ELECTRON INTEGRALS'				!
line11='  ...... END OF TWO-ELECTRON INTEGRALS .....'		!
line12='--- OPTIMIZED MCSCF MO-S ---'
line13='--- CLOSED SHELL ORBITALS ---'
VEC=' $VEC '
fEND=' $END '
entier=1							!to define the number of column
entier2=1












!************************************************************************************************************
!				~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				READLOOPS 1 ELECTRON INTEGRALS
!				~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!************************************************************************************************************


!************************************************************************************************************
!	Get the file name and echo it back to user
!************************************************************************************************************


 !WRITE (*,*) 'Please enter GAMESS output name: '
 !READ (*,*) output_name
 !WRITE (*,1001) output_name
 !1001 FORMAT (' ','The GAMESS output file name is: ', A)
!
!
!
!************************************************************************************************************
!	Open the file, and check for errors on open.
!************************************************************************************************************



 OPEN (UNIT=10,FILE=output_name,STATUS='old',form='formatted',IOSTAT=ierror)
 write(*,*) 'open',output_name,'ierror=',ierror
 openif: IF ( ierror == 0 ) THEN



!************************************************************************************************************
!	Open was ok. Read lines.
!************************************************************************************************************

	readloop: DO
		READ (10,'(a)',IOSTAT=ierror) output	!Get to next line



!************************************************************************************************************
!	The program check for line1
!************************************************************************************************************



		IF ( output(1:50) .eq. line1(1:50) ) THEN	!line1
			READ (10,1000) nbasisfunctions
			WRITE(*,*) 'nbasisfunctions',nbasisfunctions
				allocate(indice(nbasisfunctions))
				allocate(energy(nbasisfunctions))
				allocate(geo(nbasisfunctions))
				allocate(LCAO(nbasisfunctions,nbasisfunctions))
				allocate(integralOM(nbasisfunctions,nbasisfunctions))
				allocate(integralBNH(nbasisfunctions,nbasisfunctions))
				allocate(integralKE(nbasisfunctions,nbasisfunctions))
				allocate(integralZAM(nbasisfunctions,nbasisfunctions))
				allocate(integralX(nbasisfunctions,nbasisfunctions))
				allocate(integralY(nbasisfunctions,nbasisfunctions))
				allocate(integralZ(nbasisfunctions,nbasisfunctions))
				allocate(integral2E(nbasisfunctions,nbasisfunctions,nbasisfunctions,nbasisfunctions))
				allocate(LCAO2(nbasisfunctions,nbasisfunctions))
					LCAO=0.d0
					integralOM=0.d0
					integralBNH=0.d0
					integralKE=0.d0
					integralZAM=0.d0
					integralX=0.d0
					integralY=0.d0
					integralZ=0.d0
					integral2E=0.d0	
					LCAO2=0.d0		
		END IF





!************************************************************************************************************
!	The program check for line2
!************************************************************************************************************



		IF ( output(1:50) .eq. line2(1:50) ) THEN	!line2
		!WRITE(*,*) nbasisfunctions,output,int(nbasisfunctions/5),MOD(nbasisfunctions,5)
		DO i1 = 1,int(nbasisfunctions/5) + 1
			WRITE (*,*)"test i1=",i1
			READ (10,9008)
			READ (10,9028)indice(entier:entier+4)
			READ (10,9008)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:entier+4)
			!WRITE (*,9008)
			DO i2 = 1, nbasisfunctions-((i1-1)*5)	
			IF ( i2 .LT. 5 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralOM(entier:i2,i2)
				WRITE (*,*)"test i2=",i2
				WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralOM(entier:i2,i2)
			END IF
			IF ( i2 .GT. 4 ) THEN
				WRITE (*,*)"test i2=",i2
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralOM(entier:entier+4,i2)
				WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralOM(entier:entier+4,i2)
			END IF
			END DO
		END DO
		END IF

!WRITE(*,*) integralOM(:,:)


!************************************************************************************************************
!	The program check for line3
!************************************************************************************************************




		IF ( output(1:50) .eq. line3(1:50) ) THEN	!line3
		!WRITE(*,*) output
		DO i1 = 1,int(nbasisfunctions/5) + 1
			READ (10,9008)
			READ (10,9028)indice(entier:entier+4)
			READ (10,9008)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:entier+4)
			!WRITE (*,9008)
			DO i2 = 1, nbasisfunctions-((i1-1)*5)	
			IF ( i2 .LT. 5 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralBNH(entier:i2,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralBNH(entier:i2,i2)
			END IF
			IF ( i2 .GT. 4 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralBNH(entier:entier+4,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralBNH(entier:entier+4,i2)
			END IF
			END DO
		END DO
		END IF





!************************************************************************************************************
!	The program check for line4
!************************************************************************************************************




		IF ( output(1:50) .eq. line4(1:50) ) THEN	!line4
		!WRITE(*,*) output
		DO i1 = 1,int(nbasisfunctions/5) + 1
			READ (10,9008)
			READ (10,9028)indice(entier:entier+4)
			READ (10,9008)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:entier+4)
			!WRITE (*,9008)
			DO i2 = 1, nbasisfunctions-((i1-1)*5)	
			IF ( i2 .LT. 5 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralKE(entier:i2,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralKE(entier:i2,i2)
			END IF
			IF ( i2 .GT. 4 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralKE(entier:entier+4,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralKE(entier:entier+4,i2)
			END IF
			END DO
		END DO
		END IF




!************************************************************************************************************
!	The program check for line5
!************************************************************************************************************




		IF ( output(1:50) .eq. line5(1:50) ) THEN	!line5
		!WRITE(*,*) output
		DO i1 = 1,int(nbasisfunctions/5) + 1
			READ (10,9008)
			READ (10,9028)indice(entier:entier+4)
			READ (10,9008)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:entier+4)
			!WRITE (*,9008)
			DO i2 = 1, nbasisfunctions-((i1-1)*5)	
			IF ( i2 .LT. 5 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralZAM(entier:i2,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralZAM(entier:i2,i2)
			END IF
			IF ( i2 .GT. 4 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralZAM(entier:entier+4,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralZAM(entier:entier+4,i2)
			END IF
			END DO
		END DO
		END IF








!************************************************************************************************************
!	The program check for line6
!************************************************************************************************************



	a:IF ( output(1:50) .eq. line6(1:50) ) THEN	
                READ (10,'(a)',IOSTAT=ierror) output2



!************************************************************************************************************
!	When the line6 has been read, the LCAO coefficients are read in a number of column defined by 'nbasisfunctions'
!************************************************************************************************************



		b:IF ( nbasisfunctions .GT. 4 ) THEN
		DO i1 = 1, int(nbasisfunctions/5)
			READ (10,9008)
			READ (10,9028)indice(entier:entier+4)
			READ (10,9068)energy(entier:entier+4)
			READ (10,9078)geo(entier:entier+4)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:entier+4)
			!WRITE (*,9068)energy(entier:entier+4)
			!WRITE (*,9078)geo (entier:entier+4)
			DO i2 = 1, nbasisfunctions
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,LCAO(entier:entier+4,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,LCAO(entier:entier+4,i2)
			END DO
			entier=entier+5
		END DO
		END IF b
		c:IF ( MOD(nbasisfunctions,5) .GT. 0 ) THEN
			READ (10,9008)
			READ (10,9028)indice(entier:nbasisfunctions)
			READ (10,9068)energy(entier:nbasisfunctions)
			READ (10,9078)geo(entier:nbasisfunctions)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:nbasisfunctions)
			!WRITE (*,9068)energy(entier:nbasisfunctions)
			!WRITE (*,9078)geo (entier:nbasisfunctions)
			DO i2 = 1, nbasisfunctions
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,LCAO(entier:nbasisfunctions,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,LCAO(entier:nbasisfunctions,i2)
			END DO
		END IF c

	END IF a

entier=1





!************************************************************************************************************
!	The program check for line7
!************************************************************************************************************



		IF ( output(1:50) .eq. line7(1:50) ) THEN	!line7
		!WRITE(*,*) output
		READ (10,'(a)',IOSTAT=ierror) output2
		DO i1 = 1,int(nbasisfunctions/5) + 1
			READ (10,9008)
			READ (10,9028)indice(entier:entier+4)
			READ (10,9008)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:entier+4)
			!WRITE (*,9008)
			DO i2 = 1, nbasisfunctions-((i1-1)*5)	
			IF ( i2 .LT. 5 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralX(entier:i2,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralX(entier:i2,i2)
			END IF
			IF ( i2 .GT. 4 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralX(entier:entier+4,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralX(entier:entier+4,i2)
			END IF
			END DO
		END DO
		END IF





!************************************************************************************************************
!	The program check for line8
!************************************************************************************************************





		IF ( output(1:50) .eq. line8(1:50) ) THEN	!line8
		!WRITE(*,*) output
		READ (10,'(a)',IOSTAT=ierror) output2
		DO i1 = 1,int(nbasisfunctions/5) + 1
			READ (10,9008)
			READ (10,9028)indice(entier:entier+4)
			READ (10,9008)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:entier+4)
			!WRITE (*,9008)
			DO i2 = 1, nbasisfunctions-((i1-1)*5)	
			IF ( i2 .LT. 5 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralY(entier:i2,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralY(entier:i2,i2)
			END IF
			IF ( i2 .GT. 4 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralY(entier:entier+4,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralY(entier:entier+4,i2)
			END IF
			END DO
		END DO
		END IF





!************************************************************************************************************
!	The program check for line9
!************************************************************************************************************




		IF ( output(1:50) .eq. line9(1:50) ) THEN	!line9
		!WRITE(*,*) output
		READ (10,'(a)',IOSTAT=ierror) output2
		DO i1 = 1,int(nbasisfunctions/5) + 1
			READ (10,9008)
			READ (10,9028)indice(entier:entier+4)
			READ (10,9008)
			!WRITE (*,9008)
			!WRITE (*,9028)indice(entier:entier+4)
			!WRITE (*,9008)
			DO i2 = 1, nbasisfunctions-((i1-1)*5)	
			IF ( i2 .LT. 5 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralZ(entier:i2,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralZ(entier:i2,i2)
			END IF
			IF ( i2 .GT. 4 ) THEN
				READ (10,9048)ibasisfunctions,atomsymb,iatom,iorb,integralZ(entier:entier+4,i2)
				!WRITE (*,9048)ibasisfunctions,atomsymb,iatom,iorb,integralZ(entier:entier+4,i2)
			END IF
			END DO
		END DO
		END IF
		



!************************************************************************************************************
!	There is no more line to test. The IOSTAT is tested each line to determine when the file has been
!	entirely read
!************************************************************************************************************		



		IF ( ierror /= 0 ) EXIT			!Exit if not valid
		!WRITE (*,'(93a)') output		!Echo to screen (93 = line lenght)
	

END DO readloop

		



!************************************************************************************************************
!	The WHILE loop has terminated. Was it because of a READ error or
!	because of the end of the input file? 
!************************************************************************************************************



readif: IF ( ierror > 0 ) THEN 				!A READ error occured. Tell user.

	WRITE (*,1020) nlines + 1
	1020 FORMAT ('0','An error occured reading line ', I6)
	!ElSE						!The end of data was reached. Tell user.

	!WRITE (*,1030) nlines
	!1030 FORMAT ('0','End of file reached. There were ', I6, ' lines in the file.')
	END IF readif

	ELSE openif
	WRITE (*,1040) ierror
	1040 FORMAT (' ','Error opening file: IOSTAT = ', I6 )
	END IF openif




!************************************************************************************************************
!	Close file
!************************************************************************************************************


 CLOSE ( UNIT=10 )


!************************************************************************************************************
!	
!************************************************************************************************************























!************************************************************************************************************
!				~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				READLOOPS 2 ELECTRON INTEGRALS
!				~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!************************************************************************************************************
!
!************************************************************************************************************
!	Get the file name and echo it back to user
!************************************************************************************************************


 !WRITE (*,*) 'Please enter GAMESS output name: '
! READ (*,*) output_name
! WRITE (*,1001) output_name
! 1001 FORMAT (' ','The GAMESS output file name is: ', A)
!
!************************************************************************************************************
!	Open the file, and check for errors on open.
!************************************************************************************************************
 OPEN (UNIT=20,FILE=output_name2,STATUS='old',form='formatted',IOSTAT=ierror2)
IF ( ierror2 == 0 ) THEN



!************************************************************************************************************
!	Open was ok. Read lines.
!************************************************************************************************************



	DO

		READ (20,9088,IOSTAT=ierror2) IS1,JS1,KS1,LS1,NREC1,temp1,IS2,JS2,KS2,LS2,NREC2,temp2
		IF ( (ierror2 == 0) .AND. (IS1 .GT. 0) .AND. (NREC1 .GT. 0) .AND. (NREC2 .GT. 0)) THEN
                	integral2E(IS1,JS1,KS1,LS1)=temp1
                	integral2E(IS2,JS2,KS2,LS2)=temp2
			WRITE(*,9098) IS1,JS1,KS1,LS1,NREC1,integral2E(IS1,JS1,KS1,LS1)	
			WRITE(*,9098) IS2,JS2,KS2,LS2,NREC2,integral2E(IS2,JS2,KS2,LS2)	
		END IF
		IF ( (ierror2 == 0) .AND. (IS1 .GT. 0) .AND. (IS2 .EQ. 0)) THEN
                	integral2E(IS1,JS1,KS1,LS1)=temp1
			WRITE(*,9098) IS1,JS1,KS1,LS1,NREC1,integral2E(IS1,JS1,KS1,LS1)	
		END IF
                
		IF  (ierror2 .LT. 0) EXIT
		
	END DO 
else
  write(*,*) "Error when opening ",output_name2,ierror2
END IF






!************************************************************************************************************
!	Close file
!************************************************************************************************************

 CLOSE ( UNIT=20 )

!************************************************************************************************************
!	
!************************************************************************************************************




















!************************************************************************************************************
!			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!			READLOOPS LCAO COEFFICIENTS FROM .dat FILE
!			~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!************************************************************************************************************
!
!************************************************************************************************************
!	Get the file name and echo it back to user
!************************************************************************************************************


 !WRITE (*,*) 'Please enter GAMESS output name: '
! READ (*,*) output_name3
! WRITE (*,1001) output_name
! 1001 FORMAT (' ','The GAMESS output file name is: ', A)
!
!************************************************************************************************************
!	Open the file, and check for errors on open.
!************************************************************************************************************

!
!
! OPEN (UNIT=30,FILE=output_name3,STATUS='old',form='formatted',IOSTAT=ierror3)
!IF ( ierror3 == 0 ) THEN
!
!
!
!!************************************************************************************************************
!!	Open was ok. Read lines.
!!************************************************************************************************************
!
!DO
!	READ (30,'(a)',IOSTAT=ierror3) output3
!		IF ( output3(1:20) .eq. line12(1:20) ) THEN
!		DO 
!			READ (30,'(a)',IOSTAT=ierror3) output3
!			IF ( output3(1:20) .eq. VEC(1:20) ) THEN
!				DO i1=1, nbasisfunctions
!				IF ( nbasisfunctions .GT. 4 ) THEN
!					DO i2 = 1, int(nbasisfunctions/5)
!						READ (30,9108) op,ip,LCAO2(entier:entier+4,entier2)
!						WRITE (*,9108) op,ip,LCAO2(entier:entier+4,entier2)
!						entier=entier+5
!					END DO
!				END IF
!				IF ( MOD(nbasisfunctions,5) .GT. 0 ) THEN
!						READ (30,9108) op,ip,LCAO2(entier:entier+MOD(nbasisfunctions,5)-1,entier2)
!						WRITE (*,9108) op,ip,LCAO2(entier:entier+MOD(nbasisfunctions,5)-1,entier2)
!				END IF
!				entier2=entier2+1
!				END DO
!			END IF
!			!IF ( ierror3 /= 0 ) EXIT
!			IF ( output3(1:20) .eq. fEND(1:20) ) EXIT
!		END DO
!		END IF
!		IF ( ierror3 /= 0 ) EXIT
!END DO
!END IF
!
!
!
!
!
!
!
!!************************************************************************************************************
!!	Close file
!!************************************************************************************************************
!
! CLOSE ( UNIT=30 )

!************************************************************************************************************
!	
!************************************************************************************************************

write(*,*) "close unit 30"














!************************************************************************************************************
!	Formats
!************************************************************************************************************




			1000 format(1X,'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =',I5)
			9008 FORMAT(1X)
			9028 FORMAT(15X,5(4X,I4,3X))			!1,2,3,4,5,6,7
			9048 FORMAT(I5,2X,A1,2X,I1,2X,A1,1X,10F11.6)	!lcao
			9068 FORMAT(15X,10F11.4)			!energie
			9078 FORMAT(16X,10(5X,A4,2X))			!geo
			9088 FORMAT(2(4I4,F10.1,F12.9,1X))		!2 electron int
			9098 FORMAT(4I4,F10.1,F12.9,1X)			!2 electron int
			9108 FORMAT(I2,1X,I2,5E15.8)

write(*,*) "close unit 30"

end subroutine read_integrals_gamess
!

subroutine lecture_gaus_gamess(Norb,basis_file,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord)
!
!  lecture des informations concernant la base atomique
!
!	input:  
!		-Norb::		nombre d'orbitales atomiques
!
!	output:	
!		-ETA(:,:) 
!		-ZET(:,:)
!		-ICONU(:)
!		-NF:: (NS)
!		-LMMP1(:)
!
!****************************************************************************************
!****************************************************************************************

!I/O variables
integer(kind=int_4),intent(in)::Norb
 character(len=32),intent(in)::basis_file
real(kind=real_8), dimension(:,:),intent(out)::ZET,ETA,coord
integer(kind=int_4),dimension(:),intent(out)::ICONU
integer(kind=int_4),dimension(:),intent(out)::NF
integer(kind=int_4),dimension(:),intent(out)::LMNP1
integer(kind=int_4),intent(inout)::NCONS


!variables internes
integer(kind=int_4)::NGEN,NAORDS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT,i,j
integer(kind=int_4)::NDPT,IS,NS
integer(kind=int_4),dimension(NCONS)::SetCenter ! center of the contraction set
integer(kind=int_4)::counter1,ll
integer(kind=int_4), allocatable,dimension(:)::ICSU,NRCR,NC,MCONS,IGCS,MCRS,IGEN
real(kind=real_8), allocatable,dimension(:,:)::ZETtmp,ETAtmp
integer(kind=int_4)::nprimmax
real(kind=real_8),allocatable,dimension(:)::X,Y,Z,AN
 character(len=3),allocatable,dimension(:)::MTYPE
real,allocatable,dimension(:)::CHG
 character::NEANT
 character(len=21),parameter::FMT1 ='(t1,d16.8,t1,d16.8)'
 character(len=64)::read_file
 character(len=13)::pointeur
integer::nligne,reste
Integer(kind=int_4) ::ierror
 CHARACTER(len=8),dimension(:),allocatable :: ATOMNAME,ATOMNAME2

 CHARACTER(len=100)::line1,output
 CHARACTER(len=100)::linefd1,linefd2,linefd3

allocate(ATOMNAME(2),ATOMNAME2(2),X(2),Y(2),Z(2),AN(2))
allocate(MCONS(10),IGCS(10))
ZET=0.d0
ETA=0.d0
coord=0.d0
ICONU=0.d0
NF=0.d0
LMNP1=0.d0
i=1

linefd3=' ATOM      ATOMIC                      COORDINATES (BOHR)'
write(*,*)'Debut de lecture_gauss_gamess'
read_file='BASIS/'//trim(basis_file)//'/1/gamess1e.out'
read_file=trim(read_file)
open(47,file=read_file,status='old',form='formatted',IOSTAT=ierror)
open(10,file='read_gamess/H2/H2_energy.out',status='old',form='formatted',IOSTAT=ierror)
  openif: IF ( ierror == 0 ) THEN
    readloop: DO

      READ (10,'(a)',IOSTAT=ierror) output    !Get to next line
      IF ( output(1:50) .eq. linefd3(1:50) ) THEN   
          READ (10,'(a)',IOSTAT=ierror) output
          READ (10,9070) ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          write(*,*) "test",ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          coord(i,:)=(/X(i),Y(i),Z(i)/)

          i=i+1
          READ (10,9070) ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          coord(i,:)=(/X(i),Y(i),Z(i)/)


      end if
      IF ( ierror /= 0 ) EXIT
    end do readloop
  end if openif
  i=1
  MCONS(i)=33.8650000
  IGCS(i)=0.025493814541
  write(*,*) "test1111111"
  SetCenter(i)=1
  write(*,*) "test1111112"
  i=6
  MCONS(i)=33.8650000
  IGCS(i)=0.025493814541
  SetCenter(i)=2
  i=2
  MCONS(i)=5.0947900    
  IGCS(i)=0.190373108582
  SetCenter(i)=1
  i=7
  MCONS(i)=5.0947900    
  IGCS(i)=0.190373108582
  SetCenter(i)=2
  i=3
  MCONS(i)=1.1587900     
  IGCS(i)=0.852161486043
  SetCenter(i)=1
  i=8
  MCONS(i)=1.1587900     
  IGCS(i)=0.852161486043
  SetCenter(i)=2
  i=4
  MCONS(i)=0.3258400     
  IGCS(i)=1.000000000000
  SetCenter(i)=1
  i=9
  MCONS(i)=0.3258400     
  IGCS(i)=1.000000000000
  SetCenter(i)=2
  i=5
  MCONS(i)=0.1027410     
  IGCS(i)=1.000000000000
  SetCenter(i)=1
  i=10
  MCONS(i)=0.1027410     
  IGCS(i)=1.000000000000
  SetCenter(i)=2


close(47)

9070 FORMAT(1X,A8,A2,3X,F5.1,3F20.10)
end subroutine lecture_gaus_gamess
subroutine read_dimensions_gamess(norb,basis_file,NCONS,maxprim,NS)
! read key dimensions from columbus files
!****************************************************************************************
!****************************************************************************************
 character(len=32),intent(in)::basis_file
Integer(kind=int_4) ::norb,nbasisfunctions,ierror
Integer(kind=int_4),dimension(:),allocatable ::atommult
 character(len=64)::read_file
Integer(kind=int_4) :: NCONS,maxprim,NS,i
 CHARACTER(len=100)::line1,output
 CHARACTER(len=100)::output_name
 CHARACTER(len=100)::output_name2
 CHARACTER(len=100)::output_name3

 CHARACTER(len=100)::linefd1,linefd2,linefd3
real(kind=real_8),allocatable,dimension(:)::X,Y,Z,AN !Atomic number
 CHARACTER(len=8),dimension(:),allocatable :: ATOMNAME,ATOMNAME2
 CHARACTER(len=100)::test
 character(len=2) ::  shelltype
real(kind=real_8)::primexp,primcont1,primcont2,primcont3,primcont4
Integer(kind=int_4) ::  shellnumber, primnum
maxprim=6
line1=' TOTAL NUMBER OF BASIS SET SHELLS             ='     !Line to look for
linefd1=' NUMBER OF OCCUPIED ORBITALS (BETA )          ='   !wanted line is just after this one

linefd2='  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIENT(S)'
output_name='BASIS/'//trim(basis_file)//'/gamess1e.out'
output_name2='BASIS/'//trim(basis_file)//'/gamess2e.out'
output_name3='BASIS/'//trim(basis_file)//'/gamess.dat'
 !output_name='H2/H2_energy.out'	
open(10,file=output_name,status='old',form='formatted',IOSTAT=ierror)
  IF (ierror.gt.0) write(*,*) "error opening ",output_name,ierror
  write(*,*) "error opening ",output_name,ierror
  openif: IF ( ierror == 0 ) THEN


    readloop: DO
      READ (10,'(a)',IOSTAT=ierror) output    !Get to next line
      

      IF ( output(1:50) .eq. line1(1:50) ) THEN   
          READ (10,1000) nbasisfunctions
          WRITE (*,*) "nbasisfunctions=",nbasisfunctions
      end if
      IF ( output(1:50) .eq. linefd1(1:50) ) THEN   
          READ (10,1001) ns
          allocate(ATOMNAME(ns),ATOMNAME2(ns),X(ns),Y(ns),Z(ns),AN(ns),atommult(ns))
      end if

      IF ( ierror /= 0 ) EXIT
    end do readloop
  end if openif
  write(*,*) "TTTTTTTTEEEEEEESSSSSSSSSSTTTTTTTT"
norb=nbasisfunctions
close(10)
linefd3=' ATOM      ATOMIC                      COORDINATES (BOHR)'
i=1
atommult=0
NCONS=0
open(10,file='read_gamess/H2/H2_energy.out',status='old',form='formatted',IOSTAT=ierror)
  openif2: IF ( ierror == 0 ) THEN
    readloop2: DO
      READ (10,'(a)',IOSTAT=ierror) output    !Get to next line
      IF ( output(1:50) .eq. linefd3(1:50) ) THEN   
      !write(*,*) output
          READ (10,'(a)',IOSTAT=ierror) output
          READ (10,9070) ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          write(*,*) "test",ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          i=i+1
          READ (10,9070) ATOMNAME(i),test,AN(i),X(i),Y(i),Z(i)
          if (ATOMNAME(i).eq.ATOMNAME(i-1))then
            write(*,*) "SAME ATOME"
            atommult(i-1)=atommult(i-1)+1

          end if
          write(*,*) "test",ATOMNAME(i),AN(i),X(i),Y(i),Z(i)
      end if
      IF ( output(1:50) .eq. linefd2(1:50) ) THEN   
          READ (10,'(a)',IOSTAT=ierror) output
          READ (10,1002) ATOMNAME2(i)
          write(*,*) 'ATOMNAME2 = ', ATOMNAME2(i)
          READ (10,'(a)',IOSTAT=ierror) output
          !READ (10,9140) shellnumber, shelltype,primnum,primexp,primcont1!,primcont2,primcont3,primcont4
          !write(*,*) "WORKING TEST : ",shellnumber, shelltype,primnum,primexp,primcont1!,primcont2,primcont3,primcont4
          do 
            READ (10,9140) shellnumber, shelltype,primnum,primexp,primcont1!,primcont2,primcont3,primcont4
            30 if (shellnumber.gt.0) then
              NCONS=NCONS+1
              write(*,*) "WORKING TEST : ",shellnumber, shelltype,primnum,primexp,primcont1!,primcont2,primcont3,primcont4
            else

              write(*,*) "TESTTESTTEST1: ",shellnumber,shelltype,primnum,primexp,primcont1!,primcont2,primcont3,primcont4
              write(*,*) "TESTTESTTEST3: ",shellnumber,shelltype,primnum,primexp,primcont1!,primcont2,primcont3,primcont4

              READ (10,9140,IOSTAT=ierror) shellnumber, shelltype,primnum,primexp,primcont1!,primcont2,primcont3,primcont4
              if ((shellnumber.gt.0).and.ierror.eq.0) then
                NCONS=NCONS+1
                write(*,*) "WORKING TEST : ",shellnumber, shelltype,primnum,primexp,primcont1!,primcont2,primcont3,primcont4
              else
                goto 40
              end if
              write(*,*) "TESTTESTTEST2: ",output
              write(*,*) "shellnumber: ",shellnumber
              !READ (10,9140) shellnumber, shelltype,primnum,primexp,primcont1,primcont2,primcont3,primcont4
              !if (shellnumber.gt.0) go to 30
            end if
          end do
        40 continue
      end if
      IF ( ierror /= 0 ) EXIT
    end do readloop2
  end if openif2
  NCONS=NCONS*2
  write(*,*) "NCONS=",NCONS !TODO to be multiplied by the number of atom 
  !TODO do a loop for each atom type
  !TODO then sum all NCONS
close(10)

write(*,'(A37, I3)')'Number atomic orbitals to construct :',norb
1000 format(1X,'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =',I5)
1001 format(1X,'TOTAL NUMBER OF ATOMS                        =',I5)
1002 format(1X,A2)
1003 format(1X,A2)
9060 FORMAT(1X,A2,F5.1,F17.10,2F20.10) !line 6094 of fmolib.src
9070 FORMAT(1X,A8,A2,3X,F5.1,3F20.10)
9065 FORMAT(1X,A8,A2,F5.1)
9140 FORMAT(1X,I6,3X,A2,3x,I4,2x,F20.7,1x,F17.12)!,4F20.12)
end subroutine read_dimensions_gamess




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE read_gamess
