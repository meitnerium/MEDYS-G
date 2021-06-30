!Lecture des parametres d'input

Module read_input


use basics

contains
!****************************************************************************************
!****************************************************************************************
 subroutine param_debut(basis_file,nq,inR,Ne,Norb_cat,Norb_sym,FC,restriction,S,job,ni,imin,imax,di,dki,kmin,E0&
,omega,delta,theta,phy,pulsed,nper,pdt,tmin,nt,QCsoft,lp,ld,Iflg1,Iflg2)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
! Parameters that are read
 Character(len=32), Intent(out) 	              :: basis_file,job
 Integer(kind=int_4), Intent(out) 	            :: nq, inR
 Integer(kind=int_4), Intent(out) 	            :: Ne,Norb_cat,Norb_sym,FC !,opt_GS
 Integer(kind=int_4), Intent(out)	              :: nper,nt
 Integer(kind=int_4),dimension(3), Intent(out)	              :: lp
 Integer(kind=int_4),dimension(6), Intent(out)	              :: ld
 Integer(kind=int_4)			                      :: nx
 Integer(kind=int_4)		 	                      :: ny
 Integer(kind=int_4)		 	                      :: nz
 Integer(kind=int_4)		 	                      :: QCsoft,Iflg1,Iflg2 ! Iflg1 contrôle type de calcul (ex: restart ou nouveau)-Non disponible encore, Iflg2 contrôle le calcul du spectre photoelectronique (Iflg2=0 <-> no spectrum [à faire], Iflg2=1 <-> spectrum for all J, Iflg2=2 <-> spectrum for cetrain J [à faire])
 Integer(kind=int_4)		 	                      :: harmonic
 Real(kind=real_8)		 	                        :: omega,delta,theta,phy,pulsetime
 Real(kind=real_8)		 	                        :: dx,xmin,xmax ,dkx,kxmin
 Real(kind=real_8)		 	                        :: dy,ymin,ymax ,dky,kymin
 Real(kind=real_8)		 	                        :: dz,zmin,zmax ,dkz,kzmin
 Logical 				                                :: restriction,pulsed
 Real(kind=real_8), Intent(out) 	              :: S,E0,tmin
 Real(kind=real_8) 			                        :: tmax
 Real(kind=real_8)			                        :: pdt,version
!TODO : verifier ces trois parametres
 Integer(kind=int_4), Dimension(3), Intent(out)	:: ni
 Real(kind=real_8), Dimension(3), Intent(out)	  :: imin,imax,di  !,dki,kmin
 Real(kind=real_8), Dimension(3), Intent(out)	  :: dki,kmin

Namelist/grille/basis_file,nq,inR,job,QCsoft
 Namelist/planewaves/nx,xmin,ny,ymin,nz,zmin 
Namelist/drt/Ne,Norb_cat,Norb_sym,FC,S,restriction
Namelist/champ/E0,omega,delta,theta,phy,pulsed,nper,pdt,tmin,harmonic
Namelist/prim_spec/lp,ld
Namelist/control_flags/Iflg1,Iflg2

version=0.90

open(1,file='input',status='old')
read(1,nml=grille)
 read(1,nml=planewaves)
read(1,nml=drt)
read(1,nml=champ)
read(1,nml=prim_spec)
read(1,nml=control_flags)
 close(1)
basis_file=trim(basis_file)
omega=omega/dfloat(harmonic)
ni(1) = nx
if (nx .eq. 1)then
        imin(1) = 0.d0
        imax(1) = 0.d0
        di(1) = 0.d0
        dki(1) = 0.d0
        kmin(1) = 0.d0
	xmin=0.d0
	xmax=0.d0
	dx=0.d0
	dkx=0.d0
	kxmin=0.d0
else if (nx .lt. 1)then
	STOP 'nx is wrong in the input'
else
        imin(1) = xmin
        imax(1) = -xmin
        di(1) = 2.d0*dabs(xmin)/dfloat(nx-1)
	xmax=-xmin
	dx=2.d0*dabs(xmin)/dfloat(nx-1)
        dki(1) = (2.d0*pi)/(dfloat(nx)*dx)
        kmin(1) = -pi/dx
	dkx=(2.d0*pi)/(dfloat(nx)*dx)
	kxmin=-pi/dx

endif

ni(2) = ny
if (ny .eq. 1)then
        imin(2) = 0.d0
        imax(2) = 0.d0
        di(2) = 0.d0
        dki(2) = 0.d0
        kmin(2) = 0.d0
	ymax=0.d0
	ymin=0.d0
	dy=0.d0
	dky=0.d0
	kymin=0.d0
else if (ny .lt. 1)then
	STOP 'ny is wrong, check the input'
else
        imin(2) = ymin
        imax(2) = -ymin
        di(2) = 2.d0*dabs(ymin)/dfloat(ny-1)
	ymax=-ymin
	dy=2.d0*dabs(ymin)/dfloat(ny-1)
        dki(2) = (2.d0*pi)/(dfloat(ny)*dy)
        kmin(2) = -pi/dy
	dky=(2.d0*pi)/(dfloat(ny)*dy)
	kymin=-pi/dy
endif

ni(3) = nz
if (nz .eq. 1)then
        imin(3) = 0.d0
        imax(3) = 0.d0
        di(3) = 0.d0
        dki(3) = 0.d0
        kmin(3) = 0.d0
	zmax=0.d0
	zmin=0.d0
	dz=0.d0
	dkz=0.d0
	dki(3)=0.d0
	kzmin=0.d0
else if (nz .lt. 1)then
	STOP 'nz is wrong, check the input'
else
        imin(3) = zmin
        imax(3) = -zmin
        di(3) = 2.d0*dabs(zmin)/dfloat(nz-1)
	zmax=-zmin
	dz=2.d0*dabs(zmin)/dfloat(nz-1)
        dki(3) = (2.d0*pi)/(dfloat(nz)*dz)
        kmin(3) = -pi/dz
	dkz=(2.d0*pi)/(dfloat(nz)*dz)
	kzmin=-pi/dz

endif

nt=(nper*pdt)

theta=theta*pi
phy=phy*pi




pulsetime=(dfloat(nper)*2.d0*pi)/omega
pdt=pulsetime/(dfloat(nt-1))
write(*,*) "test on line 132 of readinput in all_v001"
write(*,*)' basis_file ->',basis_file
write(*,*)' job ->',job

open(50,file='Parameters.dat',status='unknown',form='formatted')

Write(50,'(A20)')'PARAMETERS SUBMITTED'
Write(50,*)
Write(50,*) 'Medys Version : ', version
Write(50,*)

Write(50,'(A15)')'PHYSICAL SYSTEM'
Write(50,'(I2,A27,F4.1)')Ne,' electrons with a spin of ',S
Write(50,*)

Write(50,'(A15)')'LASER'
Write(50,'(A11,ES9.2,A8)')'Intensity: ',E0,' a.u.'
Write(50,'(A11,ES9.2,A8)')'Intensity: ',E0**2*3.51d16,' W/cm^2 '
Write(50,'(A12,F9.2,A8,I3,A12)')'Wavelength: ',omega,' a.u. and ',nper,' optic cycles'
Write(50,'(A12,F9.2,A8,I3,A12)')'Wavelength: ',omega*219474.6305d0,' cm-1'
Write(50,*)'Wavelength: ',1.d7/(omega*219474.6305d0),' nm'
!Write(50,'(A12,F9.2,A8,I3,A12)')'Wavelength: ',1.d7/(omega*219474.6305d0),' nm'
Write(50,'(A7,F4.1,A3)')'Phase: ',delta, ' pi' 
Write(50,'(A39,F4.1,A3)')'Angle between polarisation and Z axis: ', theta/pi,' pi'
Write(50,'(A39,F4.1,A3)')'Angle between polarisation and X axis: ', phy/pi,' pi'
Write(50,*)


write(50,'(A19)') 'CONSIDERED ORBITALS'
write(50,'(A21,I3,A6,I3,A13)')'Q sub-space orbitals ',Norb_cat+FC,' with ',FC,' frozen cores'
if (restriction .eq. .TRUE.) then
	write(50,'(A23)')'With a symbolic orbital'
else
	write(50,'(A28)')'Without a symbolic orbital'
end if
Write(50,*)

write(50,'(A17)') 'SPACIAL GRID SIZE'
write(50,'(A6,F9.2,A3,F9.2,A4,E10.3,A4,I5)')'axe X ',xmin,' a ',xmax,' dx=',dx,' nx=',nx
write(50,'(A6,F9.2,A3,F9.2,A4,E10.3,A4,I5)')'axe Y ',ymin,' a ',ymax,' dy=',dy,' ny=',ny
write(50,'(A6,F9.2,A3,F9.2,A4,E10.3,A4,I5)')'axe Z ',zmin,' a ',zmax,' dz=',dz,' nz=',nz
Write(50,*)

write(50,'(A18)') 'MOMENTUM GRID SIZE'
write(50,'(A7,F9.2,A3,F9.2,A6,E10.3,A5,I4)')'axe KX ',kxmin,' a ',-kxmin,' dkx=',dkx,' nkx=',nx
write(50,'(A7,F9.2,A3,F9.2,A6,E10.3,A5,I4)')'axe KY ',kymin,' a ',-kymin,' dky=',dky,' nky=',ny
write(50,'(A7,F9.2,A3,F9.2,A6,E10.3,A5,I4)')'axe KZ ',kzmin,' a ',-kzmin,' dkz=',dkz,' nkz=',nz
Write(50,*)

write(50,'(A9)') 'TIME GRID'
write(50,'(A20,I5)')'Number of time steps',nt
Write(50,*)'Timestep',pdt,' a.u.'
!Write(50,'(A8,F9.2,A6)')'Timestep',pdt,' a.u.'
!Write(50,'(A13,F9.2,A20,F9.2,A5)')'Initial time',tmin,' a.u. and Final time',tmax,' a.u.' 

 close(50)


end subroutine param_debut

End Module read_input

