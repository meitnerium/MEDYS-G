program test_4Tung
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
use dynamique_new

Real(kind=real_8)                    :: alpha,Anl,zet1,zet2,R1,R2
Integer(kind=int_4) ::l1,l2,v1,v2,lv1,lv2,llvv,l12,dlt
Complex(kind=comp_16)   :: res,rho 

Complex(kind=comp_16)   :: temp,temp2,prefact,capB,zetC,Im_u
Real(kind=real_8)               :: z1Mod,z1Phase,z2Mod,z2Phase,zCMod,zCPhase,a
integer                               :: k,kp,ksum, k1,k2,k3,k4,p,pp,pmax,ppmax




Namelist/test4Tung/zet1,zet2,R1,R2,l1,l2,v1,v2,dlt,alpha,Anl

open(1,file='input_test4tung',status='old')
read(1,nml=test4Tung)
close(1)


call Get_GMut1DG(zet1,zet2,R1,R2,l1,l2,v1,v2,dlt,alpha,Anl,res)
write(*,*) "Res = ", res

end program test_4Tung

