subroutine Get_GMut1DG(zet1,zet2,R1,R2,l1,l2,v1,v2,dlt,alpha,Anl,res)
!***************************************************
!***************************************************
implicit none
Real(8),Intent(in)    :: alpha,Anl,zet1,zet2,R1,R2
Integer(4) ::l1,l2,v1,v2,lv1,lv2,llvv,l12,dlt
Complex(8) :: res,rho 
double precision :: pdt,pi
integer(4) :: fact
Complex(8)  :: temp,temp2,prefact,capB,zetC,Im_u
Real(8)    :: z1Mod,z1Phase,z2Mod,z2Phase,zCMod,zCPhase,a
integer(4)              :: k,kp,ksum, k1,k2,k3,k4,p,pp,pmax,ppmax

OPEN(UNIT=88, FILE='output.data',Access = 'append')
write(88,*) 'zet1 ',zet1
write(88,*) 'zet2 ',zet2
write(88,*) 'R1 ', R1
write(88,*) 'R2 ', R2
write(88,*) 'l1 ', l1
write(88,*) 'l2 ', l2
write(88,*) 'v1 ', v1
write(88,*) 'v2 ', v2
write(88,*) 'dlt ', dlt
write(88,*) 'alpha ', alpha
write(88,*) 'Anl ', Anl

Im_u=1.d0     !Def. found in Module basics : Real(kind=real_8),parameter :: pio2=pi/2.0_real_8
pi=3.141592654d0
 

prefact=sqrt(sqrt(4*zet1*zet2/pi**2 ))*cdexp(dcmplx(0.d0,-(R2+alpha)*Anl))
prefact=sqrt(pi)*sqrt( ((4*zet1)**l1)*((4*zet2)**l2) )*prefact
a=0.d0
zetC=dcmplx(zet1+zet2,0.d0)
rho=dcmplx(1.d0,0.0)
pmax=0
ppmax=0
a=1/(2.d0*dlt*pdt)
zetC=dcmplx(zet2,a)
z2Mod=cdabs(zetC)
z2Phase=1.d0
rho=Im_u*a*zetC/z2Mod**2
zetC=  -1.d0
pmax=l2
prefact=1.d0

do k=0,pmax
    do kp=0, ppmax
        p=k+kp
        if((dlt.ne.0).and.(Mod(p,2).eq.0)) then
            temp=float(comb(l2,k))*float(comb(v2,kp))*fact(p-1,2)/(sqrt(2.d0)**p) &
            *(sqrt(z2Mod)**(-(p+1))*cdexp(dcmplx(0.d0,-0.5d0*z2Phase*(p+1))))
        endif
!        temp2=0.d0
!        do k1= 0,l1
!            do k2=0,v1
!                do k3=0,pmax-k
!                    do k4=0,ppmax-kp
!                        ksum=k1+k2+k3+k4
!                        if(Mod(ksum,2).eq.0) then
!                            temp2=temp2+comb(l1,k1)*comb(v1,k2)*comb(pmax-k,k3)*comb(ppmax-kp,k4)*(rho**(pmax-k+k4))*((capB+R2-R1+alpha)**
!(l1-k1))*((capB+R2+alpha)**(v1-k2))*(capB**(pmax-k-k3))*((R2-rho*capB)**(ppmax-kp-k4))*(fact(ksum-1,2)/(dsqrt(2.d0)**ksum) )*(dsqrt(zCMod)
!)**(-(ksum+1))*cdexp(dcmplx(0.d0,-0.5d0*zCPhase*(ksum+1)))
!                           write(*,*) "temp2 = ", temp2, "res = ", res
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        res=res+temp*temp2
    enddo
enddo


zCMod=abs(zetC)
zCPhase=1.d0
capB=2.d0 
write(*,*) "capB",capB
prefact=1.d0
prefact=1.d0
write(*,*) "prefact",prefact


temp=1.d0
res=0.d0


res = prefact*res
write(88,*) 'RES = '
write(88,*) res
write(88,*) '**********************************************'
close(UNIT=88)

end subroutine Get_GMut1DG
!***************************************************
!***************************************************
integer function comb(n,k)
!***************************************************
!***************************************************
implicit none
integer(4), intent (in) :: n
integer(4), intent (in) :: k
 
comb = 2

end function





!
!!***************************************************
!!***************************************************
!!***************************************************
!subroutine Get_GMut1DG(zet1,zet2,R1,R2,l1,l2,v1,v2,dlt,alpha,Anl,res)
!!***************************************************
!!***************************************************
!Real(kind=real_8),Intent(in)                    :: alpha,Anl,zet1,zet2,R1,R2
!Integer(kind=int_4) ::l1,l2,v1,v2,lv1,lv2,llvv,l12,dlt
!Complex(kind=comp_16)   :: res,rho
!
!Complex(kind=comp_16)   :: temp,temp2,prefact,capB,zetC,Im_u
!Real(kind=real_8)               :: z1Mod,z1Phase,z2Mod,z2Phase,zCMod,zCPhase,a
!integer                               :: k,kp,ksum, k1,k2,k3,k4,p,pp,pmax,ppmax
!
!OPEN(UNIT=88, FILE='output.data',Access = 'append')
!write(88,*) 'zet1 ',zet1
!write(88,*) 'zet2 ',zet2
!write(88,*) 'R1 ', R1
!write(88,*) 'R2 ', R2
!write(88,*) 'l1 ', l1
!write(88,*) 'l2 ', l2
!write(88,*) 'v1 ', v1
!write(88,*) 'v2 ', v2
!write(88,*) 'dlt ', dlt
!write(88,*) 'alpha ', alpha
!write(88,*) 'Anl ', Anl

!Im_u=dcmplx(0.d0, 1.d0) !Def. found in Module basics : Real(kind=real_8),parameter :: pio2=pi/2.0_real_8

!
!
!prefact=sqrt(sqrt(4*zet1*zet2/pi**2 ))*cdexp(dcmplx(0.d0,-(R2+alpha)*Anl))
!prefact=sqrt(pi)*sqrt( ((4*zet1)**l1)*((4*zet2)**l2) )*prefact
!select case(dlt)
!case(0)
!a=0.d0
!zetC=dcmplx(zet1+zet2,0.d0)
!rho=dcmplx(1.d0,0.0)
!pmax=0
!ppmax=0
!case default
!a=1/(2.d0*dlt*pdt)
!zetC=dcmplx(zet2,a)
!z2Mod=abs(zetC)
!z2Phase=atan2(dimag(zetC),dreal(zetC))
!rho=Im_u*a*zetC/z2Mod**2
!zetC=  -Im_u*a*(1.d0+rho)+ zet1 ! - Im_u*a+ (a**2)zetC/z2Mod**2+zet1
!pmax=l2
!ppmax=v2
!prefact=dsqrt(a)*cdexp(dcmplx(0.d0, -pi/4))*prefact
!end select


!zCMod=abs(zetC)
!zCPhase=atan2(dimag(zetC),dreal(zetC)) !todo : verifier atan ou atan2
!capB=-dcmplx(zet1*(R2-R1 +alpha) , -Anl/2.d0)*cdexp(dcmplx(0.d0,-zCPhase))/zCMod !todo: a verifier
!write(*,*) "capB",capB
!prefact=cdexp(zetC*capB*capB)*prefact
!prefact=cdexp(-zet1*(dcmplx((alpha + R2-R1)*(alpha + R2-R1), 0.d0)))*prefact
!write(*,*) "prefact",prefact
!
!
!temp=dcmplx(1.d0,0.d0)
!res=dcmplx(0.d0,0.0)


!do k=0,pmax
!    do kp=0, ppmax
!        p=k+kp
!        if((dlt.ne.0).and.(Mod(p,2).eq.0)) then
!            temp=comb(l2,k)*comb(v2,kp)*fact( p-1,2)/(sqrt(2.d0)**p)*(sqrt(z2Mod)**(-(p+1))*cdexp(dcmplx(0.d0,-0.5d0*z2Phase*(p+1))))
!        endif
!        temp2=0.d0
!        do k1= 0,l1
!            do k2=0,v1
!                do k3=0,pmax-k
!                    do k4=0,ppmax-kp
!                        ksum=k1+k2+k3+k4
!                        if(Mod(ksum,2).eq.0) then
!                            temp2=temp2+comb(l1,k1)*comb(v1,k2)*comb(pmax-k,k3)*comb(ppmax-kp,k4)*(rho**(pmax-k+k4))*((capB+R2-R1+alpha)**
!(l1-k1))*((capB+R2+alpha)**(v1-k2))*(capB**(pmax-k-k3))*((R2-rho*capB)**(ppmax-kp-k4))*(fact(ksum-1,2)/(dsqrt(2.d0)**ksum) )*(dsqrt(zCMod)
!)**(-(ksum+1))*cdexp(dcmplx(0.d0,-0.5d0*zCPhase*(ksum+1)))
!                           write(*,*) "temp2 = ", temp2, "res = ", res
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        res=res+temp*temp2
!    enddo
!enddo
!res = prefact*res
!write(88,*) 'RES = '
!write(88,*) res
!write(88,*) '**********************************************'
!close(UNIT=88)

!end subroutine Get_GMut1DG

