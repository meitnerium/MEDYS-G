module champ1
!crée le champ
use Basics
use math
use variables
contains




!****************************************************************************************
!****************************************************************************************
function champNint( )
!****************************************************************************************
!****************************************************************************************

real(kind=real_8)::champNint
real(kind=real_8)::t!,delta

!real(kind=real_8)::cycleopt,tmax,E0,w,Tfwhm,tau,t,champ,champcontinu,enveloppe
!integer::nn
real(kind=real_8),allocatable,dimension(:)::valeur_champ

t=(tn*pdt)+tmin

! champNint=E0*Sin(w*t)*(sin((Pi*t)/(2.d0*3.d0*cycleopt))**2.d0)
 champNint=chantMich(0)
end function champNint



function Int0( )
!****************************************************************************************
!****************************************************************************************

real(kind=real_8)::Int0
real(kind=real_8)::t!,delta
!real(kind=real_8)::cycleopt,tmax,E0,w,Tfwhm,tau,t,champ,champcontinu,enveloppe
!integer::nn
real(kind=real_8),allocatable,dimension(:)::valeur_champ


t=(dfloat(tn)*pdt)+tmin


if (tn.eq.tmax) then
	Int0=0.d0
else

	Int0=champNint()*pdt
!	Int0=chantMich(1)
end if
end function Int0










!!!!!!Subroutine Michel




!****************************************************************************************
!****************************************************************************************
function chantMich(ordre)
!****************************************************************************************
!****************************************************************************************
Real(kind=real_8) ::chantMich
Real(kind=real_8) ::t
real(kind=real_8) :: nperf!,delta
integer(kind=int_4) :: nperi
integer(kind=int_4),intent(in)::ordre



!pulsed=.false.
t=tmin+ (dfloat(tn)*pdt)
!nperf=(tmax/2.d0)*w/(2.d0*pi)
!nperi=idint((tmax/2.d0)*w/(2.d0*pi))
!nper=Max((2*nperi)+nint(nperf-dfloat(nperi)),2)
!TODO change delta using input value

Select case(ordre)
	case(0)
		chantMich=EE(E0,omega,t,pulsed,nper)
	case(1)
		chantMich=AAA(E0,omega,t+pdt,pulsed,nper)-AAA(E0,omega,t,pulsed,nper)
end select

end function chantMich

!****************************************************************************************
!****************************************************************************************
Function Env(omega,t,pulsed,nper)
!****************************************************************************************
!****************************************************************************************

use basics
implicit none
real(kind=real_8) :: Env
Real(kind=real_8), intent(in) :: omega,t
integer(kind=int_4), intent(inout) :: nper
real(kind=real_8) ::tmax
logical,intent(in) :: pulsed


! tmax here denotes the time associated to maximum field amplitude
tmax=dfloat(nper)*Pi/omega
if(pulsed.eq..true.)then
	if(t.lt.2.d0*tmax)then

	Env=dsin(Pi*t/tmax/2.d0)**2.d0
	else
	Env=0.d0
	endif

else
Env=1.d0
endif


end function

!****************************************************************************************
!****************************************************************************************
Function EE(E0,omega,t,  pulsed,nper)
!****************************************************************************************
!****************************************************************************************

Use Basics

Implicit None
Real(kind=real_8) :: EE
Real(kind=real_8), intent(in) :: E0,t,omega!,delta
integer(kind=int_4), intent(inout) :: nper
logical,intent(in) :: pulsed
!TODO change this to use input value
if(pulsed.eq..true.)then
EE=Env(omega,t,pulsed,nper)*E0*sin(omega*t+delta)
else
if(t.le.dfloat(nper)*2.d0*Pi/omega)then
EE=Env(omega,t,pulsed,nper)*E0*sin(omega*t+delta)
else
EE=0.d0
endif
endif

End Function

!****************************************************************************************
!****************************************************************************************
Function AAA(E0,omega,tn, pulsed,nper)
!****************************************************************************************
!****************************************************************************************
! Attention this function is only valid for an integer number of optical cycles
! This analytical integration also only stands for more than one cycle pulses,
! (that is because of the one over n-1 factor in the expressions)
! a generalization could be done later if one is interessed in such pulses
Use Basics

Implicit None
Real(kind=real_8) :: AAA
Real(kind=real_8), intent(in) :: E0,omega,tn !delta
integer(kind=int_4), intent(in) :: nper
logical,intent(in) :: pulsed
Integer(kind=int_4) :: n
Real(kind=real_8) :: tf

n=nper
tf = dfloat(2*nper)*Pi/omega

if(pulsed.eq..true.)then

AAA = E0*(-1.d0*dcos(delta)*dcos(omega*tn)/(2.d0*omega) &
& + dsin(delta)*dsin(omega*tn)/(2.d0*omega) &
& + (dfloat(n)*dcos(delta)/(4.d0*omega))*(dcos(dfloat(n+1)/dfloat(n)*omega*tn)/dfloat(n+1)&
& + dcos(dfloat(n-1)/dfloat(n)*omega*tn)/dfloat(n-1))&
& - (dfloat(n)*dsin(delta)/(4.d0*omega))*(dsin(dfloat(n+1)/dfloat(n)*omega*tn)/dfloat(n+1)&
& + dsin(dfloat(n-1)/dfloat(n)*omega*tn)/dfloat(n-1)))

else

AAA = -1.d0*E0*dcos(tn*omega)/omega

endif

!if(tn.gt.tf)then
!AAA = 0.d0
!endif

End Function


end module champ1
