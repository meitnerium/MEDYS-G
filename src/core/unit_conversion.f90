Module unit_conversion

Use basics

Implicit None

!Private
!Public::expomanti,real 

type frequency
real(kind=real_8)::value
character(len=10)::unit_
contains
procedure:: frequency_to_real
procedure:: frequency_to_au
end type

Real(kind=real_8),parameter::AUtoSIcharge= 1.6021765314E-19! C
Real(kind=real_8),parameter::SItoAUcharge= 6.241509474153817E+18
Real(kind=real_8),parameter::AUtoSImass=9.109382616E-31! Kg 
Real(kind=real_8),parameter::SItoAUmass=1.097769236571059E+30
Real(kind=real_8),parameter::AUtoSIaction=1.0545716818E-34! J s
Real(kind=real_8),parameter::SItoAUaction=9.482522783971838E+33
Real(kind=real_8),parameter::AUtoSIlength=0.529177210818E-10! m
Real(kind=real_8),parameter::SItoAUlength=1.889726124929311E+10
Real(kind=real_8),parameter::AUtoSIenergy=4.3597441775E-18! J
Real(kind=real_8),parameter::SItoAUenergy=2.293712564973085E+17
Real(kind=real_8),parameter::AUtoSItime=2.41888432650516E-17! s
Real(kind=real_8),parameter::SItoAUtime=4.134137333655862E+16


Contains

Real(kind=real_8) function frequency_to_real (freq)
class(frequency),intent(in)::freq
frequency_to_real=freq%value
end function

type(frequency) function frequency_to_au(freq)
class(frequency),intent(in)::freq
frequency_to_au%unit_=freq%unit_
Select case(trim(freq%unit_))
 case("Hz")
 frequency_to_au%value=freq%value*SItoAUtime
 case("au")
 frequency_to_au%value=freq%value
 case default
  Stop "frequency_to_au: entered unit is not known"
end select
End function
End module unit_conversion
