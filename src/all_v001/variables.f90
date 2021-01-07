module variables
!TODO tn must no be declared here
use basics

character(len=32):: basis_file,job
Integer(kind=int_4)::norb,Norb_cat,Norb_sym,nR,nPrim,Ne,tn,nt,dimP,nper,orb_Q,nk
Integer(kind=int_8)::dimCSF,dimQ
real(kind=real_8)::pdt,cycleopt,tmax,tfwhm,tau,tmin,omega,E0,dkOpt,delta,Spin
logical::restriction,pulsed
Integer(kind=int_4) ::totPrimCount


end module variables
