Module NewSubs

use Basics

contains

!****************************************************************************************
!****************************************************************************************
subroutine CSF2State(FctP,eigenVectH0pp,StateP)
!****************************************************************************************
!****************************************************************************************
Complex(kind=comp_16),dimension(:,:),intent(in)		::fctP
Complex(kind=comp_16),Dimension(:,:),intent(in)		::eigenVectH0pp

Complex(kind=comp_16),dimension(:,:),intent(out)	::StateP

Integer(kind=int_4)					::i,j,sizeP

sizeP=size(fctP(:,1))

StateP = dcmplx(0.d0,0.d0)

do i=1,sizeP
do j=1,sizeP
	StateP(j,:)=eigenVectH0pp(i,j)*fctP(i,:)
end do
end do

end subroutine

End module
 
