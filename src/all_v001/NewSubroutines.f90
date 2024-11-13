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

subroutine vee_moTOmo(vee,moRot,norb,norb1,vee2t)
    !  
    !****************************************************************************************
    !****************************************************************************************


    use basics
    use math !chop
    !use blas95
    !use f95_precision

    implicit none
    integer(kind=int_4), intent(in) :: norb,norb1
    complex(kind = comp_16), intent(inout), dimension(:,:,:,:) :: Vee
    complex(kind = comp_16),  dimension(norb,norb,norb,norb) :: vee2
    complex(kind = comp_16),  dimension(norb,norb,norb,norb) :: vee2t
    complex(kind = comp_16), intent(in), dimension(norb,norb1) :: moRot
    complex(kind = comp_16), dimension(norb1,norb) :: moRotCon, Y !!!!! ,temp
    real(kind=real_8) :: tmp,tempo
    integer(kind=int_4) :: i,j,k,l,r,s,t,u
    complex(kind = comp_16),  dimension(norb1,norb1) :: X,X_temp !!!!!,Y
    complex(kind = comp_16),   dimension(norb,norb,norb,norb) :: temp2e

    ! Variables used for System Clock
    real(kind=real_8) :: tdeb,tfin

    ! Pour la parallelisation
    integer NTHREADS, TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
    LOGICAL :: file_exists
 
 moRotCon(:,:)=transpose(conjg(moRot(:,:)))
 do l=1,norb
            do k=1,norb
                do j=1,norb
                    do i=1,norb
                        X_temp(i,j)=Vee(i,j,k,l)
                    end do
                end do
                call gemm(moRotCon,X_temp,Y)       
                call gemm( Y,moRot, X)
                do j=1,norb1
                    do i=1,norb1
                        temp2e(i,j,k,l)=X(i,j)
                    end do
                end do
            end do
        end do
        ! TODO : add parrallel do at begining
        X=0.d0
        Y=0.d0

        !$OMP PARALLEL DO SHARED(temp2e,Vee2,Vee) PRIVATE(NTHREADS,Y,X,TID,i,j,k,l)
        do l=1,norb1
            do k=1,norb1

                do j=1,norb
                    do i=1,norb
                        X_temp(i,j)=temp2e(k,l,i,j)
                    end do
                end do

                call gemm(moRotCon,X_temp,Y)       
                call gemm( Y,moRot, X)
                do j=1,norb1
                    do i=1,norb1
                        vee2t(i,j,k,l)=X(i,j)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO


!        write(*,*)'Done.'
!        write(*,*)
        
        end subroutine


End module
 
