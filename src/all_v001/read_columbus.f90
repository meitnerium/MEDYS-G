Module read_columbus

    ! read informations from columbus files
    use basics
    use blas95
    use lapack95
    use f95_precision
use variables
contains

    !****************************************************************************************
    !****************************************************************************************
    subroutine read_dimensions(norb,basis_file)
        ! read key dimensions from columbus files
        !****************************************************************************************
        !****************************************************************************************
        Integer(kind=int_4) ::norb
        character(len=32),intent(in)::basis_file
        character(len=64)::read_file

        read_file='BASIS/'//trim(basis_file)//'/1/integral_1e_2e.dat'
        read_file=trim(read_file)

        open(40,file=read_file,status='old',form='formatted')

        read(40,'(t21,i7)')norb
        write(*,'(A37, I3)')'Number atomic orbitals to be constructed :',norb

    end subroutine read_dimensions


    !****************************************************************************************
    !****************************************************************************************
    subroutine read_integrals(SO,T1,V1,mu,mux,muy,vee,norb,basis_file)
        !read atomic integrals
        !****************************************************************************************
        !****************************************************************************************


        Implicit None
        ! Variables du bloc de lecture
        Integer(kind=int_4),intent(in) ::norb
        Integer ::i,j,k,l,num,inR,f2,ii
        character(len=7)::cable1,cable2,cable3,cable4,cable5
        integer(kind=int_4)::table1,table2,table3,table4,table5,col2,col3,col4,col5
        real(kind=real_8),dimension(:,:),intent(out)::SO,T1,V1,mu,mux,muy
        real(kind=real_8),dimension(:,:,:,:),intent(out)::vee
        real(kind=real_8)::col1
        logical::Log_SO,Log_T1,Log_V1,Log_mu,Log_mux,Log_muy
        character(len=32),intent(in)::basis_file
        character(len=64)::read_file
        character(len=4)::numeroR


        !
        !BLOC DE LECTURE DU FICHIER COLOMBUS
        !On lit SO,V1,T1,mu(z),muy,mux et Vee
        !

        SO=0.d0
        v1=0.d0
        t1=0.d0
        mu=0.d0
        mux=0.d0
        muy=0.d0
        vee=0.d0

        numeroR='1'
        read_file='BASIS/'//trim(basis_file)//'/'//trim(adjustl(numeroR))//'/integral_1e_2e.dat'
        read_file=trim(read_file)
        open(40,file=read_file,status='old',form='formatted')
        !open(50,file='essai.dat',status='unknown',form='formatted')

        Log_SO=.FALSE.
        Log_V1=.FALSE.
        Log_T1=.FALSE.
        Log_mu=.FALSE.
        Log_mux=.FALSE.
        Log_muy=.FALSE.
        write(*,*)
        write(*,*)'Reading integrals ...'
        do ii=1,9999999
            read(40,'(t4,5a7)',end=100)cable1,cable2,cable3,cable4,cable5
            !write(*,*)cable1,'|',cable2,'|',cable3,'|',cable4,'|',cable5,'|',trim(cable4)
            if (trim(cable4)=='    0' .and. trim(cable5)=='    0' .and. Log_SO==.FALSE.) then
                !write(*,*)ii
                read(cable1,*)table1
		
                num=table1
                Log_SO=.TRUE.
                read(40,'(t3,d19.12,2i4)')col1,i,j
                do k=1,num

                    read(40,'(t3,d19.12,2i4)')col1,i,j
                    !write(*,*)'test!!!!!!!!!!!',k,col1,i,j
                    !write(*,*)col1,'|',col2,'|',col3
                    SO(i,j)=col1
                    SO(j,i)=SO(i,j)
                end do
            else if (trim(cable4)=='    0' .and. trim(cable5)=='    1' .and. Log_T1==.FALSE.) then
                !write(*,*)i
                read(cable1,*)table1
		
                num=table1
                Log_T1=.TRUE.
                read(40,'(t3,d19.12,2i4)')col1,i,j
                do k=1,num
                    read(40,'(t3,d19.12,2i4)')col1,i,j
                    !write(*,*)col1,'|',col2,'|',col3
                    T1(i,j)=col1
                    T1(j,i)=T1(i,j)
                end do
            else if (trim(cable4)=='    0' .and. trim(cable5)=='    2' .and. Log_V1==.FALSE.) then
                !write(*,*)i
                read(cable1,*)table1
		
                num=table1
                Log_V1=.TRUE.
                read(40,'(t3,d19.12,2i4)')col1,i,j
                do k=1,num
                    read(40,'(t3,d19.12,2i4)')col1,i,j
                    !write(*,*)col1,'|',col2,'|',col3
                    V1(i,j)=col1
                    V1(j,i)=V1(i,j)
                end do
            !!! Modifier ici avec 1,0 pour mux et 1,1 pour muy

            !!!!!!!11

            else if (trim(cable4)=='    1' .and. trim(cable5)=='    0' .and. Log_mux==.FALSE.) then
                !write(*,*)i
                read(cable1,*)table1
		
                num=table1
                !		write(*,*)'mux'
                Log_mux=.TRUE.
                read(40,'(t3,d19.12,2i4)')col1,i,j
                do k=1,num
                    read(40,'(t3,d19.12,2i4)')col1,i,j
                    !			write(*,*)col1,'|',col2,'|',col3
                    mux(i,j)=col1
                    mux(j,i)=mux(i,j)
                end do

            else if (trim(cable4)=='    1' .and. trim(cable5)=='    1' .and. Log_muy==.FALSE.) then
                !write(*,*)i
                read(cable1,*)table1
		
                num=table1
                !		write(*,*)'muy'
                Log_muy=.TRUE.
                read(40,'(t3,d19.12,2i4)')col1,i,j
                do k=1,num
                    read(40,'(t3,d19.12,2i4)')col1,i,j
                    !			write(*,*)col1,'|',col2,'|',col3
                    muy(i,j)=col1
                    muy(j,i)=muy(i,j)
                end do

            !1111111111111!!!!!!!!!!!!!!!!!!!!!

            else if (trim(cable4)=='    1' .and. trim(cable5)=='    2' .and. Log_mu==.FALSE.) then
                !write(*,*)i
                read(cable1,*)table1
		
                num=table1
                !		write(*,*)'muz'
                Log_mu=.TRUE.
                read(40,'(t3,d19.12,2i4)')col1,i,j
                do k=1,num
                    read(40,'(t3,d19.12,2i4)')col1,i,j
                    !			write(*,*)col1,'|',col2,'|',col3
                    mu(i,j)=col1
                    mu(j,i)=mu(i,j)
                end do
            else if (trim(cable4)=='    3' .and. trim(cable5)=='    0') then
                cable4='ouhouh'
                !write(*,*)iintent(out)
                read(cable1,*)table1
                num=table1
                do f2=1,num
                    read(40,'(t3,d19.12,4i4)')col1,i,j,k,l
                    !write(*,*)col1,'|',col2,'|',col3,'|',col4,'|',col5
                    vee(i,j,k,l)=col1
                    Vee(j,i,k,l)=Vee(i,j,k,l)
                    Vee(j,i,l,k)=Vee(i,j,k,l)
                    Vee(i,j,l,k)=Vee(i,j,k,l)
                    Vee(k,l,i,j)=Vee(i,j,k,l)
                    Vee(l,k,j,i)=Vee(i,j,k,l)
                    Vee(k,l,j,i)=Vee(i,j,k,l)
                    Vee(l,k,i,j)=Vee(i,j,k,l)
                end do
            end if
        end do
100 continue
    !close(50)
    close(40)
    !
    !FIN
    !BLOC DE LECTURE DU FICHIER COLOMBUS
    !
    !
    !RECONSTRUCTION DES MATRICES
    !
    !do i=1,norb
    !do j=1,norb
    !	if (i>j) then
    !		SO(j,i)=SO(i,j)
    !		mu(j,i)=mu(i,j)
    !		muy(j,i)=muy(i,j)
    !		mux(j,i)=mux(i,j)
    !		T1(j,i)=T1(i,j)
    !		V1(j,i)=V1(i,j)
    !	end if
    !end do
    !end do
    !
    !do i=1,norb
    !do j=1,norb
    !do k=1,norb
    !do l=1,norb
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Removed because it is useless MP
    !	if (Vee(i,j,k,l).eq.0) then
    !		Vee(i,j,k,l)=Vee(j,i,l,k)
    !		Vee(i,j,k,l)=Vee(j,i,k,l)
    !		Vee(i,j,k,l)=Vee(i,j,l,k)
    !		Vee(i,j,k,l)=Vee(k,l,i,j)
    !		Vee(i,j,k,l)=Vee(l,k,j,i)
    !		Vee(i,j,k,l)=Vee(k,l,j,i)
    !		Vee(i,j,k,l)=Vee(l,k,i,j)
    !	end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !	if (Vee(i,j,k,l).ne.0) then
    !		Vee(j,i,k,l)=Vee(i,j,k,l)
    !		Vee(j,i,l,k)=Vee(i,j,k,l)
    !		Vee(i,j,l,k)=Vee(i,j,k,l)
    !		Vee(k,l,i,j)=Vee(i,j,k,l)
    !		Vee(l,k,j,i)=Vee(i,j,k,l)
    !		Vee(k,l,j,i)=Vee(i,j,k,l)
    !		Vee(l,k,i,j)=Vee(i,j,k,l)
    !	end if
    !end do
    !end do
    !end do
    !end do




    write(*,*) 'Done.'
    write(*,*)
end subroutine read_integrals


!****************************************************************************************
!****************************************************************************************
subroutine read_coeff_lcao(lcao,norb,basis_file)
    !read coefficients
    !****************************************************************************************
    !****************************************************************************************


    character(len=6)::curseur
    integer(kind=int_4),intent(in)::norb
    integer(kind=int_4)::diviseur,reste,pointeur,i,j,k,inR,nbcol
    real(kind=real_8)::coeff1,coeff2,coeff3,coeff4
    real(kind=real_8),dimension(norb,norb),intent(out)::lcao
    character(len=32),intent(in)::basis_file
    character(len=64)::read_file,numeroR

    write(*,*)'Reading LCAO coefficients ...'
    lcao=0.d0
    nbcol=4
    inR=1
    write(numeroR,'(i4)')inR
    read_file='BASIS/'//trim(basis_file)//'/'//trim(adjustl(numeroR))//'/mocoef_scf.sp'
    read_file=trim(read_file)
    open(41,file=read_file,status='old',form='formatted')
    diviseur=norb/nbcol
    reste=norb-(nbcol*diviseur)
    do i=1,99999999
        read(41,'(a6)',end=200)curseur
        if(curseur=='mocoef') then

            read(41,'(a7)')curseur !Sert juste à passer à la ligne suivante

            do j=1,norb
                pointeur=0
                do k=1,diviseur
                    read(41,'(4d20.12)')coeff1,coeff2,coeff3,coeff4
                    lcao(1+pointeur,j)=coeff1
                    lcao(2+pointeur,j)=coeff2
                    lcao(3+pointeur,j)=coeff3
                    lcao(4+pointeur,j)=coeff4
                    !write(99,'(4d20.12)')coeff1,coeff2,coeff3,coeff4
                    pointeur=pointeur+4
                end do
                select case(reste)
                    case(1)
                        read(41,'(d20.12)')coeff1
                        lcao(1+pointeur,j)=coeff1
	                                                                                !write(99,'(d20.12)')coeff1,coeff2
                    case(2)
                        read(41,'(2d20.12)')coeff1,coeff2
                        lcao(1+pointeur,j)=coeff1
                        lcao(2+pointeur,j)=coeff2
	                                                                                !write(99,'(2d20.12)')coeff1,coeff2
                    case(3)
                        read(41,'(3d20.12)')coeff1,coeff2,coeff3
                        lcao(1+pointeur,j)=coeff1
                        lcao(2+pointeur,j)=coeff2
                        lcao(3+pointeur,j)=coeff3
		                                                                !write(99,'(3d20.12)')coeff1,coeff2
                end select

            end do
		
        end if
    end do

200 continue
    close(41)



    write(*,*) 'Done.'
    write(*,*)
end subroutine read_coeff_lcao

!****************************************************************************************
subroutine read_coeff_lcao_mc(lcao,norb,basis_file,nR)
    !read coefficients
    !****************************************************************************************
    !****************************************************************************************


    character(len=6)::curseur
    integer(kind=int_4),intent(in)::norb,nR
    integer(kind=int_4)::diviseur,reste,pointeur,i,j,k,nbcol
    real(kind=real_8)::coeff1,coeff2,coeff3,coeff4
    real(kind=real_8),dimension(:,:),intent(out)::lcao
    character(len=32),intent(in)::basis_file
    character(len=64)::read_file,numeroR

    !test
    real(kind=real_8),dimension(norb)::norme

    lcao=0.d0
    nbcol=3
    read_file='BASIS/'//trim(basis_file)//'/'//trim(adjustl(numeroR))//'/mocoef_mc.sp'
    read_file=trim(read_file)
    open(41,file=read_file,status='old',form='formatted')
    diviseur=norb/nbcol
    reste=norb-(nbcol*diviseur)
    do i=1,99999999
        read(41,'(a6)',end=200)curseur
        if(curseur=='mocoef') then
            read(41,'(a7)')curseur !Sert juste à passer à la ligne suivante
            do j=1,norb
                pointeur=0
                do k=1,diviseur
                    read(41,'(3d25.15)')coeff1,coeff2,coeff3
                    lcao(1+pointeur,j)=coeff1
                    lcao(2+pointeur,j)=coeff2
                    lcao(3+pointeur,j)=coeff3
                    !write(*,'(4d20.12)')coeff1,coeff2,coeff3,coeff4
                    pointeur=pointeur+nbcol
                end do
                select case(reste)
                    case(1)
                        read(41,'(d25.15)')coeff1
                        lcao(1+pointeur,j)=coeff1
                    case(2)
                        read(41,'(2d25.15)')coeff1,coeff2
                        lcao(1+pointeur,j)=coeff1
                        lcao(2+pointeur,j)=coeff2
		                                                                !write(*,'(2d20.12)')coeff1,coeff2
                end select

            end do
		
        end if
    end do
200 continue
    close(41)

    !write(*,*)'LCAO'
    !do i=1,norb
    !	write(*,'(11(e9.2,2X))')(lcao(i,j,1),j=1,norb)
    !end do

    !do j=1,norb
    !do i=1,norb
    !	norme(j)=norme(j) + (lcao(i,j,1)*(lcao(i,j,1)))
    !end do
    !	norme(j)=dsqrt(norme(j))
    !end do
    !write(*,'(11(e9.2,2X))')(norme(i),i=1,norb)

    write(*,*) 'Coefficents LCAO lus'

end subroutine read_coeff_lcao_mc


!****************************************************************************************
!****************************************************************************************
subroutine aointTOmoint(SO,T1,V1,mux,muy,muz,vee,lcao,norb,vee2)
    ! construction des intégtrales moléculaires à partir des intégrales atomiques et des
    ! coefficients LCAO
    !****************************************************************************************
    !****************************************************************************************


    use basics
    use math !chop
    !use blas95
    !use f95_precision

    implicit none
    integer(kind=int_4), intent(in) :: norb
    real(kind=real_8), intent(inout), dimension(norb,norb) :: SO,T1,V1,mux,muy,muz
    real(kind=real_8), dimension(norb,norb)::SOtemp,T1temp,V1temp,mutemp
    real(kind=real_8), intent(inout), dimension(:,:,:,:) :: Vee
    real(kind=real_8),  dimension(norb,norb,norb,norb), intent(out) :: vee2
    real(kind=real_8), intent(in), dimension(norb,norb) :: lcao
    real(kind=real_8), dimension(norb,norb) :: lcaoCon,temp
    real(kind=real_8) :: tmp,tempo
    integer(kind=int_4) :: i,j,k,l,r,s,t,u

    ! Variables used for System Clock
    real(kind=real_8) :: tdeb,tfin

    ! Pour la parallelisation
    integer NTHREADS, TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
    LOGICAL :: file_exists

    write(*,*)'Integral transformations: AOs to MOs'
    write(*,'(A24)')'One-electron integrals'
    lcaoCon(:,:)=transpose(lcao(:,:))
    ! call cpu_time(tdeb)
    call gemm(SO,lcao,temp)
    call gemm(lcaoCon,temp,SO)

    call gemm(T1,lcao,temp)
    call gemm(lcaoCon,temp,T1)

    call gemm(V1,lcao,temp)
    call gemm(lcaoCon,temp,V1)

    call gemm(mux,lcao,temp)
    call gemm(lcaoCon,temp,mux)

    call gemm(muy,lcao,temp)
    call gemm(lcaoCon,temp,muy)

    call gemm(muz,lcao,temp)
    call gemm(lcaoCon,temp,muz)

    !Ecriture des intégrales en base moléculaire dans le fichier fort.100
    !write(100,*)'Recouvrement SO'
    !do i=1,norb
    !do j=1,norb
    !	if (i<=j .and. SO(i,j,1).ne.0) then
    !		write(100,*)SO(i,j,1),i,j
    !	end if
    !end do
    !end do

    !write(100,*)'Cinétique T1'
    !do i=1,norb
    !do j=1,norb
    !	if (i<=j .and. T1(i,j,1).ne.0) then
    !		write(100,*)T1(i,j,1),i,j
    !	end if
    !end do
    !end do

    !write(100,*)'V1'
    !do i=1,norb
    !do j=1,norb
    !	if (i<=j .and. V1(i,j,1).ne.0) then
    !		write(100,*)V1(i,j,1),i,j
    !	end if
    !end do
    !end do

    !write(100,*)'mu'
    !do i=1,norb
    !do j=1,norb
    !	if (i<=j .and. mu(i,j,1).ne.0) then
    !		write(100,*)mu(i,j,1),i,j
    !	end if
    !end do
    !end do
    !Fin écriture

    ! call cpu_time(tfin)
    !write(*,*) "durée intégrales mono-électronique:  ", tfin-tdeb

    !intégrales bi-électroniques
    ! call cpu_time(tdeb)
    open(321,file='vee2_results.txt',status='unknown')
    open(654,file='vee_results.txt',status='unknown')
    write(*,'(A24)')'Two-electron integrals'
    Vee2=0.d0
    INQUIRE(FILE="2e.bin", EXIST=file_exists)
    write(*,*) 'file_exists=',file_exists

    if (file_exists) then
        open(unit=199,file='2e.bin', form='UNFORMATTED')
        read(199)Vee2
        close(unit=199)

    else

        !!!!!!$OMP PARALLEL DO REDUCTION(+:tempo)


        !TID = OMP_GET_THREAD_NUM()
        !NTHREADS = OMP_GET_NUM_THREADS()
        !write(*,*) "THREAD No=",TID,"/",NTHREADS,"Doing i=",i

        !$OMP PARALLEL DO SHARED(Vee2,lcaoCon,Vee,lcao) PRIVATE(NTHREADS, TID,i,j,k,l,r,s,t,u,tmp)
        do l=1,norb
            write(*,*) "In aointTOmoint, l=",l,"/",norb
            do k=1,norb
                do j=1,norb
                    do i=1,norb
                        do u=1,norb
                            do t=1,norb
                                do s=1,norb
                                    do r=1,norb
                                        Vee2(i,j,k,l)=Vee2(i,j,k,l)+lcaoCon(i,r)*lcaoCon(j,s)*Vee(r,s,t,u)*lcao(t,k)*lcao(u,l) !!!à optimiser
                                    !write(321,*) Vee2(i,j,k,l)
                                    !write(654,*) Vee(r,s,t,u)
                                    enddo
                                enddo
                            enddo
                        enddo

                    !Vee2(i,j,k,l) = tempo

                    enddo
                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO

        close(321)
        close(654)
        write(*,*)'Done.'
        write(*,*)

        open(unit=199,file='2e.bin', form='UNFORMATTED')
        write(199)Vee2
        close(unit=199)
    end if
end subroutine aointTOmoint

subroutine aointTOmointn5(SO,T1,V1,mux,muy,muz,vee,lcao,norb,vee2)
    ! construction des intégtrales moléculaires à partir des intégrales atomiques et des
    ! coefficients LCAO
    !****************************************************************************************
    !****************************************************************************************


    use basics
    use math !chop
    !use blas95
    !use f95_precision

    implicit none
    integer(kind=int_4), intent(in) :: norb
    real(kind=real_8), intent(inout), dimension(norb,norb) :: SO,T1,V1,mux,muy,muz
    real(kind=real_8), dimension(norb,norb)::SOtemp,T1temp,V1temp,mutemp
    real(kind=real_8), intent(inout), dimension(:,:,:,:) :: Vee
    real(kind=real_8),  dimension(norb,norb,norb,norb), intent(out) :: vee2
    real(kind=real_8), intent(in), dimension(norb,norb) :: lcao
    real(kind=real_8), dimension(norb,norb) :: lcaoCon,temp
    real(kind=real_8) :: tmp,tempo
    integer(kind=int_4) :: i,j,k,l,r,s,t,u
    real(kind=real_8), dimension(norb,norb) :: X,Y
    real(kind=real_8),  dimension(norb,norb,norb,norb) :: temp2e

    ! Variables used for System Clock
    real(kind=real_8) :: tdeb,tfin

    ! Pour la parallelisation
    integer NTHREADS, TID, OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
    LOGICAL :: file_exists

    write(*,*)'Integral transformations: AOs to MOs with n5 algo'
    write(*,'(A24)')'One-electron integrals'
    lcaoCon(:,:)=transpose(lcao(:,:))
    ! call cpu_time(tdeb)
    call gemm(SO,lcao,temp)
    call gemm(lcaoCon,temp,SO)

    call gemm(T1,lcao,temp)
    call gemm(lcaoCon,temp,T1)

    call gemm(V1,lcao,temp)
    call gemm(lcaoCon,temp,V1)

    call gemm(mux,lcao,temp)
    call gemm(lcaoCon,temp,mux)

    call gemm(muy,lcao,temp)
    call gemm(lcaoCon,temp,muy)

    call gemm(muz,lcao,temp)
    call gemm(lcaoCon,temp,muz)



    open(321,file='vee2_results.txt',status='unknown')
    open(654,file='vee_results.txt',status='unknown')
    write(*,'(A24)')'Two-electron integrals'
    Vee2=0.d0
    INQUIRE(FILE="2e.bin", EXIST=file_exists)
    write(*,*) 'file_exists=',file_exists

    if (file_exists) then
        open(unit=199,file='2e.bin', form='UNFORMATTED')
        read(199)Vee2
        close(unit=199)

    else

        !!!!!!$OMP PARALLEL DO REDUCTION(+:tempo)


        !TID = OMP_GET_THREAD_NUM()
        !NTHREADS = OMP_GET_NUM_THREADS()
        !write(*,*) "THREAD No=",TID,"/",NTHREADS,"Doing i=",i
!!$OMP PARALLEL DO SHARED(temp2e,Vee2,lcaoCon,Vee,lcao) PRIVATE(NTHREADS, TID,i,j,k,l,X,Y)
!        do l=1,norb
!            do k=1,norb
!                do j=1,norb
!                    do i=1,norb
!                        X(i,j)=Vee(i,j,k,l)
!                    end do
!                end do
!                call gemm(lcaoCon,X,Y)      !????? les 2 indices i,j sont du cote bra. Ce sont elles que l'on transforme ici
!                call gemm(lcaoCon,Y,X)      ! ?????
!               do j=1,norb
!                    do i=1,norb
!                        temp2e(k,l,i,j)=X(i,j) !     ???? ou est-ce temp2e(i,j,k,l)=X(i,j) ????
!                    end do
!                end do
!            end do
!        end do

        !$OMP PARALLEL DO SHARED(temp2e,Vee2,lcaoCon,Vee,lcao) PRIVATE(NTHREADS, TID,i,j,k,l,X,Y)
        do l=1,norb
            do k=1,norb
                do j=1,norb
                    do i=1,norb
                        X(i,j)=Vee(i,j,k,l)
                    end do
                end do
                call gemm(lcaoCon,X,Y)       
                call gemm( Y,lcao, X)
                !call gemm(lcao,X,Y)
                !call gemm(lcao,Y,X)
                do j=1,norb
                    do i=1,norb
                        temp2e(i,j,k,l)=X(i,j)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
        X=0.d0
        Y=0.d0

        !$OMP PARALLEL DO SHARED(temp2e,Vee2,lcaoCon,Vee,lcao) PRIVATE(NTHREADS,Y,X,TID,i,j,k,l)
        do l=1,norb
            write(*,*) "In aointTOmoint, l=",l,"/",norb
            do k=1,norb

                do j=1,norb
                    do i=1,norb
                        X(i,j)=temp2e(k,l,i,j)
                    end do
                end do

                call gemm(lcaoCon,X,Y)       
                call gemm( Y,lcao, X)
                !call gemm(lcaoCon,X,Y)
                !call gemm(lcaoCon,Y,X)
                do j=1,norb
                    do i=1,norb
                        vee2(i,j,k,l)=X(i,j)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        close(321)
        close(654)
        write(*,*)'Done.'
        write(*,*)

        open(unit=199,file='2e.bin', form='UNFORMATTED')
        write(199)Vee2
        close(unit=199)
    end if
end subroutine aointTOmointn5

!****************************************************************************************
!****************************************************************************************
subroutine write_SO(SO,norb,inR,NOM)
    !  écriture des matrices SO dans le fichier NOM.dat
    !****************************************************************************************
    !****************************************************************************************

    integer(kind=int_4),intent(in)::norb,inR
    integer:: i,j
    character(len=*),intent(in)::NOM
    real(kind=real_8), intent(in), dimension(:,:) :: SO



    open(50,file=trim(adjustl(NOM))//'.dat',status='unknown',form='formatted')
    do i=1,norb
        do j=1,norb
            write(50,998) SO(i,j),i,j
        enddo
    enddo
    close(50)

998 FORMAT(E18.5,I4,I4,I4) !écriture des éléments de matrice

end subroutine write_SO

!****************************************************************************************
!****************************************************************************************
subroutine lecture_gaus(Norb,basis_file,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord)
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
    real(kind=real_8),allocatable,dimension(:)::X,Y,Z
    character(len=3),allocatable,dimension(:)::MTYPE
    real,allocatable,dimension(:)::CHG
    character::NEANT
    character(len=21),parameter::FMT1 ='(t1,d16.8,t1,d16.8)'
    character(len=64)::read_file
    character(len=13)::pointeur
    integer::nligne,reste


    ZET=0.d0
    ETA=0.d0
    coord=0.d0
    ICONU=0.d0
    NF=0.d0
    LMNP1=0.d0
    NCONS=0.d0

    write(*,*)'Début de lecture_gauss'
    read_file='BASIS/'//trim(basis_file)//'/1/argosls.sp'
    read_file=trim(read_file)
    open(47,file=read_file,status='old',form='formatted')
    do i=1,3
        read(47,*)NEANT
    end do
    read(47,*)NGEN,NS,NAORDS,NCONS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT
    !This format was working on old version of columbus, but is not working on Col7.0
    !Using * will work with both  version
    !read(47,'(t2,15i3)')NGEN,NS,NAORDS,NCONS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT
    read(47,*)NEANT
    read(47,'(t2,i3)')NDPT
    if (NDPT.ne.0) then
        do i=1,NDPT
            read(47,*)NEANT
        end do
    end if
    if (NAORDS.ne.0) then
        do i=1,NAORDS
            read(47,*)NEANT
        end do
    end if
    if (NGCS.ne.0) then
        allocate(ICSU(NGCS))
        do i=1,NGCS
            read(47,'(t2,i5)')ICSU(i)
            if (ICSU(i).ne.0) then
                do j=1,ICSU(i)
                    read(47,*)NEANT
                end do
            end if
        end do
    end if
    !write(*,*) "test on line 645 of read_columbus, NCONS=",NCONS
    allocate(NRCR(NCONS))
    allocate(ZETtmp(32,NCONS),ETAtmp(32,NCONS))
    ZETtmp=0.d0
    ETAtmp=0.d0
    do i=1,NCONS
        !DIFFÉRENCE AVEC MASTER -JN
        read(47,'(t2,3i5)')ICONU(i),LMNP1(i),NRCR(i)
        !write(99,*)'ICONU(',i,'),LMNP1(',i,'),NRCR(',i,')',ICONU(i),LMNP1(i),NRCR(i)
        if (LMNP1(i).eq.1) then
            write(99,*)' s',ICONU(i),NRCR(i)
        else if   (LMNP1(i).eq.2) then
            write(99,*)' p',ICONU(i),NRCR(i)
        end if
        do j=1,ICONU(i)
            read(47,'(t1,d16.8,d17.8)')ZETtmp(j,i),ETAtmp(j,i)
            !WRITE(99,*)'ZETtmp(',j,',',i,')',ETAtmp(j,i),ZETtmp(j,i),ETAtmp(j,i)
            WRITE(99,'(f16.10,2x,f16.10)')ZETtmp(j,i),ETAtmp(j,i)
        end do
    end do
    nprimmax=Maxval(iconu)
    ZET=ZETtmp(1:nprimmax,:)
    ETA=ETAtmp(1:nprimmax,:)
    deallocate(ETAtmp,ZETtmp)
    !Pas de partie 9 / A voir pour le futur

    allocate(x(NS),y(NS),z(NS),MTYPE(NS),NC(NS),CHG(NS),IGEN(NGEN))
    allocate(MCONS(33),IGCS(33),MCRS(33))!Ici 33 est le nombre max d'OA par atome)
    counter1=0 !compteur du nombre de contraction set
    do IS=1,NS
        read(47,'(t2,A3,2I3,F3.0)')MTYPE(IS),NF(IS),NC(IS),CHG(IS)
        !	write(*,'(t2,A3,2I3,F3.0)')MTYPE(IS),NF(IS),NC(IS),CHG(IS)
        do j=1,NC(IS)
            read(47,'(t2,3f14.8)')X(j),Y(j),Z(j) !center position
            write(99,'(t2,3f14.8)')X(j),Y(j),Z(j) !center position
            !	        write(*,'(t2,3f14.8)')X(j),Y(j),Z(j)
            coord(IS,:)=(/X(j),Y(j),Z(j)/)
        end do
        !	write(*,*)('NF')
        !	write(*,*)(NF(IS))
        IF(NC(IS) .NE. 1) THEN
            DO J = 1, NGEN
                read(47,*)IGEN(J)!A terminer
            ENDDO
        ENDIF
        DO j = 1, NF(IS)
            read(47,'(t2,2i3)')MCONS(j), IGCS(j)
            !		write(*,'(t2,2i3)')MCONS(j), IGCS(j)
            counter1=counter1+1
            !		write(*,*)'counter1',counter1
            SetCenter(counter1)=IS
        ENDDO	
        IF(NCRS .NE. 0) THEN
            read(47,*)MCRS(IS)
        !		write(*,*)MCRS(IS)
        ENDIF
    end do



!counter1=0
!Do IS=1,NS
!Do j=1,NF(IS)
!	select case(LMNP1(counter1))
!	case(1)
!		write(*,'(5I3)')IS,counter1,0,0,0
!	case(2)
!		write(*,'(5I3)')IS,counter1,1,0,0
!		write(*,'(5I3)')IS,counter1,0,1,0
!		write(*,'(5I3)')IS,counter1,0,0,1
!	case(3)
!		write(*,'(5I3)')IS,counter1,2,0,0
!		write(*,'(5I3)')IS,counter1,0,2,0
!		write(*,'(5I3)')IS,counter1,0,0,2
!		write(*,'(5I3)')IS,counter1,1,1,0
!		write(*,'(5I3)')IS,counter1,1,0,1
!		write(*,'(5I3)')IS,counter1,0,1,1
!	end select
!end do
!end do

end subroutine lecture_gaus
!!****************************************************************************************
!!****************************************************************************************
!subroutine lecture_gaus(basis_file,ETA,ZET,NCONSI,ICONU,NF,LMNP1,coord,znuc,maxprim,NSI,CMap)
!!  UPDATE DE MAI 2013 dans le share
!!  lecture des informations concernant la base atomique
!!
!!
!!	output:	
!!		-ETA(:,:) 
!!		-ZET(:,:)
!!		-ICONU(:)
!!		-NF:: (NS)
!!		-LMNP1(:)
!!
!!****************************************************************************************
!!****************************************************************************************

!!I/O variables
! character(len=32),intent(in)::basis_file
!integer(kind=int_4), intent(in) :: maxprim,NSI,NCONSI
!real(kind=real_8), dimension(maxprim,NCONSI), intent(out) :: ZET,ETA
!real(kind=real_8), dimension(NSI,3), intent(out) :: coord
!real(kind=real_8), dimension(NSI), intent(out) :: znuc
!integer(kind=int_4),dimension(NCONSI),intent(out):: ICONU,LMNP1
!integer(kind=int_4),dimension(NSI),intent(out):: NF
!integer(kind=int_4),dimension(NCONSI):: CMap ! center of the contraction set



!!variables internes
!integer(kind=int_4)::NGEN,NAORDS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT,i,j
!integer(kind=int_4)::NDPT,IS,NS,NCONS
!integer(kind=int_4)::counter1,ll
!integer(kind=int_4), allocatable,dimension(:)::ICSU,NRCR,NC,MCONS,IGCS,MCRS,IGEN
!real(kind=real_8), allocatable,dimension(:,:)::ZETtmp,ETAtmp
!integer(kind=int_4)::nprimmax
!real(kind=real_8),allocatable,dimension(:)::X,Y,Z
! character(len=3),allocatable,dimension(:)::MTYPE
!real,allocatable,dimension(:)::CHG
! character::NEANT
! character(len=21),parameter::FMT1 ='(t1,d16.8,t1,d16.8)'
! character(len=64)::read_file
! character(len=13)::pointeur
!integer::nligne,reste


!ZET=0.d0
!ETA=0.d0
!coord=0.d0
!ICONU=0.d0
!NF=0.d0
!LMNP1=0.d0
!NCONS=0.d0

!write(*,*)
!write(*,*)'Reading basis function parameters ...'
!read_file='BASIS/'//trim(basis_file)//'/1/argosls.sp'
!read_file=trim(read_file)
!open(47,file=read_file,status='old',form='formatted')
!do i=1,3
!	read(47,*)NEANT
!end do
!read(47,'(t2,15i3)')NGEN,NS,NAORDS,NCONS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT
!write(*,'(t2,i3)')NAORDS

!read(47,*)NEANT
!read(47,'(t2,i3)')NDPT
!if (NDPT.ne.0) then
!	do i=1,NDPT
!		read(47,*)NEANT
!	end do
!end if
!if (NAORDS.ne.0) then
!	do i=1,NAORDS
!		read(47,*)NEANT
!	end do
!end if
!if (NGCS.ne.0) then
!	allocate(ICSU(NGCS))
!	do i=1,NGCS
!		read(47,'(t2,i5)')ICSU(i)
!		if (ICSU(i).ne.0) then
!			do j=1,ICSU(i)
!				read(47,*)NEANT
!			end do
!		end if
!	end do
!end if
!allocate(NRCR(NCONS))
!allocate(ZETtmp(32,NCONS),ETAtmp(32,NCONS))
!ZETtmp=0.d0
!ETAtmp=0.d0
!do i=1,NCONS
!	read(47,'(t2,3i5)')ICONU(i),LMNP1(i),NRCR(i)
!!	write(*,'(t2,3i5)')ICONU(i),LMNP1(i),NRCR(i)
!	do j=1,ICONU(i)
!		read(47,'(t1,d16.8,d17.8)')ZETtmp(j,i),ETAtmp(j,i)
!!		WRITE(*,'(t1,d16.8,d17.8)')ZETtmp(j,i),ETAtmp(j,i)
!	end do	
!end do
!nprimmax=Maxval(iconu)
!ZET=ZETtmp(1:nprimmax,:)
!ETA=ETAtmp(1:nprimmax,:)
!deallocate(ETAtmp,ZETtmp)
!!Pas de partie 9 / A voir pour le futur

!allocate(x(NS),y(NS),z(NS),MTYPE(NS),NC(NS),CHG(NS),IGEN(NGEN))
!allocate(MCONS(33),IGCS(33),MCRS(33))!Ici 33 est le nombre max d'OA par atome)
!counter1=0 !compteur du nombre de contraction set
!do IS=1,NS
!	read(47,'(t2,A3,2I3,F3.0)')MTYPE(IS),NF(IS),NC(IS),CHG(IS)
!!	write(*,'(t2,A3,2I3,F3.0)')MTYPE(IS),NF(IS),NC(IS),CHG(IS)
!	do j=1,NC(IS)
!		read(47,'(t2,3f14.8)')X(j),Y(j),Z(j) !center position
!!	        write(*,'(t2,3f14.8)')X(j),Y(j),Z(j)
!		coord(IS,:)=(/X(j),Y(j),Z(j)/)
!		znuc(IS)=CHG(IS)
!	end do
!        IF(NC(IS) .NE. 1) THEN
!		DO J = 1, NGEN
!			read(47,*)IGEN(J)!A terminer
!		ENDDO
!	ENDIF
!        DO j = 1, NF(IS)
!        	read(47,'(t2,2i3)')MCONS(j), IGCS(j)
!!		write(*,'(t2,2i3)')MCONS(j), IGCS(j)
!		counter1=counter1+1
!write(*,'(A17,I2,A4,I2)')'Contraction set ',counter1,' of ',NCONS
!		CMap(counter1)=IS
!!		write(*,*)'IS',IS,'j',j
!        ENDDO	
!        IF(NCRS .NE. 0) THEN
!	        read(47,*)MCRS(IS)
!!		write(*,*)MCRS(IS)
!        ENDIF
!end do
!write(*,*)'Done.'
!write(*,*)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Prints the identification of AOs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*)'Description of the atomic orbital basis set :'
!write(*,*)
!write(*,'(A8,A11,A16)')'  AO  ',' Center ',' Contraction Set'

!counter1=0
!do i=1,NCONS
!!define number of degeneracy

!If(LMNP1(i).eq.1) then
!			ll=1
!else If(LMNP1(i).eq.2) then
!			ll=3
!else If(LMNP1(i).eq.3) then
!			ll=5
!else If(LMNP1(i).eq.4) then
!			ll=10 ! the actual number is 7
!else If(LMNP1(i).eq.5) then
!			ll=15 ! should be 9
!end if

!	do j=1,ll
!		counter1=counter1+1

!		write(*,'(I3,A4,A1,A5,A1,I2,A2,A8,I2)')counter1,' -> ',AOlabel(LMNP1(i)),AOdirection(LMNP1(i),j),' ',CMap(i),'  ',' ',i
!	end do
!end do

!write(*,*)
!do i=1,ns
!write(*,'(A9,I2,A3,A2)')'Center ',i,' : ',MTYPE(i)
!write(*,'(A6,F5.2,A5,F5.2,A5,F5.2)')'x = ',coord(i,1),' y = ',coord(i,2),' z = ',coord(i,3)
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!end subroutine lecture_gaus

!****************************************************************************************
!****************************************************************************************
! Prints Atomic Orbital label
function AOlabel(AO)
    !****************************************************************************************
    !****************************************************************************************

    use basics
    implicit none
    character(len=1) :: AOlabel
    integer(kind=int_4), intent(in) :: AO

    AOlabel='?'
    if(AO.eq.1)then
        AOlabel='s'
    endif
    if(AO.eq.2)then
        AOlabel='p'
    endif
    if(AO.eq.3)then
        AOlabel='d'
    endif
    if(AO.eq.4)then
        AOlabel='f'
    endif
    if(AO.eq.5)then
        AOlabel='g'
    endif

end function

!****************************************************************************************
!****************************************************************************************
! Prints Atomic Orbital direction
function AOdirection(AO,DIR)
    !****************************************************************************************
    !****************************************************************************************

    use basics
    implicit none
    character(len=5) :: AOdirection
    integer(kind=int_4), intent(in) :: AO,DIR
    AOdirection='?'
    if(AO.eq.1)then
        AOdirection=' '
    endif

    if(AO.eq.2)then
        if(DIR.eq.1)then
            AOdirection='z'
        endif
        if(DIR.eq.2)then
            AOdirection='x'
        endif
        if(DIR.eq.3)then
            AOdirection='y'
        endif

    endif

    if(AO.eq.3)then
        if(DIR.eq.1)then
            AOdirection='zz'
        endif
        if(DIR.eq.2)then
            AOdirection='xz'
        endif
        if(DIR.eq.3)then
            AOdirection='yz'
        endif
        if(DIR.eq.4)then
            AOdirection='xx-yy'
        endif
        if(DIR.eq.5)then
            AOdirection='xy'
        endif
    endif

    !!! orbitals over d shell are wrong ... this need to be fixed if one wishes to use that
    if(AO.eq.4)then
        if(DIR.eq.1)then
            AOdirection='xxx'
        endif
        if(DIR.eq.2)then
            AOdirection='yyy'
        endif
        if(DIR.eq.3)then
            AOdirection='zzz'
        endif
        if(DIR.eq.4)then
            AOdirection='xxy'
        endif
        if(DIR.eq.5)then
            AOdirection='xxz'
        endif
        if(DIR.eq.6)then
            AOdirection='yyx'
        endif
        if(DIR.eq.7)then
            AOdirection='yyz'
        endif
        if(DIR.eq.8)then
            AOdirection='zzx'
        endif
        if(DIR.eq.9)then
            AOdirection='zzy'
        endif
        if(DIR.eq.10)then
            AOdirection='xyz'
        endif
    endif

    if(AO.eq.5)then
        if(DIR.eq.1)then
            AOdirection='xxxx'
        endif
        if(DIR.eq.2)then
            AOdirection='yyyy'
        endif
        if(DIR.eq.3)then
            AOdirection='zzzz'
        endif
        if(DIR.eq.4)then
            AOdirection='xxxy'
        endif
        if(DIR.eq.5)then
            AOdirection='xxxz'
        endif
        if(DIR.eq.6)then
            AOdirection='yyyx'
        endif
        if(DIR.eq.7)then
            AOdirection='yyyz'
        endif
        if(DIR.eq.8)then
            AOdirection='zzzx'
        endif
        if(DIR.eq.9)then
            AOdirection='zzzy'
        endif
        if(DIR.eq.10)then
            AOdirection='xxyy'
        endif
        if(DIR.eq.11)then
            AOdirection='xxzz'
        endif
        if(DIR.eq.12)then
            AOdirection='yyzz'
        endif
        if(DIR.eq.13)then
            AOdirection='xxyz'
        endif
        if(DIR.eq.14)then
            AOdirection='yyxz'
        endif
        if(DIR.eq.15)then
            AOdirection='zzxy'
        endif
    endif
end function

!****************************************************************************************
!****************************************************************************************
subroutine lecture_NS(Ncenter,basis_file,NCONS)
    !  lecture du nombre de gaussiennes dans le fichier argosls.sp
    !****************************************************************************************
    !****************************************************************************************

    integer(kind=int_4),intent(out)::Ncenter,NCONS
    character(len=32),intent(in)::basis_file
    character(len=64)::read_file
    integer::i

    integer(kind=int_4)::NGEN,NS,NAORDS

    NCenter=0.d0
    read_file='BASIS/'//trim(basis_file)//'/1/argosls.sp'
    read_file=trim(read_file)
    open(47,file=read_file,status='old',form='formatted')
    do i=1,3
        read(47,*)
    end do
    read(47,'(t2,4i3)')NGEN,NS,NAORDS,NCONS
    NCenter=NS
    close(47)
    write(*,'(A26,I2)')'Number of atomic centers :', Ncenter
    write(*,'(A28,I3)')'Number of contraction sets :',NCONS
end subroutine lecture_NS


!****************************************************************************************
!****************************************************************************************
subroutine scan_dimension_argos(basis_file,NCONS,maxprim,NS)
    !
    !  lecture des coordonnées nucléaires (Z uniquement)
    !
    !	input: 	-basis_file::
    !
    !	output:	-NCONS
    !		-maxprim
    !		-NS
    !
    !****************************************************************************************
    !****************************************************************************************

    !  I/O variables
    character(len=32),intent(in)::basis_file
    integer(kind=int_4),intent(out)::NCONS,maxprim,NS
    !
    !  internal variables
    !variables internes
    integer(kind=int_4)::NGEN,NAORDS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT,i,j
    integer(kind=int_4)::NDPT,IS
    !integer, allocatable,dimension(:)::ICSU,ICONU,NF,NC,IGEN
    integer(kind=int_4), allocatable,dimension(:)::ICSU,ICONU,LMNP1,NRCR,NF,NC,MCONS,IGCS,MCRS,IGEN
    real(kind=real_8),allocatable,dimension(:)::X,Y,Z
    character(len=3),allocatable,dimension(:)::MTYPE
    real,allocatable,dimension(:)::CHG
    character::NEANT
    character(len=21),parameter::FMT1 ='(t1,d16.8,t1,d16.8)'
    character(len=64)::read_file,numeroR

    NCONS=0.d0
    maxprim=0.d0
    NS=0.d0

    numeroR='1'
    read_file='BASIS/'//trim(basis_file)//'/'//trim(adjustl(numeroR))//'/argosls.sp'
    read_file=trim(read_file)
    !write(*,*) "test on lien 1129 of readcolumbus read_file=",read_file
    open(47,file=read_file,status='old',form='formatted')
    do i=1,3
        read(47,*)NEANT
	                !write(*,*)"test on line 1133 of read_columbus, NEANT =",NEANT
    end do

    read(47,*)NGEN,NS,NAORDS,NCONS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT
    !read(47,'(t2,15i3)')NGEN,NS,NAORDS,NCONS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT
    read(47,*)NEANT
    read(47,'(t2,i3)')NDPT

    !write(*,*) "test on line 1138 of read_columbus, NS=",NS
    !write(*,*) "NGEN,NS,NAORDS,NCONS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT=",NGEN,NS,NAORDS,NCONS,NGCS,ITOL,ICUT,AOINTS,ONLY1E,INRM,NCRS,L1REC,L2REC,AOINT2,FSPLIT
    if (NDPT.ne.0) then
        do i=1,NDPT
            read(47,*)NEANT
        end do
    end if
    if (NAORDS.ne.0) then
        do i=1,NAORDS
            read(47,*)NEANT
        end do
    end if
    !write(*,*) "test on line 1148 of read_columbus, NGCS=",NGCS
    if (NGCS.ne.0) then
        allocate(ICSU(NGCS))
        do i=1,NGCS
            read(47,'(t2,i5)')ICSU(i)
            if (ICSU(i).ne.0) then
                do j=1,ICSU(i)
                    read(47,*)NEANT
                end do
            end if
        end do
    end if
    allocate(ICONU(NCONS))
								
    do i=1,NCONS
        read(47,'(t2,i5)')ICONU(i)
        do j=1,ICONU(i)
            read(47,*)NEANT
        end do
    end do
    maxprim=maxval(iconu)
    !Pas de partie 9 / A voir pour le futur

    allocate(x(NS),y(NS),z(NS),MTYPE(NS),NF(NS),NC(NS),CHG(NS),IGEN(NGEN))
    allocate(MCONS(33),IGCS(33),MCRS(33))!Ici 10 est le nombre max d'OA par atome)
    !counter1=0 !compteur du nombre de contraction set

    do IS=1,NS
        read(47,'(t2,A3,2I3,F3.0)')MTYPE(IS),NF(IS),NC(IS),CHG(IS)
        do j=1,NC(IS)
            read(47,'(t2,3f14.8)')X(j),Y(j),Z(j) !center position
        end do
        IF(NC(IS) .NE. 1) THEN
            DO J = 1, NGEN
                read(47,*)IGEN(J)!A terminer
            ENDDO
        ENDIF
        DO j = 1, NF(IS)
            read(47,'(t2,2i3)')MCONS(j), IGCS(j)
        ENDDO	
        IF(NCRS .NE. 0) THEN
            read(47,*)MCRS(IS)
        ENDIF
    end do


    close(47)
    !if (allocated(ICSU)) then
    deallocate(ICSU)
    !end if
    deallocate(ICONU)
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(MTYPE,NF,NC,CHG,IGEN,MCONS,IGCS,MCRS)
!deallocate(ICSU,ICONU,x,y,z,MTYPE,NF,NC,CHG,IGEN,MCONS,IGCS,MCRS)

end subroutine scan_dimension_argos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Ajouter dans un autre module !
! Routines qui affichent les fonctions de base
! Question de voir si tout se passe bien ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine buildphi(phi,lcao,eta,zet,coord,lmnp1,CMap,NS,NCONS,norb,nq,np,maxprim,xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,nx,ny,nz)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    integer(kind=int_4), intent(in) :: NS,NCONS,norb,nq,np,maxprim,nx,ny,nz
    integer(kind=int_4), dimension(NCONS), intent(in) :: LMNP1, CMap
    real(kind=real_8), dimension(NS,3), intent(in) :: coord
    real(kind=real_8), dimension(norb,norb), intent(in) :: lcao
    real(kind=real_8), dimension(maxprim,NCONS), intent(in) :: eta,zet
    real(kind=real_8), intent(in) :: xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin
    complex(kind=comp_16), dimension(np,nq), intent(out) :: phi

    real(kind=real_8) :: x,y,z,Rx,Ry,Rz
    integer(kind=int_4) :: i,j,k,l,r,s,t,u,CN,DIR
    complex(kind=comp_16), dimension(norb) :: chi
    complex(kind=comp_16) :: tempo

    do i = 1, nx
        do j = 1, ny
            do k = 1, nz

                u=0
                do s=1,NCONS
                    do t=1,2*(LMNP1(s)-1)+1
                        u=u+1
                        tempo = dcmplx(0.d0,0.d0)
                        !$OMP PARALLEL DO REDUCTION(+:tempo)
                        do l=1,maxprim
                            tempo=tempo+eta(l,s)*GTO(zet(l,s),xmin+(i-1)*dx,ymin+(j-1)*dy,zmin+(k-1)*dz,coord(CMap(s),1),coord(CMap(s),2),coord(CMap(s),3),LMNP1(s),t)
                        enddo
                        !$OMP END PARALLEL DO
                        chi(u) = tempo
                    enddo

                enddo


                do r=1,nq
                    tempo = dcmplx(0.d0)
                    !$OMP PARALLEL DO REDUCTION(+:tempo)
                    do u=1,norb
                        tempo=tempo+chi(u)*dcmplx(lcao(u,r),0.d0)
                    enddo
                    !$OMP END PARALLEL DO
                    phi((k-1)*nx*ny+(j-1)*nx+i,r) = tempo
                enddo

            enddo
        enddo
    enddo

end subroutine buildphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine buildphik(phik,lcao,eta,zet,coord,lmnp1,CMap,NS,NCONS,norb,nq,np,maxprim,xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,nx,ny,nz)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    integer(kind=int_4), intent(in) :: NS,NCONS,norb,nq,np,maxprim,nx,ny,nz
    integer(kind=int_4), dimension(NCONS), intent(in) :: LMNP1, CMap
    real(kind=real_8), dimension(NS,3), intent(in) :: coord
    real(kind=real_8), dimension(norb,norb), intent(in) :: lcao
    real(kind=real_8), dimension(maxprim,NCONS), intent(in) :: eta,zet
    real(kind=real_8), intent(in) :: xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin
    complex(kind=comp_16), dimension(np,nq), intent(out) :: phik

    real(kind=real_8) :: kx,ky,kz,Rx,Ry,Rz,zeta
    integer(kind=int_4) :: i,j,k,l,r,s,t,u,CN,DIR,nw
    complex(kind=comp_16), dimension(norb) :: chik
    complex(kind=comp_16) :: tempo

    do i = 1, nx
        do j = 1, ny
            do k = 1, nz

                u=0

                do s=1,NCONS
                    do t=1,2*(LMNP1(s)-1)+1
                        u=u+1
                        tempo = dcmplx(0.d0,0.d0)
                        !$OMP PARALLEL DO REDUCTION(+:tempo)
                        do l=1,maxprim
                            tempo=tempo+dcmplx(eta(l,s),0.d0)*GTOKS(zet(l,s),kxmin+(i-1)*dkx,kymin+(j-1)*dky,kzmin+(k-1)*dkz,coord(CMap(s),1),coord(CMap(s),2),coord(CMap(s),3),LMNP1(s),t)
                        enddo
                        !$OMP END PARALLEL DO
                        chik(u) = tempo
                    enddo

                enddo


                do r=1,nq
                    tempo=dcmplx(0.d0,0.d0)
                    !$OMP PARALLEL DO REDUCTION(+:tempo)
                    do u=1,norb
                        tempo=tempo+chik(u)*dcmplx(lcao(u,r),0.d0)
                    enddo
                    !$OMP END PARALLEL DO
                    phik((k-1)*nx*ny+(j-1)*nx+i,r) = tempo
                enddo


            enddo
        enddo
    enddo

end subroutine buildphik
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function GTOKS(zeta,kx,ky,kz,Rx,Ry,Rz,AO,DIR)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use basics
    implicit none
    complex(kind=comp_16) :: GTOKS
    real(kind=real_8), intent(in) :: zeta,kx,ky,kz,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: AO,DIR


    GTOKS = dcmplx(0.d0,0.d0)

    if(AO.eq.1)then
        GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,0)
    endif

    if(AO.eq.2)then
        !p z orbital
        if(DIR.eq.1)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,1)
        endif
        ! p x orbital
        if(DIR.eq.2)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,0)
        endif
        ! p y orbital
        if(DIR.eq.3)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,0)
        endif

    endif

    if(AO.eq.3)then
        ! d z^2 orbital
        if(DIR.eq.1)then
            GTOKS = dcmplx(2.d0,0.d0)*GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,2)-GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,0)-GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,0)
        endif
        ! d xz orbital
        if(DIR.eq.2)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,1)
        endif
        ! d yz orbital
        if(DIR.eq.3)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,1)
        endif
        ! d x^2-y^2 orbital
        if(DIR.eq.4)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,0)-GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,0)
        endif
        ! d xy orbital
        if(DIR.eq.5)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,0)
        endif

    endif

    !!! orbitals over d shell are wrong ... this need to be fixed if one wishes to use that
    if(AO.eq.4)then

        if(DIR.eq.1)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,3,0,0)
        endif
        if(DIR.eq.2)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,3,0)
        endif
        if(DIR.eq.3)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,3)
        endif
        if(DIR.eq.4)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,2,1,0)
        endif
        if(DIR.eq.5)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,1)
        endif
        if(DIR.eq.6)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,2,0)
        endif
        if(DIR.eq.7)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,1)
        endif
        if(DIR.eq.8)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,2)
        endif
        if(DIR.eq.9)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,2)
        endif
        if(DIR.eq.10)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,1)
        endif

    endif

    if(AO.eq.5)then

        if(DIR.eq.1)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,4,0,0)
        endif
        if(DIR.eq.2)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,4,0)
        endif
        if(DIR.eq.3)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,4)
        endif
        if(DIR.eq.4)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,3,1,0)
        endif
        if(DIR.eq.5)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,3,0,1)
        endif
        if(DIR.eq.6)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,3,0)
        endif
        if(DIR.eq.7)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,3,1)
        endif
        if(DIR.eq.8)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,3)
        endif
        if(DIR.eq.9)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,3)
        endif
        if(DIR.eq.10)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,2,2,0)
        endif
        if(DIR.eq.11)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,2)
        endif
        if(DIR.eq.12)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,2)
        endif
        if(DIR.eq.13)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,2,1,1)
        endif
        if(DIR.eq.14)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,2,1)
        endif
        if(DIR.eq.15)then
            GTOKS = GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,2)
        endif

    endif

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function GTO(zeta,x,y,z,Rx,Ry,Rz,AO,DIR)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use basics
    implicit none
    complex(kind=comp_16) :: GTO
    real(kind=real_8), intent(in) :: zeta,x,y,z,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: AO,DIR
    real(kind=real_8) :: arg,coef

    GTO = dcmplx(0.d0,0.d0)
    arg = -1.d0*zeta*((x-Rx)**2.d0+(y-Ry)**2.d0+(z-Rz)**2.d0)
    coef = (2.d0*zeta/Pi)**(3.d0/4.d0)

    if(AO.eq.1)then
        GTO = coef*cdexp(dcmplx(arg,0.d0))
    endif

    if(AO.eq.2)then
        !p z orbital
        if(DIR.eq.1)then
            GTO = coef*dcmplx(2.d0*dsqrt(zeta),0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        ! p x orbital
        if(DIR.eq.2)then
            GTO = coef*dcmplx(2.d0*dsqrt(zeta),0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        ! p y orbital
        if(DIR.eq.3)then
            GTO = coef*dcmplx(2.d0*dsqrt(zeta),0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif


    endif

    if(AO.eq.3)then
        ! d z^2 orbital
        if(DIR.eq.1)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**2.d0)*(dcmplx(2.d0,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))&
                &-dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))-dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0)))
        endif
        ! d xz orbital
        if(DIR.eq.2)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**2.d0)*dcmplx(x-Rx,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        ! d yz orbital
        if(DIR.eq.3)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**2.d0)*dcmplx(y-Ry,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        ! d x^2-y^2 orbital
        if(DIR.eq.4)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**2.d0)*(dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))&
                &-dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0)))
        endif
        ! d xy orbital
        if(DIR.eq.5)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**2.d0)*dcmplx(x-Rx,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif

    endif

    !!! orbitals over d shell are wrong ... this need to be fixed if one wishes to use that
    if(AO.eq.4)then

        if(DIR.eq.1)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.2)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.3)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.4)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.5)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.6)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.7)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.8)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.9)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.10)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**3.d0)*dcmplx(z-Rz,0.d0)*dcmplx(y-Rz,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif

    endif

    if(AO.eq.5)then

        if(DIR.eq.1)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.2)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.3)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.4)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.5)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.6)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.7)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.8)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(x-Rx,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.9)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.10)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.11)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.12)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.13)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(x-Rx,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.14)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(y-Ry,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(x-Rx,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif
        if(DIR.eq.15)then
            GTO = coef*(dcmplx(2.d0*dsqrt(zeta),0.d0)**4.d0)*dcmplx(x-Rx,0.d0)*dcmplx(y-Ry,0.d0)*dcmplx(z-Rz,0.d0)*dcmplx(z-Rz,0.d0)*cdexp(dcmplx(arg,0.d0))
        endif

    endif

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Hnk(kk,n,zeta)
    ! A function similar to the hermite polynomials
    ! Very useful to write Hamiltonian matrix elements in a generalized form
    ! Stops at n = 5 but n > 5 can be calculated using H_n+1 = -(kk/2zeta)H_n + d/dk (H_n)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use basics
    implicit none
    real(kind=real_8) :: Hnk
    real(kind=real_8), intent(in) :: zeta,kk
    integer(kind=int_4), intent(in) :: n

    ! default value
    Hnk = 1.d0

    if(n.eq.1)then
        Hnk = -1.d0*kk/(2.d0*zeta)
    endif

    if(n.eq.2)then
        Hnk = (kk*kk-2.d0*zeta)/((2.d0*zeta)**2.d0)
    endif

    if(n.eq.3)then
        Hnk = (-1.d0*(kk**3.d0) + 6.d0*zeta*kk )/((2.d0*zeta)**3.d0)
    endif

    if(n.eq.4)then
        Hnk = (12.d0*zeta-12.d0*zeta*kk*kk+kk**4.d0)/((2.d0*zeta)**4.d0)
    endif

    if(n.eq.5)then
        Hnk = -1.d0*(60.d0*(zeta**2.d0)*kk -20.d0*zeta*(kk**3.d0) + kk**5.d0)/((2.d0*zeta)**5.d0)
    endif

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,l,m,n)
    ! Gives the value of the overlap between a plane wave k and
    ! a general Gaussian primitive: <k|g(zeta,R,l,m,n)>
    ! zeta is the width,
    ! R is the coordinate of the atomic center
    ! and l,m and n are the pseudo quantum numbers of the GTO primitive
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use basics
    implicit none
    complex(kind=comp_16) :: GTOK
    real(kind=real_8), intent(in) :: zeta,kx,ky,kz,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: l,m,n
    real(kind=real_8) :: arg,hermite
    complex(kind=comp_16) :: coef,phase

    if(zeta.ne.0.d0)then
        arg = -1.d0*((kx**2.d0) + (ky**2.d0) + (kz**2.d0))/(4.d0*zeta)

        hermite=Hnk(kx,l,zeta)*Hnk(ky,m,zeta)*Hnk(kz,n,zeta)*((2.d0*zeta*Pi)**(-3.d0/4.d0))

        coef=dcmplx(0.d0,2.d0*sqrt(zeta))**dfloat(l+m+n)

        phase=cdexp(dcmplx(0.d0,-1.d0*kx*Rx-1.d0*ky*Ry-1.d0*kz*Rz))

        GTOK = phase*coef*dcmplx(hermite,0.d0)*cdexp(dcmplx(arg,0.d0))
    else
        GTOK = dcmplx(0.d0,0.d0)
    endif

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,l,m,n)
    ! Gives the value of x-component of the transition dipole moment
    ! between a plane wave k and a general Gaussian primitive: <k|x|g(zeta,R,l,m,n)>
    ! zeta is the width,
    ! R is the coordinate of the atomic center
    ! and l,m and n are the pseudo quantum numbers of the GTO primitive
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use basics
    implicit none
    complex(kind=comp_16) :: MUXGTO
    real(kind=real_8), intent(in) :: zeta,kx,ky,kz,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: l,m,n
    integer(kind=int_4)::l_temp

    l_temp=l+1
    MUXGTO = Rx*GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,l,m,n) + (1.d0/ (2.d0*dsqrt(zeta) ) )*GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,l_temp,m,n)

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,l,m,n)
    ! Gives the value of y-component of the transition dipole moment
    ! between a plane wave k and a general Gaussian primitive: <k|y|g(zeta,R,l,m,n)>
    ! zeta is the width,
    ! R is the coordinate of the atomic center
    ! and l,m and n are the pseudo quantum numbers of the GTO primitive
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use basics
    implicit none
    complex(kind=comp_16) :: MUYGTO
    real(kind=real_8), intent(in) :: zeta,kx,ky,kz,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: l,m,n
    integer(kind=int_4)::m_temp

    m_temp=m+1
    MUYGTO = Ry*GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,l,m,n) + (1.d0/ (2.d0*dsqrt(zeta) ) )*GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,l,m_temp,n)

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,l,m,n)
    ! Gives the value of z-component of the transition dipole moment
    ! between a plane wave k and a general Gaussian primitive: <k|z|g(zeta,R,l,m,n)>
    ! zeta is the width,
    ! R is the coordinate of the atomic center
    ! and l,m and n are the pseudo quantum numbers of the GTO primitive
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use basics
    implicit none
    complex(kind=comp_16) :: MUZGTO
    real(kind=real_8), intent(in) :: zeta,kx,ky,kz,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: l,m,n
    integer(kind=int_4)::n_temp

    n_temp=n+1
    MUZGTO = Rz*GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,l,m,n) + (1.d0/ (2.d0*dsqrt(zeta) ) )*GTOK(zeta,kx,ky,kz,Rx,Ry,Rz,l,m,n_temp)

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function MUZGTOS(zeta,kx,ky,kz,Rx,Ry,Rz,AO,DIR)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use basics
    implicit none
    complex(kind=comp_16) :: MUZGTOS
    real(kind=real_8), intent(in) :: zeta,kx,ky,kz,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: AO,DIR


    MUZGTOS = dcmplx(0.d0,0.d0)

    if(AO.eq.1)then
        MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,0)
    endif

    if(AO.eq.2)then
        !p z orbital
        if(DIR.eq.1)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,1)
        endif
        ! p x orbital
        if(DIR.eq.2)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,0)
        endif
        ! p y orbital
        if(DIR.eq.3)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,0)
        endif

    endif

    if(AO.eq.3)then
        ! d z^2 orbital
        if(DIR.eq.1)then
            MUZGTOS = dcmplx(2.d0,0.d0)*MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,2)-MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,0)-MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,0)
        endif
        ! d xz orbital
        if(DIR.eq.2)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,1)
        endif
        ! d yz orbital
        if(DIR.eq.3)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,1)
        endif
        ! d x^2-y^2 orbital
        if(DIR.eq.4)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,0)-MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,0)
        endif
        ! d xy orbital
        if(DIR.eq.5)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,0)
        endif

    endif

    !!! orbitals over d shell are wrong ... this need to be fixed if one wishes to use that
    if(AO.eq.4)then

        if(DIR.eq.1)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,0,0)
        endif
        if(DIR.eq.2)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,3,0)
        endif
        if(DIR.eq.3)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,3)
        endif
        if(DIR.eq.4)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,1,0)
        endif
        if(DIR.eq.5)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,1)
        endif
        if(DIR.eq.6)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,2,0)
        endif
        if(DIR.eq.7)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,1)
        endif
        if(DIR.eq.8)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,2)
        endif
        if(DIR.eq.9)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,2)
        endif
        if(DIR.eq.10)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,1)
        endif

    endif

    if(AO.eq.5)then

        if(DIR.eq.1)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,4,0,0)
        endif
        if(DIR.eq.2)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,4,0)
        endif
        if(DIR.eq.3)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,4)
        endif
        if(DIR.eq.4)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,1,0)
        endif
        if(DIR.eq.5)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,0,1)
        endif
        if(DIR.eq.6)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,3,0)
        endif
        if(DIR.eq.7)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,3,1)
        endif
        if(DIR.eq.8)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,3)
        endif
        if(DIR.eq.9)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,3)
        endif
        if(DIR.eq.10)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,2,0)
        endif
        if(DIR.eq.11)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,2)
        endif
        if(DIR.eq.12)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,2)
        endif
        if(DIR.eq.13)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,1,1)
        endif
        if(DIR.eq.14)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,2,1)
        endif
        if(DIR.eq.15)then
            MUZGTOS = MUZGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,2)
        endif

    endif

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function MUXGTOS(zeta,kx,ky,kz,Rx,Ry,Rz,AO,DIR)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use basics
    implicit none
    complex(kind=comp_16) :: MUXGTOS
    real(kind=real_8), intent(in) :: zeta,kx,ky,kz,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: AO,DIR


    MUXGTOS = dcmplx(0.d0,0.d0)

    if(AO.eq.1)then
        MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,0)
    endif

    if(AO.eq.2)then
        !p z orbital
        if(DIR.eq.1)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,1)
        endif
        ! p x orbital
        if(DIR.eq.2)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,0)
        endif
        ! p y orbital
        if(DIR.eq.3)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,0)
        endif

    endif

    if(AO.eq.3)then
        ! d z^2 orbital
        if(DIR.eq.1)then
            MUXGTOS = dcmplx(2.d0,0.d0)*MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,2)-MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,0)-MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,0)
        endif
        ! d xz orbital
        if(DIR.eq.2)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,1)
        endif
        ! d yz orbital
        if(DIR.eq.3)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,1)
        endif
        ! d x^2-y^2 orbital
        if(DIR.eq.4)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,0)-MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,0)
        endif
        ! d xy orbital
        if(DIR.eq.5)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,0)
        endif

    endif

    !!! orbitals over d shell are wrong ... this need to be fixed if one wishes to use that
    if(AO.eq.4)then

        if(DIR.eq.1)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,0,0)
        endif
        if(DIR.eq.2)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,3,0)
        endif
        if(DIR.eq.3)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,3)
        endif
        if(DIR.eq.4)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,1,0)
        endif
        if(DIR.eq.5)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,1)
        endif
        if(DIR.eq.6)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,2,0)
        endif
        if(DIR.eq.7)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,1)
        endif
        if(DIR.eq.8)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,2)
        endif
        if(DIR.eq.9)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,2)
        endif
        if(DIR.eq.10)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,1)
        endif

    endif

    if(AO.eq.5)then

        if(DIR.eq.1)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,4,0,0)
        endif
        if(DIR.eq.2)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,4,0)
        endif
        if(DIR.eq.3)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,4)
        endif
        if(DIR.eq.4)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,1,0)
        endif
        if(DIR.eq.5)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,0,1)
        endif
        if(DIR.eq.6)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,3,0)
        endif
        if(DIR.eq.7)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,3,1)
        endif
        if(DIR.eq.8)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,3)
        endif
        if(DIR.eq.9)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,3)
        endif
        if(DIR.eq.10)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,2,0)
        endif
        if(DIR.eq.11)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,2)
        endif
        if(DIR.eq.12)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,2)
        endif
        if(DIR.eq.13)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,1,1)
        endif
        if(DIR.eq.14)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,2,1)
        endif
        if(DIR.eq.15)then
            MUXGTOS = MUXGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,2)
        endif

    endif

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function MUYGTOS(zeta,kx,ky,kz,Rx,Ry,Rz,AO,DIR)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use basics
    implicit none
    complex(kind=comp_16) :: MUYGTOS
    real(kind=real_8), intent(in) :: zeta,kx,ky,kz,Rx,Ry,Rz
    integer(kind=int_4), intent(in) :: AO,DIR


    MUYGTOS = dcmplx(0.d0,0.d0)

    if(AO.eq.1)then
        MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,0)
    endif

    if(AO.eq.2)then
        !p z orbital
        if(DIR.eq.1)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,1)
        endif
        ! p x orbital
        if(DIR.eq.2)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,0)
        endif
        ! p y orbital
        if(DIR.eq.3)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,0)
        endif

    endif

    if(AO.eq.3)then
        ! d z^2 orbital
        if(DIR.eq.1)then
            MUYGTOS = dcmplx(2.d0,0.d0)*MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,2)-MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,0)-MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,0)
        endif
        ! d xz orbital
        if(DIR.eq.2)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,1)
        endif
        ! d yz orbital
        if(DIR.eq.3)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,1)
        endif
        ! d x^2-y^2 orbital
        if(DIR.eq.4)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,0)-MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,0)
        endif
        ! d xy orbital
        if(DIR.eq.5)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,0)
        endif

    endif

    !!! orbitals over d shell are wrong ... this need to be fixed if one wishes to use that
    if(AO.eq.4)then

        if(DIR.eq.1)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,0,0)
        endif
        if(DIR.eq.2)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,3,0)
        endif
        if(DIR.eq.3)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,3)
        endif
        if(DIR.eq.4)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,1,0)
        endif
        if(DIR.eq.5)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,1)
        endif
        if(DIR.eq.6)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,2,0)
        endif
        if(DIR.eq.7)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,1)
        endif
        if(DIR.eq.8)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,2)
        endif
        if(DIR.eq.9)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,2)
        endif
        if(DIR.eq.10)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,1)
        endif

    endif

    if(AO.eq.5)then

        if(DIR.eq.1)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,4,0,0)
        endif
        if(DIR.eq.2)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,4,0)
        endif
        if(DIR.eq.3)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,0,4)
        endif
        if(DIR.eq.4)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,1,0)
        endif
        if(DIR.eq.5)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,3,0,1)
        endif
        if(DIR.eq.6)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,3,0)
        endif
        if(DIR.eq.7)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,3,1)
        endif
        if(DIR.eq.8)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,0,3)
        endif
        if(DIR.eq.9)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,1,3)
        endif
        if(DIR.eq.10)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,2,0)
        endif
        if(DIR.eq.11)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,0,2)
        endif
        if(DIR.eq.12)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,0,2,2)
        endif
        if(DIR.eq.13)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,2,1,1)
        endif
        if(DIR.eq.14)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,2,1)
        endif
        if(DIR.eq.15)then
            MUYGTOS = MUYGTO(zeta,kx,ky,kz,Rx,Ry,Rz,1,1,2)
        endif

    endif

end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine buildmuzk(muzk,lcao,eta,zet,coord,lmnp1,CMap,NS,NCONS,norb,nq,np,maxprim,xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,nx,ny,nz)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    integer(kind=int_4), intent(in) :: NS,NCONS,norb,nq,np,maxprim,nx,ny,nz
    integer(kind=int_4), dimension(NCONS), intent(in) :: LMNP1, CMap
    real(kind=real_8), dimension(NS,3), intent(in) :: coord
    real(kind=real_8), dimension(norb,norb), intent(in) :: lcao
    real(kind=real_8), dimension(maxprim,NCONS), intent(in) :: eta,zet
    real(kind=real_8), intent(in) :: xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin
    complex(kind=comp_16), dimension(np,nq), intent(out) :: muzk

    real(kind=real_8) :: kx,ky,kz,Rx,Ry,Rz,zeta
    integer(kind=int_4) :: i,j,k,l,r,s,t,u,CN,DIR,nw
    complex(kind=comp_16), dimension(norb) :: muk
    complex(kind=comp_16) :: tempo

    do i = 1, nx
        do j = 1, ny
            do k = 1, nz

                u=0
                do s=1,NCONS
                    do t=1,2*(LMNP1(s)-1)+1
                        u=u+1

                        tempo = dcmplx(0.d0,0.d0)
                        !$OMP PARALLEL DO REDUCTION(+:tempo)
                        do l=1,maxprim
                            tempo=tempo+dcmplx(eta(l,s),0.d0)**MUZGTOS(zet(l,s),kxmin+(i-1)*dkx,kymin+(j-1)*dky,kzmin+(k-1)*dkz,coord(CMap(s),1),coord(CMap(s),2),coord(CMap(s),3),LMNP1(s),t)
                        enddo
                        !$OMP END PARALLEL DO
                        muk(u) = tempo

                    enddo
                enddo



                do r=1,nq
                    tempo = dcmplx(0.d0,0.d0)
                    !$OMP PARALLEL DO REDUCTION(+:tempo)
                    do u=1,norb
                        tempo=tempo+muk(u)*dcmplx(lcao(u,r),0.d0)
                    enddo
                    !$OMP END PARALLEL DO
                    muzk((k-1)*nx*ny+(j-1)*nx+i,r) = tempo
                enddo



            enddo
        enddo
    enddo

end subroutine buildmuzk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine buildmuxk(muxk,lcao,eta,zet,coord,lmnp1,CMap,NS,NCONS,norb,nq,np,maxprim,xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,nx,ny,nz)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    integer(kind=int_4), intent(in) :: NS,NCONS,norb,nq,np,maxprim,nx,ny,nz
    integer(kind=int_4), dimension(NCONS), intent(in) :: LMNP1, CMap
    real(kind=real_8), dimension(NS,3), intent(in) :: coord
    real(kind=real_8), dimension(norb,norb), intent(in) :: lcao
    real(kind=real_8), dimension(maxprim,NCONS), intent(in) :: eta,zet
    real(kind=real_8), intent(in) :: xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin
    complex(kind=comp_16), dimension(np,nq), intent(out) :: muxk

    real(kind=real_8) :: kx,ky,kz,Rx,Ry,Rz,zeta
    integer(kind=int_4) :: i,j,k,l,r,s,t,u,CN,DIR,nw
    complex(kind=comp_16), dimension(norb) :: muk
    complex(kind=comp_16) :: tempo

    do i = 1, nx
        do j = 1, ny
            do k = 1, nz

                u=0
                do s=1,NCONS
                    do t=1,2*(LMNP1(s)-1)+1
                        u=u+1

                        tempo = dcmplx(0.d0,0.d0)
                        !$OMP PARALLEL DO REDUCTION(+:tempo)
                        do l=1,maxprim
                            tempo=tempo+dcmplx(eta(l,s),0.d0)**MUXGTOS(zet(l,s),kxmin+(i-1)*dkx,kymin+(j-1)*dky,kzmin+(k-1)*dkz,coord(CMap(s),1),coord(CMap(s),2),coord(CMap(s),3),LMNP1(s),t)
                        enddo
                        !$OMP END PARALLEL DO
                        muk(u) = tempo

                    enddo
                enddo

                do r=1,nq
                    tempo = dcmplx(0.d0,0.d0)
                    !$OMP PARALLEL DO REDUCTION(+:tempo)
                    do u=1,norb
                        tempo=tempo+muk(u)*dcmplx(lcao(u,r),0.d0)
                    enddo
                    !$OMP END PARALLEL DO
                    muxk((k-1)*nx*ny+(j-1)*nx+i,r) = tempo
                enddo

            enddo
        enddo
    enddo

end subroutine buildmuxk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine buildmuyk(muyk,lcao,eta,zet,coord,lmnp1,CMap,NS,NCONS,norb,nq,np,maxprim,xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin,nx,ny,nz)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    integer(kind=int_4), intent(in) :: NS,NCONS,norb,nq,np,maxprim,nx,ny,nz
    integer(kind=int_4), dimension(NCONS), intent(in) :: LMNP1, CMap
    real(kind=real_8), dimension(NS,3), intent(in) :: coord
    real(kind=real_8), dimension(norb,norb), intent(in) :: lcao
    real(kind=real_8), dimension(maxprim,NCONS), intent(in) :: eta,zet
    real(kind=real_8), intent(in) :: xmin,ymin,zmin,dx,dy,dz,dkx,dky,dkz,kxmin,kymin,kzmin
    complex(kind=comp_16), dimension(np,nq), intent(out) :: muyk

    real(kind=real_8) :: kx,ky,kz,Rx,Ry,Rz,zeta
    integer(kind=int_4) :: i,j,k,l,r,s,t,u,CN,DIR,nw
    complex(kind=comp_16), dimension(norb) :: muk
    complex(kind=comp_16) :: tempo

    do i = 1, nx
        do j = 1, ny
            do k = 1, nz

                u=0
                do s=1,NCONS
                    do t=1,2*(LMNP1(s)-1)+1
                        u=u+1

                        tempo = dcmplx(0.d0,0.d0)
                        !$OMP PARALLEL DO REDUCTION(+:tempo)
                        do l=1,maxprim
                            tempo=tempo+dcmplx(eta(l,s),0.d0)**MUYGTOS(zet(l,s),kxmin+(i-1)*dkx,kymin+(j-1)*dky,kzmin+(k-1)*dkz,coord(CMap(s),1),coord(CMap(s),2),coord(CMap(s),3),LMNP1(s),t)
                        enddo
                        !$OMP END PARALLEL DO
                        muk(u) = tempo

                    enddo
                enddo


                do r=1,nq
                    tempo = dcmplx(0.d0,0.d0)
                    !$OMP PARALLEL DO REDUCTION(+:tempo)
                    do u=1,norb
                        tempo=tempo+muk(u)*dcmplx(lcao(u,r),0.d0)
                    enddo
                    !$OMP END PARALLEL DO
                    muyk((k-1)*nx*ny+(j-1)*nx+i,r) = tempo
                enddo

            enddo
        enddo
    enddo

end subroutine buildmuyk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************************************************
!!****************************************************************************************
! subroutine param_debut(basis_file,nq,inR,Ne,Norb_cat,Norb_sym,FC,restriction,S,job,nx,xmin,xmax,dx,dkx,kxmin,E0,omega,delta,theta,phy,pulsed,pulsetime,nper,nt)
!!****************************************************************************************
!!****************************************************************************************

!use basics
!implicit none
!! Parameters that are read
! Character(len=32), Intent(out) :: basis_file,job
! Integer(kind=int_4), Intent(out) :: nq, inR
! Integer(kind=int_4), Intent(out) :: Ne,Norb_cat,Norb_sym,FC
! Integer(kind=int_4), Intent(out) :: nx,nper,nt
! Real(kind=real_8), Intent(out) :: dx,dkx,xmin,xmax,kxmin,omega,delta,theta,phy,pulsetime
! Logical :: restriction,pulsed
! Real(kind=real_8), Intent(out) :: S,E0


!Namelist/grille/basis_file,nq,inR,job,nx,xmin
!Namelist/drt/Ne,Norb_cat,Norb_sym,FC,S,restriction
!Namelist/champ/E0,omega,delta,theta,phy,pulsed,pulsetime,nper,nt

!open(1,file='input',status='old')
!read(1,nml=grille)
!read(1,nml=drt)
!read(1,nml=champ)
! close(1)

!basis_file=trim(basis_file)

!dx=2.d0*dabs(xmin)/dfloat(nx)
!xmax=xmin+dfloat(nx-1)*dx

!!nperf=(pulsetime)*w/(2.d0*pi)
!!nperi=idint((pulsetime)*w/(2.d0*pi))
!!nper=Max((2*nperi)+nint(nperf-dfloat(nperi)),2)
!!pulsetime=(dfloat(nper)*2.d0*pi)/w

!!tmax=dfloat(nper*2)*pi/w
!! cycleopt=(2.d0*Pi)/w
!!tmax=6.d0*cycleopt!(5.d0/4.d0)*cycleopt
!!nt=((tmax-tmin)/pdt)+1
!!Tfwhm=Pi/w
!!tau=(4.d0*Pi)/w


!write(*,*)' basis_file ->',basis_file
!write(*,*)' job ->',job
!write(*,'(A17)') 'Spatial grid size'
!write(*,'(A6,F7.2,A3,F7.2,A4,E10.3,A4,I3)')'axe X ',xmin,' a ',xmax,' dx=',dx,' nx=',nx
!!write(*,*)'Options,','dkOpt=',dkOpt

!dkx=(2.d0*pi)/(dfloat(nx)*dx)
!kxmin=-pi/dx


!write(*,'(A22)')'Grid in momentum space'
!write(*,'(A6,F8.2,A5,F8.2,A4,E10.3)')'axe K ',kxmin,' a-1 ',kxmin+dfloat(nx-1)*dkx,' dk=',dkx
!write(*,*)'S=',S,'Ne=',Ne

!!write(*,'(A17)')'Grille temporelle'
!!write(*,*)'tmin',tmin,'pdt',pdt,'w',w,'E0',E0
!end subroutine param_debut
!!****************************************************************************************
!!****************************************************************************************
!subroutine param_debut(basis_file,norb_lect,nq,inR,orb_min,orb_max)
!!****************************************************************************************
!!****************************************************************************************

!use basics
!implicit none
!! Parameters that are read
! Character(len=32), Intent(out) :: basis_file
! Integer(kind=int_4), Intent(out) :: norb_lect, nq, inR,orb_min,orb_max


!Namelist/grille/basis_file,norb_lect,nq,inR,orb_min,orb_max

!open(1,file='input',status='old')
!read(1,nml=grille)
! close(1)

!basis_file=trim(basis_file)
!end subroutine param_debut


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! for MEDYS*

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!****************************************************************************************
    !****************************************************************************************
subroutine PrimCount(NS,NF,ICONU, NPrimCount,LMNP1) !PrimCount(NS,NF,NC,ICONU, NPrimCount)
    !****************************************************************************************
!****************************************************************************************

Integer(kind=int_4), Dimension(:),Intent(in)	::ICONU
Integer(kind=int_4), Dimension(:),Intent(in)	::NF !, NC
Integer(kind=int_4), Dimension(:),Intent(in)	::LMNP1
Integer(kind=int_4), Intent(in)					:: NS 
Integer(kind=int_4), Intent(out)	::NPrimCount
Integer(kind=int_8) ::  counter1, counter2, counter3, I,J,IS,m !, K

counter1=0
counter2=0
counter3=0
NPrimCount=0
Do IS=1,NS
    
  Do J=1,NF(IS)           !NF(IS) =  Number of AO sets for each of the symmetry-distinct atoms 
    counter1=counter1+1  ! counter1 compte les classes  d'OA (AO sets)
 ! ajouter:
!     Do K = 1, NC(IS)
!
!
!  

    select case (LMNP1(counter1))  !LMNP1 = 1,2,3 correspond à s,p,d
!
      case(1)     ! J est une orbitale s centrée sur IS
 
      counter2=counter2+1  ! (le numéro de cette OA est counter2)

      Do I=1, ICONU(counter1)
        counter3=counter3+1  ! (compte les primitives) 
      enddo

!
      case(2)     ! J est une orbitale p centree sur IS
      
      do m=1,3
        counter2=counter2+1  ! (le numéro de cette OA p_m est counter2)
        Do I=1, ICONU(counter1)
          counter3=counter3+1  ! (compte les primitives) 
        enddo
      enddo
!
      case(3)     ! J est une orbitale d centree sur IS
      
      do m=1,6
        counter2=counter2+1  ! (le numéro de cette OA d_m est counter2)
        Do I=1, ICONU(counter1)
          counter3=counter3+1  ! (compte les primitives) 
        enddo
      enddo
!
!     enddo
!
 !
    end select 

  enddo
enddo

NPrimCount=counter3

end subroutine PrimCount

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!****************************************************************************************
    !****************************************************************************************
subroutine Join_Prim_Info(NS,ETA,ZET,NCONS,ICONU,NF,LMNP1,coord,cEta,cZeta,lmn_vec,prim_center,lcao,lcCG,lp,ld,nMO)
    !****************************************************************************************
!****************************************************************************************

Integer(kind=int_4), Intent(in)			::NS,nMO
Integer(kind=int_4), Intent(in)			::NCONS
Integer(kind=int_4), Dimension(:),Intent(in)	::ICONU
Integer(kind=int_4), Dimension(:),Intent(in)	::NF
Integer(kind=int_4), Dimension(:),Intent(in)	::LMNP1
!Integer(kind=int_8) Intent(in)	::totPrimCount 
Integer(kind=int_8) ::  counter1, counter2, counter3, I,J,IS,m,K
Real(kind=real_8), Dimension(:,:),Intent(inout)	::ZET,ETA,coord
Real(kind=real_8), Dimension(:,:),Intent(in)	::lcao
Real(kind=real_8),dimension(:,:),Intent(inout)	:: cEta, lcCG,prim_center
Real(kind=real_8),dimension(:),Intent(inout)	::cZeta
!Integer(kind=int_4),allocatable,dimension(:,:)::prim_lmn
Integer(kind=int_4),dimension(:,:),Intent(inout)	:: lmn_vec
Integer(kind=int_4),dimension(3)::lp
Integer(kind=int_4),dimension(6)::ld

 counter1=0
 counter2=0
 counter3=0
lmn_vec=0
 Do IS=1,NS
    
   Do J=1,NF(IS)
      counter1=counter1+1  ! counter1 compte les classes  d'OA 
      select case (LMNP1(counter1))  !LMNP1 = 1,2,3 correspond à s,p,d
!
case(1)     ! J est une orbitale s centree sur IS
 
      counter2=counter2+1  ! (le numéro de cette OA est counter2)

       Do I=1, ICONU(counter1)
       counter3=counter3+1  ! (compte les primitives)
       cEta(counter2,counter3)=ETA(I,counter1)
       cZeta(counter3)=ZET(I,counter1) !todo: verifier ZET ou ZETA
       prim_center(counter3,:)=coord(counter1,:)
       enddo

!
case(2)     ! J est une orbitale p centree sur IS
      
      do m=1,3
         counter2=counter2+1  ! (le numéro de cette OA p_m est counter2)
       Do I=1, ICONU(counter1)
       counter3=counter3+1  ! (compte les primitives) 
       cEta(counter2,counter3)=ETA(I,counter1)
       cZeta(counter3)=ZET(I,counter1) !todo: verifier ZET ou ZETA
       prim_center(counter3,:)=coord(counter1,:)
          
          select case (lp(m))
          case (1)
           lmn_vec(counter3,1)=1
          case (2)
           lmn_vec(counter3,2)=1
          case (3)
           lmn_vec(counter3,3)=1
          end select 
          
       enddo
      enddo
!
case(3)     ! J est une orbitale d centree sur IS
      
      do m=1,6
         counter2=counter2+1  ! (le numéro de cette OA d_m est counter2)
       Do I=1, ICONU(counter1)
       counter3=counter3+1  ! (compte les primitives) 
       cEta(counter2,counter3)=ETA(I,counter1)
       cZeta(counter3)=ZET(I,counter1) !todo: verifier ZET ou ZETA
       prim_center(counter3,:)=coord(counter1,:)
          
          select case (ld(m))
          case (1)
           lmn_vec(counter3,1)=2
          case (2)
           lmn_vec(counter3,2)=2
          case (3)
           lmn_vec(counter3,3)=2
         case (4)
           lmn_vec(counter3,1)=1
           lmn_vec(counter3,2)=1
          case (5)
           lmn_vec(counter3,1)=1
           lmn_vec(counter3,3)=1
          case (6)
           lmn_vec(counter3,2)=1
           lmn_vec(counter3,3)=1
          end select 
          
       enddo
      enddo

     end select 

  enddo
 enddo
call gemm(lcao,cEta,lcCG)


end subroutine Join_Prim_Info

END MODULE read_columbus
