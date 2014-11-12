
module nakajima_zwanzig_shared

    use std_types
    use helpers
    use numer_fft
    use std_lapack
    use numer_matrix

    implicit none

    public::N1_from_type
    public::N2_from_type
    public::operator_to_exc
    public::operator_from_exc
    public::superops_to_exc ! transforms the superoperator into exc picture directly

    interface superops_to_exc
        module procedure superops_to_exc_2indexed
        module procedure superops_to_exc_4indexed
    end interface

    interface superops_from_exc
        module procedure superops_from_exc_2indexed
        module procedure superops_from_exc_4indexed
    end interface

    interface operator_to_exc
        module procedure operator_to_exc2
        module procedure operator_to_exc1
    end interface

    interface operator_from_exc
        module procedure operator_from_exc2
        module procedure operator_from_exc1
    end interface

    private::operator_from_exc2
    private::operator_from_exc1
    private::operator_to_exc2
    private::operator_to_exc1

    private::SUPERINDEX_FROM_K_L
    private::K_FROM_SUPERINDEX
    private::L_FROM_SUPERINDEX

    real(dp), dimension(:), pointer, private :: eng,en2

    integer(i4b), private :: N_of_sites
    real(dp), dimension(:,:), allocatable :: SS, S1, SS_2, S1_2, He


    contains

    subroutine init_nakajima_zwanzig_shared (Ham)
      real(dp), dimension(:,:), intent(in) :: Ham

      N_of_sites = size(Ham,1)

      allocate(SS(size(Ham,1), size(Ham,1)))
      allocate(S1(size(Ham,1), size(Ham,1)))
      allocate(He(size(Ham,1), size(Ham,1)))

      call spec(Ham,SS,He)
      call inv(SS,S1)

      write(*,*) 'SS:', SS
      write(*,*) 'S1:', S1
    end subroutine init_nakajima_zwanzig_shared

    pure function SUPERINDEX_FROM_K_L(k, l, lmax) result(res)
      integer(i4b) :: res
      integer(i4b), intent(in)  :: k, l, lmax

      res = ((l) + (lmax)*((k)-1))
    end function SUPERINDEX_FROM_K_L

    pure function K_FROM_SUPERINDEX(superindex, lmax) result(res)
      integer(i4b) :: res
      integer(i4b), intent(in)  :: superindex, lmax

      res = (((superindex) - 1) / (lmax) + 1)
    end function

    pure function L_FROM_SUPERINDEX(superindex, lmax) result(res)
      integer(i4b) :: res
      integer(i4b), intent(in)  :: superindex, lmax

      res = (mod(((superindex)-1) , (lmax)) + 1)
    end function

    pure function N1_from_type(type) result(NN1)
    ! gives left operator size according to type
         character, intent(in)           :: type
         integer(i4b)                    :: NN1

         if(type == 'g') then
             NN1 = 1
         else if(type == 'O') then
             NN1 = N_of_sites
         else if(type == 'E') then
             NN1 = N_of_sites
         else if(type == '2') then
             NN1 = N_of_sites*(N_of_sites-1)/2
         else if(type == 'F') then
             NN1 = N_of_sites*(N_of_sites-1)/2
         else
             NN1 = -1
         end if

     end function

    pure function N2_from_type(type) result(NN2)
    ! gives right operator size according to type
         character, intent(in) :: type
         integer        :: NN2

         if(type == 'g') then
             NN2 = 1
         else if(type == 'O') then
             NN2 = 1
         else if(type == 'E') then
             NN2 = N_of_sites
         else if(type == '2') then
             NN2 = N_of_sites
         else if(type == 'F') then
             NN2 = N_of_sites*(N_of_sites-1)/2
         else
             NN2 = -1
         end if

     end function

     subroutine operator_to_exc2(fromto, type)
        complex(dpc), dimension(:,:), intent(inout)                 :: fromto
        character                                                    :: type

        integer                                                        :: NN1, NN2

        NN1 = N1_from_type(type)
        NN2 = N2_from_type(type)

        if(.not.(type == 'O' .or. type == 'E' .or. type == '2' .or. type == 'F')) then
             write(*,*) 'wrong type in operator_to_exc() : '
             stop
         else if(type == 'O') then
             fromto = matmul(transpose(S1),fromto)
         else if(type == 'E') then
             fromto = matmul(matmul(S1,fromto),SS)
         else if(type == '2') then
             fromto = matmul(matmul(S1_2,fromto),SS)
         else if(type == 'F') then
             fromto = matmul(matmul(S1_2,fromto),SS_2)
         end if

     end subroutine operator_to_exc2

     subroutine operator_to_exc1(fromto, type)
        complex(dpc), dimension(:), intent(inout)                   :: fromto
        character                                                   :: type

        integer(i4b)                                                :: a
        complex(dpc), dimension(size(fromto),1)                     :: rho

        do a=1, size(fromto)
            rho(a,1) = fromto(a)
        end do

        call operator_to_exc2(rho, type)

        do a=1, size(fromto)
            fromto(a) = rho(a,1)
        end do
    end subroutine operator_to_exc1

    subroutine operator_from_exc1(fromto, type)
        complex(dpc), dimension(:), intent(inout)                   :: fromto
        character                                                   :: type

        integer(i4b)                                                :: a
        complex(dpc), dimension(size(fromto),1)                     :: rho

        do a=1, size(fromto)
            rho(a,1) = fromto(a)
        end do

        call operator_from_exc2(rho, type)

        do a=1, size(fromto)
            fromto(a) = rho(a,1)
        end do
    end subroutine operator_from_exc1

     subroutine operator_from_exc2(fromto, type)
        complex(dpc), dimension(:,:), intent(inout)                 :: fromto
        character                                                    :: type

        integer                                                        :: NN1, NN2

        NN1 = N1_from_type(type)
        NN2 = N2_from_type(type)

        if(.not.(type == 'O' .or. type == 'E' .or. type == '2' .or. type == 'F')) then
             write(*,*) 'wrong type in operator_to_exc() : '
             stop
         else if(type == 'O') then
             fromto = matmul(transpose(SS),fromto)
         else if(type == 'E') then
             fromto = matmul(matmul(SS,fromto),S1)
         else if(type == '2') then
             fromto = matmul(matmul(SS_2,fromto),S1)
         else if(type == 'F') then
             fromto = matmul(matmul(SS_2,fromto),S1_2)
         end if

     end subroutine operator_from_exc2

    subroutine superops_to_exc_2indexed(fromto, type)
        complex(dpc), dimension(:,:), intent(inout)        :: fromto
        character                                            :: type

        complex(dpc), dimension(size(fromto,1),size(fromto,1)) :: Xi, Xi1
        integer                                                 :: a,b,c,d
        integer                                                :: NN1, NN2
        complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) :: rho

        NN1 = N1_from_type(type)
        NN2 = N2_from_type(type)
        Xi  = 0.0_dp
        Xi1 = 0.0_dp

        do d=1, NN2
        do c=1, NN1
                rho = 0.0_dp
                rho(c,d) = 1.0_dp

                call operator_to_exc(rho,type)

                do b=1, NN2
                do a=1, NN1
                    Xi(SUPERINDEX_FROM_K_L(a,b,NN2),SUPERINDEX_FROM_K_L(c,d,NN2)) = rho(a,b)
                end do
                end do
        end do
        end do

        call inv(Xi,Xi1)

        fromto = matmul(matmul(Xi1,fromto),Xi)

    end subroutine superops_to_exc_2indexed

    subroutine superops_to_exc_4indexed(fromto, type)
        complex(dpc), dimension(:,:,:,:), intent(inout)    :: fromto
        character, intent(in)                           :: type

        complex(dpc), dimension(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type)) :: XX
        integer(i4b) :: i,j,k,l,  first_N, second_N

        XX = 0.0_dp

        first_N  = N1_from_type(type)
        second_N = N2_from_type(type)

        if(.not.(size(fromto,1) == first_N .or. size(fromto,2) == second_N .or. size(fromto,3) == first_N .or. size(fromto,4) == second_N)) then
             write(*,*) 'wrong dimension in superops_to_exc_4indexed() : '
             stop
         end if

        do i=1, first_N
        do j=1, second_N
          call operator_to_exc(fromto(i,j,:,:),type)
        end do
        end do

        do i=1, first_N
        do j=1, second_N
          call operator_to_exc(fromto(:,:,i,j),type)
        end do
        end do

!        do k=1, first_N
!        do l=1, second_N
!
!        do i=1, first_N
!        do j=1, second_N
!
!        XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
!            fromto(i,j,k,l)
!
!        end do
!        end do
!
!        end do
!        end do
!
!        call superops_to_exc(XX, type)
!
!        do k=1, first_N
!        do l=1, second_N
!
!        do i=1, first_N
!        do j=1, second_N
!
!        fromto(i,j,k,l) = XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N))
!
!        end do
!        end do
!
!        end do
!        end do

    end subroutine superops_to_exc_4indexed

    subroutine superops_from_exc_2indexed(fromto, type)
        complex(dpc), dimension(:,:), intent(inout)        :: fromto
        character                                            :: type

        complex(dpc), dimension(size(fromto,1),size(fromto,1)) :: Xi, Xi1
        integer                                                 :: a,b,c,d
        integer                                                :: NN1, NN2
        complex(dpc), dimension(N1_from_type(type),N2_from_type(type)) :: rho

        NN1 = N1_from_type(type)
        NN2 = N2_from_type(type)
        Xi  = 0.0_dp
        Xi1 = 0.0_dp

        do d=1, NN2
        do c=1, NN1
                rho = 0.0_dp
                rho(c,d) = 1.0_dp

                call operator_to_exc(rho,type)

                do b=1, NN2
                do a=1, NN1
                    Xi(SUPERINDEX_FROM_K_L(a,b,NN2),SUPERINDEX_FROM_K_L(c,d,NN2)) = rho(a,b)
                end do
                end do
        end do
        end do

        call inv(Xi,Xi1)

        fromto = matmul(matmul(Xi,fromto),Xi1)

    end subroutine superops_from_exc_2indexed

    subroutine superops_from_exc_4indexed(fromto, type)
        complex(dpc), dimension(:,:,:,:), intent(inout)    :: fromto
        character, intent(in)                            :: type

        complex(dpc), dimension(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type)) :: XX
        integer(i4b) :: i,j,k,l,  first_N, second_N

        XX = 0.0_dp

        first_N  = N1_from_type(type)
        second_N = N2_from_type(type)

        if(.not.(size(fromto,1) == first_N .or. size(fromto,2) == second_N .or. size(fromto,3) == first_N .or. size(fromto,4) == second_N)) then
             write(*,*) 'wrong dimension in superops_to_exc_4indexed() : '
             stop
         end if

        do i=1, first_N
        do j=1, second_N
          call operator_from_exc(fromto(i,j,:,:),type)
        end do
        end do

        do i=1, first_N
        do j=1, second_N
          call operator_from_exc(fromto(:,:,i,j),type)
        end do
        end do

!        do k=1, first_N
!        do l=1, second_N
!
!        do i=1, first_N
!        do j=1, second_N
!
!        XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
!            fromto(i,j,k,l)
!
!        end do
!        end do
!
!        end do
!        end do
!
!        call superops_from_exc(XX, type)
!
!        do k=1, first_N
!        do l=1, second_N
!
!        do i=1, first_N
!        do j=1, second_N
!
!        fromto(i,j,k,l) = XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N))
!
!        end do
!        end do
!
!        end do
!        end do

    end subroutine superops_from_exc_4indexed

    subroutine superops_4indexed_to_2indexed(from, to, type)
        complex(dpc), dimension(:,:,:,:), intent(in)        :: from
        character, intent(in)                                :: type

        complex(dpc), dimension(:,:), intent(out) :: to
        integer(i4b) :: i,j,k,l,  first_N, second_N

        if(.not.(size(to,1) == size(to,2) .and. (N1_from_type(type)*N2_from_type(type) == size(to,1)) .and. &
            (size(from,1) == size(from,3)) .and. (size(from,2) == size(from,4)) .and. (size(from,1) == N1_from_type(type)) &
            .and. (size(from,2) == N2_from_type(type)) )) then

            write(*,*) "superops_4indexed_to_2indexed - size error"
            stop
        end if

        to = 0.0_dp

        first_N  = N1_from_type(type)
        second_N = N2_from_type(type)

        do k=1, first_N
        do l=1, second_N

        do i=1, first_N
        do j=1, second_N

        to(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
            from(i,j,k,l)

        end do
        end do

        end do
        end do

    end subroutine superops_4indexed_to_2indexed


    subroutine superops_2indexed_to_4indexed(from, to, type)
        complex(dpc), dimension(:,:), intent(in)        :: from
        character, intent(in)                            :: type

        complex(dpc), dimension(:,:,:,:), intent(out) :: to
        integer(i4b) :: i,j,k,l,  first_N, second_N

        if(.not.(size(from,1) == size(from,2) .and. (N1_from_type(type)*N2_from_type(type) == size(from,1)) .and. &
            (size(to,1) == size(to,3)) .and. (size(to,2) == size(to,4)) .and. (size(to,1) == N1_from_type(type)) &
            .and. (size(to,2) == N2_from_type(type)) )) then

            write(*,*) "superops_2indexed_to_4indexed - size error"
            stop
        end if

        to = 0.0_dp

        first_N  = N1_from_type(type)
        second_N = N2_from_type(type)


        do k=1, first_N
        do l=1, second_N

        do i=1, first_N
        do j=1, second_N

        to(i,j,k,l) = from(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N))

        end do
        end do

        end do
        end do

    end subroutine superops_2indexed_to_4indexed

!    subroutine superops_4indexed_to_2indexed_updiag(from, to, type)
!        complex(dpc), dimension(:,:,:,:), intent(in)        :: from
!        character, intent(in)                               :: type
!
!        complex(dpc), dimension(:,:), intent(out) :: to
!        integer(i4b) :: i,j,k,l,  first_N, second_N,  supI, supJ
!
!        if(.not.(size(to,1) == size(to,2) .and. (N1_from_type(type) == N2_from_type(type)) .and. (N1_from_type(type)*(N1_from_type(type)+1)/2 == size(to,1)) .and. &
!            (size(from,1) == size(from,3)) .and. (size(from,2) == size(from,4)) .and. (size(from,1) == N1_from_type(type)) &
!            .and. (size(from,2) == N2_from_type(type)) )) then
!
!            write(*,*) "superops_4indexed_to_2indexed_updiag - size error", size(from,1), size(from,2)
!            stop
!        end if
!
!        to = 0.0_dp
!
!        first_N  = N1_from_type(type)
!        second_N = N2_from_type(type)
!
!        supJ = 0
!
!        do l=1, first_N
!        do k=1, l
!
!        supJ = supJ + 1
!        supI = 0
!
!        do j=1, first_N
!        do i=1, j
!
!        supI = supI + 1
!
!        to(supI,supJ) = from(i,j,k,l)
!        !write(*,*) supI,supJ, 'XX', i,j,k,l
!
!        end do
!        end do
!
!        end do
!        end do
!
!    end subroutine superops_4indexed_to_2indexed_updiag


!    subroutine superops_2indexed_to_4indexed_updiag(from, to, type)
!        complex(dpc), dimension(:,:), intent(in)        :: from
!        character, intent(in)                           :: type
!
!        complex(dpc), dimension(:,:,:,:), intent(out) :: to
!        integer(i4b) :: i,j,k,l,  first_N, second_N,  supI, supJ
!
!        if(.not.(size(from,1) == size(from,2) .and. (N1_from_type(type) == N2_from_type(type)) .and. (N1_from_type(type)*(N1_from_type(type)+1)/2 == size(from,1)) .and. &
!            (size(to,1) == size(to,3)) .and. (size(to,2) == size(to,4)) .and. (size(to,1) == N1_from_type(type)) &
!            .and. (size(to,2) == N2_from_type(type)) )) then
!
!            write(*,*) "superops_2indexed_to_4indexed_updiag - size error", size(from,1), size(from,2)
!            stop
!        end if
!
!        to = 0.0_dp
!
!        first_N  = N1_from_type(type)
!        second_N = N2_from_type(type)
!
!        supJ = 0
!
!        do l=1, first_N
!        do k=1, l
!
!        supI = 0
!        supJ = supJ + 1
!
!        do j=1, first_N
!        do i=1, j
!
!        supI = supI + 1
!
!        to(i,j,k,l) = from(supI,supJ)
!
!        end do
!        end do
!
!        end do
!        end do
!
!    end subroutine superops_2indexed_to_4indexed_updiag

!    subroutine redfield_from_evops(from, to, type, tstep)
!        complex(dpc), dimension(:,:,:,:,:), intent(in)    :: from
!        complex(dpc), dimension(:,:,:,:,:), intent(out)    :: to
!        character, intent(in)                                :: type
!        real(dp), intent(in)                                :: tstep
!
!        complex(dpc), dimension(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type)) :: XX, YY, ZZ
!        complex(dpc), dimension(:,:,:), allocatable :: xxx,yyy
!        real(dp), dimension(:), allocatable :: tmatrix
!        integer(i4b) :: i,j,k,l, NNN,  first_N, second_N
!
!        XX = 0.0_dp
!
!        first_N  = N1_from_type(type)
!        second_N = N2_from_type(type)
!
!        if(.not.(size(from,1) == first_N .or. size(from,2) == second_N .or. size(from,3) == first_N .or. size(from,4) == second_N &
!                .or. size(from,1) == size(to,1) .or. size(from,2) == size(to,2) .or. size(from,3) == size(to,3) .or. size(from,4) == size(to,4) .or. size(from,5) == size(to,5) ) ) then
!             call print_error_message(-1,  'wrong dimension in superops_to_exc_4indexed() : ')
!             stop
!         end if
!
!        ALLOCATE(xxx,(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type),size(to,5)))
!        ALLOCATE(yyy,(N1_from_type(type)*N2_from_type(type), N1_from_type(type)*N2_from_type(type),size(to,5)))
!        ALLOCATE(tmatrix,(size(to,5)))
!
!        do NNN = 1, size(to,5) ! do over time
!
!        do k=1, first_N
!        do l=1, second_N
!
!        do i=1, first_N
!        do j=1, second_N
!
!        XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N)) = &
!            from(i,j,k,l,NNN)
!
!        xxx(:,:,NNN) = XX(:,:)
!
!        end do
!        end do
!
!        end do
!        end do
!
!        tmatrix(NNN) = (tstep)*(NNN-1)
!
!        end do
!
!        XX = (xxx(:,:, 2) - xxx(:,:,1))/(tstep)
!        YY = 0.0_dp
!        call spline(tmatrix,xxx,XX,YY,yyy)
!
!        do NNN = 1, size(to,5) ! do over time
!
!        ! here is the procedure itself
!        call splint_matrix(tmatrix,xxx,yyy,tmatrix(NNN),XX)
!        call splint_matrix(tmatrix,xxx,yyy,tmatrix(NNN)+tstep/10.0_dp,YY)
!
!        ZZ = -(XX - YY)/(tstep/10.0_dp)
!        ZZ = (xxx(:,:, min(NNN+1,size(to,5)) ) - xxx(:,:,NNN))/(tstep)
!        YY = transpose(conjg(XX))
!        XX = matmul(ZZ,YY)
!
!        do k=1, first_N
!        do l=1, second_N
!
!        do i=1, first_N
!        do j=1, second_N
!
!        to(i,j,k,l,NNN) = XX(SUPERINDEX_FROM_K_L(i,j,second_N),SUPERINDEX_FROM_K_L(k,l,second_N))
!
!        end do
!        end do
!
!        end do
!        end do
!
!        end do ! do over time
!
!        DEALLOCATE(xxx)
!        DEALLOCATE(yyy)
!        DEALLOCATE(tmatrix)
!
!    end subroutine redfield_from_evops
!
!    !*************************************************************
!    !  Writing out evolution operators
!    !*************************************************************
!
!    subroutine write_evolution_operators(type)
!        character, intent(in) :: type
!        integer    (i4b)        :: i,j
!        integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2, Ublock
!        character(len=4)    :: no1,no2,no3,no4
!        character(len=100)    :: name
!        character(len=50)    :: prefix
!        complex(dpc), dimension(:,:,:,:,:), pointer        :: actual_U
!
!        Ublock = 1
!
!        ! We set indices range according to block we evaluate. Because rho0 is
!        ! whole density matrix, while evolution operators are only from particular
!        ! block, offset is set between these indices.
!        if (type == '2') then
!            actual_U => evops(Ublock,Ublock)%Ufe
!            prefix = 'Evops_fe'
!        else if (type == 'E') then
!            actual_U => evops(Ublock,Ublock)%Uee
!            prefix = 'Evops_ee'
!        else if (type == 'O') then
!            actual_U => evops(Ublock,Ublock)%Ueg
!            prefix = 'Evops_eg'
!        end if
!
!        do Uelement=1,N1_from_type(type)
!        do Uelement2=1,N2_from_type(type)
!        do Utnemele=1,N1_from_type(type)
!        do Utnemele2=1,N2_from_type(type)
!
!
!        if(Uelement < 10) then
!            write(no1,'(i1)')    Uelement
!        else if (Uelement < 100) then
!            write(no1,'(i2)')    Uelement
!        else
!            write(no1,'(i3)')    Uelement
!        endif
!        if(Uelement2 < 10) then
!            write(no2,'(i1)')    Uelement2
!        else if (Uelement2 < 100) then
!            write(no2,'(i2)')    Uelement2
!        else
!            write(no2,'(i3)')    Uelement2
!        endif
!        if(Utnemele < 10) then
!            write(no3,'(i1)')    Utnemele
!        else if (Uelement2 < 100) then
!            write(no3,'(i2)')    Utnemele
!        else
!            write(no3,'(i3)')    Utnemele
!        endif
!        if(Utnemele2 < 10) then
!            write(no4,'(i1)')    Utnemele2
!        else if (Uelement2 < 100) then
!            write(no4,'(i2)')    Utnemele2
!        else
!            write(no4,'(i3)')    Utnemele2
!        endif
!
!        name = trim(prefix) // trim(no1) // '-'//trim(no2)//'--'// trim(no3) // '-'//trim(no4)//'.dat'
!
!        open(UNIT=22, FILE = trim(file_join(out_dir,trim(name))))
!
!        i = 1
!        do while (i <= size(actual_U,5))
!            write(22,*) gt(1)*dt*(i-1),' ',real(actual_U(Uelement,Uelement2,Utnemele,Utnemele2,i)),' ',aimag(actual_U(Uelement,Uelement2,Utnemele,Utnemele2,i))
!            i = i + 1
!        end do
!
!        close(UNIT=22)
!
!        end do
!        end do
!        end do
!        end do
!    end subroutine write_evolution_operators
!
!

          subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t

            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
          end subroutine init_random_seed

    !
    ! Calculate linear absorption spectrum from polarization
    !
    subroutine create_spect_abs(rwa)
        !complex(dpc), dimension(:), intent(in)      :: in
        !complex(dpc), dimension(:), intent(out)     :: out
        real(dp), intent(in)                        ::  rwa

        complex(dpc), dimension(:),   allocatable   :: sig
        complex(dpc), dimension(:,:), allocatable   :: dat
        real(dp), dimension(:), allocatable :: rs

        complex(dpc)   :: aa, dd

        integer(i4b) :: i, NFFT, io
        real(dp) :: oma, dom, dt, a,b,c,d,e,f, lastf

        NFFT = 2048*8

        allocate(sig(NFFT))
        allocate(dat(1,NFFT))
        allocate(rs(NFFT))
     !   allocate(MC_spect_abs(NFFT))
        rs          = 0.0_dp
        !out         = 0.0_dp

        sig = 0.0_dp
        dt = 0.0_dp
        dom = 0.0_dp
        !sig = in

        CALL init_random_seed()
        lastf = 0.0;

        open(unit=9214,file='rhoE01-01.dat')
        open(unit=9215,file='rhoE01-02.dat')
        do i = 1, NFFT
            read(9214, *,IOSTAT=io) a, b, c
            read(9215, *,IOSTAT=io) d, e, f
            if(io > 0) then
              exit
            end if

            if(i == 2) then
              dt = a
              write(*,*) 'dt is', dt
              dom = 2.0_dp*PI_D/(NFFT*dt)
            end if

            !if(a > 159.1) then
            !    cycle
            !end if

            aa = b - cmplx(0,1)*c
            dd = e - cmplx(0,1)*f
            !call operator_from_exc(rhotmp,'O')
            !-0.881675* aa + 0.471858*dd, 0.471858 *aa + 0.881675 *dd
            !MC_polar_1(i) = MC_polar_1(i) + (dd(k,1)*dd(j,1))*as_orfact(k,j)*evops(kb,kb)%Ueg(k,1,j,1,i) !gcohs%C(i,k+kk)
            sig(i) = aa/1000000.0
            CALL RANDOM_NUMBER(f)
            !write(*,*) f
            !lastf = f * lastf + (1-f)*f*50
            lastf  = f
            sig(i) = ((lastf-0.5)*0.001*3 + aa*(1+(lastf-0.5)*0.01)  )/1000000.0
        end do
        close(9214)
        close(9215)

        do i = 1, NFFT/2
            sig(NFFT/2+i) = conjg(sig(NFFT/2-i+2))
        end do

        !
        ! FFT a la numerical recipes
        !
        dat(1,:) = sig(:)
        call fft_row(dat,1_i4b)
        sig(1:NFFT/2-1)  = dat(1,NFFT/2+2:NFFT)
        sig(NFFT/2:NFFT) = dat(1,1:NFFT/2+1)

        open(unit=9214,file='spect_abs.dat')

        do i = 1, NFFT
            oma = (-NFFT/2 + i)*dom*Energy_internal_to_cm + rwa
            rs(i) = abs(oma)*real(sig(i))
            write(9214, *) oma, rs(i)
        end do

        close(9214)

        deallocate(rs)
        deallocate(dat)
        deallocate(sig)

        !out = rs * dt

    end subroutine create_spect_abs

   !
    ! Priorities:
    !
    ! 1  ... highest
    !
    ! 10 ... lowest
    !

    !
    ! Prints message if its level is lower or equal to the log level
    !
    subroutine print_log_message(msg,ll)
        character(len=*), intent(in)  :: msg
        integer, intent(in)           :: ll
        !if (ll <= loglevel) then
        !    if (.not.logging) return
            write(*,'(a,a)') "Info   : ",msg
        !end if

        call flush()
    end subroutine print_log_message

    !
    ! Prints message if its level is lower or equal to the log level
    !
    subroutine print_warning_message(msg,ll)
        character(len=*), intent(in)  :: msg
        integer, intent(in)           :: ll
        !if (ll <= loglevel) then
        !    if (.not.logging) return
            write(*,'(a,a)') "Warning: ",msg
        !end if

        call flush()
    end subroutine print_warning_message


    !
    ! Prints error
    !
    subroutine print_error_message(err,errmsg)
        integer, intent(in)           :: err
        character(len=*), intent(in)  :: errmsg
        if (err == 0) return
        !if (loglevel > 0) then
        !    if (.not.logging) return
            write(*,'(a,a)') "Error  : ",errmsg
            stop
        !end if

        call flush()
    end subroutine print_error_message


 end module nakajima_zwanzig_shared
