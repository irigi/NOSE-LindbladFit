module module_goft
    use std_types
    use helpers
    use nakajima_zwanzig_shared
!    use resources_montecarlo

    implicit none

    complex(dpc), dimension(:,:), allocatable :: goft, coft, hoft

    contains

    subroutine init_goft(NSys, steps, lambda, tau, temp, dt)
        integer(i4b), intent(in) :: Nsys, steps
        real(dp), dimension(:), intent(in)   :: lambda, tau, temp
        real(dp), intent(in)   :: dt
        character(len=256) buff

        integer(i4b) :: i, j

        allocate(goft(Nsys, steps))
        allocate(hoft(Nsys, steps))
        allocate(coft(Nsys, steps))
        goft = 0.0
        hoft = 0.0
        coft = 0.0

        do i = 1, Nsys
            call brownian(dt,temp(i),tau(i),lambda(i),goft(i,:),coft(i,:),hoft(i,:))

            write(buff,'(A I1 A)') "goft",i,".dat"
            open(UNIT=1444, FILE = trim(buff))
            do j = 1, steps
                write(1444,*) dt*j," ", real(goft(i,j))," ", aimag(goft(i,j))
            end do
            close(1444)

            write(buff,'(A I1 A)') "coft",i,".dat"
            open(UNIT=1444, FILE = trim(buff))
            do j = 1, steps
                write(1444,*) dt*j," ", real(coft(i,j))," ", aimag(coft(i,j))
            end do
            close(1444)
        end do
    end subroutine init_goft

    !
    ! Complex tanh
    !
    function ctanh(x) result(nav)
        complex(dpc), intent(in)    :: x
        complex(dpc)                :: nav

        nav = (exp(x)-exp(-x))/(exp(x)+exp(-x))
    end function ctanh

    !
    ! Dimensionless brownian overdamped CC(T,x),
    !    x = beta hbar Lambda/2
    !    T = 2 t / beta hbar
    !    C(Lambda, ll, t) = Lambda ll CC(T,x)
    !
    !    CC(T,x) = (cot(x)-i) e^(-Tx) + sum_n=1^\infty e^(-n pi T) (1/(x + n pi) - 1/(x - n pi))
    !
    function dimensionless_CC(T,x) result(CC)
        real(dp), intent(in)        :: x, T
        complex(dpc)                :: CC

        real(dp)                    :: y, diff
        integer(i4b)                :: nearest_pole,i
        character(len=256)          :: buff

        ! main idea of this function is to evaluate smooth CC(T,x) by
        ! proper handling of summation of the Matsubara sum by switching
        ! between Laurent expansion around singularities of cot(x),
        ! which are regularized by proper terms of the sum, and
        ! evaluating cot(x) exactly outside the singularities

        ! we shift x to the nearest singularity
        nearest_pole = INT((x + PI_D/2) / (PI_D))
        y = x - nearest_pole*PI_D

        ! if assures that the message is written only once - at time = 0
        if(T < 1e-5) then
            write(buff,'(a,f6.3,a)') 'brownian cf taken in ',x/PI_D,', where 1 = nearest singularity'

            if(x/PI_D < 0.8) then
                call print_log_message(trim(buff),7)
            else
                call print_warning_message(trim(buff),5)
            end if
        end if

        ! if we are +-1 around pole, Laurent expansion works well, otherwise
        ! we are far enough and we don't need it
        if(abs(y) < 0.1 .and. nearest_pole > 0) then

            if(T < 1e-5) then
                call print_warning_message('interpolating cf around pole',5)
            end if

            ! series of e^(-T*y) cot y - 1/y
            CC = 0
            CC = CC + (y**3) * (15*(T**4) - 60*(T**2) - 8.0_dp)/360.0_dp
            CC = CC + (y**2) * (T/3 - (T**3)/6)
            CC = CC + (y**1) * (-1.0_dp/3.0_dp + T**2/2.0_dp)
            CC = CC - T

            ! second partial fraction (+1/(y + 2 pi n)) is added exactly
            CC = CC + 1.0_dp/(y + 2.0_dp*PI_D*nearest_pole)

            ! imaginary part is added exactly
            CC = CC - exp(-T*y)*cmplx(0,1,dpc)

            ! whole previous CC is multiplied by exp(-pi n T)
            CC = CC * exp(- PI_D*nearest_pole*T)

        else

            CC =      exp(-T*x)*(1.0_dp/tan(x) - cmplx(0,1,dpc))

            if(nearest_pole > 0) then
                CC = CC + exp(- PI_D*nearest_pole*T)*                       &
                            (1.0_dp/(x+nearest_pole*PI_D)-1.0_dp/(x-nearest_pole*PI_D))
            end if
        end if

        ! the rest of Matsubara sum is added, except the nearest pole
        do i = 1, 2000
            if(i == nearest_pole) then
                cycle
            end if

            diff = exp(- PI_D*i*T)*(1.0_dp/(x+i*PI_D)-1.0_dp/(x-i*PI_D))
            CC = CC + diff

            if(diff/abs(CC) < 0.001_dp) then
                exit
            end if
        end do

        if(T < 1e-5) then
            write(buff,'(a,i4)') 'Matsubara sum to ',i-1
            call print_log_message(trim(buff),7)
        end if

    end function dimensionless_CC

    subroutine brownian(dt,temp,tau,lambda,ggt,cct,hht,ADD)
        real(dp), intent(in)                    :: dt, temp, tau, lambda
        complex(dpc), dimension(:), intent(out) :: ggt
        complex(dpc), dimension(:), intent(out) :: cct,hht
        character(len=*), intent(in), optional  :: ADD

        complex(dpc), dimension(size(ggt))      :: cct_tmp,hht_tmp,ggt_tmp
        real(dp)                                :: BH, LLambda,t
        integer(i4b)                            :: Ntt, i

        !
        ! Set evaluation parameters
        !
        LLambda     = 1.0_dp/tau

        ! dt ... elementary grid step

        Ntt = size(ggt)

        BH = (0.6582120_dp/8.617385d-5)/temp  ! hbar/kT

        do i=1, Ntt
            t = (i-1)*dt

            cct_tmp(i) = lambda*LLambda*dimensionless_CC(2.0*t/BH, BH*LLambda/2.0)
            if(i > 1) then
                hht_tmp(i) = hht_tmp(i-1) + dt*cct_tmp(i)
                ggt_tmp(i) = ggt_tmp(i-1) + dt*hht_tmp(i)
            else
                hht_tmp(i) = dt*cct_tmp(i)
                ggt_tmp(i) = dt*hht_tmp(i)
            end if
        end do

!       open(unit=11,file='/home/olsij4am/prace/nose-debug.dat')
!       open(unit=12,file='/home/olsij4am/prace/nose-debug2.dat')

        ! write to global functions
        if (.not. present(ADD)) then
            ggt = 0.0_dp
            cct = 0.0_dp
            hht = 0.0_dp
        end if

        do i=1,Ntt
 !          write(11,*) i*dt,real(cct_tmp(i)),real(hht_tmp(i)),real(ggt_tmp(i))
  !         write(12,*) i*dt,aimag(cct_tmp(i)),aimag(hht_tmp(i)),aimag(ggt_tmp(i))

            cct(i) = cct_tmp(i) + cct(i)
            hht(i) = hht_tmp(i) + hht(i)
            ggt(i) = ggt_tmp(i) + ggt(i)
        end do

!       close(11)
!       close(12)
!       write(*,*) 'debug functions written'
!       stop

    end subroutine brownian

end module module_goft
