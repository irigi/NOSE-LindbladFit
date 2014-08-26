
module lindblad_fit_module
    use std_types
    use helpers
    use nakajima_zwanzig_shared

    implicit none

    ! declarations
    real(dp), dimension(:,:), allocatable  :: Ham
    complex(dpc), allocatable:: rho1(:,:), prhox1(:,:), prhodx1(:,:)
    !complex(dpc), allocatable:: rho2(:,:), prhox2(:,:), prhodx2(:,:)
    complex(dpc), dimension(:,:,:,:,:), allocatable   :: Ueg, Uee, Ugg
    complex(dpc), dimension(:,:), private, allocatable   :: A
    complex(dpc), dimension(:), private, allocatable     :: B
    integer(i4b) :: Nl, STEPS
    integer(i4b), private :: Lr1, Lr2, Ls1, Ls2, Lbasis, LLr1, LLr2, LLs1, LLs2, LLbasis
    real(dp) :: timeStep = 0
    character(len=64), parameter, private :: external_dir = "external", config_filename = "config.prm", directory = "."

    ! basis multiplier
    real(dp), parameter :: lindblad_basis_multiplier = 0.001
    integer(i4b), parameter :: Nbasis = 99

    public::do_lindblad_fit_work

    private::init_lindblad_fit
    private::clean_lindblad_fit
    private::read_evops
    private::open_files
    private::close_files
    private::random_test
    private::read_config_file
    private::Lmult1
    private::propagate1
    !private::Lmult2
    !private::propagate2

    private::read_Nsys
    private::ind

    public ::indices_to_superindex
    public ::superindex_to_indices

    contains

    !
    ! Do all the simulations within the module
    !
    subroutine do_lindblad_fit_work()
        integer(i4b) :: r,s,k,l,i


        call init_lindblad_fit()
        call init_nakajima_zwanzig_shared(Ham)

        !write(*,*) indices_to_superindex(1,2,2,1,351)
        !call superindex_to_indices(indices_to_superindex(1,2,2,1,351), r,s,k,l,i)
        !write(*,*) r,s,k,l,i
        !write(*,*) 1,2,2,1,351

        call open_files('E')
        call read_evops('E', Uee)
        call close_files('E')

        call create_hustomatrix(A,B)


    end subroutine do_lindblad_fit_work

    subroutine superindex_to_indices(super, i,j,k,l, tind)
      integer(i4b), intent(in)  :: super
      integer(i4b), intent(out) :: i, j, k, l, tind

      integer(i4b)              :: s

      s = super - 1

      i = mod(s, Nl) + 1
      s = (s - (i-1)) / Nl

      j = mod(s, Nl) + 1
      s = (s - (j-1)) / Nl

      k = mod(s, Nl) + 1
      s = (s - (k-1)) / Nl

      l = mod(s, Nl) + 1
      s = (s - (l-1)) / Nl

      tind = s + 1


    end subroutine superindex_to_indices

    integer(i4b) function indices_to_superindex(i,j,k,l, tind) result(super)
      integer(i4b), intent(in) :: i,j,k,l, tind

      super = (i-1) + (j-1)*Nl + (k-1)*Nl*Nl + (l-1)*Nl*Nl*Nl + (tind-1)*Nl*Nl*Nl*Nl + 1
    end function indices_to_superindex

    subroutine create_hustomatrix(A,B)
      complex(dpc), dimension(:,:), intent(out) :: A
      complex(dpc), dimension(:), intent(out)   :: B

      integer(i4b) :: i,j,k,l, tind, dummy, super1, super2

      complex(dpc), dimension(:, :,:,:,:), allocatable :: calc

      write(*,*) size(A,1), size(A,2), size(B,1)

      A = 0.0_dp
      B = 0.0_dp

      write(*,*) 'allocating', Nl*Nl*Nl*Nl*Nbasis*Nl*Nl*Nl*Nl
      allocate(calc(Nl*Nl*Nl*Nl*Nbasis,Nl,Nl,Nl,Nl))
      calc = 0.0_dp

      do tind=1,STEPS
        write(*,*) tind, 'of', STEPS
        call cache_lindblad_basis(calc, tind)

        ! cycles over basis-functions
        do Lr1=1, Nl
        do Ls1=1, Nl
        do Lr2=1, Nl
        do Ls2=1, Nl
        do Lbasis=1,Nbasis

        ! cycles over basis-functions
        do LLbasis=1,Nbasis
        do LLs2=1, Nl
        do LLr2=1, Nl
        do LLs1=1, Nl
        do LLr1=1, Nl

         super1 = indices_to_superindex(Lr1,Ls1,Lr2,Ls2,Lbasis)
         super2 = indices_to_superindex(LLr1,LLs1,LLr2,LLs2,LLbasis)
         !write(*,*) super1, super2

        ! cycles over "false-time"
        do i=1, Nl
        do j=1, Nl
        do k=1, Nl
        do l=1, Nl
          A(super1,super2) = A(super1,super2) + calc(super1,i,j,k,l)*conjg(calc(super2,i,j,k,l))

          if(super1 == 1) then
            B(super2) = B(super2) + Uee(i,j,k,l,tind)*conjg(calc(super2,i,j,k,l))
          end if
        end do
        end do
        end do
        end do

        end do
        end do
        end do
        end do
        end do

        end do
        end do
        end do
        end do
        end do
      ! over time
      end do

      deallocate(calc)

    end subroutine create_hustomatrix

    subroutine cache_lindblad_basis(calc, tind)
      complex(dpc), dimension(:, :,:,:,:), intent(inout) :: calc

      integer(i4b), intent(in) :: tind
      integer(i4b)             :: i,j,k,l,dummy
      real(dp)                 :: time

      time = tind*timeStep

      if(tind == 1) then
        calc = 0.0_dp

        ! cycles over basis-functions
        do Lr1=1, Nl
        do Ls1=1, Nl
        do Lr2=1, Nl
        do Ls2=1, Nl
        do Lbasis=1,Nbasis

        ! cycles over false time
        do i=1, Nl
        do j=1, Nl
          calc(indices_to_superindex(Lr1,Ls1,Lr2,Ls2,Lbasis),i,j,i,j) = 1.0_dp
        end do
        end do

        end do
        end do
        end do
        end do
        end do
      end if

      ! cycles over basis-functions
      do Lr1=1, Nl
      do Ls1=1, Nl
      do Lr2=1, Nl
      do Ls2=1, Nl
      do Lbasis=1,Nbasis

      ! copy to make time step
      do i=1, Nl
      do j=1, Nl
        rho1(:,:) = calc(indices_to_superindex(Lr1,Ls1,Lr2,Ls2,Lbasis),i,j,:,:)

        do dummy=1,10
          call propagate1(0.1_dp*timeStep,time + dummy*0.1*timeStep)
        end do

        calc(indices_to_superindex(Lr1,Ls1,Lr2,Ls2,Lbasis),i,j,:,:) = rho1(:,:)
      end do
      end do

      end do
      end do
      end do
      end do
      end do
    end subroutine cache_lindblad_basis

    subroutine init_lindblad_fit()
        STEPS = 2000

        Nl = read_Nsys()

        allocate(Ueg(Nl, 1,Nl, 1,STEPS) )
        allocate(Uee(Nl,Nl,Nl,Nl,STEPS) )
        allocate(Ugg(1,1,1,1,STEPS) )

        allocate(Ham(Nl,Nl))

        allocate(rho1(Nl,Nl))
        allocate(prhox1(Nl,Nl))
        allocate(prhodx1(Nl,Nl))

!        allocate(rho2(Nl,Nl))
!        allocate(prhox2(Nl,Nl))
!        allocate(prhodx2(Nl,Nl))

        allocate(A(Nl*Nl*Nl*Nl*Nbasis,Nl*Nl*Nl*Nl*Nbasis))
        allocate(B(Nl*Nl*Nl*Nl*Nbasis))

        Ueg = 0.0_dp
        Uee = 0.0_dp
        Ugg = 0.0_dp

        Ham = 0.0_dp

        rho1    = 0.0_dp
        prhox1  = 0.0_dp
        prhodx1 = 0.0_dp

        !rho2     = 0.0_dp
        !prhox2   = 0.0_dp
        !prhodx2  = 0.0_dp

        A = 0.0_dp
        B = 0.0_dp


        call read_config_file(Nl)

        write(*,*) ';Ham', Ham

    end subroutine init_lindblad_fit

    subroutine clean_lindblad_fit()
    end subroutine clean_lindblad_fit

    subroutine Lmult1 (tt, rhoin, result)
      complex(dpc), intent(in)   :: rhoin(:,:)
      real(dp), intent(in)       :: tt
      complex(dpc), intent(out)  :: result(:,:)
      integer(i4b) :: i
      complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)

      result = 0.0_dp

      result = result - iconst*MATMUL(Ham, rhoin) + iconst*MATMUL(rhoin, Ham)

      ! global Lr1, Lr2, Ls1, Ls1 must be set to desired basis element!

      ! La rho Lb+
      result(Lr1,Lr2) = result(Lr1,Lr2) + rhoin(Ls1,Ls2)     * lindblad_basis_multiplier
      if(Lr2 == Lr1) then
        do i=1,Nl
          ! - Lb+ La rho / 2
          result(Ls2,i) = result(Ls2,i) - rhoin(Ls1,i) / 2   * lindblad_basis_multiplier

          ! - rho Lb+ La / 2
          result(i,Ls1) = result(i,Ls1) - rhoin(i,Ls2) / 2   * lindblad_basis_multiplier
        end do
      end if

      result = result * arbitrary_basis_element(tt)

    end subroutine Lmult1

    real(dp) function arbitrary_basis_element(t) result(nav)
      real(dp), intent(in) :: t

      real(dp) :: omega
      integer(i4b) :: md
      integer(i4b) :: dv

      dv = (Nbasis-1) / 3
      md = mod(Nbasis-1, 3)
      omega = dv/1000

      if(md == 0) then
        nav = sin(omega*t)
      elseif(md == 1) then
        nav = cos(omega*t)
      else
        nav = exp(-omega*t)
      end if

    end function arbitrary_basis_element

    subroutine propagate1(dt,t)
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: t

      call Lmult1(t,rho1,prhodx1)
      call ode_rk4(rho1,prhodx1,t,dt,prhox1,Lmult1)

      rho1 = prhox1
    end subroutine propagate1

!    subroutine Lmult2 (tt, rhoin, result)
!      complex(dpc), intent(in)   :: rhoin(:,:)
!      real(dp), intent(in)       :: tt ! this is a dummy parameter to satisfy ode_rk4 function template
!      complex(dpc), intent(out)  :: result(:,:)
!      integer(i4b) :: i
!      complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)
!
!      result = 0.0_dp
!
!      ! we calculate this in exciton basis
!      result = result - iconst*MATMUL(Ham, rhoin) + iconst*MATMUL(rhoin, Ham)
!
!      ! global Lr1, Lr2, Ls1, Ls1 must be set to desired basis element!
!
!      ! LLa rho LLb+
!      result(LLr1,LLr2) = result(LLr1,LLr2) + rhoin(LLs1,LLs2)  * lindblad_basis_multiplier
!      if(LLr2 == LLr1) then
!        do i=1,Nl
!          ! - LLb+ LLa rho / 2
!          result(LLs2,i) = result(LLs2,i) - rhoin(LLs1,i) / 2   * lindblad_basis_multiplier
!
!          ! - rho LLb+ LLa / 2
!          result(i,LLs1) = result(i,LLs1) - rhoin(i,LLs2) / 2   * lindblad_basis_multiplier
!        end do
!      end if
!
!    end subroutine Lmult2
!
!    subroutine propagate2(dt)
!      real(dp), intent(in) :: dt
!
!      real(dp) :: t
!      t = dt ! this is a dummy parameter to satisfy ode_rk4 function template
!
!      call Lmult2(t,rho2,prhodx2)
!      call ode_rk4(rho2,prhodx2,t,dt,prhox2,Lmult2)
!
!      rho2 = prhox2
!    end subroutine propagate2

    pure function ind(i,j,k,l,read) result(res)
      integer(i4b) :: res
      integer(i4b), intent(in) :: i,j,k,l
      logical, intent(in) :: read

      res = 22 + i+Nl*(j-1)+Nl*Nl*(k-1)+Nl*Nl*Nl*(l-1)

      if(read) then
        res = res + 10000
      end if

    end function ind

    subroutine read_evops(type, actual_U)
      character, intent(in)      :: type
      complex(dpc), dimension(:,:,:,:,:), intent(out)     :: actual_U

      integer (i4b)       :: i, file_ios
      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2

      real(dp)            :: a, b, time

      do i=1,size(actual_U,5)
      do Uelement=1,N1_from_type(type)
      do Uelement2=1,N2_from_type(type)
      do Utnemele=1,N1_from_type(type)
      do Utnemele2=1,N2_from_type(type)

        read(ind(Uelement,Uelement2,Utnemele,Utnemele2,.false.),*, IOSTAT=file_ios) time, a, b
        actual_U(Uelement,Uelement2,Utnemele,Utnemele2,i) = a + b * cmplx(0,1)

        if(i == 1 .and. Uelement == 1 .and. Uelement2 == 1 .and. Utnemele == 1 .and. Utnemele2 == 1) then
          timeStep = time
        elseif(i == 2 .and. Uelement == 1 .and. Uelement2 == 1 .and. Utnemele == 1 .and. Utnemele2 == 1) then
          timeStep = time - timeStep
          write(*,*) 'timeStep', timeStep
        end if

      end do
      end do
      end do
      end do
      end do
    end subroutine read_evops

    subroutine open_files(type)
      character, intent(in) :: type

      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2, Ublock
      character(len=4)    :: no1,no2,no3,no4
      character(len=100)  :: name
      character(len=50)   :: prefix


      Ublock = 1

      ! We set indices range according to block we evaluate. Because rho0 is
      ! whole density matrix, while evolution operators are only from particular
      ! block, offset is set between these indices.
      if (type == '2') then
        !actual_U => evops(Ublock,Ublock)%Ufe
        prefix = 'Evops_fe'
      else if (type == 'E') then
        !actual_U => evops(Ublock,Ublock)%Uee
        prefix = 'Evops_ee'
      else if (type == 'O') then
        !actual_U => evops(Ublock,Ublock)%Ueg
        prefix = 'Evops_eg'
      end if

      do Uelement=1,N1_from_type(type)
      do Uelement2=1,N2_from_type(type)
      do Utnemele=1,N1_from_type(type)
      do Utnemele2=1,N2_from_type(type)


      if(Uelement < 10) then
        write(no1,'(i1)')   Uelement
      else if (Uelement < 100) then
        write(no1,'(i2)')   Uelement
      else
        write(no1,'(i3)')   Uelement
      endif
      if(Uelement2 < 10) then
        write(no2,'(i1)')   Uelement2
      else if (Uelement2 < 100) then
        write(no2,'(i2)')   Uelement2
      else
        write(no2,'(i3)')   Uelement2
      endif
      if(Utnemele < 10) then
        write(no3,'(i1)')   Utnemele
      else if (Uelement2 < 100) then
        write(no3,'(i2)')   Utnemele
      else
        write(no3,'(i3)')   Utnemele
      endif
      if(Utnemele2 < 10) then
        write(no4,'(i1)')   Utnemele2
      else if (Uelement2 < 100) then
        write(no4,'(i2)')   Utnemele2
      else
        write(no4,'(i3)')   Utnemele2
      endif

      name = trim(prefix) // trim(no1) // '-'//trim(no2)//'--'// trim(no3) // '-'//trim(no4)//'.dat'
      open(UNIT=ind(Uelement,Uelement2,Utnemele,Utnemele2,.false.), FILE = trim(name), STATUS='OLD', ACTION='READ')

      end do
      end do
      end do
      end do
    end subroutine open_files

    subroutine close_files(type)
      character, intent(in) :: type
      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2

      do Uelement=1,N1_from_type(type)
      do Uelement2=1,N2_from_type(type)
      do Utnemele=1,N1_from_type(type)
      do Utnemele2=1,N2_from_type(type)

      close(UNIT=ind(Uelement,Uelement2,Utnemele,Utnemele2,.false.))

      end do
      end do
      end do
      end do
    end subroutine close_files

    subroutine random_test()
      real(dp)     :: r, rr, rrr
      integer(i4b) :: i,j

      write(*,*) 'Montecarlo'

      call init_random_seed()

      call random_number(rr)
      call random_number(rrr)

      i = 0
      j = 0
      do while (rr < 2)
         i = i + 1
         call random_number(r)

         if (r == rr) then
           write(*,*) "perioda nalezena po ", i
           stop
         end if

         if (mod(i,100000000) == 0) then
           i = 0
           j = j + 1
           write(*,*) j
         end if
      end do

      stop
    end subroutine random_test

    integer(i4b) function read_Nsys() result(Nsys)
        character(len=256)           :: buff = ""
        real(dp)                     :: svalue
        integer(i4b)                 :: i = 0

        open(unit=32,file=trim(trim(directory)//'/'//trim(config_filename) ) , err=32, status='old')

        Nsys = 0

        do while(i == 0)
          read(32, *, iostat=i) buff

        !global_en, global_lambda, global_temp, global_tau

          if(trim(adjustl(buff)) == 'systemSize') then
            read(32, *, iostat=i) svalue
            write(*,*) buff, int(svalue)
            Nsys = int(svalue)

          else
            !write(*,*) 'unrecognised:', buff, value
          end if
        end do

        close(32)

        return

32      write(*,*) "couldn't read the supplementary config file"
        stop
    end function read_Nsys

    subroutine read_config_file(Nsys)
        integer(i4b), intent(in)     :: Nsys
        character(len=256)           :: buff = ""
        real(dp), dimension(Nsys)    :: value
        real(dp)                     :: svalue
        integer(i4b)                 :: i = 0, j

        open(unit=32,file=trim(trim(directory)//'/'//trim(config_filename) ) , err=32, status='old')

        do while(i == 0)
          read(32, *, iostat=i) buff


          if(trim(adjustl(buff)) == 'hamiltonian') then
            do j=1,Nsys
              read(32, *, iostat=i) value
              write(*,*) buff, j, value
              Ham(j,1:Nsys) = value(1:Nsys)/Energy_internal_to_cm
            end do

            buff = "-"

          else
            write(*,*) 'unrecognised:', buff, value
          end if
        end do

        close(32)

        return

32      write(*,*) "couldn't read the supplementary config file"
        stop
    end subroutine read_config_file


end module lindblad_fit_module

