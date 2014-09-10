
module lindblad_fit_module
    use std_types
    use helpers
    use nakajima_zwanzig_shared
    use numer_matrix

    implicit none
    private

    ! declarations
    real(dp), dimension(:,:), allocatable  :: Ham
    complex(dpc), allocatable:: rho1(:,:), prhox1(:,:), prhodx1(:,:)
    complex(dpc), dimension(:,:,:,:,:), allocatable      :: Evops, Devops
    complex(dpc), dimension(:,:), allocatable            :: Evops2, Devops2, AAA, VT, U
    real(dp), dimension(:), allocatable                  :: eigval
    integer(i4b)         :: Nl1, Nl2, Nl
    character, parameter :: type = 'E'
    real(dp)             :: timeStep = 0
    character(len=64), parameter, private :: external_dir = "external", config_filename = "config.prm", directory = "."

    ! basis multiplier
    real(dp), parameter :: lindblad_basis_multiplier = 0.001
    integer(i4b), public :: Nbasis = 99, STEPS = 500

    logical :: to_exciton_at_output = .false.

    public::do_lindblad_fit_work
    public::indices_to_superindex
    public::superindex_to_indices
    public::only_convert_to_exciton
    public::rates_to_evops

    contains

    !
    ! Do all the simulations within the module
    !
    subroutine do_lindblad_fit_work()
        integer(i4b) :: r,s,k,l,i

        call init_lindblad_fit()

        !write(*,*) 'READING EXTERNAL EVOPS'
        call flush()
        call open_files('r')
        call read_evops()
        call close_files('r')

        call evops_derivative()
        call calculate_coeff()

        !write(*,*) 'OUTPUTTING EVOPS'
        call flush()
        call open_files('w')
        call write_evops()
        call close_files('w')

        !write(*,*) 'OUTPUTTING EVOPS'
        call flush()
        call open_files('D')
        call write_devops()
        call close_files('D')

    end subroutine do_lindblad_fit_work

    subroutine calculate_coeff()
        integer(i4b) :: tind, i

        do tind=1,STEPS
          call superops_4indexed_to_2indexed(Evops(:,:,:,:,tind),Evops2,type)
          call superops_4indexed_to_2indexed(DEvops(:,:,:,:,tind),DEvops2,type)

          call svd(Evops2,U,EIGVAL,VT)

          !write(*,*) EIGVAL

          if(EIGVAL(1) <= 0) then
            write(*,*) 'terrible error in SVD'
            stop
          end if

          ! reciprocal eigvals
          do i=1, size(EIGVAL)
            if(EIGVAL(i)/EIGVAL(1) > 1e-7*size(EIGVAL)) then
              EIGVAL(i) = 1.0_dp / EIGVAL(i)
            else
              EIGVAL(i) = 0.0_dp
            end if
          end do

          !write(*,*) EIGVAL

          AAA = 0.0_dp
          do i=1, size(AAA,1)
            AAA(i,i) = EIGVAL(i)
          end do

          AAA = matmul(conjg(transpose(VT)),matmul(AAA, transpose(conjg(U))))

          AAA = matmul(dEvops2,AAA)

          call superops_2indexed_to_4indexed(AAA,Evops(:,:,:,:,tind),type)


        end do
    end subroutine calculate_coeff

    subroutine only_convert_to_exciton()
        integer(i4b) :: r

        call init_lindblad_fit()

        write(*,*) 'READING EXTERNAL EVOPS'
        call flush()
        call open_files('r')
        call read_evops()
        call close_files('r')

        do r=1,STEPS
          call superops_to_exc_4indexed(Evops(:,:,:,:,r),type)
        end do

        write(*,*) 'OUTPUTTING EVOPS'
        call flush()
        call open_files('w')
        call write_evops()
        call close_files('w')

    end subroutine only_convert_to_exciton

    subroutine superindex_to_indices(super, i,j,k,l, tind)
      integer(i4b), intent(in)  :: super
      integer(i4b), intent(out) :: i, j, k, l, tind

      integer(i4b)              :: s

      s = super - 1

      i = mod(s, Nl1) + 1
      s = (s - (i-1)) / Nl1

      j = mod(s, Nl2) + 1
      s = (s - (j-1)) / Nl2

      k = mod(s, Nl1) + 1
      s = (s - (k-1)) / Nl1

      l = mod(s, Nl2) + 1
      s = (s - (l-1)) / Nl2

      tind = s + 1


    end subroutine superindex_to_indices

    integer(i4b) function indices_to_superindex(i,j,k,l, tind) result(super)
      integer(i4b), intent(in) :: i,j,k,l, tind

      super = (i-1) + (j-1)*Nl1 + (k-1)*Nl1*Nl2 + (l-1)*Nl1*Nl2*Nl1 + (tind-1)*Nl1*Nl2*Nl1*Nl2 + 1
    end function indices_to_superindex

    subroutine init_lindblad_fit()
        Nl = read_Nsys()

        allocate(Ham(Nl,Nl))

        Ham = 0.0_dp

        call read_config_file(Nl)

        write(*,*) ';Ham', Ham

        call init_nakajima_zwanzig_shared(Ham)
        Nl1 = N1_from_type(type)
        Nl2 = N2_from_type(type)


        allocate(rho1(Nl1,Nl2))
        allocate(prhox1(Nl1,Nl2))
        allocate(prhodx1(Nl1,Nl2))

        allocate(Evops(Nl1,Nl2,Nl1,Nl2,STEPS) )
        allocate(DEvops(Nl1,Nl2,Nl1,Nl2,STEPS) )

        allocate(Evops2(Nl1*Nl2,Nl1*Nl2) )
        allocate(DEvops2(Nl1*Nl2,Nl1*Nl2) )
        allocate(AAA(Nl1*Nl2,Nl1*Nl2) )
        allocate(VT(Nl1*Nl2,Nl1*Nl2) )
        allocate(U(Nl1*Nl2,Nl1*Nl2) )
        allocate(eigval(Nl1*Nl2) )

        rho1    = 0.0_dp
        prhox1  = 0.0_dp
        prhodx1 = 0.0_dp

        Evops = 0.0_dp
        DEvops = 0.0_dp

        Evops2 = 0.0_dp
        DEvops2 = 0.0_dp
        AAA = 0.0_dp
        VT = 0.0_dp
        U = 0.0_dp
        eigval = 0.0_dp

    end subroutine init_lindblad_fit

    subroutine clean_lindblad_fit()
    end subroutine clean_lindblad_fit

    subroutine Lmult1 (t, rhoin, res)
      complex(dpc), intent(in)   :: rhoin(:,:)
      real(dp), intent(in)       :: t
      complex(dpc), intent(out)  :: res(:,:)
      complex(dpc), parameter:: iconst = dcmplx(0.0, 1.0)

      res = 0.0_dp

      if(type == 'E') then
        res = res - iconst*MATMUL(Ham, rhoin) + iconst*MATMUL(rhoin, Ham)
      elseif(type == 'O') then
        res = res - iconst*MATMUL(Ham, rhoin) !+ iconst*MATMUL(rhoin, Ham)
      end if

    end subroutine Lmult1

    subroutine propagate1(dt,t)
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: t

      call Lmult1(t,rho1,prhodx1)
      call ode_rk4(rho1,prhodx1,t,dt,prhox1,Lmult1)

      rho1 = prhox1
    end subroutine propagate1

    pure function ind(i,j,k,l,code) result(res)
      integer(i4b) :: res
      integer(i4b), intent(in) :: i,j,k,l
      character, intent(in) :: code

      res = 22 + i+Nl1*(j-1)+Nl1*Nl2*(k-1)+Nl1*Nl2*Nl1*(l-1)

      if(code == 'r') then
        res = res + 10000
      elseif(code == 'D') then
        res = res + 20000
      end if

    end function ind

    subroutine evops_derivative()
      integer(i4b)        :: i,j,k,l, tind

      do i=1,Nl1
      do j=1,Nl2
      do k=1,Nl1
      do l=1,Nl2
      do tind=2,size(Evops,5)-1

        DEvops(i,j,k,l,tind) = (Evops(i,j,k,l,tind+1) - Evops(i,j,k,l,tind-1)) / timeStep / 2

      end do

      DEvops(i,j,k,l,1) = DEvops(i,j,k,l,2) - (DEvops(i,j,k,l,3)-DEvops(i,j,k,l,2))
      DEvops(i,j,k,l,size(Evops,5)) = DEvops(i,j,k,l,size(Evops,5)-1) + (DEvops(i,j,k,l,size(Evops,5)-1)-DEvops(i,j,k,l,size(Evops,5)-2))

      end do
      end do
      end do
      end do

    end subroutine evops_derivative

    subroutine read_evops()
      integer (i4b)       :: i, file_ios
      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2

      real(dp)            :: a, b, time

      do i=1,size(Evops,5)
      do Uelement=1,Nl1
      do Uelement2=1,Nl2
      do Utnemele=1,Nl1
      do Utnemele2=1,Nl2

        read(ind(Uelement,Uelement2,Utnemele,Utnemele2,'r'),*, IOSTAT=file_ios) time, a, b
        Evops(Uelement,Uelement2,Utnemele,Utnemele2,i) = a + b * cmplx(0,1)

        ! read timestep from the Evops
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

    subroutine write_evops()
      integer (i4b)       :: i, file_ios
      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2

      do i=1,size(Evops,5)
      do Uelement=1,Nl1
      do Uelement2=1,Nl2
      do Utnemele=1,Nl1
      do Utnemele2=1,Nl2


        write(ind(Uelement,Uelement2,Utnemele,Utnemele2,'w'),*, IOSTAT=file_ios) timeStep*(i-1),              &
                        real(Evops(Uelement,Uelement2,Utnemele,Utnemele2,i)),                                 &
                        aimag(Evops(Uelement,Uelement2,Utnemele,Utnemele2,i))

      end do
      end do
      end do
      end do
      end do
    end subroutine write_evops

    subroutine write_devops()
      integer (i4b)       :: i, file_ios
      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2

      do i=1,size(dEvops,5)
      do Uelement=1,Nl1
      do Uelement2=1,Nl2
      do Utnemele=1,Nl1
      do Utnemele2=1,Nl2


        write(ind(Uelement,Uelement2,Utnemele,Utnemele2,'D'),*, IOSTAT=file_ios) timeStep*(i-1),              &
                        real(dEvops(Uelement,Uelement2,Utnemele,Utnemele2,i)),                                 &
                        aimag(dEvops(Uelement,Uelement2,Utnemele,Utnemele2,i))

      end do
      end do
      end do
      end do
      end do
    end subroutine write_devops

    subroutine open_files(code)
      character, intent(in) :: code

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

      do Uelement=1,Nl1
      do Uelement2=1,Nl2
      do Utnemele=1,Nl1
      do Utnemele2=1,Nl2


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

      if(code == 'r') then

        name = trim(prefix) // trim(no1) // '-'//trim(no2)//'--'// trim(no3) // '-'//trim(no4)//'.dat'
        open(UNIT=ind(Uelement,Uelement2,Utnemele,Utnemele2,code), FILE = trim(name), STATUS='OLD', ACTION='READ')

      elseif(code == 'w') then

        name = trim(prefix) // trim(no1) // '-'//trim(no2)//'--'// trim(no3) // '-'//trim(no4)//'-out.dat'
        open(UNIT=ind(Uelement,Uelement2,Utnemele,Utnemele2,code), FILE = trim(name))

      elseif(code == 'D') then
        if (type == '2') then
          prefix = 'Diss_fe'
        else if (type == 'E') then
          prefix = 'Diss_ee'
        else if (type == 'O') then
          prefix = 'Diss_eg'
        end if

        name = trim(prefix) // trim(no1) // '-'//trim(no2)//'--'// trim(no3) // '-'//trim(no4)//'-out.dat'
        open(UNIT=ind(Uelement,Uelement2,Utnemele,Utnemele2,code), FILE = trim(name))

      end if

      end do
      end do
      end do
      end do
    end subroutine open_files

    subroutine close_files(code)
      character, intent(in) :: code

      integer(i4b)        :: Uelement, Uelement2,Utnemele,Utnemele2

      do Uelement=1,Nl1
      do Uelement2=1,Nl2
      do Utnemele=1,Nl1
      do Utnemele2=1,Nl2

      close(UNIT=ind(Uelement,Uelement2,Utnemele,Utnemele2,code))

      end do
      end do
      end do
      end do
    end subroutine close_files

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


    ! testing and debugging routine
    subroutine rates_to_evops()
        integer(i4b) :: tind
        real(dp)     :: time

        real(dp), parameter :: &
                    AA = 1.0/200,    &
                    BB = 1.0/300,    &
                    CC = 1.0/100

        STEPS = 5000
        call init_lindblad_fit()

        !write(*,*) 'READING EXTERNAL EVOPS'
        !call flush()
        !call open_files('r')
        !call read_evops()
        !call close_files('r')
        timeStep = 0.65

        Evops = 0.0_dp

        do tind = 1, STEPS

        time = timeStep*(tind-1)

        Evops(1,1,1,1,tind) = (BB + AA*exp((-AA - BB)*time))/(AA + BB)
        Evops(2,2,1,1,tind) = -((AA*(-1 + exp((-AA - BB)*time)))/(AA + BB))

        Evops(1,1,2,2,tind) = -((BB*(-1 + exp((-AA - BB)*time)))/(AA + BB))
        Evops(2,2,2,2,tind) = (AA + BB*exp((-AA - BB)*time))/(AA + BB)

        Evops(1,2,1,2,tind) = exp(-CC*time)
        Evops(2,1,2,1,tind) = exp(-CC*time)

        end do

        write(*,*) 'OUTPUTTING EVOPS'
        call flush()
        call open_files('w')
        call write_evops()
        call close_files('w')

    end subroutine rates_to_evops

end module lindblad_fit_module

