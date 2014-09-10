program LindbladFit
  use lindblad_fit_module
  use std_types
  use numer_matrix
  use helpers

  implicit none

  complex(dpc), dimension(:,:), allocatable :: AA,BB,CC
  real(dp), dimension(:), allocatable :: DD
  integer(i4b) :: M,N

  integer :: num_args, ix
  character(len=12), dimension(:), allocatable :: args

  num_args = command_argument_count()
  allocate(args(num_args))

  do ix = 1, num_args
      call get_command_argument(ix,args(ix))
  end do

  if(num_args /= 2) then
    call rates_to_evops()
    !call only_convert_to_exciton()
    write(*,*) 'wrong number of arguments'
    stop
  else
    read( args(1), '(i10)' ) Nbasis
    read( args(2), '(i10)' ) STEPS
  end if

  call do_lindblad_fit_work()
  stop

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! debug code follows                         !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(num_args /= 2) then
    write(*,*) 'wrong number of arguments'
    stop
  else
    read( args(1), '(i10)' ) M
    read( args(2), '(i10)' ) N
  end if

  allocate(AA(M,N))
  allocate(BB(M,M))
  allocate(CC(N,N))
  allocate(DD(min(M,N)))

  AA = 1.0_dp
  !BB = 0.0_dp
  !BB(1,1) = 1.0_dp
  !BB(2,2) = 1.0_dp

  call svd(AA,BB,DD,CC)
  write(*,*) real(BB)
  write(*,*)
  write(*,*) real(DD)
  write(*,*)
  write(*,*) real(CC)
  stop


end program LindbladFit
