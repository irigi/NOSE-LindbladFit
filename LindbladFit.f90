program LindbladFit
  use lindblad_fit_module
  use std_types
  use numer_matrix
  use helpers

  implicit none

  real(dp), dimension(:,:), allocatable :: AA,BB,CC
  real(dp), dimension(:), allocatable :: DD
  integer(i4b) :: M,N

  integer :: num_args, ix
  character(len=12), dimension(:), allocatable :: args

  num_args = command_argument_count()
  allocate(args(num_args))

  do ix = 1, num_args
      call get_command_argument(ix,args(ix))
  end do

  call do_lindblad_fit_work()
  stop

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
  BB = 0.0_dp
  BB(1,1) = 1.0_dp
  BB(2,2) = 1.0_dp

  call svd(AA,BB,DD,CC)
  stop


end program LindbladFit
