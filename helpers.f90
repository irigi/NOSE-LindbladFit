module helpers
    use std_types

    implicit none

    !integer, dimension(:,:), allocatable :: std_type_kr

    interface trace
      module procedure trace_real
      module procedure trace_cmplx
    end interface

    interface ode_rk4
       module procedure rk4_sp
       module procedure rk4_dp
       module procedure rk4_dpc
       module procedure rk4_dpc_mat
       module procedure rk4_dpc_mat3
    end interface

    contains


    !
    ! matrix trace
    !
    real(dp) function trace_real(A) result(tr)
      real(dp), dimension(:,:), intent(in)    :: A
      integer(i4b)                            :: i

      if(size(A,1) /= size(A,2)) then
          write(*,*) "non-square matrix in trace"
          stop
      end if

      tr = 0.0_dp

      do i=1,size(A,1)
          tr = tr + A(i,i)
      end do
    end function trace_real

    complex(dpc) function trace_cmplx(A) result(tr)
      complex(dpc), dimension(:,:), intent(in)    :: A
      integer(i4b)                                :: i

      if(size(A,1) /= size(A,2)) then
          write(*,*) "non-square matrix in trace"
          stop
      end if

      tr = 0.0_dp

      do i=1,size(A,1)
          tr = tr + A(i,i)
      end do
    end function trace_cmplx

    function factorial (n) result (res)

      implicit none
      integer(i4b), intent (in) :: n
      real(dp) :: res
      integer(i4b) :: i

      !res = product ((/(i, i = 1, n)/))
      res = 1
      do i=2,n
          res  = res*i
      end do

    end function factorial

    !
    ! Dir name manipulation
    !
    function file_join(fdname1,fdname2) result (fdname3)
        character(len=*), intent(in) :: fdname1, fdname2
        character(len=256) :: fdname3

        if (len_trim(fdname1)>0) then
            fdname3 = trim(fdname1)//"/"//trim(fdname2)
        else
            fdname3 = trim(fdname2)
        end if

    end function file_join



    ! double precision complex matrix version
    subroutine rk4_dpc_mat(y,dydx,x,h,yout,derivs)
      use std_types !; use nrutil, only : assert_eq
      implicit none
      complex(DPC), dimension(:,:), intent(IN)  :: y,dydx
      real(DP),                     intent(IN)  :: x,h
      complex(DPC), dimension(:,:), intent(OUT) :: yout
      interface
         subroutine derivs(x,y,dydx)
           use std_types
           implicit none
           real(DP), intent(IN) :: x
           complex(DPC), dimension(:,:), intent(IN)  :: y
           complex(DPC), dimension(:,:), intent(OUT) :: dydx
         end subroutine derivs
      end interface
      integer(I4B) :: ndum
      real(DP) :: h6,hh,xh
      complex(DPC), dimension(size(y,1),size(y,2)) :: dym,dyt,yt
!      ndum=assert_eq(size(y),size(dydx),size(yout),'rk4_dpc')
      hh=h*0.5_dp
      h6=h/6.0_dp
      xh=x+hh
      yt=y+hh*dydx
      call derivs(xh,yt,dyt)
      yt=y+hh*dyt
      call derivs(xh,yt,dym)
      yt=y+h*dym
      dym=dyt+dym
      call derivs(x+h,yt,dyt)
      yout=y+h6*(dydx+dyt+2.0_dp*dym)
    end subroutine rk4_dpc_mat

    ! double precision complex matrix with 3 indices version
    subroutine rk4_dpc_mat3(y,dydx,x,h,yout,derivs)
      use std_types !; use nrutil, only : assert_eq
      implicit none
      complex(DPC), dimension(:,:,:), intent(IN)  :: y,dydx
      real(DP),                     intent(IN)  :: x,h
      complex(DPC), dimension(:,:,:), intent(OUT) :: yout
      interface
         subroutine derivs(x,y,dydx)
           use std_types
           implicit none
           real(DP), intent(IN) :: x
           complex(DPC), dimension(:,:,:), intent(IN)  :: y
           complex(DPC), dimension(:,:,:), intent(OUT) :: dydx
         end subroutine derivs
      end interface
      integer(I4B) :: ndum
      real(DP) :: h6,hh,xh
      complex(DPC), dimension(size(y,1),size(y,2),size(y,3)) :: dym,dyt,yt
!      ndum=assert_eq(size(y),size(dydx),size(yout),'rk4_dpc')
      hh=h*0.5_dp
      h6=h/6.0_dp
      xh=x+hh
      yt=y+hh*dydx
      call derivs(xh,yt,dyt)
      yt=y+hh*dyt
      call derivs(xh,yt,dym)
      yt=y+h*dym
      dym=dyt+dym
      call derivs(x+h,yt,dyt)
      yout=y+h6*(dydx+dyt+2.0_dp*dym)
    end subroutine rk4_dpc_mat3

!     double precision complex version
    subroutine rk4_dpc(y,dydx,x,h,yout,derivs)
      use std_types !; use nrutil, only : assert_eq
      implicit none
      complex(DPC), dimension(:), intent(IN)  :: y,dydx
      real(DP),                   intent(IN)  :: x,h
      complex(DPC), dimension(:), intent(OUT) :: yout
      interface
         subroutine derivs(x,y,dydx)
           use std_types
           implicit none
           real(DP), intent(IN) :: x
           complex(DPC), dimension(:), intent(IN)  :: y
           complex(DPC), dimension(:), intent(OUT) :: dydx
         end subroutine derivs
      end interface
      integer(I4B) :: ndum
      real(DP) :: h6,hh,xh
      complex(DPC), dimension(size(y)) :: dym,dyt,yt
!      ndum=assert_eq(size(y),size(dydx),size(yout),'rk4_dpc')
      hh=h*0.5_sp
      h6=h/6.0_sp
      xh=x+hh
      yt=y+hh*dydx
      call derivs(xh,yt,dyt)
      yt=y+hh*dyt
      call derivs(xh,yt,dym)
      yt=y+h*dym
      dym=dyt+dym
      call derivs(x+h,yt,dyt)
      yout=y+h6*(dydx+dyt+2.0_dp*dym)
    end subroutine rk4_dpc

!     double precision version
    subroutine rk4_dp(y,dydx,x,h,yout,derivs)
      use std_types !; use nrutil, only : assert_eq
      implicit none
      real(DP), dimension(:), intent(IN)  :: y,dydx
      real(DP),               intent(IN)  :: x,h
      real(DP), dimension(:), intent(OUT) :: yout
      interface
         subroutine derivs(x,y,dydx)
           use std_types
           implicit none
           real(DP), intent(IN) :: x
           real(DP), dimension(:), intent(IN)  :: y
           real(DP), dimension(:), intent(OUT) :: dydx
         end subroutine derivs
      end interface
      integer(I4B) :: ndum
      real(DP) :: h6,hh,xh
      real(DP), dimension(size(y)) :: dym,dyt,yt
!      ndum=assert_eq(size(y),size(dydx),size(yout),'rk4_dpc')
      hh=h*0.5_sp
      h6=h/6.0_sp
      xh=x+hh
      yt=y+hh*dydx
      call derivs(xh,yt,dyt)
      yt=y+hh*dyt
      call derivs(xh,yt,dym)
      yt=y+h*dym
      dym=dyt+dym
      call derivs(x+h,yt,dyt)
      yout=y+h6*(dydx+dyt+2.0_dp*dym)
    end subroutine rk4_dp

!     single precision version
    subroutine rk4_sp(y,dydx,x,h,yout,derivs)
      use std_types !; use nrutil, only : assert_eq
      implicit none
      real(SP), dimension(:), intent(IN) :: y,dydx
      real(SP), intent(IN) :: x,h
      real(SP), dimension(:), intent(OUT) :: yout
      interface
         subroutine derivs(x,y,dydx)
          use std_types
           implicit none
           real(SP), intent(IN) :: x
           real(SP), dimension(:), intent(IN) :: y
           real(SP), dimension(:), intent(OUT) :: dydx
         end subroutine derivs
      end interface
      integer(I4B) :: ndum
      real(SP) :: h6,hh,xh
      real(SP), dimension(size(y)) :: dym,dyt,yt
!      ndum=assert_eq(size(y),size(dydx),size(yout),'rk4')
      hh=h*0.5_sp
      h6=h/6.0_sp
      xh=x+hh
      yt=y+hh*dydx
      call derivs(xh,yt,dyt)
      yt=y+hh*dyt
      call derivs(xh,yt,dym)
      yt=y+h*dym
      dym=dyt+dym
      call derivs(x+h,yt,dyt)
      yout=y+h6*(dydx+dyt+2.0_sp*dym)
    end subroutine rk4_sp
end module helpers
