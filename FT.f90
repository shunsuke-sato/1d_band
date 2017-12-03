program main
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0,1d0)
  integer,parameter :: nt = 420000, nw = 400
  real(8),parameter :: fs = 0.02418d0
  real(8),parameter :: ev = 27.2114d0
  real(8),parameter :: ww_f = 100d0/ev, dw = ww_f/nw
  real(8) :: tt(0:nt),jt(0:nt),window(0:nt),dt,xx
  real(8) :: tpulse
  integer :: it,iw
  real(8) :: ww
  complex(8) :: zjw
  real(8) :: f1,f2,f3

  tpulse = 96.66d0/fs
  
  do it = 0,nt
     read(*,*)tt(it),f1,jt(it)
  end do
  dt = tt(1) - tt(0)

  window = 0d0
  do it = 0,nt
!     if(tt(it) < tpulse)then
!        window(it) = 1d0
!     else
!        xx = (tt(it)-tpulse)/(tt(nt)-tpulse)
!        window(it) = 1d0 - 3d0*xx**2 + 2d0*xx**3
!     end if
     if(tt(it) < tpulse)then
        window(it) = sin(pi*tt(it)/tpulse)**4
     end if
  end do


  do iw = 0,nw
     ww = dw*iw
     zjw = 0d0
     do it = 0,nt
        zjw = zjw + exp(zi*ww*tt(it))*window(it)*jt(it)
     end do
     zjw = zjw*dt

     write(*,"(999e26.16e3)")ww,abs(ww*zjw)**2

  end do
  
end program main
  
