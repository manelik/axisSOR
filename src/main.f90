

program main

  implicit none

  real*8, dimension(:,:), allocatable :: a,b,c,d,e,f,u
  real*8, dimension(:,:), allocatable :: x,y

  real*8 :: Lx,Ly

  real*8 :: delta_x,delta_y,rojac, smallpi
  integer :: Jx,Jy

  integer :: m,k

  integer :: nx,ny

  smallpi=acos(-1.d0)

  print*, 'dimension x'
  read(*,*) Jx
  print*, 'dimension Y'
  read(*,*) Jy
  print*, 'intervalx'
  read(*,*) delta_x
  print*, 'intervaly'
  read(*,*) delta_y
  print*, 'Blobs in x direction'
  read(*,*) nx
  print*, 'Blobs in y direction'
  read(*,*) ny

  allocate( a(Jx,Jy),b(Jx,Jy),c(Jx,Jy),d(Jx,Jy),e(Jx,Jy),&
       &f(Jx,Jy),u(Jx,Jy),x(Jx,Jy),y(Jx,Jy))
  
  do k=1,jy
     do m=1,jx
        x(m,k)= delta_x*(-.5d0+dble(m))
        y(m,k)= delta_y*(-.5d0+dble(k))
     end do
  end do

  Lx = x(Jx,Jy)
  Ly = y(Jx,Jy)

!  nx = 4
!  ny = 2
  u=0.0

!  u = 1.0
  u(1:Jx,1)=0.d0!log(x(1:J,1)**2+y(1:J,1)**2)
  u(1:Jx,Jy)=0.d0!log(x(1:J,J)**2+y(1:J,J)**2)
  u(1,1:Jy)=sin(ny*smallpi/Ly*y(1,1:Jy))!log(x(1,1:J)**2+y(1,1:J)**2)
!  u(Jx,1:Jy)=y(Jx,1:Jy)/Ly!log(x(J,1:J)**2+y(J,1:J)**2)

  a=delta_y/delta_x+0.5d0*delta_y/x
  b=delta_y/delta_x-0.5d0*delta_y/x
  c=delta_x/delta_y
  d=delta_x/delta_y
  e=-2.d0*(delta_y/delta_x+delta_x/delta_y)

  f= -delta_x*delta_y* (( ((1.d0+2.d0*dble(nx))/Lx)**2+(dble(ny)/Ly)**2)*smallpi**2*&
       &cos((1.d0+2.d0*dble(nx))*smallpi/Lx*x) +&
       &  ((1.d0+2.d0*dble(nx))/Lx)*smallpi*sin((1.d0+2.d0*dble(nx))*smallpi/Lx*x) )*sin(dble(ny)*smallpi/Ly*y)
!  u= sin(nx*smallpi/Lx*x)*sin(ny*smallpi/Ly*y)


  rojac=( delta_y/delta_x*cos(smallpi/dble(Jx)) &
       &+ delta_x/delta_y*cos(smallpi/dble(Jy)))/&
       &(delta_y/delta_x+delta_x/delta_y)   !1.d0!1.d0-2.d0*3.14159d0/J

  call sor(a,b,c,d,e,f,u,rojac,Jx,Jy)

  open(unit=666,file='out.dat',status='replace')
  
  do k=1,jy
     do m=1,jx
        write(666,*) x(m,k),y(m,k),u(m,k)
     end do
  end do

  close(666)

end program main
