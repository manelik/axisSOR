

SUBROUTINE sor(a,b,c,d,e,f,u,rjac,nx,ny,ghost,eqsym,rsym,zsym,delta_x,delta_y,x,y,rob_bc)

! Include modules


  IMPLICIT NONE
! Taken from Numerical Recipies and adapted to actually work.
! ie general mxn resolutions, not assuming n,m are odd
  
  INTEGER, intent(in) :: nx,ny,ghost,rsym,zsym
  LOGICAL, intent(in) :: eqsym,rob_bc
  REAL*8, DIMENSION(1-ghost:nx,1-ghost:ny), INTENT(IN) :: a,b,c,d,e,f,x,y
  REAL*8, DIMENSION(1-ghost:nx,1-ghost:ny), INTENT(INOUT) :: u
  REAL*8, INTENT(IN) :: rjac,delta_x,delta_y
  INTEGER, PARAMETER :: MAXITS=1000000
  REAL*8, PARAMETER :: EPS=1.0d-5
  !  Successive overrelaxation solution of equation (19.5.25) with Chebyshev acceleration. a, b,
  !  c, d, e, and f are input as the coefficients of the equation, each dimensioned to the grid
  !  size J ×J. u is input as the initial guess to the solution, usually zero, and returns with the
  !  final value. rjac is input as the spectral radius of the Jacobi iteration, or an estimate of
  !  it. Double precision is a good idea for J bigger than about 25.
  REAL*8, DIMENSION(1-ghost:nx,1-ghost:ny) :: resid
  INTEGER :: jxmax,jxm1,jxm2,jxm3,n,k
  INTEGER :: jymax,jym1,jym2,jym3
  REAL*8 :: anorm,anormf,omega
! for the moment no validation of array sizes is done since we dont use nr libs
! I assume the user will send equal sized arrays 
  jxmax=nx
  jxm1=jxmax-1
  jxm2=jxmax-2
  jxm3=jxmax-3
  jymax=ny
  jym1=jymax-1
  jym2=jymax-2
  jym3=jymax-3
  anormf=sum(abs(f(1:jxm1,1:jym1)))
  if (anormf==0.d0) anormf= 1.d-5 !default tolerance in case of no f input
  !  Compute initial norm of residual and terminate iteration when norm has been reduced by a
  !  factor EPS. This computation assumes initial u is zero.
  omega=1.0
  do n=1,MAXITS
!     First do the even-even and odd-odd squares of the grid, i.e., the red squares of the
!     checkerboard:
     resid( 1:jxm1:2,1:jym1:2)=&
          a(1:jxm1:2,1:jym1:2)*u(2:jxmax:2,1: jym1:2)+&
          b(1:jxm1:2,1:jym1:2)*u(0: jxm2:2,1: jym1:2)+&
          c(1:jxm1:2,1:jym1:2)*u(1: jxm1:2,2:jymax:2)+&
          d(1:jxm1:2,1:jym1:2)*u(1: jxm1:2,0: jym2:2)+&
          e(1:jxm1:2,1:jym1:2)*u(1: jxm1:2,1: jym1:2)-f(1:jxm1:2,1:jym1:2)
     u(1:jxm1:2,1:jym1:2)=u(1:jxm1:2,1:jym1:2)-omega*&
          resid(1:jxm1:2,1:jym1:2)/e(1:jxm1:2,1:jym1:2)
     resid( 2:jxm1:2,2:jym1:2)=&
          a(2:jxm1:2,2:jym1:2)*u(3:jxmax:2,2: jym1:2)+&
          b(2:jxm1:2,2:jym1:2)*u(1: jxm2:2,2: jym1:2)+&
          c(2:jxm1:2,2:jym1:2)*u(2: jxm1:2,3:jymax:2)+&
          d(2:jxm1:2,2:jym1:2)*u(2: jxm1:2,1: jym2:2)+&
          e(2:jxm1:2,2:jym1:2)*u(2: jxm1:2,2: jym1:2)-f(2:jxm1:2,2:jym1:2)
     u(2:jxm1:2,2:jym1:2)=u(2:jxm1:2,2:jym1:2)-omega*&
          resid(2:jxm1:2,2:jym1:2)/e(2:jxm1:2,2:jym1:2)
     omega=merge(1.0d0/(1.0d0-0.5d0*rjac**2), &
          1.0d0/(1.0d0-0.25d0*rjac**2*omega), n == 1)
     !Now do even-odd and odd-even squares of the grid, i.e., the black squares of the checkerboard:
     resid( 2:jxm1:2,1:jym1:2)=&
          a(2:jxm1:2,1:jym1:2)*u(3:jxmax:2,1: jym1:2)+&
          b(2:jxm1:2,1:jym1:2)*u(1: jxm2:2,1: jym1:2)+&
          c(2:jxm1:2,1:jym1:2)*u(2: jxm1:2,2:jymax:2)+&
          d(2:jxm1:2,1:jym1:2)*u(2: jxm1:2,0: jym2:2)+&
          e(2:jxm1:2,1:jym1:2)*u(2: jxm1:2,1: jym1:2)-f(2:jxm1:2,1:jym1:2)
     u(2:jxm1:2,1:jym1:2)=u(2:jxm1:2,1:jym1:2)-omega*&
          resid(2:jxm1:2,1:jym1:2)/e(2:jxm1:2,1:jym1:2)
     resid( 1:jxm1:2,2:jym1:2)=&
          a(1:jxm1:2,2:jym1:2)*u(2:jxmax:2,2: jym1:2)+&
          b(1:jxm1:2,2:jym1:2)*u(0: jxm2:2,2: jym1:2)+&
          c(1:jxm1:2,2:jym1:2)*u(1: jxm1:2,3:jymax:2)+&
          d(1:jxm1:2,2:jym1:2)*u(1: jxm1:2,1: jym2:2)+&
          e(1:jxm1:2,2:jym1:2)*u(1: jxm1:2,2: jym1:2)-f(1:jxm1:2,2:jym1:2)
     u(1:jxm1:2,2:jym1:2)=u(1:jxm1:2,2:jym1:2)-omega*&
          resid(1:jxm1:2,2:jym1:2)/e(1:jxm1:2,2:jym1:2)

     !apply boundary conditions, for the moment just standard robin
     ! u_inf=1, decays as 1/R
     if (rob_bc)then
        !do whole boundary, corners will be fixed later
        ! first r boundary
!        print*, 'applying r bc'
        u(jxmax,:)= (2.d0*delta_x*x(jxmax,:)/(x(jxmax,:)**2+y(jxmax,:)**2) &
             &- u(jxmax-2,:)+ 4.d0*u(jxmax-1,:))&
             /(2.d0*delta_x*x(jxmax,:)/(x(jxmax,:)**2+y(jxmax,:)**2)+3.d0)
        
        !upper boundary
!        print*, 'applying upper bc'
        u(:,jymax)= (2.d0*delta_y*y(:,jymax)/(x(:,jymax)**2+y(:,jymax)**2) &
             &- u(:,jymax-2)+ 4.d0*u(:,jymax-1))&
             /(2.d0*delta_y*y(:,jymax)/(x(:,jymax)**2+y(:,jymax)**2)+3.d0)
        !lower boundary
!        print*, 'applying lower bc'
        u(:,0)= (2.d0*delta_y*y(:,0)/(x(:,0)**2+y(:,0)**2) &
             &+ u(:,2)- 4.d0*u(:,1))&
             /(2.d0*delta_y*y(:,0)/(x(:,0)**2+y(:,0)**2)-3.d0)

     end if

     !apply symmetries

     do k = 1,ghost
        u(1-k,:)=dble(rsym)*u(k,:)
        if (eqsym) u(:,1-k)=dble(zsym)*u(:,k)
     end do

     ! sync when it has to be synched

     omega=1.0d0/(1.0d0-0.25d0*rjac**2*omega)
     anorm=sum(abs(resid(1:jxm1,1:jym1)))
     if (anorm < EPS*anormf) exit
     if (mod(n,1000)==0) print *, n,omega,anorm
  end do
  if (n > MAXITS) then !call nrerror(’MAXITS exceeded in sor’)
     print*, 'MAXITS exceeded in sor'
  else
     print*, 'SOR acheved the desired tolerance after ',n, 'iterations'
  end if

END SUBROUTINE sor
