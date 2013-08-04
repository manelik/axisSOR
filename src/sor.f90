

SUBROUTINE sor(a,b,c,d,e,f,u,rjac,nx,ny)
!  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  
  INTEGER, intent(in) :: nx,ny
  REAL*8, DIMENSION(nx,ny), INTENT(IN) :: a,b,c,d,e,f
  REAL*8, DIMENSION(nx,ny), INTENT(INOUT) :: u
  REAL*8, INTENT(IN) :: rjac
  INTEGER, PARAMETER :: MAXITS=1000000
  REAL*8, PARAMETER :: EPS=1.0d-5
  !  Successive overrelaxation solution of equation (19.5.25) with Chebyshev acceleration. a, b,
  !  c, d, e, and f are input as the coefficients of the equation, each dimensioned to the grid
  !  size J ×J. u is input as the initial guess to the solution, usually zero, and returns with the
  !  final value. rjac is input as the spectral radius of the Jacobi iteration, or an estimate of
  !  it. Double precision is a good idea for J bigger than about 25.
  REAL*8, DIMENSION(size(a,1),size(a,1)) :: resid
  INTEGER :: jxmax,jxm1,jxm2,jxm3,n
  INTEGER :: jymax,jym1,jym2,jym3
  REAL*8 :: anorm,anormf,omega
!  jmax=assert_eq((/size(a,1),size(a,2),size(b,1),size(b,2), &
!       size(c,1),size(c,2),size(d,1),size(d,2),size(e,1), &
!       size(e,2),size(f,1),size(f,2),size(u,1),size(u,2)/),’sor’)
!  for the moment no validation of array sizes is done we dont use nr libs
  jxmax=size(a,1)
  jxm1=jxmax-1
  jxm2=jxmax-2
  jxm3=jxmax-3
  jymax=size(a,2)
  jym1=jymax-1
  jym2=jymax-2
  jym3=jymax-3
  anormf=sum(abs(f(2:jxm1,2:jym1)))
  !  Compute initial norm of residual and terminate iteration when norm has been reduced by a
  !  factor EPS. This computation assumes initial u is zero.
  omega=1.0
  do n=1,MAXITS
!     First do the even-even and odd-odd squares of the grid, i.e., the red squares of the
!     checkerboard:
     resid( 2:jxm1:2,2:jym1:2)=&
          a(2:jxm1:2,2:jym1:2)*u(3:jxmax:2,2: jym1:2)+&
          b(2:jxm1:2,2:jym1:2)*u(1: jxm2:2,2: jym1:2)+&
          c(2:jxm1:2,2:jym1:2)*u(2: jxm1:2,3:jymax:2)+&
          d(2:jxm1:2,2:jym1:2)*u(2: jxm1:2,1: jym2:2)+&
          e(2:jxm1:2,2:jym1:2)*u(2: jxm1:2,2: jym1:2)-f(2:jxm1:2,2:jym1:2)
     u(2:jxm1:2,2:jym1:2)=u(2:jxm1:2,2:jym1:2)-omega*&
          resid(2:jxm1:2,2:jym1:2)/e(2:jxm1:2,2:jym1:2)
     resid( 3:jxm1:2,3:jym1:2)=&
          a(3:jxm1:2,3:jym1:2)*u(4:jxmax:2,3: jym1:2)+&
          b(3:jxm1:2,3:jym1:2)*u(2: jxm2:2,3: jym1:2)+&
          c(3:jxm1:2,3:jym1:2)*u(3: jxm1:2,4:jymax:2)+&
          d(3:jxm1:2,3:jym1:2)*u(3: jxm1:2,2: jym2:2)+&
          e(3:jxm1:2,3:jym1:2)*u(3: jxm1:2,3: jym1:2)-f(3:jxm1:2,3:jym1:2)
     u(3:jxm1:2,3:jym1:2)=u(3:jxm1:2,3:jym1:2)-omega*&
          resid(3:jxm1:2,3:jym1:2)/e(3:jxm1:2,3:jym1:2)
     omega=merge(1.0d0/(1.0d0-0.5d0*rjac**2), &
          1.0d0/(1.0d0-0.25d0*rjac**2*omega), n == 1)
     !Now do even-odd and odd-even squares of the grid, i.e., the black squares of the checkerboard:
     resid( 3:jxm1:2,2:jym1:2)=&
          a(3:jxm1:2,2:jym1:2)*u(4:jxmax:2,2: jym1:2)+&
          b(3:jxm1:2,2:jym1:2)*u(2: jxm2:2,2: jym1:2)+&
          c(3:jxm1:2,2:jym1:2)*u(3: jxm1:2,3:jymax:2)+&
          d(3:jxm1:2,2:jym1:2)*u(3: jxm1:2,1: jym2:2)+&
          e(3:jxm1:2,2:jym1:2)*u(3: jxm1:2,2: jym1:2)-f(3:jxm1:2,2:jym1:2)
     u(3:jxm1:2,2:jym1:2)=u(3:jxm1:2,2:jym1:2)-omega*&
          resid(3:jxm1:2,2:jym1:2)/e(3:jxm1:2,2:jym1:2)
     resid( 2:jxm1:2,3:jym1:2)=&
          a(2:jxm1:2,3:jym1:2)*u(3:jxmax:2,3: jym1:2)+&
          b(2:jxm1:2,3:jym1:2)*u(1: jxm2:2,3: jym1:2)+&
          c(2:jxm1:2,3:jym1:2)*u(2: jxm1:2,4:jymax:2)+&
          d(2:jxm1:2,3:jym1:2)*u(2: jxm1:2,2: jym2:2)+&
          e(2:jxm1:2,3:jym1:2)*u(2: jxm1:2,3: jym1:2)-f(2:jxm1:2,3:jym1:2)
     u(2:jxm1:2,3:jym1:2)=u(2:jxm1:2,3:jym1:2)-omega*&
          resid(2:jxm1:2,3:jym1:2)/e(2:jxm1:2,3:jym1:2)
     omega=1.0d0/(1.0d0-0.25d0*rjac**2*omega)
     anorm=sum(abs(resid(2:jxm1,2:jym1)))
     if (anorm < EPS*anormf) exit
     if (mod(n,1000)==0) print *, n,omega,anorm
  end do
  if (n > MAXITS) then !call nrerror(’MAXITS exceeded in sor’)
     print*, 'MAXITS exceeded in sor'
  else
     print*, 'SOR acheved the desired tolerance after ',n, 'iterations'
  end if

END SUBROUTINE sor
