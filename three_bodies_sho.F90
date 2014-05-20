 module utilities_sho
 implicit none
 contains

 Subroutine Fix_cor(xf,yf,pxf,pyf,x,y,vx,vy,px,py,xp,yp,pxp,pyp,t,ome,eps)
 Implicit none
 real (kind=8) :: c, s
 real (kind=8) :: xc, yc, pxc, pyc
 real (kind=8), intent(in) :: xf, yf, pxf, pyf, t, ome, eps
 real (kind=8), intent(out) :: x,y,vx,vy,px,py,xp,yp,pxp,pyp
 c = COS(ome*t)
 s = SIN(ome*t)
 x = c*xf + s*yf
 y = -s*xf + c*yf
 vx = c*pxf + s*pyf - s*xf + c*yf
 vy = -s*pxf + c*pyf - c*xf - s*yf
 end Subroutine Fix_cor
 ! --------------------- !
 Subroutine Cor_fix(x,y,vx,vy,xf,yf,pxf,pyf,t,ome)
 Implicit none
 real (kind=8) :: c, s
 real (kind=8), intent(in) :: x,y,vx,vy,t,ome
 real (kind=8), intent(out) :: xf,yf,pxf,pyf
 c = COS(ome*t)
 s = SIN(ome*t)
 xf = c*x - s*y
 yf = s*x + c*y
 pxf = c*vx - s*vy - s*x - c*y
 pyf = s*vx + c*vy + c*x - s*y
 end Subroutine Cor_fix
 ! --------------------- !
 real*8 function h(x,y,tau,px,py,ptau,ome,eps)
 Implicit none
 real (kind=8), intent(in) :: x,y,tau,px,py,ptau,ome,eps
 real (kind=8) ::  c, s, r1, r2

 c=cos(ome*tau)
 s=sin(ome*tau)

 r1 = sqrt((x+eps*c)**2+(y+eps*s)**2)
 r2 = sqrt((x-(1.d0-eps)*c)**2+(y-s*(1.d0-eps))**2)

 h = px*px/2.d0+ py*py/2.d0 + ptau
 h = h-(1.d0-eps)/r1-eps/r2

 !h = h-0.5d0*(1.d0-eps+eps**2)         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 return
 end function h
! --------------------- !


 subroutine sym2(x,y,tau,px,py,ptau,dt,ome,eps)
 Implicit none
 real (kind=8), intent(in) :: dt,ome,eps
 real (kind=8), intent(inout) :: x,y,tau,px,py,ptau
 real (kind=8) :: xnew, ynew, taunew, pxnew, pynew, ptaunew, fx, fy, ftau, fxnew, fynew, ftaunew

 call f(x,y,tau,fx,fy,ftau,ome,eps)

 xnew = x + px*dt + fx*dt**2/2.d0
 ynew = y + py*dt + fy*dt**2/2.d0
 taunew = tau + dt

 call f(xnew,ynew,taunew,fxnew,fynew,ftaunew,ome,eps)

 pxnew = px + dt*(fx+fxnew)/2.d0
 pynew = py + dt*(fy+fynew)/2.d0
 ptaunew = ptau + dt*(ftau+ftaunew)/2.d0

 x = xnew
 y = ynew
 tau = taunew
 px = pxnew
 py = pynew
 ptau = ptaunew

 contains
 subroutine f(x,y,tau,fx,fy,ftau,ome,eps)
 real (kind=8) :: c, s, r1, r2
 real (kind=8), intent(in) :: x, y, tau, ome, eps
 real (kind=8), intent(out) :: fx, fy, ftau

 c = cos(ome*tau)
 s = sin(ome*tau)

 r1 = sqrt((x+eps*c)**2+(y+eps*s)**2)
 r2 = sqrt((x-(1.d0-eps)*c)**2+(y-s*(1.d0-eps))**2)

 fx = -(1-eps)*(x+eps*c)/r1**3-eps*(x-c*(1-eps))/r2**3
 fy = -(1-eps)*(y+eps*s)/r1**3-eps*(y-s*(1-eps))/r2**3

 ftau = -(1.d0-eps)*ome*eps*(-s*(x+eps*c)+c*(y+eps*s))/r1**3
 ftau = ftau-eps*(1.d0-eps)*ome*(s*(x-(1.d0-eps)*c)-c*(y-(1.d0-eps)*s))/r2**3

 end subroutine f
 end subroutine sym2


 subroutine sym4(x,y,tau,px,py,ptau,dt,ome,eps)
 Implicit none
 real (kind=8), intent(inout) :: x,y,tau,px,py,ptau
 real (kind=8), intent(in) :: ome,eps,dt
 real (kind=8) ::  sq2, alpha, beta, dt1, dt2

 sq2 = 2.d0**(1.d0/3.d0)
 alpha = 1.d0/(2.d0-sq2)
 beta = sq2/(2.d0-sq2)
 dt1 = dt*alpha
 dt2 = -dt*beta
 call sym2(x,y,tau,px,py,ptau,dt1,ome,eps)
 call sym2(x,y,tau,px,py,ptau,dt2,ome,eps)
 call sym2(x,y,tau,px,py,ptau,dt1,ome,eps)

 end subroutine sym4

 end module utilities_sho



 Program threebodies_sho
 use utilities_sho
 Implicit none
 integer, parameter :: N_step = 1000 ! passi per periodo
 integer, parameter :: N_periods = 2 ! periodi

 logical, parameter :: chaos = .false.

 integer :: k, k_max, m_per
 real (kind=8), parameter :: pi = 3.141592653589793
 real (kind=8), parameter :: ome = 1.d0
 real (kind=8), parameter :: eps = 9.54d-4 ! 0.002 caotico

 real (kind=8) :: Jac, T_per, xcor0, ycor0, tau0, vxcor0, vycor0, ptau0, t, Hrc
 real (kind=8) :: xf0, yf0, pxf0, pyf0, x, y, vx, vy, px, py, xp, yp, pxp, pyp
 real (kind=8) :: xf, yf, tau, pxf, pyf, ptau, A_step, dt, E0, E0r, E0r_quad, errh, Ermed, E, Er

 T_per=2.d0*pi/ome


 if(chaos .eqv. .true.) then
  xcor0 = 0.56d0
 else
  xcor0 = 0.55d0
 endif
 ycor0 = 0.d0
 tau0 = 0.d0 ! sistema corotante
 vxcor0 = 0.d0 ! velocità
 vycor0 = 0.d0 ! iniziali
 ptau0 = 1.d0 ! sistema corotante
 t = 0.d0
 Jac = 3.07d0
 vycor0 = sqrt(2.d0*(1.d0-eps)/abs(xcor0+eps)+(2.d0*eps)/abs(xcor0-1.d0+eps)+xcor0**2.d0+ycor0**2.d0-Jac-vxcor0**2.d0)
 
 print*, '*** vy0 eps ***'
 print*, vycor0, vxcor0, ycor0, xcor0

 xf0 = 0.d0
 yf0 = 0.d0
 pxf0 = 0.d0
 pyf0 = 0.d0

 call Cor_fix(xcor0,ycor0,vxcor0,vycor0,xf0,yf0,pxf0,pyf0,t,ome)
 print*, xf0,yf0,pxf0,pyf0

 xf = xf0
 yf = yf0
 tau = tau0
 pxf = pxf0
 pyf = pyf0
 ptau = ptau0
 WRITE(7,*) xcor0,vxcor0,ycor0,vycor0

 A_step = Dfloat(N_step)
 dt = T_per/A_step
 k_max = N_periods*N_step

 E0 = H(xf,yf,tau,pxf,pyf,ptau,ome,eps) ! H sis.fisso

 t = 0.d0
 errh = 0.d0
 x = 0.d0
 y = 0.d0
 vx = 0.d0
 vy = 0.d0
 px = 0.d0
 py = 0.d0
 xp = 0.d0
 yp = 0.d0
 pxp = 0.d0
 pyp = 0.d0

 do k=1,k_max
  call sym4(xf,yf,tau,pxf,pyf,ptau,dt,ome,eps)
  t = t+dt
  E = H(xf,yf,tau,pxf,pyf,ptau,ome,eps)
  call Fix_cor(xf,yf,pxf,pyf, x,y,vx,vy,px,py, xp,yp,pxp,pyp,t,ome,eps)
  errh = abs(E-E0)
  WRITE(7,*) x,vx,y,vy
  WRITE(12,*) t, errh
  write(13,*) xf,yf,pxf,pyf
 enddo


 end program threebodies_sho

