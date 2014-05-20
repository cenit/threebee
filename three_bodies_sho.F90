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

 px = vx-y
 py = vy+x

 xc = 0.5d0 - eps
 yc = sqrt(3.d0)/2.d0

 pxc = -yc
 pyc = xc

 xp = x-xc
 yp = y-yc

 pxp = px-pxc
 pyp = py-pyc

 end Subroutine Fix_cor



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
 !write(*,*) '++++++++++++++++'
 !write(*,'(8g12.4)') t,c,s,pxf,pyf

 end Subroutine Cor_fix







 REAL*8 FUNCTION Hr__(x,y,vx,vy,t,eps)
 Implicit none

 real (kind=8), intent(in) :: x,y,vx,vy,t,eps
 real (kind=8) ::  xc, yc, r1, r2, T_cin, V_eff

 !-----------------------------
 !
 ! Energia nel sistema corotante H= (vx^2+vy^2)/2 + V_eff
 !
 ! Input coordnate e velocità sistema comoving x,y,vx,vy
 !
 !-----------------------------

 xc = 0.5d0 - eps
 yc = sqrt(3.d0)/2.d0

 r1 = sqrt((x+eps)**2 + y**2)
 r2 = sqrt((x-1.d0+eps)**2 + y**2)

 T_cin = (vx*vx)/2.d0 + (vy*vy)/2.d0
 V_eff = -(1.d0-eps)/r1 - eps/r2 - 0.5d0*(x**2+y**2)

 !---------------------
 ! Si rende V_eff nullo all'equilibrio
 !---------------------
 V_eff = V_eff + 1.d0 + 0.5d0 * (xc**2+yc**2)
 Hr__ = T_cin + V_eff

 write(15,'(8g12.4)') t,x,y,vx,vy,T_cin,V_eff,Hr__

 return
 end FUNCTION Hr__


 REAL*8 FUNCTION Hr(xf,yf,pxf,pyf,t,ome,eps)
 Implicit none
 real (kind=8), intent(in) :: xf,yf,pxf,pyf,t,ome,eps
 real (kind=8) ::  xc, yc, c, s, x, y, vx, vy, r1, r2, T_cin, V_eff

 !-------------
 ! Energia nel sistema corotante H= (vx^2+vy^2)/2 + V_eff
 !
 ! Input coordinate e velocità sistema fisso: xf,yf,pxf,pyf
 !-------------

 xc = 0.5d0-eps
 yc = sqrt(3.d0)/2.d0

 c = COS(OME*T)
 s = SIN(OME*T)

 x = c*xf + s*yf
 y = -s*xf + c*yf

 vx = c*pxf + s*pyf - s*xf + c*yf
 vy = -s*pxf + c*pyf - c*xf - s*yf

 r1 = SQRT((x+eps)**2 + y**2)
 r2 = SQRT((x-1.d0+eps)**2 + y**2)

 T_cin = (vx*vx)/2.d0 + (vy*vy)/2.d0

 V_eff = -(1.d0-eps)/r1 - eps/r2 - 0.5d0*(x**2+y**2)

 !---------------------
 ! Si rende V_eff nullo all'equilibrio
 !---------------------

 V_eff = V_eff + 1.d0 + 0.5d0 * (xc**2+yc**2)
 Hr = T_cin + V_eff

 write(15,*) ' '
 write(15,'(8g12.4)') t, x,y,vx,vy, T_cin, V_eff, Hr

 return
 end FUNCTION Hr




 REAL*8 function Hr_quad__(xp,yp,vx,vy,t,eps)
 Implicit none
 real (kind=8), intent(in) :: xp,yp,vx,vy,t,eps
 real (kind=8) ::  alfa, T_cin, V_eff

 !-----------------------
 ! Energia nel sistema comoving traslato approx. quadratica
 !
 ! coordinate traslate xp,yp,vx,vy
 !
 ! Hr_quad= (vx^2+vy^2)/2 + V_eff(x,y)-Veff(xc,yc)
 !-----------------------

 alfa = 3.d0 * sqrt(3.d0)/4.d0 * (1.d0-2.d0*eps)
 T_cin = (vx*vx)/2.d0 + (vy*vy)/2.D0
 V_eff = -3.d0/8.d0 * xp**2 - 9.d0/8.d0 * yp**2 - alfa*xp*yp
 Hr_quad__ = T_cin + V_eff

 write(15,*) '---- quad ---'
 write(15,'(8g12.4)') t, xp,yp, vx,vy , T_cin,V_eff , Hr_quad__
 write(15,*) ' '
 return
 end function Hr_quad__





 REAL*8 FUNCTION Hr_quad(xf,yf,pxf,pyf,t,ome,eps)
 Implicit none
 real (kind=8), intent(in) :: xf,yf,pxf,pyf,t,ome,eps
 real (kind=8) :: xc, yc, c, s, x, y, vx, vy, xp, yp, alfa, T_cin, V_eff

 xc = 0.5d0-eps
 yc = sqrt(3.d0)/2.d0

 c = COS(ome*t)
 s = SIN(ome*t)

 x = c*xf + s*yf
 y = -s*xf + c*yf
 vx = c*pxf + s*pyf - s*xf + c*yf
 vy = -s*pxf + c*pyf - c*xf - s*yf

 !-------------
 ! Traslazione posizione di equilibrio
 !-------------
 xp = x-xc
 yp = y-yc

 alfa = 3.d0*sqrt(3.d0)/4.d0 * (1.d0-2.d0*eps)
 T_cin = (vx*vx)/2.d0 + (vy*vy)/2.d0
 V_eff = -3.d0/8.d0 * xp**2 - 9.d0/8.d0 * yp**2 - alfa*xp*yp
 Hr_quad = T_cin + V_eff

 write(15,*) ' +++ '
 write(15,'(8g12.4)') t,xp,yp,vx,vy,T_cin,V_eff,Hr_quad

 return
 end FUNCTION Hr_quad





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

 h = h-0.5d0*(1.d0-eps+eps**2)

 return
 end function h





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
 integer, parameter :: N_periods = 2000 ! periodi

 logical, parameter :: chaos = .false.

 integer :: k, k_max, m_per
 real (kind=8), parameter :: pi = 3.141592653589793
 real (kind=8), parameter :: ome = 1.d0
 real (kind=8), parameter :: eps = 9.54d-4 ! 0.002 caotico

 real (kind=8) :: Jac, T_per, xc, yc, x_Sole, y_Sole, x_Giove, y_Giove, x0, y0, tau0, vx0, vy0, ptau0, t, Hrc, vvy0, yold
 real (kind=8) :: xf0, yf0, pxf0, pyf0, x, y, vx, vy, px, py, xp, yp, pxp, pyp
 real (kind=8) :: xf, yf, tau, pxf, pyf, ptau, A_step, dt, E0, E0r, E0r_quad, errh, Ermed, E, Er
 real (kind=8) :: Er__, errh_r, errh_r__, Er_quad, Er_quad__, er_q, er_q__, dist0, dist, erm


 !--------------------------------------------------------------
 ! Sistema fisso coordinate
 ! --------------------------
 !
 ! xf, yf, pxf, pyf Nota che vxf=pxf, vyf=pyf
 !
 !
 ! Sistema corontante verso orario
 ! -------------------------------
 !
 ! x, y, vx, vy momenti px= vx+y, py=vy-x
 !
 ! Posizione equilibrio sistema corotamte
 !
 ! xc= 0.5 - mu yc=sqrt(3)/2 vxc=vyc=0 pxc=yc, pyc=-xc
 !
 ! Coordinate nel sistema corotante traslato nel punto fisso
 ! ---------------------------------------------------------
 !
 ! xp=x-xc, yc=y-yc, vxp=vx, vyp=vy, pxp=px-yc, pyp=py+xc
 !
 ! sotto poniamo eps=mu= m2/(m1+m2)
 !
 ! Variabli scalate: distanza 1(sole) 2 (giove) r21=1, tempo
 !
 ! Periodo rotazione attorno a centro di massa T=1
 !
 ! Integrazione sistema fisso con simplettico ordine 4 per hamiltoniano
 !
 ! Hex= H(xf,yf,pxf,pyf,tau) + p_tau
 !
 ! dove H=(pxf^2+pyf^2)/2-(1-mu)/r1-mu/r2
 !
 !-----------------------------------------------------------


 open(1,file='fort.1')
 open(2,file='fort.2')
 open(3,file='fort.3')
 open(4,file='fort.4')
 open(7,file='fort.7')
 open(15,file='fort.15')
 open(20,file='fort.20')
 open(21,file='fort.21')


 T_per=2.d0*pi/ome


 !Tra 0.54 e 0.3 regione interessante da esplorare
 !
 !Condizioni iniziali sistema corotante verso antiorario

 xc = 0.5d0 - eps ! punto equilibio
 yc = sqrt(3.d0)/2.d0 ! di Lagrange
 x_Sole = -eps
 y_Sole = 0.d0
 x_Giove = 1.d0-eps
 y_Giove = 0.d0

 Hrc = -1.d0 - 0.5d0 * (xc**2 + yc**2)
 if(chaos .eqv. .true.) then
  x0 = 0.56d0
 else
  x0 = 0.55d0
 endif
 y0 = 0.d0
 tau0 = 0.d0 ! sistema coritante
 vx0 = 0.d0 ! velocità
 vy0 = 0.d0 ! iniziali
 ptau0 = 1.d0 ! sistema corotante
 t = 0.d0
 Jac = 3.07d0

 vvy0 = x0**2 + 2.d0*(1.d0-eps)/abs(x0+eps) + 2.d0*eps/abs(x0-1.d0+eps) - Jac
 vvy0 = vvy0 - vx0**2


 write(*,*) vvy0
 if(vvy0.lt.0) stop
 vy0 = sqrt(vvy0)
 write(*,*) '*** vy0 eps ***', vy0, eps, x0, vvy0


 write(*,*) x0,y0,vx0,vy0
 yold = y0

 write(7,'(2I8,8g13.5)') 0,0,0.d0,x0,vx0,y0,vy0,0.d0,0.d0
 write(*,'(6g15.6)') x0,y0,vx0,vy0

 xf0 = 0.d0
 yf0 = 0.d0
 pxf0 = 0.d0
 pyf0 = 0.d0

 call Cor_fix(x0,y0,vx0,vy0,xf0,yf0,pxf0,pyf0,t,ome)

 write(*,'(6g15.6)') xf0,yf0,pxf0,pyf0

 xf = xf0
 yf = yf0
 tau = tau0

 pxf = pxf0
 pyf = pyf0
 ptau = ptau0

 A_step = Dfloat(N_step)
 dt = T_per/A_step
 k_max = N_periods*N_step

 E0 = H(xf,yf,tau,pxf,pyf,ptau,ome,eps) ! H sis.fisso
 E0r = Hr(xf,yf,pxf,pyf,t,ome,eps) ! H sis corot.
 E0r_quad = Hr_quad(xf,yf,pxf,pyf,t,ome,eps) ! H cor. app. quad.

 t = 0.d0
 errh = 0.d0
 Ermed = 0.d0
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
  ermed = ermed+errh

  if(mod(k,N_step).eq.0) then !------->
   m_per = int(k/Dfloat(N_step))

   Er = Hr(xf,yf,pxf,pyf,t,ome,eps)
   Er__ = Hr__(x,y,vx,vy,t,eps)
   errh_r = abs(Er-E0r)
   errh_r__ = abs(Er__-E0r)

   Er_quad = Hr_quad(xf,yf,pxf,pyf,t,ome,eps)
   Er_quad__ = Hr_quad__(xp,yp,vx, vy,t,eps)
   er_q = abs(Er_quad-E0r_quad)
   er_q__ = abs(Er_quad__-E0r_quad)


   write(21,*) ' '
   write(21,'(I5,8g13.5)') m_per,E, errh
   write(21,'(I5,8g13.5)') m_per,Er,Er__,errh_r,errh_r__
   write(21,'(I5,8g13.5)') m_per,Er_quad,Er_quad__,er_q,er_q__


   write(20,'(I5,8g13.5)') m_per,errh,errh_r__,er_q__,errh/abs(E0),errh_r__/abs(E0r),er_q__/abs(E0r_quad)

   dist0 = sqrt((x-x0)**2 + (vx-vx0)**2)
   write(1,'(I5,8g14.5)')m_per,xf,pxf,yf,pyf,errh
   write(2,'(I5,4g13.5,I10)')m_per,xp,vx, yp,vy,k
   write(3,'(I5,8g13.5)')m_per,xp,pxp,yp,pyp
   write(7,'(2I8,8g13.5)')m_per,k,t,x,vx,y,vy,dist0,errh

   if(mod(m_per,100).eq.0) write(*,*) m_per
  endif ! <------------------

  !-------------------
  ! Sezione di Poincaré
  !-------------------
  ! Calcolo sezione di Poincaré anziché ad ogni
  ! periodo. La stessa cosa si può fare in
  ! nella subroutine Usr-Pap
  ! Introrre un flagper fare cond. iniziali
  ! orbite ed oerbite con rumore tutte calcolate
  ! sulla sezione di Poincaré anziché sulla
  ! mappa in un periodo
  !
  !-----------------------
  if(y*yold.lt.0.and.vy.gt.0) then

   write(4,'(I8,8g13.5)')k,x,y,vx,vy
   dist = sqrt(x**2+vx**2)
   if(dist.ge.100) then
    write(*,*) ' **** instability***',k,k/1000.d0, dist
    stop
   endif
  endif
  yold = y

  if(errh.ge.1.0d-4) then
   write(*,*)' divergenza ', k,m_per, errh
   stop
  endif



 enddo

 erm = ermed/Dfloat(k_max)
 write(*,*)' Err_med ', N_step, erm


 !-----------------------
 !
 ! Reversibility test
 !
 !----------------------
 !
 !do k=1,k_max
 ! call sym4(xf,yf,tau,pxf,pyf,ptau,-dt)
 ! E = H(xf,yf,tau,pxf,pyf,ptau)
 ! errh = abs(E-E0)
 ! ermed = ermed+errh
 ! t = t-dt
 ! call Fix_cor(xf,yf,pxf,pyf, x,y,vx,vy,px,py, xp,yp,pxp,pyp,t)
 ! ! call Transf(xf,yf,pxf,pyf, x,y,vx,vy,px,py, xp,yp,pxp,pyp,t)
 ! if(mod(k,N_step).eq.0) then
 !  dist0 = sqrt((x-x0)**2 + (vx-vx0)**2)
 !  dist = sqrt((xf-xf0)**2 + (yf-yf0)**2 + (pxf-pxf0)**2 + (pyf-pyf0)**2)
 !  m_per = k/N_step
 !  if(mod(m_per,100).eq.0) write(*,*) m_per
 !  write(7,'(2I8,8g13.5)') m_per,k,t,x,vx,y,vy,dist0,errh
 ! endif
 !enddo
 !
 !write(*,*) ' -- rev test --'
 !write(*,'(5g12.5,4g10.2)') t,x,x0,y,y0,dist0,dist,errh
 !erm = ermed/(2.d0*k_max)
 !write(*,*) ' err medio H ',erm

 end program threebodies_sho

