 !#define USA_OPENMP

 module parameters
 implicit none

 integer, parameter ::  N_step = 1000 ! numero passi per periodo
 integer, parameter ::  N_periods = 10 ! numero periodi
 integer, parameter ::  N_points = 1000 ! numero punti
 integer, parameter ::  N_dim = 6 ! numero dimensioni

 logical, parameter :: chaos = .false.  !  x0 = 0.55 stabile, x0 = 0.56 caotica

 real (kind=8), parameter :: ome = 1.d0
 real (kind=8), parameter :: mu = 9.54d-4 ! mu = m_G/m_S
 real (kind=8), parameter :: eps = 1.d0 ! eps ampiezza rumore
 !real (kind=8), parameter :: seme = 7255347848.d0 ! Parametri Monte-carlo
 real (kind=8), parameter :: stocas_fisso=0.031364026959689346589003385d0

 end module parameters







 module constants
 implicit none
 real (kind=8), parameter :: pi = 3.141592653589793
 end module constants






 module utilities_threebodies
 use constants
 implicit none
 contains




 pure real*8 function stocas()
 use parameters
 !real (kind=8), intent(in) :: IX
 !real (kind=8), parameter :: a = 16807.d0
 !real (kind=8), parameter :: b = 2147483647.d0
 !real (kind=8) :: p

 !p = mod(a*IX,b)
 !stocas=p/b

 !p = mod(a*seme,b)
 !stocas=p/b

 !stocas=67353735.d0/2147483647.d0

 stocas=stocas_fisso

 return
 end function stocas






 real (kind=8) function f_cor(x)
 Implicit none
 real (kind=8), intent(in) :: x

 f_cor=cos(2.d0*pi*x)
 return
 end function f_cor






 subroutine user_map (X,N_s,noise,ome,eps)
 Implicit none
 integer, parameter ::  N_dim=6 ! numero dimensioni
 real (kind=8), dimension(1:N_dim), intent(inout) :: X
 integer, intent(in) :: N_s
 real (kind=8), intent(in) :: noise,ome,eps
 integer :: k
 real (kind=8) :: xf,pxf,yf,pyf,tau,ptau
 real (kind=8) :: T_per,A_step,dt

 xf = X(1)
 pxf = X(2)
 yf = X(3)
 pyf = X(4)
 tau = X(5)
 ptau = X(6)

 T_per = 2.d0*pi/ome

 A_step = dble(N_s)
 dt = T_per/A_step

 do k=1,N_s
  call sym4(xf,yf,tau,pxf,pyf,ptau,dt,ome,eps)
 enddo

 X(1) = xf + eps*noise
 X(2) = pxf
 X(3) = yf
 X(4) = pyf
 X(5) = tau
 X(6) = ptau

 end subroutine user_map






 subroutine Cor_fix(x,y,vx,vy,xf,yf,pxf,pyf,t,ome)
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

 end subroutine Cor_fix






 subroutine Fix_cor(xf,yf,pxf,pyf,x,y,vx,vy,t,ome)
 Implicit none
 real (kind=8) :: c, s
 real (kind=8), intent(in) :: xf,yf,pxf,pyf,t,ome
 real (kind=8), intent(out) :: x,y,vx,vy

 c = COS(ome*t)
 s = SIN(ome*t)
 x = c*xf + s*yf
 y = -s*xf + c*yf
 vx = c*pxf + s*pyf - s*xf + c*yf
 vy = -s*pxf + c*pyf - c*xf - s*yf

 end subroutine Fix_cor






 subroutine sym2(x,y,tau,px,py,ptau,dt,ome,eps)
 Implicit none
 real (kind=8), intent(inout) :: x,y,tau,px,py,ptau
 real (kind=8), intent(in) :: dt,ome,eps
 real (kind=8) :: xnew,ynew,taunew,pxnew,pynew,ptaunew,fx,fy,ftau,fxnew,fynew,ftaunew

 call f(x,y,tau,fx,fy,ftau,ome,eps)
 xnew = x + px*dt + fx*dt**2/2.d0
 ynew = y + py*dt + fy*dt**2/2.d0
 taunew = tau + dt

 call f(xnew,ynew,taunew,fxnew,fynew,ftaunew,ome,eps)
 pxnew = px + dt*(fx+fxnew)/2.d0
 pynew = py + dt*(fy+fynew)/2.d0
 ptaunew = ptau + dt*(ftau+ftaunew)/2.d0

 x=xnew
 y=ynew
 tau=taunew
 px=pxnew
 py=pynew
 ptau=ptaunew

 contains
 subroutine f(x,y,tau,fx,fy,ftau,ome,eps)
 real (kind=8) :: c, s
 real (kind=8) :: r1,r2
 real (kind=8), intent(in) :: x,y,tau,ome,eps
 real (kind=8), intent(out) :: fx,fy,ftau

 c = cos(ome*tau)
 s = sin(ome*tau)

 r1 = sqrt((x + eps*c)**2 + (y + eps*s)**2)
 r2 = sqrt((x - (1.d0-eps)*c)**2 + (y - s*(1.d0 - eps))**2)

 fx = -(1.d0 - eps) * (x+eps*c) / r1**3 - eps * (x - c*(1.d0-eps)) / r2**3
 fy = -(1.d0 - eps) * (y+eps*s) / r1**3 - eps * (y - s*(1.d0-eps)) / r2**3
 ftau = -(1.d0 - eps) * ome * eps * (-s * (x + eps*c) + c * (y + eps*s)) / r1**3
 ftau = ftau - eps * (1.d0 - eps) * ome * (s * (x - (1.d0-eps)*c) - c * (y - (1.d0 - eps)*s)) / r2**3
 end subroutine f
 end subroutine sym2





 subroutine sym4(x,y,tau,px,py,ptau,dt,ome,eps)
 Implicit none
 real (kind=8), intent(inout) :: x,y,tau,px,py,ptau
 real (kind=8), intent(in) :: dt,ome,eps
 real (kind=8) :: sq2, alpha, beta, dt1, dt2
 sq2 = 2.d0**(1.d0/3.d0)
 alpha = 1.d0/(2.d0-sq2)
 beta = sq2/(2.d0-sq2)
 dt1 = dt*alpha
 dt2 = -dt*beta
 call sym2(x,y,tau,px,py,ptau,dt1,ome,eps)
 call sym2(x,y,tau,px,py,ptau,dt2,ome,eps)
 call sym2(x,y,tau,px,py,ptau,dt1,ome,eps)
 end subroutine sym4





 end module utilities_threebodies







 program threebodies
#if defined(USA_OPENMP)
 use omp_lib
#endif
 use utilities_threebodies
 use constants
 use parameters
 implicit none


 integer :: AllocateStatus, DeAllocateStatus
 real (kind=8), dimension(:, :, :), allocatable :: x, x_eps
 real (kind=8), dimension(:), allocatable :: x_3corpi
 real (kind=8), dimension(:, :), allocatable :: xx_no_noise, xx_noise, xx_i

 real (kind=8), dimension(:), allocatable :: Cor_eps
 real (kind=8), dimension(:), allocatable :: Fid_eps
 real (kind=8), dimension(:), allocatable :: Av_f, Av_mu_f
 real (kind=8), dimension(:), allocatable :: Av_mu_eps_f

 real (kind=8) :: T_per
 real (kind=8) :: Jac
 real (kind=8) :: xc, yc, x_Sole, y_Sole, x_Giove, y_Giove
 real (kind=8) :: xcor0, ycor0, vxcor0, vycor0, t
 real (kind=8) :: tau0, ptau0
 real (kind=8) :: xf0, yf0, pxf0, pyf0
 real (kind=8) :: xf, yf, pxf, pyf
 real (kind=8) :: xx,yy,vvx,vvy
 real (kind=8) :: stocastic_noise

 Integer :: n, np, j

#if defined(USA_OPENMP)
 Integer :: rank, nthreads
#endif 




 allocate (x(N_dim,N_periods+1,N_points), STAT = AllocateStatus)
 allocate (x_eps(N_dim,N_periods+1,N_points), STAT = AllocateStatus)
 allocate (x_3corpi(N_dim), STAT = AllocateStatus)
 allocate (xx_no_noise(N_dim,N_points), STAT = AllocateStatus)
 allocate (xx_noise(N_dim,N_points), STAT = AllocateStatus)
 allocate (xx_i(N_dim,N_points), STAT = AllocateStatus)
 allocate (Cor_eps(N_periods+1), STAT = AllocateStatus)
 allocate (Fid_eps(N_periods+1), STAT = AllocateStatus)
 allocate (Av_f(N_periods+1), STAT = AllocateStatus)
 allocate (Av_mu_f(N_periods+1), STAT = AllocateStatus)
 allocate (Av_mu_eps_f(N_periods+1), STAT = AllocateStatus)

 IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

 !Inizializzazione a zero di tutti gli array
 x = 0.d0
 x_eps = 0.d0
 x_3corpi = 0.d0
 xx_no_noise = 0.d0
 xx_noise = 0.d0
 xx_i = 0.d0

 !Inizializzazione a zero di tutte le variabili
 Cor_eps = 0.d0
 Fid_eps = 0.d0
 Av_f = 0.d0
 Av_mu_f = 0.d0
 Av_mu_eps_f = 0.d0
 Jac = 0.d0
 T_per = 0.d0
 xcor0 = 0.d0
 ycor0 = 0.d0
 vxcor0 = 0.d0
 vycor0 = 0.d0
 t = 0.d0
 tau0 = 0.d0
 ptau0 = 0.d0
 xf0 = 0.d0
 yf0 = 0.d0
 pxf0 = 0.d0
 pyf0 = 0.d0
 xf = 0.d0
 yf = 0.d0
 pxf = 0.d0
 pyf = 0.d0
 xx = 0.d0
 yy = 0.d0
 vvx = 0.d0
 vvy = 0.d0
 stocastic_noise = 0.d0

 !Inizializzazione a zero di tutti i contatori dei cicli loop
 n = 0
 np = 0
 j = 0



 open(4,file='Correlation.dat')
 open(9,file='Initial_conditions.dat')
 open(19,file='Evoluted.dat')
 open(20,file='Evoluted_with_noise.dat')
 open(50,file='Initial_data.dat')


 T_per = (2.d0 * pi) / ome
 Jac = 3.07d0

 !-----------------------------------------------
 ! Condizioni iniziali sistema corotante
 !-----------------------------------------------
 if(chaos .eqv. .true.) then
  xcor0 = 0.56d0
 else
  xcor0 = 0.55d0
 endif
 ycor0 = 0.d0
 vxcor0 = 0.d0
 t = 0.d0

 !------------------------------------------
 ! Condizioni iniziali sistema fisso
 !------------------------------------------
 tau0=0.d0
 ptau0=1.d0


 !------------ fine parametri, inizio verifiche e algoritmo
 vycor0 = xcor0**2 + ycor0**2 + 2.d0*(1.d0-mu)/abs(xcor0+mu) + 2.d0*mu/abs(xcor0-1.d0+mu) - Jac - vxcor0**2
 vycor0 = sqrt(vycor0)

 write(*,*) 'vy0 mu :', vycor0, mu
 write(*,*) " x0   ,  vycor0  ,  JAC "
 write(*,*)  xcor0 ,  vycor0  ,  Jac

 call Cor_fix(xcor0,ycor0,vxcor0,vycor0,xf0,yf0,pxf0,pyf0,t,ome)

 write(*,*) "  x0  ,  y0   ,  vx0   ,  vy0 "
 write(*,*)  xcor0 , ycor0 , vxcor0 , vycor0
 write(*,*) " tau0 , ptau0 "
 write(*,*)   tau0 , ptau0
 write(*,*) " xf0  , yf0   ,  pxf0  ,  pyf0 "
 write(*,*)   xf0  , yf0   ,  pxf0  ,  pyf0
 write(*,*) "  t   , ome "
 write(*,*)    t   , ome

 stocastic_noise = 1.d0 - 2.d0 * stocas()
 write(*,*) " stocastic noise "
 write(*,*)   stocastic_noise


 Cor_eps = 0.d0
 Fid_eps = 0.d0
 Av_f = 0.d0
 Av_mu_f = 0.d0
 Av_mu_eps_f = 0.d0


 x(1,1,1) = xf0
 x(2,1,1) = pxf0
 x(3,1,1) = yf0
 x(4,1,1) = pyf0
 x(5,1,1) = tau0
 x(6,1,1) = ptau0

 x_eps(1,1,1) = xf0
 x_eps(2,1,1) = pxf0
 x_eps(3,1,1) = yf0
 x_eps(4,1,1) = pyf0
 x_eps(5,1,1) = tau0
 x_eps(6,1,1) = ptau0



 x_3corpi(1) = xf0
 x_3corpi(2) = pxf0
 x_3corpi(3) = yf0
 x_3corpi(4) = pyf0
 x_3corpi(5) = tau0
 x_3corpi(6) = ptau0

 t = 0.d0



 do j=2,N_points ! loop speciale su tutti i punti allo step n = 0

  call user_map(X_3corpi,N_step,0.d0,ome,eps)

  x(1,1,j) = X_3corpi(1)
  x(2,1,j) = X_3corpi(2)
  x(3,1,j) = X_3corpi(3)
  x(4,1,j) = X_3corpi(4)
  x(5,1,j) = X_3corpi(5)
  x(6,1,j) = X_3corpi(6)

  x_eps(1,1,j) = X_3corpi(1)
  x_eps(2,1,j) = X_3corpi(2)
  x_eps(3,1,j) = X_3corpi(3)
  x_eps(4,1,j) = X_3corpi(4)
  x_eps(5,1,j) = X_3corpi(5)
  x_eps(6,1,j) = X_3corpi(6)

 enddo



 do j=1,N_points

  xf = x(1,1,j)
  pxf = x(2,1,j)
  yf = x(3,1,j)
  pyf = x(4,1,j)
  t = j * T_per

  call Fix_cor(xf,yf,pxf,pyf,xx,yy,vvx,vvy,t,ome)

  xx_i(1,j) = xx
  xx_i(2,j) = vvx
  xx_i(3,j) = yy
  xx_i(4,j) = vvy

 enddo


 do j=1,N_points
  write(50,'(I5,10g13.5)') j, x(1,1,j), x(3,1,j), x(2,1,j), x(4,1,j), x(5,1,j)/(2.d0*pi), x(6,1,j)
  write(9,'(I5,10g13.5)') j, xx_i(1,j), xx_i(3,j), xx_i(2,j), xx_i(4,j)
 enddo




#if defined (USA_OPENMP)
 !$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED(x,x_eps,xx_noise,xx_no_noise)
 rank = omp_get_thread_num()
 nthreads = omp_get_num_threads()
 if (rank==0) then
  write(*,*) "Running on ", nthreads, " threads"
 endif
 !$OMP DO
#endif
 do j=1,N_points ! Loop sui punti, secondario

  do n=1,N_periods ! Iterazione nel tempo: loop principale (N_periods e' il massimo num di iterazioni temporali, n=0 inizializzato prima)

   X_3corpi(1) = x(1,n,j)
   X_3corpi(2) = x(2,n,j)
   X_3corpi(3) = x(3,n,j)
   X_3corpi(4) = x(4,n,j)
   X_3corpi(5) = x(5,n,j)
   X_3corpi(6) = x(6,n,j)

   call user_map(X_3corpi,N_step,0.d0,ome,eps)

   x(1,n+1,j) = X_3corpi(1)
   x(2,n+1,j) = X_3corpi(2)
   x(3,n+1,j) = X_3corpi(3)
   x(4,n+1,j) = X_3corpi(4)
   x(5,n+1,j) = X_3corpi(5)
   x(6,n+1,j) = X_3corpi(6)

   X_3corpi(1) = x_eps(1,n,j)
   X_3corpi(2) = x_eps(2,n,j)
   X_3corpi(3) = x_eps(3,n,j)
   X_3corpi(4) = x_eps(4,n,j)
   X_3corpi(5) = x_eps(5,n,j)
   X_3corpi(6) = x_eps(6,n,j)


   call user_map(X_3corpi,N_step,stocastic_noise,ome,eps)

   x_eps(1,n+1,j) = X_3corpi(1)
   x_eps(2,n+1,j) = X_3corpi(2)
   x_eps(3,n+1,j) = X_3corpi(3)
   x_eps(4,n+1,j) = X_3corpi(4)
   x_eps(5,n+1,j) = X_3corpi(5)
   x_eps(6,n+1,j) = X_3corpi(6)

  enddo
 enddo
#if defined (USA_OPENMP)
 !$OMP ENDDO
 !$OMP END PARALLEL
#endif





 do j=1,N_points ! Loop sui punti, secondario

  xf  = x(1,N_periods+1,j)
  pxf = x(2,N_periods+1,j)
  yf  = x(3,N_periods+1,j)
  pyf = x(4,N_periods+1,j)
  t   = (N_periods)*T_per

  call Fix_cor(xf,yf,pxf,pyf,xx,yy,vvx,vvy,t,ome)

  xx_no_noise(1,j) = xx
  xx_no_noise(2,j) = vvx
  xx_no_noise(3,j) = yy
  xx_no_noise(4,j) = vvy

  xf  = x_eps(1,N_periods+1,j)
  pxf = x_eps(2,N_periods+1,j)
  yf  = x_eps(3,N_periods+1,j)
  pyf = x_eps(4,N_periods+1,j)
  t   = (N_periods)*T_per

  call Fix_cor(xf,yf,pxf,pyf,xx,yy,vvx,vvy,t,ome)

  xx_noise(1,j) = xx
  xx_noise(2,j) = vvx
  xx_noise(3,j) = yy
  xx_noise(4,j) = vvy


  write(19,'(I5,10g13.5 )') j, xx_no_noise(1,j), xx_no_noise(3,j), xx_no_noise(2,j), xx_no_noise(4,j)
  write(20,'(I5,10g13.5 )') j, xx_noise(1,j), xx_noise(3,j), xx_noise(2,j), xx_noise(4,j)
 enddo




 do n=1,N_periods+1 ! Iterazione nel tempo: loop principale (N_periods e' il massimo num di iterazioni temporali)

  do j=1,N_points
   Cor_eps(n) = Cor_eps(n) + f_cor(x_eps(1,n,j)) * f_cor(x_eps(1,1,j))  ! attenzione: f_cor e' di x_eps(1,1,j) e non x_eps(1,n,j)!
   Fid_eps(n) = Fid_eps(n) + f_cor(x_eps(1,n,j)) * f_cor(x(1,n,j))

   Av_f(n) = Av_f(n) + f_cor(x(1,1,j))  ! attenzione: f_cor e' di x(1,1,j) e non x(1,n,j)!
   Av_mu_f(n) = Av_mu_f(n) + f_cor(x(1,n,j))
   Av_mu_eps_f(n) = Av_mu_eps_f(n) + f_cor(x_eps(1,n,j))
  enddo

  Cor_eps(n) = Cor_eps(n) / N_points
  Fid_eps(n) = Fid_eps(n) / N_points

  Av_f(n) = Av_f(n) / N_points
  Av_mu_f(n) = Av_mu_f(n) / N_points
  Av_mu_eps_f(n) = Av_mu_eps_f(n) / N_points

  Cor_eps(n) = Cor_eps(n) - Av_f(n) * Av_mu_eps_f(n)
  Fid_eps(n) = Fid_eps(n) - Av_mu_f(n) * Av_mu_eps_f(n)

 enddo   ! chiude ciclo n=1,N_periods



 do n=1,N_periods+1
  write(4,'(I8,2g14.4)') n-1, Cor_eps(n), Fid_eps(n)
 enddo

 deallocate (x, STAT = DeAllocateStatus)
 deallocate (x_eps, STAT = DeAllocateStatus)
 deallocate (x_3corpi, STAT = DeAllocateStatus)
 deallocate (xx_no_noise, STAT = DeAllocateStatus)
 deallocate (xx_noise, STAT = DeAllocateStatus)
 deallocate (xx_i, STAT = DeAllocateStatus)
 deallocate (Cor_eps, STAT = DeAllocateStatus)
 deallocate (Fid_eps, STAT = DeAllocateStatus)
 deallocate (Av_f, STAT = DeAllocateStatus)
 deallocate (Av_mu_f, STAT = DeAllocateStatus)
 deallocate (Av_mu_eps_f, STAT = DeAllocateStatus)

 IF (DeAllocateStatus /= 0) STOP "*** Problem when deallocating memory ***"



 end program threebodies


