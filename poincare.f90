
     module data
      REAL*8 :: OME = 1.d0
      REAL*8 :: MU = 0.000954d0 !0.001d0
     end module


      PROGRAM MAIN
      use data
      IMPLICIT NONE
      REAL*8  :: dt , duepi, k_max_r8, dist
      REAL*8  :: ptau , ptau0 , px , px0 , py , py0 , vx0,vy0,t , tau , tau0 , n_step_r8
      REAL*8  :: t_per, x , x0 , y , y0, k_r8 , m_per , r10,r20,c0,s0
      REAL*8  :: xp0, yp0, pxp0, pyp0, J0, xcor0, ycor0, vxcor0, vycor0, Jac0
      REAL*8  ::   xold, yold, vxold, vyold,distcor
      REAL*8  :: xcor1, ycor1, vxcor1, vycor1, r11, r21,taucor
      REAL*8,allocatable::xin(:),vxin(:)
      INTEGER :: k , k_max  , n_step, n_period, i, n, r,s,m

      duepi = 8.d0*DATAN(1.D0)
      t_per = duepi/OME
      n_step = 1000 
      n_period = 1000 


      OPEN (UNIT=15,FILE='PROVA.txt')        
      tau0 = 0.d0
      ptau0 = 1.d0       
      Jac0 = 3.07d0

      c0 = COS(OME*tau0)
      s0 = SIN(OME*tau0)


      xcor0 =  0.55d0
      ycor0 =  0.0d0
      
      r10 = SQRT((xcor0+MU)**2+ycor0**2)
      r20 = SQRT((xcor0-1.d0+MU)**2+ycor0**2) 
 
      vxcor0 = 0.d0 !vxin(m) 
       vycor0 = sqrt(2.d0*(1.d0-MU)/r10 + (2.d0*MU)/r20 + xcor0**2.d0 + ycor0**2.d0 - Jac0 - vxcor0**2.d0)
      

     ! ANTIORARIO DA COROT A FIX
     
     x0 = c0*xcor0 - s0*ycor0
     y0 = s0*xcor0 + c0*ycor0
     px0=   -s0*xcor0 - c0*ycor0 + c0*vxcor0 - vycor0*s0 
     py0 =  c0*xcor0 - s0*ycor0 + s0*vxcor0 + c0*vycor0     
     
     ! re-declaration of variables
     
      x = x0
      y = y0
      tau = tau0
      px = px0
      py = py0
      ptau = ptau0
      n_step_r8 = real(n_step)
      dt = t_per/n_step_r8
       
       
	xold = xcor0
        yold = ycor0
        vxold = vxcor0
        vyold = vycor0    
       
      k_max = n_period*n_step
      i = 1
      k = 0
      t = 0.d0
       
      !!!!!!!!!!!!!!!
      ! forward integration
      !!!!!!!!!!!!!!!
      DO k = 1 , k_max
         k_r8 = real(k)
         CALL SYM4(x,y,tau,px,py,ptau,dt)
         t = k_r8*dt
        
        
        ! TO create the POINCARE SECTION
        
        xcor1 = x*cos(t)+y*sin(t)
        ycor1 = -x*sin(t)+y*cos(t)
        vxcor1 = px*cos(t)+py*sin(t)-x*sin(t)+y*cos(t)
        vycor1 = -px*sin(t)+py*cos(t)-x*cos(t)-y*sin(t)
        
        taucor = ycor1/(yold-ycor1)
        
        
        
        distcor=sqrt(xcor1**2 + ycor1**2 + vxcor1**2 + vycor1**2)
        if(distcor.gt. 8.d0) then
        exit
        end if
        
        
        
        if (ycor1*yold.LT.0.AND.vyold+(vycor1-vyold)*taucor.GT.0) then
        write(15,*) xold+(xcor1-xold)*taucor, vxold+(vxcor1-vxold)*taucor, x,px
        end if
	xold = xcor1
        yold = ycor1
        vxold = vxcor1
        vyold = vycor1	    
      ENDDO
      END PROGRAM MAIN
 
 
 
      SUBROUTINE F(X,Y,Tau,Fx,Fy,Ftau)
      use data
      IMPLICIT NONE
      REAL*8 :: c, Ftau , Fx , Fy , r1 , r2 , s , Tau , X , Y
       
      c = COS(OME*Tau)
      s = SIN(OME*Tau)
 
      r1 = SQRT((X+MU*c)**2+(Y+MU*s)**2) 
      r2 = SQRT((X-(1.d0-MU)*c)**2+(Y-s*(1.d0-MU))**2)  
 
      Fx = -((1.d0-MU)*(X+MU*c))/r1**3 - (MU*(X-c*(1.d0-MU)))/r2**3
      Fy = -((1.d0-MU)*(Y+MU*s))/r1**3 - (MU*(Y-s*(1.d0-MU)))/r2**3
 
      Ftau =-((1.d0-MU)*OME*MU*(-s*(X+MU*c)+c*(Y+MU*s)))/r1**3
      
      Ftau = Ftau - (MU*(1.d0-MU)*OME*(s*(X-(1.d0-MU)*c)-c*(Y-(1.d0-MU)*s)))/r2**3
 
      END
  
      SUBROUTINE SYM2(X,Y,Tau,Px,Py,Ptau,Dt)
      IMPLICIT NONE
      REAL*8 :: Dt , ftau , ftaunew , fx , fxnew , fy , fynew , Ptau
      REAL*8 :: ptaunew , Px , pxnew , Py , pynew , Tau , taunew , X
      REAL*8 :: xnew , Y , ynew
 
      CALL F(X,Y,Tau,fx,fy,ftau)
 
 
      xnew = X + Px*Dt + 0.5d0*fx*(Dt**2)
      ynew = Y + Py*Dt + 0.5d0*fy*(Dt**2) 
      taunew = Tau + Dt
 
      CALL F(xnew,ynew,taunew,fxnew,fynew,ftaunew)
      pxnew = Px + 0.5d0*Dt*(fx+fxnew)
      pynew = Py + 0.5d0*Dt*(fy+fynew)
      ptaunew = Ptau + 0.5d0*Dt*(ftau+ftaunew)
 
      X = xnew
      Y = ynew
      Tau = taunew
      Px = pxnew
      Py = pynew
      Ptau = ptaunew
 
      END

      SUBROUTINE SYM4(X,Y,Tau,Px,Py,Ptau,Dt)
      IMPLICIT NONE
      REAL*8 :: alpha , beta , Dt , dt1 , dt2 , Ptau , Px , Py , sq2
      REAL*8 :: Tau , X , Y
      sq2 = 2.d0**(1.D0/3.D0)
      alpha = 1.D0/(2.d0-sq2)
      beta = sq2/(2.d0-sq2)
      dt1 = Dt*alpha
      dt2 = -Dt*beta
      CALL SYM2(X,Y,Tau,Px,Py,Ptau,dt1)
      CALL SYM2(X,Y,Tau,Px,Py,Ptau,dt2)
      CALL SYM2(X,Y,Tau,Px,Py,Ptau,dt1)
      END
