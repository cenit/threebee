      Implicit  real*8(A-H,O-Z)
      real*8 seme, stocas
      Dimension x_0(1000000), y_0(1000000),  tau_0(1000000)
      Dimension   x(1000000),    y(1000000),   tau(1000000)
      Dimension x_eps(1000000),y_eps(1000000),tau_eps(1000000)

      Dimension px_0(1000000),  py_0(1000000),  ptau_0(1000000)
      Dimension   px(1000000),    py(1000000),    ptau(1000000)
      Dimension px_eps(1000000),py_eps(1000000)
     *     ,ptau_eps(1000000)

c-------------
c     coordinate e momenti nel sistema fisso
c-------------

      Dimension  x_v(0:1000),y_v(0:1000)
      Dimension  Cor(0:1000), Cor_eps(0:1000)
      Dimension  Fid(0:1000), Fid_eps(0:1000)
      Dimension  Av_f(0:1000), Av_mu_f(0:1000)
      Dimension  Av_mu_eps_f(0:1000)
      Dimension  X_3_corpi(10), a_3_corpi(10)
      Real*8 mu , eps , Jac
      Integer user_map

      common ome,mu
c---------------
c     Cor  correlazio
c     Fid  fidelity
c     Av_f           media di f con Lebesgue 
c     Av_mu_f        media di f con misura mu
c     Av_mu_eps_f    media di f son misura stazionaria (rumore) 
c---------------
c     
c     Standard Map 
c     
c     --------------

      open(4,file='Corr.dat')

c---- Parametri

      duepi=8*datan(1.d0)

      ome=1
      T_per=duepi/ome

      mu=0.000954D0             ! mu  = m_G/m_S
      eps= 1.E0                 ! eps ampiezza rumore


      N_step=1000               ! numero  passi per periodo
      N_periods=10              ! numero periodi
c---- 
      N_max=N_periods


      xc=0.5d0-eps              !  punto equilibio
      yc=sqrt(3.d0)/2.d0        !  di Lagrange satellite

      x_Sole=-eps
      y_Sole=0
      x_Giove=1-eps
      y_Giove=0


c-----------------------------------------
c     Costante Jacobi Jac = - 2 H  
c     H_c   hamiltoniana . sit. corot.
c     Jac valore scelto costante Jacobi
c-------------------------------------------

      H_c= -1 -  0.5d0* (xc**2+yc**2)
      Jac_c=-2*H_c
      Jac=3.07D0
c-----------------------------------------------
c     Condizioni iniziali sistema corotante 
c-----------------------------------------------
      x0=0.55D0                 !  0.56  caotico    x=0.55 regolare
      y0=0
      vx0=0
c-----
c     calcolo di vy0
c-----

      vvy0= x0**2+2*(1-mu)/abs(x0+mu)+2*mu/abs(x0-1+mu)-Jac

      vvy0= vvy0-vx0**2
      if(vvy0.lt.0)  then
         write(*,*) '**** vy0**2 < 0: cambiare cond. iniz. ',  vvy0
         stop
      endif
      vy0=sqrt (vvy0)
      write(*,*)'*** vy0   mu ***',  vy0,mu,x0, vvy0, Jac
      t=0
c------------------------------------------
c     Condizioni iniziali sistema fisso
c------------------------------------------
      call Cor_fix(x0,y0,vx0,vy0,xf0,yf0,pxf0,pyf0,t)

      tau0=0.d0
      ptau0=1.d0

c----------
c     Parametri Monte-carlo
c----------
      N_max=N_periods           !   iterazioni


      M_max=1E3                 !   numero totale dei punti


      A_max=M_max
      seme=7255347848.d0



      a_3_corpi(1)= N_step
      a_3_corpi(2)= N_periods

      zero=0
c---------------------------
c     Punti iniziali
c     
c-----------------------------
c     
c     Condizioni iniziali
c-----------------------------
      n=0

      do n=0,N_max              !  -------------------------------------->
c-----
c     loop principale iterazioni mappa
         write(*,*) ' *****  n =  *********', n
c-----



         if (n.eq.0) then       ! ------->   Caso  n=0


            x_0(1)    =xf0
            y_0(1)    =yf0
            tau_0(1)  =tau0
            px_0(1)   =pxf0
            py_0(1)   =pyf0
            ptau_0(1) =ptau0

            x(1)   = x_0(1)
            y(1)   = y_0(1)
            tau(1) = tau_0(1)
            px(1)  = px_0(1)
            py(1)  = py_0(1)
            ptau(1)= ptau_0(1)

            x_eps(1)   = x_0(1)
            y_eps(1)   = y_0(1)
            tau_eps(1) = tau_0(1)
            px_eps(1)  = px_0(1)
            py_eps(1)  = py_0(1)
            ptau_eps(1)= ptau_0(1)




            X_3_corpi(1)= x_0(1)
            X_3_corpi(2)= px_0(1)
            X_3_corpi(3)= y_0(1)
            X_3_corpi(4)= py_0(1)
            X_3_corpi(5)= tau_0(1)
            X_3_corpi(6)= ptau_0(1)

            J=0
            t=0

            do j=1,M_max        !  --->   Inizia iterazione sui punti

               ii=user_map (X_3_corpi,a_3_corpi,zero,j+1)

               x_0(j+1)   =X_3_corpi(1)
               px_0(j+1)  =X_3_corpi(2)
               y_0(j+1)   =X_3_corpi(3)
               py_0(j+1)  =X_3_corpi(4)
               tau_0(j+1) =X_3_corpi(5)
               ptau_0(j+1)=X_3_corpi(6)

               x(j+1)   =x_0(j+1)
               px(j+1)  =px_0(j+1)
               y(j+1)   =y_0(j+1)
               py(j+1)  =py_0(j+1)
               tau(j+1) =tau_0(j+1)
               ptau(j+1)=ptau_0(j+1)


               x_eps(j+1)   =x_0(j+1)
               px_eps(j+1)  =px_0(j+1)
               y_eps(j+1)   =y_0(j+1)
               py_eps(j+1)  =py_0(j+1)
               tau_eps(j+1) =tau_0(j+1)
               ptau_eps(j+1)=ptau_0(j+1)

               xf=x_0(j)
               yf=y_0(j)
               pxf=px_0(j)
               pyf=py_0(j)

               t=j*T_per
               call Fix_cor(xf,yf,pxf,pyf, xx,yy,vvx,vvy,t)
               a1=x_0(j)
               a2=y_0(j)
               a3=px_0(j)
               a4=py_0(j)
               a5=tau_0(j)/duepi
               a6=ptau_0(j)
               write(50,'(I3,10g13.5 )') j,a1,a2,a3,a4,a5,a6
               write(9,'(I3,10g13.5 )') j,  xx,yy,vvx,vvy

            enddo               !   <---

         else                   ! -------- Caso n>0

c     write(2,*) ' *****  n =  *********', n

            do j=1,M_max        ! --->
c---------
c     mappa deterministica 
c---------
               X_3_corpi(1)= x_0(j)
               X_3_corpi(2)= px_0(j)
               X_3_corpi(3)= y_0(j)
               X_3_corpi(4)= py_0(j)
               X_3_corpi(5)= tau_0(j)
               X_3_corpi(6)= ptau_0(j)

               t=0
               do np=1,n
                  ii=user_map (X_3_corpi,a_3_corpi,zero,np)
               enddo

               x(j)    =  X_3_corpi(1)
               px(j)   =  X_3_corpi(2)
               y(j)    =  X_3_corpi(3)
               py(j)   =  X_3_corpi(4)
               tau(j)  =  X_3_corpi(5)
               ptau(j) =  X_3_corpi(6)
c---------
c     mappa stocastica
c--------

               X_3_corpi(1)= x_0(j)
               X_3_corpi(2)= px_0(j)
               X_3_corpi(3)= y_0(j)
               X_3_corpi(4)= py_0(j)
               X_3_corpi(5)= tau_0(j)
               X_3_corpi(6)= ptau_0(j)


               do  np=1,n
                  xi=1-2*stocas(seme)
                  ii=user_map (X_3_corpi,a_3_corpi,xi,np)
               enddo

               x_eps(j)   = X_3_corpi(1)
               px_eps(j  )= X_3_corpi(2)
               y_eps(j)   = X_3_corpi(3)
               py_eps(j)  = X_3_corpi(4)
               tau_eps(j) = X_3_corpi(5)
               ptau_eps(j)= X_3_corpi(6)



               if(n.eq.N_max.and.j.le.2000)  then
c     write(10+n,'(10g13.5)')   x(j),y(j), px(j), py(j),   
c     *      x_eps(j),y_eps(j), px_eps(j), py_eps(j) 
                  xf=x(j)
                  yf=y(j)
                  pxf=px(j)
                  pyf=py(j)
                  t=n*T_per
                  call Fix_cor(xf,yf,pxf,pyf, xx,yy,vvx,vvy,t)
                  write(19,'(I3,10g13.5 )') j,  xx,yy,vvx,vvy
                  xf=x_eps(j)
                  yf=y_eps(j)
                  pxf=px_eps(j)
                  pyf=py_eps(j)
                  t=n*T_per
                  call Fix_cor(xf,yf,pxf,pyf, xx,yy,vvx,vvy,t)
                  write(20,'(I3,10g13.5 )') j,  xx,yy,vvx,vvy
               endif


            enddo               ! <---

         endif                  !  <------------






c------------------------------
c     
c     Calcolo correlazione
c     
c-----------------------------


         Cor_eps(n)=0
         Fid_eps(n)=0

         Av_f(n)=0
         Av_mu_f(n)=0
         Av_mu_eps_f(n)=0

         do j=0,M_max
c     write(3,'(I5,8g12.4)') J, x_0(j), y_0(j),
c     *                 x(j),y(j), x_eps(j), y_eps(j)

            Cor_eps(n)= Cor_eps(n)+ f_cor(x_eps(j),y_eps(j) )
     *           * f_cor(x_0(j),y_0(j) )
            Fid_eps(n)= Fid_eps(n)+ f_cor(x_eps(j),y_eps(j) )
     *           * f_cor(x(j),y(j))
c-----
            Av_f(n)=Av_f(n)+ f_cor(x_0(j),y_0(j) )
            Av_mu_f(n)=Av_mu_f(n)+ f_cor(x(j), x(j))
            Av_mu_eps_f(n)=Av_mu_eps_f(n)+ f_cor(x_eps(j),y_eps(j) )
c     write(3,'(I5,4g13.5)') j, Cor_eps(n), Fid_eps(n)
         enddo
         Cor_eps(n)= Cor_eps(n)/A_max
         Fid_eps(n)= Fid_eps(n)/A_max

         Av_f(n)= Av_f(n)/A_max
         Av_mu_f(n)= Av_mu_f(n)/A_max
         Av_mu_eps_f(n)= Av_mu_eps_f(n)/A_max

c     write(*,'(5g11.3)')  Av_f(n),Av_mu_f(n), Av_mu_eps_f(n)

         Cor_eps(n)= Cor_eps(n)  - Av_f(n)*Av_mu_eps_f(n)
         Fid_eps(n)= Fid_eps(n)  - Av_mu_f(n)*Av_mu_eps_f(n)
c     write(3,'(2I5,4g13.5)') n, j, Cor_eps(n), Fid_eps(n)


      enddo                     ! <---------------------------------------



      do n=0,N_max
         write(4,'(I8,2g14.4)')  n, Cor_eps(n), Fid_eps(n)
      enddo


      end


      real*8 function stocas(IX)
      REAL*8 IX
      IX=DMOD(16807.D0*IX,2147483647.D0)
      stocas=IX/2147483647.D0
      RETURN
      END

      real*8 function fmod(x)
      REAL*8  X
      if(x.ge.0) then
         ix=x
         x=x-ix
      else
         ix=x
         x=x-ix
         x=x+1
      endif
      fmod=x
      return
      end

      real*8 function f_cor(x,y)
      Implicit  real*8(A-H,O-Z)
      duepi=8*datan(1.d0)
      f_cor=cos(duepi*x)
      return
      end

c======================================================




      Integer  Function user_map (X, a,xi, n)
*     
*     lowercase name is MANDATORY
*     
      implicit real*8 (a-h, o-z)
      real*8 X(10),a(10),epsilon,gamma, mu

      common ome,mu
c=======================================================
C     
C     Calcolo standard map e cat map in forma
C     
C     utile per GIOTTO.
C     
C     Calcolo mappa inversa per n<0
C     
c=======================================================
*     epsilon=a(1)
*     
*     ich=a(2)
*     
*     ich>0  standard map     ich<=0  cat map
*     
*     Scegliere la finestra   0 <x<5   -.5 <y<.5
*     
*     Punto elletico nell'origine. Nota segno - in eps sin(x)
*     
c-------------------------------------------------------

      N_step= a(1)
      N_periods =a(2)



*--------------------------------------------------
*     this is only to avoid compiler's complaints:
*--------------------------------------------------
      user_map = 1

      xf=  X(1)
      pxf= X(2)
      yf=  X(3)
      pyf= X(4)
      tau= X(5)
      ptau=X(6)
c     write(7,*) '************* n ************', n 
      if(n.eq.2) then
c     write(9,'(I6,9g13.5)')n-1,xf,pxf,yf,pyf,tau/duepi,ptau,
c     *         t/duepi 
      endif

      duepi=8*datan(1.d0)
      ome=1
      T_per=duepi/ome

      t=0
      A_step=Dfloat(N_step)
      dt=T_per/A_step
      k_max=N_step
c     write(7,*) ,k_max, N_step
      do k=1,k_max
         call sym4(xf,yf,tau,pxf,pyf,ptau,dt)
c     if(mod(k,100).eq.0) write (7,'(2I5,8g12.4)') 
c     *     n-1,k, xf,pxf,yf,pyf,tau/duepi,ptau        
         t=t+dt
      enddo

      X(1)=  xf +eps*xi
      X(2)=  pxf
      X(3)=  yf
      X(4)=  pyf
      X(5)=  tau
      X(6)=  ptau

c     write(9,'(I6,9g13.5)')n-1,xf,pxf,yf,pyf,tau/duepi,ptau,
c     *         t/duepi 

      return
      end


      Subroutine Cor_fix(x,y,vx,vy,xf,yf,pxf,pyf,t)
      Implicit real*8 (A-H,O-Z)
      common ome,eps

      c = COS(ome*t)
      s = SIN(ome*t)

      xf =  c*x - s*y
      yf =  s*x + c*y

      pxf =  c*vx - s*vy - s*x - c*y
      pyf =  s*vx + c*vy + c*x - s*y


      return
      end

      Subroutine Fix_cor(xf,yf,pxf,pyf,x,y,vx,vy,t)
      Implicit real*8 (A-H,O-Z)
      common ome,eps

      c = COS(ome*t)
      s = SIN(ome*t)

      x =  c*xf + s*yf
      y = -s*xf + c*yf

      vx =  c*pxf + s*pyf - s*xf + c*yf
      vy = -s*pxf + c*pyf - c*xf - s*yf

      px=vx-y
      py=vy+x

      xc=0.5d0-eps
      yc=sqrt(3.d0)/2.d0

      pxc=-yc
      pyc= xc

      xp=x-xc
      yp=y-yc

      pxp=px-pxc
      pyp=py-pyc

      return
      end


      subroutine  f(x,y,tau,fx,fy,ftau)
      Implicit real*8 (A-H,O-Z)
      common ome,eps

      c=cos(ome*tau)
      s=sin(ome*tau)

      r1= sqrt((x+eps*c)**2+(y+eps*s)**2)
      r2= sqrt((x-(1-eps)*c)**2+(y-s*(1-eps))**2)


      fx=-(1-eps)*(x+eps*c)/r1**3-eps*(x-c*(1-eps))/r2**3
      fy=-(1-eps)*(y+eps*s)/r1**3-eps*(y-s*(1-eps))/r2**3

      ftau= -(1-eps)*ome*eps*(-s*(x+eps*c)+c*(y+eps*s))/r1**3
      ftaut=eps*(1-eps)*ome*(s*(x-(1-eps)*c)-c*(y-(1-eps)*s))
      ftau=ftau-ftaut/r2**3



      return
      end



      subroutine sym2(x,y,tau,px,py,ptau,dt)
      Implicit real*8 (A-H,O-Z)

      call f(x,y,tau,fx,fy,ftau)


      xnew= x+ px*dt +    fx*dt**2/2.d0

      ynew= y+ py*dt +    fy*dt**2/2.d0

      taunew= tau+ dt

      call f(xnew,ynew,taunew,fxnew,fynew,ftaunew)
      pxnew= px+ dt*(fx+fxnew )/2.d0
      pynew= py+ dt*(fy+fynew )/2.d0
      ptaunew= ptau+ dt*(ftau+ftaunew )/2.d0

      x=xnew
      y=ynew
      tau=taunew
      px=pxnew
      py=pynew
      ptau=ptaunew

      end


      subroutine sym4(x,y,tau,px,py,ptau,dt)
      Implicit real*8 (A-H,O-Z)
      sq2=2**(1.d0/3.d0)
      alpha= 1.d0/(2-sq2)
      beta= sq2/(2-sq2)
      dt1= dt*alpha
      dt2=-dt*beta
      call sym2(x,y,tau,px,py,ptau,dt1)
      call sym2(x,y,tau,px,py,ptau,dt2)
      call sym2(x,y,tau,px,py,ptau,dt1)
      return
      end

