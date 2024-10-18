c*************************************************************************
c                        KEP_UNIV.F
c*************************************************************************
c subroutine for solving kepler's equation using universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalor)
c                 q            ==>  Periastron
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c             Output:
c                 s             ==> Universal variable
c                 c0,sc1,s2c2,s3c3      ==>  s^i*c_i's from p171-172
c                                       (real scalors)
c                 iflg          ==>  =0 if converged; !=0 if not
c

      SUBROUTINE KEP_UNIV(DT,Q,MU,ALPHA,S,C0,SC1,S2C2,S3C3)
     
      USE DATA

      IMPLICIT NONE 

c...  Inputs: 
      real*8 dt,q,mu,alpha

c...  Outputs:
      real*8 s,c0,sc1,s2c2,s3c3
      integer iflg,itnew

c...  Internals:
      real*8 st,fo,fn

c----
c...  Executable code 

        call kepu_guess(dt,q,alpha,mu,s)
         
        st = s
c..     store initial guess for possible use later in
c..     laguerre's method, in case newton's method fails.

        itnew = 6
        call kepu_new(itnew,s,dt,q,mu,alpha,c0,sc1,s2c2,s3c3,iflg)
        if(iflg.ne.0) then
           call kepu_fchk(dt,q,mu,alpha,st,fo)
           call kepu_fchk(dt,q,mu,alpha,s,fn)
           if(abs(fo).lt.abs(fn)) then
               s = st 
           endif
           call kepu_lag(s,dt,q,mu,alpha,c0,sc1,s2c2,s3c3,iflg)
           if (iflg.ne.0) then
             call kepu_guess(dt,q,alpha,mu,s)
             itnew = 500
             call kepu_new(itnew,s,dt,q,mu,alpha,c0,sc1,s2c2,s3c3,iflg)
             if (iflg.ne.0) then
               WRITE(SD,*)'epic fail !'
               stop
             end if
           end if
        endif

        return
        end    ! kep_univ
c----------------------------------------------------------------------
c*************************************************************************
c                        KEPU_NEW.F
c*************************************************************************
c subroutine for solving kepler's equation in universal variables.
c using NEWTON'S METHOD
c
c             Input:
c                 itnew         ==> Max it's
c                 s             ==>  inital value of universal variable
c                 dt            ==>  time step (real scalor)
c                 q             ==>  Periastron (real scalor)
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c             Output:
c                 s             ==>  final value of universal variable
c                 fp            ==>  f' from p170  
c                                       (real scalors)
c                 c0,c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflgn          ==>  =0 if converged; !=0 if not
c
c  H Beust 03 dec 2014 : add itnew parameter

      subroutine kepu_new(itnew,s,dt,q,mu,alpha,c0,sc1,s2c2,s3c3,iflgn)

      IMPLICIT NONE 

c...  Inputs: 
      real*8 s,dt,q,mu,alpha
      integer itnew

c...  Outputs:
      real*8 fp,c0,sc1,s2c2,s3c3
      integer iflgn

c...  Internals:
      real*8, PARAMETER :: DANBYB = 1.0d-13
      integer nc
      real*8 x,ds,c1,c2,c3
      real*8 f,fpp,fppp,fdt

c----
c...  Executable code 

      do nc=0,itnew
         x = s*s*alpha
         call kepu_stumpff(x,c0,c1,c2,c3)
         sc1 = c1*s 
         s2c2 = c2*s*s 
         s3c3 = c3*s*s*s
         f = q*sc1 + mu*s3c3 - dt
         fp = q*c0 + mu*s2c2
         fpp = (-q*alpha + mu)*sc1
         fppp = (- q*alpha + mu)*c0
         ds = - f/fp
         ds = - f/(fp + ds*fpp/2.0d0)
         ds = -f/(fp + ds*fpp/2.0d0 + ds*ds*fppp/6.0d0)
         s = s + ds
         fdt = f/dt

c..      quartic convergence
         if ((fdt*fdt.lt.DANBYB*DANBYB).or.(abs(ds).lt.DANBYB)) then 
             iflgn = 0
             return
         endif
c...     newton's method succeeded

        enddo

c..     newton's method failed
        iflgn = 1
        return

        end  ! kepu_new


c----------------------------------------------------------------------
c subroutine for the calculation of stumpff functions
c
c             Input:
c                 x             ==>  argument
c             Output:
c                 c0,c1,c2,c3   ==>  c's from p171-172
c                                       (real scalors)

      subroutine kepu_stumpff(x,c0,c1,c2,c3)

      IMPLICIT NONE 

c...  Inputs: 
      real*8 x

c...  Outputs:
      real*8 c0,c1,c2,c3

c...  Internals:
      integer n,i
      real*8 xm

c----
c...  Executable code 

      n = 0
      xm = 0.1d0
      do while(abs(x).ge.xm)
         n = n + 1
         x = x/4.0d0
      enddo

      c2 = (1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x/182.d0)
     &       /132.d0)/90.d0)/56.d0)/30.d0)/12.d0)/2.d0
      c3 = (1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x/210.d0)
     &       /156.d0)/110.d0)/72.d0)/42.d0)/20.d0)/6.d0
      c1 = 1.d0 - x*c3
      c0 = 1.d0 - x*c2

      if(n.ne.0) then
         do i=n,1,-1
            c3 = (c2 + c0*c3)/4.d0
            c2 = c1*c1/2.d0
            c1 = c0*c1
            c0 = 2.d0*c0*c0 - 1.d0
            x = x * 4.d0
          enddo
       endif

       return
       end     !   kepu_stumpff

c------------------------------------------------------------------
c*************************************************************************
c Initial guess for solving kepler's equation using universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalor)
c                 q            ==>  periastron
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)

c             Output:
c                 s             ==>  initial guess for the value of 
c                                    universal variable
c

      subroutine kepu_guess(dt,q,alpha,mu,s)

      IMPLICIT NONE 

c...  Inputs: 
      real*8 dt,q,mu,alpha

c...  Inputs and Outputs:
      real*8 s

c...  Internals:
      integer iflg
      real*8 y,sy,cy,sigma,es
      real*8 x,a
      real*8 en,ec,e

c----
c...  Executable code 

        if (alpha.gt.0.0) then 
c...       find initial guess for elliptic motion

            if( dt/q .le. 0.4d0)  then
              s = dt/q
	      return
            else
              a = mu/alpha
              en = sqrt(mu/(a*a*a))
              e = 1.0d0 - q/a
              y = en*dt
              sigma = dsign(1.d0,e*sin(y))
              x = y + sigma*0.85d0*e
              s = x/sqrt(alpha)
	    endif

        else
c...       find initial guess for hyperbolic motion.
	   call kepu_p3solve(dt,q,mu,alpha,s,iflg)
	   if(iflg.ne.0) then
	      s = dt/q
	   endif
        endif

        return
        end     !   kepu_guess

c-------------------------------------------------------------------
c Returns the real root of cubic often found in solving kepler
c problem in universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalar)
c                 q            ==>  Periastron
c                 mu            ==>  Reduced mass of system (real scalar)
c                 alpha         ==>  -Twice the binding energy (real scalar)
c             Output:
c                 s             ==>  solution of cubic eqn for the  
c                                    universal variable
c                 iflg          ==>  success flag ( = 0 if O.K.) (integer)
c

      subroutine kepu_p3solve(dt,q,mu,alpha,s,iflg)

      IMPLICIT NONE 

c...  Inputs: 
      real*8 dt,q,mu,alpha

c...  Outputs:
      integer iflg
      real*8 s

c...  Internals:
      real*8 denom,a0,a1,qs,r,sq2,sq,p1,p2

c----
c...  Executable code 

	denom = - alpha*q/6.d0
	a1 = q/denom
	a0 =-dt/denom

	qs = a1/3.d0
	r = -a0/2.d0
	sq2 = qs**3 + r**2

	if( sq2 .ge. 0.d0) then
	   sq = sqrt(sq2)

	   if ((r+sq) .le. 0.d0) then
	      p1 =  -(-(r + sq))**(1.d0/3.d0)
	   else
	      p1 = (r + sq)**(1.d0/3.d0)
	   endif
	   if ((r-sq) .le. 0.d0) then
	      p2 =  -(-(r - sq))**(1.d0/3.d0)
	   else
	      p2 = (r - sq)**(1.d0/3.d0)
	   endif

	   iflg = 0
	   s = p1 + p2

	else
	   iflg = 1
	   s = 0.d0
	endif

        return
        end     !   kepu_p3solve
c-------------------------------------------------------------------
c*************************************************************************
c                        KEPU_FCHK.F
c*************************************************************************
c Returns the value of the function f of which we are trying to find the root
c in universal variables.
c
c             Input:
c                 dt            ==>  time step (real scalar)
c                 q            ==>  Periastron
c                 mu            ==>  Reduced mass of system (real scalar)
c                 alpha         ==>  Twice the binding energy (real scalar)
c                 s             ==>  Approx. root of f 
c             Output:
c                 f             ==>  function value ( = 0 if O.K.) (integer)
c

      subroutine kepu_fchk(dt,q,mu,alpha,s,f)

      IMPLICIT NONE 

c...  Inputs: 
      real*8 dt,q,mu,alpha,s

c...  Outputs:
      real*8 f

c...  Internals:
      real*8  x,c0,c1,c2,c3

c----
c...  Executable code 

        x=s*s*alpha
        call kepu_stumpff(x,c0,c1,c2,c3)
        f = q*s*c1 + mu*s*s*s*c3 - dt

        return
        end     !   kepu_fchk
c-------------------------------------------------------------------
c subroutine for solving kepler's equation in universal variables.
c using LAGUERRE'S METHOD
c
c             Input:
c                 s             ==>  inital value of universal variable
c                 dt            ==>  time step (real scalor)
c                 q            ==>  Periastron
c                 mu            ==>  Reduced mass of system (real scalor)
c                 alpha         ==>  energy (real scalor)
c             Output:
c                 s             ==>  final value of universal variable
c                 c1,c2,c3      ==>  c's from p171-172
c                                       (real scalors)
c                 iflgn          ==>  =0 if converged; !=0 if not
c
c  H. Beust 03 dec 2014 : fixed bug : -40.0 -> -q in fpp


      subroutine kepu_lag(s,dt,q,mu,alpha,c0,sc1,s2c2,s3c3,iflg)

      implicit none

c...  Inputs: 
      real*8 s,dt,q,mu,alpha

c...  Outputs:
      real*8 c0,sc1,s2c2,s3c3
      integer iflg

c...  Internals:
      integer, PARAMETER :: NLAG2 = 400
      real*8, PARAMETER :: DANBYB = 1.0d-13
      integer nc,ncmax
      real*8 ln
      real*8 x,fp,fpp,ds,f,c1,c2,c3
      real*8 fdt

      integer NTMP
      parameter(NTMP=NLAG2+1)

c----
c...  Executable code 

c...    To get close approch needed to take lots of iterations if alpha<0
        if(alpha.lt.0.0) then
           ncmax = NLAG2
        else
           ncmax = NLAG2
        endif

        ln = 5.0d0
c...    start laguere's method
        do nc =0,ncmax
           x = s*s*alpha
           call kepu_stumpff(x,c0,c1,c2,c3)
           sc1 = c1*s 
           s2c2 = c2*s*s 
           s3c3 = c3*s*s*s
           f = q*sc1 + mu*s3c3 - dt
           fp = q*c0 + mu*s2c2
           fpp = (-q*alpha + mu)*sc1
           ds = - ln*f/(fp + dsign(1.d0,fp)*sqrt(abs((ln - 1.0)*
     &       (ln - 1.0)*fp*fp - (ln - 1.0)*ln*f*fpp)))
           s = s + ds

           fdt = f/dt

c..        quartic convergence
           if ((fdt*fdt.lt.DANBYB*DANBYB).or.(abs(ds).lt.DANBYB)) then 
             iflg = 0
             return
           endif
c...      Laguerre's method succeeded
        enddo

        iflg = 2

        return

        end    !    kepu_lag
c-----------------------------------------------------------------------
