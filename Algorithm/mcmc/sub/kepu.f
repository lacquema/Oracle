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
C     
C -----------------------------------------------------------------------------
C     Compute MAP = p(data|orb)*prior(orb)
C                = prior(orb)*exp(-Chi2/2)*Product(1/sqrt(j^2+j0^2))
C     To avoid over/underflow problem, this function computes FMAP = ln(MAP)  
C-----------------------------------------------------------------------------
C
        REAL*8 FUNCTION COMPUTE_MAP(P,CHI2)

        USE DATA

        IMPLICIT NONE

        REAL*8, DIMENSION(NPLA) :: MTOT         ! Dynamical masses 
        REAL*8, DIMENSION(NPAR) :: P            ! Parameters
        REAL*8 ::       FMAP,           ! MAP probability
     &                  P1,P2,P3,P4,P5,P6,P7, ! Parameters
     &                  Q,              ! Periastron
     &                  UP2,            ! 1/Period^2
     &                  SI2,            ! sin(i/2)^2, then sin(i)^2
     &                  MUF,LMUF,       ! Prior condition & ln
     &                  SIGJV,SIGJV2,SIGV2P, ! Velocity Jitter
     &                  MSTAR,          ! Stellar mass
     &                  CHI2            ! CHI2
        INTEGER*4 ::    I,J,IV,DKAL
        
c...     MAP = exp(-Chi2/2)

        FMAP = -0.5d0*CHI2
        IF (MULTIPLA) MSTAR = EXP(P(NPAR))
        DO I = 1,NPLA
          DKAL = NEL*(I-1)
          P1 = P(DKAL+1)
          P2 = P(DKAL+2)
          P3 = P(DKAL+3)
          P4 = P(DKAL+4)
          P5 = P(DKAL+5)
          P6 = P(DKAL+6)
          P7 = P(DKAL+7)
          UP2 = EXP(-2.d0*P6)  !  1/Period^2
          Q = EXP(P1)          ! Period
          SI2 = (P4*P4+P5*P5) ! sin(i/2)^2
          SI2 = 4.d0*SI2*(1.d0-SI2)  ! sin(i)^2
          MTOT(I) = (Q**3)*DPI*DPI*UP2 ! n^2*q^3
          
          FMAP = FMAP+0.5d0*LOG(SI2*UP2)-P1 ! *sin(i)/P -> ln
c                                       -P1 => -ln(q) (*1/q, log prior)
        END DO

c...  Taking into account priors on masses in FMAP (+ln(prior(m)))    
        DO J = 0,NPRIOR            
           MUF = SUM(MPRIOR(J)%BCOF(1:NPLA)*MTOT(1:NPLA))
           IF (MULTIPLA) MUF = MUF+MPRIOR(J)%BCOF(0)*MSTAR
           CALL USER_PRIOR_MAP(J,FMAP,MUF)
        END DO

        IF (RADVEL.AND.(JITNUM.EQ.1)) THEN
           SIGJV = EXP(P(NEL*NPLA+2))
           SIGJV2 = SIGJV*SIGJV
c...  Taking into account prior on Stellar Jitter in FMAP (*1/(j+j0)) -> ln
           FMAP = FMAP-LOG(SIGJV+SIGV0)
c...  Taking into account effect of Stallar Jitter in p(data|model)
c...      = Product(1/sqrt(j^2+sig_i^2)) => -sum(1/2*ln(j^2+sig_i^2)
c...               => -sum(1/2*ln(1+j^2/sig_i^2)) (normalization)  
           FMAP = FMAP
     &        -0.5d0*SUM(LOG(1.d0+SIGJV2/STAR%SIGV2(1:STAR%NDATVR)))

        END IF
        
        COMPUTE_MAP = FMAP

        END

C
C -----------------------------------------------------------------------------
C       Just compute Chi^2, universal veriables version
C -----------------------------------------------------------------------------
C
        SUBROUTINE CHI2SEUL(P,CHI2)

        USE DATA

        IMPLICIT NONE

        REAL*8, DIMENSION(NPLA) ::
     &                  QQ,             ! Periastron
     &                  TP,             ! Tps de passage au periastre
     &                  NN,             ! Maen motion wrt q
     &                  GM,             ! Dynamical masses
     &                  EXC,            ! eccentricities
     &                  CW,SW,CP,SP,    ! cos,sin (omega+Omega,omega-Omega)
     &                  CI2,SI2,        ! cos^2(i/2), sin^2(i/2)
     &                  COM,SOM,CO,SO,  ! cos,sin (omega,Omega)
     &                  CI,SI,          ! cos(i), sin(i)
     &                  MUF,            ! Planet masses
     &                  MFRAC,          ! Fractional masses
     &                  AMX,AMY,        ! Partial radial velocity amplitudes 
     &                  AMPX,AMPY,      ! Same x MFRAC
     &                  POSX,POSY       ! Fitted positions
        REAL*8, DIMENSION(2,NPLA) ::
     &                  E1,E2           ! Vecteurs de la base propre
        REAL*8 ::       P(NPAR),        ! Parameters
     &                  V0,             ! Offset velocity
     &                  SIGJV,          ! Velocity Jitter
     &                  VRAD,           ! Radial velocity
     &                  SIGMA,          ! Cumulative mass
     &                  DX,DY,          ! Astrometric differences
     &                  DV,             ! Radial velocity difference   
     &                  CHI2            ! CHI2
        INTEGER*4       K,I
        
        CHI2 = 0.d0
        SIGJV = 0.d0
        CALL ELEMENTS(P(1:NPAR),NN,QQ,EXC,CW,SW,CI2,SI2,CP,SP,
     &                                CI,SI,COM,SOM,CO,SO,TP,MUF)
        MFRAC = 0.d0
        IF (MULTIPLA) THEN
           SIGMA = EXP(P(NPAR))
c...  Compute fractional masses
           DO I=1,NPLA
              SIGMA = SIGMA+MUF(I)
              GM(I) = SIGMA
              MFRAC(I) = MUF(I)/SIGMA
           END DO
        END IF
c...  Define vectors & amplitudes (formulas differ depending on whether
c...                             there are radial velocities or not)        
        IF (RADVEL) THEN
           IF (ISDATA(2)) THEN              
              V0 = P(NEL*NPLA+1)
              IF (JITNUM.EQ.1) SIGJV = EXP(P(NEL*NPLA+2))
           END IF
           E1(1,:) = COM*CO-SOM*CI*SO
           E1(2,:) = COM*SO+SOM*CI*CO
           E2(1,:) = -SOM*CO-COM*CI*SO
           E2(2,:) = -SOM*SO+COM*CO*CI
c           AMS = -NN*A*SI*SOM
c           AMC = +NN*A*SI*EXQ*COM
c           AMPS = AMS*MFRAC
c           AMPC = AMC*MFRAC
           AMX = -SIGMA*SOM
           AMY = SQRT(QQ*GM*(1.d0+EXC))*COM
           AMPX = AMX*MFRAC
           AMPY = AMY*MFRAC
        ELSE
           E1(1,:) = CI2*CW+SI2*CP
           E1(2,:) = CI2*SW-SI2*SP
           E2(1,:) = -CI2*SW-SI2*SP
           E2(2,:) = CI2*CW-SI2*CP           
        END IF

c...  Contribution of HCI (type 1 data) to Chi2        
        DO I = 1,NPLA
          DO K = 1,PLA(I)%NDATAS
             CALL POSFITS(PLA(I)%TAS(K),GM,
     &                    NN,QQ,EXC,E1,E2,TP,MFRAC,POSX,POSY)
             DX = (PLA(I)%X(K)-POSX(I))*PLA(I)%SIGXM1(K)
             DY = (PLA(I)%Y(K)-POSY(I))*PLA(I)%SIGYM1(K)       
             CHI2 = CHI2 +
     &      (DX*DX+DY*DY-2.d0*PLA(I)%RHOXY(K)*DX*DY)*PLA(I)%UFAC(K)
          END DO
        END DO

c...  Contribution of RV's (type 2 & 3 data) to Chi2.
c...        Note use of jitter (=0 if no jitter)
        IF (RADVEL) THEN
           IF (ISDATA(2)) THEN
              DO K = 1,STAR%NDATVR
                 CALL VFITS(STAR%TVR(K),GM,QQ,NN,V0,
     &                                      AMPX,AMPY,EXC,TP,VRAD)
                 DV = STAR%V(K)-VRAD
                 CHI2 = CHI2 + DV*DV/(STAR%SIGV2(K)+SIGJV*SIGJV)           
              END DO
           END IF
           IF (ISDATA(3)) THEN              
              DO I = 1,NPLA
                 DO K = 1,PLA(I)%NDATVR
                    CALL VFITSP(I,PLA(I)%TVR(K),GM,QQ,NN,AMPX,AMPY,
     &                              AMX(I),AMY(I),EXC,TP,VRAD)
                 END DO
                 DV = PLA(I)%V(K)-VRAD
                 CHI2 = CHI2 + DV*DV*PLA(I)%SIGVM2(K)
              END DO           
           END IF
        END IF
        
        END
      
C
C -----------------------------------------------------------------------------
C       Calculating orbital elements from a single MCMC orbital model
C -----------------------------------------------------------------------------
C

        SUBROUTINE ELEMENTS(P,NQ,Q,EXC,CW,SW,CI2,SI2,CP,SP,
     &                                CI,SI,COM,SOM,CO,SO,TP,MUF)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    DKAL,I
        REAL*8 ::       P(NPAR)         ! Parameters
        REAL*8, DIMENSION(NPLA) ::
     &       Q,                     ! Periastrons
     &       NQ,                    ! Mean motions wrt q
     &       EXC,                   ! Eccentricities
     &       CI,SI,                 ! cos,sin(inclinations)
     &       CW,SW,                 ! cos,sin(w=Omega+omega)
     &       CP,SP,                 ! cos,sin(phi=Oomega-Omega)
     &       COM,SOM,               ! cos,sin(omega's)
     &       SI2,CI2,               ! sin(i/2)^2, cos(i/2)^2
     &       CO,SO,                 ! cos,sin(Omega's)             
     &       TP,                    ! Times for periastron passages
     &       MUF                    ! Planetary masses
        REAL*8 ::
     &       PER,UP,                ! Period wrt q 
     &       MTOT,                  ! Total mass
     &       SII,                   ! sin(i/2)
     &       P1,P2,P3,P4,P5,P6,P7,  ! Intermédiaires
     &       S,                     ! Universal variable at T0
     &       ALPHA,XX,              ! -Energy
     &       SC1,S3C3,              ! Stumpff functions
     &       C0,C1,C2,C3,           !  "     "       "
     &       TT                     ! Time
        
        IF (MULTIPLA) THEN
           MTOT = EXP(P(NPAR))
        ELSE
           MTOT = STAR%MASS
        END IF
        DO I = 1,NPLA
           DKAL = NEL*(I-1)
           P1 = P(DKAL+1)
           P2 = P(DKAL+2)
           P3 = P(DKAL+3)
           P4 = P(DKAL+4)
           P5 = P(DKAL+5)
           P6 = P(DKAL+6)
           P7 = P(DKAL+7)
           Q(I) = EXP(P1)
           UP = EXP(-P6)
           NQ(I) = DPI*UP
           PER = 1.d0/UP
           MUF(I) = NQ(I)*NQ(I)*Q(I)**3 ! Dynamical mass of orbit #i
           EXC(I) = SQRT(P2*P2+P3*P3) ! e
           CW(I) = P2/EXC(I)                 ! w = Omega+omega 
           SW(I) = P3/EXC(I)
           SI2(I) = P4*P4+P5*P5               ! sin(i/2)^2 
           CI2(I) = 1.d0-SI2(I)               ! cos(i/2)^2
           SII = SQRT(SI2(I))      ! sin(i/2)
           IF (RADVEL) THEN
              CI(I) = 1.d0-2.d0*SI2(I)              ! cos(i)
              SI(I) = 2.d0*SII*SQRT(CI2(I))       ! sin(i)
              CO(I) = P4/SII  ! cos(Omega)
              SO(I) = P5/SII  ! sin(Omega)
              COM(I) = CW(I)*CO(I)+SW(I)*SO(I) ! omega = w-Omega
              SOM(I) = SW(I)*CO(I)-CW(I)*SO(I)         
              CP(I) = 0.d0
              SP(I) = 0.d0
           ELSE
              CP(I) = P4/SII   ! cos(phi)  phi=omega-Omega
              SP(I) = P5/SII   ! sin(phi)
              CI(I) = 0.d0
              SI(I) = 0.d0
              CO(I) = 0.d0
              SO(I) = 0.d0
              COM(I) = 0.d0
              SOM(I) = 0.d0
           END IF
           S = P7       ! Universal variable at T=T0
           ALPHA = MUF(I)/Q(I)*(1.d0-EXC(I))
           XX = S*S*ALPHA
           CALL KEPU_STUMPFF(XX,C0,C1,C2,C3)
           SC1 = C1*S 
           S3C3 = C3*S*S*S
           TT = MUF(I)*S3C3+Q(I)*SC1
           TP(I) = STAR%T0-TT
        END DO
        DO I = NPLA,2,-1
          MUF(I) = MUF(I)-MUF(I-1)           
        END DO
        IF (MULTIPLA) MUF(1) = MUF(1)-MTOT
        
        END    

C
C -----------------------------------------------------------------------------
C       Computation of priors ratio
C -----------------------------------------------------------------------------
C

        SUBROUTINE RATIO_PRIORS(P,PTRY,RAPQ,BAD)

        USE DATA

        IMPLICIT NONE
        
        LOGICAL ::      BAD          ! True if unreal solution
        INTEGER*4 ::    I,J,DKAL,K
        REAL*8, DIMENSION(NPAR) :: P, ! Parameters
     &                  PTRY         ! Trial model        
        REAL*8, DIMENSION(NPLA) :: MTOT,MTOTT ! Masses
        REAl*8 ::
     &       RAPQ,              ! q(x|x')/q(x'|q) (eq.42 Ford 06)
     &       P1,P2,P3,P4,P5,P6,P7, ! Intermédiaires
     &       PT1,PT2,PT3,PT4,PT5,PT6,PT7, !
     &       S,ST,              ! Variables universelles
     &       XX,ALPHA,FAC,      ! Intermédiaire
     &       DSDT,DSDTT,        ! ds/dt, ds'/dt
     &       C0,C1,C2,C3,S2C2,  ! Fonctions de Stumpff
     &       EXC,EXCT,          ! excentricite
     &       CW,SW,CWT,SWT,     ! w = Omega+omega
     &       SI2,               ! sin(i/2)^2
     &       Q,QT,              ! Periastrons   
     &       MUF,MUFT,          ! Fractional masses
     &       MSTAR, MSTART,     ! Stellar mass 
     &       SIGJV,SIGJVT,      ! Velocity Jitters
     &       UPQ,UPQT           ! Inverses de périodes

c     For each planet,
c           J/J'=(q/q')^2*(dsdt/dsdt')
c     Jabobian J (orbite->params) = 1/2*e*sin(i)*dsdt
c  Prior/planet (default) = sin(i)*1/q*1/P   (prior log in a & in P)
c     => Jacobian / Prior = 1/2*q*P*e*dsdt= (J/P)_0
c           N.B. : sin(i) vanishes in Jacobian / Prior
c     If prior u(mu) (e.g. Gaussian) then
c               mu = sum(b_i*M_i) = sum(b_i*n_i^2*q_i^3)
c     d(mu)/d(ln(q_i))=3*M_i*b_i   d(mu)/d(ln(P_i))=-2*M_i*b_i
c     Jacobian (mu<-q) = 3*(M_i*b_i)/q_i si autres variables inchangées
c     [Jacobian (mu<-P) =-2*(M_i*b_i)/P_i si autres variables inchangées]
c      With mu<->q :
c     Jacobian / Prior = (J/P)_0*prior(q)/J(q->mu)*1/prior(mu)
c                      = (J/P)_0*1/q*(q/(3*M_i*b_i))*1/u(mu)     
c     RAPQ = (J/p)/(J'/p')=(q/q')*(P/P')*(e/e')*(dsdt/dsdt')
c                                *(M'_i/M_i)*(u(mu')/u(mu))
c     Then RAPQ = RAPQ * J(jitter)/J'(jitter)*prior(jitter')/prior(jitter)
c                = RAPQ * (jitter'/jitter)*(jitter+j0)/(jitter'+j0)
c       Then, eq.(15) from Ford (2006). p(d|orb)) propto 1/(j^2+s_i^2)^(1/2)
c     => p(d|orb')|p(d|orb) = Product((j^2+s_i^2)/(j'^2+s_i^2))^(1/2)
c                              *exp(-(Chi2'-Chi2)/2)        
c
c        
        BAD = .FALSE.
        RAPQ = 1.d0
        IF (MULTIPLA) THEN
           MSTAR = EXP(P(NPAR))
           MSTART = EXP(PTRY(NPAR))
        END IF
        DO I = 1,NPLA
          DKAL = NEL*(I-1)
          P1 = P(DKAL+1)
          P2 = P(DKAL+2)
          P3 = P(DKAL+3)
          P4 = P(DKAL+4)
          P5 = P(DKAL+5)
          P6 = P(DKAL+6)
          P7 = P(DKAL+7)
          PT1 = PTRY(DKAL+1)
          PT2 = PTRY(DKAL+2)
          PT3 = PTRY(DKAL+3)
          PT4 = PTRY(DKAL+4)
          PT5 = PTRY(DKAL+5)
          PT6 = PTRY(DKAL+6)
          PT7 = PTRY(DKAL+7)
          Q = EXP(P1)
          QT = EXP(PT1)
          UPQ = EXP(-P6)
          UPQT = EXP(-PT6)
          S = P7
          ST = PT7
          MTOT(I) = Q**3*(DPI*UPQ)**2
          MTOTT(I) = QT**3*(DPI*UPQT)**2          
          EXC = SQRT(P2*P2+P3*P3)
          EXCT = SQRT(PT2*PT2+PT3*PT3)
          CW = P2/EXC
          SW = P3/EXC
          CWT = PT2/EXCT
          SWT = PT3/EXCT
c     BAD = BAD.OR.(VPRIOR(1)%BOUND(1)*UPT.GT.1.d0)
c     &             .OR.(VPRIOR(1)%BOUND(2)*UPT.LT.1.d0)
          BAD = BAD.OR.(QT.LT.VPRIOR(2)%BOUND(1))
     &             .OR.(QT.GT.VPRIOR(2)%BOUND(2))
          BAD = BAD.OR.(EXCT.GT.VPRIOR(3)%BOUND(2)).OR.
     &             (EXCT.LT.VPRIOR(3)%BOUND(1))
          SI2 = (PT4*PT4+PT5*PT5) ! sin(i/2)^2
          BAD = BAD.OR.(SI2.GT.1.d0)

          ALPHA = MTOT(I)/Q*(1.d0-EXC)
          XX = S*S*ALPHA
          CALL KEPU_STUMPFF(XX,C0,C1,C2,C3)
          S2C2 = C2*S*S 
          DSDT = 1.d0/(MTOT(I)*S2C2+Q*C0)

          ALPHA = MTOTT(I)/QT*(1.d0-EXCT)
          XX = ST*ST*ALPHA
          CALL KEPU_STUMPFF(XX,C0,C1,C2,C3)
          S2C2 = C2*ST*ST 
          DSDTT = 1.d0/(MTOTT(I)*S2C2+QT*C0)

          RAPQ = RAPQ*(Q/QT)*(UPQT/UPQ)*(EXC/EXCT)*(DSDT/DSDTT)
       END DO
c...  At this point RAPQ =
c...    Product((e/e')*(q/q')*(P/P')*(dsdt/dsdt'),planets)
c     Now consider contribution of mass priors.
c     For each prior u(mu)=u(sum(b_j*M_j))
c              RAPQ = RAPQ*u(mu')/u(mu)*product(M'_i/M_i), planets)
c                                      *m'_0/m_0 (if concerned)
       DO J = 0,NPRIOR            
           MUF = SUM(MPRIOR(J)%BCOF(1:NPLA)*MTOT(1:NPLA))
           MUFT = SUM(MPRIOR(J)%BCOF(1:NPLA)*MTOTT(1:NPLA))
           IF (MULTIPLA) THEN
              MUF = MUF+MPRIOR(J)%BCOF(0)*MSTAR
              MUFT = MUFT+MPRIOR(J)%BCOF(0)*MSTART
           END IF
c           print*,'jj',j,sngl(muf),sngl(muft)
           CALL USER_PRIOR_RATIO(J,FAC,MUF,MUFT)
           RAPQ = RAPQ*FAC
           DO I = 1,NPLA
              IF (MPRIOR(J)%BCOF(I).NE.0.d0)
     &              RAPQ = RAPQ*(MTOTT(I)/MTOT(I))
           END DO
           IF (MULTIPLA.AND.MPRIOR(J)%BCOF(0).NE.0.d0)
     &              RAPQ = RAPQ*(MSTART/MSTAR)      
        END DO
       
        IF (ISDATA(2).AND.(JITNUM.EQ.1)) THEN
c...  Contribution of velocity Jitter
           SIGJV = EXP(P(NEL*NPLA+2))
           SIGJVT = EXP(PTRY(NEL*NPLA+2))
           RAPQ = RAPQ*(SIGJVT/SIGJV)*(SIGJV+SIGV0)/(SIGJVT+SIGV0)
c...           Modified Jeffrey's prior for sigv propto 1/(sigv+sigv0)        
c...                  J/J'=Jeffrey'/Jeffrey=(sigv+sigv0)/(sigv'+sigv0)        
c...           Then variable is ln(sigv') => J/J'*sigv'/sigv
           S = 1.d0
           DO K = 1,STAR%NDATVR
              S = S*(STAR%SIGV2(K)+SIGJV*SIGJV)
     &             /(STAR%SIGV2(K)+SIGJVT*SIGJVT)
           END DO   
c...  Eq.(15) From Ford (2006). For each radv data, multiply ratio by
c...              sqrt(sigvk^2+sigv^2)/sqrt(sigvk^2+sigv'^2)        
c...          (resulting from Eqs. (5) and (12))
           RAPQ = RAPQ*SQRT(S)      
        END IF

        END         



C -----------------------------------------------------------------------------
C       Fit of the velocity without derivatives
C -----------------------------------------------------------------------------
C     

        SUBROUTINE VFITS(TT,GM,Q,NQ,V0,AMPX,AMPY,EXC,TP,VRAD)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 :: I
        REAL*8, DIMENSION(NPLA) ::
     &                  TP,             ! Time for periastron passage
     &                  Q,              ! Periastron
     &                  NQ,             ! Mean motion wrt q
     &                  GM,             ! Dynamical masses
     &                  EXC,            ! Eccentricity
     &                  AMPX,AMPY       ! Partial amplitudes
        REAL*8 ::       TT,DT,          ! Time
     &                  VRAD,           ! Vitesse radiale
     &                  UEXC,           ! 1-e
     &                  ALPHA,          ! -2*Energy
     &                  V0,             ! Offset velocity
     &                  PER,            ! Period
     &                  S,DSDT,         ! Universal variable & derivative
     &                  C0,SC1,S2C2,S3C3 ! Stumpff functions
c
c     For each planet, AMPX = -GM*mfrac*sin(omega)
c                      AMPY = +sqrt(q*GM*(1+e))*mfrac*cos(omega)         
c     vrad = V0+sum((ampx(i)*s(i)*c1(i)+ampy(i)*c0(i))*(ds/dt(i)))
c             ds/dt(i) = 1/r(i) = 1/(GM(i)*s(i)^2*c2(i)+q(i)*c0(i))
        
        VRAD = V0
        DO I = 1,NPLA
          DT = TT-TP(I)
          UEXC = 1.d0-EXC(I)
          ALPHA = GM(I)/Q(I)*UEXC
          IF (EXC(I).LT.1.d0) THEN
            PER = (DPI/NQ(I))/(UEXC*SQRT(UEXC))
            DT = DT-FLOOR(DT/PER+0.5d0)*PER
          END IF
          CALL KEP_UNIV(DT,Q(I),GM(I),ALPHA,S,C0,SC1,S2C2,S3C3)
          DSDT = 1.d0/(GM(I)*S2C2+Q(I)*C0) ! ds/dt = 1/r     
          VRAD = VRAD+(AMPX(I)*SC1+AMPY(I)*C0)*DSDT
        END DO

        END
C
C
C -----------------------------------------------------------------------------
C       Fit of relative velocity of a planet / star without derivatives
C -----------------------------------------------------------------------------
C     

        SUBROUTINE VFITSP(IPLA,TT,GM,Q,N,AMPX,AMPY,AMX,AMY,EXC,TP,VRAD)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    IPLA,            ! PLanet to be considered
     &                  I
        REAL*8, DIMENSION(NPLA) ::
     &                  TP,             ! Time for periastron passage
     &                  Q,              ! Periastron
     &                  N,              ! Mean motion wrt q
     &                  GM,             ! Dynamical masses
     &                  EXC,            ! Eccentricity
     &                  AMPX,AMPY       ! Partial amplitudes
        REAL*8 ::       TT,DT,          ! Time
     &                  VRAD,           ! Vitesse radiale
     &                  AMX,AMY,        ! Amplitudes without mass
     &                  UEXC,           ! 1-e
     &                  ALPHA,          ! -2*Energy
     &                  V0,             ! Offset velocity
     &                  PER,            ! Period
     &                  S,DSDT,         ! Universal variable & derivative
     &                  C0,SC1,S2C2,S3C3 ! Stumpff functions
c
c     For each planet, AMPX = -GM*mfrac*sin(omega)
c                      AMPY = +sqrt(q*GM*(1+e))*mfrac*cos(omega)         
c     vrad = -sum((ampx(i)*s(i)*c1(i)+ampy(i)*c0(i))*(ds/dt(i)),i=1..ipla-1)
c             -(amx(ipla)*s(ipla)*c1(ipla)+amy(ipla)*c0(pla))*(ds/dt(ipla))
c     ds/dt(i) = 1/r(i) = 1/(GM(i)*s(i)^2*c2(i)+q(i)*c0(i))
c                AMX, AMY = Same as AMPX,AMPY without mfrac
        VRAD = 0.d0
        DO I = 1,IPLA-1
          DT = TT-TP(I)
          UEXC = 1.d0-EXC(I)
          ALPHA = GM(I)/Q(I)*UEXC
          IF (EXC(I).LT.1.d0) THEN
            PER = (DPI/N(I))/(UEXC*SQRT(UEXC))
            DT = DT-FLOOR(DT/PER+0.5d0)*PER
          END IF
          CALL KEP_UNIV(DT,Q(I),GM(I),ALPHA,S,C0,SC1,S2C2,S3C3)
          DSDT = 1.d0/(GM(I)*S2C2+Q(I)*C0) ! ds/dt = 1/r     
          VRAD = VRAD-(AMPX(I)*SC1+AMPY(I)*C0)*DSDT
        END DO
        DT = TT-TP(IPLA)
        UEXC = 1.d0-EXC(IPLA)
        ALPHA = GM(IPLA)/Q(IPLA)*UEXC
        IF (EXC(IPLA).LT.1.d0) THEN
           PER = (DPI/N(IPLA))/(UEXC*SQRT(UEXC))
           DT = DT-FLOOR(DT/PER+0.5d0)*PER
        END IF
        CALL KEP_UNIV(DT,Q(IPLA),GM(IPLA),ALPHA,S,C0,SC1,S2C2,S3C3)
        DSDT = 1.d0/(GM(IPLA)*S2C2+Q(IPLA)*C0) ! ds/dt = 1/r     
        VRAD = VRAD-(AMX*SC1+AMY*C0)*DSDT
        END

C
C -----------------------------------------------------------------------------
C       Fit of the position without derivatives
C -----------------------------------------------------------------------------
C

        SUBROUTINE POSFITS(TT,GM,N,QQ,EXC,E1,E2,TP,MFRAC,POSX,POSY)

        USE DATA

        IMPLICIT NONE

        REAl*8, DIMENSION(NPLA) ::
     &       JACX,JACY,         ! Jacobi positions
     &       QQ,                ! Periastron
     &       TP,                ! Time for periastron passage
     &       GM,                ! Dynamical masses
     &       MFRAC,             ! Fractional planetary masses
     &       N,                 ! Mean motion wrt q
     &       POSX,POSY,         ! Heliocentric positions
     &       EXC                ! Eccentricities
        REAL*8, DIMENSION(2,NPLA) :: E1,E2 ! Base vectors for orbits 
        REAL*8 ::
     &       TT,DT,             ! Time
     &       PER,               ! Orbital period
     &       UEXC,              ! 1-e
     &       ALPHA,             ! -2*Energie (s=u/sqrt(|alpha|))
     &       C0,SC1,S2C2,S3C3,  ! s^i*c_i(alpha*s^2), i=1,2,3
     &       S,                 ! Variable universelle
     &       GX,GY,             ! Partial center of mass position
     &       RX,RY              ! Coordonnes dans la base propre
        INTEGER*4 ::    I
        
        GX = 0.d0
        GY = 0.d0
        DO I = 1,NPLA
           DT = TT-TP(I)
           UEXC = 1.d0-EXC(I)
           ALPHA = GM(I)/QQ(I)*UEXC
           IF (EXC(I).LT.1.d0) THEN
              PER = (DPI/N(I))/(UEXC*SQRT(UEXC))
              DT = DT-FLOOR(DT/PER+0.5d0)*PER
           END IF
           CALL KEP_UNIV(DT,QQ(I),GM(I),ALPHA,S,C0,SC1,S2C2,S3C3)
           RX = QQ(I)-GM(I)*S2C2
           RY = SQRT(QQ(I)*GM(I)*(1.d0+EXC(I)))*SC1
           JACX(I) = RX*E1(1,I)+RY*E2(1,I)
           JACY(I) = RX*E1(2,I)+RY*E2(2,I)
           POSX(I) = JACX(I)+GX
           POSY(I) = JACY(I)+GY
           IF (MULTIPLA) THEN
             GX = GX+MFRAC(I)*JACX(I)
             GY = GY+MFRAC(I)*JACY(I)
           END IF
        END DO

        END
C
C
C -----------------------------------------------------------------------------
C     Routine to enter initial data
C -----------------------------------------------------------------------------
C

        SUBROUTINE INIT_DATA(N,P,COV,FILES,NEW)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 :: N,I,DKAL
        REAL*8, DIMENSION(N) :: P ! Initial Parameters
        REAL*8, DIMENSION(N,N) :: COV ! Initial covariances
        REAL*8 ::
     &       M0,                ! Anomalie moyenne de reference
     &       MASS,              ! Dynamical mass
     &       S,                 ! Universal variable   
     &       HH2,               ! -2xEnergie = alpha de kepu
     &       C0,SC1,S2C2,S3C3   ! Fonctions de Stumpff

        LOGICAL ::      OK,NEW
        CHARACTER*(*), DIMENSION(NFIL) :: FILES ! #1 = output, #2 = dump,
c                                                #3 = data

        STAR%T0 = 0.d0
        STAR%SIGJV = 0.d0
        CALL READ_DATAFILE(FILES(3))
        NFREE = 2*SUM(PLA%NDATAS)+STAR%NDATVR
     &          +SUM(PLA%NDATVR)+2*STAR%NDATAS-NPAR  ! # of degrees od freedom

        IF (NEW) THEN
           WRITE(SD,*)'Give reference time for data (JD) : '
           READ(5,*)STAR%T0
           IF (ISDATA(2)) THEN
c...  Initialize offset velocity & jitter
              WRITE(SD,*)
     &              'Give initial offset velocity and jitter in km/s '
              READ(5,*)STAR%V0,STAR%SIGJV
              STAR%SIGJV = STAR%SIGJV*MPS
              STAR%V0 = STAR%V0*MPS
           END IF             
        END IF

c...  Enter orbits and configure parameters for Levenberg-Marquardt        
        IF (NEW) THEN
           CALL ENTER_ORBITS(NPLA)
           DO I = 1,NPLA
              DKAL = NEL*(I-1)
              P(DKAL+1) = LOG(PLA(I)%Q) 
              P(DKAL+2) = PLA(I)%EXC*PLA(I)%CW
              P(DKAL+3) = PLA(I)%EXC*PLA(I)%SW
              P(DKAL+4) = PLA(I)%TI*PLA(I)%CP
              P(DKAL+5) = PLA(I)%TI*PLA(I)%SP
              IF (RADVEL) THEN
                 P(DKAL+4) = PLA(I)%TI*PLA(I)%CO
                 P(DKAL+5) = PLA(I)%TI*PLA(I)%SO
              END IF
              P(DKAL+6) = LOG(PLA(I)%PER)
              MASS = (DPI/PLA(I)%PER)**2*PLA(I)%Q**3
              HH2 = -(PLA(I)%EXC-1.d0)*MASS/PLA(I)%Q
c...                                hh2 = -2xEnergie = alpha pour kepu          
              M0 =  STAR%T0-PLA(I)%TP
              CALL KEP_UNIV(M0,PLA(I)%Q,MASS,HH2,S,C0,SC1,S2C2,S3C3)
              P(DKAL+7) = S
           END DO  
           IF (MULTIPLA) P(NPAR) = LOG(STAR%MASS)
           IF (ISDATA(2)) P(NPAR-1) = STAR%V0
           IF (ISDATA(2).AND.(JITNUM.EQ.1)) THEN
              P(NPAR-1) = LOG(STAR%MASS)
c...                     NPAR-1 because no Jitter at first round             
              P(NPAR-2) = STAR%V0
           END IF
c           WRITE(SD,*)'Range of periods : (days)'
c           READ(5,*)VPRIOR(1)%BOUND(1:2)
           WRITE(SD,*)'Range of periastra : (au)'
           READ(5,*)VPRIOR(2)%BOUND(1:2)
           WRITE(SD,*)'Range of eccentricities : '
           READ(5,*)VPRIOR(3)%BOUND(1:2)   
           VPRIOR(1:NLIM)%MEAN = 0.d0
           VPRIOR(1:NLIM)%SDEV = 0.d0
           VPRIOR(1:2)%TYP = 0
           VPRIOR(3)%TYP = 3
        END IF

        END


C      
C -----------------------------------------------------------------------------
C       Fit of the position and derivatives (for Levenberg-Marquardt)
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQPOSFIT(N1,TT,P,POSX,POSY,DPOSX,DPOSY)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    N1               ! Dimension of parameter space 
        REAL*8, DIMENSION(N1) :: P       ! Parameters 
        REAL*8, DIMENSION(NPLA) :: POSX,POSY ! Heliocentric positions
        REAL*8, DIMENSION(N1,NPLA) :: DPOSX,DPOSY ! Derivatives / parameters
        REAL*8, DIMENSION(NPLA) :: JACX,JACY ! Jacobi pos
        REAL*8, DIMENSION(7,NPLA) :: DJACX,DJACY ! Jac. deriv.
        REAL*8, DIMENSION(NPLA) :: MDYN,PER ! Dynamical masses & Periods / q
        REAL*8, DIMENSION(NPLA) :: MFRAC ! Fractional masses
        INTEGER*4 ::    I,J,L,DKAL ! Indexes
        REAL*8 ::
     &       TT,DT,             ! Temps
     &       MTOT,MSTAR,        ! Total mass up to planet #i
     &       Z,Q,               ! Periastron & ln()
     &       TP,                ! Tps de passage au periastre
     &       N,                 ! Mean motion / q
     &       K,H,               ! e*[cos(w),sin(w)]
     &       QQ,PP,             ! tan(i/2)*cos(phi),tan(i/2)*sin(phi)
     &       DJACXDZ,DJACXDLP,DJACXDS0,DJACXDK,DJACXDH,
     &       DJACXDPP,DJACXDQQ,DJACXDQ,DJACXDPER,
     &       DJACYDZ,DJACYDLP,DJACYDS0,DJACYDK,DJACYDH,
     &       DJACYDPP,DJACYDQQ,DJACYDQ,DJACYDPER,   ! All derivatives
     &       S,S0,              ! Universal variable (at T and T0)
     &       TI,T2,              ! tan(i/2) + square
     &       EXC,               ! eccentricity
     &       CW,SW,CP,SP,CO,SO, ! cos,sin (w,phi,O)
     &       CI,SI,COM,SOM,     ! cos,sin (i,Om)
     &       ALPHA,XX,UFA,      ! Energy
     &       CI2,SI2,           ! cos(i/2)^2, sin(i/2)^2
     &       E1(2),E2(2),       ! Vecteurs de la base propre (X,Y)
     &       U,CU,SU,           ! Anomalie excentrique
     &       FACT,              ! 1/(1-e*cos(u))
     &       RX,RY,             ! Coordonnes dans la base propre
     &       C0,C1,C2,C3,SC1,S2C2,S3C3,    ! Stumpff functions
     &       DTPDS0,DTPDALPHA,DTPDM,DTPDQ, ! Derivatives of TP
     &       DSDTP,DSDQ,DSDALPHA,DSDE,DSDPER,DSDM, ! Derivatives of S
     &       DOMDK,DOMDH,DODP,DODQ,DIDP,DIDQ,DOMDP,DOMDQ, ! Intern derivs. 
     &       DWDK,DWDH,DPHIDP,DPHIDQ, ! Derivees interm.
     &       DCI2DP,DCI2DQ,DSI2DP,DSI2DQ, !  
     &       DMFDM,DMFDMM,DUDH,DUDK, !
     &       DE1DK(2),DE1DH(2),DE1DP(2),DE1DQ(2), ! Derivees de E1
     &       DE2DK(2),DE2DH(2),DE2DP(2),DE2DQ(2), ! Derivees de E2
     &       DEDK,DEDH,         ! Derivees de Exc
     &       DRXDQ,DRXDS,DRXDE,DRXDS0,DRXDTP, ! Derivees de RX (internal)
     &       DRXDPER,DRXDALPHA,DRXDK,DRXDH,DRXDPP,DRXDQQ,
     &       DRYDQ,DRYDS,DRYDE,DRYDS0,DRYDTP, ! Derivees de RY (internal)
     &       DRYDPER,DRYDALPHA,DRYDK,DRYDH,DRYDPP,DRYDQQ
        
        POSX = 0.d0
        POSY = 0.d0
        DPOSX = 0.d0
        DPOSY = 0.d0
        IF (MULTIPLA) MSTAR = EXP(P(N1))
        MTOT = MSTAR
        DO I = 1,NPLA
          DKAL = NEL*(I-1)
          Z = P(DKAL+1)
          K = P(DKAL+2)
          H = P(DKAL+3)
          QQ = P(DKAL+4)
          PP = P(DKAL+5)
          PER(I) = EXP(P(DKAL+6))
          S0 = P(DKAL+7)
          N = DPI/PER(I)         ! Mean motion / q
          Q = EXP(Z)            ! Semi-major axis
          MDYN(I) = N*N*(Q**3)    ! Dynamical mass
          IF (MULTIPLA) THEN
             MFRAC(I) = (MDYN(I)-MTOT)/MDYN(I) ! Mfrac(i) = (M(i)-M(i-1))/M(i)
             MTOT = MDYN(I)
          END IF
          EXC = SQRT(K*K+H*H) ! e
          CW = K/EXC            ! cos(w)
          SW = H/EXC                   ! sin(w)
          ALPHA = MDYN(I)/Q*(1.d0-EXC) ! -2*Energy
          XX = S0*S0*ALPHA
          CALL KEPU_STUMPFF(XX,C0,C1,C2,C3)
          SC1 = C1*S0 
          S2C2 = C2*S0*S0
          S3C3 = C3*S0*S0*S0
          DT = MDYN(I)*S3C3+Q*SC1  ! Kepler's equation
          TP = STAR%T0-DT
          DTPDS0 = -MDYN(I)*S2C2-Q*C0
          DTPDALPHA = -(0.5d0/ALPHA)*(MDYN(I)*(S0*S2C2-3.d0*S3C3)
     &          +Q*(S0*C0-SC1))
          DTPDM = -S3C3
          DTPDQ = -SC1

          DT = TT-TP
          CALL KEP_UNIV(DT,Q,MDYN(I),ALPHA,S,C0,SC1,S2C2,S3C3)
          RX = Q-MDYN(I)*S2C2
          UFA = SQRT(Q*MDYN(I)*(1.d0+EXC))
          RY = UFA*SC1
          FACT = 1.d0/(MDYN(I)*S2C2+Q*C0)
          DSDTP = -FACT
          DSDQ = -SC1*FACT
          DSDM = -S3C3*FACT
          DSDALPHA = -0.5d0*FACT/ALPHA*
     &                  (MDYN(I)*(S*S2C2-3.d0*S3C3)+Q*(S*C0-SC1))
          DSDE = -DSDALPHA*MDYN(I)/Q
          DSDQ = 2.d0*DSDALPHA*ALPHA/Q+DSDQ
          DSDPER = -2.d0*DSDALPHA*ALPHA/PER(I)
          DSDQ = DSDQ+DSDM*3.d0*MDYN(I)/Q
          DSDPER = DSDPER-DSDM*2.d0*MDYN(I)/PER(I)
          
          DRXDQ = 1.d0-3.d0*MDYN(I)*S2C2/Q
          DRXDS = -MDYN(I)*SC1
          DRXDALPHA = -0.5d0*MDYN(I)/ALPHA*(S*SC1-2.d0*S2C2)
          DRXDPER = +2.d0*MDYN(I)*S2C2/PER(I)
          
          DRYDQ = 2.d0*RY/Q   
          DRYDPER = -RY/PER(I)
          DRYDE = 0.5d0*RY/(1.d0+EXC)
          DRYDS = UFA*C0
          DRYDALPHA = 0.5d0*UFA/ALPHA*(S*C0-SC1)
        
          DRXDTP = DRXDS*DSDTP
          DRXDQ = DRXDQ+DRXDS*DSDQ+2.d0*DRXDALPHA*ALPHA/Q
          DRXDPER = DRXDPER+DRXDS*DSDPER-2.d0*DRXDALPHA*ALPHA/PER(I) 
          DRXDE = -DRXDALPHA*MDYN(I)/Q+DRXDS*DSDE

          DRYDTP = DRYDS*DSDTP
          DRYDQ = DRYDQ+DRYDS*DSDQ+2.d0*DRYDALPHA*ALPHA/Q
          DRYDPER = DRYDPER+DRYDS*DSDPER-2.d0*DRYDALPHA*ALPHA/PER(I)
          DRYDE = DRYDE+DRYDS*DSDE-DRYDALPHA*MDYN(I)/Q
          
          DRXDQ = DRXDQ+DRXDTP*(DTPDQ+3.d0*DTPDM*MDYN(I)/Q)
          DRXDPER = DRXDPER-2.d0*DRXDTP*DTPDM*MDYN(I)/PER(I)
          DRXDQ = DRXDQ+2.d0*DRXDTP*DTPDALPHA*ALPHA/Q
          DRXDPER = DRXDPER-2.d0*DRXDTP*DTPDALPHA*ALPHA/PER(I)
          DRXDE = DRXDE-DRXDTP*DTPDALPHA*MDYN(I)/Q
          DRXDS0 = DRXDTP*DTPDS0
          
          DRYDQ = DRYDQ+DRYDTP*(DTPDQ+3.d0*DTPDM*MDYN(I)/Q)
          DRYDPER = DRYDPER-2.d0*DRYDTP*DTPDM*MDYN(I)/PER(I)         
          DRYDQ = DRYDQ+2.d0*DRYDTP*DTPDALPHA*ALPHA/Q
          DRYDPER = DRYDPER-2.d0*DRYDTP*DTPDALPHA*ALPHA/PER(I)
          DRYDE = DRYDE-DRYDTP*DTPDALPHA*MDYN(I)/Q
          DRYDS0 = DRYDTP*DTPDS0

          DEDK = CW
          DEDH = SW
          DWDK = -SW/EXC
          DWDH = CW/EXC

          T2 = QQ*QQ+PP*PP        
          TI = SQRT(T2)        ! Tan(i/2)
          IF (RADVEL) THEN
             SI = 2.d0*TI/(1.d0+T2)      ! sin(i) 
             CI = (1.d0-T2)/(1.d0+T2)  ! cos(i)
             CO = QQ/TI         ! cos(Omega)
             SO = PP/TI         ! sin(Omega)
             COM = CW*CO+SW*SO  ! omega = w-Omega
             SOM = SW*CO-CW*SO  !
          
             DOMDK = -SW/EXC ! dw/dk
             DOMDH = CW/EXC  ! dw/dh
             DODQ = -PP/T2
             DODP = QQ/T2
             DIDP = 2.d0*SO/(1.d0+T2)
             DIDQ = 2.d0*CO/(1.d0+T2)
             DOMDQ = PP/T2 ! -dO/dq
             DOMDP = -QQ/T2 ! -dO/dp
             E1 = (/ COM*CO-SOM*CI*SO , COM*SO+SOM*CI*CO /)
             E2 = (/ -SOM*CO-COM*CI*SO, -SOM*SO+COM*CO*CI /)

             DE1DK = (/ E2(1)*DOMDK, E2(2)*DOMDK /)
             DE1DH = (/ E2(1)*DOMDH, E2(2)*DOMDH /)
             DE2DK = (/ -E1(1)*DOMDK, -E1(2)*DOMDK /)
             DE2DH = (/ -E1(1)*DOMDH, -E1(2)*DOMDH /)     

             DE1DP = (/ E2(1)*DOMDP-E1(2)*DODP+SOM*SI*SO*DIDP, 
     &                E2(2)*DOMDP+E1(1)*DODP-SOM*CO*SI*DIDP /)
             DE1DQ = (/ E2(1)*DOMDQ-E1(2)*DODQ+SOM*SI*SO*DIDQ, 
     &                E2(2)*DOMDQ+E1(1)*DODQ-SOM*CO*SI*DIDQ /)
             DE2DP = (/ -E1(1)*DOMDP-E2(2)*DODP+COM*SI*SO*DIDP, 
     &               -E1(2)*DOMDP+E2(1)*DODP-COM*CO*SI*DIDP /)
             DE2DQ = (/ -E1(1)*DOMDQ-E2(2)*DODQ+COM*SI*SO*DIDQ, 
     &               -E1(2)*DOMDQ+E2(1)*DODQ-COM*CO*SI*DIDQ /)       
          ELSE
             CI2 = 1.d0/(1.d0+T2)   ! cos(i/2)^2
             SI2 = T2*CI2       ! sin(i/2)^2
             CP = QQ/TI         ! cos(phi)
             SP = PP/TI         ! sin(phi)
             DPHIDQ = -PP/T2
             DPHIDP = QQ/T2
             DCI2DP = -2.d0*CI2*CI2*PP
             DCI2DQ = -2.d0*CI2*CI2*QQ
             DSI2DP = -DCI2DP
             DSI2DQ = -DCI2DQ
                
             E1 = (/ CI2*CW+SI2*CP, CI2*SW-SI2*SP /)
             E2 = (/ -CI2*SW-SI2*SP, CI2*CW-SI2*CP /)

             DE1DK = (/ -CI2*SW*DWDK, CI2*CW*DWDK /)
             DE1DH = (/ -CI2*SW*DWDH, CI2*CW*DWDH /)
             DE2DK = (/ -CI2*CW*DWDK, -CI2*SW*DWDK /)
             DE2DH = (/ -CI2*CW*DWDH, -CI2*SW*DWDH /)
       
             DE1DP = (/ DCI2DP*CW+DSI2DP*CP-SI2*SP*DPHIDP,
     &             DCI2DP*SW-DSI2DP*SP-SI2*CP*DPHIDP /)
             DE1DQ = (/ DCI2DQ*CW+DSI2DQ*CP-SI2*SP*DPHIDQ,
     &             DCI2DQ*SW-DSI2DQ*SP-SI2*CP*DPHIDQ /)
             DE2DP = (/ -DCI2DP*SW-DSI2DP*SP-SI2*CP*DPHIDP , 
     &             DCI2DP*CW-DSI2DP*CP+SI2*SP*DPHIDP /)
             DE2DQ = (/ -DCI2DQ*SW-DSI2DQ*SP-SI2*CP*DPHIDQ , 
     &             DCI2DQ*CW-DSI2DQ*CP+SI2*SP*DPHIDQ /)
          END IF

          JACX(I) = RX*E1(1)+RY*E2(1)
          JACY(I) = RX*E1(2)+RY*E2(2)
          
          DRXDK = DRXDE*DEDK
          DRXDH = DRXDE*DEDH
          DRYDK = DRYDE*DEDK
          DRYDH = DRYDE*DEDH
          DRXDPP = 0.d0
          DRXDQQ = 0.d0
          DRYDPP = 0.d0
          DRYDQQ = 0.d0           

          DJACXDQ = DRXDQ*E1(1)+DRYDQ*E2(1)
          DJACXDK = DRXDK*E1(1)+DRYDK*E2(1)+RX*DE1DK(1)+RY*DE2DK(1)
          DJACXDH = DRXDH*E1(1)+DRYDH*E2(1)+RX*DE1DH(1)+RY*DE2DH(1)
          DJACXDQQ = DRXDQQ*E1(1)+DRYDQQ*E2(1)+RX*DE1DQ(1)+RY*DE2DQ(1)
          DJACXDPP = DRXDPP*E1(1)+DRYDPP*E2(1)+RX*DE1DP(1)+RY*DE2DP(1)
          DJACXDPER = DRXDPER*E1(1)+DRYDPER*E2(1)
          DJACXDS0 = DRXDS0*E1(1)+DRYDS0*E2(1)
          
          DJACYDQ = DRXDQ*E1(2)+DRYDQ*E2(2)
          DJACYDK = DRXDK*E1(2)+DRYDK*E2(2)+RX*DE1DK(2)+RY*DE2DK(2)
          DJACYDH = DRXDH*E1(2)+DRYDH*E2(2)+RX*DE1DH(2)+RY*DE2DH(2)
          DJACYDQQ = DRXDQQ*E1(2)+DRYDQQ*E2(2)+RX*DE1DQ(2)+RY*DE2DQ(2)
          DJACYDPP = DRXDPP*E1(2)+DRYDPP*E2(2)+RX*DE1DP(2)+RY*DE2DP(2)
          DJACYDPER = DRXDPER*E1(2)+DRYDPER*E2(2)
          DJACYDS0 = DRXDS0*E1(2)+DRYDS0*E2(2)

          DJACXDZ = Q*DJACXDQ                                  
          DJACYDZ = Q*DJACYDQ 
          DJACXDLP = PER(I)*DJACXDPER
          DJACYDLP = PER(I)*DJACYDPER

          DJACX(:,I) =
     & (/ DJACXDZ,DJACXDK,DJACXDH,DJACXDQQ,DJACXDPP,DJACXDLP,DJACXDS0 /)
          DJACY(:,I) =
     & (/ DJACYDZ,DJACYDK,DJACYDH,DJACYDQQ,DJACYDPP,DJACYDLP,DJACYDS0 /)
        END DO
       
c        print*,'mdyn',sngl(mdyn),'mfrac',sngl(mfrac)
        DO I = 1,NPLA
c... Basic formula Jacobi coords => helio positions
           POSX(I) = JACX(I)+SUM(MFRAC(1:I-1)*JACX(1:I-1))
           POSY(I) = JACY(I)+SUM(MFRAC(1:I-1)*JACY(1:I-1))
c...  Derivatives
           DKAL = NEL*(I-1)
           DPOSX(DKAL+1:DKAL+7,I) = DJACX(:,I)
           DPOSY(DKAL+1:DKAL+7,I) = DJACY(:,I)
           DO J = 1,I-1
              DKAL = NEL*(J-1)   ! Derivatives relative to JACX/Y(J)
              DPOSX(DKAL+1:DKAL+7,I) = DPOSX(DKAL+1:DKAL+7,I)
     &                +MFRAC(J)*DJACX(:,J)
              DPOSY(DKAL+1:DKAL+7,I) = DPOSY(DKAL+1:DKAL+7,I)
     &                +MFRAC(J)*DJACY(:,J)
              IF (J.GT.1) THEN   ! Derivatives relative to MFRAC(J)
                 MTOT = MDYN(J-1)  ! MFRAC(J) = (MDYN(J)-MTOT)/MDYN(J)
              ELSE
                 MTOT = MSTAR
              END IF
              DMFDM = MTOT/MDYN(J) ! d(MFRAC(J))/d(ln(MDYN(J))
              DPOSX(DKAL+1,I) = DPOSX(DKAL+1,I)+3.d0*DMFDM*JACX(J)
c                                 d(ln(MDYN(J))/d(ln(q(J))=3.
              DPOSY(DKAL+1,I) = DPOSY(DKAL+1,I)+3.d0*DMFDM*JACY(J)
c                                 => d(POSX,POSY)/d(ln(q(J))
              DPOSX(DKAL+6,I) = DPOSX(DKAL+6,I)-2.d0*DMFDM*JACX(J)
c                                 d(ln(MDYN(J))/d(ln(per(J))=-2.
              DPOSY(DKAL+6,I) = DPOSY(DKAL+6,I)-2.d0*DMFDM*JACY(J)
c                                 => d(POSX,POSY)/d(ln(per(J))
              IF (J.GT.1) THEN
                 DPOSX(DKAL-6,I) = DPOSX(DKAL-6,I)-3.d0*DMFDM*JACX(J)
c                               d(ln(MDYN(J-1))/d(ln(q(J-1))=3.
                 DPOSY(DKAL-6,I) = DPOSY(DKAL-6,I)-3.d0*DMFDM*JACY(J)
c                               => d(POSX,POSY)/d(ln(q(J-1))
                 DPOSX(DKAL-1,I) = DPOSX(DKAL-1,I)+2.d0*DMFDM*JACX(J)
c                               d(ln(MDYN(J-1))/d(ln(per(J-1))=-2.
                 DPOSY(DKAL-1,I) = DPOSY(DKAL-1,I)+2.d0*DMFDM*JACY(J)
c                               => d(POSX,POSY)/d(ln(per(J-1))
              ELSE
                 DPOSX(N1,I) = DPOSX(N1,I)-DMFDM*JACX(J)
c                                     d(ln(MDYN(J-1))/d(ln(MSTAR)) = 1
                 DPOSY(N1,I) = DPOSY(N1,I)-DMFDM*JACY(J)
c                                     => d(POSX,POSY)/d(ln(MSTAR)
              END IF
           END DO
        END DO

        END    

C      
C -----------------------------------------------------------------------------
C     Planet/star radial velocity fit and derivatives (for Levenberg-Marquardt)
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQVFITP(N1,TT,P,VRADP,DVRADP)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    N1              ! Dimension of parameter space 
        REAL*8, DIMENSION(N1) :: P       ! Parameters 
        REAL*8, DIMENSION(NPLA) :: VRADP ! Radial velocities
        REAL*8, DIMENSION(N1,NPLA) :: DVRADP ! Derivatives / parameters
        REAL*8, DIMENSION(NPLA) :: MDYN,PER ! Dynamical masses & Periods 
        REAL*8, DIMENSION(NPLA) :: MFRAC ! Fractional masses
        REAL*8, DIMENSION(NPLA) :: VJAC ! Jacobi velocities
        REAL*8, DIMENSION(7,NPLA) :: DVJAC ! Jac. deriv.
        INTEGER*4 ::    I,J,DKAL
        REAL*8 ::
     &       TT,DT,             ! Temps
     &       MTOT,MSTAR,        ! Total mass up to planet #i
     &       Z,Q,               ! Periastron & ln()
     &       TP,                ! Tps de passage au periastre
     &       N,                 ! Mean motion / q
     &       K,H,               ! e*[cos(w),sin(w)]
     &       QQ,PP,             ! tan(i/2)*cos(phi),tan(i/2)*sin(phi)
     &       S,S0,              ! Universal variable (at T and T0)
     &       TI,T2,             ! tan(i/2) + square
     &       EXC,               ! eccentricity
     &       CW,SW,CP,SP,CO,SO, ! cos,sin (w,phi,O)
     &       CI,SI,COM,SOM,     ! cos,sin (i,Om)
     &       ALPHA,XX,UFA,      ! Energy
     &       VX,VY,             ! Velocity aling eigen axes
     &       U,CU,SU,           ! Anomalie excentrique
     &       FACT,              ! 1/(1-e*cos(u))
     &       RX,RY,             ! Coordonnes dans la base propre
     &       C0,C1,C2,C3,SC1,S2C2,S3C3,    ! Stumpff functions
     &       DTPDS0,DTPDALPHA,DTPDM,DTPDQ, ! Derivatives of TP
     &       DSDTP,DSDQ,DSDALPHA,DSDE,DSDM, ! Derivatives of S
     &       DOMDK,DOMDH,DODP,DODQ,DIDP,DIDQ,DOMDP,DOMDQ, ! Intern derivs. 
     &       DWDK,DWDH,DPHIDP,DPHIDQ, ! Derivees interm.
     &       DALPHADE,DALPHADM,DALPHADQ, !  
     &       DMFDM,DMFDMM,DUDH,DUDK, !
     &       DEDK,DEDH,         ! Derivees de Exc
     &       DVXDQ,DVXDS,DVXDE,DVXDS0,DVXDTP,DVXDM, ! Derivees de VX (internal)
     &       DVXDPER,DVXDALPHA,DVXDK,DVXDH,DVXDPP,DVXDQQ,
     &       DVYDQ,DVYDS,DVYDE,DVYDS0,DVYDTP,DVYDM, ! Derivees de VY (internal)
     &       DVYDPER,DVYDALPHA,DVYDK,DVYDH,DVYDPP,DVYDQQ,
     &       DVJACDQ,DVJACDZ,DVJACDK,DVJACDH,DVJACDQQ,DVJACDPER,
     &       DVJACDPP,DVJACDLP,DVJACDS0   ! Derivatives of VJAC   

        MSTAR = EXP(P(N1))
        MTOT = MSTAR       
        DO I = 1,NPLA
          DKAL = NEL*(I-1)
          Z = P(DKAL+1)
          K = P(DKAL+2)
          H = P(DKAL+3)
          QQ = P(DKAL+4)
          PP = P(DKAL+5)
          PER(I) = EXP(P(DKAL+6))
          S0 = P(DKAL+7)
          N = DPI/PER(I)         ! Mean motion / q
          Q = EXP(Z)            ! Semi-major axis
          MDYN(I) = N*N*(Q**3)    ! Dynamical mass
          IF (MULTIPLA) THEN
             MFRAC(I) = (MDYN(I)-MTOT)/MDYN(I) ! Mfrac(i) = (M(i)-M(i-1))/M(i)
             MTOT = MDYN(I)
          END IF
          EXC = SQRT(K*K+H*H) ! e
          CW = K/EXC                   ! cos(w)
          SW = H/EXC                   ! sin(w)
          ALPHA = MDYN(I)/Q*(1.d0-EXC) ! -2*Energy
          XX = S0*S0*ALPHA
          CALL KEPU_STUMPFF(XX,C0,C1,C2,C3)
          SC1 = C1*S0 
          S2C2 = C2*S0*S0
          S3C3 = C3*S0*S0*S0
          DT = MDYN(I)*S3C3+Q*SC1  ! Kepler's equation
          TP = STAR%T0-DT
          DTPDS0 = -MDYN(I)*S2C2-Q*C0
          DTPDALPHA = -(0.5d0/ALPHA)*(MDYN(I)*(S0*S2C2-3.d0*S3C3)
     &          +Q*(S0*C0-SC1))
          DTPDM = -S3C3
          DTPDQ = -SC1

          DT = TT-TP
          CALL KEP_UNIV(DT,Q,MDYN(I),ALPHA,S,C0,SC1,S2C2,S3C3)
          UFA = SQRT(Q*MDYN(I)*(1.d0+EXC))
          FACT = 1.d0/(MDYN(I)*S2C2+Q*C0) ! FACT = ds/dt = 1/r
          VX = -MDYN(I)*SC1*FACT
          VY = UFA*C0*FACT
          
          DSDTP = -FACT
          DSDQ = -SC1*FACT
          DSDM = -S3C3*FACT
          DSDALPHA = -0.5d0*FACT/ALPHA*
     &                  (MDYN(I)*(S*S2C2-3.d0*S3C3)+Q*(S*C0-SC1))
          DALPHADM = ALPHA/MDYN(I)
          DALPHADE = -MDYN(I)/Q
          DALPHADQ = -ALPHA/Q
          DSDE = DSDALPHA*DALPHADE
          DSDQ = DSDALPHA*DALPHADQ+DSDQ
          DSDM = DSDALPHA*DALPHADM+DSDM

          DVXDALPHA = 0.5d0*VX*((-SC1+S*C0)/SC1
     &          -MDYN(I)*(EXC*S*SC1-2.d0*S2C2)*FACT)/ALPHA
          DVXDS = VX*(C0/SC1-EXC*MDYN(I)*SC1*FACT)          
          DVXDQ = -C0*VX*FACT ! d/dFACT(=>q)
          DVXDM =  Q*C0*VX*FACT/MDYN(I)

          DVXDM = DVXDM+DVXDS*DSDM+DVXDALPHA*DALPHADM
          DVXDQ = DVXDQ+DVXDS*DSDQ+DVXDALPHA*DALPHADQ

          DVXDTP = DVXDS*DSDTP
          DVXDE = DVXDS*DSDE+DVXDALPHA*DALPHADE
          DVXDQ = DVXDQ+3.d0*DVXDM*MDYN(I)/Q
          DVXDPER = -2.d0*DVXDM*MDYN(I)/PER(I)

          DVXDQ = DVXDQ+DVXDTP*(DTPDQ+3.d0*DTPDM*MDYN(I)/Q)
          DVXDPER = DVXDPER-2.d0*DVXDTP*DTPDM*MDYN(I)/PER(I)
          DVXDQ = DVXDQ+2.d0*DVXDTP*DTPDALPHA*ALPHA/Q
          DVXDPER = DVXDPER-2.d0*DVXDTP*DTPDALPHA*ALPHA/PER(I)
          DVXDE = DVXDE-DVXDTP*DTPDALPHA*MDYN(I)/Q
          DVXDS0 = DVXDTP*DTPDS0          
          
          DVYDALPHA = -0.5d0*VY*(S*SC1/C0
     &          +MDYN(I)*(EXC*S*SC1-2.d0*S2C2)*FACT/ALPHA)          
          DVYDS = -VY*(ALPHA*SC1/C0+EXC*MDYN(I)*SC1*FACT)
          DVYDQ = 0.5d0*VY*(MDYN(I)*S2C2-Q*C0)*FACT/Q
          DVYDM = 0.5d0*VY*(-MDYN(I)*S2C2+Q*C0)*FACT/MTOT

          DVYDM = DVYDM+DVYDS*DSDM+DVYDALPHA*DALPHADM
          DVYDQ = DVYDQ+DVYDS*DSDQ+DVYDALPHA*DALPHADQ

          DVYDTP = DVYDS*DSDTP
          DVYDE = 0.5d0*VY/(1.d0+EXC)+DVYDS*DSDE+DVYDALPHA*DALPHADE
          DVYDQ = DVYDQ+3.d0*DVYDM*MDYN(I)/Q
          DVYDPER = -2.d0*DVYDM*MDYN(I)/PER(I)

          DVYDQ = DVYDQ+DVYDTP*(DTPDQ+3.d0*DTPDM*MDYN(I)/Q)
          DVYDPER = DVYDPER-2.d0*DVYDTP*DTPDM*MDYN(I)/PER(I)
          DVYDQ = DVYDQ+2.d0*DVYDTP*DTPDALPHA*ALPHA/Q
          DVYDPER = DVYDPER-2.d0*DVYDTP*DTPDALPHA*ALPHA/PER(I)
          DVYDE = DVYDE-DVYDTP*DTPDALPHA*MDYN(I)/Q
          DVYDS0 = DVYDTP*DTPDS0
          
          T2 = QQ*QQ+PP*PP        
          TI = SQRT(T2)        ! Tan(i/2)
          SI = 2.d0*TI/(1.d0+T2)      ! sin(i) 
          CI = (1.d0-T2)/(1.d0+T2)  ! cos(i)
          CO = QQ/TI         ! cos(Omega)
          SO = PP/TI         ! sin(Omega)
          COM = CW*CO+SW*SO  ! omega = w-Omega
          SOM = SW*CO-CW*SO  !
          VJAC(I) = (VX*SOM+VY*COM)*SI
          
          DEDK = CW
          DEDH = SW
          DOMDK = -SW/EXC       ! = dw/dk (OM = W-O)
          DOMDH = CW/EXC  ! = dw/dh
          DODQ = -PP/T2
          DODP = QQ/T2
          DOMDQ = PP/T2         ! -dO/dq
          DOMDP = -QQ/T2        ! -dO/dp
          DIDP = 2.d0*SO/(1.d0+T2)
          DIDQ = 2.d0*CO/(1.d0+T2)
          
          DVXDK = DVXDE*DEDK
          DVXDH = DVXDE*DEDH
          DVYDK = DVYDE*DEDK
          DVYDH = DVYDE*DEDH

          DVJACDQ = (DVXDQ*SOM+DVYDQ*COM)*SI
          DVJACDK = (DVXDK*SOM+DVYDK*COM+(VX*COM-VY*SOM)*DOMDK)*SI
          DVJACDH = (DVXDH*SOM+DVYDH*COM+(VX*COM-VY*SOM)*DOMDH)*SI          
          DVJACDPP = (VX*COM-VY*SOM)*DOMDP*SI+(VX*SOM+VY*COM)*CI*DIDP
          DVJACDQQ = (VX*COM-VY*SOM)*DOMDQ*SI+(VX*SOM+VY*COM)*CI*DIDQ
          DVJACDPER = (DVXDPER*SOM+DVYDPER*COM)*SI
          DVJACDS0 = (DVXDS0*SOM+DVYDS0*COM)*SI

          DVJACDZ = DVJACDQ*Q
          DVJACDLP = DVJACDPER*PER(I)

          DVJAC(:,I) =
     & (/ DVJACDZ,DVJACDK,DVJACDH,DVJACDQQ,DVJACDPP,DVJACDLP,DVJACDS0 /) 
          
        END DO
        
        DVRADP = 0.d0
c...  Radial velocities of planets relative to star
        DO I = 1,NPLA
c... Basic formula Jacobi coords => helio vels. (*(-1) for radial vels.)
           VRADP(I) = -VJAC(I)-SUM(MFRAC(1:I-1)*VJAC(1:I-1))
c...  Derivatives
           DKAL = NEL*(I-1)
           DVRADP(DKAL+1:DKAL+7,I) = -DVJAC(:,I)
           DO J = 1,I-1
              DKAL = NEL*(J-1)   ! Derivatives relative to VJAC(J)
              DVRADP(DKAL+1:DKAL+7,I) = DVRADP(DKAL+1:DKAL+7,I)
     &                -MFRAC(J)*DVJAC(:,J)
              IF (J.GT.1) THEN   ! Derivatives relative to MFRAC(J)
                 MTOT = MDYN(J-1)  ! MFRAC(J) = (MDYN(J)-MTOT)/MDYN(J)
              ELSE
                 MTOT = MSTAR
              END IF
              DMFDM = MTOT/MDYN(J) ! d(MFRAC(J))/d(ln(MDYN(J))
              DMFDMM = DMFDM*PER(J)*PER(J)
              DVRADP(DKAL+1,I) = DVRADP(DKAL+1,I)-3.d0*DMFDM*VJAC(J)
c                        d(ln(MDYN(J))/d(ln(a(J))=3.=> d(VRADP)/d(ln(a(J))
              DVRADP(DKAL+6,I) = DVRADP(DKAL+6,I)
     &                             -2.d0*DMFDMM*VJAC(J)*P(DKAL+6)
c              d(ln(MDYN(J))/d(P6(J))=2*P6/(P6^2+P7^2)=2*PER(J)^2*P6(J)
c                                                  => d(VRADP)/d(P6(J))
              DVRADP(DKAL+7,I) = DVRADP(DKAL+7,I)
     &                             -2.d0*DMFDMM*VJAC(J)*P(DKAL+7)
c              d(ln(MDYN(J))/d(P7(J))=2*P7/(P6^2+P7^2)=2*PER(J)^2*P7(J)
c                                                  => d(POSX,POSY)/d(P7(J)         
              IF (J.GT.1) THEN
                 DMFDMM = DMFDM*PER(J-1)*PER(J-1)
                 DVRADP(DKAL-6,I) = DVRADP(DKAL-6,I)+3.d0*DMFDM*VJAC(J)
c                  d(ln(MDYN(J-1))/d(ln(a(J-1))=3.=> d(VRADP)/d(ln(a(J-1))
                 DVRADP(DKAL-1,I) = DVRADP(DKAL-1,I)
     &                      +2.d0*DMFDMM*VJAC(J)*P(DKAL-1)
c            d(ln(MDYN(J-1))/d(P6(J-1))=2*P6/(P6^2+P7^2)=2*PER(J-1)^2*P6(J-1)
c                                              => d(VRADP)/d(P6(J-1))
                 DVRADP(DKAL,I) = DVRADP(DKAL,I)
     &                      +2.d0*DMFDMM*VJAC(J)*P(DKAL)
c            d(ln(MDYN(J-1))/d(P7(J-1))=2*P7/(P6^2+P7^2)=2*PER(J-1)^2*P7(J-1)
c                                              => d(VRADP)/d(P7(J-1))
              ELSE
                 DVRADP(N1,I) = DVRADP(N1,I)+DMFDM*VJAC(J)
c                d(ln(MDYN(J-1))/d(ln(MSTAR)) = 1 => d(VRADP)/d(ln(MSTAR)
              END IF
           END DO
        END DO

        END    
C
C -----------------------------------------------------------------------------
C       Radial velocity fit and derivatives (For Leveneberg-Marquardt)
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQVFIT(N1,TT,P,VRAD,DVRAD)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    N1              ! Dimension of parameter space 
        REAL*8, DIMENSION(NPLA) :: MDYN,PER ! Dynamical masses & Periods 
        REAL*8, DIMENSION(NPLA) :: MFRAC ! Fractional masses
        REAL*8, DIMENSION(NPLA) :: VJAC ! Jacobi velocities
        REAL*8, DIMENSION(7,NPLA) :: DVJAC ! Jac. deriv.
        REAL*8, DIMENSION(N1) :: P,DVRAD ! Parameters & derivatives
        INTEGER*4 ::    I,J,DKAL
        REAL*8 ::
     &       TT,DT,             ! Temps
     &       VRAD,              ! Radial velocity
     &       MTOT,MSTAR,        ! Total mass up to planet #i
     &       Z,Q,               ! Periastron & ln()
     &       TP,                ! Tps de passage au periastre
     &       N,                 ! Mean motion / q
     &       K,H,               ! e*[cos(w),sin(w)]
     &       QQ,PP,             ! tan(i/2)*cos(phi),tan(i/2)*sin(phi)
     &       S,S0,              ! Universal variable (at T and T0)
     &       TI,T2,             ! tan(i/2) + square
     &       EXC,               ! eccentricity
     &       CW,SW,CP,SP,CO,SO, ! cos,sin (w,phi,O)
     &       CI,SI,COM,SOM,     ! cos,sin (i,Om)
     &       ALPHA,XX,UFA,      ! Energy
     &       VX,VY,             ! Velocity aling eigen axes
     &       U,CU,SU,           ! Anomalie excentrique
     &       FACT,              ! 1/(1-e*cos(u))
     &       RX,RY,             ! Coordonnes dans la base propre
     &       C0,C1,C2,C3,SC1,S2C2,S3C3,    ! Stumpff functions
     &       DTPDS0,DTPDALPHA,DTPDM,DTPDQ, ! Derivatives of TP
     &       DSDTP,DSDQ,DSDALPHA,DSDE,DSDM, ! Derivatives of S
     &       DOMDK,DOMDH,DODP,DODQ,DIDP,DIDQ,DOMDP,DOMDQ, ! Intern derivs. 
     &       DWDK,DWDH,DPHIDP,DPHIDQ, ! Derivees interm.
     &       DALPHADE,DALPHADM,DALPHADQ, !  
     &       DMFDM,DMFDMM,DUDH,DUDK, !
     &       DEDK,DEDH,         ! Derivees de Exc
     &       DVXDQ,DVXDS,DVXDE,DVXDS0,DVXDTP,DVXDM, ! Derivees de VX (internal)
     &       DVXDPER,DVXDALPHA,DVXDK,DVXDH,DVXDPP,DVXDQQ,
     &       DVYDQ,DVYDS,DVYDE,DVYDS0,DVYDTP,DVYDM, ! Derivees de VY (internal)
     &       DVYDPER,DVYDALPHA,DVYDK,DVYDH,DVYDPP,DVYDQQ,
     &       DVJACDQ,DVJACDZ,DVJACDK,DVJACDH,DVJACDQQ,DVJACDPER,
     &       DVJACDPP,DVJACDLP,DVJACDS0   ! Derivatives of VJAC   

        VRAD = P(NEL*NPLA+1)    ! Offset velocity
        DVRAD = 0.d0
        MSTAR = EXP(P(N1))
        MTOT = MSTAR       
        DO I = 1,NPLA
          DKAL = NEL*(I-1)
          Z = P(DKAL+1)
          K = P(DKAL+2)
          H = P(DKAL+3)
          QQ = P(DKAL+4)
          PP = P(DKAL+5)
          PER(I) = EXP(P(DKAL+6))
          S0 = P(DKAL+7)
          N = DPI/PER(I)         ! Mean motion / q
          Q = EXP(Z)            ! Semi-major axis
          MDYN(I) = N*N*(Q**3)    ! Dynamical mass
          IF (MULTIPLA) THEN
             MFRAC(I) = (MDYN(I)-MTOT)/MDYN(I) ! Mfrac(i) = (M(i)-M(i-1))/M(i)
             MTOT = MDYN(I)
          END IF
          EXC = SQRT(K*K+H*H) ! e
          CW = K/EXC                   ! cos(w)
          SW = H/EXC                   ! sin(w)
          ALPHA = MDYN(I)/Q*(1.d0-EXC) ! -2*Energy
          XX = S0*S0*ALPHA
          CALL KEPU_STUMPFF(XX,C0,C1,C2,C3)
          SC1 = C1*S0 
          S2C2 = C2*S0*S0
          S3C3 = C3*S0*S0*S0
          DT = MDYN(I)*S3C3+Q*SC1  ! Kepler's equation
          TP = STAR%T0-DT
          DTPDS0 = -PLA(I)%MDYN*S2C2-PLA(I)%Q*C0
          DTPDALPHA = -(0.5d0/ALPHA)*(MDYN(I)*(S0*S2C2-3.d0*S3C3)
     &          +Q*(S0*C0-SC1))
          DTPDM = -S3C3
          DTPDQ = -SC1

          DT = TT-TP
          CALL KEP_UNIV(DT,Q,MDYN(I),ALPHA,S,C0,SC1,S2C2,S3C3)
          UFA = SQRT(Q*MDYN(I)*(1.d0+EXC))
          FACT = 1.d0/(MDYN(I)*S2C2+Q*C0) ! FACT = ds/dt = 1/r
          VX = -MDYN(I)*SC1*FACT
          VY = UFA*C0*FACT
          
          DSDTP = -FACT
          DSDQ = -SC1*FACT
          DSDM = -S3C3*FACT
          DSDALPHA = -0.5d0*FACT/ALPHA*
     &                  (MDYN(I)*(S*S2C2-3.d0*S3C3)+Q*(S*C0-SC1))
          DALPHADM = ALPHA/MDYN(I)
          DALPHADE = -MDYN(I)/Q
          DALPHADQ = -ALPHA/Q
          DSDE = DSDALPHA*DALPHADE
          DSDQ = DSDALPHA*DALPHADQ+DSDQ
          DSDM = DSDALPHA*DALPHADM+DSDM

          DVXDALPHA = 0.5d0*VX*((-SC1+S*C0)/SC1
     &          -MDYN(I)*(EXC*S*SC1-2.d0*S2C2)*FACT)/ALPHA
          DVXDS = VX*(C0/SC1-EXC*MDYN(I)*SC1*FACT)          
          DVXDQ = -C0*VX*FACT ! d/dFACT(=>q)
          DVXDM =  Q*C0*VX*FACT/MDYN(I)

          DVXDM = DVXDM+DVXDS*DSDM+DVXDALPHA*DALPHADM
          DVXDQ = DVXDQ+DVXDS*DSDQ+DVXDALPHA*DALPHADQ

          DVXDTP = DVXDS*DSDTP
          DVXDE = DVXDS*DSDE+DVXDALPHA*DALPHADE
          DVXDQ = DVXDQ+3.d0*DVXDM*MDYN(I)/Q
          DVXDPER = -2.d0*DVXDM*MDYN(I)/PER(I)

          DVXDQ = DVXDQ+DVXDTP*(DTPDQ+3.d0*DTPDM*MDYN(I)/Q)
          DVXDPER = DVXDPER-2.d0*DVXDTP*DTPDM*MDYN(I)/PER(I)
          DVXDQ = DVXDQ+2.d0*DVXDTP*DTPDALPHA*ALPHA/Q
          DVXDPER = DVXDPER-2.d0*DVXDTP*DTPDALPHA*ALPHA/PER(I)
          DVXDE = DVXDE-DVXDTP*DTPDALPHA*MDYN(I)/Q
          DVXDS0 = DVXDTP*DTPDS0
          
          DVYDALPHA = -0.5d0*VY*(S*SC1/C0
     &          +MDYN(I)*(EXC*S*SC1-2.d0*S2C2)*FACT/ALPHA)          
          DVYDS = -VY*(ALPHA*SC1/C0+EXC*MDYN(I)*SC1*FACT)
          DVYDQ = 0.5d0*VY*(MDYN(I)*S2C2-Q*C0)*FACT/Q
          DVYDM = 0.5d0*VY*(-MDYN(I)*S2C2+Q*C0)*FACT/MTOT

          DVYDM = DVYDM+DVYDS*DSDM+DVYDALPHA*DALPHADM
          DVYDQ = DVYDQ+DVYDS*DSDQ+DVYDALPHA*DALPHADQ

          DVYDTP = DVYDS*DSDTP
          DVYDE = 0.5d0*VY/(1.d0+EXC)+DVYDS*DSDE+DVYDALPHA*DALPHADE
          DVYDQ = DVYDQ+3.d0*DVYDM*MDYN(I)/Q
          DVYDPER = -2.d0*DVYDM*MDYN(I)/PER(I)

          DVYDQ = DVYDQ+DVYDTP*(DTPDQ+3.d0*DTPDM*MDYN(I)/Q)
          DVYDPER = DVYDPER-2.d0*DVYDTP*DTPDM*MDYN(I)/PER(I)
          DVYDQ = DVYDQ+2.d0*DVYDTP*DTPDALPHA*ALPHA/Q
          DVYDPER = DVYDPER-2.d0*DVYDTP*DTPDALPHA*ALPHA/PER(I)
          DVYDE = DVYDE-DVYDTP*DTPDALPHA*MDYN(I)/Q
          DVYDS0 = DVYDTP*DTPDS0

          T2 = QQ*QQ+PP*PP        
          TI = SQRT(T2)        ! Tan(i/2)
          SI = 2.d0*TI/(1.d0+T2)      ! sin(i) 
          CI = (1.d0-T2)/(1.d0+T2)  ! cos(i)
          CO = QQ/TI         ! cos(Omega)
          SO = PP/TI         ! sin(Omega)
          COM = CW*CO+SW*SO  ! omega = w-Omega
          SOM = SW*CO-CW*SO  !
          VJAC(I) = (VX*SOM+VY*COM)*SI
          
          DEDK = CW
          DEDH = SW
          DOMDK = -SW/EXC       ! = dw/dk (OM = W-O)
          DOMDH = CW/EXC  ! = dw/dh
          DODQ = -PP/T2
          DODP = QQ/T2
          DOMDQ = PP/T2         ! -dO/dq
          DOMDP = -QQ/T2        ! -dO/dp
          DIDP = 2.d0*SO/(1.d0+T2)
          DIDQ = 2.d0*CO/(1.d0+T2)
          
          DVXDK = DVXDE*DEDK
          DVXDH = DVXDE*DEDH
          DVYDK = DVYDE*DEDK
          DVYDH = DVYDE*DEDH

          DVJACDQ = (DVXDQ*SOM+DVYDQ*COM)*SI
          DVJACDK = (DVXDK*SOM+DVYDK*COM+(VX*COM-VY*SOM)*DOMDK)*SI
          DVJACDH = (DVXDH*SOM+DVYDH*COM+(VX*COM-VY*SOM)*DOMDH)*SI          
          DVJACDPP = (VX*COM-VY*SOM)*DOMDP*SI+(VX*SOM+VY*COM)*CI*DIDP
          DVJACDQQ = (VX*COM-VY*SOM)*DOMDQ*SI+(VX*SOM+VY*COM)*CI*DIDQ
          DVJACDPER = (DVXDPER*SOM+DVYDPER*COM)*SI
          DVJACDS0 = (DVXDS0*SOM+DVYDS0*COM)*SI

          DVJACDZ = DVJACDQ*Q
          DVJACDLP = DVJACDPER*PER(I)

          DVJAC(:,I) =
     & (/ DVJACDZ,DVJACDK,DVJACDH,DVJACDQQ,DVJACDPP,DVJACDLP,DVJACDS0 /) 
          
        END DO

c... Basic formula Jacobi coords => radial velocity
        VRAD = VRAD+SUM(MFRAC(1:NPLA)*VJAC(1:NPLA))
        DVRAD = 0.d0
        DVRAD(NEL*NPLA+1) = 1.d0
        DO J = 1,NPLA
           DKAL = NEL*(J-1)
           DVRAD(DKAL+1:DKAL+7) = MFRAC(J)*DVJAC(:,J)
           IF (J.GT.1) THEN   ! Derivatives relative to MFRAC(J)
              MTOT = MDYN(J-1)  ! MFRAC(J) = (MDYN(J)-MTOT)/MDYN(J)
           ELSE
              MTOT = MSTAR
           END IF
           DMFDM = MTOT/MDYN(J) ! d(MFRAC(J))/d(ln(MDYN(J))
           DMFDMM = DMFDM*PER(J)*PER(J)
           DVRAD(DKAL+1) = DVRAD(DKAL+1)+3.d0*DMFDM*VJAC(J)
c                d(ln(MDYN(J))/d(ln(a(J))=3. => d(VRAD)/d(ln(a(J))
           DVRAD(DKAL+6) = DVRAD(DKAL+6)+2.d0*DMFDMM*VJAC(J)*P(DKAL+6)
c              d(ln(MDYN(J))/d(P6(J))=2*P6/(P6^2+P7^2)=2*PER(J)^2*P6(J)
c                                           => d(VRAD)/d(P6(J))
           DVRAD(DKAL+7) = DVRAD(DKAL+7)+2.d0*DMFDMM*VJAC(J)*P(DKAL+7)
c              d(ln(MDYN(J))/d(P7(J))=2*P7/(P6^2+P7^2)=2*PER(J)^2*P7(J)
c                                           => d(VRAD)/d(P7(J))   
           IF (J.GT.1) THEN
              DMFDMM = DMFDM*PER(J-1)*PER(J-1)
              DVRAD(DKAL-6) = DVRAD(DKAL-6)-3.d0*DMFDM*VJAC(J)
c           d(ln(MDYN(J-1))/d(ln(a(J-1))=3. => d(POSX,POSY)/d(ln(a(J-1))
              DVRAD(DKAL-1) = DVRAD(DKAL-1)
     &                      -2.d0*DMFDMM*VJAC(J)*P(DKAL-1)
c           d(ln(MDYN(J-1))/d(P6(J-1))=2*P6/(P6^2+P7^2)=2*PER(J-1)^2*P6(J-1)
c                                           => d(VRAD)/d(P6(J-1))
              DVRAD(DKAL) = DVRAD(DKAL)
     &                      -2.d0*DMFDMM*VJAC(J)*P(DKAL)
c           d(ln(MDYN(J-1))/d(P7(J-1))=2*P7/(P6^2+P7^2)=2*PER(J-1)^2*P7(J-1)
c                                           => d(VRAD)/d(P7(J-1))              
           ELSE
              DVRAD(N1) = DVRAD(N1)-DMFDM*VJAC(J)
c                                     d(ln(MDYN(J-1))/d(ln(MSTAR)) = 1
           END IF
        END DO        
        END    

      
C
C -----------------------------------------------------------------------------
C       Computation of orbital elements and uncertainties from LM solution
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQELEMENTS(P,COV)

        USE DATA

        IMPLICIT NONE

        REAL*8, DIMENSION(NPAR) ::
     &                  DEDP,DWDP,DIDP,DMUDP,DPERDP,DTPDP
        INTEGER*4 ::    I,J,DKAL
        REAL*8 ::
     &       P(NPAR),                   ! Parameters
     &       COV(NPAR,NPAR),            ! Covariances
     &       NN,DNN,                    ! Moyen mouvement
     &       K,H,                       ! e*(cos(w),sin(w))        
     &       S,                         ! Universal variable at T0
     &       ALPHA,XX,T2,               ! Auxiliary variables
     &       C0,C1,C2,C3,SC1,S2C2,S3C3, ! Stumpff functions
     &       TT,                        ! T-T0 (Kepler's equation)   
     &       DTTDS,DTTDA,DTTDM,DTTDQ,   ! Partial derivatives of TT
     &       QQ,PP,                     ! tan(i/2)*(cos(phi),sin(phi))
     &       ERR_COMB                   ! Fonction erreur combinée
        
        STAR%SIGJV = 0.d0
        IF (MULTIPLA) THEN
           STAR%MASS = EXP(P(NPAR))
           STAR%DMASS =  SQRT(ABS(COV(NPLA,NPLA)))*STAR%MASS
        END IF
        IF (RADVEL) THEN
           STAR%V0 = P(NEL*NPLA+1)
           STAR%DV0 = SQRT(ABS(COV(NEL*NPLA+1,NEL*NPLA+1)))
           IF (JITNUM.EQ.1) THEN
              STAR%SIGJV = EXP(P(NEL*NPLA+2))
              STAR%DSIGJV = SQRT(ABS(COV(NEL*NPLA+2,NEL*NPLA+2)))
     &                                                *STAR%SIGJV
           END IF
        END IF
        DO I = 1,NPLA
           DKAL = NEL*(I-1)
           PLA(I)%Q = EXP(P(DKAL+1))
           PLA(I)%DQ = SQRT(ABS(COV(DKAL+1,DKAL+1)))*PLA(I)%Q
           PLA(I)%PER = EXP(P(DKAL+6))
           PLA(I)%DPER = SQRT(ABS(COV(DKAL+6,DKAL+6)))*PLA(I)%PER
           NN = DPI/PLA(I)%PER
           S = P(DKAL+7)
           K = P(DKAL+2)
           H = P(DKAL+3)
           QQ = P(DKAL+4)
           PP = P(DKAL+5)
           PLA(I)%EXC = SQRT(K*K+H*H)
           PLA(I)%A = PLA(I)%Q/(1.d0-PLA(I)%EXC)
           PLA(I)%CW = K/PLA(I)%EXC ! cos(w=omega+Omega)
           PLA(I)%SW = H/PLA(I)%EXC          ! sin(w)
           PLA(I)%W = ATAN2(PLA(I)%SW,PLA(I)%CW)*DR
           DEDP(1:NPAR) = 0.d0
           DWDP(1:NPAR) = 0.d0
           DWDP(DKAL+2) = -PLA(I)%SW/PLA(I)%EXC
           DWDP(DKAL+3) = PLA(I)%CW/PLA(I)%EXC
           DEDP(DKAL+2) = PLA(I)%CW
           DEDP(DKAL+3) = PLA(I)%SW
           PLA(I)%DEXC = ERR_COMB(NPAR,DEDP(1:NPAR),COV(1:NPAR,1:NPAR))         
           PLA(I)%DW = ERR_COMB(NPAR,DWDP(1:NPAR),COV(1:NPAR,1:NPAR))*DR
           PLA(I)%MDYN = NN*NN*PLA(I)%Q**3
           DMUDP = 0.d0
           DMUDP(DKAL+1) = 3.d0*PLA(I)%MDYN   ! dm/m = 3*dq/q-2*dP/P
           DMUDP(DKAL+6) = -2.d0*PLA(I)%MDYN
           PLA(I)%DMDYN =
     &              ERR_COMB(NPAR,DMUDP(1:NPAR),COV(1:NPAR,1:NPAR))           
           IF (I.EQ.1) THEN
             IF (NPLA.GT.1) THEN
                PLA(I)%MU = PLA(I)%MDYN-STAR%MASS
                DMUDP(NPAR) = -STAR%MASS   ! dmu = dmass = -d(lnmass)*mass 
                PLA(I)%DMU =
     &                ERR_COMB(NPAR,DMUDP(1:NPAR),COV(1:NPAR,1:NPAR))
             ELSE
                STAR%MASS = PLA(I)%MDYN    ! Only one orbit
                STAR%DMASS =
     &                ERR_COMB(NPAR,DMUDP(1:NPAR),COV(1:NPAR,1:NPAR))
             END IF
           ELSE
             PLA(I)%MU = PLA(I)%MDYN-PLA(I-1)%MDYN
             DMUDP(DKAL+1-NEL) = -3.d0*PLA(I-1)%MDYN 
             DMUDP(DKAL+6-NEL) = +2.d0*PLA(I-1)%MDYN
             PLA(I)%DMU =
     &                ERR_COMB(NPAR,DMUDP(1:NPAR),COV(1:NPAR,1:NPAR))
           END IF
           DMUDP(NPAR) = 0.d0
           DMUDP(DKAL+1-NEL) = 0.d0
           DMUDP(DKAL+6-NEL) = 0.d0           

           ALPHA = PLA(I)%MDYN/PLA(I)%Q*(1.d0-PLA(I)%EXC)
           XX = S*S*ALPHA
           CALL KEPU_STUMPFF(XX,C0,C1,C2,C3)
           SC1 = C1*S
           S2C2 = C2*S*S
           S3C3 = C3*S*S*S
           TT = PLA(I)%MDYN*S3C3+PLA(I)%Q*SC1 ! Kepler's equation
           DTTDS = PLA(I)%MDYN*S2C2+PLA(I)%Q*C0
           DTTDA = (0.5d0/ALPHA)*(PLA(I)%MDYN*(S*S2C2-3.d0*S3C3)
     &          +PLA(I)%Q*(S*C0-SC1))
           DTTDM = S3C3
           DTTDQ = SC1
           PLA(I)%TP = STAR%T0-TT         
           DTPDP(1:NPAR) = 0.d0
           DTPDP(DKAL+1) = -DTTDQ*PLA(I)%Q
           DTPDP(1:NPAR) = DTPDP(1:NPAR)-DTTDM*DMUDP(1:NPAR) 
           DTPDP(DKAL+7) = DTPDP(DKAL+7)-DTTDS
           DTPDP(DKAL+1) = DTPDP(DKAL+1)+DTTDA*ALPHA
           DTPDP(1:NPAR) = DTPDP(1:NPAR)
     &          -DTTDA*ALPHA/PLA(I)%MDYN*DMUDP(1:NPAR)
     &          +DTTDA*PLA(I)%MDYN/PLA(I)%Q*DEDP(1:NPAR)
           PLA(I)%DTP =
     &             ERR_COMB(NPAR,DTPDP(1:NPAR),COV(1:NPAR,1:NPAR))

           T2 = QQ*QQ+PP*PP        
           PLA(I)%TI = SQRT(T2)        ! Tan(i/2)
           PLA(I)%INC = 2.d0*ATAN(PLA(I)%TI)*DR
           DWDP(1:NPAR) = 0.d0
           DIDP(1:NPAR) = 0.d0
           DWDP(DKAL+4) = -PP/T2
           DWDP(DKAL+5) = QQ/T2

           IF (RADVEL) THEN
              PLA(I)%CO = QQ/PLA(I)%TI ! cos(Omega)
              PLA(I)%SO = PP/PLA(I)%TI         ! sin(Omega)
              PLA(I)%O = ATAN2(PLA(I)%SO,PLA(I)%CO)*DR
              DIDP(4) =  2.d0*PLA(I)%CO/(1.d0+T2)
              DIDP(5) =  2.d0*PLA(I)%SO/(1.d0+T2)
              PLA(I)%DOO = ERR_COMB(NPAR,DWDP(1:NPAR),
     &                               COV(1:NPAR,1:NPAR))*DR
              PLA(I)%OM = PLA(I)%W-PLA(I)%O
              PLA(I)%PHI = PLA(I)%OM-PLA(I)%O
              PLA(I)%CP = COS(PLA(I)%PHI/DR)
              PLA(I)%SP = SIN(PLA(I)%PHI/DR)
              PLA(I)%DINC = ERR_COMB(NPAR,DIDP(1:NPAR),
     &                               COV(1:NPAR,1:NPAR))*DR
              PLA(I)%COM = COS(PLA(I)%OM/DR)
              PLA(I)%SOM = SIN(PLA(I)%OM/DR)

              DWDP(1:NPAR) = 0.d0
              DWDP(DKAL+2) = -PLA(I)%SW/PLA(I)%EXC
              DWDP(DKAL+3) = PLA(I)%CW/PLA(I)%EXC
              DWDP(DKAL+4) = PP/T2
              DWDP(DKAL+5) = -QQ/T2
              PLA(I)%DOM = ERR_COMB(NPAR,
     &                        DWDP(1:NPAR),COV(1:NPAR,1:NPAR))*DR
              DWDP(1:NPAR) = 0.d0
              DWDP(DKAL+2) = -PLA(I)%SW/PLA(I)%EXC
              DWDP(DKAL+3) = PLA(I)%CW/PLA(I)%EXC           
              DWDP(DKAL+4) = 2.d0*PP/T2
              DWDP(DKAL+5) = -2.d0*QQ/T2
              PLA(I)%DPHI = ERR_COMB(NPAR,
     &                        DWDP(1:NPAR),COV(1:NPAR,1:NPAR))*DR
           ELSE
              PLA(I)%CP = QQ/PLA(I)%TI         ! cos(omega-Omega)
              PLA(I)%SP = PP/PLA(I)%TI ! sin(omega-Omega)
              PLA(I)%PHI = ATAN2(PLA(I)%SP,PLA(I)%CP)*DR  
              DIDP(4) =  2.d0*PLA(I)%CP/(1.d0+T2)
              DIDP(5) =  2.d0*PLA(I)%SP/(1.d0+T2)
              PLA(I)%DPHI = ERR_COMB(NPAR,DWDP(1:NPAR),
     &                               COV(1:NPAR,1:NPAR))*DR
              PLA(I)%OM = 0.5d0*(PLA(I)%W+PLA(I)%PHI)
              PLA(I)%O = 0.5d0*(PLA(I)%W-PLA(I)%PHI)
              PLA(I)%CO = COS(PLA(I)%O/DR)
              PLA(I)%SO = SIN(PLA(I)%O/DR)
              PLA(I)%COM = COS(PLA(I)%OM/DR)
              PLA(I)%SOM = SIN(PLA(I)%OM/DR)
              PLA(I)%CO = COS(PLA(I)%O/DR)
              PLA(I)%SO = SIN(PLA(I)%O/DR)

              DWDP(1:NPAR) = 0.d0
              DWDP(DKAL+2) = -0.5d0*PLA(I)%SW/PLA(I)%EXC
              DWDP(DKAL+3) = 0.5d0*PLA(I)%CW/PLA(I)%EXC           
              DWDP(DKAL+4) = -0.5d0*PP/T2
              DWDP(DKAL+5) = 0.5d0*QQ/T2
              PLA(I)%DOM = ERR_COMB(NPAR,
     &                        DWDP(1:NPAR),COV(1:NPAR,1:NPAR))*DR
              DWDP(1:NPAR) = 0.d0
              DWDP(DKAL+2) = -0.5d0*PLA(I)%SW/PLA(I)%EXC
              DWDP(DKAL+3) = 0.5d0*PLA(I)%CW/PLA(I)%EXC     
              DWDP(DKAL+4) = 0.5d0*PP/T2
              DWDP(DKAL+5) = -0.5d0*QQ/T2
              PLA(I)%DOO = ERR_COMB(NPAR,
     &                        DWDP(1:NPAR),COV(1:NPAR,1:NPAR))*DR
           END IF
        END DO

        END      
      
C
C -----------------------------------------------------------------------------
C       Checking validity of current model in Levenberg-Marquardt
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQCHECK(N,P,TEST)

        USE DATA

        IMPLICIT NONE

        INTEGER*4       N               ! Dimension of parameter space

        REAL*8, DIMENSION(NPLA) :: MDYN ! Planet masses 
        REAL*8, DIMENSION(N) :: P    ! Parameters
        INTEGER*4 ::    I,J,DKAL
        REAL*8 ::       MTOT,MSTAR,MFRAC,! Total mass (partial)
     &                  Z,Q,            ! Periaston & ln()
     &                  TP,             ! Tps de passage au periastre
     &                  NN,              ! Moyen mouvement / q
     &                  K,H,            ! e*[cos(w),sin(w)]/sqrt(1-e^2)
     &                  QQ,PP,          ! tan(i/2)*[cos(O),sin(O)]
     &                  PER,            ! Period / q
     &                  S0,             ! e/sqrt(1-e^2) + carre
     &                  EXC             ! excentricite

        LOGICAL ::      TEST
        
        TEST = .FALSE.
c       print*,'avant',sngl(p)
        IF (MULTIPLA) THEN
           MSTAR = EXP(P(N))
        ELSE
           MSTAR = STAR%MASS
        END IF
        MTOT = MSTAR
        DO I = 1,NPLA
          DKAL = NEL*(I-1)
          Z = P(DKAL+1)
          K = P(DKAL+2)
          H = P(DKAL+3)
          QQ = P(DKAL+4)
          PP = P(DKAL+5)
          Q = EXP(Z)            ! Periastron       
          PER = EXP(P(DKAL+6))
          S0 = P(DKAL+7)
          NN = DPI/PER          ! Mean motion
          MDYN(I) = NN*NN*(Q**3)    ! Dynamical mass
          IF (MULTIPLA) THEN
             MFRAC = (MDYN(I)-MTOT)/MDYN(I) ! Mfrac = (M(i)-M(i-1))/M(i)
             MTOT = MDYN(I)
c             TEST = TEST.OR.(MFRAC.LE.0.d0)
c             IF (TEST) print*,'Bahhh !'
          END IF
          EXC = SQRT(K*K+H*H)
          TEST = TEST.OR.(EXC.GT.VPRIOR(3)%BOUND(2))
c          if (test) print*,'test exc',test,exc,k,h
c          TEST = TEST.OR.(PER.GT.VPRIOR(1)%BOUND(2))
c     &               .OR.(PER.LT.VPRIOR(1)%BOUND(1))
c          if (test) print*,'test per',test,per
          TEST = TEST.OR.(Q.GT.VPRIOR(2)%BOUND(2))
     &               .OR.(Q.LT.VPRIOR(2)%BOUND(1))
c          if (test) print*,'test a',test,a
        END DO
        IF (MULTIPLA)
     &          TEST = TEST.OR.(MSTAR.LT.MTINY).OR.(MSTAR.GT.MHUGE)
c        if (test) print*,'masses',mtiny,mstar,mhuge,a,per,exc
        DO I = 1,N
          IF (P(I).NE.P(I)) TEST = .TRUE.
        END DO
c        if (test) stop
        END    
C
      
C -----------------------------------------------------------------------------
C       Restart from new solution in L-M
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQRENEW(N,P)

        USE DATA

        IMPLICIT NONE

        INTEGER*4       N               ! Dimension of parameter space

        REAL*8, DIMENSION(NPAR) :: P    ! Parameters
        INTEGER*4 ::    I,DKAL
        REAL*8 ::       R,      ! Random number
     &                  MASS,   ! Dynamical mass   
     &                  M0,     ! Mean anomaly at T0
     &                  HH2,    ! Energy
     &                  C0,SC1,S2C2,S3C3, ! Stumpff functions        
     &                  S       ! Universal variable

        DO I = 1,NPLA
           CALL RANDOM_NUMBER(R)
           PLA(I)%O = 360.d0*R
           CALL RANDOM_NUMBER(R)
           PLA(I)%OM = 360.d0*R
           PLA(I)%W = PLA(I)%OM+PLA(I)%O
           PLA(I)%PHI = PLA(I)%OM-PLA(I)%O
           PLA(I)%CW = COS(PLA(I)%W/DR)
           PLA(I)%SW = SIN(PLA(I)%W/DR)
           PLA(I)%CP = COS(PLA(I)%PHI/DR)
           PLA(I)%SP = SIN(PLA(I)%PHI/DR)
           PLA(I)%COM = COS(PLA(I)%OM/DR)
           PLA(I)%SOM = SIN(PLA(I)%OM/DR)              
           PLA(I)%CO = COS(PLA(I)%O/DR)
           PLA(I)%SO = SIN(PLA(I)%O/DR)   
           DKAL = NEL*(I-1)
           P(DKAL+1) = LOG(PLA(I)%Q) 
           P(DKAL+2) = PLA(I)%EXC*PLA(I)%CW
           P(DKAL+3) = PLA(I)%EXC*PLA(I)%SW
           P(DKAL+4) = PLA(I)%TI*PLA(I)%CP
           P(DKAL+5) = PLA(I)%TI*PLA(I)%SP
           IF (RADVEL) THEN
              P(DKAL+4) = PLA(I)%TI*PLA(I)%CO
              P(DKAL+5) = PLA(I)%TI*PLA(I)%SO
           END IF
           P(DKAL+4) = PLA(I)%TI*PLA(I)%CO
           P(DKAL+5) = PLA(I)%TI*PLA(I)%SO
           P(DKAL+6) = LOG(PLA(I)%PER)
           CALL RANDOM_NUMBER(R)
           MASS = (DPI/PLA(I)%PER)**2*PLA(I)%Q**3
           HH2 = -(PLA(I)%EXC-1.d0)*MASS/PLA(I)%Q
c...                                hh2 = -2xEnergie = alpha pour kepu          
           M0 =  (STAR%T0-PLA(I)%TP)*(2.d0*R-1.d0)
           CALL KEP_UNIV(M0,PLA(I)%Q,MASS,HH2,S,C0,SC1,S2C2,S3C3)
           P(DKAL+7) = S
        END DO
        IF (RADVEL) THEN
           P(NEL*NPLA+1) = STAR%V0
           IF (JITNUM.EQ.1) P(NEL*NPLA+2) = LOG(STAR%SIGJV)
        END IF
        IF (MULTIPLA) P(NPAR) = LOG(STAR%MASS)        

        END    
        
C     
C -----------------------------------------------------------------------------
C    Calculation of FMAP & derivatives (HCI part only)
C        MAP = exp(-Chi2/2)*f(parameters) (f = priors)
C     => Replace Chi2 with Chi2-2*ln(f(parameters)) = -2*ln(MAP) = Chi2+FMAP    
C-----------------------------------------------------------------------------
C

        SUBROUTINE MRQMAPASTROM(N,P,FMAP,DFMAP,D2FMAP)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::   N               ! Dimension of parameter space
        REAL*8, DIMENSION(N,N) :: D2FMAP ! Second derivatives
        REAL*8, DIMENSION(N) :: DFMAP        ! First derivatives
        REAL*8, DIMENSION(NPLA) :: MDYN         ! Dynamical masses 
        REAL*8, DIMENSION(N) :: P            ! Parameters
c        REAL*8, DIMENSION(N) :: DLMUF        ! dln(mass)/dpars
c        REAL*8, DIMENSION(N,N) :: D2LMUF  ! d2ln(mass)/dpars^2
        REAL*8 ::       FMAP,           ! Merit function
     &                  MSTAR,          ! Stellar mass
     &                  QQ,PP,          ! tan(i/2)*[cos(O),sin(O)]
     &                  FM,DFM,D2FM,    ! Mass function & derivatives
     &                  MUF,LMUF,       ! Mass prior & lna
     &                  T2,T,           ! tan(i/2)
     &                  SI,CI,          ! sin(i),cos(i)
     &                  FF1,FF2,        ! Auxiliary (second der.)
     &                  Q,PER,NN,       ! Periastron, Period/q, mean-motion/q
     &                  LNQ,LNP,        ! ln(periastron, period / q)
     &                  CO,SO           ! cos,sin (Omega-omega)
        REAL*8, DIMENSION(:), ALLOCATABLE :: DLMUF    ! dln(mass)/dpars
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: D2LMUF  ! d2ln(mass)/dpars^2
        INTEGER*4 ::    I,J,K,KM,IM          
c
        ALLOCATE(DLMUF(N))
        ALLOCATE(D2LMUF(N,N))
        IF (MULTIPLA) MSTAR = EXP(P(N))

        DO I = 1,NPLA
          IM = NEL*(I-1)
          LNQ = P(IM+1)
          PER = EXP(P(IM+6))
          Q = EXP(LNQ)
          NN = DPI/PER
          MDYN(I) = NN*NN*(Q**3)           
        END DO

c...  Taking into account priors on masses in FMAP (-2*ln(prior(m)))    
        DLMUF = 0.d0
        D2LMUF = 0.d0
        DO J = 0,NPRIOR            
           MUF = SUM(MPRIOR(J)%BCOF(1:NPLA)*MDYN(1:NPLA))
           IF (MULTIPLA) MUF = MUF+MPRIOR(J)%BCOF(0)*MSTAR
           CALL FMAP_PRIOR(J,MUF,FM,DFM,D2FM)
           IF (MULTIPLA) THEN
              FF1 = MPRIOR(J)%BCOF(0)*MSTAR/MUF
              DLMUF(N) =  FF1
              D2LMUF(N,N) = FF1*(1.d0-FF1)
           END IF
           DO I = 1,NPLA
              IM = NEL*(I-1)
              FF1 = MPRIOR(J)%BCOF(I)*MDYN(I)/MUF
              DLMUF(IM+1) = 3.d0*FF1
              DLMUF(IM+6) = -2.d0*FF1
              FF2 = FF1*(1.d0-FF1)
              D2LMUF(IM+1,IM+1) = 9.d0*FF2
              D2LMUF(IM+1,IM+6) = -6.d0*FF2
              D2LMUF(IM+6,IM+6) = 4.d0*FF2
              D2LMUF(IM+6,IM+1) = D2LMUF(IM+1,IM+6)
              DO K = 1,NPLA
                 IF (K.NE.I) THEN
                    FF2 = FF1*MPRIOR(J)%BCOF(K)*MDYN(K)/MUF
                    KM = NEL*(K-1)
                    D2LMUF(IM+1,KM+1) = -9.d0*FF2
                    D2LMUF(IM+1,KM+6) = 6.d0*FF2
                    D2LMUF(IM+6,KM+1) = 6.d0*FF2
                    D2LMUF(IM+6,KM+6) = -4.d0*FF2
                 END IF
              END DO
              IF (MULTIPLA) THEN
                 FF2 = FF1*MPRIOR(J)%BCOF(0)*MSTAR/MUF
                 D2LMUF(IM+1,N) = -3.d0*FF2
                 D2LMUF(IM+6,N) = 2.d0*FF2
                 D2LMUF(N,IM+1) = -3.d0*FF2
                 D2LMUF(N,IM+6) = 2.d0*FF2
              END IF
           END DO
           FMAP = FMAP+FM
           DFMAP = DFMAP+DFM*DLMUF
           DO I = 1,N
              D2FMAP(:,I) = D2FMAP(:,I)+DFM*D2LMUF(:,I)
     &                                 +D2FM*DLMUF(I)*DLMUF(:)
           END DO
        END DO

        DO I = 1,NPLA
          IM = NEL*(I-1)
          LNQ = P(IM+1)
          LNP = P(IM+6)
          Q = EXP(LNQ)
          NN = DPI/EXP(LNP)
          QQ = P(IM+4)
          PP = P(IM+5)
          T2 = QQ*QQ+PP*PP    
          T = SQRT(T2)        ! Tan(i/2)
          CI = (1.d0-T2)/(1.d0+T2)  ! cos(i)
          SI = 2.d0*T/(1.d0+T2)      ! sin(i)           
          CO = QQ/T
          SO = PP/T

c...  Taking into account Period and periastron prior in FMAP
c...  prior(P) = 1/P (logarithmic prior) => -2*ln(prior(P))=+2*ln(P)          
c...  Idem for periastron
          FMAP = FMAP+2.d0*(LNP+LNQ)
          DFMAP(IM+1) = DFMAP(IM+1)+2.d0
          DFMAP(IM+6) = DFMAP(IM+6)+2.d0
c...  Taking into account inc prior in FMAP (-2*ln(prior(i))=-2*ln(sin(i)))
c...        and related derivatives
          FMAP = FMAP-2.d0*LOG(SI)
          DFMAP(IM+4) = DFMAP(IM+4)-2.d0*CO*CI/T
          DFMAP(IM+5) = DFMAP(IM+5)-2.d0*SO*CI/T
          FF1 = CI/T2
c          FF2 = 2.d0*CI*CI/T2-4.d0/(T2*(1.d0+T2)*(1.d0+T2))
          FF2 = (2.d0*T2*T2-4.d0*T2-2.d0)/(T2*(1.d0+T2)**2)
          D2FMAP(IM+4,IM+4) = -2.d0*(FF2*CO*CO+FF1)
          D2FMAP(IM+4,IM+5) = D2FMAP(IM+4,IM+5)-2.d0*FF2*CO*SO
          D2FMAP(IM+5,IM+4) = -2.d0*FF2*CO*SO
          D2FMAP(IM+5,IM+5) = -2.d0*(FF2*SO*SO+FF1)
        END DO
        DEALLOCATE(DLMUF)
        DEALLOCATE(D2LMUF)
        
        END
      
       
