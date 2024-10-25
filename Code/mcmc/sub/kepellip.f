C
C       Routines for Keplerian elliptical orbits       
C     
C -----------------------------------------------------------------------------
C    Solving Kepler's equation by quartic Newton-Raphson method                
C -----------------------------------------------------------------------------
C                                                                              
        SUBROUTINE KEPLER0(M,E,U)           

        IMPLICIT NONE

        INTEGER*4, PARAMETER :: SD = 6
        REAL*8          M,      ! Mean anomaly
     &                  AM,S,SM,SME,   !
     &                  EPSI,          ! Accuracy   
     &                  E,             ! Eccentricity
     &                  U,ECU,ESU,     ! Eccentric anomaly
     &                  F,F1,F2,F3,    ! Intermediaires   
     &                  DELA           ! Correction        
        PARAMETER(      EPSI=1d-12        )                
        REAL*8, PARAMETER :: PI = 3.14159265358979323846d0
        REAL*8, PARAMETER :: DPI = 2.d0*PI
        LOGICAL         OK,QUAD  
        INTEGER         I                                 
        INTEGER*8 ::    I8
c     
        
c        QUAD = .TRUE.
        AM = M-NINT(M/DPI,KIND(I8))*DPI
        S = 1.d0
        IF (AM.LT.0.d0) S = -1.d0
        AM = ABS(AM)
c... Starter from Esmaelzadeh & Ghadiri, Int. Jour. Comput. Appl. 89, 975 (2014)
        IF (AM.LT.0.25d0) THEN
           SM = SIN(AM)
           SME = SIN(AM+E)
           U = AM+E*SM/(1.d0-SME+SM)
        ELSE IF (AM.LT.2.d0) THEN
           U = AM+E
        ELSE
           U = AM+E*(PI-AM)/(1.d0+E)
        END IF
c...
c     U = AM+E*SIN(AM)
        OK = .FALSE.
        I=0         
        DO WHILE (.NOT.OK)
           ESU = E*SIN(U)
           ECU = E*COS(U)
           F = U-ESU-AM
           F1 = 1.d0-ECU
           F2 = ESU        
           F3 = ECU
           DELA = -F/F1
c          IF (QUAD) THEN
            DELA = -F/(F1+0.5d0*DELA*F2)
            DELA = -F/(F1+0.5d0*DELA*F2+DELA*DELA*F3/6.d0)
c          END IF
          U = U+DELA
          I = I+1
          OK = ((ABS(DELA).LT.EPSI).OR.(I.GE.500))
        END DO
        U = U*S
        if (i.gt.300) then
           WRITE(SD,*)'glup',i,e,m,u,f
           stop
        end if
        
        END

C
C -----------------------------------------------------------------------------
C    Converting true longitude at reference time to time for periastron
C -----------------------------------------------------------------------------
C                                                                              
        SUBROUTINE LAMBDA_TO_TP(T0,N,EXC,EXQ,CL,SL,CW,SW,TP)           

        IMPLICIT NONE

        REAL*8 ::       T0,            ! Reference time
     &                  N,             ! Mean motion
     &                  EXC,EXQ,       ! Eccentricity & sqrt(1-e^2)
     &                  CL,SL,         ! Cos & Sin (mean longitude)
     &                  CW,SW,         ! Cos & Sin (long. of periastron)
     &                  CV,SV,         ! Cos & Sin (true anomaly)
     &                  CU,SU,         ! Cos & Sin (eccentric anomaly)   
     &                  M0,            ! Mean anomaly
     &                  TP             ! Time of periastron (output)        

        CV = CW*CL+SW*SL                ! v = lambda-w
        SV = SL*CW-CL*SW                ! 
        CU = (EXC+CV)/(1.d0+EXC*CV)     ! cos(u)
        SU = EXQ*SV/(1.d0+EXC*CV)       ! sin(u)
        M0 = ATAN2(SU,CU)-EXC*SU        ! Kepler's equation
        TP = T0-M0/N
        
        END 

C
C -----------------------------------------------------------------------------
C    Converting time for periastron to true longitude at reference time
C -----------------------------------------------------------------------------
C                                                                              
        SUBROUTINE TP_TO_LAMBDA(T0,TP,N,EXC,EXQ,CW,SW,CL,SL)           

        IMPLICIT NONE

        REAL*8 ::       T0,            ! Reference time
     &                  TP,            ! Time of periastron (output) 
     &                  N,             ! Mean motion
     &                  EXC,EXQ,       ! Eccentricity & sqrt(1-e^2)
     &                  CW,SW,         ! Cos & Sin (long. of periastron)
     &                  CV,SV,         ! Cos & Sin (true anomaly)
     &                  U,CU,SU,         ! Cos & Sin (eccentric anomaly)   
     &                  M0,            ! Mean anomaly
     &                  CL,SL          ! Cos & Sin (mean longitude) (output)       
        M0 =  (T0-TP)*N
        CALL KEPLER0(M0,EXC,U)          ! Kepler's equation
        CU = COS(U)
        SU = SIN(U)
        CV = (CU-EXC)/(1.d0-EXC*CU)     ! True anomaly
        SV = EXQ*SU/(1.d0-EXC*CU)
        CL = CW*CV-SW*SV
        SL = CW*SV+SW*CV                !  lambda = Omega+omega+v = w+v
        
        CV = CW*CL+SW*SL                ! v = lambda-w
        SV = SL*CW-CL*SW                ! 
        CU = (EXC+CV)/(1.d0+EXC*CV)     ! cos(u)
        SU = EXQ*SV/(1.d0+EXC*CV)       ! sin(u)
        M0 = ATAN2(SU,CU)-EXC*SU        ! Kepler's equation
        TP = T0-M0/N

        END       


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
     &                  A,              ! Semi-major axis
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
          UP2 = P7*P7+P6*P6  !  1/Period^2
          A = EXP(P1)
          SI2 = (P4*P4+P5*P5)*SQRT(1.d0+P2*P2+P3*P3) ! sin(i/2)^2
          SI2 = 4.d0*SI2*(1.d0-SI2)  ! sin(i)^2
          MTOT(I) = A**3*DPI*DPI*UP2 ! n^2*a^3
          
          FMAP = FMAP+0.5d0*LOG(SI2*UP2)-P1 ! *sin(i)/P -> ln
c                                       -P1 => -ln(a) (*1/a, log prior)
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

C -----------------------------------------------------------------------------
C       Just compute Chi^2
C -----------------------------------------------------------------------------
C
        SUBROUTINE CHI2SEUL(P,CHI2)

        USE DATA

        IMPLICIT NONE

        REAL*8, DIMENSION(NPLA) ::
     &                  A,              ! 1/2 grand axe
     &                  TP,             ! Tps de passage au periastre
     &                  NN,             ! Moyen mouvement
     &                  EXC,EXQ,        ! eccentricity + sqrt(1-e^2)
     &                  CW,SW,CP,SP,    ! cos,sin (omega+Omega,omega-Omega)
     &                  CI2,SI2,        ! cos^2(i/2), sin^2(i/2)
     &                  COM,SOM,CO,SO,  ! cos,sin (omega,Omega)
     &                  CI,SI,          ! cos(i), sin(i)
     &                  MUF,            ! Planet masses
     &                  MFRAC,          ! Fractional masses
     &                  AMC,AMS,        ! Partial radial velocity amplitudes 
     &                  AMPC,AMPS,      ! Same x MFRAC
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
        INTEGER*4 ::    K,I
        
        CHI2 = 0.d0
        SIGJV = 0.d0
        CALL ELEMENTS(P(1:NPAR),NN,A,EXC,EXQ,CW,SW,CI2,SI2,CP,SP,
     &                                  CI,SI,COM,SOM,CO,SO,TP,MUF)
        MFRAC = 0.d0
        IF (MULTIPLA) THEN
           SIGMA = EXP(P(NPAR))
c...  Compute fractional masses
           DO I=1,NPLA
              SIGMA = SIGMA+MUF(I)
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
           AMS = -NN*A*SI*SOM
           AMC = +NN*A*SI*EXQ*COM
           AMPS = AMS*MFRAC
           AMPC = AMC*MFRAC
        ELSE
           E1(1,:) = CI2*CW+SI2*CP
           E1(2,:) = CI2*SW-SI2*SP
           E2(1,:) = -CI2*SW-SI2*SP
           E2(2,:) = CI2*CW-SI2*CP           
        END IF

c...  Contribution of HCI (type 1 data) to Chi2        
        DO I = 1,NPLA
          DO K = 1,PLA(I)%NDATAS
             CALL POSFITS(PLA(I)%TAS(K),
     &                    NN,A,EXC,EXQ,E1,E2,TP,MFRAC,POSX,POSY)
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
                 CALL VFITS(STAR%TVR(K),NN,V0,AMPS,AMPC,EXC,TP,VRAD)
                 DV = STAR%V(K)-VRAD
                 CHI2 = CHI2 + DV*DV/(STAR%SIGV2(K)+SIGJV*SIGJV)           
              END DO
           END IF
           IF (ISDATA(3)) THEN              
              DO I = 1,NPLA
                 DO K = 1,PLA(I)%NDATVR
                    CALL VFITSP(I,PLA(I)%TVR(K),NN,AMPS,AMPC,
     &                              AMS(I),AMC(I),EXC,TP,VRAD)
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

        SUBROUTINE ELEMENTS(P,NN,A,EXC,EXQ,CW,SW,CI2,SI2,CP,SP,
     &                                  CI,SI,COM,SOM,CO,SO,TP,MUF)
      
        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    DKAL,I
        REAL*8 ::       P(NPAR)         ! Parameters
        REAL*8, DIMENSION(NPLA) ::
     &                  A,              ! Semi-major axes
     &                  NN,             ! Mean motions
     &                  EXC,EXQ,        ! Eccentricities
     &                  CI,SI,          ! cos,sin(inclinations)
     &                  CW,SW,          ! cos,sin(w=Omega+omega)
     &                  CP,SP,          ! cos,sin(phi=Oomega-Omega)
     &                  COM,SOM,        ! cos,sin(omega's)
     &                  SI2,CI2,        ! sin(i/2)^2, cos(i/2)^2
     &                  CO,SO,          ! cos,sin(Omega's)             
     &                  TP,             ! Times for periastron passages
     &                  MUF             ! Planetary masses
        REAL*8 ::       MTOT,           ! Total mass (up to planet #i)
     &                  PER,UP,UP2,     ! Period 
     &                  SQEXQ,          ! SQRT(EXQ) = (1-e^2)^(1/4)
     &                  S2,S,SII,       ! Intermediaires
     &                  P1,P2,P3,P4,P5,P6,P7, ! Intermédiaires
     &                  M0,             ! Anomalie moyenne de reference
     &                  CL,SL,          ! Lambda = Omega+omega+v
     &                  CV,SV,          ! v = anomalie vraie
     &                  CU,SU           ! Anomalie excentrique 
        
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
           A(I) = EXP(P1)
           UP2 = P7*P7+P6*P6
           UP = SQRT(UP2)
           NN(I) = DPI*UP
           PER = 1.d0/UP
           MUF(I) = NN(I)*NN(I)*A(I)**3 ! Dynamical mass of orbit #i
           S2 = (P2*P2+P3*P3)
           S = SQRT(S2)        ! e/sqrt(1-e^2)
           EXC(I) = S/SQRT(1.d0+S2)
           EXQ(I) = EXC(I)/S
           SQEXQ = SQRT(EXQ(I))
           CW(I) = P2/S                 ! w = Omega+omega 
           SW(I) = P3/S 
           SI2(I) = (P4*P4+P5*P5)/EXQ(I)         ! sin(i/2)^2 
           CI2(I) = 1.d0-SI2(I)               ! cos(i/2)^2
           SII = SQRT(SI2(I))      ! sin(i/2)
           IF (RADVEL) THEN
              CI(I) = 1.d0-2.d0*SI2(I)              ! cos(i)
              SI(I) = 2.d0*SII*SQRT(CI2(I))       ! sin(i)
              CO(I) = P4/SII/SQEXQ ! cos(Omega)
              SO(I) = P5/SII/SQEXQ ! sin(Omega)
              COM(I) = CW(I)*CO(I)+SW(I)*SO(I) ! omega = w-Omega
              SOM(I) = SW(I)*CO(I)-CW(I)*SO(I)         
              CP(I) = 0.d0
              SP(I) = 0.d0
           ELSE
              CP(I) = P4/SII/SQEXQ ! cos(phi)  phi=omega-Omega
              SP(I) = P5/SII/SQEXQ ! sin(phi)
              CI(I) = 0.d0
              SI(I) = 0.d0
              CO(I) = 0.d0
              SO(I) = 0.d0
              COM(I) = 0.d0
              SOM(I) = 0.d0
           END IF
           CL = PER*P6       ! lambda = Omega+omega+v0 = w+v0
           SL = PER*P7

           CALL LAMBDA_TO_TP(STAR%T0,NN(I),EXC(I),EXQ(I),
     &           CL,SL,CW(I),SW(I),TP(I))
        END DO
        DO I = NPLA,2,-1
          MUF(I) = MUF(I)-MUF(I-1)           
        END DO
        IF (MULTIPLA) MUF(1) = MUF(1)-MTOT
        
        END    
C
C-----------------------------------------------------------------------------
C       Computation of priors ratio
C -----------------------------------------------------------------------------
C

        SUBROUTINE RATIO_PRIORS(P,PTRY,RAPQ,BAD)

        USE DATA

        IMPLICIT NONE
     
        LOGICAL ::      BAD           ! True if unreal solution
        INTEGER*4 ::    I,J,DKAL,K
        REAL*8, DIMENSION(NPAR) :: P, ! Parameters
     &                  PTRY         ! Trial model
        REAL*8, DIMENSION(NPLA) :: MTOT,MTOTT ! Masses
        REAl*8 ::
     &       RAPQ,              ! q(x|x')/q(x'|x) (eq.42 Ford 06) (output)
     &       S2,S,              ! Intermediaires
     &       P1,P2,P3,P4,P5,P6,P7, ! Intermédiaires
     &       PT1,PT2,PT3,PT4,PT5,PT6,PT7, !
     &       EXC,EXQ,EXCT,EXQT, ! excentricite
     &       CW,SW,CWT,SWT,     ! w = Omega+omega
     &       CL,SL,             ! Lambda = w+v = Omega+omega+v
     &       SI2,               ! sin(i/2)^2
     &       CV,SV,             ! V = anomalie vraie
     &       UP,UPT,            ! Inverses de périodes
     &       SIGJV,SIGJVT,      ! Velocity Jitters
     &       A, AT,             ! Semi-major axes
     &       MUF,MUFT,          ! Fractional masses & ln
     &       MSTAR, MSTART,     ! Stellar mass 
     &       FAC,FACT           ! Jacobiens

c     For each planet,
c  Jabobian J (orbite->params) = 1/2*e*sin(i)/((1-e^2)^3/2*P^3*a)*dnu/dM 
c     = 1/2*e*sin(i)*(1+e*cos(nu))^2/((1-e^2)^3*P^3*a)
c  Prior/planet (default) = sin(i)*1/a*1/P   (prior log in a & in P)
c     => Jacobian / Prior = 1/2*e*(1+e*cos(nu))^2/((1-e^2)^3*P^2)= (J/P)_0
c           N.B. : sin(i) and a vanished in Jacobian / Prior
c     If prior u(mu) (e.g. Gaussian) then
c               mu = sum(b_i*M_i) = sum(b_i*n_i^2*a_i^3)
c     d(mu)/d(ln(a_i))=3*M_i*b_i   d(mu)/d(ln(P_i))=-2*M_i*b_i
c     Jacobian (mu<-a) = 3*(M_i*b_i)/a_i si autres variables inchangées
c     [Jacobian (mu<-P) =-2*(M_i*b_i)/P_i si autres variables inchangées]
c      With mu<->a :
c     Jacobian / Prior = (J/P)_0*prior(a)/J(a->mu)*1/prior(mu)
c                      = (J/P)_0*1/a*(a/(3*M_i*b_i))*1/u(mu)     
c     RAPQ = (J/p)/(J'/p')=(e/e')*(FAC/FAC')*((1-e'^2)^3/(1-e^2)^3)
c                      *(P'^2/P^2)*(M'_i/M_i)*(u(mu')/u(mu))
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
          A = EXP(P1)
          AT = EXP(PT1)
          UP = P7*P7+P6*P6
          UPT = PT7*PT7+PT6*PT6
          MTOT(I) = A**3*DPI*DPI*UP
          MTOTT(I) = AT**3*DPI*DPI*UPT          
          UP = SQRT(UP)
          UPT = SQRT(UPT)
          BAD = BAD.OR.(VPRIOR(1)%BOUND(1)*UPT.GT.1.d0)
     &             .OR.(VPRIOR(1)%BOUND(2)*UPT.LT.1.d0)
          BAD = BAD.OR.(AT.LT.VPRIOR(2)%BOUND(1))
     &             .OR.(AT.GT.VPRIOR(2)%BOUND(2))

          S2 = P2*P2+P3*P3
          S = SQRT(S2)        ! e/sqrt(1-e^2)
          EXC = S/SQRT(1.d0+S2)
          EXQ = EXC/S
          CW = P2/S                 ! w = Omega+omega 
          SW = P3/S       

          S2 = PT2*PT2+PT3*PT3
          S = SQRT(S2)        ! e/sqrt(1-e^2)
          EXCT = S/SQRT(1.d0+S2)
          EXQT = EXCT/S
          CWT = PT2/S                 ! w = Omega+omega 
          SWT = PT3/S    
          BAD = BAD.OR.(EXCT.GT.VPRIOR(3)%BOUND(2)).OR.
     &             (EXCT.LT.VPRIOR(3)%BOUND(1)).OR.(EXCT.GT.1.d0)

          SI2 = (PT4*PT4+PT5*PT5)/EXQT ! sin(i/2)^2
          BAD = BAD.OR.(SI2.GT.1.d0)

          CL = P6/UP
          SL = P7/UP
          CV = CL*CW+SL*SW    ! v = lambda-w = (Omega+omega+v)-(Omega+omega)  
          FAC = (1.d0+EXC*CV)
          FAC = FAC*FAC

          CL = PT6/UPT
          SL = PT7/UPT
          CV = CL*CWT+SL*SWT 
          FACT = (1.d0+EXCT*CV)
          FACT = FACT*FACT              

          RAPQ = RAPQ*EXC/EXCT       ! e/e'
          S = EXQT/EXQ                  ! sqrt(1-e'^2)/sqrt(1-e^2)
          S = S*S
          S = S*S*S  !  *(EXQT/EXQ)          ! (1-e^'2)/(1-e^2)^(6/2)
          RAPQ = RAPQ*(FAC/FACT)*S 
          S = UP/UPT            ! (P'/P)
          S = S*S                 ! (P'/P)^2 si prior ln(P)
c       S = S*S*S               ! (P'/P)^3 si prior P
c        S = S**(10.d0/3.d0)  ! (P'/P)^(10/3) si prior a
          RAPQ = RAPQ*S
       END DO
c...  At this point RAPQ =
c...    Product((e/e')*(FAC/FAC')*((1-e'^2)^3/(1-e^2)^3)*(P'^2/P^2),planets)
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


C
C -----------------------------------------------------------------------------
C       Fit of the velocity without derivatives
C -----------------------------------------------------------------------------
C

        SUBROUTINE VFITS(TT,N,V0,AMPS,AMPC,EXC,TP,VRAD)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 :: I
        REAL*8, DIMENSION(NPLA) ::
     &                  TP,             ! Time for periastron passage
     &                  N,              ! Mean motion
     &                  EXC,            ! Eccentricity
     &                  AMPC,AMPS       ! Partial amplitudes
        REAL*8 ::       TT,             ! Time
     &                  VRAD,           ! Vitesse radiale
     &                  M,              ! Anomalie moyenne
     &                  V0,             ! Offset velocity
     &                  U,CU,SU        ! Anomalie excentrique
c
c     For each planet, AMPS = n*a*sin(i)*mfrac*sin(omega)
c                      AMPC = -n*a*sin(i)*mfrac*sqrt(1-e^2)*cos(omega)         
c     vrad = V0+sum((amps(i)*sin(u)+ampc(i)*cos(u))/(1-e(i)*cos(u))
        
        VRAD = V0
        DO I = 1,NPLA
           M = N(I)*(TT-TP(I))
           CALL KEPLER0(M,EXC(I),U)
           CU = COS(U)
           SU = SIN(U)
           VRAD = VRAD+(AMPS(I)*SU+AMPC(I)*CU)/(1.d0-EXC(I)*CU)
        END DO
        END

C
C -----------------------------------------------------------------------------
C       Fit of relative velocity of a planet / star without derivatives
C -----------------------------------------------------------------------------
C

        SUBROUTINE VFITSP(IPLA,TT,N,AMPS,AMPC,AMS,AMC,EXC,TP,VRAD)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    IPLA,            ! PLanet to be considered
     &                  I
        REAL*8, DIMENSION(NPLA) ::
     &                  TP,             ! Time for periastron passage
     &                  N,              ! Mean motion
     &                  EXC,            ! Eccentricity
     &                  AMPC,AMPS       ! Partial amplitudes
        REAL*8 ::       TT,             ! Time
     &                  VRAD,           ! Vitesse radiale
     &                  AMS,AMC,        ! Amplitudes
     &                  M,              ! Anomalie moyenne
     &                  V0,             ! Offset velocity
     &                  U,CU,SU         ! Anomalie excentrique
c
c     For each planet, AMPS = n*a*sin(i)*mfrac*sin(omega)
c                      AMPC = -n*a*sin(i)*mfrac*sqrt(1-e^2)*cos(omega)         
c                AMS, AMC = Same as AMPS,AMPC without mfrac
c     vrad = -sum((amps(i)*sin(u)+ampc(i)*cos(u))/(1-e(i)*cos(u)),i=1..ipla-1
c            -ams(ipla)*sin(u)+amc(ipla)*cos(u))/(1-e(ipla)*cos(u)),       
        VRAD = 0.d0
        DO I = 1,IPLA-1
           M = N(I)*(TT-TP(I))
           CALL KEPLER0(M,EXC(I),U)
           CU = COS(U)
           SU = SIN(U)
           VRAD = VRAD-(AMPS(I)*SU+AMPC(I)*CU)/(1.d0-EXC(I)*CU)
        END DO
        M = N(IPLA)*(TT-TP(IPLA))
        CALL KEPLER0(M,EXC(IPLA),U)
        CU = COS(U)
        SU = SIN(U)
        VRAD = VRAD-(AMS*SU+AMC*CU)/(1.d0-EXC(IPLA)*CU)

        END
C
C -----------------------------------------------------------------------------
C       Fit of the position without derivatives
C -----------------------------------------------------------------------------
C

        SUBROUTINE POSFITS(TT,N,A,EXC,EXQ,E1,E2,TP,MFRAC,POSX,POSY)

        USE DATA

        IMPLICIT NONE

        REAl*8, DIMENSION(NPLA) ::
     &                  JACX,JACY,      ! Jacobi positions
     &                  A,              ! Semi-major axis
     &                  TP,             ! Time for periastron passage
     &                  MFRAC,          ! Fractional planetary masses
     &                  N,              ! Mean motion
     &                  POSX,POSY,      ! Heliocentric positions
     &                  EXC,EXQ         ! Eccentricity + sqrt(1-e^2)
        REAL*8, DIMENSION(2,NPLA) :: E1,E2 ! Base vectors for orbits 
        REAL*8 ::       TT,             ! Time
     &                  M,M0,           ! Anomalie moyenne
     &                  U,CU,SU,        ! Anomalie excentrique
     &                  GX,GY,          ! Partial center of mass position
     &                  RX,RY           ! Coordonnes dans la base propre
        INTEGER*4 ::    I
        
        GX = 0.d0
        GY = 0.d0
        DO I = 1,NPLA
           M = N(I)*(TT-TP(I))
           CALL KEPLER0(M,EXC(I),U)
           CU = COS(U)
           SU = SIN(U)

           RX = A(I)*(CU-EXC(I))
           RY = A(I)*EXQ(I)*SU
           JACX(I) = RX*E1(1,I)+RY*E2(1,I)
           JACY(I) = RX*E1(2,I)+RY*E2(2,I)
c           POSX(I) = JACX(I)+SUM(MFRAC(1:I-1)*JACX(1:I-1))
c           POSY(I) = JACY(I)+SUM(MFRAC(1:I-1)*JACY(1:I-1))
           POSX(I) = JACX(I)+GX
           POSY(I) = JACY(I)+GY
c           print*,'jacs',sngl(tt),i,sngl(m),sngl(tp(i))
           IF (MULTIPLA) THEN
             GX = GX+MFRAC(I)*JACX(I)
             GY = GY+MFRAC(I)*JACY(I)
           END IF
        END DO

        END

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
        REAL*8 :: LAMBDA,CL,SL    ! Mean longitude at T0
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
              P(DKAL+1) = LOG(PLA(I)%A) 
              P(DKAL+2) = PLA(I)%EXC*PLA(I)%CW/PLA(I)%EXQ
              P(DKAL+3) = PLA(I)%EXC*PLA(I)%SW/PLA(I)%EXQ              
              P(DKAL+4) = PLA(I)%TI*PLA(I)%CP
              P(DKAL+5) = PLA(I)%TI*PLA(I)%SP
              IF (RADVEL) THEN
                 P(DKAL+4) = PLA(I)%TI*PLA(I)%CO
                 P(DKAL+5) = PLA(I)%TI*PLA(I)%SO
              END IF
              LAMBDA = DPI/PLA(I)%PER*(STAR%T0-PLA(I)%TP)+PLA(I)%W/DR
              CL = COS(LAMBDA)
              SL = SIN(LAMBDA)
              P(DKAL+6) = CL/PLA(I)%PER
              P(DKAL+7) = SL/PLA(I)%PER
           END DO  
           IF (MULTIPLA) P(NPAR) = LOG(STAR%MASS)
           IF (ISDATA(2)) P(NPAR-1) = STAR%V0
           IF (ISDATA(2).AND.(JITNUM.EQ.1)) THEN
              P(NPAR-1) = LOG(STAR%MASS)
c...                     NPAR-1 because no Jitter at first round             
              P(NPAR-2) = STAR%V0
           END IF
           WRITE(SD,*)'Range of periods : (days)'
           READ(5,*)VPRIOR(1)%BOUND(1:2)
           WRITE(SD,*)'Range of semi-major axes : (days)'
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
        REAL*8, DIMENSION(NPLA) :: MDYN,PER ! Dynamical masses & Periods 
        REAL*8, DIMENSION(NPLA) :: MFRAC ! Fractional masses
        INTEGER*4 ::    I,J,L,DKAL ! Indexes
        REAL*8 ::       TT,             ! Temps
     &                  MTOT,MSTAR,     ! Total mass up to planet #i
     &                  Z,A,            ! 1/2 grand axe et ln()
     &                  TP,             ! Tps de passage au periastre
     &                  N,              ! Mean motion
     &                  K,H,            ! e*[cos(w),sin(w)]/sqrt(1-e^2)
     &                  QQ,PP,          ! tan(i/2)*cos(phi),tan(i/2)*sin(phi)
     &                  DJACXDZ,DJACXDP6,DJACXDP7,
     &                  DJACXDK,DJACXDH,DJACXDP,DJACXDQ,
     &                  DJACYDZ,DJACYDP6,DJACYDP7,
     &                  DJACYDK,DJACYDH,
     &                  DJACYDP,DJACYDQ,! All derivatives
     &                  M,M0,           ! Anomalie moyenne
     &                  TI,S,T2,S2,     ! tan(i/2), e/sqrt(1-e^2) + carres
     &                  EXC,EXQ,EXQ3,   ! excentricite + sqrt(1-e^2)
     &                  CW,SW,CP,SP,    ! cos,sin (w,phi)
     &                  CO,SO,          ! cos,sin (O)
     &                  CI,SI,COM,SOM,  ! cos,sin (i,Om)
     &                  CL,SL,CM,SM,       ! cos,sin (lambda,M0)
     &                  CI2,SI2,        ! cos(i/2)^2, sin(i/2)^2
     &                  E1(2),E2(2),    ! Vecteurs de la base propre (X,Y)
     &                  U,CU,SU,        ! Anomalie excentrique
     &                  FACT,           ! 1/(1-e*cos(u))
     &                  RX,RY,          ! Coordonnes dans la base propre
     &                  DOMDK,DOMDH,DODP,DODQ,    ! Derivees interm.
     &                  DIDP,DIDQ,DOMDP,DOMDQ,    !  
     &                  DWDK,DWDH,DPHIDP,DPHIDQ, ! Derivees interm.
     &                  DCI2DP,DCI2DQ,DSI2DP,DSI2DQ, !  
     &                  DMFDM,DMFDMM,DUDH,DUDK,    !
     &                  DE1DK(2),DE1DH(2),DE1DP(2),DE1DQ(2), ! Derivees de E1
     &                  DE2DK(2),DE2DH(2),DE2DP(2),DE2DQ(2),! Derivees de E2
     &                  DEDK,DEDH,DEXQDK,DEXQDH,! Derivees de Exc,Exq
     &                  DUDN,DUDE,DUDP6,DUDP7,! Derivees de U
     &                  NDTPDK,NDTPDH,NDTPDP6,NDTPDP7, ! 
     &                  DRXDA,DRXDU,DRXDE, ! Derivees de RX (interm.)
     &                  DRXDP6,DRXDP7,DRXDK,DRXDH,DRXDP,DRXDQ,DNDP6,
     &                  DRYDA,DRYDU,DRYDE, ! Derivees de RY (interm.)
     &                  DRYDP6,DRYDP7,DRYDK,DRYDH,DRYDP,DRYDQ,DNDP7,
     &                  DJACXDA,DJACXDPER, ! Derivees interm. de JACX
     &                  DJACYDA,DJACYDPER ! Derivees interm. de JACY
        
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
          PER(I) = 1.d0/SQRT(P(DKAL+6)**2+P(DKAL+7)**2)
          CL = PER(I)*P(DKAL+6)
          SL = PER(I)*P(DKAL+7)
          N = DPI/PER(I)         ! Mean motion
          A = EXP(Z)            ! Semi-major axis
          MDYN(I) = N*N*(A**3)    ! Dynamical mass
          IF (MULTIPLA) THEN
             MFRAC(I) = (MDYN(I)-MTOT)/MDYN(I) ! Mfrac(i) = (M(i)-M(i-1))/M(i)
             MTOT = MDYN(I)
          END IF
          S2 = K*K+H*H
          S = SQRT(S2)        ! e/sqrt(1-e^2)
          EXC = S/SQRT(1.d0+S2)
          EXQ = EXC/S
          CW = K/S                   ! cos(w)
          SW = H/S                   ! sin(w)
          CM = CW*CL+SW*SL ! M0 = lambda-w
          SM = CW*SL-SW*CL
          M0 = ATAN2(SM,CM)
          TP = -M0/N+STAR%T0 ! Tp = T0-M0/NN
          M = N*(TT-TP)         ! Mean anomaly

          CALL KEPLER0(M,EXC,U)
          CU = COS(U)
          SU = SIN(U)
          FACT = 1.d0/(1.d0-EXC*CU) 
          EXQ3 = EXQ*EXQ*EXQ

          DNDP6 = DPI*CL
          DNDP7 = DPI*SL
          
          DWDK = -H/S2
          DWDH = K/S2
          NDTPDK = DWDK
          NDTPDH = DWDH
          NDTPDP6 = +PER(I)*SL+M0*PER(I)*CL
          NDTPDP7 = -PER(I)*CL+M0*PER(I)*SL                    

          DEDK = EXQ3*CW
          DEDH = EXQ3*SW
          DUDN = FACT*(TT-TP)
          DUDP6 = DUDN*DNDP6
          DUDP7 = DUDN*DNDP7

          DUDE = FACT*SU
          DUDK = -NDTPDK*FACT
          DUDH = -NDTPDH*FACT
          DUDP6 = DUDP6-NDTPDP6*FACT
          DUDP7 = DUDP7-NDTPDP7*FACT                   

          RX = A*(CU-EXC)
          RY = A*EXQ*SU
          
          T2 = QQ*QQ+PP*PP        
          TI = SQRT(T2)        ! Tan(i/2)
          IF (RADVEL) THEN
             SI = 2.d0*TI/(1.d0+T2)      ! sin(i) 
             CI = (1.d0-T2)/(1.d0+T2)  ! cos(i)
             CO = QQ/TI         ! cos(Omega)
             SO = PP/TI         ! sin(Omega)
             COM = CW*CO+SW*SO  ! omega = w-Omega
             SOM = SW*CO-CW*SO  !
          
             DOMDK = -H/S2 ! dw/dk
             DOMDH = K/S2  ! dw/dh
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

c          print*,'jac',sngl(tt),i,sngl(m),sngl(tp)
          
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
          
          DRXDA = CU-EXC
          DRXDU = -A*SU
          DRXDE = -A
          DRXDP6 = DRXDU*DUDP6
          DRXDP7 = DRXDU*DUDP7
          DRXDK = (DRXDU*DUDE+DRXDE)*DEDK+DRXDU*DUDK
          DRXDH = (DRXDU*DUDE+DRXDE)*DEDH+DRXDU*DUDH
          DRXDP = 0.d0
          DRXDQ = 0.d0

          DRYDA = EXQ*SU
          DRYDE = -A*SU*EXC/EXQ
          DRYDU = A*EXQ*CU
          DRYDP6 = DRYDU*DUDP6
          DRYDP7 = DRYDU*DUDP7
          DRYDK = (DRYDU*DUDE+DRYDE)*DEDK+DRYDU*DUDK
          DRYDH = (DRYDU*DUDE+DRYDE)*DEDH+DRYDU*DUDH
          DRYDP = 0.d0
          DRYDQ = 0.d0

          DJACXDA = DRXDA*E1(1)+DRYDA*E2(1)
          DJACXDK = DRXDK*E1(1)+DRYDK*E2(1)+RX*DE1DK(1)+RY*DE2DK(1)
          DJACXDH = DRXDH*E1(1)+DRYDH*E2(1)+RX*DE1DH(1)+RY*DE2DH(1)
          DJACXDQ = DRXDQ*E1(1)+DRYDQ*E2(1)+RX*DE1DQ(1)+RY*DE2DQ(1)
          DJACXDP = DRXDP*E1(1)+DRYDP*E2(1)+RX*DE1DP(1)+RY*DE2DP(1)
          DJACXDP6 = DRXDP6*E1(1)+DRYDP6*E2(1)
          DJACXDP7 = DRXDP7*E1(1)+DRYDP7*E2(1)
          
          DJACYDA = DRXDA*E1(2)+DRYDA*E2(2)
          DJACYDK = DRXDK*E1(2)+DRYDK*E2(2)+RX*DE1DK(2)+RY*DE2DK(2)
          DJACYDH = DRXDH*E1(2)+DRYDH*E2(2)+RX*DE1DH(2)+RY*DE2DH(2)
          DJACYDQ = DRXDQ*E1(2)+DRYDQ*E2(2)+RX*DE1DQ(2)+RY*DE2DQ(2)
          DJACYDP = DRXDP*E1(2)+DRYDP*E2(2)+RX*DE1DP(2)+RY*DE2DP(2)
          DJACYDP6 = DRXDP6*E1(2)+DRYDP6*E2(2)
          DJACYDP7 = DRXDP7*E1(2)+DRYDP7*E2(2)
          
          DJACXDZ = A*DJACXDA                                  
          DJACYDZ = A*DJACYDA 
          
          DJACX(:,I) =
     &  (/ DJACXDZ,DJACXDK,DJACXDH,DJACXDQ,DJACXDP,DJACXDP6,DJACXDP7 /)
          DJACY(:,I) =
     &  (/ DJACYDZ,DJACYDK,DJACYDH,DJACYDQ,DJACYDP,DJACYDP6,DJACYDP7 /)
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
              DMFDMM = DMFDM*PER(J)*PER(J)
              DPOSX(DKAL+1,I) = DPOSX(DKAL+1,I)+3.d0*DMFDM*JACX(J)
c                                 d(ln(MDYN(J))/d(ln(a(J))=3.
              DPOSY(DKAL+1,I) = DPOSY(DKAL+1,I)+3.d0*DMFDM*JACY(J)
c                                 => d(POSX,POSY)/d(ln(a(J))
              DPOSX(DKAL+6,I) = DPOSX(DKAL+6,I)
     &             +2.d0*DMFDMM*JACX(J)*P(DKAL+6)
c              d(ln(MDYN(J))/d(P6(J))=2*P6/(P6^2+P7^2)=2*PER(J)^2*P6(J)
              DPOSY(DKAL+6,I) = DPOSY(DKAL+6,I)
     &             +2.d0*DMFDMM*JACY(J)*P(DKAL+6)
c              => d(POSX,POSY)/d(P6(J))
              DPOSX(DKAL+7,I) = DPOSX(DKAL+7,I)
     &             +2.d0*DMFDMM*JACX(J)*P(DKAL+7)
c              d(ln(MDYN(J))/d(P7(J))=2*P7/(P6^2+P7^2)=2*PER(J)^2*P7(J)
              DPOSY(DKAL+7,I) = DPOSY(DKAL+7,I)
     &             +2.d0*DMFDMM*JACY(J)*P(DKAL+7)
c              => d(POSX,POSY)/d(P7(J)         
              IF (J.GT.1) THEN
                 DMFDMM = DMFDM*PER(J-1)*PER(J-1)
                 DPOSX(DKAL-6,I) = DPOSX(DKAL-6,I)-3.d0*DMFDM*JACX(J)
c                               d(ln(MDYN(J-1))/d(ln(a(J-1))=3.
                 DPOSY(DKAL-6,I) = DPOSY(DKAL-6,I)-3.d0*DMFDM*JACY(J)
c                               => d(POSX,POSY)/d(ln(a(J-1))
                 DPOSX(DKAL-1,I) = DPOSX(DKAL-1,I)
     &                      -2.d0*DMFDMM*JACX(J)*P(DKAL-1)
c            d(ln(MDYN(J-1))/d(P6(J-1))=2*P6/(P6^2+P7^2)=2*PER(J-1)^2*P6(J-1)
                 DPOSY(DKAL-1,I) = DPOSY(DKAL-1,I)
     &                      -2.d0*DMFDMM*JACY(J)*P(DKAL-1)
c              => d(POSX,POSY)/d(P6(J-1))
                 DPOSX(DKAL,I) = DPOSX(DKAL,I)
     &                      -2.d0*DMFDMM*JACX(J)*P(DKAL)
c            d(ln(MDYN(J-1))/d(P7(J-1))=2*P7/(P6^2+P7^2)=2*PER(J-1)^2*P7(J-1)
                 DPOSY(DKAL,I) = DPOSY(DKAL,I)
     &                     -2.d0*DMFDMM*JACY(J)*P(DKAL)
c               => d(POSX,POSY)/d(P7(J-1))                
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
        REAL*8 ::       TT,             ! Time
     &                  MTOT,MSTAR,     ! Total mass (partial)
     &                  Z,A,AMP,           ! 1/2 gd axe, amplitude et ln()
     &                  TP,             ! Tps de passage au periastre
     &                  N,              ! Moyen mouvement
     &                  K,H,            ! e*[cos(w),sin(w)]/sqrt(1-e^2)
     &                  QQ,PP,          ! tan(i/2)*[cos(O),sin(O)]
     &                  VRAD,           ! Vitesse radiale
     &                  TI2,TII,        ! tan(i/2)
     &                  SI,             ! sin(i)
     &                  M,M0,           ! Anomalie moyenne
     &                  S,S2,           ! e/sqrt(1-e^2) + carre
     &                  EXC,EXQ,EXQ3,   ! excentricite + sqrt(1-e^2)
     &                  COM,SOM,CO,SO,  ! cos,sin (omega,Omega)
     &                  CL,SL,CM,SM,    ! cos,sin (lambda,M0)
     &                  CW,SW,          ! cos, sin (w=omega+Omega)
     &                  U,CU,SU,        ! Anomalie excentrique
     &                  FACT,FACT2,     ! 1/(1-e*cos(u))     
     &                  RX,RY,          ! Coordonnes dans la base propre
     &                  VITZ,           ! Vitesse en Z
     &                  DOMDK,DOMDH,DOMDQ,DOMDP, ! Derivees interm.     
     &                  DWDK,DWDH, ! Derivees interm.
     &                  DFACTU,DFACTE,     ! Derivees interm.
     &                  DMFDM,DMFDMM,DUDH,DUDK,    !
     &                  DEDK,DEDH,         ! Derivees de Exc
     &                  DUDN,DUDE,DUDP6,DUDP7,! Derivees de U
     &                  NDTPDK,NDTPDH,NDTPDP6,NDTPDP7, ! 
     &                  DVZDU,DVZDE,DVZDOM,DVZDAMP,
     &                  DNDA,DNDMUF,DNDP6,DNDP7,    ! Dérivées de N
     &                  DAMPDZ,DSI,DAMPDQ,DAMPDP,DAMPDN, ! Dérivées de AMP
     &                  DVJACP6,DVJACZ,DVJACP7,
     &                  DVJACH,DVJACK,DVJACP,DVJACQ ! Dérivées de VJAC        
        
        MSTAR = EXP(P(N1))
        MTOT = MSTAR       
        DO I = 1,NPLA
          DKAL = NEL*(I-1)
          Z = P(DKAL+1)
          K = P(DKAL+2)
          H = P(DKAL+3)
          QQ = P(DKAL+4)
          PP = P(DKAL+5)
          PER(I) = 1.d0/SQRT(P(DKAL+6)**2+P(DKAL+7)**2)
          CL = PER(I)*P(DKAL+6)
          SL = PER(I)*P(DKAL+7)
          N = DPI/PER(I)         ! Mean motion
          A = EXP(Z)            ! Semi-major axis
          MDYN(I) = N*N*(A**3)    ! Dynamical mass
          IF (MULTIPLA) THEN
             MFRAC(I) = (MDYN(I)-MTOT)/MDYN(I) ! Mfrac(i) = (M(i)-M(i-1))/M(i)
             MTOT = MDYN(I)
          END IF
          S2 = K*K+H*H
          S = SQRT(S2)        ! e/sqrt(1-e^2)
          EXC = S/SQRT(1.d0+S2)
          EXQ = EXC/S
          CW = K/S                   ! cos(w)
          SW = H/S                   ! sin(w)
          CM = CW*CL+SW*SL ! M0 = lambda-w
          SM = CW*SL-SW*CL
          M0 = ATAN2(SM,CM)
          TP = -M0/N+STAR%T0 ! Tp = T0-M0/NN
          M = N*(TT-TP)         ! Mean anomaly
          
          CALL KEPLER0(M,EXC,U)
          CU = COS(U)
          SU = SIN(U)
          FACT = 1.d0/(1.d0-EXC*CU) 
          EXQ3 = EXQ*EXQ*EXQ
          TI2 = QQ*QQ+PP*PP        
          TII = SQRT(TI2)        ! Tan(i/2)
c          CI = (1.d0-T2)/(1.d0+T2)  ! cos(i)
          SI = 2.d0*TII/(1.d0+TI2)      ! sin(i) 
          CO = QQ/TII                    ! cos(Omega)
          SO = PP/TII                    ! sin(Omega)
          COM = CW*CO+SW*SO     ! omega = w-Omega
          SOM = SW*CO-CW*SO              !
          AMP = -N*A*SI
          VJAC(I) = AMP*FACT*(SU*SOM-EXQ*CU*COM)
          
          DNDP6 = DPI*CL
          DNDP7 = DPI*SL
          
          DWDK = -H/S2
          DWDH = K/S2
          NDTPDK = DWDK
          NDTPDH = DWDH
          NDTPDP6 = +PER(I)*SL+M0*PER(I)*CL
          NDTPDP7 = -PER(I)*CL+M0*PER(I)*SL                    

          DEDK = EXQ3*CW
          DEDH = EXQ3*SW
          DUDN = FACT*(TT-TP)
          DUDP6 = DUDN*DNDP6
          DUDP7 = DUDN*DNDP7          
          
          DAMPDZ = AMP
          DAMPDN = AMP/N 
          DSI = (1.d0-TI2)/(TI2*(1.d0+TI2))
          DAMPDQ = AMP*DSI*QQ
          DAMPDP = AMP*DSI*PP

          FACT2 = FACT*FACT
          DFACTU = -EXC*SU*FACT2
          DFACTE = CU*FACT2

          DUDE = FACT*SU
          DUDK = -NDTPDK*FACT
          DUDH = -NDTPDH*FACT
          DUDP6 = DUDP6-NDTPDP6*FACT
          DUDP7 = DUDP7-NDTPDP7*FACT     

          DOMDK = DWDK ! -H/S2
          DOMDH = DWDH ! K/S2
          DOMDQ = PP/TI2 ! -dO/dK
          DOMDP = -QQ/TI2 ! -dO/dH
c          DNDA = -1.5d0*N/A
c          DNDMUF = +0.5d0*N/MTOT
  
          DVZDU = VJAC(I)/FACT*DFACTU+AMP*FACT*(CU*SOM+EXQ*SU*COM)
          DVZDE = DVZDU*DUDE+VJAC(I)/FACT*DFACTE+AMP*FACT*EXC/EXQ*CU*COM
          DVZDOM = AMP*FACT*(SU*COM+EXQ*CU*SOM)
          DVZDAMP = VJAC(I)/AMP
          DVJACP6 = DVZDU*DUDP6+DVZDAMP*DAMPDN*DNDP6
          DVJACP7 = DVZDU*DUDP7+DVZDAMP*DAMPDN*DNDP7
          DVJACZ = DVZDAMP*DAMPDZ
          DVJACK = DVZDE*DEDK+DVZDOM*DOMDK+DVZDU*DUDK
          DVJACH = DVZDE*DEDH+DVZDOM*DOMDH+DVZDU*DUDH
          DVJACQ = DVZDAMP*DAMPDQ+DVZDOM*DOMDQ
          DVJACP = DVZDAMP*DAMPDP+DVZDOM*DOMDP

          DVJAC(:,I) =
     &          (/ DVJACZ,DVJACK,DVJACH,DVJACQ,DVJACP,DVJACP6,DVJACP7 /)        

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
C       Radial velocity fit and derivatives (For Lveneberg-Marquardt)
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQVFIT(N1,TT,P,VRAD,DVRAD)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    N1              ! Dimension of parameter space 
        REAL*8, DIMENSION(NPLA) :: VJAC ! Jacobi velocities
        REAL*8, DIMENSION(7,NPLA) :: DVJAC ! Jac. deriv.
        REAL*8, DIMENSION(NPLA) :: MDYN,PER ! Dynamical masses & periods 
        REAL*8, DIMENSION(NPLA) :: MFRAC ! Fractional masses
        REAL*8, DIMENSION(N1) :: P,DVRAD ! Parameters & derivatives
        INTEGER*4 ::    I,J,DKAL
        REAL*8 ::       TT,             ! Time
     &                  MTOT,MSTAR,     ! Total mass (partial)
     &                  Z,A,AMP,           ! 1/2 gd axe, amplitude et ln()
     &                  TP,             ! Tps de passage au periastre
     &                  N,              ! Moyen mouvement
     &                  K,H,            ! e*[cos(w),sin(w)]/sqrt(1-e^2)
     &                  QQ,PP,          ! tan(i/2)*[cos(O),sin(O)]
     &                  VRAD,           ! Vitesse radiale
     &                  TI2,TII,        ! tan(i/2)
     &                  SI,             ! sin(i)
     &                  M,M0,           ! Anomalie moyenne
     &                  S,S2,           ! e/sqrt(1-e^2) + carre
     &                  EXC,EXQ,EXQ3,   ! excentricite + sqrt(1-e^2)
     &                  COM,SOM,CO,SO,  ! cos,sin (omega,Omega)
     &                  CL,SL,CM,SM,    ! cos,sin (lambda,M0)
     &                  CW,SW,          ! cos, sin (w=omega+Omega)
     &                  U,CU,SU,        ! Anomalie excentrique
     &                  FACT,FACT2,     ! 1/(1-e*cos(u))     
     &                  RX,RY,          ! Coordonnes dans la base propre
     &                  VITZ,           ! Vitesse en Z
     &                  DOMDK,DOMDH,DOMDQ,DOMDP, ! Derivees interm.     
     &                  DWDK,DWDH, ! Derivees interm.
     &                  DFACTU,DFACTE,     ! Derivees interm.
     &                  DMFDM,DMFDMM,DUDH,DUDK,    !
     &                  DEDK,DEDH,         ! Derivees de Exc
     &                  DUDN,DUDE,DUDP6,DUDP7,! Derivees de U
     &                  NDTPDK,NDTPDH,NDTPDP6,NDTPDP7, ! 
     &                  DVZDU,DVZDE,DVZDOM,DVZDAMP,
     &                  DNDA,DNDMUF,DNDP6,DNDP7,    ! Dérivées de N
     &                  DAMPDZ,DSI,DAMPDQ,DAMPDP,DAMPDN, ! Dérivées de AMP
     &                  DVJACP6,DVJACZ,DVJACP7,
     &                  DVJACH,DVJACK,DVJACP,DVJACQ ! Dérivées de VJAC        
        
        VRAD = P(NEL*NPLA+1)  ! Offset velocity
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
          PER(I) = 1.d0/SQRT(P(DKAL+6)**2+P(DKAL+7)**2)
          CL = PER(I)*P(DKAL+6)
          SL = PER(I)*P(DKAL+7)
          N = DPI/PER(I)         ! Mean motion
          A = EXP(Z)            ! Semi-major axis
          MDYN(I) = N*N*(A**3)    ! Dynamical mass
          IF (MULTIPLA) THEN
             MFRAC(I) = (MDYN(I)-MTOT)/MDYN(I) ! Mfrac(i) = (M(i)-M(i-1))/M(i)
             MTOT = MDYN(I)
          END IF
          S2 = K*K+H*H
          S = SQRT(S2)        ! e/sqrt(1-e^2)
          EXC = S/SQRT(1.d0+S2)
          EXQ = EXC/S
          CW = K/S                   ! cos(w)
          SW = H/S                   ! sin(w)
          CM = CW*CL+SW*SL ! M0 = lambda-w
          SM = CW*SL-SW*CL
          M0 = ATAN2(SM,CM)
          TP = -M0/N+STAR%T0 ! Tp = T0-M0/NN
          M = N*(TT-TP)         ! Mean anomaly
          
          CALL KEPLER0(M,EXC,U)
          CU = COS(U)
          SU = SIN(U)
          FACT = 1.d0/(1.d0-EXC*CU) 
          EXQ3 = EXQ*EXQ*EXQ
          TI2 = QQ*QQ+PP*PP        
          TII = SQRT(TI2)        ! Tan(i/2)
c          CI = (1.d0-T2)/(1.d0+T2)  ! cos(i)
          SI = 2.d0*TII/(1.d0+TI2)      ! sin(i) 
          CO = QQ/TII                    ! cos(Omega)
          SO = PP/TII                    ! sin(Omega)
          COM = CW*CO+SW*SO     ! omega = w-Omega
          SOM = SW*CO-CW*SO              !
          AMP = -N*A*SI
          VJAC(I) = AMP*FACT*(SU*SOM-EXQ*CU*COM)
          
          DNDP6 = DPI*CL
          DNDP7 = DPI*SL
          
          DWDK = -H/S2
          DWDH = K/S2
          NDTPDK = DWDK
          NDTPDH = DWDH
          NDTPDP6 = +PER(I)*SL+M0*PER(I)*CL
          NDTPDP7 = -PER(I)*CL+M0*PER(I)*SL                    

          DEDK = EXQ3*CW
          DEDH = EXQ3*SW
          DUDN = FACT*(TT-TP)
          DUDP6 = DUDN*DNDP6
          DUDP7 = DUDN*DNDP7          
          
          DAMPDZ = AMP
          DAMPDN = AMP/N 
          DSI = (1.d0-TI2)/(TI2*(1.d0+TI2))
          DAMPDQ = AMP*DSI*QQ
          DAMPDP = AMP*DSI*PP

          FACT2 = FACT*FACT
          DFACTU = -EXC*SU*FACT2
          DFACTE = CU*FACT2

          DUDE = FACT*SU
          DUDK = -NDTPDK*FACT
          DUDH = -NDTPDH*FACT
          DUDP6 = DUDP6-NDTPDP6*FACT
          DUDP7 = DUDP7-NDTPDP7*FACT     

          DOMDK = DWDK ! -H/S2
          DOMDH = DWDH ! K/S2
          DOMDQ = PP/TI2 ! -dO/dK
          DOMDP = -QQ/TI2 ! -dO/dH
c          DNDA = -1.5d0*N/A
c          DNDMUF = +0.5d0*N/MTOT
  
          DVZDU = VJAC(I)/FACT*DFACTU+AMP*FACT*(CU*SOM+EXQ*SU*COM)
          DVZDE = DVZDU*DUDE+VJAC(I)/FACT*DFACTE+AMP*FACT*EXC/EXQ*CU*COM
          DVZDOM = AMP*FACT*(SU*COM+EXQ*CU*SOM)
          DVZDAMP = VJAC(I)/AMP
          DVJACP6 = DVZDU*DUDP6+DVZDAMP*DAMPDN*DNDP6
          DVJACP7 = DVZDU*DUDP7+DVZDAMP*DAMPDN*DNDP7
          DVJACZ = DVZDAMP*DAMPDZ
          DVJACK = DVZDE*DEDK+DVZDOM*DOMDK+DVZDU*DUDK
          DVJACH = DVZDE*DEDH+DVZDOM*DOMDH+DVZDU*DUDH
          DVJACQ = DVZDAMP*DAMPDQ+DVZDOM*DOMDQ
          DVJACP = DVZDAMP*DAMPDP+DVZDOM*DOMDP

          DVJAC(:,I) =
     &          (/ DVJACZ,DVJACK,DVJACH,DVJACQ,DVJACP,DVJACP6,DVJACP7 /)        

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
     &                  DEDP,DWDP,DIDP,DMUDP,DNNDP,DPERDP,DLDP,DTPDP
        INTEGER*4 ::    I,J,DKAL
        REAL*8 ::       P(NPAR), ! Paramètres
     &                  COV(NPAR,NPAR), ! Covariances
     &                  NN,DNN,         ! Moyen mouvement
     &                  K,H,            ! e*(cos(w),sin(w))/sqrt(1-e^2)
     &                  SL,CL,          ! sin & cos(real longitude)
     &                  M0,CM,SM,       ! Mean anomaly
     &                  S2,S,T2,           ! Intermediaire
     &                  EXQ3,           ! Excentricite^3
     &                  SIGMA,          ! Cumulative mass
     &                  QQ,PP,          ! tan(i/2)*(cos(phi),sin(phi))
     &                  ERR_COMB        ! Fonction erreur combinée
        
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
           PLA(I)%A = EXP(P(DKAL+1))
           PLA(I)%DA = SQRT(COV(DKAL+1,DKAL+1))*PLA(I)%A
           PLA(I)%PER = 1.d0/SQRT(P(DKAL+6)**2+P(DKAL+7)**2)
           NN = DPI/PLA(I)%PER
           CL = P(DKAL+6)*PLA(I)%PER
           SL = P(DKAL+7)*PLA(I)%PER
           K = P(DKAL+2)
           H = P(DKAL+3)
           S2 = K*K+H*H
           S = SQRT(S2)         ! e/sqrt(1-e^2)
           PLA(I)%EXC = S/SQRT(1.d0+S2)
           PLA(I)%EXQ = PLA(I)%EXC/S          ! sqrt(1-e^2)
           EXQ3 = PLA(I)%EXQ**3
           PLA(I)%CW = K/S             ! cos(w=omega+Omega)
           PLA(I)%SW = H/S          ! sin(w)

           PLA(I)%W = ATAN2(PLA(I)%SW,PLA(I)%CW)*DR
           CM = PLA(I)%CW*CL+PLA(I)%SW*SL ! M0 = lambda-w
           SM = PLA(I)%CW*SL-PLA(I)%SW*CL
           M0 = ATAN2(SM,CM)
           PLA(I)%TP = -M0/NN+STAR%T0 ! Tp = T0-M0/NN
           DEDP(1:NPAR) = 0.d0
           DWDP(1:NPAR) = 0.d0
           DWDP(DKAL+2) = -H/S2
           DWDP(DKAL+3) = K/S2
           DEDP(DKAL+2) = EXQ3*PLA(I)%CW
           DEDP(DKAL+3) = EXQ3*PLA(I)%SW
           PLA(I)%DEXC = ERR_COMB(NPAR,DEDP(1:NPAR),COV(1:NPAR,1:NPAR))         
           PLA(I)%DW = ERR_COMB(NPAR,DWDP(1:NPAR),COV(1:NPAR,1:NPAR))*DR
           DPERDP(1:NPAR) = 0.d0
           DPERDP(DKAL+6) = -PLA(I)%PER**2*CL
           DPERDP(DKAL+7) = -PLA(I)%PER**2*SL
           PLA(I)%DPER =
     &             ERR_COMB(NPAR,DPERDP(1:NPAR),COV(1:NPAR,1:NPAR))
           DLDP(1:NPAR) = 0.d0
           DLDP(DKAL+6) = -PLA(I)%PER*SL
           DLDP(DKAL+7) = +PLA(I)%PER*CL
           DTPDP(1:NPAR) = 0.d0
           DTPDP(1:NPAR) = -DPERDP(1:NPAR)/DPI*M0
     &           +(DWDP(1:NPAR)-DLDP(1:NPAR))/NN
           PLA(I)%DTP =
     &             ERR_COMB(NPAR,DTPDP(1:NPAR),COV(1:NPAR,1:NPAR))           

           PLA(I)%MDYN = NN*NN*PLA(I)%A**3
           DMUDP = 0.d0
           DMUDP(DKAL+1) = 3.d0*PLA(I)%MDYN ! dm/m = 3*da/a-2*dP/P
           DMUDP(DKAL+6) = 2.d0*DPI*DPI*P(DKAL+6)*PLA(I)%A**3 ! d(m)/d(p6,p7)
           DMUDP(DKAL+7) = 2.d0*DPI*DPI*P(DKAL+7)*PLA(I)%A**3
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
             DMUDP(DKAL+6-NEL) =
     &             -2.d0*DPI*DPI*P(DKAL+6-NEL)*PLA(I-1)%A**3
             DMUDP(DKAL+7-NEL) =
     &             -2.d0*DPI*DPI*P(DKAL+7-NEL)*PLA(I-1)%A**3
             PLA(I)%DMU =
     &                ERR_COMB(NPAR,DMUDP(1:NPAR),COV(1:NPAR,1:NPAR))
           END IF

           QQ = P(DKAL+4)
           PP = P(DKAL+5)
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
              DWDP(DKAL+2) = -H/S2
              DWDP(DKAL+3) = K/S2           
              DWDP(DKAL+4) = PP/T2
              DWDP(DKAL+5) = -QQ/T2
              PLA(I)%DOM = ERR_COMB(NPAR,
     &                        DWDP(1:NPAR),COV(1:NPAR,1:NPAR))*DR
              DWDP(1:NPAR) = 0.d0
              DWDP(DKAL+2) = -H/S2
              DWDP(DKAL+3) = K/S2           
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
              DWDP(DKAL+2) = -0.5d0*H/S2
              DWDP(DKAL+3) = 0.5d0*K/S2           
              DWDP(DKAL+4) = -0.5d0*PP/T2
              DWDP(DKAL+5) = 0.5d0*QQ/T2
              PLA(I)%DOM = ERR_COMB(NPAR,
     &                        DWDP(1:NPAR),COV(1:NPAR,1:NPAR))*DR
              DWDP(1:NPAR) = 0.d0
              DWDP(DKAL+2) = -0.5d0*H/S2
              DWDP(DKAL+3) = 0.5d0*K/S2           
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
     &                  Z,A,            ! 1/2 gd axe et ln()
     &                  TP,             ! Tps de passage au periastre
     &                  NN,              ! Moyen mouvement
     &                  K,H,            ! e*[cos(w),sin(w)]/sqrt(1-e^2)
     &                  QQ,PP,          ! tan(i/2)*[cos(O),sin(O)]
     &                  PER,            ! Period
     &                  S,S2,           ! e/sqrt(1-e^2) + carre
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
          A = EXP(Z)            ! 1/2 gd axe       
          PER = 1.d0/SQRT(P(DKAL+6)*P(DKAL+6)+P(DKAL+7)*P(DKAL+7))
          NN = DPI/PER         ! Mean motion
          MDYN(I) = NN*NN*(A**3)    ! Dynamical mass
          IF (MULTIPLA) THEN
             MFRAC = (MDYN(I)-MTOT)/MDYN(I) ! Mfrac = (M(i)-M(i-1))/M(i)
             MTOT = MDYN(I)
c             TEST = TEST.OR.(MFRAC.LE.0.d0)
c             IF (TEST) print*,'Bahhh !'
          END IF
          S2 = K*K+H*H
          S = SQRT(S2)        ! e/sqrt(1-e^2)
          EXC = S/SQRT(1.d0+S2)
          TEST = TEST.OR.(EXC.GT.VPRIOR(3)%BOUND(2))
c          if (test) print*,'test exc',test,exc,k,h
          TEST = TEST.OR.(PER.GT.VPRIOR(1)%BOUND(2))
     &               .OR.(PER.LT.VPRIOR(1)%BOUND(1))
c          if (test) print*,'test per',test,per
          TEST = TEST.OR.(A.GT.VPRIOR(2)%BOUND(2))
     &               .OR.(A.LT.VPRIOR(2)%BOUND(1))
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
        REAL*8 ::       R,              ! Random number
     &                  LAMBDA          ! Mean longitude at T0
        
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
           P(DKAL+1) = LOG(PLA(I)%A) 
           P(DKAL+2) = PLA(I)%EXC*PLA(I)%CW/PLA(I)%EXQ
           P(DKAL+3) = PLA(I)%EXC*PLA(I)%SW/PLA(I)%EXQ              
           P(DKAL+4) = PLA(I)%TI*PLA(I)%CP
           P(DKAL+5) = PLA(I)%TI*PLA(I)%SP
           IF (RADVEL) THEN
              P(DKAL+4) = PLA(I)%TI*PLA(I)%CO
              P(DKAL+5) = PLA(I)%TI*PLA(I)%SO
           END IF
           CALL RANDOM_NUMBER(R)
           LAMBDA = DPI*R
           PLA(I)%TP = (PLA(I)%W/DR-LAMBDA)*PLA(I)%PER/DPI+STAR%T0
           P(DKAL+6) = COS(LAMBDA)/PLA(I)%PER
           P(DKAL+7) = SIN(LAMBDA)/PLA(I)%PER
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
        REAL*8, DIMENSION(NPLA) :: PER2         ! Periods^2
        REAL*8, DIMENSION(N) :: P            ! Parameters
        REAL*8, DIMENSION(N) :: DLMUF        ! dln(mass)/dpars
        REAL*8, DIMENSION(N,N) :: D2LMUF  ! d2ln(mass)/dpars^2
        REAL*8 ::       FMAP,           ! Merit function
     &                  MSTAR,          ! Stellar mass
     &                  QQ,PP,          ! tan(i/2)*[cos(O),sin(O)]
     &                  FM,DFM,D2FM,    ! Mass function & derivatives
     &                  MUF,LMUF,       ! Mass prior & lna
     &                  T2,T,           ! tan(i/2)
     &                  SI,CI,          ! sin(i),cos(i)
     &                  FF1,FF2,        ! Auxiliary (second der.)
     &                  A,PER,NN,       ! Semi-major axis, Period, mean-motion
     &                  LNA,LNP,        ! ln(semi-major axis)
     &                  P6,P7,          ! cos(lambda),sin(lambda)/P
     &                  CO,SO           ! cos,sin (Omega-omega)

        INTEGER*4 ::    I,J,K,KM,IM          
c
        IF (MULTIPLA) MSTAR = EXP(P(N))

        DO I = 1,NPLA
          IM = NEL*(I-1)
          LNA = P(IM+1)
          P6 = P(IM+6)
          P7 = P(IM+7)
          PER2(I) = 1.d0/(P6*P6+P7*P7)
          A = EXP(LNA)
          NN = DPI/SQRT(PER2(I))
          MDYN(I) = NN*NN*(A**3)           
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
              DLMUF(IM+6) = 2.d0*FF1*PER2(I)*P(IM+6)
              DLMUF(IM+7) = 2.d0*FF1*PER2(I)*P(IM+7)
              FF2 = FF1*(1.d0-FF1)
              D2LMUF(IM+1,IM+1) = 9.d0*FF2
              D2LMUF(IM+1,IM+6) = 6.d0*FF2*PER2(I)*P(IM+6)
              D2LMUF(IM+1,IM+7) = 6.d0*FF2*PER2(I)*P(IM+7)
              D2LMUF(IM+6,IM+6) = 2.d0*FF1*PER2(I)*
     &             (1.d0-2.d0*FF1*PER2(I)*P(IM+6)**2)
              D2LMUF(IM+7,IM+7) = 2.d0*FF1*PER2(I)*
     &             (1.d0-2.d0*FF1*PER2(I)*P(IM+7)**2)              
              D2LMUF(IM+6,IM+7) = -4.d0*FF1*FF1*PER2(I)*PER2(I)*
     &             P(IM+6)*P(IM+7)
              D2LMUF(IM+7,IM+6) = D2LMUF(IM+6,IM+7)
              D2LMUF(IM+6,IM+1) = D2LMUF(IM+1,IM+6)
              D2LMUF(IM+7,IM+1) = D2LMUF(IM+1,IM+7)
              DO K = 1,NPLA
                 IF (K.NE.I) THEN
                    FF2 = FF1*MPRIOR(J)%BCOF(K)*MDYN(K)/MUF
                    KM = NEL*(K-1)
                    D2LMUF(IM+1,KM+1) = -9.d0*FF2
                    D2LMUF(IM+1,KM+6) = -6.d0*PER2(K)*P(KM+6)
                    D2LMUF(IM+1,KM+7) = -6.d0*PER2(K)*P(KM+7)
                    D2LMUF(IM+6,KM+1) = -6.d0*FF2*PER2(I)*P(IM+6)
                    D2LMUF(IM+7,KM+1) = -6.d0*FF2*PER2(I)*P(IM+7)
                    D2LMUF(IM+6,KM+6) = -4.d0*FF2*PER2(I)*PER2(K)*
     &                                              P(IM+6)*P(KM+6)
                    D2LMUF(IM+7,KM+7) = -4.d0*FF2*PER2(I)*PER2(K)*
     &                                              P(IM+7)*P(KM+7)
                    D2LMUF(IM+6,KM+7) = -4.d0*FF2*PER2(I)*PER2(K)*
     &                                              P(IM+6)*P(KM+7)
                    D2LMUF(IM+7,KM+6) = -4.d0*FF2*PER2(I)*PER2(K)*
     &                                              P(IM+7)*P(KM+6)
                 END IF
              END DO
              IF (MULTIPLA) THEN
                 FF2 = FF1*MPRIOR(J)%BCOF(0)*MSTAR/MUF
                 D2LMUF(IM+1,N) = -3.d0*FF2
                 D2LMUF(IM+6,N) = -2.d0*FF2*PER2(I)*P(IM+6)
                 D2LMUF(IM+7,N) = -2.d0*FF2*PER2(I)*P(IM+7)
                 D2LMUF(N,IM+1) = -3.d0*FF2
                 D2LMUF(N,IM+6) = -2.d0*FF2*PER2(I)*P(IM+6)
                 D2LMUF(N,IM+7) = -2.d0*FF2*PER2(I)*P(IM+7)
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
          LNA = P(IM+1)
          P6 = P(IM+6)
          P7 = P(IM+7)
          LNP = 0.5d0*LOG(PER2(I))
          A = EXP(LNA)
          NN = DPI/SQRT(PER2(I))
          QQ = P(IM+4)
          PP = P(IM+5)
          T2 = QQ*QQ+PP*PP    
          T = SQRT(T2)        ! Tan(i/2)
          CI = (1.d0-T2)/(1.d0+T2)  ! cos(i)
          SI = 2.d0*T/(1.d0+T2)      ! sin(i)           
          CO = QQ/T
          SO = PP/T

c...  Taking into account Period and semi-major axis in FMAP
c...  prior(P) = 1/P (logarithmic prior) => -2*ln(prior(P))=+2*ln(P)          
c...  Idem for semi-major axis
          FMAP = FMAP+2.d0*(LNP+LNA)
          DFMAP(IM+1) = DFMAP(IM+1)+2.d0
          DFMAP(IM+6) = DFMAP(IM+6)-2.d0*P6*PER2(I)
          DFMAP(IM+7) = DFMAP(IM+7)-2.d0*P7*PER2(I)
          D2FMAP(IM+6,IM+6) = D2FMAP(IM+6,IM+6)
     &                    +2.d0*(P6*P6-P7*P7)*PER2(I)*PER2(I)
          D2FMAP(IM+7,IM+7) = D2FMAP(IM+6,IM+6)
     &                    +2.d0*(-P6*P6+P7*P7)*PER2(I)*PER2(I)          
          D2FMAP(IM+6,IM+7) = D2FMAP(IM+6,IM+7)
     &                    +4.d0*P6*P7*PER2(I)*PER2(I)
          D2FMAP(IM+7,IM+6) = D2FMAP(IM+7,IM+6)
     &                    +4.d0*P6*P7*PER2(I)*PER2(I)
c...  Taking into account inc prior in FMAP (-2*ln(prior(i))=-2*ln(sin(i)))
c...        and related derivatives
          FMAP = FMAP-2.d0*LOG(SI)
          DFMAP(IM+4) = DFMAP(IM+4)-2.d0*CO*CI/T
          DFMAP(IM+5) = DFMAP(IM+5)-2.d0*SO*CI/T
          FF1 = CI/T2
          FF2 = 2.d0*CI*CI/T2-4.d0/(T2*(1.d0+T2)*(1.d0+T2))
          D2FMAP(IM+4,IM+4) = -2.d0*(FF2*CO*CO+FF1)
          D2FMAP(IM+4,IM+5) = -2.d0*FF2*CO*SO
          D2FMAP(IM+5,IM+4) = -2.d0*FF2*CO*SO
          D2FMAP(IM+5,IM+5) = -2.d0*(FF2*SO*SO+FF1)
        END DO
        END
      
       
