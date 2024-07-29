C
C        Common utilitary routines
C
C-----------------------------------------------------------------------
C      Upcase conversion
C-----------------------------------------------------------------------
C
        SUBROUTINE UPCASE(CH)

        IMPLICIT NONE

        CHARACTER*(*), INTENT(INOUT) :: CH       ! Chain
        INTEGER       ::        I                ! Index

        DO I=1,LEN(CH)
          IF ((IACHAR(CH(I:I)).GE.IACHAR('a')).AND.
     &        (IACHAR(CH(I:I)).LE.IACHAR('z'))) THEN
            CH(I:I)=ACHAR(IACHAR(CH(I:I))-IACHAR('a')+IACHAR('A'))
          END IF
        END DO
        END SUBROUTINE UPCASE

C     
C-----------------------------------------------------------------
C         Generateur aleatoire Gaussien normalise
C------------------------------------------------------------------
C
        SUBROUTINE GASDEV(S)

        IMPLICIT NONE

        REAL*8          S,R,V(2),S1,S2   ! Parametres internes
        INTEGER*4       J                ! aleatoire 1 ou 2


        CALL RANDOM_NUMBER(S1)
        J=1+FLOOR(2.d0*S1)
        R=2.d0
        DO WHILE(R.GE.1.d0)
          CALL RANDOM_NUMBER(S1)
          CALL RANDOM_NUMBER(S2)
          V(1)=2.d0*S1-1.d0
          V(2)=2.d0*S2-1.d0
          R=V(1)*V(1)+V(2)*V(2)
        END DO
        S = V(J)*SQRT(-2.d0*LOG(R)/R)
        END
C
C -----------------------------------------------------------------------------
C       Drawing a random permutation
C -----------------------------------------------------------------------------
C

        SUBROUTINE PERMUTATION(N,L)

        IMPLICIT NONE

        INTEGER*4 ::    N,           ! Dimension of parameter space  
     &                  K,           ! Index of permutation
     &                  L(N),        ! Permutation
     &                  J,C,KL       ! Indexes
        REAL*8 :: R

        DO J = 1,N
          L(J) = J
        END DO        
        DO J = N,2,-1
          CALL RANDOM_NUMBER(R)          
          KL = FLOOR(J*R)+1
          C = L(J)
          L(J) = L(KL)
          L(KL) = C
        END DO
        END

c        INTEGER*4 ::    N,           ! Dimension of parameter space  
c     &                  K,           ! Index of permutation
c     &                  L(N),        ! Permutation
c     &                  J,KJ,KL,C    ! Indexes
c        KJ = K
c          DO J = 1,N
c          L(J) = J
c        END DO        
c        DO J = 2,N
c          KJ = KJ*(J-1)
c          KL = MOD(KJ,J)+1
c          C = L(J)
c          L(J) = L(KL)
c          L(KL) = C
c        END DO
c        END

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
C       Gaussian elimination
C -----------------------------------------------------------------------------
C

        SUBROUTINE GAUSS(N,A,B,AM,X,TEST)

        IMPLICIT NONE

        INTEGER*4, PARAMETER :: SD = 6
        INTEGER*4 ::    N       ! Dimension of parameter space

        REAL*8 ::       A(N,N),         ! Matrix
     &                  AM(N,N),        ! Inverse matrix
     &                  B(N),           ! Right hand vector
     &                  X(N),           ! Solution vector
     &                  M               ! Running value
        INTEGER*4 ::    I,J,K           ! Indexes
        LOGICAL ::      TEST

        TEST=.TRUE.
        AM = 0.d0
        DO I=1,N
          AM(I,I)=1.d0
        END DO

        DO I=1,N-1
          M=ABS(A(I,I))
          J=I
          DO K=I+1,N
            IF (ABS(A(K,I)).GT.M) THEN
              M=ABS(A(K,I))
              J=K
            END IF
          END DO
          if (m.eq.0.) THEN
            WRITE(SD,*)'Singular matrix !!'
            TEST=.FALSE.
          END IF
          IF (J.NE.I) THEN
            DO K=I,N
              M=A(I,K)
              A(I,K)=A(J,K)
              A(J,K)=M
            END DO
            DO K=1,N
              M=AM(I,K)
              AM(I,K)=AM(J,K)
              AM(J,K)=M
            END DO
            M=B(I)
            B(I)=B(J)
            B(J)=M
          END IF
          DO J=I+1,N
            M=-A(J,I)/A(I,I)
            A(J,I)=0.d0
c            A(J,I+1:N) = A(J,I+1:N)+M*A(I,J+1:N)
c            AM(J,1:N) = AM(J,1:N)+M*AM(I,1:N)
           DO K=I+1,N
              A(J,K)=A(J,K)+M*A(I,K)
            END DO
            DO K=1,N
              AM(J,K)=AM(J,K)+M*AM(I,K)
            END DO
            B(J)=B(J)+M*B(I)
          END DO
        END DO        

        DO I=N,1,-1
          X(I)=B(I)
          DO J=I+1,N
            X(I)=X(I)-A(I,J)*X(J)
c            AM(I,1:N) = AM(I,1:N)-A(I,J)*AM(J,1:N)
            DO K=1,N
              AM(I,K)=AM(I,K)-A(I,J)*AM(J,K)
            END DO
          END DO
          M=1.d0/A(I,I)
          X(I)=X(I)*M
c          AM(I,1:N) = AM(I,1:N)*M
          DO K=1,N
            AM(I,K)=AM(I,K)*M
          END DO            
        END DO

        END
                     

c------------------------------------------------------------
c  To find eigenvalues of a symmetric NxN matrix
c------------------------------------------------------------

      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

      IMPLICIT NONE

      INTEGER :: N,NP,NROT
      INTEGER, PARAMETER :: NMAX = 500
      REAL*8 :: A(NP,NP),D(NP),V(NP,NP)
      INTEGER I,IP,IQ,J
      REAL*8 :: C,G,H,S,SM,T,TAU,THETA,TRESH,B(NMAX),Z(NMAX)
      V = 0.0d0
      DO IP = 1,N
        V(IP,IP) = 1.d0
      END DO
      DO IP = 1,N
        B(IP) = A(IP,IP)
        D(IP) = B(IP)
        Z(IP) = 0.D0
      END DO
      NROT = 0
      DO I = 1,50
        SM = 0.D0
        DO IP = 1,N-1
          DO IQ = IP+1,N
            SM = SM+ABS(A(IP,IQ))
          END DO
        END DO
        IF(SM.EQ.0.D0)RETURN
        IF(I.LT.4)THEN
          TRESH = 0.2D0*SM/N**2
        ELSE
          TRESH = 0.D0
        ENDIF
        DO IP = 1,N-1
          DO IQ = IP+1,N
            G = 100.d0*ABS(A(IP,IQ))
            IF ((I.GT.4).AND.
     &         (ABS(D(IP))+G.EQ.ABS(D(IP))).AND.
     &         (ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ) = 0.D0
            ELSE IF (ABS(A(IP,IQ)).GT.TRESH)THEN
              H = D(IQ)-D(IP)
              IF (ABS(H)+G.EQ.ABS(H))THEN
                T = A(IP,IQ)/H
              ELSE
                THETA = 0.5d0*H/A(IP,IQ)
                T = 1.d0/(ABS(THETA)+SQRT(1.d0+THETA**2))
                IF(THETA.LT.0.d0) T = -T
              ENDIF
              C = 1.d0/SQRT(1.d0+T**2)
              S = T*C
              TAU = S/(1.D0+C)
              H = T*A(IP,IQ)
              Z(IP) = Z(IP)-H
              Z(IQ) = Z(IQ)+H
              D(IP) = D(IP)-H
              D(IQ) = D(IQ)+H
              A(IP,IQ) = 0.D0
              DO J = 1,IP-1
                G = A(J,IP)
                H = A(J,IQ)
                A(J,IP) = G-S*(H+G*TAU)
                A(J,IQ) = H+S*(G-H*TAU)
              END DO
              DO J = IP+1,IQ-1
                G = A(IP,J)
                H = A(J,IQ)
                A(IP,J) = G-S*(H+G*TAU)
                A(J,IQ) = H+S*(G-H*TAU)
              END DO
              DO J = IQ+1,N
                G = A(IP,J)
                H = A(IQ,J)
                A(IP,J) = G-S*(H+G*TAU)
                A(IQ,J) = H+S*(G-H*TAU)
              END DO
              DO J = 1,N
                G = V(J,IP)
                H = V(J,IQ)
                V(J,IP) = G-S*(H+G*TAU)
                V(J,IQ) = H+S*(G-H*TAU)
              END DO
              NROT = NROT+1
            ENDIF
          END DO
        END DO
        DO IP = 1,N
          B(IP) = B(IP)+Z(IP)
          D(IP) = B(IP)
          Z(IP) = 0.D0
        END DO
      END DO
      STOP 'too many iterations in Jacobi'
      RETURN
      END
 
C
C -----------------------------------------------------------------------------
C    Returns the error of a combined variable
C -----------------------------------------------------------------------------
C

        FUNCTION ERR_COMB(N,DFDV,COV)

        IMPLICIT NONE     

        INTEGER*4       N
	REAL*8		DFDV(N),COV(N,N),VF,
     &                  ERR_COMB
        INTEGER*4       J,I
        
        VF = 0.d0
        DO J = 1,N
          DO I = 1,N
            VF = VF+COV(I,J)*DFDV(I)*DFDV(J)
          END DO
        END DO
        ERR_COMB = SQRT(ABS(VF))

        END


C
C-----------------------------------------------------------------------------
C    Converts a date to Julian Day
C     Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
C        4th ed., Duffet-Smith and Zwart, 2011.
C-----------------------------------------------------------------------------
C

        REAL*8 FUNCTION DATE_TO_JD(JJ,MM,YY)

        IMPLICIT NONE

        INTEGER*4 ::   MM,YY    ! Month, Year
        INTEGER*4 ::   MMP,YYP
        REAL*8 ::      JJ       ! Day
        INTEGER*4 :: A,B,C,D
c        REAL*8 ::    DATE_TO_JD

        IF ((MM.EQ.1).OR.(MM.EQ.2)) THEN
           YYP = YY-1
           MMP = MM+12
        ELSE
           YYP = YY
           MMP = MM
        END IF        

c... This checks where we are in relation to October 15, 1582, the beginning
c... of the Gregorian calendar.
        IF ((YY.LT.1582).OR.
     &       ((YY.EQ.1582).AND.(MM.LT.10)).OR.
     &       ((YY.EQ.1582).AND.(MM.EQ.10).AND.(JJ.LT.15))) THEN
           B = 0
c...         before start of Gregorian calendar
        ELSE
c...        after start of Gregorian calendar
           A = FLOOR(0.01d0*DBLE(YYP))
           B = 2 - A + FLOOR(0.25d0*DBLE(A))
        END IF        

        IF (YYP.LT.0) THEN
          C = FLOOR(365.25d0*DBLE(YYP)-0.75d0)
        ELSE
          C = FLOOR(365.25d0*DBLE(YYP))
        END IF        
        D = FLOOR(30.6001d0*DBLE(MMP+1))
        DATE_TO_JD = DBLE(B+C+D+JJ)+1720994.5d0

        END
C
C-----------------------------------------------------------------------------
C    Converts a Julian Day to a date
C     Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
C        4th ed., Duffet-Smith and Zwart, 2011.
C-----------------------------------------------------------------------------
C

        SUBROUTINE JD_TO_DATE(JD,JJ,MM,YY)

        IMPLICIT NONE

        INTEGER*4 ::   MM,YY    ! Month, Year
        REAL*8 ::      JJ,JD    ! Day, date
        INTEGER*4 :: A,B,C,D,I,E,F,G

        I = INT(JD)
        F = JD-I
        A = FLOOR((DBLE(I)-1867216.25d0)/36524.25d0)

        IF (I.GT.2299160) THEN
           B = I+1+A-FLOOR(0.25d0*DBLE(A))
        ELSE
           B = I
        END IF
        C = B+1524
        D = FLOOR((DBLE(C)-122.1d0)/365.25d0)
        E = FLOOR(365.25d0*D)
        G = FLOOR(DBLE(C-E)/30.6001d0)
        JJ = C-E+F-FLOOR(30.6001d0*DBLE(G))

        IF (G.LT.13.5d0) THEN
           MM = G-1
        ELSE
           MM = G-13
        END IF
        IF (MM.GT.2.5d0) THEN
           YY = D-4716
        ELSE
           YY = D-4715
        END IF

        END
     

C-----------------------------------------------------------------------------
C    Entering inital guess for orbits
C-----------------------------------------------------------------------------
C

        SUBROUTINE ENTER_ORBITS(N)

        USE DATA
        
        IMPLICIT NONE

        INTEGER*4 ::    I,N     ! # of planet
        REAL*8 ::       SIGMA,  ! Cumulative mass
     &                  NN      ! Mean motion                 
        CHARACTER*2 ::  UNIT

c        FORMAT(f6.1,a2,(f))
        WRITE(SD,'(a)')
     & 'Enter initial orbital guess for each planet under the form'
        WRITE(SD,'(a)')' mass unit: ms ou mj (solar or jupiter mass)'
        WRITE(SD,'(a)')' m unit a e i Omega omega tp, where'
        WRITE(SD,'(a)')'   m = Mass'
        WRITE(SD,'(a)')'   unit = ms ou mj (Solar or Jupiter mass)'
        WRITE(SD,'(a)')'   a = Semi major axis (au)'
        WRITE(SD,'(a)')'   e = Eccentricity'
        WRITE(SD,'(a)')'   i = Inclination (deg)'
        WRITE(SD,'(a)')
     &    '   Omega = Long. of periastron (deg, relative to North)' 
        WRITE(SD,'(a)')'   omega = Argument of periastron (deg)'
        WRITE(SD,'(a)')'   tp = Time for periastron (JD-offset)'
        SIGMA = STAR%MASS
        DO I = 1,N
           WRITE(SD,*)'Give input for planet #',I
           READ(5,*)PLA(I)%MU,UNIT,PLA(I)%A,PLA(I)%EXC,PLA(I)%INC,
     &                        PLA(I)%O,PLA(I)%OM,PLA(I)%TP
           CALL UPCASE(UNIT)
           IF (UNIT.EQ.'MS') PLA(I)%MU = PLA(I)%MU*SMAS
           IF (UNIT.EQ.'MJ') PLA(I)%MU = PLA(I)%MU*MJUP
           SIGMA = SIGMA+PLA(I)%MU
           PLA(I)%MDYN = SIGMA
           NN = SQRT(SIGMA/PLA(I)%A**3)
           PLA(I)%PER = DPI/NN
           PLA(I)%EXQ = SQRT(1.d0-PLA(I)%EXC**2)
           PLA(I)%TI = TAN(0.5d0*PLA(I)%INC/DR)
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
        END DO  
        END            

C
C-----------------------------------------------------------------------------
C    Reading user given prior parameters (for masses)
C-----------------------------------------------------------------------------
C

        SUBROUTINE ENTER_PRIOR(IV)

        USE DATA
        
        IMPLICIT NONE

        INTEGER*4 :: IV,I         ! # of variable
        CHARACTER*2 ::   UNIT     ! Unity (Msun,Mjup...)
        REAL*8 :: FAC

c        FORMAT((f),a2,(f))

        ALLOCATE(MPRIOR(IV)%ACOF(0:NPLA))
        ALLOCATE(MPRIOR(IV)%BCOF(0:NPLA))      
        MPRIOR(IV)%MEAN = 0.d0
        MPRIOR(IV)%SDEV = 0.d0
        MPRIOR(IV)%BOUND = 0.d0
        WRITE(SD,*)'Prior on masses : Enter list of coefficients'
        WRITE(SD,*)'        a_0*m0+a_1*m1+...+a_n*mn'       
        READ(5,*)MPRIOR(IV)%ACOF(0:NPLA)
        DO I = 0,NPLA-1
           MPRIOR(IV)%BCOF(I) = MPRIOR(IV)%ACOF(I)-MPRIOR(IV)%ACOF(I+1)
        END DO
        MPRIOR(IV)%BCOF(NPLA) = MPRIOR(IV)%ACOF(NPLA)
        WRITE(SD,*)' Enter type of prior :'
        WRITE(SD,*)'   0 = Logarithmic (default)'
        WRITE(SD,*)'   1 = Gaussian'
        WRITE(SD,*)'   2 = Log-Normal'
        WRITE(SD,*)'   3 = Linear'
        WRITE(SD,*)'   4 = Fixed'
        READ(5,*)MPRIOR(IV)%TYP
        SELECT CASE (MPRIOR(IV)%TYP)
        CASE(1) ! Normal
           WRITE(SD,*)
     &           'Enter mean and standard deviation + unit (ms or mj) '
          READ(5,*)MPRIOR(IV)%MEAN,MPRIOR(IV)%SDEV,UNIT
          CALL CONVERT_UNIT(UNIT,FAC)
          MPRIOR(IV)%MEAN = MPRIOR(IV)%MEAN*FAC ! Convert into the right unit
          MPRIOR(IV)%SDEV = MPRIOR(IV)%MEAN*FAC
        CASE(2)                  ! Log-normal
          WRITE(SD,*)
     &           'Enter mean and standard deviation + unit (ms or mj) '
          READ(5,'(a)')MPRIOR(IV)%MEAN,MPRIOR(IV)%SDEV,UNIT
          CALL CONVERT_UNIT(UNIT,FAC)
          MPRIOR(IV)%SDEV = MPRIOR(IV)%SDEV/MPRIOR(IV)%MEAN ! ln(1±s/x)~s/x
          MPRIOR(IV)%MEAN = LOG(MPRIOR(IV)%MEAN*FAC)
        CASE(3) ! Linear
          WRITE(SD,*)'Enter bounds + unit (ms or mj) '
          READ(5,'(a)')MPRIOR(IV)%BOUND(1:2),UNIT
          CALL CONVERT_UNIT(UNIT,FAC)
          MPRIOR(IV)%BOUND(1:2) = MPRIOR(IV)%BOUND(1:2)*FAC
        CASE(4) ! Fixed
          WRITE(SD,*)'Enter value + unit (ms or mj) '
          READ(5,'(a)')MPRIOR(IV)%MEAN,UNIT
          CALL CONVERT_UNIT(UNIT,FAC)
          MPRIOR(IV)%MEAN = MPRIOR(IV)%MEAN*FAC ! Convert into the right unit
          MPRIOR(IV)%SDEV = 0.d0
        END SELECT

c...    Internal routine : treatment of the mass unit character
        
        CONTAINS

          SUBROUTINE CONVERT_UNIT(UNT,FC)

          IMPLICIT NONE

          CHARACTER*2 :: UNT
          REAL*8 :: FC
          
          CALL UPCASE(UNIT)
          IF (UNT.EQ.'MS') FC = SMAS
          IF (UNT.EQ.'MJ') FC = MJUP        

          END SUBROUTINE CONVERT_UNIT

        END
      
C
C -----------------------------------------------------------------------------
C     Routine to initiate the simulation
C -----------------------------------------------------------------------------
C

        SUBROUTINE INIT_SIMU(CNEW,NMOD,FILES,NEW,CORR)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    I,J,CNEW,ERROR,NMOD
        INTEGER*4, DIMENSION(4) :: ISDATAN
        LOGICAL ::      OK,NEW,CORR
        CHARACTER*3 ::  UNIT
        CHARACTER*80, DIMENSION(4) :: FILES ! #1 = output, #2 = dump,
c                                                #3--NFIL = data
        
        WRITE(SD0,*)'Choice :'
        WRITE(SD0,*)'  1 = Continuation of an interruped computation'
        WRITE(SD0,*)'  2 = New MCMC with starting least squares'
        WRITE(SD0,*)'  3 = New MCMC without starting least squares'
        WRITE(SD0,*)'  4 = Least square only with no MCMC'
        WRITE(SD0,*)'  5 = Generate simulated data'
        READ(5,*)CNEW
        NEW = (CNEW.GT.1)

        DATATYP = 0
        FILES = ' '
        IF (NEW) THEN
           OPEN(SDL,FILE=LOGFILE,STATUS='UNKNOWN')
        ELSE
           OPEN(SDL,FILE=LOGFILE,STATUS='OLD',POSITION='APPEND')
        END IF
        
        IF (NEW) THEN
           WRITE(SD,*)
     &'Give types of data present (1=present, 0=absent, ex: 1 1 0 0):'
           WRITE(SD,*)'  Number #1 => relative astrometric data'
           WRITE(SD,*)'  Number #2 => stellar radial velocity data'
           WRITE(SD,*)'  Number #3 => relative radial velocity data'           
           WRITE(SD,*)'  Number #4 => absolute stellar astrometric data'
           READ(5,*)ISDATAN(1:NDATATYP)
           DO I = 1,4
              ISDATA(I) = (ISDATAN(I).EQ.1)
           END DO
           RADVEL = (ISDATA(2).OR.ISDATA(3))
           WRITE(SD,*)'Number of planets :'
           READ(5,*)NPLA
           MULTIPLA = (NPLA.GT.1)
           MULTIPLA = (MULTIPLA.OR.RADVEL)
c...  MULTIPLA = .TRUE. whenever individual masses will have to be fitted.
c...  This is the case if more than one orbit is considered (NPLA>1)
c...  or when radial velocity data are present.           
           ALLOCATE(PLA(NPLA))
           WRITE(SD,*)'Name of astrometric data file ?'
           READ(5,'(a)')FILES(3)
           IF (RADVEL) THEN 
              WRITE(SD,*)'Name of radial velocity data file ?'
              READ(5,'(a)')FILES(4)
           END IF
           WRITE(SD,*)'Precision ?'
           READ(5,*)EPS
           WRITE(SD,*)'Choose data line format :'
           WRITE(SD,*)
     & ' (Dec,RA,Sep in mas, PA in deg, Corr=correlation coefficient)'
           WRITE(SD,*)'   1 = ipl Day Month Year Dec RA dDec dRA Corr'
           WRITE(SD,*)'   2 = ipl Day Month Year Dec RA dDec dRA'           
           WRITE(SD,*)'   3 = ipl Jul.day Dec RA dDec dRA Corr'
           WRITE(SD,*)'   4 = ipl Jul.day Dec RA dDec dRA'
           WRITE(SD,*)'   5 = ipl Day Month Year Sep PA dSep dPA Corr'
           WRITE(SD,*)'   6 = ipl Day Month Year Sep PA dSep dPA'           
           WRITE(SD,*)'   7 = ipl Jul.day Sep PA dSep dPA Corr'
           WRITE(SD,*)'   8 = ipl Jul.day Sep PA dSep dPA'
           READ(5,*)DATATYP(1)
           CORR = (MOD(DATATYP(1),2).NE.0)
           IF (RADVEL) THEN              
              WRITE(SD,*)'Choose RV data line format :'
              WRITE(SD,*)' (RV in km/s)'
              WRITE(SD,*)'   1 = Day Month Year RV dRV'
              WRITE(SD,*)'   2 = Jul.day RV dRV'
              READ(5,*)DATATYP(2)
              WRITE(SD,*)'Should we consider a jitter (1=yes, 0=no) ?'
              READ(5,*)JITNUM
           END IF   
           WRITE(SD,*)
     &      'Distance of system (pc) or parallax (mas) (specify unit) ?'
           READ(5,*)STAR%DIST,UNIT
           CALL UPCASE(UNIT)
           IF (TRIM(UNIT).EQ.'PC') STAR%PARX = 1.d3/STAR%DIST
           IF (UNIT.EQ.'MAS') THEN
              STAR%PARX = STAR%DIST
              STAR%DIST = 1.d3/STAR%PARX
           END IF
           WRITE(SD,*)'Mass of central body (Msun) (Initial guess) ?'
           READ(5,*)STAR%MASS
           STAR%MASS = STAR%MASS*SMAS   ! Convert mass into the right unit
           STAR%DMASS = 0.d0
           IF (RADVEL) STAR%V0 = 0.d0
           IF (CNEW.LT.5) THEN
              WRITE(SD,*)'Name of output file ?'
              READ(5,'(a)')FILES(1)
              FILES(1) = TRIM(FILES(1))//'.dat'
              WRITE(SD,*)'Name of dump file ?'
              READ(5,'(a)')FILES(2)
              WRITE(SD,*)'Dump frequency ?'
              READ(5,*)FREQDUMP 
              IF (MULTIPLA) THEN
                 NPAR = NEL*NPLA+1  ! Store central mass into the last variable
              ELSE
                 NPAR = NEL*NPLA
              END IF
              IF (RADVEL) NPAR = NPAR+1+JITNUM
              NMOD = 0
c...  If MULTIPLA=F then NPAR = NEL*NPLA (variables = params * nb of orbits)
c...  If MULTIPLA=T but RADVEL=F then NPAR=NEL*NPAR+1
c...       Variable #NPAR = ln(central mass)
c...  If RADVEL=T and no jitter NPAR = NEL*NPAR+2
c...       Variables #NPAR-1 = V0, #NPAR = ln(central mass)
c...  If RADVEL=T and jitter  NPAR = NEL*NPAR+3
c...       Variables #NPAR-2 = V0, #NPAR-1=ln(Jitter), #NPAR = ln(central mass)
              ALLOCATE(PSTART(NPAR+2))
              ALLOCATE(PSAV(NPAR+2,NMAX))
c...  NPAR+2 => #NPAR+1 = CHI2, #NPAR+2 = FMAP in output
              WRITE(SD,*)'Number of mass priors'
              READ(5,*)NPRIOR
              IF (MULTIPLA) THEN
                 NPRIOR = NPRIOR+NPLA
              ELSE
                 NPRIOR = NPRIOR-1
              END IF
c...  Numbering of mass priors starts at #0. 
              IF (NPRIOR.GE.0) THEN
                 ALLOCATE(MPRIOR(0:NPRIOR))
                 IF (MULTIPLA) THEN
                    J = NPLA+1
                 ELSE
                    J = 0
                 END IF
                 DO I = J,NPRIOR            
                    CALL ENTER_PRIOR(I)
                 END DO
              END IF              
c...  Priors #0..NPLA tell that individual masses must be positive
              IF (MULTIPLA) THEN
                 DO I = 0,NPLA
                    ALLOCATE(MPRIOR(I)%ACOF(0:NPLA))
                    ALLOCATE(MPRIOR(I)%BCOF(0:NPLA))      
                    MPRIOR(I)%TYP = 6
                    MPRIOR(I)%MEAN = 0.d0
                    MPRIOR(I)%SDEV = MTINY
                    MPRIOR(I)%BOUND = 0.d0
                    MPRIOR(I)%ACOF = 0.d0
                    MPRIOR(I)%BCOF = 0.d0
                    MPRIOR(I)%ACOF(I) = 1.d0
                    IF (I.GT.0) MPRIOR(I)%BCOF(I-1) = -1.d0
                    MPRIOR(I)%BCOF(I) = 1.d0
                 END DO   
              END IF
           END IF
        ELSE
           WRITE(SD,*)'Name of dump file to read ?'
           READ(5,'(a)')FILES(2)
           WRITE(SD,*)'Dump frequency ?'
           READ(5,*)FREQDUMP  
           CALL READ_DUMP0(NMOD,FILES)
           CORR = (MOD(DATATYP(1),2).NE.0)
        END IF

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
        LOGICAL ::      OK,NEW,CORR
        CHARACTER*80, DIMENSION(4) :: FILES ! #1 = output, #2 = dump,
c                                                #3--NFIL = data

        CALL READ_ASTROM(1,FILES(3))
        NFREE = 2*SUM(PLA%NDATAS)-NPAR ! # of degrees od freedom

        IF (NEW) THEN
           WRITE(SD,*)'Give reference time for data (JD) : '
           READ(5,*)STAR%T0
           STAR%V0 = 0.d0
           STAR%SIGJV = 0.d0
           IF (RADVEL) THEN
              CALL READ_RV(2,FILES(4))
c...  Initialize offset velocity & jitter
              WRITE(SD,*)
     &              'Give initial offset velocity and jitter in km/s '
              READ(5,*)STAR%V0,STAR%SIGJV
              STAR%SIGJV = STAR%SIGJV*MPS
              STAR%V0 = STAR%V0*MPS
           END IF             
        END IF
        IF (RADVEL) NFREE = 2*SUM(PLA%NDATAS)+STAR%NDATVR-NPAR
        
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
           IF (RADVEL) P(NPAR-1) = STAR%V0
           IF (RADVEL.AND.(JITNUM.EQ.1)) THEN
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
C     Routine to enter initial data
C -----------------------------------------------------------------------------
C

        SUBROUTINE SIMUL_DATA(IDAT,FILES)

        USE DATA

        IMPLICIT NONE

        INTEGER*4, PARAMETER :: NMW = 100
        INTEGER*4, DIMENSION(NMW) :: IPLW ! Read planets numbers
        REAL*8, DIMENSION(NMW) :: JDW     ! Read dates
        REAL*8, DIMENSION(2,NPLA) ::
     &                  E1,E2                ! Vecteurs de la base propre
        REAL*8 ::      DATE_TO_JD            ! Conversion function
        REAL*8, DIMENSION(NPLA) :: CI,       ! Cos(i)
     &                             MFRAC,    ! Fractional masses
     &                             NN,       ! Mean motions
     &                             POSX,POSY ! Heliocentic positions
        INTEGER*4, DIMENSION(NPLA) :: IPL    ! Curent indexes
        INTEGER*4 :: J,I,ERROR,MM,YR,IPLA,NDAT
        INTEGER*4 :: IDAT
        REAL*8 :: LAMBDA,CL,SL,DERR,JJ,JD,S 
        LOGICAL ::      OK,CORR
        CHARACTER*80, DIMENSION(4) :: FILES ! #1 = output, #2 = dump,
c                                                #3--NFIL = data

        CALL ENTER_ORBITS(NPLA)
        CI(:) = (1.d0-PLA(:)%TI*PLA(:)%TI)/(1.d0+PLA(:)%TI*PLA(:)%TI)
        E1(1,:) = PLA(:)%COM*PLA(:)%CO-PLA(:)%SOM*CI(:)*PLA(:)%SO
        E1(2,:) = PLA(:)%COM*PLA(:)%SO+PLA(:)%SOM*CI(:)*PLA(:)%CO
        E2(1,:) = -PLA(:)%SOM*PLA(:)%CO-PLA(:)%COM*CI(:)*PLA(:)%SO
        E2(2,:) = -PLA(:)%SOM*PLA(:)%SO+PLA(:)%COM*PLA(:)%CO*CI(:)
        MFRAC(:) = PLA(:)%MU/PLA(:)%MDYN
        NN(:) = DPI/PLA(:)%PER
        
        OK = .TRUE.
        PLA%NDATAS = 0
        SELECT CASE(DATATYP(IDAT))
        CASE(3,4,7,8)
           WRITE(SD,*)'Give JD offset'
           READ(5,*)STAR%JDOFFSET        
        END SELECT
        WRITE(SD,*)'Give error in mas'
        READ(5,*)DERR
        NDAT = 0
        CORR = (MOD(DATATYP(IDAT),2).NE.0)
        DO WHILE(OK)
           SELECT CASE (DATATYP(IDAT))
           CASE(1,2,5,6)
              READ(5,*,IOSTAT=ERROR)IPLA,JJ,MM,YR
           CASE(3,4,7,8)
              READ(5,*,IOSTAT=ERROR)IPLA,JD
           END SELECT
           OK=(ERROR.EQ.0)
           IF (OK) THEN   
              NDAT = NDAT+1
              PLA(IPLA)%NDATAS = PLA(IPLA)%NDATAS+1
              IPLW(NDAT) = IPLA
              SELECT CASE (DATATYP(IDAT))
              CASE(1,2,5,6)
c... Conversion JJ-MM-YR -> JD
                 JD = DATE_TO_JD(JJ,MM,YR)-STAR%JDOFFSET
              CASE(3,4,7,8)
                 JD = JD   ! +STAR%JDOFFSET
              END SELECT  
              JDW(NDAT) = JD
           END IF
        END DO
        DO I = 1,NPLA
           J = PLA(I)%NDATAS
           ALLOCATE(PLA(I)%TAS(J))
           ALLOCATE(PLA(I)%X(J))
           ALLOCATE(PLA(I)%SIGX(J))
           ALLOCATE(PLA(I)%Y(J))
           ALLOCATE(PLA(I)%SIGY(J))
           ALLOCATE(PLA(I)%RHOXY(J))
           ALLOCATE(PLA(I)%SIGXM1(J))
           ALLOCATE(PLA(I)%SIGYM1(J))
           ALLOCATE(PLA(I)%SIGXM2(J))
           ALLOCATE(PLA(I)%SIGYM2(J))
           ALLOCATE(PLA(I)%UFAC(J))
        END DO
        IPL = 0
        DO J = 1,NDAT
           IPLA = IPLW(J)
           IPL(IPLA) = IPL(IPLA)+1
           I = IPL(IPLA)
           PLA(IPLA)%TAS(I) = JDW(J) ! Time stored in Jul days-JDoffset
           CALL POSFITS(JDW(J),NN,PLA(:)%A,PLA(:)%EXC,PLA(:)%EXQ,
     &                E1,E2,PLA(:)%TP,MFRAC(:),POSX,POSY)
           CALL GASDEV(S)
c           print*,s
           PLA(IPLA)%X(I) = POSX(IPLA)+S*DERR/STAR%PARX              
c           print*,s,pla(ipla)%x(i),posx(ipla)
           CALL GASDEV(S)
           PLA(IPLA)%Y(I) = POSY(IPLA)+S*DERR/STAR%PARX
           PLA(IPLA)%SIGX(I) = DERR/STAR%PARX
           PLA(IPLA)%SIGY(I) = DERR/STAR%PARX
           PLA(IPLA)%SIGXM1(I) = 1.d0/PLA(IPLA)%SIGX(I)
           PLA(IPLA)%SIGYM1(I) = 1.d0/PLA(IPLA)%SIGY(I)
           PLA(IPLA)%SIGXM2(I) = PLA(IPLA)%SIGXM1(I)**2
           PLA(IPLA)%SIGYM2(I) = PLA(IPLA)%SIGYM1(I)**2
           IF (CORR) THEN
              CALL RANDOM_NUMBER(S)
              PLA(IPLA)%RHOXY(I) = 2.d0*S-1.d0
           ELSE
              PLA(IPLA)%RHOXY(I) = 0.d0
           END IF
           PLA(IPLA)%UFAC(I) =
     &          1.d0/(1.d0-PLA(IPLA)%RHOXY(I)**2)
        END DO
        CALL WRITE_ASTROM(IDAT,FILES(3))

        END

C     
C -----------------------------------------------------------------------------
C       Writing a dump file
C -----------------------------------------------------------------------------
C
        SUBROUTINE WRITE_DUMP(NMOD,COUNTOK,LONG,NGEL,INFER,
     &           PROBINFER,CHN,DIRS,BETA0,FILES)

        USE DATA

        IMPLICIT NONE
     
        CHARACTER*(*), DIMENSION(NFIL) :: FILES
        TYPE(CHAIN), DIMENSION(NCH) :: CHN ! The chains
        INTEGER*4 ::    C,I,    ! Indexes
     &                  COUNTOK,     ! Counter
     &                  ERROR,
     &                  NMOD               ! Number of models   
        INTEGER*8 ::    LONG,
     &                  NGEL         ! Time for Gelman-Rubin
        LOGICAL ::      INFER        ! Inference flag
        REAl*8 ::       PROBINFER,   ! Probability of retaining a model
     &                  DIRS(NPAR,NPAR), ! Directions
     &                  BETA0(NPAR)  ! Scale parameter vector
        
        OPEN(18,FILE=FILES(2),STATUS='UNKNOWN')
        WRITE(18,*,IOSTAT=ERROR)ISDATA(1:NDATATYP)
        WRITE(18,*,IOSTAT=ERROR)NPAR,NPLA,NPRIOR,JITNUM
        DO I=1,NLIM        
           WRITE(18,*,IOSTAT=ERROR)VPRIOR(I)%BOUND(1:2)
        END DO
        IF (NPRIOR.GE.0) THEN
           DO I=0,NPRIOR
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%TYP
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%ACOF(0:NPLA)
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%BCOF(0:NPLA)           
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%MEAN,MPRIOR(I)%SDEV
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%BOUND
           END DO
        END IF
        WRITE(18,*,IOSTAT=ERROR)EPS,STAR%DIST,STAR%PARX
        WRITE(18,*,IOSTAT=ERROR)STAR%MASS,STAR%T0,STAR%V0,STAR%SIGJV
        WRITE(18,'(a)',IOSTAT=ERROR)FILES(1)
        DO I = 2,NFIL
          WRITE(18,'(a)',IOSTAT=ERROR)FILES(I)
        END DO
        WRITE(18,*,IOSTAT=ERROR)DATATYP(1:NDATATYP)
        WRITE(18,*,IOSTAT=ERROR)PSTART(1:(NPAR+2))
        WRITE(18,*,IOSTAT=ERROR)NMOD
        WRITE(18,*,IOSTAT=ERROR)DIRS(1:NPAR,1:NPAR)
        WRITE(18,*,IOSTAT=ERROR)BETA0(1:NPAR)
        DO C = 1,NCH
           WRITE(18,*,IOSTAT=ERROR)CHN(C)%MOY(1:NPAR)
           WRITE(18,*,IOSTAT=ERROR)CHN(C)%BETA(1:NPAR)
           WRITE(18,*,IOSTAT=ERROR)CHN(C)%P(1:NPAR)           
           WRITE(18,*,IOSTAT=ERROR)CHN(C)%COV(1:NPAR,1:NPAR)
           WRITE(18,*,IOSTAT=ERROR)CHN(C)%CHI2
        END DO
        DO I = 1,NMAX
           WRITE(18,*,IOSTAT=ERROR)PSAV(1:NPAR+2,I)
        END DO
        WRITE(18,*,IOSTAT=ERROR)LONG,NGEL
        WRITE(18,*,IOSTAT=ERROR)INFER
        WRITE(18,*,IOSTAT=ERROR)PROBINFER
        WRITE(18,*,IOSTAT=ERROR)COUNTOK
        CLOSE(18)  

        END
C     
C -----------------------------------------------------------------------------
C       Reading a dump file (begin)
C -----------------------------------------------------------------------------
C
        SUBROUTINE READ_DUMP0(NMOD,FILES)

        USE DATA

        IMPLICIT NONE
     
        CHARACTER*(*), DIMENSION(NFIL) :: FILES
        INTEGER*4 ::    ERROR,I,NP1,
     &                  NMOD             ! Number of models   

        OPEN(18,FILE=FILES(2),STATUS='UNKNOWN')
        READ(18,*,IOSTAT=ERROR)ISDATA(1:NDATATYP)
        RADVEL = (ISDATA(2).OR.ISDATA(3))
        READ(18,*,IOSTAT=ERROR)NPAR,NPLA,NPRIOR,JITNUM
        MULTIPLA = (NPLA.GT.1)
        MULTIPLA = (MULTIPLA.OR.RADVEL)

        DO I = 1,NLIM
           READ(18,*,IOSTAT=ERROR)VPRIOR(I)%BOUND(1:2)
           VPRIOR(I)%MEAN = 0.0d0
           VPRIOR(I)%SDEV = 0.d0
        END DO
        VPRIOR(1:2)%TYP = 0
        VPRIOR(3)%TYP = 3
        ALLOCATE(PSTART(NPAR+2))
        ALLOCATE(PLA(NPLA))
        ALLOCATE(PSAV(NPAR+2,NMAX))
        IF (NPRIOR.GE.0) THEN
           ALLOCATE(MPRIOR(0:NPRIOR))
           DO I = 0,NPRIOR
              ALLOCATE(MPRIOR(I)%ACOF(0:NPLA))
              ALLOCATE(MPRIOR(I)%BCOF(0:NPLA))     
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%TYP
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%ACOF(0:NPLA)
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%BCOF(0:NPLA)           
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%MEAN,MPRIOR(I)%SDEV
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%BOUND(1:2)
           END DO
        END IF
        READ(18,*,IOSTAT=ERROR)EPS,STAR%DIST,STAR%PARX
        READ(18,*,IOSTAT=ERROR)STAR%MASS,STAR%T0,STAR%V0,STAR%SIGJV
        READ(18,'(a)',IOSTAT=ERROR)FILES(1)
        DO I = 2,NFIL
          READ(18,'(a)',IOSTAT=ERROR)FILES(I)
        END DO
        READ(18,*,IOSTAT=ERROR)DATATYP(1:NDATATYP) 
        READ(18,*)PSTART(1:(NPAR+2))
        READ(18,*)NMOD
        END
        
C     
C -----------------------------------------------------------------------------
C       Reading a dump file (end)
C -----------------------------------------------------------------------------
C
        SUBROUTINE READ_DUMP(NMOD,COUNTOK,LONG,NGEL,INFER,
     &     PROBINFER,CHN,DIRS,BETA0)

        USE DATA

        IMPLICIT NONE
        TYPE(CHAIN), DIMENSION(NCH) :: CHN ! The chains     
        INTEGER*4 ::    C,I,         ! Indexes
     &                  COUNTOK,     ! Counter
     &                  NMOD               ! Number of models   
        INTEGER*8 ::    LONG,
     &                  NGEL         ! Time for Gelman-Rubin
        LOGICAL ::      INFER        ! Inference flag
        REAl*8 ::       PROBINFER,   ! Probability of retaining a model
     &                  DIRS(NPAR,NPAR),  ! Directions
     &                  BETA0(NPAR)   ! Scale parameter vector
        
        READ(18,*)DIRS(1:NPAR,1:NPAR)
        READ(18,*)BETA0(1:NPAR)
        DO C = 1,NCH
           READ(18,*)CHN(C)%MOY(1:NPAR)
           READ(18,*)CHN(C)%BETA(1:NPAR)
           READ(18,*)CHN(C)%P(1:NPAR)           
           READ(18,*)CHN(C)%COV(1:NPAR,1:NPAR)
           READ(18,*)CHN(C)%CHI2
        END DO
        DO I = 1,NMAX
           READ(18,*)PSAV(1:NPAR+2,I)   ! Sauvetage inversé
        END DO
        READ(18,*)LONG,NGEL
        READ(18,*)INFER
        READ(18,*)PROBINFER
        READ(18,*)COUNTOK
        CLOSE(18)        

        END
        
C
C -----------------------------------------------------------------------------
C       Ecrire un ficher resultat (doublement des solutions (O,w))
C -----------------------------------------------------------------------------
C
        SUBROUTINE WRITE_DISTRIB_DOUBLING(NMOD,DEV)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    NMOD,         ! Number of models 
     &                  NSAV          ! Number of pars. to store per planet
        CHARACTER*(*) :: DEV          ! 
        REAL*8, DIMENSION(:), ALLOCATABLE ::
     &                  NN,           ! Moyen mouvement
     &                  CI2,SI2       ! cos^2(i/2), sin^2(i/2)
        REAL*8 ::       CHI2,         ! Chi2
     &                  FMAP,         ! MAP
     &                  SIGMA         ! Cumulative mass
        REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: DATA4
        INTEGER*4       I,K

        NSAV = NEL+4
        ALLOCATE(DATA4(2*(NMOD+1),NSAV,NPLA))
        ALLOCATE(NN(NPLA))
        ALLOCATE(CI2(NPLA))
        ALLOCATE(SI2(NPLA))
        CALL ELEMENTS(PSTART(1:NPAR),NN,PLA%A,PLA%EXC,PLA%EXQ,
     &       PLA%CW,PLA%SW,CI2,SI2,PLA%CP,PLA%SP,PLA%TP,PLA%MU)
        IF (MULTIPLA) STAR%MASS = EXP(PSTART(NPAR))
        
        PLA%PER = DPI/NN
        PLA%W = ATAN2(PLA%SW,PLA%CW)
        PLA%PHI = ATAN2(PLA%SP,PLA%CP)
        PLA%OM = 0.5d0*(PLA%W+PLA%PHI)
        PLA%O = 0.5d0*(PLA%W-PLA%PHI)
        PLA%INC = 2.d0*ATAN2(SQRT(SI2),SQRT(CI2))
        CHI2 = PSTART(NPAR+1)
        FMAP = PSTART(NPAR+2)
        DO K = 1,NPLA
           DATA4(1,1:NSAV,K) = (/ SNGL(PLA(K)%A),SNGL(PLA(K)%PER),
     &       SNGL(PLA(K)%EXC),SNGL(PLA(K)%OM),SNGL(PLA(K)%INC),
     &       SNGL(PLA(K)%O),SNGL(PLA(K)%TP),SNGL(PLA(K)%MU/MJUP),
     &       SNGL(STAR%MASS/SMAS),SNGL(CHI2),SNGL(FMAP) /)
        END DO
        PLA%OM = MOD(PLA%OM+PI+PI,DPI)-PI
        PLA%O = MOD(PLA%O+PI+PI,DPI)-PI
        DO K = 1,NPLA
           DATA4(2,1:NSAV,K) = (/ SNGL(PLA(K)%A),SNGL(PLA(K)%PER),
     &       SNGL(PLA(K)%EXC),SNGL(PLA(K)%OM),SNGL(PLA(K)%INC),
     &       SNGL(PLA(K)%O),SNGL(PLA(K)%TP),SNGL(PLA(K)%MU/MJUP),
     &       SNGL(STAR%MASS/SMAS),SNGL(CHI2),SNGL(FMAP) /)
        END DO
        DO I = 1,NMOD
          CALL ELEMENTS(PSAV(1:NPAR,I),NN,PLA%A,PLA%EXC,PLA%EXQ,
     &       PLA%CW,PLA%SW,CI2,SI2,PLA%CP,PLA%SP,PLA%TP,PLA%MU)
          IF (MULTIPLA) STAR%MASS = EXP(PSAV(NPAR,I))
          PLA%PER = DPI/NN
          PLA%W = ATAN2(PLA%SW,PLA%CW)
          PLA%PHI = ATAN2(PLA%SP,PLA%CP)
          PLA%OM = 0.5d0*(PLA%W+PLA%PHI)
          PLA%O = 0.5d0*(PLA%W-PLA%PHI)
          PLA%INC = 2.d0*ATAN2(SQRT(SI2),SQRT(CI2))
          CHI2 = PSAV(NPAR+1,I)
          FMAP = PSAV(NPAR+2,I)
          DO K = 1,NPLA
             DATA4(2*I+1,1:NSAV,K) = (/ SNGL(PLA(K)%A),SNGL(PLA(K)%PER),
     &       SNGL(PLA(K)%EXC),SNGL(PLA(K)%OM),SNGL(PLA(K)%INC),
     &       SNGL(PLA(K)%O),SNGL(PLA(K)%TP),SNGL(PLA(K)%MU/MJUP),
     &       SNGL(STAR%MASS/SMAS),SNGL(CHI2),SNGL(FMAP) /)
          END DO
          PLA%OM = MOD(PLA%OM+PI+PI,DPI)-PI
          PLA%O = MOD(PLA%O+PI+PI,DPI)-PI
          DO K = 1,NPLA
             DATA4(2*I+2,1:NSAV,K) = (/ SNGL(PLA(K)%A),SNGL(PLA(K)%PER),
     &       SNGL(PLA(K)%EXC),SNGL(PLA(K)%OM),SNGL(PLA(K)%INC),
     &       SNGL(PLA(K)%O),SNGL(PLA(K)%TP),SNGL(PLA(K)%MU/MJUP),
     &       SNGL(STAR%MASS/SMAS),SNGL(CHI2),SNGL(FMAP) /)
          END DO
        END DO
c...      Save everything into output file
        OPEN(18,FILE=DEV,STATUS='UNKNOWN')
        WRITE(18,*)2*(NMOD+1),NSAV,NPLA
        WRITE(18,*)DATA4
        CLOSE(18)
        DEALLOCATE(DATA4)
        DEALLOCATE(NN)
        DEALLOCATE(CI2)
        DEALLOCATE(SI2)
        END

C
C -----------------------------------------------------------------------------
C       Ecrire un ficher resultat
C -----------------------------------------------------------------------------
C
        SUBROUTINE WRITE_DISTRIB(NMOD,DEV)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    NMOD,         ! Number of models 
     &                  NSAV          ! Number of pars. to store per planet
        CHARACTER*(*) :: DEV          ! 
        REAL*8, DIMENSION(:), ALLOCATABLE ::
     &                  NN,           ! Moyen mouvement
     &                  CI,SI         ! cos(i), sin(i)
        REAL*8 ::       CHI2,         ! Chi2
     &                  FMAP,         ! MAP
     &                  SIGMA         ! Cumulative mass
        REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: DATA4
        INTEGER*4       I,K

        STAR%SIGJV = 0.d0
        NSAV = NEL+6
        ALLOCATE(DATA4(NMOD+1,NSAV,NPLA))
        ALLOCATE(NN(NPLA))
        ALLOCATE(CI(NPLA))
        ALLOCATE(SI(NPLA))
        CALL ELEMENTS(PSTART(1:NPAR),NN,PLA%A,PLA%EXC,PLA%EXQ,
     &       CI,SI,PLA%COM,PLA%SOM,PLA%CO,PLA%SO,PLA%TP,PLA%MU)

        STAR%MASS = EXP(PSTART(NPAR))
        STAR%V0 = PSTART(NEL*NPLA+1)
        IF (JITNUM.EQ.1) STAR%SIGJV = EXP(PSTART(NEL*NPLA+2))
        PLA%PER = DPI/NN
        PLA%OM = ATAN2(PLA%SOM,PLA%COM)
        PLA%O = ATAN2(PLA%SO,PLA%CO)
        PLA%INC = ATAN2(SI,CI)
        CHI2 = PSTART(NPAR+1)
        FMAP = PSTART(NPAR+2)
        DO K = 1,NPLA
           DATA4(1,1:NSAV,K) = (/ SNGL(PLA(K)%A),SNGL(PLA(K)%PER),
     &       SNGL(PLA(K)%EXC),SNGL(PLA(K)%OM),SNGL(PLA(K)%INC),
     &          SNGL(PLA(K)%O),SNGL(PLA(K)%TP),
     &          SNGL(PLA(K)%MU/MJUP),SNGL(STAR%V0/MPS),
     &          SNGL(STAR%SIGJV/MPS),SNGL(STAR%MASS/SMAS),
     &          SNGL(CHI2),SNGL(FMAP) /)
        END DO
        DO I = 1,NMOD
          CALL ELEMENTS(PSAV(1:NPAR,I),NN,PLA%A,PLA%EXC,PLA%EXQ,
     &       CI,SI,PLA%COM,PLA%SOM,PLA%CO,PLA%SO,PLA%TP,PLA%MU)
          STAR%MASS = EXP(PSAV(NPAR,I))
          STAR%V0 = PSAV(NEL*NPLA+1,I)
          IF (JITNUM.EQ.1) STAR%SIGJV = EXP(PSAV(NEL*NPLA+2,I))
          PLA%PER = DPI/NN
          PLA%OM = ATAN2(PLA%SOM,PLA%COM)
          PLA%O = ATAN2(PLA%SO,PLA%CO)
          PLA%INC = ATAN2(SI,CI)
          CHI2 = PSAV(NPAR+1,I)
          FMAP = PSAV(NPAR+2,I)
          DO K = 1,NPLA
             DATA4(I+1,1:NSAV,K) = (/ SNGL(PLA(K)%A),SNGL(PLA(K)%PER),
     &       SNGL(PLA(K)%EXC),SNGL(PLA(K)%OM),SNGL(PLA(K)%INC),
     &          SNGL(PLA(K)%O),SNGL(PLA(K)%TP),
     &          SNGL(PLA(K)%MU/MJUP),SNGL(STAR%V0/MPS),
     &          SNGL(STAR%SIGJV/MPS),SNGL(STAR%MASS/SMAS),
     &          SNGL(CHI2),SNGL(FMAP) /)
          END DO
        END DO
c...      Save everything into output file
        OPEN(18,FILE=DEV,STATUS='UNKNOWN')
        WRITE(18,*)NMOD+1,NSAV,NPLA
        WRITE(18,*)DATA4
        CLOSE(18)
        DEALLOCATE(DATA4)
        DEALLOCATE(NN)
        DEALLOCATE(CI)
        DEALLOCATE(SI)
        END


      
C
C-----------------------------------------------------------------------------
C    Displays a solution on screen
C-----------------------------------------------------------------------------
C

        SUBROUTINE DISPLAY_SOLUTION()

        USE DATA
        
        IMPLICIT NONE

        INTEGER*4 :: I          ! Planet index

 1      FORMAT(a25,' = ',f15.6,'   +/- ',f12.6,(a))

        DO I = 1,NPLA
           WRITE(SD,*)'-------------------'
           WRITE(SD,*)'Planet #',I, ' :'
           WRITE(SD,1)'Semi major axis ',PLA(I)%A,ABS(PLA(I)%DA),' AU'
           WRITE(SD,1)'Period ',PLA(I)%PER,ABS(PLA(I)%DPER),' days'
           WRITE(SD,1)'Eccentricity ',PLA(I)%EXC,ABS(PLA(I)%DEXC)
           IF (RADVEL) THEN
              WRITE(SD,1)'Argument of periastron (omega) ',
     &             PLA(I)%OM,ABS(PLA(I)%DOM),' °'
              WRITE(SD,1)'inclination ',
     &              PLA(I)%INC,ABS(PLA(I)%DINC),' °'
              WRITE(SD,1)'Longitude of ascending node (Omega) ',
     &             PLA(I)%O,ABS(PLA(I)%DOO),' °'
           ELSE
              WRITE(SD,1)'Argument of periastron (omega) ',
     &                     PLA(I)%OM,ABS(PLA(I)%DOM),' ° (mod 180°)'
              WRITE(SD,1)'inclination ',
     &              PLA(I)%INC,ABS(PLA(I)%DINC),' °'
              WRITE(SD,1)'Longitude of ascending node (Omega) ',
     &             PLA(I)%O,ABS(PLA(I)%DOO),' ° (mod 180°)'
           END IF
              
           WRITE(SD,1)'Time of periastron ',PLA(I)%TP,ABS(PLA(I)%DTP)
           WRITE(SD,1)'Dynamical mass of orbit ',PLA(I)%MDYN/SMAS,
     &          ABS(PLA(I)%DMDYN/SMAS),' Msun'
           IF (MULTIPLA) WRITE(SD,1)'Mass :',PLA(I)%MU/MJUP,
     &                ABS(PLA(I)%DMU)/MJUP,' Mjup'
        END DO
        WRITE(SD,*)'-------------------'
        IF (MULTIPLA) WRITE(SD,1)'Stellar mass ',STAR%MASS/SMAS,
     &       STAR%DMASS/SMAS,' Msun'
        IF (RADVEL) THEN
           WRITE(SD,1)'Offset velocity ',STAR%V0/MPS,
     &       STAR%DV0/MPS,' km/s'
           IF (JITNUM.EQ.1) WRITE(SD,1)'Velocity Jitter',STAR%SIGJV/MPS,
     &         STAR%DSIGJV/MPS,' km/s'
        END IF

        END 

c....
c     WRITe(SD,*)STAR%MASS/SMAS,STAR%V0/MPS,STAR%SIGJV/MPS
c          DO I = 1,NPLA
c             WRITE(SD,*)PLA(I)%A,PLA(I)%PER,PLA(I)%EXC,PLA(I)%OM,
c     &            PLA(I)%INC,PLA(I)%O,
c     &            PLA(I)%TP,PLA(I)%MU/MJUP
c          END DO
        
C
C-----------------------------------------------------------------------------
C    Exiting simulation
C-----------------------------------------------------------------------------
C

        SUBROUTINE UTIL_EXIT()

        USE DATA
        
        IMPLICIT NONE

        WRITE(SD,*)'--------------------'
        WRITE(SD,*)'Computation finished'
        WRITE(SD,*)'--------------------'
        CLOSE(SDL)

        END
