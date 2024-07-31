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
