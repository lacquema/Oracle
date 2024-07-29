C
C      Levenberg-Marquardt common routines      
C      
C -----------------------------------------------------------------------------
C       Search for minimal chi2 by Levenberg-Marquardt method
C -----------------------------------------------------------------------------
C
        SUBROUTINE MRQMIN(N,P,COV,CHI2,FMAP)

        USE DATA
        
        IMPLICIT NONE

        REAL*8, PARAMETER :: LAMBDA0 = 1d-3
        INTEGER*4       N       ! Dimension of parameter space

        REAL*8          ALPHA(N,N),     ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  COV(N,N),       ! Covariance matrix
     &                  SIGP(N),        ! Variance for parameters
     &                  DP(N),          ! Delta-Parameters
     &                  P(N),           ! Parameters
     &                  CHI2,FMAP,      ! Chi2 (total), FMAP
     &                  LAMBDA          ! Parameterr of method
        INTEGER*4       I,ITS,iv,jv
        LOGICAL         OK,TEST

        LAMBDA = LAMBDA0
        CALL MRQCOF(N,P,ALPHA,BETA,CHI2,FMAP)
        WRITE(SD,*)'Initial chi2 = ',chi2
        OK = .FALSE.
        ITS = 0
        DO WHILE(.NOT.OK)        
           CALL MRQSTEP(N,LAMBDA,P,ALPHA,BETA,CHI2,FMAP,OK)
           IF ((MOD(ITS,1000).EQ.0).or.(lambda.lt.1d-4)) THEN
            WRITE(SD,*)'        ----> ',sngl(CHI2),' <----',sngl(lambda)
         END IF
           ITS = ITS+1
           if (chi2.ne.chi2) stop
        END DO
        WRITE(SD,*)'End of Levenberg-Marquardt'
        CALL GAUSS(N,ALPHA,BETA,COV,DP,TEST)
        DO I=1,N
          SIGP(I)=SQRT(COV(I,I))
        END DO
        END

C
C -----------------------------------------------------------------------------
C       One step of Levenberg-Marquardt method
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQSTEP(N,LAMBDA,P,ALPHA,BETA,CHI2,FMAP,OK)

        USE DATA

        IMPLICIT NONE

        REAL*8, PARAMETER :: LAMBDA0 = 1d-3
        INTEGER*4       N       ! Dimension of parameter space
        REAL*8          A(N,N),ATRY(N,N),! Trial matrix
     &                  ALPHA(N,N),     ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  B(N),BTRY(N),   ! Trial vector
     &                  COV(N,N),       ! Covariance matrix (useless here)
     &                  DP(N),          ! Delta-Parameters
     &                  PTRY(N),        ! Trial new parameters
     &                  P(N),           ! Parameter set
     &                  CHI2,FMAP,      ! Chi2,MAP (CHi2 = Real Chi2 + FMAP)
     &                  CTRY,           ! Trial Chi2
     &                  LAMBDA          ! parameter for mrq method

        LOGICAL         OK,TEST,chk     ! Test for stopping
        INTEGER*4       I,J

        A(1:N,1:N) = ALPHA(1:N,1:N)
        B(1:N) = BETA(1:N)
c...  Multiply diagonal by 1+lambda
        DO I=1,N
           A(I,I) = A(I,I)*(1.0d0+LAMBDA)
        END DO        
c        WRITE(SD,*)'a',sngl(a),'b',sngl(b)
        CALL GAUSS(N,A,B,COV,DP,TEST)
c...   Trial step
        
        PTRY(1:N) = P(1:N)+DP(1:N)
c        do i=1,n
c          print*,i,'dp',sngl(dp(i)),'p',sngl(p(i)),sngl(ptry(i))
c        end do
c...  Check if trial step is still in the authorized bounds
        CALL MRQCHECK(N,PTRY,CHK)
        
c...  Compute coefs and trial Chi2. If error then assume bad step
        IF (.NOT.CHK) THEN
          CALL MRQCOF(N,PTRY,ATRY,BTRY,CTRY,FMAP)        
          IF (CTRY.NE.CTRY) CHK = .TRUE. 
        END IF
c...  We are sill in valid region
        IF (.NOT.CHK) THEN
c...    Degraded CHi2. Refuse step and increase lambda
          IF (CTRY.GT.CHI2) THEN
            IF (LAMBDA.EQ.0.d0) THEN
              LAMBDA = LAMBDA0
            ELSE
              LAMBDA = 10.d0*LAMBDA
            END IF
            OK=.FALSE.
c...    Improved Chi2. Accept step and decrease lambda
          ELSE
            P(1:N) = PTRY(1:N)
            BETA(1:N) = BTRY(1:N)
            ALPHA(1:N,1:N) = ATRY(1:N,1:N)
            LAMBDA = 0.1d0*LAMBDA
            OK = (ABS(CHI2-CTRY).LT.EPS)
            CHI2 = CTRY
         END IF
c... Case where the trial step is out of bounds or erroneous. Restart
        ELSE  
          OK = .FALSE.
          CALL MRQRENEW(N,P)
          CALL MRQCOF(N,P,ALPHA,BETA,CHI2,FMAP) 
          LAMBDA = LAMBDA0
c          CHK = .TRUE.
c          DO WHILE(CHK)
c            A(1:N,1:N) = ALPHA(1:N,1:N)
c            B(1:N) = BETA(1:N)
c...  Multiply diagonal by 1+lambda
c            DO I=1,N
c              A(I,I) = A(I,I)*(1.0d0+LAMBDA)
c            END DO        
c        WRITE(SD,*)'a',sngl(a),'b',sngl(b)
c            CALL GAUSS(N,A,B,COV,DP,TEST)
c            PTRY(1:N) = P(1:N)+DP(1:N)
c            CALL MRQCHECK(N,PTRY,CHK)
c            IF (CHK) LAMBDA = LAMBDA*10.d0
c          END DO
        END IF
          
        END

C
C -----------------------------------------------------------------------------
C       Computation of HCI contribution to Levenberg-Marquardt coefficients
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQASTROM(N,P,ALPHA,BETA,CHI2)

        USE DATA

        IMPLICIT NONE

        INTEGER*4       N               ! Dimension of parameter space
        
        REAL*8, DIMENSION(NPAR,NPLA) :: DXF,DYF
        REAL*8, DIMENSION(NPLA) :: XF,YF
        REAL*8 ::       ALPHA(N,N),     ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  P(N),           ! The parameters
     &                  CHI2,           ! Chi2
     &                  F1X,F1Y,F2XX,F2YY,F2XY, ! Intermediaires
     &                  SXM1,SYM1,RHO,UF,!
     &                  DX,DY           ! Normalized differences   
        INTEGER*4       I,J,K,L,SJ,iv,jv  

c        REAL*8 ::       PB(NPAR),
c     &       xfz(npla),yfz(npla),xfq(npla),yfq(npla),es
        
c        pb(1:n) = p(1:n)
        DXF = 0.d0
        DYF = 0.d0
c        es = 1.d-8
        DO L = 1,NPLA
           DO I = 1,PLA(L)%NDATAS
c         do iv=1,n
c         pb(iv) = p(iv)+es
c         call mrqposfit(n,pla(l)%tas(i),pb(1:n),xfz,yfz,
c     &            dxf(1:n,1:npla),dyf(1:n,1:npla))
c         pb(iv) = p(iv)-es
c         call mrqposfit(n,pla(l)%tas(i),pb(1:n),xfq,yfq,
c     &            dxf(1:n,1:npla),dyf(1:n,1:npla))

              CALL MRQPOSFIT(N,PLA(L)%TAS(I),P(1:N),
     &                      XF(1:NPLA),YF(1:NPLA),
     &                      DXF(1:N,1:NPLA),DYF(1:N,1:NPLA))
c              WRITE(SD,*)iv,'x', ! sngl(xfz),sngl(xfq),
c     &         sngl((xfz-xfq)/(2.*es)),
c     &                 sngl(dxf(iv,:))
c              WRITE(SD,*)iv,'y', ! sngl(yfz),sngl(yfq),
c     &         sngl((yfz-yfq)/(2.*es)),
c     &                 sngl(dyf(iv,:))
c         end do
c        stop

c         print*,sngl(pla(l)%x(i)/star%dist*1d3),
c     &            sngl(pla(l)%y(i)/star%dist*1d3),
c     &            sngl(xf(l)/star%dist*1d3),
c     &            sngl(yf(l)/star%dist*1d3)    
             
            SXM1 = PLA(L)%SIGXM1(I)
            SYM1 = PLA(L)%SIGYM1(I)
            RHO = PLA(L)%RHOXY(I)
            UF = PLA(L)%UFAC(I)
            DX = (PLA(L)%X(I)-XF(L))*SXM1
            DY = (PLA(L)%Y(I)-YF(L))*SYM1       
            CHI2 = CHI2 + (DX*DX+DY*DY-2.d0*RHO*DX*DY)*UF
            F1X = (DX-RHO*DY)*SXM1*UF  
            F1Y = (DY-RHO*DX)*SYM1*UF
            F2XX = PLA(L)%SIGXM2(I)*UF
            F2YY = PLA(L)%SIGYM2(I)*UF
            F2XY = -RHO*SXM1*SYM1*UF
            BETA(1:N) = BETA(1:N)+F1X*DXF(1:N,L)+F1Y*DYF(1:N,L)
            DO J=1,N
                ALPHA(J:N,J) = ALPHA(J:N,J)
     &             +F2XX*DXF(J,L)*DXF(J:N,L)+F2YY*DYF(J,L)*DYF(J:N,L)
     &             +F2XY*(DXF(J,L)*DYF(J:N,L)+DXF(J:N,L)*DYF(J,L))
            END DO        
          END DO
        END DO
c      WRITE(SD,*)chi2,'r'
c        stop

        DO J=2,N
          ALPHA(1:J-1,J) = ALPHA(J,1:J-1)
        END DO 

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
          DFMAP(IM+1) = 2.d0
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
C       Radial velocity fit and derivatives
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
           P(DKAL+4) = PLA(I)%TI*PLA(I)%CO
           P(DKAL+5) = PLA(I)%TI*PLA(I)%SO
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
C     
C -----------------------------------------------------------------------------
C       Calculation of Chi2 and coefficients for Levenberg-Marquardt
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQCOF(N,P,ALPHA,BETA,CHI2,FMAP)

        USE DATA

        IMPLICIT NONE

        INTEGER*4       N               ! Dimension of parameter space

        REAL*8, DIMENSION(NPAR,NPLA) :: DXF,DYF
        REAL*8, DIMENSION(NPAR) :: DVF
        REAL*8, DIMENSION(NPLA) :: XF,YF
        REAL*8, DIMENSION(NPAR,NPAR) :: D2FMAP  ! Second derivatives / MAP
        REAL*8, DIMENSION(NPAR) :: DFMAP        ! First derivatives /MAP
        REAL*8 ::       ALPHA(N,N), ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  P(N),           ! The parameters
     &                  CHI2,           ! Chi2
     &                  FMAP,           ! Map merit function (added to Chi2)
     &                  F1X,F1Y,F2XX,F2YY,F2XY,F1,F2, ! Intermediaires
     &                  SIGJV,          ! Velocity Jitter
     &                  SIGV2P,            ! Incertitude^2 augmentée
     &                  FS,
     &                  SXM1,SYM1,RHO,UF,  !
     &                  VF,                ! Fitted velocity
     &                  DX,DY,DV           ! Ecarts normalisés   
        
        
c        REAL*8 ::       PB(NPAR),es,
c     &                  xfz(npla),yfz(npla),xfq(npla),yfq(npla), 
c     &                  fmappp,fmapmp,fmappm,fmapmm,
c     &                  vfq,vfz,fmapz,fmapq ! Vitesse radiale
        INTEGER*4 ::     I,J,K,L,SJ,iv,jv          
        LOGICAL ::       JITTER
c     
        CHI2 = 0.d0       
        ALPHA(1:N,1:N) = 0.d0
        BETA(1:N) = 0.d0
        JITTER = RADVEL.AND.(JITNUM.EQ.1).AND.(N.EQ.NPAR)
        IF (RADVEL) CALL MRQRADVEL(N,JITTER,P(1:N),
     &                     ALPHA(1:N,1:N),BETA(1:N),CHI2)

        CALL MRQASTROM(N,P,ALPHA(1:N,1:N),BETA(1:N),CHI2)

c...  Map merit function & derivatives
c        pb(1:n) = p(1:n)
c           es=1d-5
c           do iv = 1,n
c             do jv = 1,n
c                pb = p
c                if (jv.ne.iv) then
c                  pb(iv)=p(iv)+es
c                  pb(jv)=p(jv)+es
c                  call mrqmap(pb,fmappp,dfmap,d2fmap)
c                  pb(iv) = p(iv)-es
c                  pb(jv) = p(jv)+es
c                  call mrqmap(pb,fmapmp,dfmap,d2fmap)
c                  pb(iv)=p(iv)+es
c                  pb(jv)=p(jv)-es
c                  call mrqmap(pb,fmappm,dfmap,d2fmap)
c                  pb(iv) = p(iv)-es
c                  pb(jv) = p(jv)-es
c                  call mrqmap(pb,fmapmm,dfmap,d2fmap)
c                  CALL MRQMAP(P,FMAP,DFMAP,D2FMAP)
c                if (d2fmap(iv,jv).ne.0.d0) write(sd,*)iv,jv,
c     &               sngl(0.25d0*(fmappp-fmapmp-fmappm+fmapmm)/(es*es)),
c     &               sngl(d2fmap(iv,jv))
c             else
c               pb(iv) =p(iv)+es 
c               call mrqmap(pb,fmapz,dfmap,d2fmap)
c               pb(iv)=p(iv)-es
c               call mrqmap(pb,fmapq,dfmap,d2fmap)
c               call mrqmap(p,fmap,dfmap,d2fmap)
c               if (d2fmap(iv,iv).ne.0.d0) write(sd,*)iv,iv,
c     &              sngl((fmapq+fmapz-2.d0*fmap)/(es*es)),
c     &               sngl(d2fmap(iv,iv))               
c             end if
c            end do
c          end do
c            stop
        CALL MRQMAP(JITTER,N,P(1:N),FMAP,DFMAP(1:N),D2FMAP(1:N,1:N))
        CHI2 = CHI2+FMAP        ! Her we take into account MAP in Chi2
        BETA(1:N) = BETA(1:N)-0.5d0*DFMAP(1:N) 
        ALPHA(1:N,1:N) = ALPHA(1:N,1:N)+0.5d0*D2FMAP(1:N,1:N)
        END    

C
C     
C -----------------------------------------------------------------------------
C    Calculation of FMAP & derivatives
C        MAP = exp(-Chi2/2)*f(parameters) (f = priors)
C     => Replace Chi2 with Chi2-2*ln(f(parameters)) = -2*ln(MAP) = Chi2+FMAP    
C-----------------------------------------------------------------------------
C

        SUBROUTINE MRQMAP(JITTER,N,P,FMAP,DFMAP,D2FMAP)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 :: N                  ! Dimension of parameter space
        REAL*8, DIMENSION(N,N) :: D2FMAP ! Second derivatives
        REAL*8, DIMENSION(N) :: DFMAP   ! First derivatives
        REAL*8, DIMENSION(N) :: P       ! Parameters
        REAL*8 ::       FMAP,           ! Merit function
     &                  SIGJV,SIGJV2,   ! Jitter,Jitter^2
     &                  SIGV2P          ! Jitter^2+Error_i^2

        INTEGER*4 ::    I,SJ          
        LOGICAL ::      JITTER
c     
        IF (JITTER) THEN
          SJ = NEL*NPLA+2
          SIGJV = EXP(P(SJ))        ! Velocity Jitter
          SIGJV2 = SIGJV*SIGJV
        END IF

        FMAP = 0.d0
        DFMAP = 0.d0
        D2FMAP = 0.d0

        CALL MRQMAPASTROM(N,P,FMAP,DFMAP(1:N),D2FMAP(1:N,1:N))
        
c... Taking into account prior on Stellar Jitter in FMAP (+2*ln(j+j0))
        IF (JITTER) THEN
          FMAP = FMAP+2.d0*LOG(SIGJV+SIGV0)
          DFMAP(SJ) = DFMAP(SJ)+2.d0*SIGJV/(SIGJV+SIGV0)
          D2FMAP(SJ,SJ) = D2FMAP(SJ,SJ)
     &                    +2.d0*SIGJV*SIGV0/(SIGJV+SIGV0)**2
c...  Taking into account effect of Stellar Jitter in p(data|model)
c...        = Product(cte*1/sqrt(j^2+sig_i^2)) => +2*1/2*ln(j^2+sig_i^2)        
          DO I = 1,STAR%NDATVR
            SIGV2P = STAR%SIGV2(I)+SIGJV2
            FMAP = FMAP+LOG(1.d0+SIGJV2/STAR%SIGV2(I))
c            FMAP = FMAP+LOG(SIGJV2+STAR%SIGV2(I))
            DFMAP(SJ) = DFMAP(SJ)+2.d0*SIGJV2/SIGV2P  !!!!!!
            D2FMAP(SJ,SJ) = D2FMAP(SJ,SJ)
     &                  +4.d0*STAR%SIGV2(I)*SIGJV2/(SIGV2P*SIGV2P)
          END DO
        END IF

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
      
