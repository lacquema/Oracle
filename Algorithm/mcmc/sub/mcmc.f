C
C   Common MCMC routines 
C      
C
C -----------------------------------------------------------------------------
C       General MCMC routine (sequential)
C -----------------------------------------------------------------------------
C

        SUBROUTINE MCMC(NMOD,NEW,FILES)

        USE DATA

        IMPLICIT NONE
     
        CHARACTER*(*),DIMENSION(NFIL) :: FILES

        TYPE(CHAIN), DIMENSION(NCH) :: CHN
        TYPE(CHAIN) ::  CH0,CHC      ! Starting & current chain  
        INTEGER*4 ::    ERROR,
     &                  C,IU,JU,I,K, ! Indexes
     &                  COUNTOK,     ! Counter
     &                  LL(NPAR),    ! Permutation of integers
     &                  NMOD         ! # of models
        INTEGER*8 ::    LONG,
     &                  NGEL        ! Time for Gelman-Rubin
        LOGICAL ::      NEW,         ! True if new computation
     &                  ACC,         ! Acceptance flag
     &                  VERBOSE,     ! Flag for print
     &                  INFER,       ! Inference flag
     &                  INFCRIT,     ! Instantaneous inference flag
     &                  OK           ! Convergence flag
        REAl*8 ::       PROBINFER,   ! Probability of retaining a model
     &                  DIRS(NPAR,NPAR), ! Directions
     &                  DP(NPAR)     ! Departure from model

        CALL INIT_MCMC(PROBINFER,CHN,CH0,DIRS,
     &              NMOD,COUNTOK,LONG,NGEL,INFER,NEW)
        ALLOCATE(CHC%P(NPAR))
        ALLOCATE(CHC%MOY(NPAR))
        ALLOCATE(CHC%BETA(NPAR))        
        ALLOCATE(CHC%COV(NPAR,NPAR))
        
        DO WHILE ((COUNTOK.LT.5).OR.(NMOD.LT.NMAX))  
          CALL PERMUTATION(NPAR,LL)          
          DO C = 1,NCH
            DO I = 1,NPAR
              CALL METRO_STEP(CHN(C)%P,DIRS(1:NPAR,LL(I)),
     &                        CHN(C)%CHI2,CHN(C)%BETA(LL(I)),ACC)       
            END DO
            DP = CHN(C)%P-CHN(C)%MOY
            CHN(C)%MOY = CHN(C)%MOY+DP/(LONG+1_8)
            DO IU=1,NPAR          
              DO JU=1,IU
                IF (LONG.EQ.1_8) THEN
                  CHN(C)%COV(IU,JU) = DP(IU)*DP(JU)
                ELSE
                  CHN(C)%COV(IU,JU) = (LONG-1_8)*CHN(C)%COV(IU,JU)/LONG 
     &                           + DP(IU)*DP(JU)/(LONG+1_8)
                END IF
                CHN(C)%COV(JU,IU) = CHN(C)%COV(IU,JU)
              END DO
            END DO            
          END DO
          LONG = LONG+1_8

c....  Gelman-Rubin statistics

          IF (LONG.EQ.NGEL) THEN
            CALL GELMAN_RUBIN(PROBINFER,CHN,NMOD,LONG,INFER,OK,FILES)

            IF (.NOT.OK) THEN
              NGEL = LONG+100_8
              COUNTOK = 0
            ELSE
              IF (COUNTOK.LT.10) THEN
                NGEL = NINT(LONG*1.01d0,KIND=8)
              ELSE
                NGEL = MIN(LONG+100_8,NINT(LONG*1.01d0,KIND=8))
              END IF
              COUNTOK = COUNTOK+1
              WRITE(SD,*)'countok = ',countok
            END IF
          END IF

c...  Compute new directions          

          IF ((NDIAG.GT.0).AND.(MOD(LONG,NDIAG).EQ.0)) THEN
             CALL CHANGE_DIRS(LONG,CHN,DIRS)
          END IF
          
c......... Resume current status => dump file

          IF (MOD(LONG,FREQDUMP).LE.1) THEN
            WRITE(SD,'(a)',ADVANCE='NO')'Dump file created...'
            CALL  WRITE_DUMP(NMOD,COUNTOK,LONG,NGEL,INFER,
     &         PROBINFER,CHN,DIRS,CH0%BETA,FILES)
            WRITE(SD,*)'ok'
          END IF

        END DO

        DO C = 1,NCH
           DEALLOCATE(CHN(C)%P)
           DEALLOCATE(CHN(C)%MOY)
           DEALLOCATE(CHN(C)%BETA)
           DEALLOCATE(CHN(C)%COV)
        END DO
        DEALLOCATE(CH0%P)
        DEALLOCATE(CH0%MOY)
        DEALLOCATE(CH0%BETA)        
        DEALLOCATE(CH0%COV)
        DEALLOCATE(CHC%P)
        DEALLOCATE(CHC%MOY)
        DEALLOCATE(CHC%BETA)        
        DEALLOCATE(CHC%COV)                
        END


C     
C -----------------------------------------------------------------------------
C       General MCMC routine (parallel version OPEN-MP)
C -----------------------------------------------------------------------------
C
      
        SUBROUTINE MCMCP(NMOD,NEW,FILES)

        USE DATA

        IMPLICIT NONE

        CHARACTER*(*),DIMENSION(NFIL) :: FILES

        TYPE(CHAIN), DIMENSION(NCH) :: CHN
        TYPE(CHAIN) ::  CH0          ! Starting chain   
        INTEGER*4 ::    ERROR,
     &                  C,IU,JU,I,J,K, ! Indexes
     &                  COUNTOK,     ! Counter
     &                  LL(NPAR,200),! Permutation of integers
     &                  LLC(NPAR),
     &                  NMOD         ! # of models
        INTEGER*8 ::    LONG,LONGC,
     &                  NGROUP,JGROUP,! Number of grouped steps (=NGEL..)
     &                  NGEL         ! Time for Gelman-Rubin
        LOGICAL ::      NEW,         ! True if new computation
     &                  ACC,         ! Acceptance flag
     &                  VERBOSE,     ! Flag for print
     &                  INFER,       ! Inference flag
     &                  INFCRIT,     ! Instantaneous inference flag
     &                  OK           ! Convergence flag
        REAl*8 ::       PROBINFER,   ! Probability of retaining a model
     &                  DIRS(NPAR,NPAR), ! Directions
     &                  DIRC(NPAR,NPAR), ! Directions
     &                  DMOY(NPAR),! Mean of chains
     &                  DP(NPAR) ! Departure from model
        REAL*8 :: CHCMOY(NPAR),CHCBETA(NPAR),CHCP(NPAR) ! Current partial chain
        REAL*8 :: CHCCOV(NPAR,NPAR),CHCCHI2

        CALL INIT_MCMC(PROBINFER,CHN,CH0,DIRS,
     &              NMOD,COUNTOK,LONG,NGEL,INFER,NEW)
        NGROUP = NGEL
        DO WHILE ((COUNTOK.LT.5).OR.(NMOD.LT.NMAX))
           JGROUP = NGEL/NGROUP
           DO J = 1,JGROUP
              DO K = 1,NGROUP
                 CALL PERMUTATION(NPAR,LL(1:NPAR,K))
              END DO
!$omp parallel do
!$omp& default(shared)
!$omp& private(i,k,iu,ju,acc,dirc,dmoy,llc,dp)
!$omp& private(chcmoy,chcbeta,chcp,chccov,chcchi2)
              DO C = 1,NCH
                 CHCMOY = 0.d0
                 CHCCOV = 0.d0
                 CHCP = CHN(C)%P
                 CHCCHI2 = CHN(C)%CHI2
                 DO K = 1,NGROUP
                    LLC(1:NPAR) = LL(1:NPAR,K)
                    DO I = 1,NPAR
                       DIRC(1:NPAR,I) = DIRS(1:NPAR,LLC(I))
                       CHCBETA(I) = CHN(C)%BETA(LLC(I))
                    END DO
                    DO I = 1,NPAR
                       CALL METRO_STEP(CHCP,DIRC(1:NPAR,I),
     &                      CHCCHI2,CHCBETA(I),ACC)       
                    END DO
                    DP = CHCP-CHCMOY
                    CHCMOY = CHCMOY+DP/K
                    IF (K.GT.1) THEN
                       DO IU=1,NPAR 
                          DO JU=1,IU
                             CHCCOV(IU,JU) = (K-2)*CHCCOV(IU,JU)/(K-1)
     &                            + DP(IU)*DP(JU)/K
                             CHCCOV(JU,IU) = CHCCOV(IU,JU)
                          END DO
                       END DO                
                    END IF
                 END DO 
                 CHN(C)%P = CHCP
                 CHN(C)%CHI2 = CHCCHI2
                 DMOY = CHCMOY-CHN(C)%MOY
                 CHN(C)%MOY = CHN(C)%MOY+NGROUP*DMOY/(NGROUP+LONG)
                 DO IU=1,NPAR  
                    DO JU=1,IU
                       CHN(C)%COV(IU,JU) = ((LONG-1_8)*CHN(C)%COV(IU,JU)
     &                    + (NGROUP-1_8)*CHCCOV(IU,JU)
     &                    + DMOY(IU)*DMOY(JU)*NGROUP*LONG/(NGROUP+LONG))
     &                    /(NGROUP+LONG-1_8)
                       CHN(C)%COV(JU,IU) = CHN(C)%COV(IU,JU)
                    END DO
                 END DO
              END DO
!$omp end parallel do
              LONG = LONG+NGROUP
           END DO

c....  Gelman-Rubin statistics

           CALL GELMAN_RUBIN(PROBINFER,CHN,NMOD,LONG,INFER,OK,FILES)

c.......
           IF (.NOT.OK) THEN
              NGEL = 100_8
              COUNTOK = 0
           ELSE
              IF (COUNTOK.LT.10) THEN
                 NGEL = NINT(LONG*1d-2/NGROUP,KIND=8)*NGROUP
              ELSE
                 NGEL = MIN(NGROUP,NINT(LONG*1d-2/NGROUP,KIND=8)*NGROUP)
              END IF
              COUNTOK = COUNTOK+1
              WRITE(SD,*)'countok = ',countok
           END IF

c...  Compute new directions

           IF ((NDIAG.GT.0).AND.(MOD(LONG,NDIAG).LT.10)) THEN
             CALL CHANGE_DIRS(LONG,CHN,DIRS)
           END IF

c......... Resume current status => dump file

          IF (MOD(LONG,FREQDUMP).LT.10) THEN
            WRITE(SD,'(a)',ADVANCE='NO')'Dump file created...'
            CALL  WRITE_DUMP(NMOD,COUNTOK,LONG,NGEL,INFER,
     &           PROBINFER,CHN,DIRS,CH0%BETA,FILES)
            WRITE(SD,*)'ok'
          END IF

        END DO

        DO C = 1,NCH
           DEALLOCATE(CHN(C)%P)
           DEALLOCATE(CHN(C)%MOY)
           DEALLOCATE(CHN(C)%BETA)
           DEALLOCATE(CHN(C)%COV)
        END DO
        DEALLOCATE(CH0%P)
        DEALLOCATE(CH0%MOY)
        DEALLOCATE(CH0%BETA)        
        DEALLOCATE(CH0%COV)
        
        END
C     
C -----------------------------------------------------------------------------
C       Initialization of MCMC
C -----------------------------------------------------------------------------
C
      
        SUBROUTINE INIT_MCMC(PROBINFER,CHN,CH0,DIRS,
     &              NMOD,COUNTOK,LONG,NGEL,INFER,NEW)

        USE DATA

        IMPLICIT NONE

        TYPE(CHAIN), DIMENSION(NCH) :: CHN
        TYPE(CHAIN)     CH0          ! Starting chain   
        INTEGER*4 ::    C,I,         ! Indexes
     &                  COUNTOK,     ! Counter
     &                  NMOD         ! # of models
        INTEGER*8 ::    LONG,
     &                  NGEL         ! Time for Gelman-Rubin
        LOGICAL ::      INFER,       ! Inference flag
     &                  FAIL,        ! True if serach beta failed
     &                  NEW          ! True if new simulation
        REAl*8 ::       PROBINFER,   ! Probability of retaining a model
     &                  DIRS(NPAR,NPAR)  ! Directions

        write(sd,*)'init_mcmc'
        DO C = 1,NCH
           ALLOCATE(CHN(C)%P(NPAR))
           ALLOCATE(CHN(C)%MOY(NPAR))
           ALLOCATE(CHN(C)%BETA(NPAR))
           ALLOCATE(CHN(C)%COV(NPAR,NPAR))
        END DO
        ALLOCATE(CH0%P(NPAR))
        ALLOCATE(CH0%MOY(NPAR))
        ALLOCATE(CH0%BETA(NPAR))        
        ALLOCATE(CH0%COV(NPAR,NPAR))
        
        IF (NEW) THEN
           FAIL = .TRUE.
           DO WHILE(FAIL)
              DIRS = 0.d0
              DO I=1,NPAR
                 DIRS(I,I) = 1.d0
              END DO
              CH0%BETA = 0.5d0
              CH0%P = PSTART(1:NPAR)
              CH0%MOY = PSTART(1:NPAR)
              CH0%COV = 0.d0
              CALL CHI2SEUL(CH0%P,CH0%CHI2)
c              print*,'chi2seul',ch0%chi2
              CALL SEARCH_BETA(CH0%P,DIRS(1:NPAR,1:NPAR),CH0%BETA,FAIL)
           END DO
           WRITE(SD,*)'beta initiation successful ! '
           WRITE(SD,*)'beta values = ',sngl(ch0%beta)
           DO C=1,NCH
              CHN(C) = CH0             
           END DO
           LONG = 1_8
           INFER = .FALSE.
           NGEL = 100_8
           COUNTOK = 0
           NMOD = 0
           WRITE(SD,*)'Starting MCMC...'
           WRITE(SD,*)sngl(chn(1:nch)%chi2)
        ELSE
           CALL READ_DUMP(NMOD,COUNTOK,LONG,NGEL,INFER,
     &        PROBINFER,CHN,DIRS,CH0%BETA)
           CH0%P = PSTART(1:NPAR)
           CH0%MOY = PSTART(1:NPAR)
           CH0%COV = 0.d0
           CALL CHI2SEUL(CH0%P,CH0%CHI2)
           WRITE(SD,*)"Here we go again, folks !"
        END IF

        END  
C     
C -----------------------------------------------------------------------------
C       Gelman-Rubin statistics & inference 
C -----------------------------------------------------------------------------
C
      
        SUBROUTINE GELMAN_RUBIN(PROBINFER,CHN,NMOD,LONG,INFER,OK,FILES)

        USE DATA

        IMPLICIT NONE

        CHARACTER*(*),DIMENSION(NFIL) :: FILES

        TYPE(CHAIN), DIMENSION(NCH) :: CHN
        INTEGER*4 ::    C,I,IU,JU,   ! Indexes
     &                  NMOD         ! # of models
        INTEGER*8 ::    LONG
        LOGICAL ::      INFER,       ! Inference flag
     &                  VERBOSE,     ! Flag for print
     &                  INFCRIT,     ! Instantaneous inference criterion 
     &                  OK           ! True if convergence reached
        REAL*8 ::       PROBINFER,   ! Probability of retaining a model
     &                  COMPUTE_MAP, ! Function => MAP
     &                  R,W,B,VARP,TS,     ! TEsts for Gelman-Rubin
     &                  RR,RS,          ! Random number
     &                  MOYY(NPAR)         ! Global mean over chains
      
        OK = .TRUE.
        VERBOSE = (MOD(LONG,VERBFREQ).LT.2)
        INFCRIT = .TRUE.
        DO IU = 1,NPAR
           MOYY(IU) = 0.d0
           DO I = 1,NCH              
              MOYY(IU) = MOYY(IU)+CHN(I)%MOY(IU)
           END DO
           MOYY(IU) = MOYY(IU)/NCH
           W = 0.d0
           B = 0.d0
           DO I = 1,NCH
              B = B+(CHN(I)%MOY(IU)-MOYY(IU))**2
              W = W+CHN(I)%COV(IU,IU)
           END DO
           W = W/NCH
           B = B*LONG/(NCH-1)
           VARP = (LONG-1_8)*W/LONG+B/LONG
           R = SQRT(VARP/W)  ! Gelman-Rubin Criterion
           TS = LONG*NCH*MIN(VARP/B,1.d0)
           IF (VERBOSE) WRITE(SD,*)iu,'r',sngl(r),'ts',sngl(ts)
           OK = OK.AND.((R.LT.RCONV).AND.(TS.GT.TSCONV))
           INFCRIT = INFCRIT.AND.((R.LT.RINF).AND.(TS.GT.TSINF))
        END DO
        IF (VERBOSE) WRITE(SD,*)'nmod = ',nmod,
     &          'long = ',long,' chi2 = ',sngl(chn(1:nch)%chi2)
        IF (INFCRIT.AND.(.NOT.INFER)) THEN
           INFER = .TRUE.
           PROBINFER = MAX(NINFER*100.d0/(4.d0*LONG),MINPROB)
           WRITE(SD,*)
     &          'Inference now, sampling probability = ',probinfer
        END IF
        
c...  Inférence quand le critère est atteint

        IF (INFCRIT) THEN
           CALL RANDOM_NUMBER(RR)
           IF (RR.LT.PROBINFER) THEN
              DO C = 1,NCH
                 IF (NMOD.LT.NMAX) THEN
                    NMOD = NMOD+1
                    I = NMOD
                    if (mod(nmod,1000).eq.0)
     &                    WRITE(SD,*)'modele ',nmod
                 ELSE
                    CALL RANDOM_NUMBER(RS)
                    I = FLOOR(RS*NINFER)*NCH+C
                 END IF
                 PSAV(1:NPAR,I) = CHN(C)%P
                 PSAV(NPAR+1,I) = CHN(C)%CHI2
                 PSAV(NPAR+2,I) = COMPUTE_MAP(CHN(C)%P,CHN(C)%CHI2)
              END DO
           END IF
           IF ((MOD(LONG,WRITFREQ).LT.10).AND.(NMOD.EQ.NMAX)) THEN
              WRITE(SD,*)'Creating an output file...'
              CALL WRITE_DISTRIB(NMOD,FILES(1))
           END IF
        END IF

        END

C     
C -----------------------------------------------------------------------------
C       Changing the directions
C -----------------------------------------------------------------------------
C
        SUBROUTINE CHANGE_DIRS(LONG,CHN,DIRS)

        USE DATA

        IMPLICIT NONE
     
        TYPE(CHAIN), DIMENSION(NCH) :: CHN ! The chains
        INTEGER*4 ::    C,IU,JU,NROT,I      ! Indexes
        INTEGER*8 ::    LONG
        REAl*8 ::       W,B,               ! Means
     &                  COVP(NPAR,NPAR),   ! Averaged covariances
     &                  MOYY(NPAR),   ! Global means over chains
     &                  VP(NPAR),        ! Eigenvalues
     &                  DIRS(NPAR,NPAR) ! Directions
        LOGICAL ::      FAIL            ! True if search beta failed
        
c...  Computing averaged covariances

        DO IU = 1,NPAR
           MOYY(IU) = 0.d0
           DO I = 1,NCH              
              MOYY(IU) = MOYY(IU)+CHN(I)%MOY(IU)
           END DO
           MOYY(IU) = MOYY(IU)/NCH           
        END DO
c        
        DO IU = 1,NPAR
           DO JU = 1,IU
              W = 0.d0
              B = 0.d0
              DO I = 1,NCH              
                 W = W+CHN(I)%COV(IU,JU)
                 B = B+(CHN(I)%MOY(IU)-MOYY(IU))
     &                       *(CHN(I)%MOY(JU)-MOYY(JU))
              END DO
              W = W/NCH
              B = B*LONG/(NCH-1)
              COVP(IU,JU) = (LONG-1_8)*W/LONG+B/LONG
              COVP(JU,IU) = COVP(IU,JU)
           END DO
        END DO
              
c...    Computing eigenvectors of COVP matrix => new directions

        WRITE(SD,*)'******************************************'
        WRITE(SD,*)'Changing directions :'
        CALL JACOBI(COVP,NPAR,NPAR,VP,DIRS,NROT)
        WRITE(SD,'(a)',ADVANCE='NO')'beta ok..'
        DO C = 1,NCH
           CALL SEARCH_BETA(CHN(C)%P,DIRS,CHN(C)%BETA,FAIL)    
           WRITE(SD,'(I3)',ADVANCE='NO')C
        END DO
        WRITE(SD,*)' '
        WRITE(SD,*)'******************************************'

        END
      
C     
C -----------------------------------------------------------------------------
C       1 Step of Metropolis-Hastings algorithm
C -----------------------------------------------------------------------------
C

        SUBROUTINE METRO_STEP(P,DIR,CHI2,BETA,ACC)

        USE DATA

        IMPLICIT NONE
     
        INTEGER*4 ::    IU           ! # of parameter updated
        LOGICAL ::      ACC,         ! Acceptance flag
     &                  BAD          ! True if unreal solution
        REAL*8, DIMENSION(NPAR) :: P, ! Parameters
     &                  DIR,          ! Direction vector       
     &                  PTRY  ! Trial model
        REAl*8 ::       R,           ! Random number
     &                  ALPHA,       ! Acceptance probability                  
     &                  BETA,        ! Scale parameter vector
     &                  RAPQ,        ! q(x|x')/q(x'|q) (eq.42 Ford 06)
     &                  CHI2,        ! Chi^2
     &                  CHI2TRY      ! Trial CHi2    

c...   alpha = min(p(d|orb')/p(d|orb)*q(orb|orb')/q(orb'|orb), 1)     
c...        p(d|orb')|p(d|orb) = Product((j^2+s_i^2)/(j'^2+s_i^2))^(1/2)
c...                              *exp(-(Chi2'-Chi2)/2)
c...            "Product" incorporated into RAPQ         
c...      RAPQ = (Jacobian / Prior) / (Jacobian' / Prior')
c...      alpha = MAP(orbit')/MAP(orbit)*Jacobian/Jacobian'             
c...
        
C        PTRY(1:NPAR) = P(1:NPAR)
        CALL GASDEV(R)
        PTRY(1:NPAR) = P(1:NPAR)+BETA*R*DIR(1:NPAR)
C        DO IU = 1,NPAR
C          PTRY(IU) = PTRY(IU)+BETA*R*DIR(IU)
C        END DO
        BAD = .FALSE.
        
        CALL RATIO_PRIORS(P,PTRY,RAPQ,BAD)
        IF (BAD) THEN
          ALPHA = 0.d0
        ELSE
           CALL CHI2SEUL(PTRY,CHI2TRY)
          ALPHA = MIN(EXP(-0.5d0*(CHI2TRY-CHI2))*RAPQ,1.d0)
        END IF
        CALL RANDOM_NUMBER(R)
        ACC = .FALSE.
c        WRITE(SD,*)sngl(alpha),sngl(chi2try),sngl(chi2),bad,sngl(RAPQ)

        IF (R.LT.ALPHA) THEN
          P(1:NPAR) = PTRY(1:NPAR)
          CHI2 = CHI2TRY
          ACC = .TRUE.
        END IF
c        write(sd,*)sngl(alpha),acc
        END         

C
C -----------------------------------------------------------------------------
C       Routine to fix the scale parameters
C -----------------------------------------------------------------------------
C

        SUBROUTINE SEARCH_BETA(PTRY,DIRS,BETA,FAIL)

        USE DATA

        IMPLICIT NONE
     
        INTEGER*4 ::    I,           ! Index
     &                  IU,          ! # of parameter updated
     &                  LL(NPAR),    ! Permutation on integers
     &                  NS,          ! Number of steps
     &                  NSL(NPAR),   ! # of steps since last update
     &                  NACC(NPAR)   ! # of accepted steps per parameter
        LOGICAL ::      ACC(NPAR),   ! Acceptance flags
     &                  OK(NPAR),    ! Convergence flags per parameter
     &                  FAIL,        ! True if process failed => retry
     &                  SOK          ! Global convergence flag
        REAl*8 ::       P(NPAR),     ! Parameter vector
     &                  PTRY(NPAR),  ! Starting vector
     &                  DIRS(NPAR,NPAR),   ! Directions
     &                  CHI2,        ! CHi2
     &                  BETA(NPAR),  ! Scale vector
     &                  BETAM1(NPAR),BETAM2(NPAR), ! Idem
     &                  S2(NPAR),    ! Frequency parameter
     &                  PSI(NPAR),   ! Acceptance rates
     &                  R,           ! Ratio
     &                  LASTBETA     ! Saving variable

c        CALL PRIOR(N,P)
c        p(1) = 1.d0/22.581947d0
c        p(2) = 0.21907479d0*cos(132.8643*PI/180.d0)
c        p(3) = 0.21907479d0*sin(132.8643*PI/180.d0)
c        p(4) = (T0-2005.9468)*DPI*p(1)+(132.8643*PI/180.d0)
c        WRITE(SD,*)p(1:n)
c        CALL CHI2SEUL(P,CHI2)
c        WRITE(SD,*)'chi2init = ',chi2
c        DIRS = 0.d0
c        DO I=1,N
c          DIRS(I,I) = 1.d0
c        END DO
        FAIL = .FALSE.
        P(1:NPAR) = PTRY(1:NPAR)
c        WRITE(SD,*)p(1:n)        
        CALL CHI2SEUL(P,CHI2)
c        WRITE(SD,*)'chi2init = ',chi2
        NS = 0
        NSL(1:NPAR) = 0
        BETAM1(1:NPAR) = BETA(1:NPAR)
        BETAM2(1:NPAR) = BETA(1:NPAR)
        OK(1:NPAR) = .FALSE.
        S2(1:NPAR)= 2.d0
        SOK = .FALSE.
        PSI(1:NPAR) = 0.d0
        NACC(1:NPAR) = 0
        DO WHILE (.NOT.(SOK.OR.FAIL))
c          CALL RANDOM_NUMBER(R)
c          I = FLOOR(R*FACTN)
          CALL PERMUTATION(NPAR,LL)
          DO I=1,NPAR
            CALL METRO_STEP(P,DIRS(1:NPAR,LL(I)),
     &                       CHI2,BETA(LL(I)),ACC(LL(I)))
          END DO
          NS = NS+1
          DO IU = 1,NPAR
            NSL(IU) = NSL(IU)+1
            IF (ACC(IU)) NACC(IU) = NACC(IU)+1
            PSI(IU) = DBLE(NACC(IU))/DBLE(NS)
            LASTBETA = BETA(IU)
            CALL UPDATE_BETA(NSL(IU),S2(IU),PSI(IU),BETA(IU),OK(IU))
            IF (BETA(IU).NE.LASTBETA) THEN
              BETAM2(IU) = BETAM1(IU)
              BETAM1(IU) = LASTBETA
              R = (BETAM1(IU)-BETAM2(IU))*(BETA(IU)-BETAM1(IU))
              IF (R.LT.0.d0) S2(IU) = S2(IU)+1.d0
            END IF
          END DO
          SOK = .TRUE.
          FAIL = .FALSE.
          DO IU=1,NPAR
            SOK = (SOK.AND.OK(IU))
            FAIL = (FAIL.OR.(ABS(BETA(IU)).GT.BETAMAX))           
          END DO
          IF (FAIL)
     &          WRITE(SD,*)'Beta initiation process failed. Retrying...'
          
c     if (mod(ns,1).eq.0) then
c            WRITE(SD,*)ns,'beta',sngl(beta(1:npar)),ok(1:npar)
c            WRITE(SD,*)'psi',sngl(psi(1:npar))
c            print*,HUGE(r),beta(npar)
c            if (beta(npar).gt.HUGE(r)) stop
c         end if
        END DO
c        WRITE(SD,*)ns,'beta',sngl(beta(1:npar)),ok(1:npar)
c        WRITE(SD,*)'psi',sngl(psi(1:npar))
c        WRITE(SD,*)'chi2',chi2
c        WRITE(SD,*)'p',sngl(p(1:npar)),1.d0/p(1)
        END


C
C -----------------------------------------------------------------------------
C       Routine to update the beta's
C -----------------------------------------------------------------------------
C

        SUBROUTINE UPDATE_BETA(NSL,S2,PSI,BETA,OK)

        USE DATA

        IMPLICIT NONE
     
        REAL*8, PARAMETER :: PSI0 = 0.44d0
        INTEGER*4 ::    NSL          ! NUmber of steps per parameter
        LOGICAL ::      OK           ! Validity test
        REAl*8 ::       PSI,         ! Acceptance rate
     &                  BETA,        ! Scale parameter
     &                  RHO,         ! Rapport
     &                  S2           ! Frequency parameter

        RHO = PSI/PSI0
        IF ((PSI-PSI0)**2.GT.S2*PSI0*(1.d0-PSI0)/NSL) THEN
          NSL = 0
          IF (PSI.GT.0.5d0*PSI0) THEN
            BETA = BETA*RHO
          ELSE IF (PSI.GT.0.2*PSI0) THEN
            BETA = BETA*RHO*SQRT(RHO)
          ELSE IF (PSI.GT.0.1*PSI0) THEN
            BETA = BETA*RHO*RHO
          ELSE
            BETA = BETA*0.01d0
          END IF
        END IF
c        WRITE(SD,*)beta,PSI,PSI0,RHO
        OK = (ABS(RHO-1.d0).LT.1d-1)

        END


C -----------------------------------------------------------------------------
C       Taking info account user supplied priors on masses in MAP (log)
C -----------------------------------------------------------------------------
C

        SUBROUTINE USER_PRIOR_MAP(IV,FMAP,EX)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    IV              ! # of variable
        REAL*8 ::       FMAP,           ! MAP to modify
     &                  X,              ! ln(mass)
     &                  EX,             ! mass,
     &                  XBAR,SIGMA,    ! Mean and standard deviation
     &                  S               ! Auxiliary variable

        SELECT CASE (MPRIOR(IV)%TYP)
           CASE (0)               ! Logarithmic prior (default)
             X = LOG(EX)
             FMAP = FMAP-X     
           CASE (1)                ! Normal (Gaussian) prior
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = MPRIOR(IV)%SDEV
             S = (EX-XBAR)/SIGMA
             FMAP = FMAP-0.5d0*S*S
           CASE (2)           ! Log-normal prior                    
             X = LOG(EX)
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = MPRIOR(IV)%SDEV
             S = (X-XBAR)/SIGMA
             FMAP = FMAP-0.5d0*S*S-X
           CASE (4)    ! Fixed value (treated as steep gaussian 0.1%)
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = 1d-3*XBAR
             S = (EX-XBAR)/SIGMA
             FMAP = FMAP-0.5d0*S*S
           CASE (6)    ! Heaviside
             EX = EXP(X)
             SIGMA = MPRIOR(IV)%SDEV
             XBAR = MPRIOR(IV)%MEAN
             S = (EX-XBAR)/SIGMA
             FMAP = FMAP-LOG(1.d0+EXP(-S))
          END SELECT
        END
C     
C -----------------------------------------------------------------------------
C       Taking info account user supplied priors on masses in prior ratio
C -----------------------------------------------------------------------------
C

        SUBROUTINE USER_PRIOR_RATIO(IV,RAPQ,EX,EXTRY)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    IV              ! # of variable
        REAL*8 ::       RAPQ,           ! Prior ratio
     &                  X,XTRY,         ! ln(mass) & ln(mass_try)
     &                  EX,EXTRY,       ! masses,
     &                  XBAR,SIGMA,      ! Mean and standard deviation
     &                  S,ST            ! Auxiliary variables

        SELECT CASE (MPRIOR(IV)%TYP)
           CASE(0)
              X = LOG(EX)
              XTRY = LOG(EXTRY)
              EXTRY = EXP(XTRY)
              RAPQ = X/XTRY ! log prior (default)
           CASE (1)                ! Normal (Gaussian) prior
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = MPRIOR(IV)%SDEV
             S = (EX-XBAR)/SIGMA
             ST = (EXTRY-XBAR)/SIGMA
             RAPQ = EXP(0.5d0*(S*S-ST*ST))
           CASE (2)           ! Log-normal prior                    
             X = LOG(EX)
             XTRY = LOG(EXTRY)
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = MPRIOR(IV)%SDEV
             S = (X-XBAR)/SIGMA
             ST = (XTRY-XBAR)/SIGMA
             RAPQ = EXP(0.5d0*(S*S-ST*ST))*EX/EXTRY
           CASE(3)      ! Linear
              RAPQ = 1.d0
           CASE(4)               ! Fixed
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = XBAR*1d-6
             S = (EX-XBAR)/SIGMA
             ST = (EXTRY-XBAR)/SIGMA
             RAPQ = EXP(0.5d0*(S*S-ST*ST))
           CASE(6)    ! Heaviside
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = MPRIOR(IV)%SDEV
             S = (EX-XBAR)/SIGMA
             ST = (EXTRY-XBAR)/SIGMA
             IF (ST.LT.-1d3) THEN
                RAPQ = 0.d0
             ELSE                
                RAPQ = (1.d0+EXP(-S))/(1.d0+EXP(-ST)) 
             END IF
c             print*,'xb',iv,sngl(extry),
c     &             sngl(s),sngl(st),sngl(rapq)
c             if (xtry.ne.xtry) stop
          END SELECT
        END
              
C
C -----------------------------------------------------------------------------
C     Taking info account user supplied priors on masses into FMAP
C     FMAP = -2*ln(prior) (F(X), dF/DX,d2f/DX^2)
C -----------------------------------------------------------------------------
C

        SUBROUTINE FMAP_PRIOR(IV,EX,F,DF,D2F)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    IV
        REAL*8 ::       F,              ! Function to add into Map
     &                  X,EX,           ! ln(mass) & mass
     &                  DF,D2F,         ! dF/dX, d2F/dX^2
     &                  XBAR,SIGMA,     ! Mean and standard deviation
     &                  S,ES,FF         ! Auxiliary quantities

        SELECT CASE (MPRIOR(IV)%TYP)
          CASE (0)  ! Logarithmic prior (default) (~ 1/X)
             X = LOG(EX)
             F = 2.d0*X
             DF = 2.d0
             D2F = 0.
          CASE (1) ! Gaussian (Normal) prior  (~ exp(X^2/2))
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = MPRIOR(IV)%SDEV
             S = (EX-XBAR)/SIGMA
             F = S*S
             DF = 2.d0*S*EX/SIGMA
             D2F = 2.d0*EX*(2.d0*EX-XBAR)/(SIGMA*SIGMA)
          CASE (2)              ! Log-normal prior (~ exp(-lnX^2/2)/X)
             X = LOG(EX)
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = MPRIOR(IV)%SDEV
             S = (X-XBAR)/SIGMA
             F = S*S+2.d0*X
             DF = 2.d0*S/SIGMA+2.d0
             D2F = 2.d0/(SIGMA*SIGMA)
          CASE (3)               ! Linear prior (constant)
             F = 0.d0
             DF = 0.d0
             D2F = 0.d0
          CASE (4)      ! Fixed value (treated as steep gaussian 1d-6)
             XBAR = MPRIOR(IV)%MEAN
             SIGMA = 1.d-6*XBAR
             S = (EX-XBAR)/SIGMA
             F = S*S
             DF = 2.d0*S*EX/SIGMA
             D2F = 2.d0*EX*(2.d0*EX-XBAR)/(SIGMA*SIGMA)
          CASE(6)       ! Heaviside (1/(1+exp(-k*X)))
             SIGMA = MPRIOR(IV)%SDEV
             XBAR = MPRIOR(IV)%MEAN
             S = (EX-XBAR)/SIGMA
             ES = EXP(-S)
             F = 2.d0*LOG(1.d0+ES)
             FF = (1.d0+ES)*SIGMA
             DF = -2.d0*EX*ES/FF
             D2F = 2.d0*(EX-FF)*EX*ES/(FF*FF)
          END SELECT
        END

      
C
C -----------------------------------------------------------------------------
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
        REAl*8 ::       RAPQ,   ! q(x|x')/q(x'|x) (eq.42 Ford 06) (output)
     &                  R,           ! Random number
     &                  S2,S,        ! Intermediaires
     &                  P1,P2,P3,P4,P5,P6,P7, ! Intermédiaires
     &                  PT1,PT2,PT3,PT4,PT5,PT6,PT7, !
     &                  EXC,EXQ,EXCT,EXQT, ! excentricite
     &                  CW,SW,CWT,SWT,  ! w = Omega+omega
     &                  CL,SL,          ! Lambda = w+v = Omega+omega+v
     &                  SI2,         ! sin(i/2)^2
     &                  CV,SV,          ! V = anomalie vraie
     &                  UP,UPT,      ! Inverses de périodes
     &                  SIGJV,SIGJVT, ! Velocity Jitters
     &                  A, AT,       ! Semi-major axes
     &                  MUF,MUFT,LMUF,LMUFT,  ! Fractional masses & ln
     &                  MSTAR, MSTART,        ! Stellar mass 
     &                  FAC,FACT              ! Jacobiens

c     For each planet,
c  Jabobian J (orbite->params) = 1/2*e*sin(i)/((1-e^2)^3/2*P^3*m)*dnu/dM 
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
          BAD = BAD.OR.(EXQT.LE.0.d0).OR.
     &             (EXCT.GT.VPRIOR(3)%BOUND(2)).OR.
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
       
        IF (RADVEL.AND.(JITNUM.EQ.1)) THEN
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

c...  Contribution of RV's (type 2 data) to Chi2.
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
              SI(I) = 2.d0*SII*SQRT(1.d0-SI2(I))       ! sin(i)
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
           CL = PER*P6       ! lambda = Omega+omega+M0 = w+M0
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
C    Reads a radial velocity data file
C-----------------------------------------------------------------------------
C

        SUBROUTINE READ_RV(IDAT,FILNAM)

        USE DATA
        
        IMPLICIT NONE

        CHARACTER*(*) :: FILNAM ! File name  
        INTEGER*4 :: IDAT       ! Type of data
        INTEGER*4 :: NDAT       ! Current data
        INTEGER*4 ::   MM,YR    ! Month, Year
        REAL*8 ::      JJ,JD    ! Day, date
        REAL*8 ::      VV,DVV   ! Radial velocity
        REAL*8 ::      DATE_TO_JD ! Conversion function
        INTEGER*4 ::   I,ERROR
        LOGICAL ::     OK
        CHARACTER*80 :: LIG 
        
        OPEN(15,FILE=FILNAM,STATUS='UNKNOWN')                                   
              
        READ(15,'(a)')LIG
        READ(15,'(a)')LIG
        IF (DATATYP(IDAT).EQ.2) READ(15,*)STAR%JDOFFSET
        OK = .TRUE.    
        STAR%NDATVR = 0
        DO WHILE(OK)
           SELECT CASE (DATATYP(IDAT))
           CASE(1)
             READ(15,*,IOSTAT=ERROR)JJ,MM,YR,VV,DVV           
           CASE(2)
             READ(15,*,IOSTAT=ERROR)JD,VV,DVV    
           END SELECT
           OK=(ERROR.EQ.0)      
           IF (OK) STAR%NDATVR = STAR%NDATVR+1
        END DO
        NDAT = STAR%NDATVR
        ALLOCATE(STAR%TVR(NDAT))
        ALLOCATE(STAR%V(NDAT))
        ALLOCATE(STAR%SIGV(NDAT))
        ALLOCATE(STAR%SIGVM2(NDAT))
        ALLOCATE(STAR%SIGV2(NDAT))
        REWIND(15)
        READ(15,'(a)')LIG                                                       
        READ(15,'(a)')LIG
        IF (DATATYP(IDAT).EQ.2) READ(15,*)STAR%JDOFFSET
        STAR%NDATVR = 0
        OK = .TRUE.
        DO WHILE(OK)
           SELECT CASE (DATATYP(IDAT))
           CASE(1)
             READ(15,*,IOSTAT=ERROR)JJ,MM,YR,VV,DVV           
           CASE(2)
             READ(15,*,IOSTAT=ERROR)JD,VV,DVV    
           END SELECT
           OK=(ERROR.EQ.0)      
           IF (OK) THEN
              STAR%NDATVR = STAR%NDATVR+1
              NDAT = STAR%NDATVR
c...  Conversion JJ-MM-YR -> JD
              SELECT CASE (DATATYP(IDAT))
              CASE(1)
                 JD = DATE_TO_JD(JJ,MM,YR)-STAR%JDOFFSET 
              CASE(2)
                 JD = JD   !   +STAR%JDOFFSET
              END SELECT
              STAR%TVR(NDAT) = JD
              STAR%V(NDAT) = VV*MPS ! Conversion to AU/day
              STAR%SIGV(NDAT) = DVV*MPS 
              STAR%SIGV2(NDAT) = STAR%SIGV(NDAT)*STAR%SIGV(NDAT)
              STAR%SIGVM2(NDAT) = 1.d0/STAR%SIGV2(NDAT)
           END IF
        END DO
        CLOSE(15)

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

