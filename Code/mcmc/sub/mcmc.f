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
        INTEGER*8 ::    LONG,
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
c                 write(sd,*)'ps ',sngl(chn(c)%p)
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
        IF (VERBOSE) THEN
           WRITE(SD,*)'nmod = ',nmod,
     &          'long = ',long,' chi2 = ',sngl(chn(1:nch)%chi2)
c           WRITE(SD,*)'exc 1',
c     &          (sngl(sqrt(chn(i)%p(2)**2+chn(i)%p(3)**2)),i=1,nch)
c           WRITE(SD,*)'exc 2',
c     &      (sngl(sqrt(chn(i)%p(9)**2+chn(i)%p(10)**2)),i=1,nch) 
        end if
        
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
              IF (RADVEL) THEN
                 CALL WRITE_DISTRIB(NMOD,FILES(1))
              ELSE
                 CALL WRITE_DISTRIB_DOUBLING(NMOD,FILES(1))
              END IF          
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
        
C     PTRY(1:NPAR) = P(1:NPAR)
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
c        WRITE(SD,*)'p = ',sngl(p(1:npar))        
        CALL CHI2SEUL(P,CHI2)
c        WRITE(SD,*)'chi2init = ',sngl(chi2)
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
          
c        if (mod(ns,1).eq.0) then
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
c        print*,beta,PSI,PSI0,RHO
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

      
