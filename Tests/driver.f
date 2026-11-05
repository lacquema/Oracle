c-------------------------------------------------------------------------
c     Programme de test : evalue les routines de Code/mcmc a un jeu de
c     parametres fixe, sans hasard ni MCMC. Il lit la meme entree que
c     ASTROM_MCMCO et n'ecrit ses resultats que sur les lignes "TEST".
c     Argument : 1 = variables classiques, 2 = variables universelles.
c-------------------------------------------------------------------------
        MODULE DATA

        PUBLIC

        INCLUDE 'commons.f'

        END MODULE DATA

c-------------------------------------------------------------------------

        PROGRAM DRIVER

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    I,              ! Index
     &                  DKAL,           ! Offset of an orbit in P
     &                  L,              ! Index of a displaced point
     &                  NMOD,           ! Nombre de modeles
     &                  CNEW            ! Choice
        REAL*8, DIMENSION(:), ALLOCATABLE ::
     &                  P,              ! Parameters
     &                  PM,             ! Parameters of the MCMC
     &                  PB,             ! Displaced parameters
     &                  BETA,           ! Vector of L-M
     &                  BETAB           ! Dummy vector
        REAL*8, DIMENSION(:,:), ALLOCATABLE ::
     &                  COV,            ! Covariance matrix
     &                  ALPHA           ! Matrix of L-M
        REAL*8 ::       CHI2,           ! Chi2 from CHI2SEUL
     &                  MERIT,          ! Chi2+Fmap from MRQCOF
     &                  FMAP,           ! Fmap from MRQCOF
     &                  MPLUS,MMOINS,   ! Merit function at P+H, P-H
     &                  DUMMY,          ! Unused output
     &                  CL,SL,          ! cos,sin(longitude at T0)
     &                  RAPQ,           ! Prior ratio of the MCMC
     &                  CHIT,           ! Chi2 of a displaced point
     &                  H,              ! Finite difference step
     &                  COMPUTE_MAP     ! Function to compute MAP
        CHARACTER*80, DIMENSION(NFIL) :: FILES ! Output, dump, data
        CHARACTER*8 ::  ARG             ! Command line argument
        LOGICAL ::      NEW,            ! New simulation
     &                  BAD             ! Unreal solution

        CALL GET_COMMAND_ARGUMENT(1,ARG)
        READ(ARG,*)METHOD

        CALL INIT_SIMU(CNEW,NMOD,FILES,NEW)
        ALLOCATE(P(NPAR))
        ALLOCATE(PM(NPAR))
        ALLOCATE(PB(NPAR))
        ALLOCATE(BETA(NPAR))
        ALLOCATE(BETAB(NPAR))
        ALLOCATE(COV(NPAR,NPAR))
        ALLOCATE(ALPHA(NPAR,NPAR))
        CALL INIT_DATA(NPAR,P,COV,FILES,NEW)

c...  Same shift as ASTROM_MCMCO when there is a jitter
        IF (ISDATA(2).AND.(JITNUM.EQ.1)) THEN
           P(NPAR) = P(NPAR-1)
           P(NEL*NPLA+2) = LOG(STAR%SIGJV)
        END IF

 1      FORMAT('TEST ',A,1X,I4,1X,ES24.16)
 2      FORMAT('TEST ',A,1X,ES24.16)
        WRITE(SD,'(A,3(1X,I4),1X,L1)')
     &       'TEST META',NPAR,NPLA,NFREE,MULTIPLA

c...  Levenberg-Marquardt and the MCMC do not use the same parameters :
c...  L-M    : Par 4,5 = tan(i/2)*(cos,sin)
c...  Classical  MCMC : Par 4,5 = (1-e^2)^(1/4)*sin(i/2)*(cos,sin)
c...  Universal  MCMC : Par 4,5 = sin(i/2)*(cos,sin)
c...  Same conversion as ASTROM_MCMCO, without going through MRQELEMENTS
        PM = P
        DO I = 1,NPLA
           DKAL = NEL*(I-1)
           PM(DKAL+4:DKAL+5) = P(DKAL+4:DKAL+5)
     &                         /SQRT(1.d0+PLA(I)%TI*PLA(I)%TI)
           IF (METHOD.EQ.1) THEN
              PM(DKAL+4:DKAL+5) = PM(DKAL+4:DKAL+5)*SQRT(PLA(I)%EXQ)
c...  L-M : mean longitude / MCMC : true longitude at T0
              CALL LONGITUDE_T0(I,CL,SL)
              PM(DKAL+6) = CL/PLA(I)%PER
              PM(DKAL+7) = SL/PLA(I)%PER
           END IF
        END DO

c...  Chi2, MAP and merit function along the two code paths
        CALL CHI2SEUL(PM,CHI2)
        WRITE(SD,2)'CHI2',CHI2
        WRITE(SD,2)'MAP',COMPUTE_MAP(PM,CHI2)
        CALL MRQCOF(NPAR,P,ALPHA,BETA,MERIT,FMAP)
        WRITE(SD,2)'MERIT',MERIT
        WRITE(SD,2)'FMAP',FMAP

c...  Prior ratio of the MCMC between the point PM and displaced points
        DO I = 1,NPAR
           WRITE(SD,1)'PM',I,PM(I)
        END DO
        WRITE(SD,2)'LNPRIOR0',COMPUTE_MAP(PM,CHI2)+0.5d0*CHI2
        DO L = 1,4
           PB = PM
           DO I = 1,NPAR
              PB(I) = PM(I)+SIN(1.7d0*L*I)*MAX(5.d-4*ABS(PM(I)),1.d-7)
           END DO
           CALL RATIO_PRIORS(PM,PB,RAPQ,BAD)
           CALL CHI2SEUL(PB,CHIT)
           DO I = 1,NPAR
              WRITE(SD,1)'PT',1000*L+I,PB(I)
           END DO
           WRITE(SD,1)'RAPQ',L,RAPQ
           WRITE(SD,1)'LNPRIOR',L,COMPUTE_MAP(PB,CHIT)+0.5d0*CHIT
        END DO

c...  Gradient of Levenberg-Marquardt against finite differences :
c...  BETA = -1/2 d(Chi2+Fmap)/dP
        DO I = 1,NPAR
           H = 1.d-6*MAX(ABS(P(I)),1.d-2)
           PB = P
           PB(I) = P(I)+H
           CALL MRQCOF(NPAR,PB,ALPHA,BETAB,MPLUS,DUMMY)
           PB(I) = P(I)-H
           CALL MRQCOF(NPAR,PB,ALPHA,BETAB,MMOINS,DUMMY)
           WRITE(SD,1)'BETA',I,BETA(I)
           WRITE(SD,1)'BETAFD',I,-0.25d0*(MPLUS-MMOINS)/H
           WRITE(SD,1)'P',I,P(I)
        END DO

c...  Elements rebuilt from the parameters
        COV = 0.d0
        DO I = 1,NPAR
           COV(I,I) = (1.d-3*I)**2
        END DO
        CALL MRQELEMENTS(P,COV)
        WRITE(SD,2)'MSTAR',STAR%MASS
        WRITE(SD,2)'DMSTAR',STAR%DMASS
        WRITE(SD,2)'COVLAST',SQRT(COV(NPAR,NPAR))
        WRITE(SD,2)'PLAST',P(NPAR)
        DO I = 1,NPLA
           WRITE(SD,1)'MU',I,PLA(I)%MU
           WRITE(SD,1)'MDYN',I,PLA(I)%MDYN
        END DO

        CLOSE(SDL)
        END
