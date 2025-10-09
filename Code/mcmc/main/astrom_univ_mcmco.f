c 10h30 v3/3
        MODULE DATA
        
        PUBLIC
        
        INCLUDE '../sub/commons.f'

        END MODULE DATA

c-------------------------------------------------------------------------
       
        PROGRAM ASTROM_UNIV_MCMCO

        USE DATA

        IMPLICIT NONE

        INTEGER*4        I,DKAL,mm,yy,
     &       NMOD               ! Nombre de modeles
        REAL*8, DIMENSION(:), ALLOCATABLE :: P ! Parameters
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: COV ! Covariance matrix
	REAL*8 ::
     &       CHI2,FMAP,         ! CHI2
     &       UP,                ! 1/period
     &       MASS,              ! Dynamical mass
     &       COMPUTE_MAP,       ! Function to compute MAP
     &       SI2,               ! sin(i/2)
     &       M0,                ! Anomalie moyenne de reference
     &       HH2,               ! -2xEnergie = alpha de kepu
     &       S,                 ! niversal variable
     &       C0,SC1,S2C2,S3C3   ! Fonctions de Stumpff

        CHARACTER*80 ::  LIG
        CHARACTER*80, DIMENSION(NFIL) :: FILES ! #1 = output, #2 = dump,
c                                                #3 = data
        INTEGER*4 ::         CNEW
        LOGICAL ::          OK,NEW

        METHOD = 2
        CALL INIT_SIMU(CNEW,NMOD,FILES,NEW)
        
        IF (CNEW.EQ.5) THEN
           CALL SIMUL_DATA(1,FILES)
           CALL UTIL_EXIT()
           STOP
        END IF       

        ALLOCATE(P(NPAR))
        ALLOCATE(COV(NPAR,NPAR))  
       
        CALL INIT_DATA(NPAR,P,COV,FILES,NEW)

        IF (NEW) THEN           
           COV = 0.d0
           IF ((CNEW.EQ.2).OR.(CNEW.EQ.4)) THEN
              IF (ISDATA(2).AND.(JITNUM.EQ.1)) THEN 
c...  If RV & jitter, first perform a L-M fit with no Jitter
                 CALL MRQMIN(NPAR-1,P,COV,CHI2,FMAP)
                 WRITE(SD,*)'Chi2+Fmap =',SNGL(CHI2),
     &             'Chi2 =',SNGL(CHI2-FMAP),
     &             'Fmap =',SNGL(FMAP)         
                 CALL MRQELEMENTS(P,COV)
                 WRITE(SD,*)'q',sngl(pla(1:npla)%q)
                 WRITE(SD,*)'exc',sngl(pla(1:npla)%exc)
                 WRITE(SD,*)'i',sngl(pla(1:npla)%inc)
                 WRITE(SD,*)'tp',sngl(pla(1:npla)%tp)
c...  Then use the result of the first L-M fit as starting point for
c..   a new fit including jitter                 
                 P(NPAR) = P(NPAR-1)
                 P(NEL*NPLA+2) = LOG(STAR%SIGJV)
              END IF
              CALL MRQMIN(NPAR,P,COV,CHI2,FMAP)
              WRITE(SD,*)'Chi2+Fmap =',SNGL(CHI2),
     &             'Chi2 =',SNGL(CHI2-FMAP),
     &             'Fmap =',SNGL(FMAP)         
           END IF

           CALL MRQELEMENTS(P,COV)

           DO I = 1,NPLA
              DKAL = NEL*(I-1)
              UP = 1.d0/PLA(I)%PER
              SI2 = PLA(I)%TI/SQRT(1.d0+PLA(I)%TI*PLA(I)%TI) ! sin(i/2) 
              PSTART(DKAL+2) = PLA(I)%EXC*PLA(I)%CW ! Par 2 = e*cos(w)
              PSTART(DKAL+3) = PLA(I)%EXC*PLA(I)%SW ! Par 3 = e*sin(w)
              IF (RADVEL) THEN
                 PSTART(DKAL+4) = SI2*PLA(I)%CO ! Par 4 = sin(i/2)*cos(O)
                 PSTART(DKAL+5) = SI2*PLA(I)%SO ! Par 5 = sin(i/2)*sin(O)
              ELSE
                 PSTART(DKAL+4) = SI2*PLA(I)%CP ! Par 4 = sin(i/2)*cos(phi)
                 PSTART(DKAL+5) = SI2*PLA(I)%SP ! Par 5 = sin(i/2)*sin(phi)
              END IF

c              MASS = (DPI*UP)**2*PLA(I)%Q**3
c              HH2 = -(PLA(I)%EXC-1.d0)*MASS/PLA(I)%Q
c...                                hh2 = -2xEnergie = alpha pour kepu          
c              M0 =  STAR%T0-P(DKAL+6)
c              CALL KEP_UNIV(M0,PLA(I)%Q,MASS,HH2,S,C0,SC1,S2C2,S3C3)
              PSTART(DKAL+1) = LOG(PLA(I)%Q)           ! Par 1 = ln(q)
              PSTART(DKAL+6) = LOG(PLA(I)%PER)         ! Par 6 = ln(perq)
              PSTART(DKAL+7) = P(DKAL+7) ! Par 7=f(E)*S=f(E)*univ.var@t=T0    
           END DO
           IF (MULTIPLA) PSTART(NPAR) = LOG(STAR%MASS) ! Stellar mass
           IF (ISDATA(2)) THEN
              PSTART(NEL*NPLA+1) = STAR%V0 !  Offset velocity
              IF (JITNUM.EQ.1) PSTART(NEL*NPLA+2) = LOG(STAR%SIGJV) ! Jitter
           END IF
           
           CALL CHI2SEUL(PSTART(1:NPAR),CHI2)
c...  Store CHi2 and MAP in additional variables
           PSTART(NPAR+1) = CHI2
           PSTART(NPAR+2) = COMPUTE_MAP(PSTART(1:NPAR),CHI2)

           WRITE(SD,*)'chi2init = ',chi2,'map = ', pstart(npar+2)
           WRITE(SD,*)'q',sngl(pla(1:npla)%q)
           WRITE(SD,*)'exc',sngl(pla(1:npla)%exc)
           WRITE(SD,*)'i',sngl(pla(1:npla)%inc)
           WRITE(SD,*)'tp',sngl(pla(1:npla)%tp)
           IF (ISDATA(2).AND.(JITNUM.EQ.1))
     &                WRITE(SD,*)'jitter',sngl(star%sigjv/mps)
           WRITE(SD,*)'M*',sngl(STAR%mass/smas)
        END IF

        IF (CNEW.LT.4) THEN
          CALL MCMCP(NMOD,NEW,FILES)
          IF (RADVEL) THEN
             CALL WRITE_DISTRIB(NMOD,FILES(1))
          ELSE
             CALL WRITE_DISTRIB_DOUBLING(NMOD,FILES(1))
          END IF          
        ELSE
          CALL DISPLAY_SOLUTION()
        END IF
       
        CALL UTIL_EXIT()
        END
