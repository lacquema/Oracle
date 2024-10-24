c 10h30 v3/3
        MODULE DATA
        
        PUBLIC
        
        INCLUDE '../sub/commons.f'

        END MODULE DATA

c-------------------------------------------------------------------------
       
        PROGRAM ASTROM_MCMCO

        USE DATA

        IMPLICIT NONE

        INTEGER*4        I,DKAL,mm,yy,
     &       METHOD / 1 / ,     ! Identification of the code 
     &       NMOD               ! Nombre de modeles
        REAL*8, DIMENSION(:), ALLOCATABLE :: P ! Parameters
        REAL*8, DIMENSION(:,:), ALLOCATABLE :: COV ! Covariance matrix
	REAL*8 ::
     &       CHI2,FMAP,         ! CHI2
     &       UP,                ! 1/period
     &       COMPUTE_MAP,       ! Function to compute MAP
     &       SI2,               ! sin(i/2)
     &       CL,SL              ! Lambda = Omega+omega+v,

        CHARACTER*80 ::  LIG
        CHARACTER*80, DIMENSION(NFIL) :: FILES ! #1 = output, #2 = dump,
c                                                #3 = data
        INTEGER*4 ::         CNEW
        LOGICAL ::          OK,NEW

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
                 WRITE(SD,*)'a',sngl(pla(1:npla)%a)
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
              PSTART(DKAL+2) = PLA(I)%EXC/PLA(I)%EXQ*PLA(I)%CW
c                                            Par 2 = e*cos(w)/sqrt(1-e^2)
              PSTART(DKAL+3) = PLA(I)%EXC/PLA(I)%EXQ*PLA(I)%SW
c                                            Par 3 = e*sin(w)/sqrt(1-e^2)
              IF (RADVEL) THEN
                 PSTART(DKAL+4) = SQRT(PLA(I)%EXQ)*SI2*PLA(I)%CO
c                                     Par 4 = sin(i/2)*cos(O)*(1-e^2)^(1/4)
                 PSTART(DKAL+5) = SQRT(PLA(I)%EXQ)*SI2*PLA(I)%SO
c                                     Par 5 = sin(i/2)*sin(O)*(1-e^2)^(1/4)
              ELSE
                 PSTART(DKAL+4) = SQRT(PLA(I)%EXQ)*SI2*PLA(I)%CP
c                                     Par 4 = sin(i/2)*cos(phi)*(1-e^2)^(1/4)
                 PSTART(DKAL+5) = SQRT(PLA(I)%EXQ)*SI2*PLA(I)%SP
c                                     Par 5 = sin(i/2)*sin(phi)*(1-e^2)^(1/4)
              END IF
              CALL TP_TO_LAMBDA(STAR%T0,PLA(I)%TP,DPI*UP,PLA(I)%EXC,
     &              PLA(I)%EXQ,PLA(I)%CW,PLA(I)%SW,CL,SL)              
              PSTART(DKAL+6) = UP*CL    ! Par 6 = cos(lambda)/P
              PSTART(DKAL+7) = UP*SL    ! Par 7 = sin(lambda)/P
              PSTART(DKAL+1) = LOG(PLA(I)%A)   ! Par 1 = ln(a)
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
           WRITE(SD,*)'a',sngl(pla(1:npla)%a)
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
          CALL DISPLAY_SOLUTION(METHOD)
        END IF
       
        CALL UTIL_EXIT()
        END
