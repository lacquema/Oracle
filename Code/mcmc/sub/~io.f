C
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
        WRITE(SD,'(a)')' m unit a/q e i Omega omega tp, where'
        WRITE(SD,'(a)')'   m = Mass'
        WRITE(SD,'(a)')'   unit = ms ou mj (Solar or Jupiter mass)'
        WRITE(SD,'(a)')'   a/q = Semi major axis / Periastron (au)'
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
           IF (UNIT.EQ.'MS') PLA(I)%MUNIT = SMAS
           IF (UNIT.EQ.'MJ') PLA(I)%MUNIT = MJUP
           PLA(I)%MU = PLA(I)%MU*PLA(I)%MUNIT
           SIGMA = SIGMA+PLA(I)%MU
           PLA(I)%MDYN = SIGMA
           SELECT CASE (METHOD)
           CASE(1) ! Elliptic variables
              PLA(I)%Q = PLA(I)%A*(1.d0-PLA(I)%EXC)
           CASE(2)              ! Universal variables
              PLA(I)%Q = PLA(I)%A
           END SELECT
           NN = SQRT(SIGMA/PLA(I)%A**3)
           PLA(I)%PER = DPI/NN
           IF (PLA(I)%EXC.LT.1.d0) PLA(I)%EXQ = SQRT(1.d0-PLA(I)%EXC**2)
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
          MPRIOR(IV)%MUNIT = FAC
          MPRIOR(IV)%MEAN = MPRIOR(IV)%MEAN*FAC ! Convert into the right unit
          MPRIOR(IV)%SDEV = MPRIOR(IV)%SDEV*FAC
        CASE(2)                  ! Log-normal
          WRITE(SD,*)
     &           'Enter mean and standard deviation + unit (ms or mj) '
          READ(5,'(a)')MPRIOR(IV)%MEAN,MPRIOR(IV)%SDEV,UNIT
          CALL CONVERT_UNIT(UNIT,FAC)
          MPRIOR(IV)%MUNIT = FAC
          MPRIOR(IV)%SDEV = MPRIOR(IV)%SDEV/MPRIOR(IV)%MEAN ! ln(1±s/x)~s/x
          MPRIOR(IV)%MEAN = LOG(MPRIOR(IV)%MEAN*FAC)
        CASE(3) ! Linear
          WRITE(SD,*)'Enter bounds + unit (ms or mj) '
          READ(5,'(a)')MPRIOR(IV)%BOUND(1:2),UNIT
          CALL CONVERT_UNIT(UNIT,FAC)
          MPRIOR(IV)%MUNIT = FAC
          MPRIOR(IV)%BOUND(1:2) = MPRIOR(IV)%BOUND(1:2)*FAC
        CASE(4) ! Fixed
          WRITE(SD,*)'Enter value + unit (ms or mj) '
          READ(5,'(a)')MPRIOR(IV)%MEAN,UNIT
          CALL CONVERT_UNIT(UNIT,FAC)
          MPRIOR(IV)%MUNIT = FAC
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

        SUBROUTINE INIT_SIMU(CNEW,NMOD,FILES,NEW)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    I,J,CNEW,ERROR,NMOD
        INTEGER*4, DIMENSION(4) :: ISDATAN
        LOGICAL ::      OK,NEW
        CHARACTER*3 ::  UNIT
        CHARACTER*(*), DIMENSION(NFIL) :: FILES ! #1 = output, #2 = dump,
c                                                #3 = data

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
           MULTIPLA = (MULTIPLA.OR.ISDATA(2).OR.ISDATA(4))
c...  MULTIPLA = .TRUE. whenever individual masses will have to be fitted.
c...  This is the case if more than one orbit is considered (NPLA>1)
c...  or when information on stellar proper motion is present (Type 2 & 4 data)
           ALLOCATE(PLA(NPLA))
           WRITE(SD,*)'Name of data file ?'
           READ(5,'(a)')FILES(3)
           WRITE(SD,*)'Precision ?'
           READ(5,*)EPS
           WRITE(SD,*)'Choose data format options :'
           WRITE(SD,*)
     &'DDMMYR/JD (1/2), (DEC,RA)/(SEP,PA) (1/2), corr. coefs (1/0)' 
           READ(5,*)DATEINPUT,XYFORMAT,CORRNUM
           CORR = (CORRNUM.EQ.1)
           DATATYP(1) = 4*(XYFORMAT-1)+2*(DATEINPUT-1)+2-CORRNUM
           DATATYP(2) = DATEINPUT
           DATATYP(3) = DATEINPUT
           DATATYP(4) = DATATYP(1)
           IF (ISDATA(2)) THEN
              WRITE(SD,*)'Should we consider a jitter (1=yes, 0=no) ?'
              READ(5,*)JITNUM
           END IF                 

c...  Recall chosen input formats :           
c...  Case of Relative astrometry (HCI)
           WRITE(SD,*)'Chosen HCI data (type 1 data) line format :'
           SELECT CASE (DATATYP(1))
           CASE(1)
              WRITE(SD,*)
     &'#Planet Day Month Year Dec(mas) RA(mas) dDec(mas) dRA(mas) Corr'
           CASE(2)
              WRITE(SD,*)
     &'#Planet Day Month Year Dec(mas) RA(mas) dDec(mas) dRA(mas)'           
           CASE(3)
              WRITE(SD,*)
     &'#Planet JD(-offset) Dec(mas) RA(mas) dDec(mas) dRA(mas) Corr'
           CASE(4)
              WRITE(SD,*)
     &'#Planet JD(-offset) Dec(mas) RA(mas) dDec(mas) dRA(mas)'
           CASE(5)
              WRITE(SD,*)
     &'#Planet Day Month Year Sep(mas) PA(deg) dSep(mas) dPA(deg) Corr'
           CASE(6)
              WRITE(SD,*)
     &'#Planet Day Month Year Sep(mas) PA(deg) dSep(mas) dPA(deg)'           
           CASE(7)
              WRITE(SD,*)
     &'#Planet ipl JD(-offest) Sep(mas) PA(deg) dSep(mas) dPA(deg) Corr'
           CASE(8)
              WRITE(SD,*)
     &'#Planet ipl JD(-offset) Sep(mas) PA(deg) dSep(mas) dPA(deg)'
           END SELECT              
c...  Case of stellar (absolute) radial velocities'
           IF (ISDATA(2)) THEN
              WRITE(SD,*)
     &'Chosen stellar RV data (type 2 data) line format :'
              SELECT CASE (DATATYP(2))
              CASE(1)
                 WRITE(SD,*)'Day Month Year RV(km/s) dRV(km/s)'
              CASE(2)
                 WRITE(SD,*)'JD(-offset) RV(km/s) dRV(km/s)'
              END SELECT
           END IF
c...  Case of relative (planet/star) radial velocities'
           IF (ISDATA(3)) THEN
              WRITE(SD,*)
     &'Chosen relative planet RV data (type 3 data) line format :'
              SELECT CASE (DATATYP(3))
              CASE(1)
                 WRITE(SD,*)'#Planet Day Month Year RV(km/s) dRV(km/s)'
              CASE(2)
                 WRITE(SD,*)'#Planet JD(-offset) RV(km/s) dRV(km/s)'
              END SELECT
           END IF
c...   Case of absolute astrometry data (cf. Gaia)
           IF (ISDATA(4)) THEN
              WRITE(SD,*)
     &'Chosen absolute astrometry data (type 4 data) line format :'
              SELECT CASE (DATATYP(4))
              CASE(1)
                 WRITE(SD,*)
     &'Day Month Year Dec(mas) RA(mas) dDec(mas) dRA(mas) Corr'
              CASE(2)
                 WRITE(SD,*)
     &'Day Month Year Dec(mas) RA(mas) dDec(mas) dRA(mas)'           
              CASE(3)
                 WRITE(SD,*)
     &'JD(-offset) Dec(mas) RA(mas) dDec(mas) dRA(mas) Corr'
              CASE(4)
                 WRITE(SD,*)
     &'JD(-offset) Dec(mas) RA(mas) dDec(mas) dRA(mas)'
              CASE(5)
                 WRITE(SD,*)
     &'Day Month Year Sep(mas) PA(deg) dSep(mas) dPA(deg) Corr'
              CASE(6)
                 WRITE(SD,*)
     &'Day Month Year Sep(mas) PA(deg) dSep(mas) dPA(deg)'           
              CASE(7)
                 WRITE(SD,*)
     &'JD(-offest) Sep(mas) PA(deg) dSep(mas) dPA(deg) Corr'
              CASE(8)
                 WRITE(SD,*)
     &'JD(-offset) Sep(mas) PA(deg) dSep(mas) dPA(deg)'
              END SELECT     
           END IF
c... Distance or parallax
           WRITE(SD,*)
     &      'Distance of system (pc) or parallax (mas) (specify unit) ?'
           READ(5,*)STAR%DIST,UNIT
           CALL UPCASE(UNIT)
           IF (TRIM(UNIT).EQ.'PC') STAR%PARX = 1.d3/STAR%DIST
           IF (UNIT.EQ.'MAS') THEN
              STAR%PARX = STAR%DIST
              STAR%DIST = 1.d3/STAR%PARX
           END IF
c... Central mass
           WRITE(SD,*)'Mass of central body (Msun) (Initial guess) ?'
           READ(5,*)STAR%MASS
           STAR%MASS = STAR%MASS*SMAS ! Convert mass into the right unit
           STAR%DMASS = 0.d0
           STAR%V0 = 0.d0

           IF (CNEW.LT.5) THEN
              WRITE(SD,*)'Name of output file ?'
              READ(5,'(a)')FILES(1)
            !   FILES(1) = TRIM(FILES(1))//'.dat'
              FILES(1) = TRIM(FILES(1))
              WRITE(SD,*)'Name of dump file ?'
              READ(5,'(a)')FILES(2)
              WRITE(SD,*)'Dump frequency ?'
              READ(5,*)FREQDUMP 
              IF (MULTIPLA) THEN
                 NPAR = NEL*NPLA+1  ! Store central mass into the last variable
              ELSE
                 NPAR = NEL*NPLA
              END IF
              IF (ISDATA(2)) NPAR = NPAR+1+JITNUM
              NMOD = 0
c...  If MULTIPLA=F then NPAR = NEL*NPLA (variables = params * nb of orbits)
c...  If MULTIPLA=T but ISDATA(2)=F then NPAR=NEL*NPAR+1
c...       Variable #NPAR = ln(central mass)
c...  If ISDATA(2)=T and no jitter NPAR = NEL*NPAR+2
c...       Variables #NPAR-1 = V0, #NPAR = ln(central mass)
c...  If ISDATA(2)=T and jitter  NPAR = NEL*NPAR+3
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
                    MPRIOR(I)%MUNIT = 1.d0
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
cxx           CORR = (MOD(DATATYP(1),2).NE.0)
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
        LOGICAL ::      OK
        CHARACTER*(*), DIMENSION(NFIL) :: FILES ! #1 = output, #2 = dump,
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
        WRITE(18,*,IOSTAT=ERROR)ISDATA(1:NDATATYP),CORR
        WRITE(18,*,IOSTAT=ERROR)NPAR,NPLA,NPRIOR,JITNUM,
     &                          DATEINPUT,XYFORMAT,CORRNUM
        DO I=1,NLIM        
           WRITE(18,*,IOSTAT=ERROR)VPRIOR(I)%BOUND(1:2)
        END DO
        WRITE(18,*,IOSTAT=ERROR)PLA%MUNIT
        IF (NPRIOR.GE.0) THEN
           DO I=0,NPRIOR
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%TYP
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%ACOF(0:NPLA)
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%BCOF(0:NPLA)           
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%MUNIT,MPRIOR(I)%MEAN,
     &                                                   MPRIOR(I)%SDEV
              WRITE(18,*,IOSTAT=ERROR)MPRIOR(I)%BOUND
           END DO
        END IF
        WRITE(18,*,IOSTAT=ERROR)EPS,STAR%DIST,STAR%PARX
        WRITE(18,*,IOSTAT=ERROR)STAR%MASS,STAR%T0,STAR%V0,STAR%SIGJV
        DO I = 1,NFIL
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
        READ(18,*,IOSTAT=ERROR)ISDATA(1:NDATATYP),CORR
        RADVEL = (ISDATA(2).OR.ISDATA(3))
        READ(18,*,IOSTAT=ERROR)NPAR,NPLA,NPRIOR,JITNUM,
     &                         DATEINPUT,XYFORMAT,CORRNUM
        MULTIPLA = (NPLA.GT.1)
        MULTIPLA = (MULTIPLA.OR.ISDATA(2).OR.ISDATA(4))

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
        READ(18,*,IOSTAT=ERROR)PLA%MUNIT
        IF (NPRIOR.GE.0) THEN
           ALLOCATE(MPRIOR(0:NPRIOR))
           DO I = 0,NPRIOR
              ALLOCATE(MPRIOR(I)%ACOF(0:NPLA))
              ALLOCATE(MPRIOR(I)%BCOF(0:NPLA))     
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%TYP
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%ACOF(0:NPLA)
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%BCOF(0:NPLA)           
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%MUNIT,MPRIOR(I)%MEAN,
     &                                                 MPRIOR(I)%SDEV
              READ(18,*,IOSTAT=ERROR)MPRIOR(I)%BOUND(1:2)
           END DO
        END IF
        READ(18,*,IOSTAT=ERROR)EPS,STAR%DIST,STAR%PARX
        READ(18,*,IOSTAT=ERROR)STAR%MASS,STAR%T0,STAR%V0,STAR%SIGJV
        DO I = 1,NFIL
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
           SELECT CASE (METHOD)
           CASE(1)              ! Elliptic variables
             WRITE(SD,1)'Semi major axis ',PLA(I)%A,ABS(PLA(I)%DA),' AU'
             WRITE(SD,1)'Period ',PLA(I)%PER,ABS(PLA(I)%DPER),' days'
           CASE(2)
             WRITE(SD,1)'Periastron ',PLA(I)%Q,ABS(PLA(I)%DQ),' AU'
           END SELECT
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

        
C
C-----------------------------------------------------------------------------
C    Writes an astrometric data file
C-----------------------------------------------------------------------------
C

        SUBROUTINE WRITE_ASTROM(IDAT,FILNAM)

        USE DATA
        
        IMPLICIT NONE

        CHARACTER*(*) :: FILNAM ! File name  
        INTEGER*4 :: IDAT       ! Type of data
        INTEGER*4 :: NDAT       ! Current data
        INTEGER*4 ::   MM,YR    ! Month, Year
        REAL*8 ::      PX       ! Parallax
        REAL*8 ::      JJ,JD    ! Day, date
        REAL*8 ::      PA,DPA   ! Position angle
        REAL*8 ::      SEP,DSEP ! Separation
        REAL*8 ::      COV,RHO  ! Covariance, coef. correlation
        REAL*8 ::      XX,DXX   ! Declination
        REAL*8 ::      YY,DYY   ! Right Ascension
        INTEGER*4 ::   I,J,ERROR
        LOGICAL ::     OK

 7      FORMAT(i1,1x,f5.2,1x,i2,1x,i4,5(1x,f9.3))        
 8      FORMAT(I1,1x,f11.3,5(1x,f9.3))

        OPEN(15,FILE=FILNAM,STATUS='UNKNOWN')                                   
        SELECT CASE(DATATYP(IDAT))
        CASE(1,2,5,6)
           WRITE(15,'(a)')
     &    'I JJ MM YY   DEC(mas)  RA(mas)  DDec(mas)  DRa(mas) Corr'           
           WRITE(15,'(a)')
     &    '--------------------------------------------------------'
        CASE(3,4,7,8)
           WRITE(15,'(a)')
     &    'I JD     DEC(mas)  RA(mas)  DDec(mas)  DRa(mas) Corr'
           WRITE(15,'(a)')
     &    '----------------------------------------------------'
           WRITE(15,*)STAR%JDOFFSET        
        END SELECT
        PX = STAR%PARX
        DO I = 1,NPLA
           DO J = 1,PLA(I)%NDATAS
              SELECT CASE (DATATYP(IDAT))
              CASE(5,6,7,8)
c...  Conversion (dec,RA) -> (SEP,PA)
                 XX = PLA(I)%X(J)
                 YY = PLA(I)%Y(J)                 
                 SEP = SQRT(XX*XX+YY*YY)
                 PA = ATAN2(YY,XX)
                 DXX = PLA(I)%SIGX(J)
                 DYY = PLA(I)%SIGY(J)
                 COV = PLA(I)%RHOXY(J)*DXX*DYY
                 DSEP = SQRT(DXX*DXX*XX*XX+DYY*DYY*YY*YY+2.d0*COV*XX*YY)
     &                           /SEP
                 DPA = SQRT(DXX*DXX*YY*YY+DYY*DYY*XX*XX-2.d0*COV*XX*YY)
     &                           /(SEP*SEP)
                 RHO = (XX*YY*(DYY*DYY-DXX*DXX)+(XX*XX-YY*YY)*COV)
     &                /(SEP**3*DSEP*DPA)
              END SELECT
              SELECT CASE (DATATYP(IDAT))
              CASE(1,2)
                 CALL JD_TO_DATE(PLA(I)%TAS(J)+STAR%JDOFFSET,JJ,MM,YR)
                 WRITE(15,7,IOSTAT=ERROR)I,JJ,MM,YR,
     &                 PLA(I)%X(J)*PX,PLA(I)%Y(J)*PX,
     &                 PLA(I)%SIGX(J)*PX,PLA(I)%SIGY(J)*PX,
     &                 PLA(I)%RHOXY(J)       
              CASE(3,4)
                 WRITE(15,8,IOSTAT=ERROR)I,PLA(I)%TAS(J),
     &                 PLA(I)%X(J)*PX,PLA(I)%Y(J)*PX,
     &                 PLA(I)%SIGX(J)*PX,PLA(I)%SIGY(J)*PX,
     &                 PLA(I)%RHOXY(J)       
              CASE(5,6)
                 CALL JD_TO_DATE(PLA(I)%TAS(J)+STAR%JDOFFSET,JJ,MM,YR)
                 WRITE(15,7,IOSTAT=ERROR)I,JJ,MM,YR,
     &                 SEP*PX,PA*DR,DSEP*PX,DPA*DR,RHO    
              CASE(7,8)
                 WRITE(15,8,IOSTAT=ERROR)I,PLA(I)%TAS(J),
     &                SEP*PX,PA*DR,DSEP*PX,DPA*DR,RHO
              END SELECT
           END DO
        END DO
        CLOSE(15)

        END

C     
C
C-----------------------------------------------------------------------------
C    Reads the whole data file (4 types of data)
C-----------------------------------------------------------------------------
C

        SUBROUTINE READ_DATAFILE(FILNAM)

        USE DATA
        
        IMPLICIT NONE

        CHARACTER*(*) :: FILNAM ! File name  
        INTEGER*4 :: IPLA       ! Planet index
        INTEGER*4 :: IDAT       ! Type of data
        INTEGER*4 :: NDAT       ! Current data
        REAL*8 ::      JD       ! Julian Day
        REAL*8 ::      VV,DVV   ! Radial velocity
        REAL*8 ::      XX,DXX   ! Declination
        REAL*8 ::      YY,DYY   ! Right Ascension
        REAL*8 ::      RHO      ! coef. correlation
        INTEGER*4 ::   I,ERROR
        LOGICAL ::     OK,OKF
        CHARACTER*80 :: LIG 
        
        OPEN(15,FILE=FILNAM,STATUS='UNKNOWN')                                   

        OKF = .TRUE.
        STAR%JDOFFSET = 0
        READ(15,*)STAR%JDOFFSET
        READ(15,'(a)')LIG
        PLA%NDATAS = 0
        PLA%NDATVR = 0
        STAR%NDATVR = 0
        STAR%NDATAS = 0
c... First reading to determine sizes
        DO WHILE(OKF)
           READ(15,*,IOSTAT=ERROR)IDAT
           OKF = (ERROR.EQ.0)
           IF (OKF) THEN              
              READ(15,'(a)')LIG
              READ(15,'(a)')LIG
              OK = .TRUE.
              DO WHILE(OK)
                 CALL READ_INPUTLINE(IDAT,IPLA,JD,XX,YY,DXX,DYY,
     &                                              RHO,VV,DVV,OK)
                 IF (OK) THEN
                    SELECT CASE(IDAT)
                    CASE(1)
                       PLA(IPLA)%NDATAS = PLA(IPLA)%NDATAS+1
                    CASE(2)
                       STAR%NDATVR = STAR%NDATVR+1
                    CASE(3)
                       PLA(IPLA)%NDATVR = PLA(IPLA)%NDATVR+1
                    CASE(4)
                       STAR%NDATAS = STAR%NDATAS+1
                    END SELECT
                 END IF   
              END DO
           END IF
        END DO
c...  Allocate data arrays to the exact needed size
        DO I = 1,NPLA
           NDAT = PLA(I)%NDATAS
           IF (NDAT.NE.0) THEN
              ALLOCATE(PLA(I)%TAS(NDAT))
              ALLOCATE(PLA(I)%X(NDAT))
              ALLOCATE(PLA(I)%SIGX(NDAT))
              ALLOCATE(PLA(I)%Y(NDAT))
              ALLOCATE(PLA(I)%SIGY(NDAT))
              ALLOCATE(PLA(I)%RHOXY(NDAT))
              ALLOCATE(PLA(I)%SIGXM1(NDAT))
              ALLOCATE(PLA(I)%SIGYM1(NDAT))
              ALLOCATE(PLA(I)%SIGXM2(NDAT))
              ALLOCATE(PLA(I)%SIGYM2(NDAT))
              ALLOCATE(PLA(I)%UFAC(NDAT))
           END IF
           NDAT = PLA(I)%NDATVR
           IF (NDAT.NE.0) THEN
              ALLOCATE(PLA(I)%TVR(NDAT))
              ALLOCATE(PLA(I)%V(NDAT))
              ALLOCATE(PLA(I)%SIGV(NDAT))
              ALLOCATE(PLA(I)%SIGVM2(NDAT))
              ALLOCATE(PLA(I)%SIGV2(NDAT))
           END IF
        END DO
        NDAT = STAR%NDATAS
        IF (NDAT.NE.0) THEN
           ALLOCATE(STAR%TAS(NDAT))
           ALLOCATE(STAR%X(NDAT))
           ALLOCATE(STAR%SIGX(NDAT))
           ALLOCATE(STAR%Y(NDAT))
           ALLOCATE(STAR%SIGY(NDAT))
           ALLOCATE(STAR%RHOXY(NDAT))
           ALLOCATE(STAR%SIGXM1(NDAT))
           ALLOCATE(STAR%SIGYM1(NDAT))
           ALLOCATE(STAR%SIGXM2(NDAT))
           ALLOCATE(STAR%SIGYM2(NDAT))
           ALLOCATE(STAR%UFAC(NDAT))
        END IF
        NDAT = STAR%NDATVR
        IF (NDAT.NE.0) THEN
           ALLOCATE(STAR%TVR(NDAT))
           ALLOCATE(STAR%V(NDAT))
           ALLOCATE(STAR%SIGV(NDAT))
           ALLOCATE(STAR%SIGVM2(NDAT))
           ALLOCATE(STAR%SIGV2(NDAT))
        END IF

c...  Rewind file and reread to store data
        REWIND(15)
        READ(15,*)STAR%JDOFFSET
        READ(15,'(a)')LIG
        OKF = .TRUE.
        PLA%NDATAS = 0
        PLA%NDATVR = 0
        STAR%NDATVR = 0
        STAR%NDATAS = 0   
        DO WHILE(OKF)
           READ(15,*,IOSTAT=ERROR)IDAT
           OKF = (ERROR.EQ.0)
           IF (OKF) THEN              
              READ(15,'(a)')LIG
              READ(15,'(a)')LIG
              OK = .TRUE.
              DO WHILE(OK)
                 CALL READ_INPUTLINE(IDAT,IPLA,JD,XX,YY,DXX,DYY,RHO,
     &                                                     VV,DVV,OK)
                 IF (OK) THEN
                    SELECT CASE(IDAT)
                    CASE(1)              
                       PLA(IPLA)%NDATAS = PLA(IPLA)%NDATAS+1
                       NDAT = PLA(IPLA)%NDATAS
                       PLA(IPLA)%TAS(NDAT) = JD ! Time stored in JD-offset 
                       PLA(IPLA)%X(NDAT) = XX/STAR%PARX  
                       PLA(IPLA)%SIGX(NDAT) = DXX/STAR%PARX
                       PLA(IPLA)%Y(NDAT) = YY/STAR%PARX
                       PLA(IPLA)%SIGY(NDAT) = DYY/STAR%PARX
                       PLA(IPLA)%RHOXY(NDAT) = RHO
                       PLA(IPLA)%SIGXM1(NDAT) =
     &                           1.d0/PLA(IPLA)%SIGX(NDAT)
                       PLA(IPLA)%SIGYM1(NDAT) =
     &                           1.d0/PLA(IPLA)%SIGY(NDAT)
                       PLA(IPLA)%SIGXM2(NDAT) =
     &                           PLA(IPLA)%SIGXM1(NDAT)**2
                       PLA(IPLA)%SIGYM2(NDAT) =
     &                           PLA(IPLA)%SIGYM1(NDAT)**2
                       PLA(IPLA)%UFAC(NDAT) =
     &                           1.d0/(1.d0-PLA(IPLA)%RHOXY(NDAT)**2)
                    CASE(2)
                       STAR%NDATVR = STAR%NDATVR+1
                       NDAT = STAR%NDATVR
                       STAR%TVR(NDAT) = JD
                       STAR%V(NDAT) = VV*MPS ! Conversion to AU/day
                       STAR%SIGV(NDAT) = DVV*MPS 
                       STAR%SIGV2(NDAT) = STAR%SIGV(NDAT)**2
                       STAR%SIGVM2(NDAT) = 1.d0/STAR%SIGV2(NDAT)
                    CASE(3)
                       PLA(IPLA)%NDATVR = PLA(IPLA)%NDATVR+1
                       NDAT = PLA(IPLA)%NDATVR
                       PLA(IPLA)%TVR(NDAT) = JD
                       PLA(IPLA)%V(NDAT) = VV*MPS ! Conversion to AU/day
                       PLA(IPLA)%SIGV(NDAT) = DVV*MPS 
                       PLA(IPLA)%SIGV2(NDAT) = PLA(IPLA)%SIGV(NDAT)**2
                       PLA(IPLA)%SIGVM2(NDAT) =
     &                            1.d0/PLA(IPLA)%SIGV2(NDAT)
                    CASE(4)
                       STAR%NDATAS = STAR%NDATAS+1
                       NDAT = STAR%NDATAS
                       STAR%TAS(NDAT) = JD ! Time stored in JD-offset 
                       STAR%X(NDAT) = XX/STAR%PARX  
                       STAR%SIGX(NDAT) = DXX/STAR%PARX
                       STAR%Y(NDAT) = YY/STAR%PARX
                       STAR%SIGY(NDAT) = DYY/STAR%PARX
                       STAR%RHOXY(NDAT) = RHO
                       STAR%SIGXM1(NDAT) = 1.d0/STAR%SIGX(NDAT)
                       STAR%SIGYM1(NDAT) = 1.d0/STAR%SIGY(NDAT)
                       STAR%SIGXM2(NDAT) = STAR%SIGXM1(NDAT)**2
                       STAR%SIGYM2(NDAT) = STAR%SIGYM1(NDAT)**2
                       STAR%UFAC(NDAT) =
     &                           1.d0/(1.d0-STAR%RHOXY(NDAT)**2)       
                    END SELECT
                 END IF   
              END DO
           END IF
        END DO

        CLOSE(15)        

c....... Internal routine : reading an input line.................
        
        CONTAINS

        SUBROUTINE READ_INPUTLINE(IDAT,IPLA,JD,XX,YY,DXX,DYY,
     &                                          RHO,VV,DVV,OK)

          IMPLICIT NONE

          LOGICAL :: OK         ! True if succesful line reading
          INTEGER*4 :: IPLA     ! Planet index
          INTEGER*4 :: IDAT     ! Type of data
          INTEGER*4 ::   MM,YR  ! Month, Year
          REAL*8 ::      JJ,JD  ! Day, date
          REAL*8 ::      PA,DPA,CPA,SPA,VPA ! Position angle
          REAL*8 ::      SEP,DSEP,VSEP ! Separation
          REAL*8 ::      CV,RHO ! Covariance, coef. correlation
          REAL*8 ::      VV,DVV   ! Radial velocity
          REAL*8 ::      XX,DXX ! Declination
          REAL*8 ::      YY,DYY ! Right Ascension
          REAL*8 ::      DATE_TO_JD ! Conversion function
          INTEGER*4 :: ERROR
          
          RHO = 0.d0
          IPLA = 0
          SELECT CASE(IDAT)
          CASE(1)                  ! Relative astrometry data line
             SELECT CASE (DATATYP(1))
             CASE(1)
                READ(15,*,IOSTAT=ERROR)IPLA,JJ,MM,YR,XX,YY,DXX,DYY,RHO
             CASE(2)
                READ(15,*,IOSTAT=ERROR)IPLA,JJ,MM,YR,XX,YY,DXX,DYY    
             CASE(3)
                READ(15,*,IOSTAT=ERROR)IPLA,JD,XX,YY,DXX,DYY,RHO           
             CASE(4)
                READ(15,*,IOSTAT=ERROR)IPLA,JD,XX,YY,DXX,DYY    
             CASE(5)
                READ(15,*,IOSTAT=ERROR)IPLA,JJ,MM,YR,SEP,PA,DSEP,DPA,RHO 
             CASE(6)
                READ(15,*,IOSTAT=ERROR)IPLA,JJ,MM,YR,SEP,PA,DSEP,DPA    
             CASE(7)
                READ(15,*,IOSTAT=ERROR)IPLA,JD,SEP,PA,DSEP,DPA,RHO           
             CASE(8)
                READ(15,*,IOSTAT=ERROR)IPLA,JD,SEP,PA,DSEP,DPA    
             END SELECT
          CASE(2)               ! Stellar radial velocity data
             SELECT CASE (DATATYP(2))
             CASE(1)
                READ(15,*,IOSTAT=ERROR)JJ,MM,YR,VV,DVV           
             CASE(2)
                READ(15,*,IOSTAT=ERROR)JD,VV,DVV    
             END SELECT
          CASE(3)   ! Relative radial velocity data
             SELECT CASE (DATATYP(3))
             CASE(1)
                READ(15,*,IOSTAT=ERROR)IPLA,JJ,MM,YR,VV,DVV           
             CASE(2)
                READ(15,*,IOSTAT=ERROR)IPLA,JD,VV,DVV    
             END SELECT
          CASE(4)
             SELECT CASE (DATATYP(4))
             CASE(1)
                READ(15,*,IOSTAT=ERROR)JJ,MM,YR,XX,YY,DXX,DYY,RHO
             CASE(2)
                READ(15,*,IOSTAT=ERROR)JJ,MM,YR,XX,YY,DXX,DYY    
             CASE(3)
                READ(15,*,IOSTAT=ERROR)JD,XX,YY,DXX,DYY,RHO           
             CASE(4)
                READ(15,*,IOSTAT=ERROR)JD,XX,YY,DXX,DYY    
             CASE(5)
                READ(15,*,IOSTAT=ERROR)JJ,MM,YR,SEP,PA,DSEP,DPA,RHO 
             CASE(6)
                READ(15,*,IOSTAT=ERROR)JJ,MM,YR,SEP,PA,DSEP,DPA    
             CASE(7)
                READ(15,*,IOSTAT=ERROR)JD,SEP,PA,DSEP,DPA,RHO           
             CASE(8)
                READ(15,*,IOSTAT=ERROR)JD,SEP,PA,DSEP,DPA    
             END SELECT
          END SELECT
          OK = (ERROR.EQ.0)      
          IF (OK) THEN
             IF (XYFORMAT.EQ.2) THEN
c... Conversion (sep,PA) -> (dec,RA)
                PA = PA/DR
                DPA = DPA/DR
                CPA = COS(PA)
                SPA = SIN(PA)
                XX = SEP*CPA
                YY = SEP*SPA
                VSEP = DSEP*DSEP
                VPA = DPA*DPA
                CV = RHO*DSEP*DPA
                DXX = SQRT(CPA*CPA*VSEP+SEP*SEP*SPA*SPA*VPA
     &               -2.d0*SEP*CPA*SPA*CV)                 
                DYY = SQRT(SPA*SPA*VSEP+SEP*SEP*CPA*CPA*VPA
     &               +2.d0*SEP*CPA*SPA*CV)
                RHO = (CPA*SPA*VSEP-SEP*SPA*CPA*VPA
     &               +(CPA*CPA-SPA*SPA)*SEP*CV)/(DXX*DYY)                
             END IF
             IF (DATEINPUT.EQ.1) THEN
c... Conversion JJ-MM-YR -> JD
                JD = DATE_TO_JD(JJ,MM,YR)-STAR%JDOFFSET
             END IF
          END IF
          END SUBROUTINE READ_INPUTLINE

        END

C
C -----------------------------------------------------------------------------
C       Writing an output file (with (O,w) solutions doubling)
C -----------------------------------------------------------------------------
C
        SUBROUTINE WRITE_DISTRIB_DOUBLING(NMOD,DEV)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    NMOD,         ! Number of models 
     &                  NSAV,         ! Number of pars. to store per planet
     &                  LOG_TO_INT    ! Conversion function
        CHARACTER*(*) :: DEV    ! 
        REAL*8, DIMENSION(:), ALLOCATABLE ::
     &                  NN,           ! Mean motion (wrt/q if univ. var.)
     &                  O,OM,         ! Omega's ± Pi
     &                  CI,SI,        ! cos(i),sin(i)
     &                  CI2,SI2       ! cos^2(i/2), sin^2(i/2)

        REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: DATA4
        INTEGER*4       I,J,DNPLA,IC

        DNPLA = 2**NPLA
        NSAV = NEL+4
        ALLOCATE(DATA4(DNPLA*(NMOD+1),NSAV,NPLA))
        ALLOCATE(NN(NPLA))
        ALLOCATE(CI2(NPLA))
        ALLOCATE(SI2(NPLA))
        ALLOCATE(CI(NPLA))
        ALLOCATE(SI(NPLA))
        ALLOCATE(OM(NPLA))
        ALLOCATE(O(NPLA))

        CALL STORE_SOLUTION(0,PSTART,NN,CI2,SI2,CI,SI,OM,O)

        DO I = 1,NMOD
           CALL STORE_SOLUTION(DNPLA*I,PSAV(:,I),
     &                               NN,CI2,SI2,CI,SI,OM,O)
        END DO
c...      Save everything into output file
        OPEN(18,FILE=DEV,STATUS='UNKNOWN')
        WRITE(18,*)METHOD,DATEINPUT,NPLA,LOG_TO_INT(MULTIPLA),
     &        LOG_TO_INT(RADVEL),JITNUM,NPRIOR,NLIM
        WRITE(18,*)LOG_TO_INT(ISDATA(1)),LOG_TO_INT(ISDATA(2)),
     &             LOG_TO_INT(ISDATA(3)),LOG_TO_INT(ISDATA(4))   
        WRITE(18,*)STAR%JDOFFSET,STAR%DIST,STAR%PARX
        WRITE(18,*)STAR%NDATAS,STAR%NDATVR
        DO I = 1,STAR%NDATAS
           WRITE(18,*)STAR%TAS(I),STAR%X(I),STAR%Y(I),
     &                            STAR%SIGX(I),STAR%SIGY(I)
        END DO
        DO I = 1,STAR%NDATVR
           WRITE(18,*)STAR%TVR(I),STAR%V(I),STAR%SIGV(I)
        END DO
        WRITE(18,*)PLA%NDATAS
        WRITE(18,*)PLA%NDATVR
        DO I = 1,NPLA
           DO J = 1,PLA(I)%NDATAS
              WRITE(18,*)PLA(I)%TAS(J),PLA(I)%X(J),PLA(I)%Y(J),
     &             PLA(I)%SIGX(J),PLA(I)%SIGY(J),PLA(I)%RHOXY(J)
           END DO
           DO J = 1,PLA(I)%NDATVR
              WRITE(18,*)PLA(I)%TVR(J),PLA(I)%V(J),PLA(I)%SIGV(J)
           END DO
        END DO
        WRITE(18,*)PLA%MUNIT        
c... Only store relevant priors #n+1..nprior / #0..n just tell masses are >0
        DO I = NPLA+1,NPRIOR
           WRITE(18,*)MPRIOR(I)%TYP
           WRITE(18,*)MPRIOR(I)%ACOF(0:NPLA)
           WRITE(18,*)MPRIOR(I)%BCOF(0:NPLA)
           WRITE(18,*)MPRIOR(I)%MUNIT,MPRIOR(I)%MEAN,MPRIOR(I)%SDEV,
     &                                               MPRIOR(I)%BOUND
        END DO           
        DO I = 1,NLIM
           WRITE(18,*)VPRIOR(I)%BOUND(1:2)
        END DO

        WRITE(18,*)DNPLA*(NMOD+1),NSAV,NPLA
        WRITE(18,*)DATA4
        CLOSE(18)
        DEALLOCATE(DATA4)
        DEALLOCATE(NN)
        DEALLOCATE(CI2)
        DEALLOCATE(SI2)
        DEALLOCATE(CI)
        DEALLOCATE(SI)
        DEALLOCATE(OM)
        DEALLOCATE(O)

c....... Internal routine : storing 1 solution .................
        
        CONTAINS

        SUBROUTINE STORE_SOLUTION(DK,P,NN,CI2,SI2,CI,SI,OM,O)

          IMPLICIT NONE

          INTEGER*4 :: DK       ! Position in DATA array
          REAL*8 ::    P(NPAR)  ! Parameters of the solution
          REAL*8, DIMENSION(NPLA) ::
     &                  NN,           ! Mean motion (wrt/q if univ. var.)
     &                  O,OM,         ! Omega's ± Pi
     &                  CI,SI,        ! cos(i),sin(i)
     &                  CI2,SI2       ! cos^2(i/2), sin^2(i/2)
          REAL*8 ::     CHI2,         ! Chi2
     &                  FMAP         ! MAP

          INTEGER*4 :: IC,K
          
          CALL ELEMENTS(P,NN,PLA%A,PLA%EXC,PLA%EXQ,
     &         PLA%CW,PLA%SW,CI2,SI2,PLA%CP,PLA%SP,CI,SI,
     &         PLA%COM,PLA%SOM,PLA%CO,PLA%SO,PLA%TP,PLA%MU)
c...  Note : In the case of universal variables,
c...           the periastron is stored here in pla%a 
          IF (MULTIPLA) STAR%MASS = EXP(P(NPAR))
          PLA%PER = DPI/NN
          PLA%W = ATAN2(PLA%SW,PLA%CW)
          PLA%PHI = ATAN2(PLA%SP,PLA%CP)
          PLA%OM = 0.5d0*(PLA%W+PLA%PHI)
          PLA%O = 0.5d0*(PLA%W-PLA%PHI)
          PLA%INC = 2.d0*ATAN2(SQRT(SI2),SQRT(CI2))
          CHI2 = P(NPAR+1)
          FMAP = P(NPAR+2)
          DO IC = 1,DNPLA
             OM(1:NPLA) = PLA%OM
             O(1:NPLA) = PLA%O
             DO K = 1,NPLA
                IF (BTEST(IC-1,K-1)) THEN
                   OM(K) = MOD(PLA(K)%OM+DPI,DPI)-PI
                   O(K) = MOD(PLA(K)%O+DPI,DPI)-PI
                ELSE
                   OM(K) = PLA(K)%OM
                   O(K) = PLA(K)%O
                END IF
                DATA4(DK+IC,1:NSAV,K) = (/ SNGL(PLA(K)%A),
     &               SNGL(PLA(K)%PER),SNGL(PLA(K)%EXC),SNGL(OM(K)),
     &               SNGL(PLA(K)%INC),SNGL(O(K)),SNGL(PLA(K)%TP),
     &               SNGL(PLA(K)%MU/PLA(K)%MUNIT),SNGL(STAR%MASS/SMAS),
     &               SNGL(CHI2),SNGL(FMAP) /)
             END DO
          END DO

          END SUBROUTINE STORE_SOLUTION
        
        END

C
C -----------------------------------------------------------------------------
C       Writing an output file
C -----------------------------------------------------------------------------
C
        SUBROUTINE WRITE_DISTRIB(NMOD,DEV)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    NMOD,         ! Number of models 
     &                  NSAV,         ! Number of pars. to store per planet
     &                  LOG_TO_INT    ! Conversion function
        CHARACTER*(*) :: DEV    ! 
        REAL*8, DIMENSION(:), ALLOCATABLE ::
     &                  NN,           ! Mean motion (wrt/q if univ. var.)
     &                  CI,SI,        ! cos(i),sin(i)
     &                  CI2,SI2       ! cos^2(i/2), sin^2(i/2)
        REAL*8 ::       CHI2,         ! Chi2
     &                  FMAP,         ! MAP
     &                  SIGMA         ! Cumulative mass
        REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: DATA4
        INTEGER*4       I,J

        STAR%SIGJV = 0.d0
        NSAV = NEL+6
        ALLOCATE(DATA4(NMOD+1,NSAV,NPLA))
        ALLOCATE(NN(NPLA))
        ALLOCATE(CI(NPLA))
        ALLOCATE(SI(NPLA))
        ALLOCATE(CI2(NPLA))
        ALLOCATE(SI2(NPLA))
c...  Note : In the case of universal variables,
c...           the periastron is stored here in pla%a 
        CALL STORE_SOLUTION(0,PSTART,NN,CI2,SI2,CI,SI)
C
        DO I = 1,NMOD
           CALL STORE_SOLUTION(I,PSAV(:,I),NN,CI2,SI2,CI,SI)
        END DO
c...      Save everything into output file
        OPEN(18,FILE=DEV,STATUS='UNKNOWN')
        WRITE(18,*)METHOD,DATEINPUT,NPLA,LOG_TO_INT(MULTIPLA),
     &        LOG_TO_INT(RADVEL),JITNUM,NPRIOR,NLIM
        WRITE(18,*)LOG_TO_INT(ISDATA(1)),LOG_TO_INT(ISDATA(2)),
     &             LOG_TO_INT(ISDATA(3)),LOG_TO_INT(ISDATA(4))   
        WRITE(18,*)STAR%JDOFFSET,STAR%DIST,STAR%PARX
        WRITE(18,*)STAR%NDATAS,STAR%NDATVR
        DO I = 1,STAR%NDATAS
           WRITE(18,*)STAR%TAS(I),STAR%X(I),STAR%Y(I),
     &                            STAR%SIGX(I),STAR%SIGY(I)
        END DO
        DO I = 1,STAR%NDATVR
           WRITE(18,*)STAR%TVR(I),STAR%V(I),STAR%SIGV(I)
        END DO
        WRITE(18,*)PLA%NDATAS
        WRITE(18,*)PLA%NDATVR
        DO I = 1,NPLA
           DO J = 1,PLA(I)%NDATAS
              WRITE(18,*)PLA(I)%TAS(J),PLA(I)%X(J),PLA(I)%Y(J),
     &             PLA(I)%SIGX(J),PLA(I)%SIGY(J),PLA(I)%RHOXY(J)
           END DO
           DO J = 1,PLA(I)%NDATVR
              WRITE(18,*)PLA(I)%TVR(J),PLA(I)%V(J),PLA(I)%SIGV(J)
           END DO
        END DO
        WRITE(18,*)PLA%MUNIT        
c... Only store relevant priors #n+1..nprior / #0..n just tell masses are >0
        DO I = NPLA+1,NPRIOR
           WRITE(18,*)MPRIOR(I)%TYP
           WRITE(18,*)MPRIOR(I)%ACOF(0:NPLA)
           WRITE(18,*)MPRIOR(I)%BCOF(0:NPLA)
           WRITE(18,*)MPRIOR(I)%MUNIT,MPRIOR(I)%MEAN,MPRIOR(I)%SDEV,
     &                                               MPRIOR(I)%BOUND
        END DO           
        DO I = 1,NLIM
           WRITE(18,*)VPRIOR(I)%BOUND(1:2)
        END DO

        WRITE(18,*)NMOD+1,NSAV,NPLA
        WRITE(18,*)DATA4
        CLOSE(18)
        
        DEALLOCATE(DATA4)
        DEALLOCATE(NN)
        DEALLOCATE(CI)
        DEALLOCATE(SI)
        DEALLOCATE(CI2)
        DEALLOCATE(SI2)

c....... Internal routine : storing 1 solution .................
        
        CONTAINS

        SUBROUTINE STORE_SOLUTION(DK,P,NN,CI2,SI2,CI,SI)

          IMPLICIT NONE

          INTEGER*4 :: DK       ! Position in DATA array
          REAL*8 ::    P(NPAR)  ! Parameters of the solution
          REAL*8, DIMENSION(NPLA) ::
     &                  NN,           ! Mean motion (wrt/q if univ. var.)
     &                  CI,SI,        ! cos(i),sin(i)
     &                  CI2,SI2       ! cos^2(i/2), sin^2(i/2)
          REAL*8 ::     CHI2,         ! Chi2
     &                  FMAP         ! MAP

          INTEGER*4 :: IC,K
          
          CALL ELEMENTS(P,NN,PLA%A,PLA%EXC,PLA%EXQ,
     &         PLA%CW,PLA%SW,CI2,SI2,PLA%CP,PLA%SP,CI,SI,
     &         PLA%COM,PLA%SOM,PLA%CO,PLA%SO,PLA%TP,PLA%MU)

          STAR%MASS = EXP(P(NPAR))
          STAR%V0 = P(NEL*NPLA+1)
          IF (JITNUM.EQ.1) STAR%SIGJV = EXP(P(NEL*NPLA+2))

          PLA%PER = DPI/NN
          PLA%OM = ATAN2(PLA%SOM,PLA%COM)
          PLA%O = ATAN2(PLA%SO,PLA%CO)
          PLA%INC = ATAN2(SI,CI)
          CHI2 = P(NPAR+1)
          FMAP = P(NPAR+2)

          DO K = 1,NPLA
             DATA4(DK+1,1:NSAV,K) = (/ SNGL(PLA(K)%A),SNGL(PLA(K)%PER),
     &            SNGL(PLA(K)%EXC),SNGL(PLA(K)%OM),SNGL(PLA(K)%INC),
     &            SNGL(PLA(K)%O),SNGL(PLA(K)%TP),
     &            SNGL(PLA(K)%MU/PLA(K)%MUNIT),SNGL(STAR%V0/MPS),
     &            SNGL(STAR%SIGJV/MPS),SNGL(STAR%MASS/SMAS),
     &            SNGL(CHI2),SNGL(FMAP) /)
          END DO

          END SUBROUTINE STORE_SOLUTION

        END


      
       
