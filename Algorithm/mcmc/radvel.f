C      
C        Radial velocity utilitary routines
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

