C
C-----------------------------------------------------------------------------
C    Reads an astrometric data file
C-----------------------------------------------------------------------------
C

        SUBROUTINE READ_ASTROM(IDAT,FILNAM)

        USE DATA
        
        IMPLICIT NONE

        CHARACTER*(*) :: FILNAM ! File name  
        INTEGER*4 :: IPLA       ! Planet index
        INTEGER*4 :: IDAT       ! Type of data
        INTEGER*4 :: NDAT       ! Current data
        INTEGER*4 ::   MM,YR    ! Month, Year
        REAL*8 ::      JJ,JD    ! Day, date
        REAL*8 ::      PA,DPA,CPA,SPA,VPA ! Position angle
        REAL*8 ::      SEP,DSEP,VSEP ! Separation
        REAL*8 ::      COV,RHO  ! Covariance, coef. correlation
        REAL*8 ::      XX,DXX   ! Declination
        REAL*8 ::      YY,DYY   ! Right Ascension
        REAL*8 ::      DATE_TO_JD ! Conversion function
c        REAL*8, DIMENSION(NDATMAX) :: X,Y,T ! Astrometric data
c        REAL*8, DIMENSION(NDATMAX) :: SIGX,SIGY,SIGXM2,SIGYM2,
c     &                                SIGXM1,SIGYM1,RHOXY,UFAC   
        INTEGER*4 ::   I,ERROR
        LOGICAL ::     OK
        CHARACTER*80 :: LIG 
        
        OPEN(15,FILE=FILNAM,STATUS='UNKNOWN')                                   
      
        READ(15,'(a)')LIG                                                       
        READ(15,'(a)')LIG                                                       
        SELECT CASE(DATATYP(IDAT))
        CASE(3,4,7,8)
           READ(15,*)STAR%JDOFFSET        
        END SELECT
        OK = .TRUE.
        PLA%NDATAS = 0
        DO WHILE(OK)
           SELECT CASE (DATATYP(IDAT))
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
           IF (MOD(DATATYP(1),2).EQ.0) RHO=0.d0
           OK=(ERROR.EQ.0)
           IF (OK) PLA(IPLA)%NDATAS = PLA(IPLA)%NDATAS+1
        END DO
        DO I = 1,NPLA
           NDAT = PLA(I)%NDATAS
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
        END DO
        REWIND(15)
        READ(15,'(a)')LIG                                                       
        READ(15,'(a)')LIG                                                       
        STAR%JDOFFSET = 0.d0
        SELECT CASE(DATATYP(IDAT)) 
        CASE(3,4,7,8)
           READ(15,*)STAR%JDOFFSET        
        END SELECT
        PLA%NDATAS = 0
        OK = .TRUE.                                                             
        DO WHILE(OK)
           SELECT CASE (DATATYP(IDAT))
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
           IF (MOD(DATATYP(1),2).EQ.0) RHO=0.d0
           OK=(ERROR.EQ.0)
           IF (OK) THEN   
              PLA(IPLA)%NDATAS = PLA(IPLA)%NDATAS+1
              NDAT = PLA(IPLA)%NDATAS
              SELECT CASE (DATATYP(IDAT))
              CASE(5,6,7,8)
c... Conversion (sep,PA) -> (dec,RA)
                 PA = PA/DR
                 DPA = DPA/DR
                 CPA = COS(PA)
                 SPA = SIN(PA)
                 XX = SEP*CPA
                 YY = SEP*SPA
                 VSEP = DSEP*DSEP
                 VPA = DPA*DPA
                 COV = RHO*DSEP*DPA
                 DXX = SQRT(CPA*CPA*VSEP+SEP*SEP*SPA*SPA*VPA
     &                       -2.d0*SEP*CPA*SPA*COV)                 
                 DYY = SQRT(SPA*SPA*VSEP+SEP*SEP*CPA*CPA*VPA
     &                   +2.d0*SEP*CPA*SPA*COV)
                 RHO = (CPA*SPA*VSEP-SEP*SPA*CPA*VPA
     &                   +(CPA*CPA-SPA*SPA)*SEP*COV)/(DXX*DYY)
              END SELECT
              SELECT CASE (DATATYP(IDAT))
              CASE(1,2,5,6)
c... Conversion JJ-MM-YR -> JD
                 JD = DATE_TO_JD(JJ,MM,YR)-STAR%JDOFFSET
              CASE(3,4,7,8)
                 JD = JD   ! +STAR%JDOFFSET
              END SELECT
              
              PLA(IPLA)%TAS(NDAT) = JD ! Time stored in Jul days-JDoffset 
              PLA(IPLA)%X(NDAT) = XX/STAR%PARX  
              PLA(IPLA)%SIGX(NDAT) = DXX/STAR%PARX
              PLA(IPLA)%Y(NDAT) = YY/STAR%PARX
              PLA(IPLA)%SIGY(NDAT) = DYY/STAR%PARX
              PLA(IPLA)%RHOXY(NDAT) = RHO
              PLA(IPLA)%SIGXM1(NDAT) = 1.d0/PLA(IPLA)%SIGX(NDAT)
              PLA(IPLA)%SIGYM1(NDAT) = 1.d0/PLA(IPLA)%SIGY(NDAT)
              PLA(IPLA)%SIGXM2(NDAT) = PLA(IPLA)%SIGXM1(NDAT)**2
              PLA(IPLA)%SIGYM2(NDAT) = PLA(IPLA)%SIGYM1(NDAT)**2
              PLA(IPLA)%UFAC(NDAT) =
     &                1.d0/(1.d0-PLA(IPLA)%RHOXY(NDAT)**2)
           END IF
        END DO        
        CLOSE(15)        
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

C     
