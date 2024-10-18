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

