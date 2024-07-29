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
C       Computation of RV contribution to Levenberg-Marquardt coefficients
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQRADVEL(N,JITTER,P,ALPHA,BETA,CHI2)

        USE DATA

        IMPLICIT NONE

        INTEGER*4       N               ! Dimension of parameter space
        
        REAL*8, DIMENSION(NPAR) :: DVF
        REAL*8 ::       ALPHA(N,N),     ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  P(N),           ! The parameters
     &                  CHI2,           ! Chi2
     &                  F1,F2,FS,       ! Intermediaires
     &                  SIGJV,          ! Velocity Jitter
     &                  SIGV2P,         ! Incertitude^2 augmentée
     &                  VF,                ! Fitted velocity
     &                  DV              ! Normalized difference   
        INTEGER*4 ::      I,J,K,L,SJ,iv,jv  
        LOGICAL ::        JITTER
        
c        REAL*8 ::       PB(NPAR),es,vfq,vfz ! Vitesse radiale

        SJ = NEL*NPLA+2
c        pb(1:n) = p(1:n)
        DVF(1:N) = 0.d0

        IF (JITTER) THEN
           SIGJV = EXP(P(SJ))   ! Velocity Jitter
        ELSE
           SIGJV = 0.d0
        END IF
c     es = 1.d-8

        DO I = 1,STAR%NDATVR
c     do iv=1,n
c          pb(iv) = p(iv)+es
c          call mrqvfit(n,star%tvr(i),pb(1:n),vfz,dvf(1:n))
c         pb(iv) = p(iv)-es
c          call mrqvfit(n,star%tvr(i),pb(1:n),vfq,dvf(1:n))
           CALL MRQVFIT(N,STAR%TVR(I),P(1:N),VF,DVF(1:N))
c       WRITE(SD,*)iv,sngl(vfz),sngl(vfq),sngl((vfz-vfq)/(2.d0*es)),
c     &                 sngl(dvf(iv))
c         end do
c         stop
           DV = STAR%V(I)-VF
           SIGV2P = STAR%SIGV2(I)+SIGJV*SIGJV
           CHI2 = CHI2+DV*DV/SIGV2P
c     WRITE(SD,*)(v(i)-vf)*AU/YEAR,SIG(I)*AU/YEAR,
c     &(V(I)-VF)*(V(I)-VF)/(SIG(I)*SIG(I))
c          WRITE(SD,*)sngl(x(i)),sngl(xf),sngl(y(i)),sngl(yf),
c     &           sngl(sigm2x(i)),sngl(sigm2y(i)),sngl(chi2)
          
           F1 = DV/SIGV2P
           FS = SIGJV/SIGV2P
           BETA(1:N) = BETA(1:N)+F1*DVF(1:N)
           DO J=1,N
             F2 = DVF(J)/SIGV2P
             ALPHA(J:N,J) = ALPHA(J:N,J)+F2*DVF(J:N)
             IF (JITTER.AND.(SJ.GT.J))
     &             ALPHA(SJ,J) = ALPHA(SJ,J)+2.d0*DV*FS*FS*DVF(J)
           END DO  
           IF (JITTER) THEN
              ALPHA(N,SJ) = ALPHA(N,SJ)+2.d0*DV*FS*FS*DVF(N)
              BETA(SJ) = BETA(SJ)+FS*FS*DV*DV
              ALPHA(SJ,SJ) = ALPHA(SJ,SJ)
     &                      +FS*FS*DV*DV*(4.d0*FS*SIGJV-2.d0)
           END IF
        END DO

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
