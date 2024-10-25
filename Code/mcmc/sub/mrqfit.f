C
C      Levenberg-Marquardt common routines      
C      
C -----------------------------------------------------------------------------
C       Search for minimal chi2 by Levenberg-Marquardt method
C -----------------------------------------------------------------------------
C
        SUBROUTINE MRQMIN(N,P,COV,CHI2,FMAP)

        USE DATA

        IMPLICIT NONE

        REAL*8, PARAMETER :: LAMBDA0 = 1d-3
        INTEGER*4       N       ! Dimension of parameter space

        REAL*8          ALPHA(N,N),     ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  COV(N,N),       ! Covariance matrix
     &                  SIGP(N),        ! Variance for parameters
     &                  DP(N),          ! Delta-Parameters
     &                  P(N),           ! Parameters
     &                  CHI2,FMAP,      ! Chi2 (total), FMAP
     &                  LAMBDA          ! Parameterr of method
        INTEGER*4       I,ITS,iv,jv
        LOGICAL         OK,TEST

        LAMBDA = LAMBDA0
        CALL MRQCOF(N,P,ALPHA,BETA,CHI2,FMAP)
        WRITE(SD,*)'Initial chi2 = ',chi2
        OK = .FALSE.
        ITS = 0
        DO WHILE(.NOT.OK)        
           CALL MRQSTEP(N,LAMBDA,P,ALPHA,BETA,CHI2,FMAP,OK)
           IF ((MOD(ITS,1000).EQ.0).or.(lambda.lt.1d-4)) THEN
            WRITE(SD,*)'        ----> ',sngl(CHI2),' <----',sngl(lambda)
         END IF
           ITS = ITS+1
           if (chi2.ne.chi2) stop
        END DO
        WRITE(SD,*)'End of Levenberg-Marquardt'
        CALL GAUSS(N,ALPHA,BETA,COV,DP,TEST)
        DO I=1,N
          SIGP(I)=SQRT(COV(I,I))
        END DO
        END

C
C -----------------------------------------------------------------------------
C       One step of Levenberg-Marquardt method
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQSTEP(N,LAMBDA,P,ALPHA,BETA,CHI2,FMAP,OK)

        USE DATA

        IMPLICIT NONE

        REAL*8, PARAMETER :: LAMBDA0 = 1d-3
        INTEGER*4       N       ! Dimension of parameter space
        REAL*8          A(N,N),ATRY(N,N),! Trial matrix
     &                  ALPHA(N,N),     ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  B(N),BTRY(N),   ! Trial vector
     &                  COV(N,N),       ! Covariance matrix (useless here)
     &                  DP(N),          ! Delta-Parameters
     &                  PTRY(N),        ! Trial new parameters
     &                  P(N),           ! Parameter set
     &                  CHI2,FMAP,      ! Chi2,MAP (CHi2 = Real Chi2 + FMAP)
     &                  CTRY,           ! Trial Chi2
     &                  LAMBDA          ! parameter for mrq method

        LOGICAL         OK,TEST,chk     ! Test for stopping
        INTEGER*4       I,J

        A(1:N,1:N) = ALPHA(1:N,1:N)
        B(1:N) = BETA(1:N)
c...  Multiply diagonal by 1+lambda
        DO I=1,N
           A(I,I) = A(I,I)*(1.0d0+LAMBDA)
        END DO        
c        WRITE(SD,*)'a',sngl(a),'b',sngl(b)
        CALL GAUSS(N,A,B,COV,DP,TEST)
c...  Trial step
        
        PTRY(1:N) = P(1:N)+DP(1:N)
c        do i=1,n
c          print*,i,'dp',sngl(dp(i)),'p',sngl(p(i)),sngl(ptry(i))
c        end do
c...  Check if trial step is still in the authorized bounds
        CALL MRQCHECK(N,PTRY,CHK)
        
c...  Compute coefs and trial Chi2. If error then assume bad step
        IF (.NOT.CHK) THEN
          CALL MRQCOF(N,PTRY,ATRY,BTRY,CTRY,FMAP)        
          IF (CTRY.NE.CTRY) CHK = .TRUE. 
        END IF
c...  We are sill in valid region
        IF (.NOT.CHK) THEN
c...    Degraded CHi2. Refuse step and increase lambda
          IF (CTRY.GT.CHI2) THEN
            IF (LAMBDA.EQ.0.d0) THEN
              LAMBDA = LAMBDA0
            ELSE
              LAMBDA = 10.d0*LAMBDA
            END IF
            OK=.FALSE.
c...    Improved Chi2. Accept step and decrease lambda
          ELSE
            P(1:N) = PTRY(1:N)
            BETA(1:N) = BTRY(1:N)
            ALPHA(1:N,1:N) = ATRY(1:N,1:N)
            LAMBDA = 0.1d0*LAMBDA
            OK = (ABS(CHI2-CTRY).LT.EPS)
            CHI2 = CTRY
         END IF
c... Case where the trial step is out of bounds or erroneous. Restart
        ELSE  
          OK = .FALSE.
          CALL MRQRENEW(N,P)
          CALL MRQCOF(N,P,ALPHA,BETA,CHI2,FMAP) 
          LAMBDA = LAMBDA0
c          CHK = .TRUE.
c          DO WHILE(CHK)
c            A(1:N,1:N) = ALPHA(1:N,1:N)
c            B(1:N) = BETA(1:N)
c...  Multiply diagonal by 1+lambda
c            DO I=1,N
c              A(I,I) = A(I,I)*(1.0d0+LAMBDA)
c            END DO        
c        WRITE(SD,*)'a',sngl(a),'b',sngl(b)
c            CALL GAUSS(N,A,B,COV,DP,TEST)
c            PTRY(1:N) = P(1:N)+DP(1:N)
c            CALL MRQCHECK(N,PTRY,CHK)
c            IF (CHK) LAMBDA = LAMBDA*10.d0
c          END DO
        END IF
          
        END

C
C -----------------------------------------------------------------------------
C       Computation of HCI contribution to Levenberg-Marquardt coefficients
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQASTROM(N,P,ALPHA,BETA,CHI2)

        USE DATA

        IMPLICIT NONE

        INTEGER*4       N               ! Dimension of parameter space
        
        REAL*8, DIMENSION(NPAR,NPLA) :: DXF,DYF
        REAL*8, DIMENSION(NPLA) :: XF,YF
        REAL*8 ::       ALPHA(N,N),     ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  P(N),           ! The parameters
     &                  CHI2,           ! Chi2
     &                  F1X,F1Y,F2XX,F2YY,F2XY, ! Intermediaires
     &                  SXM1,SYM1,RHO,UF,!
     &                  DX,DY           ! Normalized differences   
        INTEGER*4       I,J,K,L,SJ,iv,jv  

c        REAL*8 ::       PB(NPAR),
c     &       xfz(npla),yfz(npla),xfq(npla),yfq(npla),es
        
c        pb(1:n) = p(1:n)
        DXF = 0.d0
        DYF = 0.d0
c        es = 1.d-8
        DO L = 1,NPLA
           DO I = 1,PLA(L)%NDATAS
c         do iv=1,n
c         pb(iv) = p(iv)+es
c         call mrqposfit(n,pla(l)%tas(i),pb(1:n),xfz,yfz,
c     &            dxf(1:n,1:npla),dyf(1:n,1:npla))
c         pb(iv) = p(iv)-es
c         call mrqposfit(n,pla(l)%tas(i),pb(1:n),xfq,yfq,
c     &            dxf(1:n,1:npla),dyf(1:n,1:npla))

              CALL MRQPOSFIT(N,PLA(L)%TAS(I),P(1:N),
     &                      XF(1:NPLA),YF(1:NPLA),
     &                      DXF(1:N,1:NPLA),DYF(1:N,1:NPLA))
c              WRITE(SD,*)iv,'x', ! sngl(xfz),sngl(xfq),
c     &         sngl((xfz-xfq)/(2.*es)),
c     &                 sngl(dxf(iv,:))
c              WRITE(SD,*)iv,'y', ! sngl(yfz),sngl(yfq),
c     &         sngl((yfz-yfq)/(2.*es)),
c     &                 sngl(dyf(iv,:))
c         end do
c        stop

c         print*,sngl(pla(l)%x(i)/star%dist*1d3),
c     &            sngl(pla(l)%y(i)/star%dist*1d3),
c     &            sngl(xf(l)/star%dist*1d3),
c     &            sngl(yf(l)/star%dist*1d3)    
             
            SXM1 = PLA(L)%SIGXM1(I)
            SYM1 = PLA(L)%SIGYM1(I)
            RHO = PLA(L)%RHOXY(I)
            UF = PLA(L)%UFAC(I)
            DX = (PLA(L)%X(I)-XF(L))*SXM1
            DY = (PLA(L)%Y(I)-YF(L))*SYM1       
            CHI2 = CHI2 + (DX*DX+DY*DY-2.d0*RHO*DX*DY)*UF
            F1X = (DX-RHO*DY)*SXM1*UF  
            F1Y = (DY-RHO*DX)*SYM1*UF
            F2XX = PLA(L)%SIGXM2(I)*UF
            F2YY = PLA(L)%SIGYM2(I)*UF
            F2XY = -RHO*SXM1*SYM1*UF
            BETA(1:N) = BETA(1:N)+F1X*DXF(1:N,L)+F1Y*DYF(1:N,L)
            DO J=1,N
                ALPHA(J:N,J) = ALPHA(J:N,J)
     &             +F2XX*DXF(J,L)*DXF(J:N,L)+F2YY*DYF(J,L)*DYF(J:N,L)
     &             +F2XY*(DXF(J,L)*DYF(J:N,L)+DXF(J:N,L)*DYF(J,L))
            END DO        
          END DO
        END DO
c      WRITE(SD,*)chi2,'r'
c        stop

        DO J=2,N
          ALPHA(1:J-1,J) = ALPHA(J,1:J-1)
        END DO 

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
C  Computation of Relative RV contribution to Levenberg-Marquardt coefficients
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQRADVELP(N,P,ALPHA,BETA,CHI2)

        USE DATA

        IMPLICIT NONE

        INTEGER*4       N               ! Dimension of parameter space
        
        REAL*8, DIMENSION(NPAR,NPLA) :: DVF
        REAL*8, DIMENSION(NPLA) :: VF
        REAL*8 ::       ALPHA(N,N),     ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  P(N),           ! The parameters
     &                  CHI2,           ! Chi2
     &                  F1,F2,          ! Intermediaires
     &                  DV              ! Normalized difference   
        INTEGER*4       I,J,K,L,SJ,iv,jv  

c        REAL*8 ::       PB(NPAR),vfz(npla),vfq(npla),es
        
c        pb(1:n) = p(1:n)
        DVF = 0.d0
c        es = 1.d-8
        DO L = 1,NPLA
           DO I = 1,PLA(L)%NDATVR
c         do iv=1,n
c         pb(iv) = p(iv)+es
c         call mrqvfitp(n,pla(l)%tvr(i),pb(1:n),vfz,dvf(1:n,1:npla))
c         pb(iv) = p(iv)-es
c         call mrqvfitp(n,pla(l)%tvr(i),pb(1:n),vfq,dvq(1:n,1:npla))

              CALL MRQVFITP(N,PLA(L)%TVR(I),P(1:N),
     &                      VF(1:NPLA),DVF(1:N,1:NPLA))

c     WRITE(SD,*)iv,'v', ! sngl(vfz),sngl(vfq),
c     &         sngl((vfz-vfq)/(2.*es)),sngl(dvf(iv,:))
c         end do
c        stop

              DV = PLA(L)%V(I)-VF(L)
              CHI2 = CHI2 + DV*DV*PLA(L)%SIGVM2(I)          
              F1 = DV*PLA(L)%SIGVM2(I)
              BETA(1:N) = BETA(1:N)+F1*DVF(1:N,L)
              DO J = 1,N
                 F2 = DVF(J,L)*PLA(L)%SIGVM2(I)
                 ALPHA(J:N,J) = ALPHA(J:N,J)+F2*DVF(J:N,L)
              END DO        
           END DO
        END DO
c      WRITE(SD,*)chi2,'r'
c        stop

        DO J = 2,N
          ALPHA(1:J-1,J) = ALPHA(J,1:J-1)
        END DO 

        END 
C
C     
C -----------------------------------------------------------------------------
C       Calculation of Chi2 and coefficients for Levenberg-Marquardt
C -----------------------------------------------------------------------------
C

        SUBROUTINE MRQCOF(N,P,ALPHA,BETA,CHI2,FMAP)

        USE DATA

        IMPLICIT NONE

        INTEGER*4       N               ! Dimension of parameter space

        REAL*8, DIMENSION(NPAR,NPLA) :: DXF,DYF
        REAL*8, DIMENSION(NPAR) :: DVF
        REAL*8, DIMENSION(NPLA) :: XF,YF
        REAL*8, DIMENSION(NPAR,NPAR) :: D2FMAP  ! Second derivatives / MAP
        REAL*8, DIMENSION(NPAR) :: DFMAP        ! First derivatives /MAP
        REAL*8 ::       ALPHA(N,N), ! Matrix coefficients
     &                  BETA(N),        ! Vector
     &                  P(N),           ! The parameters
     &                  CHI2,           ! Chi2
     &                  FMAP,           ! Map merit function (added to Chi2)
     &                  F1X,F1Y,F2XX,F2YY,F2XY,F1,F2, ! Intermediaires
     &                  SIGJV,          ! Velocity Jitter
     &                  SIGV2P,            ! Incertitude^2 augmentée
     &                  FS,
     &                  SXM1,SYM1,RHO,UF,  !
     &                  VF,                ! Fitted velocity
     &                  DX,DY,DV           ! Ecarts normalisés   
        
        
        REAL*8 ::       PB(NPAR),es,
     &                  xfz(npla),yfz(npla),xfq(npla),yfq(npla), 
     &                  fmappp,fmapmp,fmappm,fmapmm,
     &                  vfq,vfz,fmapz,fmapq ! Vitesse radiale
        INTEGER*4 ::     I,J,K,L,SJ,iv,jv          
        LOGICAL ::       JITTER
c     
        CHI2 = 0.d0       
        ALPHA(1:N,1:N) = 0.d0
        BETA(1:N) = 0.d0
        JITTER = ISDATA(2).AND.(JITNUM.EQ.1).AND.(N.EQ.NPAR)
        IF (ISDATA(2)) CALL MRQRADVEL(N,JITTER,P(1:N),
     &                     ALPHA(1:N,1:N),BETA(1:N),CHI2)
        IF (ISDATA(3)) CALL MRQRADVELP(N,P(1:N),
     &                     ALPHA(1:N,1:N),BETA(1:N),CHI2)
        
        CALL MRQASTROM(N,P,ALPHA(1:N,1:N),BETA(1:N),CHI2)

c...  Map merit function & derivatives
c        pb(1:n) = p(1:n)
c           es=1d-6
c           do iv = 1,n
c             do jv = 1,n
c                pb = p
c                if (jv.ne.iv) then
c                  pb(iv)=p(iv)+es
c                  pb(jv)=p(jv)+es
c           call mrqmap(jitter,n,pb,fmappp,dfmap(1:n),d2fmap(1:n,1:n))
c                  pb(iv) = p(iv)-es
c                  pb(jv) = p(jv)+es
c           call mrqmap(jitter,n,pb,fmapmp,dfmap(1:n),d2fmap(1:n,1:n))
c                  pb(iv)=p(iv)+es
c                  pb(jv)=p(jv)-es
c           call mrqmap(jitter,n,pb,fmappm,dfmap(1:n),d2fmap(1:n,1:n))
c                  pb(iv) = p(iv)-es
c                  pb(jv) = p(jv)-es
c           call mrqmap(jitter,n,pb,fmapmm,dfmap(1:n),d2fmap(1:n,1:n))
c                  CALL MRQMAP(JITTER,N,P(1:N),FMAP,DFMAP(1:N),
c     &                    D2FMAP(1:N,1:N))
c                if (d2fmap(iv,jv).ne.0.d0) write(sd,*)iv,jv,
c     &               sngl(0.25d0*(fmappp-fmapmp-fmappm+fmapmm)/(es*es)),
c     &               sngl(d2fmap(iv,jv))
c             else
c               pb(iv) =p(iv)+es 
c               call mrqmap(jitter,n,pb,fmapz,dfmap(1:n),d2fmap(1:n,1:n))
c               pb(iv)=p(iv)-es
c               call mrqmap(jitter,n,pb,fmapq,dfmap(1:n),d2fmap(1:n,1:n))
c               call mrqmap(jitter,n,p(1:n),fmap,dfmap(1:n),
c     &               d2fmap(1:N,1:N))
c               if (d2fmap(iv,iv).ne.0.d0) write(sd,*)iv,iv,
c     &              sngl((fmapq+fmapz-2.d0*fmap)/(es*es)),
c     &               sngl(d2fmap(iv,iv))               
c             end if
c            end do
c          end do
c            stop
        CALL MRQMAP(JITTER,N,P(1:N),FMAP,DFMAP(1:N),D2FMAP(1:N,1:N))
        CHI2 = CHI2+FMAP        ! Her we take into account MAP in Chi2
        BETA(1:N) = BETA(1:N)-0.5d0*DFMAP(1:N) 
        ALPHA(1:N,1:N) = ALPHA(1:N,1:N)+0.5d0*D2FMAP(1:N,1:N)
        END    

C
C     
C -----------------------------------------------------------------------------
C    Calculation of FMAP & derivatives
C        MAP = exp(-Chi2/2)*f(parameters) (f = priors)
C     => Replace Chi2 with Chi2-2*ln(f(parameters)) = -2*ln(MAP) = Chi2+FMAP    
C-----------------------------------------------------------------------------
C

        SUBROUTINE MRQMAP(JITTER,N,P,FMAP,DFMAP,D2FMAP)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 :: N                  ! Dimension of parameter space
        REAL*8, DIMENSION(N,N) :: D2FMAP ! Second derivatives
        REAL*8, DIMENSION(N) :: DFMAP   ! First derivatives
        REAL*8, DIMENSION(N) :: P       ! Parameters
        REAL*8 ::       FMAP,           ! Merit function
     &                  SIGJV,SIGJV2,   ! Jitter,Jitter^2
     &                  SIGV2P          ! Jitter^2+Error_i^2

        INTEGER*4 ::    I,SJ          
        LOGICAL ::      JITTER
c     
        IF (JITTER) THEN
          SJ = NEL*NPLA+2
          SIGJV = EXP(P(SJ))        ! Velocity Jitter
          SIGJV2 = SIGJV*SIGJV
        END IF

        FMAP = 0.d0
        DFMAP = 0.d0
        D2FMAP = 0.d0

        CALL MRQMAPASTROM(N,P,FMAP,DFMAP(1:N),D2FMAP(1:N,1:N))
             
c... Taking into account prior on Stellar Jitter in FMAP (+2*ln(j+j0))
        IF (JITTER) THEN
          FMAP = FMAP+2.d0*LOG(SIGJV+SIGV0)
          DFMAP(SJ) = DFMAP(SJ)+2.d0*SIGJV/(SIGJV+SIGV0)
          D2FMAP(SJ,SJ) = D2FMAP(SJ,SJ)
     &                    +2.d0*SIGJV*SIGV0/(SIGJV+SIGV0)**2
c...  Taking into account effect of Stellar Jitter in p(data|model)
c...        = Product(cte*1/sqrt(j^2+sig_i^2)) => +2*1/2*ln(j^2+sig_i^2)        
          DO I = 1,STAR%NDATVR
            SIGV2P = STAR%SIGV2(I)+SIGJV2
            FMAP = FMAP+LOG(1.d0+SIGJV2/STAR%SIGV2(I))
c            FMAP = FMAP+LOG(SIGJV2+STAR%SIGV2(I))
            DFMAP(SJ) = DFMAP(SJ)+2.d0*SIGJV2/SIGV2P  !!!!!!
            D2FMAP(SJ,SJ) = D2FMAP(SJ,SJ)
     &                  +4.d0*STAR%SIGV2(I)*SIGJV2/(SIGV2P*SIGV2P)
          END DO
        END IF

        END


