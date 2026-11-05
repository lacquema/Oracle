c-------------------------------------------------------------------------
c     Programme de test, variables classiques : longitude vraie a T0.
c     TP_TO_LAMBDA n'existe que dans kepellip.f.
c-------------------------------------------------------------------------
        SUBROUTINE LONGITUDE_T0(I,CL,SL)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    I               ! Index of the orbit
        REAL*8 ::       CL,SL           ! cos,sin(longitude at T0)

        CALL TP_TO_LAMBDA(STAR%T0,PLA(I)%TP,DPI/PLA(I)%PER,
     &             PLA(I)%EXC,PLA(I)%EXQ,PLA(I)%CW,PLA(I)%SW,CL,SL)

        END
