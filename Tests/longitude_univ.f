c-------------------------------------------------------------------------
c     Programme de test, variables universelles : routine vide, car ces
c     variables n'utilisent pas de longitude a T0.
c-------------------------------------------------------------------------
        SUBROUTINE LONGITUDE_T0(I,CL,SL)

        IMPLICIT NONE

        INTEGER*4 ::    I               ! Index of the orbit
        REAL*8 ::       CL,SL           ! cos,sin(longitude at T0)

        CL = 0.d0
        SL = 0.d0

        END
