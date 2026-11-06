        SUBROUTINE LONGITUDE_T0(I,CL,SL)

        USE DATA

        IMPLICIT NONE

        INTEGER*4 ::    I               ! Index of the orbit
        REAL*8 ::       CL,SL,          ! cos,sin(true longitude at T0)
     &                  MASS,HH2,M0,S,C0,SC1,S2C2,S3C3,RX,RY,RR

        MASS = (DPI/PLA(I)%PER)**2*PLA(I)%Q**3
        HH2 = -(PLA(I)%EXC-1.d0)*MASS/PLA(I)%Q
        M0 = STAR%T0-PLA(I)%TP
        CALL KEP_UNIV(M0,PLA(I)%Q,MASS,HH2,S,C0,SC1,S2C2,S3C3)
        RX = PLA(I)%Q-MASS*S2C2
        RY = SQRT(PLA(I)%Q*MASS*(1.d0+PLA(I)%EXC))*SC1
        RR = SQRT(RX*RX+RY*RY)
        CL = (RX*PLA(I)%CW-RY*PLA(I)%SW)/RR
        SL = (RX*PLA(I)%SW+RY*PLA(I)%CW)/RR

        END
