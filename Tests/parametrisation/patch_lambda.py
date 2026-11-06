"""Variante : parametre 6,7 = (cos lambda, sin lambda)/Pq, lambda = longitude vraie a T0."""
import sys
def sub(path, old, new, count=1):
    s = open(path, newline='').read()
    nl = '\r\n' if '\r\n' in s else '\n'
    o, n = old.replace('\n', nl), new.replace('\n', nl)
    assert s.count(o) == count, (path, s.count(o), old[:70])
    s = s.replace(o, n)
    open(path, 'w', newline='').write(s)
L, M = sys.argv[1], sys.argv[2]
k = L + '/kepu.f'
# COMPUTE_MAP : 1/Pq^2
sub(k, "          UP2 = EXP(-2.d0*P6)  !  1/Period^2\n", "          UP2 = P6*P6+P7*P7  !  1/Period^2 (variante lambda)\n")
# ELEMENTS
sub(k, "     &       TT                     ! Time\n", "     &       TT,                    ! Time\n     &       CNU,SNU                ! cos,sin(nu) at T0\n        LOGICAL :: OKNU\n")
sub(k, "           UP = EXP(-P6)\n", "           UP = SQRT(P6*P6+P7*P7)\n")
sub(k, "           CALL ECC_FUNC(EXC(I),FE,DFE)\n           S = P7/FE  ! S=Universal variable at T=T0 / P7 = S*f(e)\n",
       "           CNU = (P6*CW(I)+P7*SW(I))/UP ! nu = lambda-w\n           SNU = (P7*CW(I)-P6*SW(I))/UP\n           CALL NU_TO_S(CNU,SNU,EXC(I),Q(I),MUF(I),S,OKNU)\n")
# RATIO_PRIORS
sub(k, "          UPQ = EXP(-P6)\n          UPQT = EXP(-PT6)\n", "          UPQ = SQRT(P6*P6+P7*P7)\n          UPQT = SQRT(PT6*PT6+PT7*PT7)\n")
sub(k, "     &       UPQ,UPQT           ! Inverses de périodes\n", "     &       UPQ,UPQT,          ! Inverses de périodes\n     &       CNU,CNUT,FNU,FNUT  ! cos(nu) at T0 & Jacobian factors\n")
sub(k, "          CALL ECC_FUNC(EXC,FE,DFE)\n          CALL ECC_FUNC(EXCT,FET,DFET)\n          S = P7/FE     ! P7 = S*f(E)\n          ST = PT7/FET\n", "")
old = s0 = None
s = open(k, newline='').read()
nl = '\r\n' if '\r\n' in s else '\n'
a = s.index("          ALPHA = MTOT(I)/Q*(1.d0-EXC)" + nl + "          XX = S*S*ALPHA")
b = s.index("(EXC/EXCT)*(FE/FET)*(DSDT/DSDTT)", a); b = s.index(nl, b) + len(nl)
new = nl.join([
"c...  Variante lambda : J/prior = e*(1+e*cos(nu))^2/((1+e)^(3/2)*Pq^3)",
"          CNU = (P6*CW+P7*SW)/UPQ",
"          CNUT = (PT6*CWT+PT7*SWT)/UPQT",
"          BAD = BAD.OR.(1.d0+EXCT*CNUT.LE.0.d0)",
"          FNU = (1.d0+EXC*CNU)**2/(1.d0+EXC)**1.5d0",
"          FNUT = (1.d0+EXCT*CNUT)**2/(1.d0+EXCT)**1.5d0",
"          RAPQ = RAPQ*(EXC/EXCT)*(FNU/FNUT)*(UPQ/UPQT)**3", ""])
s = s[:a] + new + s[b:]
s += nl.join(["", "C", "C  Variante lambda : variable universelle s a partir de l'anomalie vraie",
"C", "        SUBROUTINE NU_TO_S(CNU,SNU,EXC,Q,MU,S,OK)", "", "        IMPLICIT NONE", "",
"        REAL*8 ::       CNU,SNU,        ! cos,sin(true anomaly)",
"     &                  EXC,Q,MU,       ! Eccentricity, periastron, mass",
"     &                  S,              ! Universal variable (output)",
"     &                  NU,ALPHA,X      ! Internals",
"        LOGICAL ::      OK", "",
"        OK = .TRUE.", "        NU = ATAN2(SNU,CNU)", "        ALPHA = MU/Q*(1.d0-EXC)",
"        IF (EXC.LT.1.d0-1.d-10) THEN",
"           X = 2.d0*ATAN2(SQRT(1.d0-EXC)*SIN(0.5d0*NU),",
"     &                    SQRT(1.d0+EXC)*COS(0.5d0*NU))",
"           S = X/SQRT(ALPHA)",
"        ELSE IF (EXC.GT.1.d0+1.d-10) THEN",
"           OK = (1.d0+EXC*CNU.GT.0.d0)",
"           X = SQRT((EXC-1.d0)/(EXC+1.d0))*TAN(0.5d0*NU)",
"           X = MIN(MAX(X,-1.d0+1.d-15),1.d0-1.d-15)",
"           S = 2.d0*ATANH(X)/SQRT(-ALPHA)",
"        ELSE",
"           S = TAN(0.5d0*NU)*SQRT(2.d0*Q/MU)",
"        END IF", "        END", ""])
open(k, 'w', newline='').write(s)
# Programme principal : PSTART(6,7) depuis tp
m = M + '/astrom_univ_mcmco.f'
sub(m, "        LOGICAL ::          OK,NEW\n", "        LOGICAL ::          OK,NEW\n        REAL*8 :: RX,RY,RR,CL,SL\n")
sub(m, "              PSTART(DKAL+6) = LOG(PLA(I)%PER)         ! Par 6 = ln(perq)\n              PSTART(DKAL+7) = P(DKAL+7) ! Par 7=f(E)*S=f(E)*univ.var@t=T0    \n",
"""c...  Variante lambda : Par 6,7 = (cos,sin)(lambda)/Pq, lambda = w+nu(T0)
              MASS = (DPI*UP)**2*PLA(I)%Q**3
              HH2 = -(PLA(I)%EXC-1.d0)*MASS/PLA(I)%Q
              M0 = STAR%T0-PLA(I)%TP
              CALL KEP_UNIV(M0,PLA(I)%Q,MASS,HH2,S,C0,SC1,S2C2,S3C3)
              RX = PLA(I)%Q-MASS*S2C2
              RY = SQRT(PLA(I)%Q*MASS*(1.d0+PLA(I)%EXC))*SC1
              RR = SQRT(RX*RX+RY*RY)
              CL = (RX*PLA(I)%CW-RY*PLA(I)%SW)/RR
              SL = (RX*PLA(I)%SW+RY*PLA(I)%CW)/RR
              PSTART(DKAL+6) = CL*UP
              PSTART(DKAL+7) = SL*UP
""")
print('ok')
