
### --- Packages --- ###
from PyQt6.QtWidgets import QWidget
import numpy as np
import Utils as ut
import random as rd


class SelectOrbitsClass(QWidget):
    def __init__(self, NbBodies, NbOrbits, P, a, e, i, w, W, tp, m, Mdyn, Chi2, map, NbSelectOrbits, NbPtsEllipse, StarDist):
        super().__init__()

        # Random selection of orbits
        self.NbParams = 10
        self.SelectP, self.Selecta, self.Selecte, self.Selecti, self.Selectw, self.SelectW, self.Selecttp, self.Selectm, self.SelectMdyn, self.SelectChi2 = [np.zeros((NbBodies, NbSelectOrbits)) for k in range(self.NbParams)]

        for j in range(NbBodies):
            for k in range(NbSelectOrbits):
                indexRd = rd.randint(0, NbOrbits-1)
                self.SelectP[j][k] = P[j][indexRd]
                self.Selecta[j][k] = a[j][indexRd]
                self.Selecte[j][k] = e[j][indexRd]
                self.Selecti[j][k] = i[j][indexRd]
                self.Selectw[j][k] = w[j][indexRd]
                self.SelectW[j][k] = W[j][indexRd]
                self.Selecttp[j][k] =  tp[j][indexRd]
                self.Selectm[j][k] = m[j][indexRd]
                self.SelectMdyn[j][k] = Mdyn[j][indexRd]
                self.SelectChi2[j][k] = Chi2[j][indexRd]

        self.SelectParams = [NbBodies, NbSelectOrbits, self.SelectP, self.Selecta, self.Selecte, self.Selecti, self.Selectw, self.SelectW, self.Selecttp, self.Selectm, self.SelectMdyn, self.SelectChi2]

        # Ellipse
        self.Selectt, self.SelectX, self.SelectY, self.SelectZ = [np.zeros((NbBodies, NbSelectOrbits, NbPtsEllipse)) for k in range(4)]
        for j in range(NbBodies):
            for k in range(NbSelectOrbits):
                self.Selectt[j][k], self.SelectX[j][k], self.SelectY[j][k], self.SelectZ[j][k] = ut.Ellipse(self.SelectP[j][k], self.Selecta[j][k], self.Selecte[j][k], self.Selecti[j][k], self.Selectw[j][k], self.SelectW[j][k], self.Selecttp[j][k], NbPtsEllipse, Time=True)

        self.SelectEllipses = [NbBodies, NbSelectOrbits, NbPtsEllipse, self.SelectP, self.Selectt, -self.SelectX/StarDist*1000, self.SelectY/StarDist*1000, self.SelectZ/StarDist*1000]
