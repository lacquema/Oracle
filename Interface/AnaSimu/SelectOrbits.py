### --- Packages --- ###
from PyQt6.QtWidgets import QWidget
import numpy as np
import Utils as ut
import random as rd

class SelectOrbitsClass(QWidget):
    def __init__(self, NbBodies, NbOrbits, P, a, e, i, w, W, tp, m, m0, V0, Jitter, Chi2, map, NbSelectOrbits, NbPtsEllipse, StarDist):
        super().__init__()

        # Number of parameters
        self.NbParams = 10
        # Initialize selected orbits parameters
        self.SelectP, self.Selecta, self.Selecte, self.Selecti, self.Selectw, self.SelectW, self.Selecttp, self.Selectm, self.Selectm0, self.SelectChi2 = [np.zeros((NbBodies, NbSelectOrbits)) for _ in range(self.NbParams)]

        # Randomly select orbits for each body
        for j in range(NbBodies):
            for k in range(NbSelectOrbits):
                indexRd = rd.randint(0, NbOrbits-1)
                self.SelectP[j][k], self.Selecta[j][k], self.Selecte[j][k], self.Selecti[j][k], self.Selectw[j][k], self.SelectW[j][k], self.Selecttp[j][k], self.Selectm[j][k], self.Selectm0[j][k], self.SelectChi2[j][k] = [param[j][indexRd] for param in [P, a, e, i, w, W, tp, m, m0, Chi2]]

        self.SelectParams = [NbBodies, NbSelectOrbits, self.SelectP, self.Selecta, self.Selecte, self.Selecti, self.Selectw, self.SelectW, self.Selecttp, self.Selectm, self.Selectm0, self.SelectChi2]

        # Ellipse calculations
        self.Selectt, self.SelectX, self.SelectY, self.SelectZ, self.SelectRa, self.SelectDec, self.SelectSep, self.SelectPa, self.SelectRV = [np.zeros((NbBodies, NbSelectOrbits, NbPtsEllipse)) for _ in range(9)]
        
        for j in range(NbBodies):
            for k in range(NbSelectOrbits):
                self.Selectt[j][k], self.SelectX[j][k], self.SelectY[j][k], self.SelectZ[j][k] = ut.Ellipse(self.SelectP[j][k], self.Selecta[j][k], self.Selecte[j][k], self.Selecti[j][k], self.Selectw[j][k], self.SelectW[j][k], self.Selecttp[j][k], NbPtsEllipse, Time=True)

                # Radial Velocity
                self.SelectRV[j][k] = np.gradient(self.SelectZ[j][k], self.Selectt[j][k]) * 1.495978707e8 / 86400.0 # in km/s
                
                # Conversion to milliarcseconds
                self.SelectRa[j][k] = -self.SelectX[j][k]/StarDist*1000
                self.SelectDec[j][k] = self.SelectY[j][k]/StarDist*1000
                self.SelectZ[j][k] = self.SelectZ[j][k]/StarDist*1000

                # Separation and position angle
                self.SelectSep[j][k] = np.sqrt(self.SelectRa[j][k]**2+self.SelectDec[j][k]**2)
                self.SelectPa[j][k] = np.rad2deg(np.arctan2(self.SelectRa[j][k], self.SelectDec[j][k]))

        self.SelectEllipses = [NbBodies, NbSelectOrbits, NbPtsEllipse, self.SelectP, self.Selectt, self.SelectRa, self.SelectDec, self.SelectZ, self.SelectSep, self.SelectPa, self.SelectRV]

