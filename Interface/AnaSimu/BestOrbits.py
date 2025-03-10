### --- Packages --- ###
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem
import numpy as np
import Utils as ut

class BestOrbitsClass(QWidget):
    def __init__(self, NbBodies, NbOrbits, P, a, e, i, w, W, tp, m, Mdyn, Chi2, map, NbPtsEllipse, StarDist, NbInputData, Corr):
        super().__init__()

        # Number of parameters
        self.NbParams = 10
        # Initialize best fit parameters for each body
        self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.BestMdyn, self.BestChi2 = [np.zeros(NbBodies) for _ in range(self.NbParams)]

        # Store input parameters
        self.Params = [P, a, e, i, w, W, tp, m, Mdyn]

        # Find best fit parameters for each body
        for j in range(NbBodies):
            self.BestChi2[j] = np.min(Chi2[j])
            IndexBestChi2 = list(Chi2[j]).index(self.BestChi2[j])  # Adjust index calculation
            self.BestChi2[j] = self.BestChi2[j] / NbInputData # Reduced Chi2
            self.BestP[j], self.Besta[j], self.Beste[j], self.Besti[j], self.Bestw[j], self.BestW[j], self.Besttp[j], self.Bestm[j], self.BestMdyn[j] = [param[j][IndexBestChi2] for param in self.Params]

        self.BestParams = [NbBodies, self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.BestMdyn, self.BestChi2]

        # Widget container
        self.Widget = QWidget()
        Layout = QVBoxLayout()
    
        # Label
        LblBestFit = QLabel('Best fit for each bodies:')
        LblBestFit.setStatusTip('Orbits parameters corresponding to the best Chi2 fit')
        Layout.addWidget(LblBestFit)

        # Table
        TblBestFit = QTableWidget()
        TblBestFit.setEnabled(False)
        TblBestFit.setStatusTip('The number of bodies is counting from the center of the system outwards')
        TblBestFit.setRowCount(NbBodies)
        TblBestFit.setColumnCount(self.NbParams)
        self.LabelParams = ['P (yr)', 'a (AU)', 'e', 'i (deg)', 'w (deg)', 'W (deg)', 'tp (MJD)', 'm (Mjup)', 'Mdyn (Msol)', 'Chi2']
        TblBestFit.setHorizontalHeaderLabels(self.LabelParams)
        TblBestFit.setFixedSize(102*self.NbParams, 25*(NbBodies+1)+8)

        # Fill table with best fit parameters
        for j in range(NbBodies):
            for k in range(self.NbParams):
                TblBestFit.setItem(j, k, QTableWidgetItem('{}'.format(np.around(self.BestParams[k+1][j], 3))))

        Layout.addWidget(TblBestFit)
        self.Widget.setLayout(Layout)

        # Ellipse calculations
        self.Bestt, self.BestX, self.BestY, self.BestZ, self.BestRa, self.BestDec, self.BestSep, self.BestPa = [np.zeros((NbBodies, NbPtsEllipse)) for _ in range(8)]
        for j in range(NbBodies):
            self.Bestt[j], self.BestX[j], self.BestY[j], self.BestZ[j] = ut.Ellipse(self.BestP[j], self.Besta[j], self.Beste[j], self.Besti[j], self.Bestw[j], self.BestW[j], self.Besttp[j], NbPtsEllipse, Time=True)
            
            # Conversion to milliarcseconds
            self.BestRa[j] = -self.BestX[j]/StarDist*1000
            self.BestDec[j] = self.BestY[j]/StarDist*1000
            self.BestZ[j] = self.BestZ[j]/StarDist*1000

            # Separation and position angle
            self.BestSep[j] = np.sqrt(self.BestRa[j]**2+self.BestDec[j]**2)
            self.BestPa[j] = np.rad2deg(np.arctan2(self.BestRa[j], self.BestDec[j]))
        
        self.BestEllipses = [NbBodies, NbPtsEllipse, self.BestP, self.Bestt, self.BestRa, self.BestDec, self.BestZ, self.BestSep, self.BestPa]

