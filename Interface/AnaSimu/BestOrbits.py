### --- Packages --- ###
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem
import numpy as np
import Utils as ut


class BestOrbitsClass(QWidget):
    def __init__(self, NbBodies, NbOrbits, P, a, e, i, w, W, tp, m, Mdyn, Chi2, map, NbPtsEllipse, StarDist, NbInputData, Corr):
        super().__init__()

        # Best fit for each bodies
        self.NbParams = 10
        self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.BestMdyn, self.BestChi2 = [np.zeros(NbBodies) for k in range(self.NbParams)]

        self.NbParamsLib = 0
        self.Params = [P, a, e, i, w, W, tp, m, Mdyn]

        for j in range(NbBodies):
            self.BestChi2[j] = np.min(Chi2[j])
            IndexBestChi2 = list(Chi2[j]).index(self.BestChi2)
            self.BestP[j] = P[j][IndexBestChi2]
            self.Besta[j] = a[j][IndexBestChi2]
            self.Beste[j] = e[j][IndexBestChi2]
            self.Besti[j] = i[j][IndexBestChi2]
            self.Bestw[j] = w[j][IndexBestChi2]
            self.BestW[j] = W[j][IndexBestChi2]
            self.Besttp[j] =  tp[j][IndexBestChi2]
            self.Bestm[j] = m[j][IndexBestChi2]
            self.BestMdyn[j] = Mdyn[j][IndexBestChi2]

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

        for j in range(NbBodies):
            for k in range(self.NbParams):
                TblBestFit.setItem(j, k, QTableWidgetItem('{}'.format(np.around(self.BestParams[k+1][j], 3))))

        Layout.addWidget(TblBestFit)
        self.Widget.setLayout(Layout)

        # Ellipse
        self.Bestt, self.BestX, self.BestY, self.BestZ = [np.zeros((NbBodies, NbPtsEllipse)) for k in range(4)]
        for j in range(NbBodies):
            self.Bestt[j], self.BestX[j], self.BestY[j], self.BestZ[j] = ut.Ellipse(self.BestP[j], self.Besta[j], self.Beste[j], self.Besti[j], self.Bestw[j], self.BestW[j], self.Besttp[j], NbPtsEllipse, Time=True)

        self.BestEllipses = [NbBodies, NbPtsEllipse, self.BestP, self.Bestt, -self.BestX/StarDist*1000, self.BestY/StarDist*1000, self.BestZ/StarDist*1000]

