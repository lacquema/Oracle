### --- Packages --- ###
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem, QHeaderView, QSizePolicy
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFontMetrics
import numpy as np
import Utils as ut

class BestOrbitClass(QWidget):
    def __init__(self, NbBodies, NbOrbits, PlanetsMassUnit, P, a, e, i, w, W, tp, m, m0, V0, Jitter, Chi2, map, NbPtsEllipse, StarDist, NbInputData):
        super().__init__()

        # Number of parameters
        self.NbParams = 10
        # Initialize best fit parameters for each body
        self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.Bestm0, self.BestChi2 = [np.zeros(NbBodies) for _ in range(self.NbParams)]

        # Store input parameters
        self.Params = [P, a, e, i, w, W, tp, m, m0, Chi2]

        # Find best fit parameters for each body
        for j in range(NbBodies):
            self.BestChi2[j] = np.min(Chi2[j])
            IndexBestChi2 = list(Chi2[j]).index(self.BestChi2[j])  # Adjust index calculation
            # self.BestChi2[j] = self.BestChi2[j] / NbInputData # Reduced Chi2
            self.BestChi2[j] = self.BestChi2[j]
            self.BestP[j], self.Besta[j], self.Beste[j], self.Besti[j], self.Bestw[j], self.BestW[j], self.Besttp[j], self.Bestm[j], self.Bestm0[j] = [param[j][IndexBestChi2] for param in self.Params[:-1]]

        self.BestParams = [NbBodies, PlanetsMassUnit, self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.Bestm0, self.BestChi2]

        # Widget container
        self.Widget = QWidget()
        Layout = QVBoxLayout()
    
        # Label
        LblBestFit = QLabel('Best fit for each orbit :')
        LblBestFit.setStatusTip('Orbits parameters corresponding to the best Chi2 fit')
        Layout.addWidget(LblBestFit)

        # Table
        TblBestFit = QTableWidget()
        TblBestFit.setStatusTip('Orbits parameters corresponding to the best Chi2 fit')
        TblBestFit.setRowCount(NbBodies)
        TblBestFit.setColumnCount(self.NbParams)
        self.LabelParams = ['P [yr]', 'a [AU]', 'e', 'i [°]', 'w [°]', 'W [°]', 'tp [MJD]', 'm ['+PlanetsMassUnit+']', 'm0 [Msun]', 'Chi2']
        TblBestFit.setHorizontalHeaderLabels(self.LabelParams)

        # Disable selection
        TblBestFit.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        # Disable editing
        TblBestFit.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        # Disable focus (no keyboard/mouse interaction)
        TblBestFit.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        # Set column width mode
        TblBestFit.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)

        # Add a stylesheet for better aesthetics with extra padding
        TblBestFit.setStyleSheet("""
            QTableWidget {
                border: none;
                font-size: 12px;
            }
            QTableWidget::item {
                padding: 8px;
            }
            QHeaderView::section {
                border: 1px solid;
                font-weight: italic;
                padding: 2px;
            }
        """)

        # Center align row numbers (vertical header)
        TblBestFit.verticalHeader().setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)

        # Keep a fixed, consistent dimension policy across analysis tables
        TblBestFit.setSizeAdjustPolicy(QTableWidget.SizeAdjustPolicy.AdjustIgnored)
        TblBestFit.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        TblBestFit.horizontalHeader().setStretchLastSection(True)

        # Fill table with best fit parameters and center align values
        for j in range(NbBodies):
            for k in range(self.NbParams):
                # Check if variance is zero for the parameter across all bodies
                if np.var(self.Params[k][j]) == 0 or self.Params[k][j][0] == float('inf'):
                    value = "/"
                else:
                    if k == 6:  # tp [MJD]
                        value = '{:.0f}'.format(np.around(self.BestParams[k+2][j], 0))
                    else:
                        value = '{}'.format(np.around(self.BestParams[k+2][j], 3))
                # value = '...'
                item = QTableWidgetItem(value)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)  # Center align text
                TblBestFit.setItem(j, k, item)

        # Ensure initial width shows all columns (max of header and values)
        table_width = (
            TblBestFit.verticalHeader().width()
            + 2 * TblBestFit.frameWidth()
            + TblBestFit.verticalScrollBar().sizeHint().width()
        )
        cell_metrics = QFontMetrics(TblBestFit.font())
        header_metrics = QFontMetrics(TblBestFit.horizontalHeader().font())
        for col in range(TblBestFit.columnCount()):
            value_width = 0
            for row in range(TblBestFit.rowCount()):
                item = TblBestFit.item(row, col)
                if item is None:
                    continue
                value_width = max(value_width, cell_metrics.horizontalAdvance(item.text()))
            header_item = TblBestFit.horizontalHeaderItem(col)
            header_text = header_item.text() if header_item is not None else ''
            header_width = header_metrics.horizontalAdvance(header_text)
            col_width = max(value_width, header_width)
            col_width += 24
            table_width += col_width
        self.TableMinWidth = table_width
        TblBestFit.setMinimumWidth(table_width)

        # Limit table height to header + filled rows only
        TblBestFit.resizeRowsToContents()
        table_height = TblBestFit.horizontalHeader().height() + 2 * TblBestFit.frameWidth()
        for row in range(TblBestFit.rowCount()):
            table_height += TblBestFit.rowHeight(row)
        TblBestFit.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        TblBestFit.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        TblBestFit.setFixedHeight(table_height)

        Layout.addWidget(TblBestFit)
        self.Widget.setLayout(Layout)
        widget_margins = Layout.contentsMargins()
        self.WidgetMinWidth = self.TableMinWidth + widget_margins.left() + widget_margins.right()
        self.Widget.setMinimumWidth(self.WidgetMinWidth)

        # Ellipse calculations
        self.Bestt, self.BestX, self.BestY, self.BestZ, self.BestRa, self.BestDec, self.BestSep, self.BestPa, self.BestRV = [np.zeros((NbBodies, NbPtsEllipse)) for _ in range(9)]
        for j in range(NbBodies):
            self.Bestt[j], self.BestX[j], self.BestY[j], self.BestZ[j] = ut.Ellipse(self.BestP[j], self.Besta[j], self.Beste[j], self.Besti[j], self.Bestw[j], self.BestW[j], self.Besttp[j], NbPtsEllipse, Time=True)

            # Radial Velocity
            self.BestRV[j] = np.gradient(self.BestZ[j], self.Bestt[j]) * 1.495978707e8 / 86400.0 # in km/s
            
            # Conversion to milliarcseconds
            self.BestRa[j] = -self.BestX[j]/StarDist*1000
            self.BestDec[j] = self.BestY[j]/StarDist*1000
            self.BestZ[j] = self.BestZ[j]/StarDist*1000

            # Separation and position angle
            self.BestSep[j] = np.sqrt(self.BestRa[j]**2+self.BestDec[j]**2)
            self.BestPa[j] = np.rad2deg(np.arctan2(self.BestRa[j], self.BestDec[j]))

        self.BestEllipse = [NbBodies, NbPtsEllipse, self.BestP, self.Bestt, self.BestRa, self.BestDec, self.BestZ, self.BestSep, self.BestPa, self.BestRV]
