### --- Packages --- ###
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem, QHeaderView, QSizePolicy
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFontMetrics
import numpy as np
import Utils as ut

class LMOrbitClass(QWidget):
    def __init__(self, NbBodies, NbOrbits, PlanetsMassUnit, P, a, e, i, w, W, tp, m, m0, V0, Jitter, Chi2, map, NbPtsEllipse, StarDist, NbInputData):
        super().__init__()

        # Number of parameters
        self.NbParams = 10
        # Initialize LM fit parameters for each body
        self.LMP, self.LMa, self.LMe, self.LMi, self.LMw, self.LMW, self.LMtp, self.LMm, self.LMm0, self.LMChi2 = [np.zeros(NbBodies) for _ in range(self.NbParams)]

        # Store input parameters
        self.Params = [P, a, e, i, w, W, tp, m, m0, Chi2]

        self.LMParams = [NbBodies, PlanetsMassUnit, self.LMP, self.LMa, self.LMe, self.LMi, self.LMw, self.LMW, self.LMtp, self.LMm, self.LMm0, self.LMChi2]

        for j in range(NbBodies):
            for lm_values, param_values in zip(self.LMParams[2:], self.Params):
                lm_values[j] = param_values[j][0]

        # Widget container
        self.Widget = QWidget()
        Layout = QVBoxLayout()
    
        # Label
        LblLMFit = QLabel('Levenberg-Marquardt fit for each orbit :')
        LblLMFit.setStatusTip('Orbits parameters corresponding to the LM fit')
        Layout.addWidget(LblLMFit)

        # Table
        TblLMFit = QTableWidget()
        TblLMFit.setStatusTip('Orbits parameters corresponding to the LM fit')
        TblLMFit.setRowCount(NbBodies)
        TblLMFit.setColumnCount(self.NbParams)
        self.LabelParams = ['P [yr]', 'a [AU]', 'e', 'i [°]', 'w [°]', 'W [°]', 'tp [MJD]', 'm ['+PlanetsMassUnit+']', 'm0 [Msun]', 'Chi2']
        TblLMFit.setHorizontalHeaderLabels(self.LabelParams)

        # Disable selection
        TblLMFit.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        # Disable editing
        TblLMFit.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        # Disable focus (no keyboard/mouse interaction)
        TblLMFit.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        # Set column width mode
        TblLMFit.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)

        # Add a stylesheet for better aesthetics with extra padding
        TblLMFit.setStyleSheet("""
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
        TblLMFit.verticalHeader().setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)

        # Keep a fixed, consistent dimension policy across analysis tables
        TblLMFit.setSizeAdjustPolicy(QTableWidget.SizeAdjustPolicy.AdjustIgnored)
        TblLMFit.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        TblLMFit.horizontalHeader().setStretchLastSection(True)

        # Fill table with LM fit parameters and center align values
        for j in range(NbBodies):
            for k in range(self.NbParams):
                # Check if variance is zero for the parameter across all bodies
                if np.var(self.Params[k][j]) == 0 or self.Params[k][j][0] == float('inf'):
                    value = "/"
                else:
                    if k == 6:  # tp [MJD]
                        value = '{:.0f}'.format(np.around(self.LMParams[k+2][j], 0))
                    else:
                        value = '{}'.format(np.around(self.LMParams[k+2][j], 3))
                # value = '...'
                item = QTableWidgetItem(value)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)  # Center align text
                TblLMFit.setItem(j, k, item)

        # Ensure initial width shows all columns (max of header and values)
        table_width = (
            TblLMFit.verticalHeader().width()
            + 2 * TblLMFit.frameWidth()
            + TblLMFit.verticalScrollBar().sizeHint().width()
        )
        cell_metrics = QFontMetrics(TblLMFit.font())
        header_metrics = QFontMetrics(TblLMFit.horizontalHeader().font())
        for col in range(TblLMFit.columnCount()):
            value_width = 0
            for row in range(TblLMFit.rowCount()):
                item = TblLMFit.item(row, col)
                if item is None:
                    continue
                value_width = max(value_width, cell_metrics.horizontalAdvance(item.text()))
            header_item = TblLMFit.horizontalHeaderItem(col)
            header_text = header_item.text() if header_item is not None else ''
            header_width = header_metrics.horizontalAdvance(header_text)
            col_width = max(value_width, header_width)
            col_width += 24
            table_width += col_width
        self.TableMinWidth = table_width
        TblLMFit.setMinimumWidth(table_width)

        # Limit table height to header + filled rows only
        TblLMFit.resizeRowsToContents()
        table_height = TblLMFit.horizontalHeader().height() + 2 * TblLMFit.frameWidth()
        for row in range(TblLMFit.rowCount()):
            table_height += TblLMFit.rowHeight(row)
        TblLMFit.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        TblLMFit.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        TblLMFit.setFixedHeight(table_height)

        Layout.addWidget(TblLMFit)
        self.Widget.setLayout(Layout)
        widget_margins = Layout.contentsMargins()
        self.WidgetMinWidth = self.TableMinWidth + widget_margins.left() + widget_margins.right()
        self.Widget.setMinimumWidth(self.WidgetMinWidth)

        # Ellipse calculations
        self.LMt, self.LMX, self.LMY, self.LMZ, self.LMRa, self.LMDec, self.LMSep, self.LMPa, self.LMRV = [np.zeros((NbBodies, NbPtsEllipse)) for _ in range(9)]
        for j in range(NbBodies):
            self.LMt[j], self.LMX[j], self.LMY[j], self.LMZ[j] = ut.Ellipse(self.LMP[j], self.LMa[j], self.LMe[j], self.LMi[j], self.LMw[j], self.LMW[j], self.LMtp[j], NbPtsEllipse, Time=True)

            # Radial Velocity
            self.LMRV[j] = np.gradient(self.LMZ[j], self.LMt[j]) * 1.495978707e8 / 86400.0 # in km/s
            
            # Conversion to milliarcseconds
            self.LMRa[j] = -self.LMX[j]/StarDist*1000
            self.LMDec[j] = self.LMY[j]/StarDist*1000
            self.LMZ[j] = self.LMZ[j]/StarDist*1000

            # Separation and position angle
            self.LMSep[j] = np.sqrt(self.LMRa[j]**2+self.LMDec[j]**2)
            self.LMPa[j] = np.rad2deg(np.arctan2(self.LMRa[j], self.LMDec[j]))

        self.LMEllipse = [NbBodies, NbPtsEllipse, self.LMP, self.LMt, self.LMRa, self.LMDec, self.LMZ, self.LMSep, self.LMPa, self.LMRV]
