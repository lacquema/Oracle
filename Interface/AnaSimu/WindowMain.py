#! /Users/lacquema/Oracle.env/bin/python3

### --- Packages --- ###

# Standard packages
import sys

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QStatusBar, QWidget, QApplication, QLabel, QGridLayout
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QGuiApplication, QCursor

# My packages
from Tools import Delimiter, SpaceView, TempoView, Conv, Hist, Hist2D, Corner, PosAtDate
from LMOrbit import LMOrbitClass
from BestOrbit import BestOrbitClass
from SelectOrbits import SelectOrbitsClass
from TransferData import TransfertSimu

### --- Main Window Generating --- ###

class WindowMainClass(QMainWindow):
    SignalCloseWindowMain = pyqtSignal()

    def __init__(self, SimuName, InputData, OutputParams, NbSelectOrbits, SystDist, NbPtsEllipse=1000):
        super().__init__()

        # Window settings
        self.setWindowTitle(f'Oracle analysis of {SimuName}')
        # self.setWindowTitle(f'Oracle analysis of ...')

        # Move the window to the top-left of the active screen (where the mouse is)
        mouse_pos = QCursor.pos()
        for screen in QGuiApplication.screens():
            if screen.geometry().contains(mouse_pos):
                geo = screen.geometry()
                self.move(geo.x(), geo.y())

        # Layout initialization
        Layout = QVBoxLayout()

        # Data
        if InputData != None:
            NbInputData = InputData['NbData']
        else:
            NbInputData = 1

        # LM orbit
        LMOrbit = LMOrbitClass(*OutputParams, NbPtsEllipse, SystDist, NbInputData)
        Layout.addWidget(LMOrbit.Widget)

        # Best orbits
        BestOrbit = BestOrbitClass(*OutputParams, NbPtsEllipse, SystDist, NbInputData)
        Layout.addWidget(BestOrbit.Widget)
        # print("Best orbits done")

        # Select orbits
        SelectOrbits = SelectOrbitsClass(*OutputParams, NbSelectOrbits, NbPtsEllipse, SystDist)
        # print("Select orbits done")

        # Separation
        Layout.addWidget(Delimiter())

        # Grid layout
        GridLayout = QGridLayout()

        # Space view
        GridLayout.addWidget(SpaceView(InputData, SelectOrbits.SelectEllipses, BestOrbit.BestEllipse, LMOrbit.LMEllipse), 0, 0, 1, 3)
        # print("Space view done")

        # Temporal study
        GridLayout.addWidget(TempoView(InputData, SelectOrbits.SelectEllipses, BestOrbit.BestEllipse, LMOrbit.LMEllipse), 1, 0)
        # print("Temporal view done")

        # Convergence of orbit parameters
        GridLayout.addWidget(Conv(OutputParams), 1, 1)
        # print("Convergence done")

        # Histogram of orbit parameters
        GridLayout.addWidget(Hist(OutputParams, BestOrbit.BestParams, LMOrbit.LMParams), 1, 2)
        # print("Histogram done")

        # Histogram 2D of orbit parameters
        GridLayout.addWidget(Hist2D(OutputParams, BestOrbit.BestParams, LMOrbit.LMParams), 2, 0)
        # print("2D Histogram done")

        # Corner plot
        GridLayout.addWidget(Corner(SelectOrbits.SelectParams, BestOrbit.BestParams, LMOrbit.LMParams), 2, 1)
        # print("Corner plot done")

        # Position at date
        GridLayout.addWidget(PosAtDate(InputData, OutputParams, SelectOrbits.SelectEllipses, BestOrbit.BestParams, BestOrbit.BestEllipse, LMOrbit.LMParams, LMOrbit.LMEllipse, SystDist), 2, 2)
        # print("Position at date done")

        # Add grid layout to layout
        Layout.addLayout(GridLayout)

        # Widget container
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

        # Initial window width: ensure both top tables are fully visible
        top_tables_width = max(LMOrbit.WidgetMinWidth, BestOrbit.WidgetMinWidth)
        layout_margins = Layout.contentsMargins()
        required_width = top_tables_width + layout_margins.left() + layout_margins.right() + 30
        self.setMinimumWidth(required_width)
        self.resize(required_width, self.sizeHint().height())

        # Status bar
        self.setStatusBar(QStatusBar(self))
        
        # Showing
        self.show()

    # Close programme when the main window is closed
    def closeEvent(self, e):
        try:
            app.closeAllWindows()
        except:
            self.SignalCloseWindowMain.emit()

if __name__ == "__main__":
    app = QApplication(sys.argv)  # Application creation

    # Simu = '/Users/lacquema/Simulations/ggtau/ggtau_Ab2_Ab1_fit/ggtau_Ab2_Ab1_fit_1/solggtauAb12.dat'
    SimuPath = '/Users/lacquema/Simulations/ggtau/ggtau_Ab2Aa_Ab1_fit/ggtau_Ab2Aa_Ab1_fit_Herve_3e/solggtaua_3e.dat'
    SimuName = SimuPath.split('/')[-1]
    # PathInputData = Simu + 'ggtau_Ab12_hci.dat'
    # PathOutputData = Simu + 'solggtauAb12.dat'
    # PathInputData = Simu + 'ggtau_data.txt'
    # PathOutputData = Simu + 'solggtaua_3f.dat'
    InputData, OutputParams = TransfertSimu(SimuPath, 2)
    SystDist = 145

    WindowMain = WindowMainClass(SimuName, InputData, OutputParams, NbSelectOrbits=10000, SystDist=SystDist)  # Main window showing
    sys.exit(app.exec())  # Application execution

