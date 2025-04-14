#! /Users/lacquema/Oracle.env/bin/python3

### --- Packages --- ###

# Standard packages
import sys

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QStatusBar, QWidget, QApplication, QLabel, QGridLayout
from PyQt6.QtCore import pyqtSignal

# My packages
from Tools import Delimiter, SpaceView, TempoView, Conv, Hist, Hist2D, Corner, PosAtDate
from TransferData import TransfertSimu
from BestOrbits import BestOrbitsClass
from SelectOrbits import SelectOrbitsClass

### --- Main Window Generating --- ###

class WindowMainClass(QMainWindow):
    SignalCloseWindowMain = pyqtSignal()

    def __init__(self, SystName, PathSimu, NbSelectOrbits, NbPtsEllipse=1000, SystDist=None):
        super().__init__()

        # Window settings
        self.setWindowTitle(f'Oracle analysis of {SystName}')

        # Layout initialization
        Layout = QVBoxLayout()

        # Data
        InputData, OutputParams = TransfertSimu(PathSimu)
        if InputData:
            NbInputData = InputData['NbData']
            SystDist = InputData['SystDist']['pc']
        else:
            NbInputData = 1
            if SystDist is None:
                print('No distance provided, please provide a distance in parsec')
                return

        # Best orbits
        BestOrbits = BestOrbitsClass(*OutputParams, NbPtsEllipse, SystDist, NbInputData)
        Layout.addWidget(BestOrbits.Widget)

        # Select orbits
        SelectOrbits = SelectOrbitsClass(*OutputParams, NbSelectOrbits, NbPtsEllipse, SystDist)

        # Separation
        Layout.addWidget(Delimiter(Title='Plots :'))

        # Grid layout
        GridLayout = QGridLayout()

        # Space view
        GridLayout.addWidget(SpaceView(InputData, SelectOrbits.SelectEllipses, BestOrbits.BestEllipses), 0, 0, 1, 3)

        # Temporal study
        GridLayout.addWidget(TempoView(InputData, SelectOrbits.SelectEllipses, BestOrbits.BestEllipses), 1, 0)

        # Convergence of orbit parameters
        GridLayout.addWidget(Conv(OutputParams), 1, 1)

        # Histogram of orbit parameters
        GridLayout.addWidget(Hist(OutputParams, BestOrbits.BestParams), 1, 2)

        # Histogram 2D of orbit parameters
        GridLayout.addWidget(Hist2D(OutputParams, BestOrbits.BestParams), 2, 0)

        # Corner plot
        GridLayout.addWidget(Corner(SelectOrbits.SelectParams, BestOrbits.BestParams), 2, 1)

        # Position at date
        GridLayout.addWidget(PosAtDate(InputData, SelectOrbits.SelectEllipses, BestOrbits.BestEllipses), 2, 2)

        # Add grid layout to layout
        Layout.addLayout(GridLayout)

        # Widget container
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

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
    Simu = '/Users/lacquema/Simulations/ggtau/ggtau_Ab2Aa_Ab1_fit/ggtau_Ab2Aa_Ab1_fit_Herve_3f/solggtaua_3f.dat'
    # PathInputData = Simu + 'ggtau_Ab12_hci.dat'
    # PathOutputData = Simu + 'solggtauAb12.dat'
    # PathInputData = Simu + 'ggtau_data.txt'
    # PathOutputData = Simu + 'solggtaua_3f.dat'
    SystDist = 145
    SystName = 'ggtau'

    WindowMain = WindowMainClass(SystName=SystName, PathSimu=Simu, NbSelectOrbits=10000, NbPtsEllipse=500, SystDist=None)  # Main window showing
    sys.exit(app.exec())  # Application execution

