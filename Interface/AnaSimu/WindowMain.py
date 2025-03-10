#! /Users/lacquema/Oracle.env/bin/python3

### --- Packages --- ###

# Standard packages
import sys

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QStatusBar, QWidget, QApplication
from PyQt6.QtCore import pyqtSignal

# My packages
from Tools import Delimiter, SpaceView, TempoView, Conv, Hist, Hist2D, Corner, PosAtDate
from TransferData import TransfertData, TransfertSimu
from BestOrbits import BestOrbitsClass
from SelectOrbits import SelectOrbitsClass

### --- Main Window Generating --- ###

class WindowMainClass(QMainWindow):
    SignalCloseWindowMain = pyqtSignal()

    def __init__(self, SystName, PathData, PathSimu, SystDist, NbSelectOrbits, NbPtsEllipse):
        super().__init__()

        # Window settings
        self.setWindowTitle(f'Oracle data analysis of {SystName}')

        # Layout initialization
        Layout = QVBoxLayout()

        # Data
        InputData = TransfertData(PathData)[1]
        OutputParams = TransfertSimu(PathSimu)

        # Best orbits
        BestOrbits = BestOrbitsClass(*OutputParams, NbPtsEllipse, SystDist, InputData[0], InputData[-2])
        Layout.addWidget(BestOrbits.Widget)

        # Select orbits
        SelectOrbits = SelectOrbitsClass(*OutputParams, NbSelectOrbits, NbPtsEllipse, SystDist)

        # Separation
        Layout.addWidget(Delimiter())

        # Space view
        Layout.addWidget(SpaceView(InputData, SelectOrbits.SelectEllipses, BestOrbits.BestEllipses))

        # Temporal study
        Layout.addWidget(TempoView(InputData, SelectOrbits.SelectEllipses, BestOrbits.BestEllipses))

        # Convergence of orbit parameters
        Layout.addWidget(Conv(OutputParams))

        # Histogram of orbit parameters
        Layout.addWidget(Hist(OutputParams, BestOrbits.BestParams))

        # Histogram 2D of orbit parameters
        Layout.addWidget(Hist2D(OutputParams, BestOrbits.BestParams))

        # Corner plot
        Layout.addWidget(Corner(SelectOrbits.SelectParams, BestOrbits.BestParams))

        # Position at date
        Layout.addWidget(PosAtDate(InputData, SelectOrbits.SelectEllipses, BestOrbits.BestEllipses))

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

    Simu = '/Users/lacquema/Simulations/ggtau/ggtau_Ab2_Ab1_fit/ggtau_Ab2_Ab1_fit_1/'
    PathInputData = Simu + 'ggtau_Ab12_hci.dat'
    PathOutputData = Simu + 'solggtauAb12.dat'
    SystDist = 145
    SystName = 'ggtau'

    WindowMain = WindowMainClass(SystName, PathInputData, PathOutputData, SystDist, 10000, 500)  # Main window showing
    sys.exit(app.exec())  # Application execution

