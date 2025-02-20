#! /Users/lacquema/Oracle.env/bin/python3

### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QStatusBar, QWidget, QApplication, QProgressBar
from PyQt6.QtCore import pyqtSignal

# My packages
from Tools import *
from TransferData import TransfertData, TransfertSimu
# from WindowMenu import LoadWindowClass
from BestOrbits import BestOrbitsClass
from SelectOrbits import SelectOrbitsClass



### --- Main Window Generating --- ###

class WindowMainClass(QMainWindow):

    SignalCloseWindowMain = pyqtSignal()

    def __init__(self, SystName, PathData, PathSimu, SystDist, NbSelectOrbits, NbPtsEllipse):
        super().__init__()

        # Window settings
        self.setWindowTitle('Oracle data analysis of {}'.format(SystName))

        # Layout intialisation
        Layout = QVBoxLayout()

        # Data
        InputData = TransfertData(PathData)[1]
        OutputParams = TransfertSimu(PathSimu)

        # Chi2 reduction
        # NbParamsLib = 0
        # for k in range(2,11):
        #     for j in range(OutputParams[0]):
        #         if np.min(OutputParams[k][j]) != np.max(OutputParams[k][j]):
        #             NbParamsLib += 1 # variable parameter
        #             print(NbParamsLib)

        # NbData = 2*InputData[0]
        # for k in range(InputData[0]):
        #     if InputData[-2][k] == 1:
        #         NbData -= 1 # Input data with corr egal 1 count for 1 and not 2 (x, y)

        # reduct = NbData - NbParamsLib
        # OutputParams[-2] = OutputParams[-2]/reduct
    
        # Best orbits
        BestOrbits = BestOrbitsClass(*OutputParams, NbPtsEllipse, SystDist, InputData[0], InputData[-2])
        BestOrbitsWidget = BestOrbits.Widget
        BestOrbitsParams = BestOrbits.BestParams
        BestOrbitsEllipses = BestOrbits.BestEllipses
        Layout.addWidget(BestOrbitsWidget)

        # Select orbits
        SelectOrbits = SelectOrbitsClass(*OutputParams, NbSelectOrbits, NbPtsEllipse, SystDist)
        SelectOrbitsParams = SelectOrbits.SelectParams
        SelectOrbitsEllipses = SelectOrbits.SelectEllipses

        # Separation
        Layout.addWidget(Delimiter())

        # Space view
        SpaceViewWidget = SpaceView(InputData, SelectOrbitsEllipses, BestOrbitsEllipses)
        Layout.addWidget(SpaceViewWidget)

        # Temporal study
        TempoWidget = TempoView(InputData, SelectOrbitsEllipses, BestOrbitsEllipses)
        Layout.addWidget(TempoWidget)

        # Convergence of orbit parameters
        ConvWidget = Conv(OutputParams)
        Layout.addWidget(ConvWidget)

        # Histogram of orbit parameters
        HistWidget = Hist(OutputParams, BestOrbitsParams)
        Layout.addWidget(HistWidget)

        # Histogram 2D of orbit parameters
        Hist2DWidget = Hist2D(OutputParams, BestOrbitsParams)
        Layout.addWidget(Hist2DWidget)

        # Corner plot
        CornerWidget = Corner(SelectOrbitsParams)
        Layout.addWidget(CornerWidget)

        # Position at date
        PosAtDateWidget = PosAtDate(InputData, SelectOrbitsEllipses, BestOrbitsEllipses)
        Layout.addWidget(PosAtDateWidget)

        # Widget container
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))
        
        # Showing
        self.show()
        # LoadWin.close()

    # Close programme when the main window are closed
    def closeEvent(self, e):
        try:
            app.closeAllWindows()
        except:
            self.SignalCloseWindowMain.emit()




if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation

    # n = '1'
    Simu = '/Users/lacquema/Simulations/ggtau/ggtau_Ab2_Ab1_fit/ggtau_Ab2_Ab1_fit_1/'
    PathInputData = Simu+'ggtau_Ab12_hci.dat'
    PathOutputData = Simu+'solggtauAb12.dat'
    SystDist = 240  
    SystName = 'ggtau'

    # LoadWin = LoadWindowClass() # Loading window showing
    # app.processEvents() # Continue the program
    WindowMain = WindowMainClass(SystName, PathInputData, PathOutputData, SystDist, 10000, 500) # Main window showing
    sys.exit(app.exec()) # Application execution

    