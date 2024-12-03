#! /Users/lacquema/Oracle.env/bin/python3
import sys
import os
sys.path.append(os.path.dirname(__file__)+'/..')

### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication, QVBoxLayout, QPushButton, QFileDialog
from PyQt6.QtCore import Qt, pyqtSignal

# My packages
# from NewSimu.TabsParamNew import *
from Parameters import *
from Utils import DelAllWidgetsBtw
from WindowMain import WindowMainClass
# from WindowMenu import LoadWindowClass


### --- Parameters Window Generating --- ###

class WindowSetAnaSimu(QMainWindow):

    SignalCloseWindowSetAnaSimu = pyqtSignal()
    ReSignalCloseWindowMain = pyqtSignal()
    
    def __init__(self):
        super().__init__()

        # Window characteristics
        self.setWindowTitle('Settings of the analysis')
        self.setMinimumWidth(600)

        # Layout initialisation
        self.Layout = QVBoxLayout()

        # Reset button
        self.BtnReset = QPushButton('Reset')
        self.BtnReset.setStatusTip('Reset all tab settings')
        self.Layout.addWidget(self.BtnReset)
        self.BtnReset.clicked.connect(self.ResetParams)
        
        self.InitWidgets()

        # Widget Container
        self.Container = QWidget()
        self.Container.setLayout(self.Layout)

        # Container
        self.setCentralWidget(self.Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))



    def InitWidgets(self):
        
        self.SimuPath = PathBrowser('Directory path', 'Path to the adjustment to analyse', 0)
        self.Layout.addWidget(self.SimuPath)
        self.SimuPath.EditPath.textChanged.connect(self.FindSettings)

        # self.InputFileName = LineEdit('Input file', 'Name of the input simulation file with extension', 'start.sh')
        # self.Layout.addWidget(self.InputFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.SimuFileName = LineEdit('Simulation file', 'Name of the simulation file to analyse with extension', '')
        # self.Layout.addWidget(self.SimuFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.DataFileName = LineEdit('Data file', 'Name of data file with extension', '')
        # self.Layout.addWidget(self.DataFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.SystDistValue = 0
        # self.SystDist = DoubleSpinBox('System distance', 'Distance from us of the studied system', 0, 0, None, 0.01)
        # self.Layout.addWidget(self.SystDist, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.SystDistUnit = ComboBox(None, 'Unit', ['pc', 'mas'])
        # self.SystDist.Layout.addWidget(self.SystDistUnit, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.BtnAutoComp = QPushButton('Autocomplete by go file')
        # self.Layout.addWidget(self.BtnAutoComp)
        # self.BtnAutoComp.clicked.connect(self.DialBrowseInputFile)

        self.Layout.addWidget(Delimiter(Title='Options:'))

        self.NbSelectOrbits = SpinBox('Number of orbits', 'Number of ramdom selected orbits to analyse', 10000, 1, None, 1)
        self.Layout.addWidget(self.NbSelectOrbits, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.NbPtsEllipse = SpinBox('Points by ellipse', 'Number of points by computed ellipse', 500, 1, None, 1)
        # self.Layout.addWidget(self.NbPtsEllipse, alignment=Qt.AlignmentFlag.AlignLeft)

        self.BtnStart = QPushButton('Analyse')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnStart.clicked.connect(self.AnalyseSimu)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)


    def FindSettings(self):
        try:  
            with open(self.SimuPath.EditPath.text()+'go_mcmco.sh', 'r') as file:
                self.SimuName = self.SimuPath.EditPath.text().split('/')[-2]
                GoFileLines = file.readlines()
                self.DataFileName = GoFileLines[7][:-1]
                self.SimuFileName = GoFileLines[12][:-1]+'.dat'
                self.SystDist = float(GoFileLines[10].split(' ')[0])
                self.SystDistUnit = GoFileLines[10].split(' ')[1]
                print(f'Do you want analyse {self.SimuName} simulation with {self.SimuFileName} solution file, {self.DataFileName} data file and a system distance of {self.SystDist} {self.SystDistUnit} ?')
        except:
            print('Simulation not found')
        
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()


    def AnalyseSimu(self):
        if len(self.SimuPath.EditPath.text()) == 0:
            print('Simulation directory path not given')
        else:
            self.OpenWinMain()


    def OpenWinMain(self):
        if self.SystDistUnit == 'pc':
            self.SystDistValue = self.SystDist
        elif self.SystDistUnit == 'mas':
            self.SystDistValue = 1000/self.SystDist
    
        self.WinMain = WindowMainClass(self.SimuPath.EditPath.text().split('/')[-2], self.SimuPath.EditPath.text()+self.DataFileName, self.SimuPath.EditPath.text()+self.SimuFileName, self.SystDistValue, self.NbSelectOrbits.SpinParam.value(), 1000)
        self.WinMain.SignalCloseWindowMain.connect(self.ReSignalCloseWindowMain.emit)
        self.WinMain.show()
        self.close()


    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        try: # In the case self.WinMain openning, dont show again WindowMenu
            if self.WinMain.isVisible() == False:
                self.SignalCloseWindowSetAnaSimu.emit() 
        except:
            self.SignalCloseWindowSetAnaSimu.emit() 
        

# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetAnaSimu()
    WindowParam.show()
    app.exec() # Application execution
