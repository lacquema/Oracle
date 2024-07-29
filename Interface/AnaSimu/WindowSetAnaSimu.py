#! /Users/lacquema/ByeGildas/bin/python3
import sys
import os

### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication, QVBoxLayout, QPushButton, QFileDialog
from PyQt6.QtCore import Qt, pyqtSignal

# My packages
# from NewSimu.TabsParamNew import *
from Parameters import *
from UtilsAnaSimu import DelAllWidgetsBtw
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
        # self.SimuPath.EditPath.textChanged.connect(self.FindInputFiles)

        self.SimuFileName = LineEdit('Adjustment file', 'Name of the adjustment solution file with extension', 'adjustment.dat')
        self.Layout.addWidget(self.SimuFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.DataFileName = LineEdit('Data file', 'Name of data file with extension', 'data.txt')
        self.Layout.addWidget(self.DataFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.SystDistValue = 0
        self.SystDist = DoubleSpinBox('System distance', 'Distance from us of the studied system', 0, 0, None, 0.01)
        self.Layout.addWidget(self.SystDist, alignment=Qt.AlignmentFlag.AlignLeft)

        self.SystDistUnit = ComboBox(None, 'Unit', ['pc', 'mas'])
        self.SystDist.Layout.addWidget(self.SystDistUnit, alignment=Qt.AlignmentFlag.AlignLeft)

        self.NbSelectOrbits = SpinBox('Selected orbits number', 'Number of ramdom selected orbits for this analysis', 10000, 1, None, 1)
        self.Layout.addWidget(self.NbSelectOrbits, alignment=Qt.AlignmentFlag.AlignLeft)

        self.NbPtsEllipse = SpinBox('Points by ellipse', 'Number of points by computed ellipse', 500, 1, None, 1)
        self.Layout.addWidget(self.NbPtsEllipse, alignment=Qt.AlignmentFlag.AlignLeft)

        self.BtnStart = QPushButton('Analyse the adjustment')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnStart.clicked.connect(self.AnalyseSimu)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()

    def AnalyseSimu(self):
        if len(self.SimuPath.EditPath.text()) == 0:
            print('Adjustement directory path not given.')
            print('Check your inputs.')
        else:
            if self.SystDist.SpinParam.value() == 0:
                print('System distance is zero.')
                print('Check your inputs.')
            else:
                self.OpenWinMain()

    def OpenWinMain(self):
        if self.SystDistUnit.ComboParam.currentText() == 'pc':
            self.SystDistValue = self.SystDist.SpinParam.value()
        elif self.SystDistUnit.ComboParam.currentText() == 'mas':
            self.SystDistValue = 1000/self.SystDist.SpinParam.value()
        try:
            self.WinMain = WindowMainClass(self.SimuPath.EditPath.text()[:-1].split('/')[-1], self.SimuPath.EditPath.text()+self.DataFileName.EditParam.text(), self.SimuPath.EditPath.text()+self.SimuFileName.EditParam.text(), self.SystDistValue, self.NbSelectOrbits.SpinParam.value(), self.NbPtsEllipse.SpinParam.value())
            self.WinMain.SignalCloseWindowMain.connect(self.ReSignalCloseWindowMain.emit)
            self.WinMain.show()
            self.close()
        except:
            print('Adjustment not found: check the directory path and the name of input files.')

    # def FindInputFiles(self):
    #     Files = os.listdir(self.SimuPath.EditPath.text()[:-1])
    #     self.SimuName.EditParam.setText(self.SimuPath.EditPath.text()[:-1].split('/')[-1])
    #     txtFiles = []
    #     datFiles = []
    #     for x in Files:
    #         if x[-4:] == '.txt':
    #             txtFiles.append(x)
    #         if x[-4:] == '.dat':
    #             datFiles.append(x)
    #     if 'logfile.dat' in datFiles: datFiles.remove('logfile.dat')
    #     if len(txtFiles) == 1:
    #         self.DataFileName.EditParam.setText(txtFiles[0])
    #     else:
    #         self.DataFileName.EditParam.setText('')
    #         print('Data file not found')
    #     if len(datFiles) == 1:
    #         self.SimuFileName.EditParam.setText(datFiles[0])
    #     else:
    #         self.SimuFileName.EditParam.setText('')
    #         print('Simulation file not found')
        

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
