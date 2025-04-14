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
from WindowWithFinder import WindowWithFinder


### --- Parameters Window Generating --- ###

class WindowSetAnaSimu(WindowWithFinder):

    SignalCloseWindowSetAnaSimu = pyqtSignal()
    ReSignalCloseWindowMain = pyqtSignal()
    
    def __init__(self):
        super().__init__()

        # Window characteristics
        self.setWindowTitle('Settings of the analysis')
        self.setMinimumWidth(1000)

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

        # # Container
        # self.setCentralWidget(self.Container)

        # Add container to the split main window
        self.Splitter.addWidget(self.Container)

        # Connect folder to edit path
        self.Finder.doubleClicked.connect(self.ChangePath)

        # Status bar
        self.setStatusBar(QStatusBar(self))

    
    def ChangePath(self):
        index = self.Finder.selectedIndexes()[0]
        info = self.Model.fileInfo(index)
        self.SimuPath.EditParam.setText(info.absoluteFilePath()+'/')




    def InitWidgets(self):
        
        self.SimuPath = LineEdit('Directory path', 'Path to the adjustment to analyse', '')
        self.Layout.addWidget(self.SimuPath)
        self.SimuPath.EditParam.textChanged.connect(self.FindSettings)

        # self.InputFileNameW = LineEdit('Start file', 'Name of the input simulation file with extension', 'start.sh')
        # self.Layout.addWidget(self.InputFileNameW, alignment=Qt.AlignmentFlag.AlignLeft)
        # self.InputFileNameW.EditParam.textChanged.connect(self.FindSettings)

        self.SimuFileNameW = LineEdit('Simulation file', 'Name of the simulation file to analyse with extension', '')
        self.Layout.addWidget(self.SimuFileNameW, alignment=Qt.AlignmentFlag.AlignLeft)

        self.CheckHeader = CheckBox('Header with data', 'Check if the simulation file has a header with data')
        self.Layout.addWidget(self.CheckHeader, alignment=Qt.AlignmentFlag.AlignLeft)
        self.CheckHeader.CheckParam.setChecked(True)
        
        # self.DataFileNameW = LineEdit('Data file', 'Name of data file with extension', '')
        # self.Layout.addWidget(self.DataFileNameW, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.SystDistW = DoubleSpinBox('System distance', 'Distance from us of the studied system', 0, 0, None)
        # self.Layout.addWidget(self.SystDistW, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.SystDistUnitW = ComboBox(None, 'Unit', ['pc', 'mas'])
        # self.SystDistW.Layout.addWidget(self.SystDistUnitW, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.HeaderCheck(self.CheckHeader.CheckParam.isChecked())

        # self.BtnAutoComp = QPushButton('Autocomplete by go file')
        # self.Layout.addWidget(self.BtnAutoComp)
        # self.BtnAutoComp.clicked.connect(self.DialBrowseInputFile)

        self.Layout.addWidget(Delimiter(Title='Options :'))

        self.NbSelectOrbits = SpinBox('Number of orbits', 'Number of ramdom selected orbits to analyse', 10000, 1, None, 1)
        self.Layout.addWidget(self.NbSelectOrbits, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.NbPtsEllipse = SpinBox('Points by ellipse', 'Number of points by computed ellipse', 500, 1, None, 1)
        # self.Layout.addWidget(self.NbPtsEllipse, alignment=Qt.AlignmentFlag.AlignLeft)

        self.BtnStart = QPushButton('Analyse')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnStart.clicked.connect(self.AnalyseSimu)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # self.EnableBtnStartOrNot()
        # self.SimuFileNameW.EditParam.textChanged.connect(self.EnableBtnStartOrNot)
        # self.DataFileNameW.EditParam.textChanged.connect(self.EnableBtnStartOrNot)
        # self.SystDistW.SpinParam.valueChanged.connect(self.EnableBtnStartOrNot)

    # def HeaderCheck(self, state):
    #     self.Header = state
    #     self.DataFileNameW.setVisible(not state)
    #     self.SystDistW.setVisible(not state)
    #     self.SystDistUnitW.setVisible(not state)


    # def EnableBtnStartOrNot(self):
    #     self.BtnStart.setEnabled(False)
    #     if self.SimuFileNameW.EditParam.text()!=0 and self.DataFileNameW.EditParam.text()!=0 and self.SystDistW.SpinParam.value()!=0:
    #         # if self.DataFileName[-4:]=='.dat' and self.SimuFileName[-4:]=='.dat' and type(self.SystDist)==float:
    #         self.BtnStart.setEnabled(True)


    def FindSettings(self):
        try:  
            with open(self.SimuPath.EditParam.text()+'start.sh', 'r') as file:
                self.SimuName = self.SimuPath.EditParam.text().split('/')[-2]
                GoFileLines = file.readlines()
                self.DataFileName = GoFileLines[8].split('/')[-1][:-1]
                self.SimuFileName = GoFileLines[13].split('/')[-1][:-1]
                self.SystDist = float(GoFileLines[11].split(' ')[0])
                self.SystDistUnit = GoFileLines[11].split(' ')[1]
                if self.DataFileName[-4:]=='.dat' and self.SimuFileName[-4:]=='.dat' and type(self.SystDist)==float:
                    # self.DataFileNameW.EditParam.setText(self.DataFileName)
                    self.SimuFileNameW.EditParam.setText(self.SimuFileName)
                    # self.SystDistW.SpinParam.setValue(self.SystDist)
                    # self.SystDistUnitW.ComboParam.setCurrentText(self.SystDistUnit)
                else:
                    print('\nAutocomplete failed')
                    self.ClearEdits()
                    # self.EnableBtnStartOrNot()
                # print(f'Do you want analyse {self.SimuPath.EditParam.text()} simulation with {self.SimuFileName} solution file, {self.DataFileName} data file and a system distance of {self.SystDist} {self.SystDistUnit} ?')
        except:
            print('\nAutocomplete failed')
            # self.EnableBtnStartOrNot()
            self.ClearEdits()


    def ClearEdits(self):
        # self.DataFileNameW.EditParam.setText('')
        self.SimuFileNameW.EditParam.setText('')
        # self.SystDistW.SpinParam.setValue(0)
        # self.SystDistUnitW.ComboParam.setCurrentIndex(0)
        
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()


    def AnalyseSimu(self):
        if len(self.SimuPath.EditParam.text()) == 0:
            print('\nSimulation directory path not given')
        else:
            try:
                self.OpenWinMain()
            except:
                print('\nThere is a problem with inputs')


    def OpenWinMain(self):
        # if self.Header:
        #     DataFile = None
            # self.SystDistValue = None
        # else:
        #     DataFile = self.SimuPath.EditParam.text()+self.DataFileNameW.EditParam.text()
            # if self.SystDistUnitW.ComboParam.currentText() == 'pc':
            #     self.SystDistValue = self.SystDistW.SpinParam.value()
            # elif self.SystDistUnitW.ComboParam.currentText() == 'mas':
            #     self.SystDistValue = 1000/self.SystDistW.SpinParam.value() # mas to pc
    
        self.WinMain = WindowMainClass(self.SimuPath.EditParam.text().split('/')[-2], 
                                       self.SimuPath.EditParam.text()+self.SimuFileNameW.EditParam.text(),
                                       self.NbSelectOrbits.SpinParam.value(), 
                                       1000,
                                       self.CheckHeader.CheckParam.isChecked())
        
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
