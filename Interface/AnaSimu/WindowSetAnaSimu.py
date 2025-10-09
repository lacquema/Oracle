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
from TransferData import HeaderDataIn, TransfertSimu


### --- Parameters Window Generating --- ###

class WindowSetAnaSimu(WindowWithFinder):

    SignalCloseWindowSetAnaSimu = pyqtSignal()
    ReSignalCloseWindowMain = pyqtSignal()
    
    def __init__(self):
        super().__init__()

        # File to save the last path used
        self.last_path_file = os.path.join(os.path.dirname(__file__), '.last_path')

        # Window characteristics
        self.setWindowTitle('Settings of the analysis')
        self.setMinimumWidth(1000)
        self.setMinimumHeight(500)

        # Layout initialisation
        self.Layout = QVBoxLayout()

        # Reset button
        self.BtnReset = QPushButton('Reset')
        self.BtnReset.setStatusTip('Reset all tab settings')
        self.Layout.addWidget(self.BtnReset)
        self.BtnReset.clicked.connect(self.ResetParams)
        
        self.InitWidgets()

        # Pré-remplir le champ avec le dernier chemin utilisé si dispo
        last_path = self.load_last_path()
        if last_path:
            # print(f'\nLast path used: {last_path}')
            self.SimuFilePathW.EditParam.setText(last_path)
            # self.check_change_path(file_path=last_path)

        # Widget Container
        self.Container = QWidget()
        self.Container.setLayout(self.Layout)

        # Add container to the split main window
        self.Splitter.addWidget(self.Container)

        # Connect folder to edit path
        self.Finder.doubleClicked.connect(self.check_change_path)

        # Status bar
        self.setStatusBar(QStatusBar(self))

    
    def check_change_path(self):
        index = self.Finder.selectedIndexes()[0]
        info = self.Model.fileInfo(index)
        if info.isDir():
            print('\nSelected file is a directory')
            self.ClearEdits()
        else:
            file_path = info.absoluteFilePath()
            if self.is_valid_file(file_path):  # Check if the file is valid
                self.SimuFilePathW.EditParam.setText(file_path)
            else:
                print('\nSelected file is not valid')
                self.ClearEdits()

    def save_last_path(self, path):
        try:
            with open(self.last_path_file, "w", encoding="utf-8") as f:
                f.write(path)
        except Exception as e:
            print(f"Erreur lors de la sauvegarde du chemin : {e}")

    def load_last_path(self):
        try:
            if os.path.exists(self.last_path_file):
                with open(self.last_path_file, "r", encoding="utf-8") as f:
                    return f.read().strip()
        except Exception as e:
            print(f"Erreur lors du chargement du chemin : {e}")
        return ""

    def is_valid_file(self, file_path):
        # Add logic to validate the file (e.g., check extension, content, etc.)
        return file_path.endswith('.dat')  # Example: only allow .dat files

    def InitWidgets(self):
        
        # self.SimuPath = LineEdit('Directory path', 'Path to the adjustment to analyse', '')
        # self.Layout.addWidget(self.SimuPath)
        # self.SimuPath.EditParam.textChanged.connect(self.FindSettings)

        self.SimuFilePathW = LineEdit('Simulation file', 'Path to the simulation file to analyse', '')
        self.Layout.addWidget(self.SimuFilePathW)

        # self.CheckHeader = CheckBox('Header with data', 'Check if the simulation file has a header with data')
        # self.Layout.addWidget(self.CheckHeader, alignment=Qt.AlignmentFlag.AlignLeft)
        # self.CheckHeader.CheckParam.setChecked(True)
        
        # self.DataFileNameW = LineEdit('Data file', 'Name of data file with extension', '')
        # self.Layout.addWidget(self.DataFileNameW, alignment=Qt.AlignmentFlag.AlignLeft)

        self.SystDistW = DoubleSpinBox('System distance', 'Distance from us of the studied system', 0, 0, None)
        self.Layout.addWidget(self.SystDistW, alignment=Qt.AlignmentFlag.AlignLeft)
        self.SystDistW.setVisible(False)

        self.SystDistUnitW = ComboBox(None, 'Unit', ['pc', 'mas'])
        self.SystDistW.Layout.addWidget(self.SystDistUnitW, alignment=Qt.AlignmentFlag.AlignLeft)

        self.UnivYNW = CheckBox('Universal variable', 'Check if the simulation is in the universal variable')
        self.Layout.addWidget(self.UnivYNW)
        self.UnivYNW.setVisible(False)

        self.SimuFilePathW.EditParam.textChanged.connect(self.HeaderOptionsVisible)



        # self.HeaderCheck(self.CheckHeader.CheckParam.isChecked())

        # self.BtnAutoComp = QPushButton('Autocomplete by go file')
        # self.Layout.addWidget(self.BtnAutoComp)
        # self.BtnAutoComp.clicked.connect(self.DialBrowseInputFile)

        self.Layout.addWidget(Delimiter(Title='Options :'))

        self.NbSelectOrbits = SpinBox('Number of orbits', 'Number of ramdom selected orbits to analyse', 10000, 1, None, 1)
        self.Layout.addWidget(self.NbSelectOrbits, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.NbPtsEllipse = SpinBox('Points by ellipse', 'Number of points by computed ellipse', 500, 1, None, 1)
        # self.Layout.addWidget(self.NbPtsEllipse, alignment=Qt.AlignmentFlag.AlignLeft)

        self.BtnStart = QPushButton('Analyse the simulation')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnStart.clicked.connect(self.AnalyseSimu)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # self.EnableBtnStartOrNot()
        # self.SimuFileNameW.EditParam.textChanged.connect(self.EnableBtnStartOrNot)
        # self.DataFileNameW.EditParam.textChanged.connect(self.EnableBtnStartOrNot)
        # self.SystDistW.SpinParam.valueChanged.connect(self.EnableBtnStartOrNot)

    def HeaderOptionsVisible(self):
        try:
            if HeaderDataIn(self.SimuFilePathW.EditParam.text()):
                self.SystDistW.setVisible(False)
                self.UnivYNW.setVisible(False)
            else:
                self.SystDistW.setVisible(True)
                self.UnivYNW.setVisible(True)
        except:
            self.SystDistW.setVisible(False)
            self.UnivYNW.setVisible(False)

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


    # def FindSettings(self):
    #     directory = self.SimuPath.EditParam.text()
    #     largest_file = None
    #     largest_size = 0

    #     # Find the largest .dat file in the directory
    #     for file_name in os.listdir(directory):
    #         if file_name.endswith('.dat'):
    #             file_path = os.path.join(directory, file_name)
    #             file_size = os.path.getsize(file_path)
    #             if file_size > largest_size:
    #                 largest_size = file_size
    #                 largest_file = file_name

    #     if largest_file:
    #         self.DataFileName = largest_file
    #         self.SimuFilePathW.EditParam.setText(self.DataFileName)
    #     else:
    #         print('\nNo .dat files found in the directory')
    #         self.ClearEdits()


    def ClearEdits(self):
        # self.DataFileNameW.EditParam.setText('')
        self.SimuFilePathW.EditParam.setText('')
        # self.SystDistW.SpinParam.setValue(0)
        # self.SystDistUnitW.ComboParam.setCurrentIndex(0)
        
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()


    def AnalyseSimu(self):
        if len(self.SimuFilePathW.EditParam.text()) == 0:
            print('\nSimulation path not given')
        else:
            try:
                self.OpenWinMain()
                self.save_last_path(self.SimuFilePathW.EditParam.text())
            except Exception as e:
                print(f'\nThere is a problem with inputs: {e}')


    def OpenWinMain(self):
        if self.UnivYNW.CheckParam.isChecked():
            self.UnivYN = 2
        else:
            self.UnivYN = 1
        PathSimu = self.SimuFilePathW.EditParam.text()
        SimuName = PathSimu.split('/')[-1]
        InputData, OutputParams = TransfertSimu(PathSimu, self.UnivYN)
        if InputData != None:
            self.SystDistValue = InputData['SystDist']['pc']
        else:
            if self.SystDistUnitW.ComboParam.currentText() == 'pc':
                self.SystDistValue = self.SystDistW.SpinParam.value()
            elif self.SystDistUnitW.ComboParam.currentText() == 'mas':
                self.SystDistValue = 1000/self.SystDistW.SpinParam.value() # mas to pc

        self.WinMain = WindowMainClass(SimuName, InputData, OutputParams,
                                       self.NbSelectOrbits.SpinParam.value(), 
                                       self.SystDistValue)
        
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
