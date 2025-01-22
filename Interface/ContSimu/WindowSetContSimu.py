#! /Users/lacquema/Oracle.env/bin/python3
import sys
import os
sys.path.append(os.path.dirname(__file__)+'/..')
import shutil
import subprocess

### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication, QVBoxLayout, QPushButton, QFileDialog
from PyQt6.QtCore import Qt, pyqtSignal

# My packages
from Parameters import *
from Utils import *


### --- Parameters Window Generating --- ###

class WindowSetContSimu(QMainWindow):

    SignalCloseWindowSetContSimu = pyqtSignal()
    
    def __init__(self):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Window characteristics
        self.setWindowTitle('Parameters to continue the simulation')
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
        
        self.SimuPath = PathBrowser('Directory path', 'Path to the directory of the simulation to continue', 0)
        self.Layout.addWidget(self.SimuPath)
        self.SimuPath.EditPath.textChanged.connect(self.FindSettings)
        # self.SimuPath.EditPath.textChanged.connect(self.ChangeStartOrder)

        # self.InputFileName = LineEdit('Continuation file', 'Name you want to give to the input continuation shell file with extension', 'inputcont.sh')
        # self.Layout.addWidget(self.InputFileName, alignment=Qt.AlignmentFlag.AlignLeft)
        # self.InputFileName.EditParam.textChanged.connect(self.ChangeStartOrder)

        # self.DumpFileName = LineEdit('Dump file', 'Name of dump file with extension', 'dump.dat')
        # self.Layout.addWidget(self.DumpFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Layout.addWidget(Delimiter(Title='Options:'))

        # self.DumpFileName.Layout.addSpacing(100)
        self.DumpFreq = SpinBox('Save frequency', 'Save frequency in the dump file [yr]', 10000000, 0, 1000000000, 100000)
        self.Layout.addWidget(self.DumpFreq, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.BtnAutoComp = QPushButton('Autocomplete by go file')
        # self.Layout.addWidget(self.BtnAutoComp)
        # self.BtnAutoComp.clicked.connect(self.DialBrowseInputFile)

        self.NbCores = SpinBox('Number of cores', 'Number of cores to be used', 8, 1, None, 1)
        self.Layout.addWidget(self.NbCores, alignment=Qt.AlignmentFlag.AlignLeft)

        self.CheckParallel = CheckBox('Parallelization', 'Parallelization of the simulation algorithm')
        self.Layout.addWidget(self.CheckParallel)

        # self.NbHours = SpinBox('Simulation duration', 'Simulation duration [hour]', 48, 1, 48, 1)
        # self.Layout.addWidget(self.NbHours, alignment=Qt.AlignmentFlag.AlignLeft)
        # self.NbHours.SpinParam.valueChanged.connect(self.ChangeStartOrder)

        # self.CheckOrder = CheckBox('Starting order :', 'If you just want to create the input file, but dont want to run the command in the terminal')
        # self.Layout.addWidget(self.CheckOrder)
        # self.CheckOrder.CheckParam.stateChanged.connect(lambda: self.StartOrder.setEnabled(self.CheckOrder.CheckParam.isChecked()))

        # self.StartOrderValue = f'oarsub -l nodes=1/core=8,walltime={self.NbHours.SpinParam.value()} --project dynapla ./inputcont.sh'
        # self.StartOrder = LineEdit(None, 'Terminal order to continue the adjustment', self.StartOrderValue)
        # self.CheckOrder.Layout.addWidget(self.StartOrder)
        # self.StartOrder.setEnabled(self.CheckOrder.CheckParam.isChecked())

        self.BtnCont = QPushButton('Continue')
        self.Layout.addWidget(self.BtnCont, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnCont.clicked.connect(self.ContinueSimu)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()
        

    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowSetContSimu.emit() 

    # def ChangeStartOrder(self):
    #     self.StartOrder.EditParam.setText(f'oarsub -l nodes=1/core=8,walltime={self.NbHours.SpinParam.value()} --project dynapla {self.SimuPath.EditPath.text()+self.InputFileName.EditParam.text()}')

    def FindSettings(self):
        try:  
            with open(self.SimuPath.EditPath.text()+'go_mcmco.sh', 'r') as file:
                self.SimuName = self.SimuPath.EditPath.text().split('/')[-2]
                GoFileLines = file.readlines()
                self.AlgoFileName = GoFileLines[3].split()[0].split('/')[-1][:-4]
                self.SimuFileName = GoFileLines[12][:-1]+'.dat'
                self.DumpFileName = GoFileLines[13][:-1]
                print(f'Do you want continue {self.SimuName} simulation from {self.DumpFileName} dump file ?')
        except:
            print('Simulation not found')

    def ContinueSimu(self):
        self.ContFilePath = self.SimuPath.EditPath.text()+'cont_mcmco.sh'
        if len(self.SimuPath.EditPath.text())==0:
            print('Simulation path not given')
        else:
            # if os.path.exists(self.ContFilePath):
            #     print('Continuation file already exists')
            # else:
                self.DoContFile()
                print('Continuation file was created')
                print('Just run it')

            # if self.CheckOrder.CheckParam.isChecked():
            #     result = subprocess.run(self.StartOrder.EditParam.text(), shell=True, capture_output=True, text=True)
            #     error = result.stderr
            #     if len(error)!=0:
            #         print(result.stderr)
            #         print('Simulation not launched but you can still launch yourself the input continuation shell file.')
            # else:
            #     print('All you have to do is launch the input continuation shell file.')



    def DoContFile(self):
        if self.CheckParallel.CheckParam.isChecked(): self.AlgoFileName += '_par'

        self.EnvPath = '/'.join(self.DirPath.split('/')[:-2])

        with open(self.ContFilePath, "w") as file:
            file.write('#! /bin/bash') # Header
            file.write('\n')
            file.write('export OMP_NUM_THREADS='+self.NbCores.SpinParam.text())
            file.write('\n')
            file.write('export STACKSIZE=1000000')  
            file.write('\n')
            file.write(self.EnvPath+'/Code/bin/'+self.AlgoFileName+' <<!')
            file.write('\n')
            file.write('1') # continuation
            file.write(' # Simulation continuation')
            file.write('\n')
            file.write(self.SimuPath.EditPath.text()+self.DumpFileName)
            file.write('\n')
            file.write(self.DumpFreq.SpinParam.text())
            file.write(' # Dump frequency')
            file.write('\n')
            file.write('!')



# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetContSimu()
    WindowParam.show()
    app.exec() # Application execution

