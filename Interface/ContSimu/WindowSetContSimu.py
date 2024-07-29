#! /Users/lacquema/ByeGildas/bin/python3
import sys
import os
import shutil
import subprocess

### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication, QVBoxLayout, QPushButton, QFileDialog
from PyQt6.QtCore import Qt, pyqtSignal

# My packages
from Parameters import *
from UtilsContSimu import *


### --- Parameters Window Generating --- ###

class WindowSetContSimu(QMainWindow):

    SignalCloseWindowSetContSimu = pyqtSignal()
    
    def __init__(self):
        super().__init__()

        # Window characteristics
        self.setWindowTitle('Settings of the running ajustement')
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

        # Container
        self.setCentralWidget(self.Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))



    def InitWidgets(self):
        
        self.SimuPath = PathBrowser('Directory path', 'Path to the directory of the adjustment to analyse', 0)
        self.Layout.addWidget(self.SimuPath)
        self.SimuPath.EditPath.textChanged.connect(self.ChangeStartOrder)

        self.InputFileName = LineEdit('Continuation file', 'Name you want to give to the input continuation shell file with extension', 'inputcont.sh')
        self.Layout.addWidget(self.InputFileName, alignment=Qt.AlignmentFlag.AlignLeft)
        self.InputFileName.EditParam.textChanged.connect(self.ChangeStartOrder)

        self.DumpFileName = LineEdit('Dump file', 'Name of dump file with extension', 'dump.dat')
        self.Layout.addWidget(self.DumpFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.DumpFileName.Layout.addSpacing(100)
        self.DumpFreq = SpinBox('Save frequency', 'Save frequency in the dump file [yr]', 10000000, 0, 1000000000, 100000)
        self.DumpFileName.Layout.addWidget(self.DumpFreq, alignment=Qt.AlignmentFlag.AlignLeft)

        self.NbHours = SpinBox('Simulation duration', 'Simulation duration [hour]', 48, 1, 48, 1)
        self.Layout.addWidget(self.NbHours, alignment=Qt.AlignmentFlag.AlignLeft)
        self.NbHours.SpinParam.valueChanged.connect(self.ChangeStartOrder)

        self.CheckOrder = CheckBox('Starting order :', 'If you just want to create the input file, but dont want to run the command in the terminal')
        self.Layout.addWidget(self.CheckOrder)
        self.CheckOrder.CheckParam.stateChanged.connect(lambda: self.StartOrder.setEnabled(self.CheckOrder.CheckParam.isChecked()))

        self.StartOrderValue = f'oarsub -l nodes=1/core=8,walltime={self.NbHours.SpinParam.value()} --project dynapla ./inputcont.sh'
        self.StartOrder = LineEdit(None, 'Terminal order to continue the adjustment', self.StartOrderValue)
        self.CheckOrder.Layout.addWidget(self.StartOrder)
        self.StartOrder.setEnabled(self.CheckOrder.CheckParam.isChecked())

        self.BtnStart = QPushButton('Continue the ajustement')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnStart.clicked.connect(self.ContinueSimulation)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()
        

    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowSetContSimu.emit() 

    def ChangeStartOrder(self):
        self.StartOrder.EditParam.setText(f'oarsub -l nodes=1/core=8,walltime={self.NbHours.SpinParam.value()} --project dynapla {self.SimuPath.EditPath.text()+self.InputFileName.EditParam.text()}')

    def ContinueSimulation(self):
        if len(self.SimuPath.EditPath.text())==0:
            print('Simulation path not given.')
            print('Check your inputs.')
        else:
            self.DoInputContShell()
            print('Input continuation shell file was created.')
            if self.CheckOrder.CheckParam.isChecked():
                result = subprocess.run(self.StartOrder.EditParam.text(), shell=True, capture_output=True, text=True)
                error = result.stderr
                if len(error)!=0:
                    print(result.stderr)
                    print('Simulation not launched but you can still launch yourself the input continuation shell file.')
            else:
                print('All you have to do is launch the input continuation shell file.')


    def DoInputContShell(self):
        with open(self.SimuPath.EditPath.text()+"InputCont.sh", "w") as file:
            file.write('#! /bin/bash\nexport OMP_NUM_THREADS=8\nexport STACKSIZE=1000000\n./astrom_mcmcop <<!') # Header
            file.write('\n')
            file.write('1') # continuation
            file.write('\n')
            file.write(self.DumpFileName.EditParam.text())
            file.write('\n')
            file.write(self.DumpFreq.SpinParam.text())
            file.write('\n')
            file.write('exit')
            file.write('\n')
            file.write('!')





        


# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetContSimu()
    WindowParam.show()
    app.exec() # Application execution

