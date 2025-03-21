#! /Users/lacquema/Oracle.env/bin/python3
import sys
import os
sys.path.append(os.path.dirname(__file__)+'/..')
import shutil
import subprocess

### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication, QVBoxLayout, QPushButton, QFileDialog, QTreeWidgetItem, QTreeView
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QIcon

# My packages
from Parameters import *
from Utils import *
from WindowWithFinder import WindowWithFinder


### --- Parameters Window Generating --- ###

class WindowSetContSimu(WindowWithFinder):

    SignalCloseWindowSetContSimu = pyqtSignal()
    
    def __init__(self):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Window characteristics
        self.setWindowTitle('Parameters to continue the simulation')
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
        
        self.SimuPath = LineEdit('Directory path', 'Path to the directory of the simulation to continue', '')
        self.Layout.addWidget(self.SimuPath)
        self.SimuPath.EditParam.textChanged.connect(self.FindSettings)
        # self.SimuPath.EditPath.textChanged.connect(self.ChangeStartOrder)

        # self.StartFileName = LineEdit('Startup file', 'Name you want to give to the input continuation shell file with extension', 'inputcont.sh')
        # self.Layout.addWidget(self.InputFileName, alignment=Qt.AlignmentFlag.AlignLeft)
        # self.InputFileName.EditParam.textChanged.connect(self.ChangeStartOrder)

        # self.DumpFileName = LineEdit('Dump file', 'Name of dump file with extension', '')
        # self.Layout.addWidget(self.DumpFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Layout.addWidget(Delimiter(Title='Options :'))

        # self.DumpFileName.Layout.addSpacing(100)
        self.DumpFreq = SpinBox('Save frequency', 'Save frequency in the dump file [yr]', 10000000, 0, 1000000000, 100000)
        self.Layout.addWidget(self.DumpFreq, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.BtnAutoComp = QPushButton('Autocomplete by go file')
        # self.Layout.addWidget(self.BtnAutoComp)
        # self.BtnAutoComp.clicked.connect(self.DialBrowseInputFile)

        # self.CheckParallel = CheckBox('Parallelization :', 'Parallelization of the simulation algorithm')
        # self.Layout.addWidget(self.CheckParallel)
        # self.CheckParallel.CheckParam.setChecked(True)

        # self.CheckParallel.Layout.setSpacing(60)
        # self.NbCores = SpinBox('Number of cores', 'Number of cores use for parallelization', 8, 1, None, 1)
        # self.CheckParallel.Layout.addWidget(self.NbCores, alignment=Qt.AlignmentFlag.AlignLeft)
        # self.CheckParallel.CheckParam.stateChanged.connect(self.NbCores.setEnabled)

        self.BtnCont = QPushButton('Create continuation file')
        self.Layout.addWidget(self.BtnCont, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnCont.clicked.connect(self.CreateContFile)

        self.Layout.addWidget(Delimiter(Title='Starting :'))

        # self.CheckOrder = CheckBox('Starting order :', 'If you want to start with bash order')
        # self.Layout.addWidget(self.CheckOrder)
        # self.CheckOrder.CheckParam.stateChanged.connect(self.CheckStartOrderChange)

        self.NbOrdersValue = 0
        self.OrdersValue = []
        with open(self.DirPath+'/../Orders.txt', 'r') as file:
            for x in file:
                self.NbOrdersValue += 1
                self.OrdersValue.append(x.replace('\n',''))

        self.ComboOrder = ComboBox('Saved orders', 'All orders saved', self.OrdersValue)
        self.Layout.addWidget(self.ComboOrder)

        self.BtnDelOrder = QPushButton('del')
        self.ComboOrder.Layout.addWidget(self.BtnDelOrder)
        self.BtnDelOrder.clicked.connect(self.DelOrder)

        self.BtnChangeOrder = QPushButton()
        self.BtnChangeOrder.setIcon(QIcon(f'{self.DirPath}/Items/arrowDown.png'))
        self.BtnChangeOrder.clicked.connect(self.ChangeOrder)
        self.ComboOrder.Layout.addWidget(self.BtnChangeOrder)

        self.StartOrder = LineEdit('Order', 'Terminal order to continue the adjustment', self.ComboOrder.ComboParam.currentText())
        self.Layout.addWidget(self.StartOrder)

        self.BtnSaveOrder = QPushButton('save')
        self.StartOrder.Layout.addWidget(self.BtnSaveOrder)
        self.BtnSaveOrder.clicked.connect(self.SaveOrder)

        self.LblGo = QLabel('./continuation.sh')
        self.StartOrder.Layout.addWidget(self.LblGo)

        self.BtnStart = QPushButton('Continue the simulation')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnStart.clicked.connect(self.Start)

        # self.CheckStartOrderChange(self.CheckOrder.CheckParam.isChecked())

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        self.BtnStart.setEnabled(False)
        self.BtnCont.setEnabled(False)

    # def CheckStartOrderChange(self, state):
    #     # self.NbHours.setEnabled(state)
    #     self.ComboOrder.setEnabled(state)
    #     self.StartOrder.setEnabled(state)
    #     self.BtnStart.setEnabled(state)

    def ChangeOrder(self):
        self.StartOrder.EditParam.setText(self.ComboOrder.ComboParam.currentText())

    def SaveOrder(self):
        with open(self.DirPath+'/../Orders.txt', 'a') as file:
            file.write('\n')
            file.write(self.StartOrder.EditParam.text())
            self.ComboOrder.ComboParam.addItem(self.StartOrder.EditParam.text())
            self.NbOrdersValue += 1
            self.ComboOrder.ComboParam.setCurrentIndex(self.NbOrdersValue-1)

    def DelOrder(self):
        index = self.ComboOrder.ComboParam.currentIndex()
        if index>2:
            lines = []
            c = 0 
            with open(self.DirPath+'/../Orders.txt', 'r') as file :
                for x in file:
                    if c != index:
                        lines.append(x.replace('\n',''))
                    c += 1
            with open(self.DirPath+'/../Orders.txt', 'w') as file :
                file.write('\n'.join(lines))
            self.ComboOrder.ComboParam.removeItem(index)
            self.NbOrdersValue -= 1
            self.ComboOrder.ComboParam.setCurrentIndex(-1)
        else:
            print('\nImpossible to remove this order')

    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()
        

    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowSetContSimu.emit() 

    # def ChangeStartOrder(self):
    #     self.StartOrder.EditParam.setText(f'oarsub -l nodes=1/core=8,walltime={self.NbHours.SpinParam.value()} --project dynapla {self.SimuPath.EditParam.text()+self.InputFileName.EditParam.text()}')

    def FindSettings(self):
        try:  
            with open(self.SimuPath.EditParam.text()+'start.sh', 'r') as file:
                self.SimuName = self.SimuPath.EditParam.text().split('/')[-2]
                GoFileLines = file.readlines()
                self.NbCores = GoFileLines[2][:-1].split('=')[-1]
                self.AlgoFileName = GoFileLines[4].split()[0].split('/')[-1]
                self.SimuFileName = GoFileLines[13][:-1]+'.dat'
                self.DumpFileName = GoFileLines[14][:-1].split('/')[-1]
                if self.AlgoFileName[-4:]=='_par':
                    print(f'\nDo you want continue the {self.SimuPath.EditParam.text()} simulation from {self.DumpFileName} dump file (in parallel with {self.NbCores} cores)?')
                else:
                    print(f'\nDo you want continue the {self.SimuPath.EditParam.text()} simulation from {self.DumpFileName} dump file (not in parallel)?')
                self.BtnStart.setEnabled(True)
                self.BtnCont.setEnabled(True)

        except:
            print('\nSimulation not found')
            self.BtnStart.setEnabled(False)
            self.BtnCont.setEnabled(False)


    def CreateContFile(self):
        self.ContFilePath = self.SimuPath.EditParam.text()+'continuation.sh'

        if len(self.SimuPath.EditParam.text())==0:
            print('\nSimulation path not given')
        else:
            self.DoContFile()
            print('\nContinuation file was created')
            print('Just run it')



    def DoContFile(self):
        self.EnvPath = '/'.join(self.DirPath.split('/')[:-2])

        with open(self.ContFilePath, "w") as file:
            file.write('#! /bin/bash') # Header
            file.write('\n')
            file.write(f'cd {self.SimuPath.EditParam.text()}')
            file.write('\n')
            # if self.CheckParallel.CheckParam.isChecked():
            # if self.AlgoFileName[-4:]=='_par':
            file.write('export OMP_NUM_THREADS='+self.NbCores) # Header
            # else:
            #     file.write('export OMP_NUM_THREADS=1')
            file.write('\n')
            file.write('export STACKSIZE=1000000')
            file.write('\n')
            file.write(self.EnvPath+'/Code/bin/'+self.AlgoFileName+' <<!')
            file.write('\n')
            file.write('1') # continuation
            file.write(' # Simulation continuation')
            file.write('\n')
            file.write(self.SimuPath.EditParam.text()+self.DumpFileName)
            file.write('\n')
            file.write(self.DumpFreq.SpinParam.text())
            file.write(' # Dump frequency')
            file.write('\n')
            file.write('exit')
            file.write('\n')
            file.write('!')


    def Start(self):
        self.GoPath = self.SimuPath.EditParam.text()+'continuation.sh'
        if not os.path.exists(self.GoPath):
            self.CreateContFile()
            if os.path.exists(self.GoPath):
                print(f'> {self.StartOrder.EditParam.text()} {self.GoPath} &')
                # subprocess.Popen(f'cd {self.SimuPath.EditParam.text()}', shell=True, text=True)
                # os.chdir(self.SimuPath.EditParam.text(), )
                subprocess.run(f'chmod +x {self.GoPath}', shell=True, text=True)
                subprocess.run(f'{self.StartOrder.EditParam.text()} {self.GoPath} &', shell=True, text=True, cwd=self.SimuPath.EditParam.text())
        else:
            print(f'> {self.StartOrder.EditParam.text()} {self.GoPath} &')
            # subprocess.Popen(f'cd {self.SimuPath.EditParam.text()}', shell=True, text=True)
            # os.chdir(self.SimuPath.EditParam.text())
            subprocess.run(f'chmod +x {self.GoPath}', shell=True, text=True)
            subprocess.run(f'{self.StartOrder.EditParam.text()} {self.GoPath} &', shell=True, text=True, cwd=self.SimuPath.EditParam.text())





# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetContSimu()
    WindowParam.show()
    app.exec() # Application execution

