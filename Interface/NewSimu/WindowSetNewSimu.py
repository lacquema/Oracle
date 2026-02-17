#! /Users/lacquema/Oracle.env/bin/python3
import sys
import os
sys.path.append(os.path.dirname(__file__)+'/..')
import shutil
import subprocess


### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication

# My packages
from Tabs import *
from WindowWithFinder import WindowWithFinder


### --- Parameters Window Generating --- ###

class WindowSetNewSimu(WindowWithFinder):

    SignalCloseWindowSetNewSimu = pyqtSignal() # initiation of the closeEvent signal
    
    def __init__(self):
        super().__init__()


        self.minimumSize = (1200, 500)
        self.resize(*self.minimumSize)

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Window characteristics
        self.setWindowTitle('Settings of the new simulation')

        # Widget Container
        self.Container = QTabWidget()

        # Tab 1
        self.TabSimuSet = TabSimuSet()
        self.Container.addTab(self.TabSimuSet, 'Simulation settings')
        self.InitInterTabConnect('TabSimuSet')

        # Tab 2
        self.TabDataSet = TabDataSet()
        self.Container.addTab(self.TabDataSet, 'Data settings')
        self.InitInterTabConnect('TabDataSet')

        # Tab 3
        self.TabGuessSet = TabGuessSet()
        self.Container.addTab(self.TabGuessSet, 'Guesses')
        self.InitInterTabConnect('TabGuessSet')

        # Tab 4
        self.TabRangeSet = TabRangeSet()
        self.Container.addTab(self.TabRangeSet, 'Ranges')
        self.InitInterTabConnect('TabRangeSet')

        # Tab 5
        self.TabStartSet = TabStartSet()
        # self.TabStartSet.BtnStart.clicked.connect(self.StartSimulation)
        self.Container.addTab(self.TabStartSet, 'Starting')
        self.InitInterTabConnect('TabStartSet')

        # Initial state of some widgets
        self.EnableUnivVarOrNot(self.TabSimuSet.CheckUnivVar.CheckParam.isChecked())
        self.ChangeNbBodies()

        # Container
        # self.setCentralWidget(self.Container)

        # Add container to the split main window
        self.Splitter.addWidget(self.Container)
        # self.Container.setMinimumWidth(1000)

        # Connect folder to edit path
        self.Finder.doubleClicked.connect(self.ChangePath)

        # Status bar
        self.setStatusBar(QStatusBar(self))


    def ChangePath(self):
        index = self.Finder.selectedIndexes()[0]
        info = self.Model.fileInfo(index)
        self.TabSimuSet.SimuPath.EditParam.setText(info.absoluteFilePath()+'/')


    def InitInterTabConnect(self, IdTab):
        if IdTab=='TabSimuSet':
            self.TabSimuSet.BtnReset.clicked.connect(lambda: self.ReInit('TabSimuSet'))
            self.TabSimuSet.CheckLM.stateChanged.connect(self.StartBtnAvailableOrNot)
            self.TabSimuSet.CheckMCMC.stateChanged.connect(self.StartBtnAvailableOrNot)
            self.TabSimuSet.CheckUnivVar.CheckParam.stateChanged.connect(self.EnableUnivVarOrNot)
        elif IdTab=='TabDataSet':
            self.TabDataSet.BtnReset.clicked.connect(lambda: self.ReInit('TabDataSet'))
            self.TabDataSet.NbBodies.SpinParam.valueChanged.connect(self.ChangeNbBodies)
        elif IdTab=='TabGuessSet':
            self.TabGuessSet.BtnReset.clicked.connect(lambda: self.ReInit('TabGuessSet'))
        elif IdTab=='TabRangeSet':
            self.TabRangeSet.BtnReset.clicked.connect(lambda: self.ReInit('TabRangeSet'))
        elif IdTab=='TabStartSet':
            self.TabStartSet.BtnReset.clicked.connect(lambda: self.ReInit('TabStartSet'))
            self.TabStartSet.BtnCreate.clicked.connect(self.CreateInputFiles)
            self.TabStartSet.BtnStart.clicked.connect(self.Start)
        
    def ReInit(self, IdTab):
        self.InitInterTabConnect(IdTab)
        self.EnableUnivVarOrNot(self.TabSimuSet.CheckUnivVar.CheckParam.isChecked())
        self.ChangeNbBodies()
        self.StartBtnAvailableOrNot()
        if IdTab=='TabGuessSet':
            for i in range(self.TabGuessSet.TablePriors.rowCount()):
                for j in range(self.TabGuessSet.TablePriors.columnCount()):
                    self.TabGuessSet.TablePriors.setItem(i, j, QTableWidgetItem('0.'))
                    self.TabGuessSet.TablePriors.item(i, j).setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        if IdTab=='TabDataSet':
            for x in self.TabRangeSet.ListPriorMassId:
                self.TabRangeSet.DelThisPriorMass(x)


    def EnableUnivVarOrNot(self, state):
        self.TabRangeSet.PMin.setEnabled(not state)
        self.TabRangeSet.aMin.setEnabled(not state)
        self.TabRangeSet.PeriMin.setEnabled(state)
        if not state:
            self.TabRangeSet.eMin.SpinParam.setMaximum(1)
            self.TabRangeSet.eMax.SpinParam.setMaximum(1)
            self.TabGuessSet.LabelParams[1] = 'a [AU]'
        else:
            self.TabRangeSet.eMin.SpinParam.setMaximum(10)
            self.TabRangeSet.eMax.SpinParam.setMaximum(10)
            self.TabGuessSet.LabelParams[1] = 'q [AU]'

        self.TabGuessSet.TablePriors.setColumnCount(len(self.TabGuessSet.LabelParams))
        self.TabGuessSet.TablePriors.setHorizontalHeaderLabels(self.TabGuessSet.LabelParams)

    def InitNbBodies(self, NewValue):
        self.TabDataSet.NbBodiesValue = NewValue
        self.TabRangeSet.NbBodiesValue = NewValue
        self.TabGuessSet.NbBodiesValue = NewValue

    def ChangeNbBodies(self):
        # Update
        # print("Old value:", self.TabDataSet.NbBodiesValue
        #       , "\nNew value:", self.TabDataSet.NbBodies.SpinParam.value())
        OldValue = self.TabDataSet.NbBodiesValue
        NewValue = self.TabDataSet.NbBodies.SpinParam.value()
        # Change the number of bodies in the different tabs
        self.InitNbBodies(NewValue)
        # Change the text of the number of orbits
        NbOrbitsValue = NewValue-1
        if  NbOrbitsValue==1: self.TabDataSet.NbOrbits.setText(f'\t=>    {NbOrbitsValue} Orbit')
        else: self.TabDataSet.NbOrbits.setText(f'\t=>    {NbOrbitsValue} Orbits')
        # Change the number of row in the priors table

        if self.TabGuessSet.TablePriors.rowCount() < NbOrbitsValue:
            print("rowcount inf a new value")
            self.TabGuessSet.TablePriors.setRowCount(NbOrbitsValue)
            for n in range(len(self.TabGuessSet.LabelParams)):
                self.TabGuessSet.TablePriors.setItem(NbOrbitsValue-1, n, QTableWidgetItem('0.'))
                self.TabGuessSet.TablePriors.item(NbOrbitsValue-1, n).setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        else:
            self.TabGuessSet.TablePriors.setRowCount(NbOrbitsValue)

        if NewValue > OldValue:
            for x in self.TabRangeSet.ListPriorMass:
                x.Layout.insertWidget(x.Layout.indexOf(x.Distrib)-1, QLabel('+'), alignment=Qt.AlignmentFlag.AlignLeft)
                x.Layout.insertWidget(x.Layout.indexOf(x.Distrib)-1, SpinBox(None, 'Coefficient of mass', 0, 0, 1, 1), alignment=Qt.AlignmentFlag.AlignLeft)
                x.Layout.insertWidget(x.Layout.indexOf(x.Distrib)-1, QLabel('m'+str(NewValue)), alignment=Qt.AlignmentFlag.AlignLeft)

        elif NewValue < OldValue:
            for x in self.TabRangeSet.ListPriorMass:
                for i in range(3):
                    x.Layout.removeWidget(x.Layout.itemAt(x.Layout.indexOf(x.Distrib)-2).widget())

    # def EnableOrNotPriorJitter(self):
    #     if self.TabDataSet.RelRV.CheckParam.isChecked() or self.TabDataSet.AbsRV.CheckParam.isChecked():
    #         self.TabGuessSet.CheckJitter.setEnabled(True)
    #     else:
    #         self.TabGuessSet.CheckJitter.setEnabled(False)
    #         self.TabGuessSet.CheckJitter.CheckParam.setChecked(False)


    def StartBtnAvailableOrNot(self):
        if self.TabSimuSet.CheckLM.isChecked() or self.TabSimuSet.CheckMCMC.isChecked():
            self.TabStartSet.BtnStart.setEnabled(True)
            self.TabStartSet.BtnCreate.setEnabled(True)
        else:
            self.TabStartSet.BtnStart.setEnabled(False)
            self.TabStartSet.BtnCreate.setEnabled(False)

    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowSetNewSimu.emit() 


    def Start(self):
        self.GoPath = self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()+'/start.sh'
        if not os.path.exists(self.GoPath):
            self.CreateInputFiles()
            if os.path.exists(self.GoPath):
                print(f'{self.TabStartSet.StartOrder.EditParam.text()} {self.GoPath} &')
                # subprocess.run(f'cd {self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()}', shell=True, text=True)
                subprocess.run(f'chmod +x {self.GoPath}', shell=True, text=True)
                subprocess.run(f'{self.TabStartSet.StartOrder.EditParam.text()} {self.GoPath} &', shell=True, text=True, cwd=self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text())
        else:
            print(f'{self.TabStartSet.StartOrder.EditParam.text()} {self.GoPath} &')
            # subprocess.run(f'cd {self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()}', shell=True, text=True)
            subprocess.run(f'chmod +x {self.GoPath}', shell=True, text=True)
            subprocess.run(f'{self.TabStartSet.StartOrder.EditParam.text()} {self.GoPath} &', shell=True, text=True, cwd=self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text())

    # def ChangeStartOrder(self):
    #     self.TabStartSet.StartOrder.EditParam.setText(f'oarsub -l nodes=1/core={self.TabStartSet.NbCores.SpinParam.value()},walltime={self.TabStartSet.NbHours.SpinParam.value()} --project dynapla {self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()}/{self.TabSimuSet.InputFileName.EditParam.text()}')
    
    def CreateInputFiles(self):
        if len(self.TabSimuSet.SimuPath.EditParam.text())==0:
            print('\nSimulation path not given')
        else:
            if os.path.exists(self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()):
                print('\nThis simulation already exists')
            else: 
                os.makedirs(self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text())
                if len(self.TabDataSet.PathData.EditPath.text())==0:
                    os.rmdir(self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text())
                    print('\nData file not given')
                else:
                    if self.TabDataSet.FormatRelAstro.ComboParam.currentIndex() == 0:
                        os.rmdir(self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text())
                        print('\nAstrometric data format not given')
                    else:
                        print(f'{self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()}/ directory was created.')
                        # subprocess.run(f'cd {self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()}', shell=True, text=True)
                        shutil.copy(self.TabDataSet.PathData.EditPath.text(), self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()+'/'+self.TabDataSet.PathData.EditPath.text().split('/')[-1])                    
                        print('\nData file was copied')
                        self.DoGoFile()
                        print('\nGo file was created')
                        print('Just run it')
                

    def DoGoFile(self):
        self.NbBodiesValue = self.TabDataSet.NbBodies.SpinParam.value()
        self.AlgoFileName = 'astrom'
        if self.TabSimuSet.CheckUnivVar.CheckParam.isChecked(): self.AlgoFileName += '_univ'
        self.AlgoFileName += '_mcmco'
        if self.TabStartSet.CheckParallel.CheckParam.isChecked(): self.AlgoFileName += '_par'

        self.EnvPath = '/'.join(self.DirPath.split('/')[:-2])

        self.SimuDir = self.TabSimuSet.SimuPath.EditParam.text()+self.TabSimuSet.SimuName.EditParam.text()+'/'
        
        with open(self.SimuDir+'start.sh', 'w') as file:
            file.write('#! /bin/bash')
            file.write('\n')
            file.write(f'cd {self.SimuDir}')
            file.write('\n')
            if self.TabStartSet.CheckParallel.CheckParam.isChecked():
                file.write('export OMP_NUM_THREADS='+self.TabStartSet.NbCores.SpinParam.text()) # Header
            else:
                file.write('export OMP_NUM_THREADS=1')
            file.write('\n')
            file.write('export STACKSIZE=1000000')
            file.write('\n')
            file.write(self.EnvPath+'/Code/bin/'+self.AlgoFileName+' <<!')
            file.write('\n')
            if self.TabSimuSet.CheckLM.isChecked() and self.TabSimuSet.CheckMCMC.isChecked(): # Choice
                file.write('2')
                file.write(' # New simulation LM and MCMC')
            elif not self.TabSimuSet.CheckLM.isChecked() and self.TabSimuSet.CheckMCMC.isChecked():
                file.write('3')
                file.write(' # New simulation LM only')
            elif self.TabSimuSet.CheckLM.isChecked() and not self.TabSimuSet.CheckMCMC.isChecked():
                file.write('4')
                file.write(' # New simulation MCMC only')
            file.write('\n')
            if self.TabDataSet.CheckRelAstro.CheckParam.isChecked(): file.write('1 ') # Data
            else: file.write('0 ')
            if self.TabDataSet.AbsRV.CheckParam.isChecked(): file.write('1 ')
            else: file.write('0 ')
            if self.TabDataSet.RelRV.CheckParam.isChecked(): file.write('1 ')
            else: file.write('0 ')
            if self.TabDataSet.CheckAbsAstro.CheckParam.isChecked(): file.write('1')
            else: file.write('0')
            file.write(' # Type of data (RelAstro? AbsRV? RelRV? AbsAstro?)')
            file.write('\n')
            file.write(str(self.NbBodiesValue-1)) # Number of orbits
            file.write(' # Number of orbits')
            file.write('\n')
            if self.TabDataSet.CheckAbsAstro.CheckParam.isChecked() or self.TabDataSet.CheckRelAstro.CheckParam.isChecked() or self.TabDataSet.AbsRV.CheckParam.isChecked() or self.TabDataSet.RelRV.CheckParam.isChecked(): 
                file.write(self.SimuDir+self.TabDataSet.PathData.EditPath.text().split('/')[-1])
                # file.write(' # Data file')
                file.write('\n')
            file.write('1d-'+str(self.TabSimuSet.Precision.SpinParam.value())) # Precision of simulation
            file.write(' # Precision')
            file.write('\n')
            file.write(str(self.TabDataSet.FormatDate.ComboParam.currentIndex()+1))
            file.write(' '+str(self.TabDataSet.FormatRelAstro.ComboParam.currentIndex()))
            if self.TabDataSet.CheckCorrCoef.CheckParam.isChecked(): file.write(' 1')
            else: file.write(' 0')
            file.write(' # Format data (1=DDMMYYYY/2=JD 1=(DEC,RA)/2=(SEP,PA) CorrCoeff?')
            file.write('\n')
            if self.TabDataSet.AbsRV.CheckParam.isChecked():
                if self.TabDataSet.CheckJitter.CheckParam.isChecked(): 
                    file.write('1')
                else: 
                    file.write('0')
                file.write(' # Jitter?')
                file.write('\n')
            file.write(self.TabDataSet.SystDist.SpinParam.text()+' '+self.TabDataSet.SystDistUnit.ComboParam.currentText())
            file.write(' # Distance')
            file.write('\n')
            file.write(self.TabGuessSet.FirstGuessCenterMass.SpinParam.text())
            file.write(' # First guess of center mass [Msun]')
            file.write('\n')
            # file.write(self.SimuDir+self.TabSimuSet.OutFileName.EditParam.text()+'.dat')
            file.write('results.dat')
            # file.write(' # Result file')
            file.write('\n')
            # file.write(self.SimuDir+self.TabSimuSet.DumpFileName.EditParam.text()+'.dat')
            file.write('dump.dat')
            # file.write(' # Dump file')
            file.write('\n')
            file.write(self.TabSimuSet.DumpFreq.SpinParam.text())
            file.write(' # Dump frequency')
            file.write('\n')
            c=0 # Number of mass prior
            for i in range(len(self.TabRangeSet.ListPriorMass)):
                DistribIndex = self.TabRangeSet.ListPriorMass[i].Layout.itemAt(3*self.NbBodiesValue+1).widget().ComboParam.currentIndex()
                if DistribIndex != 0: c+=1
            file.write(str(c))
            file.write(' # Number of masses prior')
            file.write('\n')
            for i in range(len(self.TabRangeSet.ListPriorMass)):
                DistribIndex = self.TabRangeSet.ListPriorMass[i].Layout.itemAt(3*self.NbBodiesValue+1).widget().ComboParam.currentIndex()
                if DistribIndex != 0:
                    for j in range(self.NbBodiesValue):
                        file.write(self.TabRangeSet.ListPriorMass[i].Layout.itemAt(3*j+1).widget().SpinParam.text()+' ')
                    file.write(' # Coefficients')
                    file.write('\n')
                    file.write(str(DistribIndex))
                    file.write(' # Distribution (1=Normal, 2=Log, 3=Uniform, 4=Fixed)')
                    file.write('\n')
                    if DistribIndex == 1 or DistribIndex == 2:
                        file.write(self.TabRangeSet.ListPriorMass[i].Mean.SpinParam.text())
                        file.write(' ')
                        file.write(self.TabRangeSet.ListPriorMass[i].SD.SpinParam.text())
                        file.write(' ')
                    if DistribIndex == 3:
                        file.write(self.TabRangeSet.ListPriorMass[i].Min.SpinParam.text())
                        file.write(' ')
                        file.write(self.TabRangeSet.ListPriorMass[i].Max.SpinParam.text())
                        file.write(' ')
                    if DistribIndex == 4:
                        file.write(self.TabRangeSet.ListPriorMass[i].Value.SpinParam.text())
                        file.write(' ')
                    if self.TabRangeSet.ListPriorMass[i].PriorUnit.ComboParam.currentText() == 'Msun':
                        file.write('ms')
                    elif self.TabRangeSet.ListPriorMass[i].PriorUnit.ComboParam.currentText() == 'Mjup':
                        file.write('mj')
                    file.write(self.TabRangeSet.ListPriorMass[i].PriorUnit.ComboParam.currentText())
                    file.write(' # Distribution parameters')
                    file.write('\n')
            file.write(self.TabGuessSet.RefTime.MJDWidget.SpinParam.text())
            file.write(' # Reference of time')
            file.write('\n')
            if self.TabDataSet.CheckJitter.CheckParam.isChecked(): 
                file.write(self.TabDataSet.Jitter.SpinParam.text()+' '+self.TabDataSet.V0.SpinParam.text())
                file.write(' # Jitter and V0 [km/s]')
                file.write('\n')
            for i in range(0, self.NbBodiesValue-1):
                for j in range(len(self.TabGuessSet.LabelParams)):
                    if float(self.TabGuessSet.TablePriors.item(i, j).text()) == 0:
                        file.write('0.0001')
                    else:
                        file.write(self.TabGuessSet.TablePriors.item(i, j).text())
                    if j == 0: 
                        if self.TabGuessSet.FirstGuessOtherMassUnit.ComboParam.currentIndex() == 0:
                            file.write(' mj')
                        else:
                            file.write(' ms')
                    file.write(' ')
                if not self.TabSimuSet.CheckUnivVar.CheckParam.isChecked(): 
                    file.write(f' # First guess of orbit {i} parameters (m[Mjup] a[AU] e i[deg] Om[deg] om[deg] tp[MJD])')                
                else: 
                    file.write(f' # First guess of orbit {i} parameters (m[Mjup] q[AU] e i[deg] Om[deg] om[deg] tp[MJD])')
                file.write('\n')
            if not self.TabSimuSet.CheckUnivVar.CheckParam.isChecked():
                file.write(self.TabRangeSet.PMin.SpinParam.text()+' '+self.TabRangeSet.PMax.SpinParam.text())
                file.write(' # Range of permited period')
                file.write('\n')
                file.write(self.TabRangeSet.aMin.SpinParam.text()+' '+self.TabRangeSet.aMax.SpinParam.text())
                file.write(' # Range of permited semi-major axis')
                file.write('\n')
            else:
                file.write(self.TabRangeSet.PeriMin.SpinParam.text()+' '+self.TabRangeSet.PeriMax.SpinParam.text())
                file.write(' # Range of permited periastron')
                file.write('\n')
            file.write(self.TabRangeSet.eMin.SpinParam.text()+' '+self.TabRangeSet.eMax.SpinParam.text())
            file.write(' # Range of permited eccentricity')
            file.write('\n')
            file.write('!')

            



            


            





# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetNewSimu()
    # WindowParam.show()
    WindowParam.show()
    app.exec() # Application execution
