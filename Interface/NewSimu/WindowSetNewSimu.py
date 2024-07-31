#! /var/guix/profiles/per-user/lacquema/Oracle/bin/python3
import sys
import os
import shutil
import subprocess
sys.path.append(os.path.dirname(__file__))


### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication

# My packages
from Tabs import *


### --- Parameters Window Generating --- ###

class WindowSetNewSimu(QMainWindow):

    SignalCloseWindowSetNewSimu = pyqtSignal() # initiation of the closeEvent signal
    
    def __init__(self):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Window characteristics
        self.setWindowTitle('Settings of the new simulation')

        # Widget Container
        self.Container = QTabWidget()

        # Tab 1
        self.Tab1 = TabSimuSet()
        self.Container.addTab(self.Tab1, 'Simulation settings')
        self.Tab1.SimuPath.EditPath.textChanged.connect(self.ChangeStartOrder)
        self.Tab1.SimuName.EditParam.textChanged.connect(self.ChangeStartOrder)

        # Tab 2
        self.Tab2 = TabDataSet()
        self.Container.addTab(self.Tab2, 'Data settings')
        self.Tab2.RelRV.CheckData.stateChanged.connect(self.EnableOrNotPriorJitter)
        self.Tab2.AbsRV.CheckData.stateChanged.connect(self.EnableOrNotPriorJitter)

        # Tab 3
        self.Tab3 = TabPriorSet()
        self.Container.addTab(self.Tab3, 'Priors')
        self.Tab3.BtnReset.clicked.connect(self.EnableOrNotPriorJitter)

        # Tab 4
        self.Tab4 = TabStartSet()
        self.Tab4.BtnStart.clicked.connect(self.StartSimulation)
        self.Container.addTab(self.Tab4, 'Starting')
        self.Tab4.NbHours.SpinParam.valueChanged.connect(self.ChangeStartOrder)

        # Container
        self.setCentralWidget(self.Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))

    def EnableOrNotPriorJitter(self):
        if self.Tab2.RelRV.CheckData.isChecked() or self.Tab2.AbsRV.CheckData.isChecked():
            self.Tab3.CheckJitter.setEnabled(True)
        else:
            self.Tab3.CheckJitter.setEnabled(False)
            self.Tab3.CheckJitter.setChecked(False)


    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowSetNewSimu.emit() 

    def ChangeStartOrder(self):
        self.Tab4.StartOrder.EditParam.setText(f'oarsub -l nodes=1/core=8,walltime={self.Tab4.NbHours.SpinParam.value()} --project dynapla {self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text()}/{self.Tab1.InputFileName.EditParam.text()}')
    
    def StartSimulation(self):
        if len(self.Tab1.SimuPath.EditPath.text())==0:
            print('Simulation path not given.')
            print('Check your inputs.')
            print('\n')
        else:
            if os.path.exists(self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text()):
                print('This directory already exists.')
                print('Check your inputs.')
                print('\n')
            else: 
                os.makedirs(self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text())
                if len(self.Tab2.PathData.EditPath.text())==0:
                    os.rmdir(self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text())
                    print('Data file not given.')
                    print('Check your inputs.')
                    print('\n')
                else:
                    if (self.Tab2.RelAstro.CheckData.isChecked() and self.Tab2.RelAstro.FormatData.currentIndex()==0) \
                        or (self.Tab2.AbsAstro.CheckData.isChecked() and self.Tab2.AbsAstro.FormatData.currentIndex()==0) \
                        or (self.Tab2.RelRV.CheckData.isChecked() and self.Tab2.RelRV.FormatData.currentIndex()==0) \
                        or (self.Tab2.AbsRV.CheckData.isChecked() and self.Tab2.AbsRV.FormatData.currentIndex()==0):
                        os.rmdir(self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text())
                        print('Data format not given.')
                        print('Check your inputs.')
                        print('\n')
                    else:
                        print(f'{self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text()}/ directory was created.')
                        shutil.copy(self.Tab2.PathData.EditPath.text(), self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text()+'/'+self.Tab2.DataFileName.EditParam.text())                    
                        print('Data file was copied.')
                        self.DoInputShell()
                        print('Input shell file was created.')
                        if self.Tab4.CheckOrder.CheckParam.isChecked():
                            command = 'cd '+self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text()+';chmod u+x '+self.Tab1.InputFileName.EditParam.text()+';'+self.Tab4.StartOrder.EditParam.text()
                            print(command)
                            result = subprocess.run(command, shell=True, text=True)
                            print('Simulation launched')
                            print('\n')
                            error = result.stderr
                            if type(error)!= type(None):
                                print(result.stderr)
                                print('Simulation not launched but you can still launch yourself the input shell file created in the desired directory.\n')
                                print('\n')
                        else:
                            print('All you have to do is launch the input shell file created in the desired directory.')
                            print('\n')

                    
    


    def DoInputShell(self):
        with open(self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text()+'/'+self.Tab1.InputFileName.EditParam.text(), "w") as file:
            file.write('#! /bin/bash\nexport OMP_NUM_THREADS=8\nexport STACKSIZE=1000000\n'+self.DirPath+'/../../Algorithm/bin/astrom_mcmcop <<!') # Header
            file.write('\n')
            if self.Tab4.CheckLM.isChecked() and self.Tab4.CheckMCMC.isChecked(): # Choice
                file.write('2')
                file.write(' ! New simulation LM and MCMC')
            elif not self.Tab4.CheckLM.isChecked() and self.Tab4.CheckMCMC.isChecked():
                file.write('3')
                file.write(' ! New simulation LM only')
            elif self.Tab4.CheckLM.isChecked() and not self.Tab4.CheckMCMC.isChecked():
                file.write('4')
                file.write(' ! New simulation MCMC only')
            file.write('\n')
            if self.Tab2.RelAstro.CheckData.isChecked(): file.write('1 ') # Data
            else: file.write('0 ')
            if self.Tab2.AbsRV.CheckData.isChecked(): file.write('1 ')
            else: file.write('0 ')
            if self.Tab2.RelRV.CheckData.isChecked(): file.write('1 ')
            else: file.write('0 ')
            if self.Tab2.AbsAstro.CheckData.isChecked(): file.write('1')
            else: file.write('0')
            file.write(' ! Type of data')
            file.write('\n')
            file.write(str(self.Tab3.NbOrbitsValue)) # Number of orbits
            file.write(' ! Number of orbits')
            file.write('\n')
            if self.Tab2.AbsAstro.CheckData.isChecked() or self.Tab2.RelAstro.CheckData.isChecked() or self.Tab2.AbsRV.CheckData.isChecked() or self.Tab2.RelRV.CheckData.isChecked(): 
                file.write(self.Tab2.DataFileName.EditParam.text())
                file.write(' ! Data file')
                file.write('\n')
            file.write('1d-'+str(self.Tab1.Precision.SpinParam.value())) # Precision of simulation
            file.write(' ! Precision')
            file.write('\n')
            if self.Tab2.RelAstro.CheckData.isChecked(): # Format of data
                file.write(str(self.Tab2.RelAstro.FormatData.currentIndex()))
                file.write('\n')
            if self.Tab2.AbsRV.CheckData.isChecked(): 
                file.write(str(self.Tab2.AbsRV.FormatData.currentIndex()))
                file.write('\n')
            if self.Tab2.RelRV.CheckData.isChecked(): 
                file.write(str(self.Tab2.RelRV.FormatData.currentIndex()))
                file.write('\n')
            if self.Tab2.AbsAstro.CheckData.isChecked(): 
                file.write(str(self.Tab2.AbsAstro.FormatData.currentIndex()))
                file.write('\n')
            if self.Tab3.CheckJitter.isChecked(): 
                file.write('1')
                file.write(' ! Jitter')
            else: 
                file.write('0')
                file.write(' ! No Jitter')
            file.write('\n')
            file.write(self.Tab3.SystDist.SpinParam.text()+' '+self.Tab3.SystDistUnit.ComboParam.currentText())
            file.write(' ! Distance')
            file.write('\n')
            file.write(str(float(self.Tab3.TablePriors.item(0,0).text())*0.000954588))
            file.write(' ! First guess of center mass (ms)')
            file.write('\n')
            file.write(self.Tab1.OutFileName.EditParam.text())
            file.write(' ! Result file')
            file.write('\n')
            file.write(self.Tab1.DumpFileName.EditParam.text())
            file.write(' ! Dump file')
            file.write('\n')
            file.write(self.Tab1.DumpFreq.SpinParam.text())
            file.write(' ! Dump frequency')
            file.write('\n')
            c=0 # Number of mass prior
            for i in range(len(self.Tab3.ListPriorMass)):
                DistribIndex = self.Tab3.ListPriorMass[i].Layout.itemAt(3*self.Tab3.NbBodies.SpinParam.value()+1).widget().ComboParam.currentIndex()
                if DistribIndex != 0: c+=1
            file.write(str(c))
            file.write(' ! Number of masses prior')
            file.write('\n')
            for i in range(len(self.Tab3.ListPriorMass)):
                DistribIndex = self.Tab3.ListPriorMass[i].Layout.itemAt(3*self.Tab3.NbBodies.SpinParam.value()+1).widget().ComboParam.currentIndex()
                if DistribIndex != 0:
                    for j in range(self.Tab3.NbBodies.SpinParam.value()):
                        file.write(self.Tab3.ListPriorMass[i].Layout.itemAt(3*j+1).widget().SpinParam.text()+' ')
                    file.write(' ! Coefficients')
                    file.write('\n')
                    file.write(str(DistribIndex))
                    file.write(' ! Distribution (1:Normal, 2:Log, 3:Uniform, 4:Fixed)')
                    file.write('\n')
                    if DistribIndex == 1 or DistribIndex == 2:
                        file.write(self.Tab3.ListPriorMass[i].Mean.SpinParam.text())
                        file.write(' ')
                        file.write(self.Tab3.ListPriorMass[i].SD.SpinParam.text())
                        file.write(' ')
                    if DistribIndex == 3:
                        file.write(self.Tab3.ListPriorMass[i].Min.SpinParam.text())
                        file.write(' ')
                        file.write(self.Tab3.ListPriorMass[i].Max.SpinParam.text())
                        file.write(' ')
                    if DistribIndex == 4:
                        file.write(self.Tab3.ListPriorMass[i].Value.SpinParam.text())
                        file.write(' ')
                    file.write(self.Tab3.ListPriorMass[i].PriorUnit.ComboParam.currentText())
                    file.write(' ! Distribution parameters')
                    file.write('\n')
            file.write(self.Tab3.RefTime.SpinParam.text())
            file.write(' ! Reference of time')
            file.write('\n')
            if self.Tab3.CheckJitter.isChecked(): 
                file.write(self.Tab3.Jitter.SpinParam.text()+' '+self.Tab3.V0.SpinParam.text())
                file.write(' ! Initial VO and Jitter')
                file.write('\n')
            for i in range(1, self.Tab3.NbBodies.SpinParam.value()):
                for j in range(len(self.Tab3.LabelParams)):
                    file.write(self.Tab3.TablePriors.item(i, j).text())
                    if j == 0: file.write(' mj')
                    file.write(' ')
                file.write(' ! First guess of orbit parameters (m[mj] a[AU] e i[deg] Om[deg] om[deg] tp[MJD])')
                file.write('\n')
            file.write(self.Tab3.PMin.SpinParam.text()+' '+self.Tab3.PMax.SpinParam.text())
            file.write(' ! Range of permited period')
            file.write('\n')
            file.write(self.Tab3.aMin.SpinParam.text()+' '+self.Tab3.aMax.SpinParam.text())
            file.write(' ! Range of permited half major axis')
            file.write('\n')
            file.write(self.Tab3.eMin.SpinParam.text()+' '+self.Tab3.eMax.SpinParam.text())
            file.write(' ! Range of permited eccentricity')
            file.write('\n')
            file.write('exit')
            file.write('\n')
            file.write('!')
            
            



            


            





# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetNewSimu()
    WindowParam.show()
    app.exec() # Application execution
