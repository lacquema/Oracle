#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import sys
import os

# PyQt packages
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtWidgets import QGridLayout, QTableWidgetItem, QTableWidget, QFileDialog, QTabWidget, QScrollArea, QMainWindow, QVBoxLayout, QHBoxLayout, QPushButton, QWidget, QStatusBar, QApplication, QProgressBar, QLabel, QCheckBox
from PyQt6.QtGui import QIcon, QFont

# My packages
from Data import DataClass
from Parameters import *
from PriorMass import PriorMassClass
from Utils import DelAllWidgetsBtw



class GeneralTab(QWidget):
    
    def __init__(self):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Layout initialisation
        self.Layout = QVBoxLayout()

        # Reset button
        self.BtnReset = QPushButton('Reset')
        self.BtnReset.setStatusTip('Reset tab settings')
        self.Layout.addWidget(self.BtnReset)
        self.BtnReset.clicked.connect(self.ResetParams)
        
        self.InitWidgets()

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        
        self.setLayout(self.Layout)

    # Widgets initialisations
    def InitWidgets(self):
        return
    
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()


    def ValidatedIfIn(self, WidgetEditing, ListValidCharacters):
        TextNew = WidgetEditing.text()
        c = 0
        for i in range(len(TextNew)):
            if TextNew[i] not in  ListValidCharacters:
                c=+1
                break
        if c != 0:
            WidgetEditing.setText(self.TextOld)
            WidgetEditing.setCursorPosition(i)
        else:
            self.TextOld = TextNew
    


class TabSimuSet(GeneralTab):
    
    def __init__(self):
        super().__init__()

    def InitWidgets(self):

        
        self.SimuPath = LineEdit('Path', 'Path where create the adjustment folder', '')
        self.Layout.addWidget(self.SimuPath)

        self.SimuPath.Layout.addSpacing(20)

        self.SimuName = LineEdit('New folder', 'Name you want to give to the adjustment folder', '')
        self.SimuPath.Layout.addWidget(self.SimuName)

        # self.InputFileName = LineEdit('Start file', 'Name you want to give to the start file', 'go_mcmco')
        # self.Layout.addWidget(self.InputFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Layout.addWidget(Delimiter(Title='Adjustment :'))

        self.WidgetH = QWidget()
        self.LayoutH = QHBoxLayout()
        self.CheckLM = QCheckBox('Levenberg-Marquardt')
        self.CheckLM.setChecked(True)
        self.LayoutH.addWidget(self.CheckLM, alignment=Qt.AlignmentFlag.AlignLeft)
        self.LayoutH.setSpacing(100)
        self.CheckMCMC = QCheckBox('MCMC')
        self.CheckMCMC.setChecked(True)
        self.LayoutH.addWidget(self.CheckMCMC, alignment=Qt.AlignmentFlag.AlignLeft)
        self.WidgetH.setLayout(self.LayoutH)
        self.Layout.addWidget(self.WidgetH, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Precision = SpinBox('Precision order', 'Order of adjustment precision in powers of 10', 7, 0, 10, 1)
        self.Layout.addWidget(self.Precision, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Layout.addWidget(Delimiter(Title='Outputs :'))

        # self.DumpFileName = LineEdit('Dump file', 'Name you want to give to the dump file', 'dump')
        # self.Layout.addWidget(self.DumpFileName, alignment=Qt.AlignmentFlag.AlignLeft)
        # self.DumpFileName.Layout.addSpacing(50)
        self.DumpFreq = SpinBox('Save frequency', 'Number of iterations between two saves to dump file', 10000000, 0, 1000000000, 100000)
        self.Layout.addWidget(self.DumpFreq, alignment=Qt.AlignmentFlag.AlignLeft)

        # self.OutFileName = LineEdit('Results file', 'Name you want to give to the results file', 'simulation')
        # self.Layout.addWidget(self.OutFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.Layout.setSpacing(0)


class TabDataSet(GeneralTab):
    
    def __init__(self):
        super().__init__()


    def InitWidgets(self):

        self.PathData = PathBrowser('Path to data file', 'Path to the existing data file', 1)
        self.Layout.addWidget(self.PathData)
        # self.PathData.setEnabled(self.CheckAbsAstro.CheckData.isChecked() or self.CheckRelAstro.CheckData.isChecked() or self.AbsRV.CheckData.isChecked() or self.RelRV.CheckData.isChecked())

        # self.DataFileName = LineEdit(' --> \t Data file', 'Name you want to give to the data file', 'data.txt')
        # self.PathData.Layout.addWidget(self.DataFileName)
        # self.DataFileName.setEnabled(self.CheckAbsAstro.CheckData.isChecked() or self.CheckRelAstro.CheckData.isChecked() or self.AbsRV.CheckData.isChecked() or self.RelRV.CheckData.isChecked())

        self.FormatDate = ComboBox('Dates format', 'Format of dates', ['Day Month Year', 'MJD'])
        self.Layout.addWidget(self.FormatDate)

        self.WidgetDataType = QWidget()
        self.LayoutDataTypeH = QHBoxLayout()
        
        self.LayoutAstroV = QVBoxLayout()

        self.LayoutAstroV.addWidget(Delimiter(Title='Astrometric data :'))

        self.CheckRelAstro = CheckBox('Relative:', 'If you want use relative astrometric data')
        self.LayoutAstroV.addWidget(self.CheckRelAstro)
        # self.CheckRelAstro.CheckData.setCheckState()
        self.CheckRelAstro.CheckParam.setChecked(True)
        self.CheckRelAstro.CheckParam.clicked.connect(lambda: self.CheckRelAstro.CheckParam.setChecked(True))
        self.CheckRelAstro.Layout.setSpacing(0)
        self.CheckRelAstro.Layout.addSpacing(30)
        self.FormatRelAstro = ComboBox(None, 'Format of astrometric data', ['Format', 'Dec RA', 'Sep PA'], 0)
        self.CheckRelAstro.Layout.addWidget(self.FormatRelAstro)
        # self.CheckRelAstro.Layout.addSpacing(10)
        self.CheckCorrCoef = CheckBox('Correlation', 'If you wish to use a correlation coefficient between coordinates')
        self.CheckRelAstro.Layout.addWidget(self.CheckCorrCoef)

        self.CheckAbsAstro = CheckBox('Absolute', 'If you want use absolute astrometric data')
        self.LayoutAstroV.addWidget(self.CheckAbsAstro)
        # self.CheckAbsAstro.CheckData.stateChanged.connect(self.EnableOrNotPutDataPath)

        self.LayoutAstroV.setSpacing(0)

        self.LayoutDataTypeH.addLayout(self.LayoutAstroV)

        self.LayoutRVV = QVBoxLayout()

        self.LayoutRVV.addWidget(Delimiter(Title='Radial velocity data :'))

        self.LayoutRVV.addSpacing(10)

        self.RelRV = CheckBox('Relative', 'If you want use relative radial velocity data')
        self.LayoutRVV.addWidget(self.RelRV)
        # self.RelRV.CheckData.stateChanged.connect(self.EnableOrNotPutDataPath)

        self.LayoutRVV.addSpacing(10)

        self.AbsRV = CheckBox('Absolute', 'If you want use absolute radial velocity data')
        self.LayoutRVV.addWidget(self.AbsRV)
        # self.AbsRV.CheckData.stateChanged.connect(self.EnableOrNotPutDataPath)

        self.LayoutDataTypeH.addLayout(self.LayoutRVV)

        self.WidgetDataType.setLayout(self.LayoutDataTypeH)
        self.Layout.addWidget(self.WidgetDataType)

        self.LayoutAstroV.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.LayoutRVV.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)

    # def EnableOrNotPutDataPath(self):
    #     self.DataFileName.setEnabled(self.CheckAbsAstro.CheckData.isChecked() or self.CheckRelAstro.CheckData.isChecked() or self.AbsRV.CheckData.isChecked() or self.RelRV.CheckData.isChecked())
    #     self.PathData.setEnabled(self.CheckAbsAstro.CheckData.isChecked() or self.CheckRelAstro.CheckData.isChecked() or self.AbsRV.CheckData.isChecked() or self.RelRV.CheckData.isChecked())
    #     if not self.CheckAbsAstro.CheckData.isChecked() and not self.CheckRelAstro.CheckData.isChecked() and not self.AbsRV.CheckData.isChecked() and not self.RelRV.CheckData.isChecked():
    #         print('Adjustments need data.')


class TabPriorSet(GeneralTab):
    
    def __init__(self):
        super().__init__()


    def InitWidgets(self):

        # self.BtnBetaPic = QPushButton('Beta Pic model')
        # self.Layout.addWidget(self.BtnBetaPic)

        # self.BtnGGTau = QPushButton('GGTau model')
        # self.Layout.addWidget(self.BtnGGTau)

        self.ListPriorMass = []
        self.ListPriorMassId = []
        self.c = 0 # counter

        self.WidgetH = QWidget()
        self.LayoutH = QHBoxLayout()
        
        self.LayoutV1 = QVBoxLayout()

        self.NbBodies = SpinBox('Number of bodies', 'Number of bodies', 2, 2, 5, 1)
        self.LayoutV1.addWidget(self.NbBodies, alignment=Qt.AlignmentFlag.AlignLeft)
        self.NbBodiesValue = self.NbBodies.SpinParam.value()
        self.NbBodies.SpinParam.valueChanged.connect(self.ChangeNbBodies)

        self.NbOrbitsValue = self.NbBodiesValue-1
        self.NbOrbits = QLabel(f'=>    {self.NbOrbitsValue} Orbit')
        self.NbBodies.Layout.addWidget(self.NbOrbits)

        self.SystDist = DoubleSpinBox('System distance', 'Distance from us of the studied system', 0, 0, None, 1, 2)
        self.LayoutV1.addWidget(self.SystDist)

        self.SystDistUnit = ComboBox(None, 'Unit', ['pc', 'mas'])
        self.SystDist.Layout.addWidget(self.SystDistUnit, alignment=Qt.AlignmentFlag.AlignLeft)

        self.CheckJitter = CheckBox('Jitter :')
        self.CheckJitter.setEnabled(False)
        self.CheckJitter.CheckParam.stateChanged.connect(self.EnablePriorJitterOrNot)
        self.LayoutV1.addWidget(self.CheckJitter, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Jitter = DoubleSpinBox(None, 'Jitter [m/s]', 0, 0, 2147483647, None, 2)
        self.Jitter.setEnabled(False)
        self.CheckJitter.Layout.addWidget(self.Jitter, alignment=Qt.AlignmentFlag.AlignLeft)
        self.LblPlusV0 = QLabel('+')
        self.LblPlusV0.setEnabled(False)
        self.CheckJitter.Layout.addWidget(self.LblPlusV0)
        self.V0 = DoubleSpinBox(None, 'Systematic error V0 [m/s]', 0, 0, 2147483647, None, 2)
        self.V0.setEnabled(False)
        self.CheckJitter.Layout.addWidget(self.V0, alignment=Qt.AlignmentFlag.AlignLeft)

        self.CheckUnivVar = CheckBox('Universal variables', 'If the eccentricity of an orbit is very large or if it is poorly constrained')
        self.LayoutV1.addWidget(self.CheckUnivVar)
        self.CheckUnivVar.CheckParam.stateChanged.connect(self.EnableUnivVarOrNot)

        self.LayoutV1.addWidget(Delimiter(Title='Range of parameters :'))

        # self.LayoutV1.setFixedWidth(300)

        self.PMin = DoubleSpinBox('Period', 'Mimimum of orbits period [day]', None, 0, None, 1, 2)
        self.PMin.LblParam.setMinimumWidth(500)
        self.LayoutV1.addWidget(self.PMin)
        self.PMin.Layout.addSpacing(15)
        self.PMin.Layout.addWidget(QLabel('<->'))
        self.PMin.Layout.addSpacing(10)
        self.PMax = DoubleSpinBox(None, 'Maximum of orbits period [day]', None, 0, None, 1, 2)
        self.PMin.Layout.addWidget(self.PMax)

        self.aMin = DoubleSpinBox('Semi-major axis', 'Mimimum of orbits semi-major axis [AU]', None, 0, None, 1, 2)
        self.LayoutV1.addWidget(self.aMin)
        self.aMin.Layout.addSpacing(15)
        self.aMin.Layout.addWidget(QLabel('<->'))
        self.aMin.Layout.addSpacing(10)
        self.aMax = DoubleSpinBox(None, 'Maximum of orbits semi-major axis [AU]', None, 0, None, 1, 2)
        self.aMin.Layout.addWidget(self.aMax)

        self.PeriMin = DoubleSpinBox('Periastron', 'Mimimum of orbits periastron [AU]', None, 0, None, 1, 2)
        self.LayoutV1.addWidget(self.PeriMin)
        self.PeriMin.Layout.addSpacing(15)
        self.PeriMin.Layout.addWidget(QLabel('<->'))
        self.PeriMin.Layout.addSpacing(10)
        self.PeriMax = DoubleSpinBox(None, 'Maximum of orbits periastron [AU]', None, 0, None, 1, 2)
        self.PeriMin.Layout.addWidget(self.PeriMax)

        self.eMin = DoubleSpinBox('Eccentricity', 'Mimimum of orbits eccentricity', None, 0, 10, 0.1, 2)
        self.LayoutV1.addWidget(self.eMin)
        self.eMin.Layout.addSpacing(15)
        self.eMin.Layout.addWidget(QLabel('<->'))
        self.eMin.Layout.addSpacing(10)
        self.eMax = DoubleSpinBox(None, 'Maximum of orbits eccentricity', None, 0, 10, 0.1, 2)
        self.eMin.Layout.addWidget(self.eMax)

        self.LayoutV1.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.LayoutV1.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.LayoutV1.setSpacing(0)
        self.LayoutH.addLayout(self.LayoutV1)

        self.LayoutH.addSpacing(20)

        self.LayoutV2 = QVBoxLayout()

        self.LayoutV2.addWidget(Delimiter(Title = 'Orbits parameters first guess :'))

        self.FirstGuessCenterMass = DoubleSpinBox('Center mass', 'First guess of center mass [Msun]', 0, 0, None, 1, 2)
        self.LayoutV2.addWidget(self.FirstGuessCenterMass)

        # self.FirstGuessCenterMassUnit = ComboBox(None, 'Unit', ['Msun', 'Mjup'])
        # self.FirstGuessCenterMass.Layout.addWidget(self.FirstGuessCenterMassUnit)

        self.FirstGuessOtherMassUnit = ComboBox('Other bodies mass unit', 'Unit of mass of the other bodies', ['Mjup', 'Msun'])
        self.FirstGuessOtherMassUnit.ComboParam.currentIndexChanged.connect(self.OtherMassUnitChange)
        self.LayoutV2.addWidget(self.FirstGuessOtherMassUnit)

        self.TablePriors = QTableWidget()

        self.LabelParams = ['m [Mjup]', 'a [AU]', 'e', 'i [°]', 'w [°]', 'W [°]', 'tp [MJD]']
        self.EnableUnivVarOrNot(self.CheckUnivVar.CheckParam.isChecked())

        self.TablePriors.setStatusTip('First guess of orbits parameters of each bodies.')
        self.TablePriors.setRowCount(self.NbOrbitsValue)
        self.LayoutV2.addWidget(self.TablePriors, alignment=Qt.AlignmentFlag.AlignVCenter)
        for i in range(self.NbOrbitsValue):
            for j in range(len(self.LabelParams)):
                # if i==0 and j!=0: self.TablePriors.setItem(i, j, QTableWidgetItem('X'))
                # else: 
                self.TablePriors.setItem(i, j, QTableWidgetItem('0.'))
                self.TablePriors.item(i, j).setTextAlignment(Qt.AlignmentFlag.AlignCenter)

        # self.ContainerTest = QWidget()
        # self.LayoutTest = QHBoxLayout()

        # self.TablePriors.setCellWidget(0, 0, self.FirstGuessCenterMassUnit)
        self.TablePriors.itemChanged.connect(self.ValidationItemTbl)
        self.TablePriors.cellClicked.connect(self.SaveOldTextTbl)

        self.RefTime = DateAndMJDEdit('Time reference', 'Reference of the time to count the orbit phase, preferably in the middle of data times')
        self.LayoutV2.addWidget(self.RefTime, alignment=Qt.AlignmentFlag.AlignLeft)

        self.LayoutV2.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.LayoutH.addLayout(self.LayoutV2)

        self.WidgetH.setLayout(self.LayoutH)
        self.Layout.addWidget(self.WidgetH)

        self.ContainerMass = QWidget()
        self.LayoutMass = QVBoxLayout()

        self.BtnNewPriorMass = QPushButton('+')
        self.LayoutMass.addWidget(self.BtnNewPriorMass)
        self.BtnNewPriorMass.clicked.connect(self.AddNewPriorMass)

        self.ContainerMass.setLayout(self.LayoutMass)
        self.Layout.addWidget(self.ContainerMass)

        self.NbWidgetsWithoutPriorMass = self.LayoutMass.count()

        self.Layout.setSpacing(0)

        
    def AddNewPriorMass(self):
        self.PriorMass = PriorMassClass(self.NbBodiesValue)
        self.ListPriorMass.append(self.PriorMass)
        self.c += 1
        self.PriorMass.Id = self.c
        self.ListPriorMassId.append(self.PriorMass.Id)
        self.LayoutMass.addWidget(self.PriorMass, alignment=Qt.AlignmentFlag.AlignLeft)
        self.PriorMass.SignalDelPrior.connect(self.DelThisPriorMass)
        self.LayoutMass.setSpacing(0)

    def DelThisPriorMass(self, Id):
        index = self.ListPriorMassId.index(Id)
        self.ListPriorMass.pop(index)
        self.ListPriorMassId.pop(index)

        DelAllWidgetsBtw(self.LayoutMass, self.NbWidgetsWithoutPriorMass, self.LayoutMass.count())
        
        for i in range(len(self.ListPriorMass)):
            self.LayoutMass.addWidget(self.ListPriorMass[i], alignment=Qt.AlignmentFlag.AlignLeft)

    def ChangeNbBodies(self):
        # Change widgets
        OldValue = self.NbBodiesValue
        self.NbBodiesValue = self.NbBodies.SpinParam.value()
        self.NbOrbitsValue = self.NbBodiesValue-1
        if  self.NbOrbitsValue==1: self.NbOrbits.setText(f'=>    {self.NbOrbitsValue} Orbit')
        else: self.NbOrbits.setText(f'=>    {self.NbOrbitsValue} Orbits')

        # Change the number of row in the priors table
        self.TablePriors.setRowCount(self.NbOrbitsValue)

        if self.NbBodiesValue > OldValue:
            for n in range(len(self.LabelParams)):
                self.TablePriors.setItem(self.NbOrbitsValue-1, n, QTableWidgetItem('0.'))
                self.TablePriors.item(self.NbOrbitsValue-1, n).setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            
            for x in self.ListPriorMass:
                x.Layout.insertWidget(x.Layout.indexOf(x.Distrib)-1, QLabel('+'), alignment=Qt.AlignmentFlag.AlignLeft)
                x.Layout.insertWidget(x.Layout.indexOf(x.Distrib)-1, SpinBox(None, 'Coefficient of mass', 0, 0, 1, 1), alignment=Qt.AlignmentFlag.AlignLeft)
                x.Layout.insertWidget(x.Layout.indexOf(x.Distrib)-1, QLabel('m'+str(self.NbBodiesValue)), alignment=Qt.AlignmentFlag.AlignLeft)

        elif self.NbOrbitsValue < OldValue:
            for x in self.ListPriorMass:
                for i in range(3):
                    x.Layout.removeWidget(x.Layout.itemAt(x.Layout.indexOf(x.Distrib)-2).widget())

    def ValidationItemTbl(self):
        if self.TablePriors.currentItem()!=None:
            TextNew = self.TablePriors.currentItem().text()
            try:
                TextNew = float(TextNew)
                if self.TablePriors.currentColumn() not in [3,4,5] and TextNew<0:
                    self.TablePriors.currentItem().setText(self.TextOld)
            except:
                self.TablePriors.currentItem().setText(self.TextOld)
                
    def SaveOldTextTbl(self):
        self.TextOld = self.TablePriors.currentItem().text()

    def EnablePriorJitterOrNot(self, state):
        self.Jitter.setEnabled(state)
        self.V0.setEnabled(state)
        self.LblPlusV0.setEnabled(state)

    def EnableUnivVarOrNot(self, state):
        self.PMin.setEnabled(not state)
        self.aMin.setEnabled(not state)
        self.PeriMin.setEnabled(state)
        if not state:
            self.eMin.SpinParam.setMaximum(1)
            self.eMax.SpinParam.setMaximum(1)
            self.LabelParams[1] = 'a [AU]'
        else:
            self.eMin.SpinParam.setMaximum(10)
            self.eMax.SpinParam.setMaximum(10)
            self.LabelParams[1] = 'q [AU]'

        self.TablePriors.setColumnCount(len(self.LabelParams))
        self.TablePriors.setHorizontalHeaderLabels(self.LabelParams)

    def InputBetaPicValues(self):
        return

    def OtherMassUnitChange(self, index):
        if index==0:
            self.LabelParams[0] = 'm [Mjup]'
        elif index==1: 
            self.LabelParams[0] = 'm [Msun]'

        self.TablePriors.setColumnCount(len(self.LabelParams))
        self.TablePriors.setHorizontalHeaderLabels(self.LabelParams)


class TabStartSet(GeneralTab):
    
    def __init__(self):
        super().__init__()


    def InitWidgets(self):

        self.Layout.addWidget(Delimiter(Title='Options :'))

        self.CheckParallel = CheckBox('Parallelization', 'Parallelization of the simulation algorithm')
        self.Layout.addWidget(self.CheckParallel)
        self.CheckParallel.CheckParam.setChecked(True)

        self.CheckParallel.Layout.setSpacing(60)
        self.NbCores = SpinBox('Number of cores', 'Number of cores use for parallelization', 8, 1, None, 1)
        self.CheckParallel.Layout.addWidget(self.NbCores, alignment=Qt.AlignmentFlag.AlignLeft)
        self.CheckParallel.CheckParam.stateChanged.connect(self.CheckParallelChange)

        self.BtnCreate = QPushButton('Create startup files')
        self.Layout.addWidget(self.BtnCreate, alignment=Qt.AlignmentFlag.AlignRight)

        self.Layout.addWidget(Delimiter(Title='Order :'))

        # self.CheckOrder = CheckBox('Starting order', 'If you want to start with bash order')
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

        self.StartOrder = LineEdit('Shell order', 'Terminal order to start the adjustment', self.ComboOrder.ComboParam.currentText())
        self.Layout.addWidget(self.StartOrder)

        self.BtnSaveOrder = QPushButton('save')
        self.StartOrder.Layout.addWidget(self.BtnSaveOrder)
        self.BtnSaveOrder.clicked.connect(self.SaveOrder)

        self.LblGo = QLabel('./start.sh')
        self.StartOrder.Layout.addWidget(self.LblGo)

        self.BtnStart = QPushButton('Start the simulation')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)

        # self.CheckStartOrderChange(self.CheckOrder.CheckParam.isChecked())

    # def CheckStartOrderChange(self, state):
    #     # self.NbHours.setEnabled(state)
    #     self.ComboOrder.setEnabled(state)
    #     self.StartOrder.setEnabled(state)
    #     self.BtnStart.setEnabled(state)

    def CheckParallelChange(self, state):
        self.NbCores.setEnabled(state)
        self.NbCores.SpinParam.setValue(state*8//2)

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

        
    
