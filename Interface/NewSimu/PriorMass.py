#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# PyQt packages
from PyQt6.QtWidgets import QWidget, QHBoxLayout, QLabel, QPushButton, QCheckBox, QFileDialog
from PyQt6.QtCore import pyqtSignal, Qt

# My packages
from Parameters import *
from UtilsNewSimu import DelAllWidgetsBtw


class PriorMassClass(QWidget):

    SignalDelPrior = pyqtSignal(int)

    def __init__(self, NbBodies):
        super().__init__()

        # Identity
        self.Id = 0
        
        # Layout
        self.Layout = QHBoxLayout()

        self.NbBodies = NbBodies

        self.ButtonDel = QPushButton('-')
        self.Layout.addWidget(self.ButtonDel, alignment=Qt.AlignmentFlag.AlignLeft)
        self.ButtonDel.clicked.connect(lambda: self.SignalDelPrior.emit(self.Id))

        for i in range(1, NbBodies+1):
            self.Layout.addWidget(SpinBox(None, 'Coefficient of mass', 0, 0, 1, 1), alignment=Qt.AlignmentFlag.AlignLeft)
            self.Layout.addWidget(QLabel('m'+str(i)), alignment=Qt.AlignmentFlag.AlignLeft)

            if i != NbBodies: self.Layout.addWidget(QLabel('+'), alignment=Qt.AlignmentFlag.AlignLeft)
            else : self.Layout.addWidget(QLabel('='), alignment=Qt.AlignmentFlag.AlignLeft)
        

        self.Distrib = ComboBox(None, 'Choice of a priori mass combination', ['Distribution', 'Normal', 'Log', 'Linear', 'Fixed'])
        self.Layout.addWidget(self.Distrib, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Distrib.ComboParam.currentIndexChanged.connect(self.AddDistribParams)

        # Widget container  
        self.setLayout(self.Layout)
        


    def AddDistribParams(self):

        DelAllWidgetsBtw(self.Layout, self.Layout.indexOf(self.Distrib)+1, self.Layout.count())

        self.DistribIndex = self.Distrib.ComboParam.currentIndex()

        if self.DistribIndex == 1 or self.DistribIndex == 2: # Normal or Log
            self.Mean = DoubleSpinBox(None, 'Mean of the distribution', 0, 0, 2147483647, 0.01)
            self.Layout.addWidget(self.Mean, alignment=Qt.AlignmentFlag.AlignLeft)

            self.Layout.addWidget(QLabel('+/-'))

            self.SD = DoubleSpinBox(None, 'Standart deviation of the distribution', 0, 0, 2147483647, 0.01)
            self.Layout.addWidget(self.SD)

        elif self.DistribIndex == 3: # Uniform
            self.Min = SpinBox(None, 'Minimum of the distribution', 0, 0, 2147483647, 1)
            self.Layout.addWidget(self.Min)

            self.Layout.addWidget(QLabel('<->'))

            self.Max = SpinBox(None, 'Maximum of the distribution', 0, 0, 2147483647, 1)
            self.Layout.addWidget(self.Max)

        elif self.DistribIndex == 4: # Fixed
            self.Value = DoubleSpinBox(None, 'Value', 0, 0, 2147483647, 0.01)
            self.Layout.addWidget(self.Value)

        if self.DistribIndex != 0: 
            self.PriorUnit = ComboBox(None, 'Unit', ['ms', 'mj'])
            self.Layout.addWidget(self.PriorUnit)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignLeft)




### --- Check --- ###
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    PriorMassWidget = PriorMassClass(4)
    PriorMassWidget.show()
    app.exec() # Application execution