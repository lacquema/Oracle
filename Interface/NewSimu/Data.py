#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import sys

# PyQt packages
from PyQt6.QtWidgets import QWidget, QHBoxLayout, QLabel, QLineEdit, QComboBox, QSpinBox, QApplication, QCheckBox, QDateEdit, QComboBox, QFileDialog, QPushButton, QStatusBar
from PyQt6.QtGui import QDoubleValidator, QIntValidator
from Parameters import *


### --- Parameters Generating --- ###

class DataClass(QWidget):

    def __init__(self, DataType):
        super().__init__()

        # Layout
        self.Layout = QHBoxLayout()

        # Checkbox
        self.CheckData = QCheckBox(str(DataType))
        self.Layout.addWidget(self.CheckData)
        self.CheckData.stateChanged.connect(self.CheckStateChanged)

        self.Layout.addSpacing(50)

        # Format 
        FormatsAstro = ['ID Day Month Year Dec RA dDec dRA Corr',
                        'ID Day Month Year Dec RA dDec dRA',           
                        'ID JD Dec RA dDec dRA Corr',
                        'ID JD Dec RA dDec dRA',
                        'ID Day Month Year Sep PA dSep dPA Corr',
                        'ID Day Month Year Sep PA dSep dPA',           
                        'ID JD Sep PA dSep dPA Corr',
                        'ID JD Sep PA dSep dPA']
        
        FormatsRV = ['ID Day Month Year RV dRV',
                     'ID JD Year RV dRV']
        
        if DataType.split(' ')[1] == 'Astrometry':
            self.FormatData = QComboBox()
            self.FormatData.addItem('Format')
            self.FormatData.addItems(FormatsAstro)
            self.FormatData.setStatusTip('Format of data       Dec,RA,Sep in mas / PA in deg / Corr=correlation coefficient')

        elif DataType.split(' ')[1] == 'RV':
            self.FormatData = QComboBox()
            self.FormatData.addItem('Format')
            self.FormatData.addItems(FormatsRV)
            self.FormatData.setStatusTip('Format of data       RV in km/s')

        self.Layout.addWidget(self.FormatData)
        self.FormatData.setEnabled(False)

        self.Layout.addSpacing(50)

        # # Path of data file
        # self.PathData = PathBrowser(None, 'Path of the data file', 1)
        # self.Layout.addWidget(self.PathData)
        
        
        # Buttons behaviors
        self.FormatData.setEnabled(False)
        # self.PathData.setEnabled(False)

        # Widget container
        self.setLayout(self.Layout) # ParamsClass is directly the widget container

    def CheckStateChanged(self):
        self.FormatData.setEnabled(self.CheckData.isChecked())
        # self.PathData.setEnabled(self.CheckData.isChecked())
            





# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    DataWidget = DataClass('Relative RV')
    DataWidget.show()
    app.exec() # Application execution
