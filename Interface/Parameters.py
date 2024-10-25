#! /Users/lacquema/Oracle.env/bin/python3

### --- Packages --- ###

# Transverse packages
import sys
from Utils import *

# PyQt packages
from PyQt6.QtWidgets import QVBoxLayout, QProgressBar, QPushButton, QDateEdit, QCheckBox, QWidget, QHBoxLayout, QLabel, QLineEdit, QComboBox, QSpinBox, QApplication, QDoubleSpinBox, QFileDialog
from PyQt6.QtGui import QDoubleValidator, QIntValidator
from PyQt6.QtCore import Qt, QDate, QDateTime


### --- Parameters Generating --- ###

class GeneralParam(QWidget):

    def __init__(self, ParamName):
        super().__init__()

        # Layout
        self.Layout = QHBoxLayout()

        # Parameters label
        if ParamName != None:
            self.LblParam = QLabel("{} :".format(ParamName))
            self.Layout.setSpacing(20)
            self.Layout.addWidget(self.LblParam)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignLeft)

        # Widget container
        self.setLayout(self.Layout) # GeneralParam is directly the widget container



# General LineEdit
class LineEdit(GeneralParam):

    def __init__(self, ParamName, ParamStatus, ParamDefault):
        super().__init__(ParamName)
        self.EditParam = QLineEdit()
        if ParamDefault != None: self.EditParam.setText(str(ParamDefault))
        if ParamStatus != None: self.EditParam.setStatusTip(ParamStatus)
        self.Layout.addWidget(self.EditParam)


# LineEdit with float validator
class LineEditFloat(LineEdit):

    def __init__(self, ParamName, ParamStatus, ParamDefault):
        super().__init__(ParamName, ParamStatus, ParamDefault)
        self.EditParam.setValidator(QDoubleValidator()) # Only float accepted


# LineEdit with int validator
class LineEditInt(LineEdit):

    def __init__(self, ParamName, ParamStatus, ParamDefault):
        super().__init__(ParamName, ParamStatus, ParamDefault)
        self.EditParam.setValidator(QIntValidator()) # Only int accepted


# LineEdit with positive int validator
class LineEditIntPos(LineEditInt):

    def __init__(self, ParamName, ParamStatus, ParamDefault):
        super().__init__(ParamName, ParamStatus, ParamDefault)
        self.EditParam.textChanged.connect(self.ValidatorPos) # Only positive int accepted

    def ValidatorPos(self):
        text = self.EditParam.text()
        if len(text) == 1:
            if text[0] == '-':
                self.EditParam.setText('')
        if len(text) > 1:
            if text[0] == '-':
                self.EditParam.setText(text[1:])
                self.EditParam.setCursorPosition(0)


# LineEdit with positive float validator
class LineEditFloatPos(LineEditFloat):

    def __init__(self, ParamName, ParamStatus, ParamDefault):
        super().__init__(ParamName, ParamStatus, ParamDefault)
        self.EditParam.textChanged.connect(self.ValidatorPos) # Only positive float accepted

    def ValidatorPos(self):
        text = self.EditParam.text()
        if len(text) == 1:
            if text[0] == '-':
                self.EditParam.setText('')
        if len(text) > 1:
            if text[0] == '-':
                self.EditParam.setText(text[1:])
                self.EditParam.setCursorPosition(0)


# LineEdit with variable validator
class LineEditValidated(LineEdit):

    def __init__(self, ParamName, ParamStatus, ParamDefault, ListValidCharacters):
        super().__init__(ParamName, ParamStatus, ParamDefault)
        self.EditParam.textChanged.connect(lambda: self.Validator(ListValidCharacters)) # Only characters in ListValidCharacters accepted

    def Validator(self, ListValidCharacters):
        TextNew = self.EditParam.text()
        c = 0
        for i in range(len(TextNew)):
            if TextNew[i] not in  ListValidCharacters:
                c=+1
                break
        if c != 0:
            self.EditParam.setText(self.TextOld)
            self.EditParam.setCursorPosition(i)
        else:
            self.TextOld = TextNew




# Gerenal ComboBox
class ComboBox(GeneralParam):

    def __init__(self, ParamName, ParamStatus, ParamItems, ParamIndexDefault=None):
        super().__init__(ParamName)
        self.ComboParam = QComboBox()
        self.ComboParam.addItems(ParamItems)
        if ParamStatus != None: self.ComboParam.setStatusTip(ParamStatus)
        if ParamIndexDefault !=None: self.ComboParam.setCurrentIndex(ParamIndexDefault)
        self.Layout.addWidget(self.ComboParam)



# Gerenal SpinBox
class SpinBox(GeneralParam):
     
     def __init__(self, ParamName, ParamStatus, ParamDefault=None, ParamMin=None, ParamMax=None, ParamIncrement=None):
        super().__init__(ParamName)
        self.SpinParam = QSpinBox()
        if ParamStatus != None: self.SpinParam.setStatusTip(ParamStatus)
        if ParamMin!=None: self.SpinParam.setMinimum(ParamMin)
        else: self.SpinParam.setMinimum(-2147483647)
        if ParamMax!=None: self.SpinParam.setMaximum(ParamMax)
        else: self.SpinParam.setMaximum(2147483647)
        if ParamDefault != None: self.SpinParam.setValue(ParamDefault)
        if ParamIncrement != None: self.SpinParam.setSingleStep(ParamIncrement)
        self.Layout.addWidget(self.SpinParam)



# Float SpinBox
class DoubleSpinBox(GeneralParam):
     
     def __init__(self, ParamName, ParamStatus, ParamDefault=None, ParamMin=None, ParamMax=None, ParamIncrement=None, ParamPrecision=None):
        super().__init__(ParamName)
        self.SpinParam = QDoubleSpinBox()
        if ParamStatus != None: self.SpinParam.setStatusTip(ParamStatus)
        if ParamMin!=None: self.SpinParam.setMinimum(ParamMin)
        else: self.SpinParam.setMinimum(-2147483647)
        if ParamMax!=None: self.SpinParam.setMaximum(ParamMax)
        else: self.SpinParam.setMaximum(2147483647)
        if ParamDefault != None: self.SpinParam.setValue(ParamDefault)
        if ParamIncrement != None: self.SpinParam.setSingleStep(ParamIncrement)
        if ParamPrecision != None: self.SpinParam.setDecimals(ParamPrecision)
        self.Layout.addWidget(self.SpinParam)


# Gerenal CheckBox
class CheckBox(GeneralParam):
    def __init__(self, ParamName, ParamStatus=None):
        super().__init__(None)
        self.CheckParam = QCheckBox(ParamName)
        if ParamStatus != None: self.CheckParam.setStatusTip(ParamStatus)
        self.Layout.addWidget(self.CheckParam)



class PathBrowser(GeneralParam):
    def __init__(self, ParamName, ParamStatus, TargetType):
        super().__init__(ParamName)
        # Edit path
        self.EditPath = QLineEdit()
        self.EditPath.setMinimumWidth(200)
        if ParamStatus != None: self.EditPath.setStatusTip(ParamStatus)
        self.Layout.addWidget(self.EditPath)

        # Btn browse
        self.BtnBrowse = QPushButton('Browse')
        self.Layout.addWidget(self.BtnBrowse)

        if TargetType == 0:
            self.BtnBrowse.clicked.connect(self.DialBrowseDirectory)
        elif TargetType == 1:
            self.BtnBrowse.clicked.connect(self.DialBrowseFile)

    def DialBrowseDirectory(self):
        self.Directory = QFileDialog.getExistingDirectory(self)
        if len(self.Directory) !=0:
            self.EditPath.setText(self.Directory+'/')

    def DialBrowseFile(self):
        self.File = QFileDialog.getOpenFileName(self)[0]
        if len(self.File) !=0:
            self.EditPath.setText(self.File)


class Delimiter(QWidget):
    def __init__(self, Height=None, Width=None, Title=None):
        super().__init__()

        # Layout
        self.Layout = QVBoxLayout()

        # Delimiter
        self.Delim = QProgressBar()
        self.Layout.addWidget(self.Delim, alignment=Qt.AlignmentFlag.AlignBaseline)

        if Height!=None: self.Delim.setFixedHeight(Height)
        else: self.Delim.setFixedHeight(3)

        if Width!=None: self.Delim.setFixedWidth(Width)

        # Title
        if Title!=None: 
            self.Title = QLabel(Title)
            self.Title.setStyleSheet('QLabel{font: bold italic}')
            self.Layout.addWidget(self.Title, alignment=Qt.AlignmentFlag.AlignCenter)

        self.Layout.setSpacing(5)

        self.setLayout(self.Layout)


# Gerenal DateEdit
class DateEdit(GeneralParam):
    def __init__(self, ParamName, ParamStatus):
        super().__init__(ParamName)
        self.EditParam = QDateEdit(calendarPopup=True)
        if ParamStatus != None: self.EditParam.setStatusTip(ParamStatus)
        self.Layout.addWidget(self.EditParam)
        

class DateAndMJDEdit(GeneralParam):
    def __init__(self, ParamName, ParamStatus):
        super().__init__(ParamName)

        # Current date
        self.CurDate = QDateTime.currentDateTime().date().getDate()
        self.CurMJD = DatetoMJD(*self.CurDate)

        # Date on calendar
        self.DateWidget = DateEdit(None, ParamStatus)
        self.DateWidget.EditParam.setDate(QDate(*self.CurDate))
        self.Layout.addWidget(self.DateWidget)

        # Date on MJD
        self.MJDWidget = SpinBox(None, 'Date [MJD]', self.CurMJD, 0, None)
        self.DateWidget.Layout.addWidget(self.MJDWidget)

        # Connections
        self.DateWidget.EditParam.dateChanged.connect(self.DateChanged)
        self.MJDWidget.SpinParam.textChanged.connect(self.MJDChanged)

        self.setLayout(self.Layout)
    
    def DateChanged(self):
        Date = self.DateWidget.EditParam.date().getDate() # YY MM JJ
        MJD = DatetoMJD(*Date)
        self.MJDWidget.SpinParam.setValue(MJD)

    def MJDChanged(self):
        MJD = self.MJDWidget.SpinParam.value()
        Date = MJDtoDate(MJD)
        self.DateWidget.EditParam.setDate(QDate(*Date))
    




# Check
if __name__=="__main__":
    ParamName, ParamStatus, ParamDefault = 'Xmin (AU)','Minimum abscissa','-100'
    app = QApplication(sys.argv) # Application creation
    ParamWidget = DateAndMJDEdit('Date', 'Date of wanted observation')
    ParamWidget.show()
    app.exec() # Application execution
