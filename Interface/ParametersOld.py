#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import sys

# PyQt packages
from PyQt6.QtWidgets import QWidget, QHBoxLayout, QLabel, QLineEdit, QComboBox, QSpinBox, QApplication, QCheckBox, QDateEdit
from PyQt6.QtGui import QDoubleValidator, QIntValidator


### --- Parameters Generating --- ###

class ParamsClass(QWidget):

    def __init__(self, ParamName):
        super().__init__()

        # Layout
        self.Layout = QHBoxLayout()

        # Parameters label
        if ParamName != None:
            self.LblParam = QLabel("{} :".format(ParamName))
            self.Layout.addWidget(self.LblParam)

        # Widget container
        self.setLayout(self.Layout) # ParamsClass is directly the widget container




# General LineEdit
class LineEdit(ParamsClass):

    def __init__(self, ParamName, ParamStatus, ParamDefault):
        super().__init__(ParamName)
        self.EditParam = QLineEdit()
        self.EditParam.setText(str(ParamDefault))
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
class ComboBox(ParamsClass):
    def __init__(self, ParamName, ParamStatus, ParamItems):
        super().__init__(ParamName)
        self.ComboParam = QComboBox()
        self.ComboParam.addItems(ParamItems)
        if ParamStatus != None: self.ComboParam.setStatusTip(ParamStatus)
        self.Layout.addWidget(self.ComboParam)



# Gerenal SpinBox
class SpinBox(ParamsClass):
    def __init__(self, ParamName, ParamStatus, ParamMin=None, ParamMax=None, ParamDefault=None):
        super().__init__(ParamName)
        self.SpinParam = QSpinBox()
        if ParamMin!=None: self.SpinParam.setMinimum(ParamMin)
        if ParamMax!=None: self.SpinParam.setMaximum(ParamMax)
        if ParamDefault != None: self.SpinParam.setValue(ParamDefault)
        if ParamStatus != None: self.SpinParam.setStatusTip(ParamStatus)
        self.Layout.addWidget(self.SpinParam)


# Gerenal CheckBox
class CheckBox(ParamsClass):
    def __init__(self, ParamName, ParamStatus):
        super().__init__(None)
        self.CheckParam = QCheckBox(ParamName)
        if ParamStatus != None: self.CheckParam.setStatusTip(ParamStatus)
        self.Layout.addWidget(self.CheckParam)


# Gerenal DateEdit
class DateEdit(ParamsClass):
    def __init__(self, ParamName, ParamStatus):
        super().__init__(ParamName)
        self.EditParam = QDateEdit(calendarPopup=True)
        if ParamStatus != None: self.EditParam.setStatusTip(ParamStatus)
        self.Layout.addWidget(self.EditParam)




# Check
if __name__=="__main__":
    ParamName, ParamStatus, ParamDefault = 'Xmin (AU)','Minimum abscissa','-100'
    app = QApplication(sys.argv) # Application creation
    ParamWidget = LineEdit(ParamName, ParamStatus, ParamDefault)
    ParamWidget.show()
    app.exec() # Application execution
