#! /Users/lacquema/Oracle.env/bin/python3
import sys
import os

sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__)+'/AnaSimu')
sys.path.append(os.path.dirname(__file__)+'/NewSimu')
sys.path.append(os.path.dirname(__file__)+'/ContSimu')

### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QApplication

# My packages
from WindowMenu import WindowMenuClass
from NewSimu.WindowSetNewSimu import WindowSetNewSimu
from ContSimu.WindowSetContSimu import WindowSetContSimu
from AnaSimu.WindowSetAnaSimu import WindowSetAnaSimu


### --- Main Window Generating --- ###

class MainClass():

    def __init__(self):
        super().__init__()

        # Windows initialisation
        self.WinMenu = WindowMenuClass()
        self.WinSetNewSimu = WindowSetNewSimu()
        self.WinSetNewSimu.SignalCloseWindowSetNewSimu.connect(self.ReOpenWinLoad)
        self.WinSetContSimu = WindowSetContSimu()
        self.WinSetContSimu.SignalCloseWindowSetContSimu.connect(self.ReOpenWinLoad)
        self.WinSetAnaSimu = WindowSetAnaSimu()
        self.WinSetAnaSimu.SignalCloseWindowSetAnaSimu.connect(self.ReOpenWinLoad)
        self.WinSetAnaSimu.ReSignalCloseWindowMain.connect(self.ReOpenWinLoad)
        
        # Showing
        self.WinMenu.show()
        
        # Actions of buttons
        self.WinMenu.BtnNew.clicked.connect(self.OpenWinSetNewSimu)
        self.WinMenu.BtnContinue.clicked.connect(self.OpenWinSetContSimu)
        self.WinMenu.BtnAnalyse.clicked.connect(self.OpenWinSetAnaSimu)

    def OpenWinSetNewSimu(self):
        self.WinMenu.close()
        self.WinSetNewSimu.show()

    def OpenWinSetContSimu(self):    
        self.WinMenu.close()
        self.WinSetContSimu.show()

    def OpenWinSetAnaSimu(self):
        self.WinMenu.close()
        self.WinSetAnaSimu.show()

    def closeEvent(self, e):
        self.WinMenu.show()

    def ReOpenWinLoad(self):
        app.closeAllWindows()
        self.WinMenu.show()





if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowMain = MainClass() # Main window showing
    sys.exit(app.exec()) # Application execution

    