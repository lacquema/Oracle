#! /Users/lacquema/Oracle.env/bin/python3


### --- Packages --- ###

# Transverse packages
import sys
import os

# PyQt packages
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QPushButton, QWidget, QStatusBar, QApplication
from PyQt6.QtGui import QIcon


### --- Parameters Window Generating --- ###

class WidgetParam(QWidget):
    
    def __init__(self):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Layout
        self.Layout = QVBoxLayout()
        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)

        # Reset button
        self.BtnReset = QPushButton('Reset')
        self.BtnReset.setStatusTip('Reset')
        self.Layout.addWidget(self.BtnReset)

        # Refresh button
        self.BtnRefresh = QPushButton()
        self.BtnRefresh.setIcon(QIcon(f'{self.DirPath}/Items/arrowCircle.png'))
        self.BtnRefresh.setStatusTip('Refresh')
        self.Layout.addWidget(self.BtnRefresh)

        # Widget container
        self.setLayout(self.Layout)

       
# Check
if __name__=="__main__":
    app = QApplication(sys.argv)
    Plot = WidgetParam()
    Window = QMainWindow()
    Window.setCentralWidget(Plot)
    Window.setStatusBar(QStatusBar())
    Window.show()
    app.exec()
    
