#! /Users/lacquema/ByeGildas/bin/python3

# Transverse packages
import sys
import os

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QLabel, QStatusBar, QWidget, QApplication, QPushButton
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt 

### --- Loading Window Generating --- ###
class WindowMenuClass(QMainWindow):
    
    def __init__(self):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Window settings
        self.setWindowTitle("Oracle Launcher")
        self.setFixedSize(850, 450)

        # Background image
        self.setStyleSheet(f"background-image: url({self.DirPath}/Items/LoadingBackground.png)") 
        # button.setStyleSheet('QPushButton {background-color: #A3C1DA; color: red;}')

        # Layout intialisation
        Layout = QVBoxLayout()

        Layout.addSpacing(20)

        # Button new simulation
        self.BtnNew = QPushButton('New Simulation')
        # Front = QFont('', 20, italic=True)
        # self.BtnNew.setFont(Front)
        self.BtnNew.setFixedSize(200, 40)
        self.BtnNew.setStyleSheet('QPushButton {background-color: grey; color: white; font: italic 15px;}')
        Layout.addWidget(self.BtnNew, alignment=Qt.AlignmentFlag.AlignHCenter)

        Layout.addSpacing(20)

        # Button simulation Continue
        self.BtnContinue = QPushButton('Continue')
        # self.BtnContinue.setFont(Front)
        self.BtnContinue.setFixedSize(200, 40)
        self.BtnContinue.setStyleSheet('QPushButton {background-color: grey; color: white; font: italic 15px;}')
        Layout.addWidget(self.BtnContinue, alignment=Qt.AlignmentFlag.AlignHCenter)

        Layout.addSpacing(95)

        # Button analyse simulation
        self.BtnAnalyse = QPushButton('Analyse')
        # self.BtnAnalyse.setFont(Front)
        self.BtnAnalyse.setFixedSize(200, 40)
        self.BtnAnalyse.setStyleSheet('QPushButton {background-color: grey; color: white; font: italic 15px;}')
        Layout.addWidget(self.BtnAnalyse, alignment=Qt.AlignmentFlag.AlignHCenter)

        Layout.addSpacing(500)

        # Credits in status bar
        StatusBar = QStatusBar(self)
        StatusBar.addWidget(QLabel(' Version 0, 2024, IPAG, Antoine Lacquement'))
        self.setStatusBar(StatusBar)

        # Widget container
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

        # Showing
        self.show()

if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    LoadWin = WindowMenuClass() # Loading window showing
    sys.exit(app.exec()) # Application execution