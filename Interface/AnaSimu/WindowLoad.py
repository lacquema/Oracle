#! /Users/lacquema/ByeGildas/bin/python3

# Transverse packages
import sys
import os

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QLabel, QStatusBar, QWidget, QApplication
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt 

### --- Loading Window Generating --- ###
class LoadWindowClass(QMainWindow):
    
    def __init__(self):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Window settings
        self.setWindowTitle("...Load Window...")
        self.setFixedSize(850, 450)

        # Background image
        self.setStyleSheet(f"background-image: url({self.DirPath}/Items/LoadingBackground.png)") 

        # Layout intialisation
        Layout = QVBoxLayout()
    
        # Title
        Title = QLabel('<b>Oracle analysis data interface</b>')
        Title.setFont(QFont('Helvetica', 30))
        Layout.addWidget(Title, alignment=Qt.AlignmentFlag.AlignCenter)

        # Put the title in the top of the window
        Layout.addSpacing(370)

        # Credits in status bar
        StatusBar = QStatusBar(self)
        StatusBar.addWidget(QLabel(' Version 1.1, 2024, IPAG, Antoine Lacquement'))
        self.setStatusBar(StatusBar)

        # Widget container
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

        # Showing
        self.show()

if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    LoadWin = LoadWindowClass()
    app.processEvents() # Continue the program
    sys.exit(app.exec()) # Application execution