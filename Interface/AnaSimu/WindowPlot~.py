#! /Users/lacquema/Oracle.env/bin/python3
import sys

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QApplication, QSplitter, QMainWindow, QStatusBar

from WidgetParam import WidgetParam
from WidgetPlot import WidgetPlot

class WindowPlot(QMainWindow):

    SignalCloseWindowPlot = pyqtSignal() # initiation of the closeEvent signal

    def __init__(self, ToolName):
        super().__init__()

        # Window settings
        self.setWindowTitle(ToolName)

        # Splitter widget
        self.Splitter = QSplitter()

        # Parameters widget
        self.WidgetParam = WidgetParam()
        self.Splitter.addWidget(self.WidgetParam)

        # Plotting widget
        self.WidgetPlot = WidgetPlot()
        self.Splitter.addWidget(self.WidgetPlot)

        # Status bar
        self.setStatusBar(QStatusBar(self))

        # Container
        self.setCentralWidget(self.Splitter)


    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowPlot.emit() 


if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    Window = WindowPlot('test')
    Window.show()
    app.exec() # Application execution