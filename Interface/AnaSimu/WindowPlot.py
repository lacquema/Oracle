#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import sys

# PyQt packages
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QApplication

# Matplotlib packages
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


### --- Canvas Generating --- ###

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self):
        self.fig = Figure()
        super(MplCanvas, self).__init__(self.fig)


### --- Plot Window Generating --- ###

class WindowPlotClass(QMainWindow):

    SignalCloseWindowPlot = pyqtSignal() # initiation of the closeEvent signal

    def __init__(self, ToolName):
        super().__init__()

        # Window characteristics
        self.setWindowTitle(ToolName+': Plot')
        # self.move(self.frameGeometry().topLeft())

        # Layout
        Layout = QVBoxLayout()

        # Canvas initialisation
        self.Canvas = MplCanvas()
        
        # Plot on the Canvas
        Layout.addWidget(self.Canvas)

        # Toolbar
        Toolbar = NavigationToolbar(self.Canvas, self)
        Layout.addWidget(Toolbar)

        # Widget container
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

    # Emition of the CloseEvent signal when the plot window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowPlot.emit() 


if __name__=="__main__":
    app = QApplication(sys.argv)
    WindowPlot = WindowPlotClass('ToolName')
    WindowPlot.show()
    app.exec()