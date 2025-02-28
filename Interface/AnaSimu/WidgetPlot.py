#! /Users/lacquema/Oracle.env/bin/python3


### --- Packages --- ###

# Transverse packages
import sys

# PyQt packages
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

class WidgetPlot(QWidget):

    def __init__(self):
        super().__init__()

        # Layout
        self.Layout = QVBoxLayout()

        # Canvas initialisation
        self.Canvas = MplCanvas()

        # Toolbar
        Toolbar = NavigationToolbar(self.Canvas, self)
        self.Layout.addWidget(Toolbar)
        
        # Plot on the Canvas
        self.Layout.addWidget(self.Canvas)

        # Widget container
        self.setLayout(self.Layout)



if __name__=="__main__":
    app = QApplication(sys.argv)
    Plot = WidgetPlot()
    Window = QMainWindow()
    Window.setCentralWidget(Plot)
    Window.show()
    app.exec()