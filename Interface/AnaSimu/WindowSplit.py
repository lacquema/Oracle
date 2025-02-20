#! /Users/lacquema/Oracle.env/bin/python3
import sys

from PyQt6.QtGui import QFileSystemModel
from PyQt6.QtCore import QDir, Qt
from PyQt6.QtWidgets import QApplication, QColumnView, QSplitter, QTreeView, QMainWindow, QPushButton

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from Tools import *


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self):
        self.fig = Figure()
        super(MplCanvas, self).__init__(self.fig)



class WindowSpliter(QMainWindow):

    def __init__(self):
        super().__init__()

        # Splitter widget
        self.Splitter = QSplitter()

        
        # # Param layout initialisation
        # self.Param_Layout = QVBoxLayout()

        # # self.Param_Layout.addWidget(QPushButton('Bonjour'))
        
        # # Plot layout initialisation
        # self.Plot_Layout = QVBoxLayout()

        # # Canvas initialisation
        # self.Canvas = MplCanvas()
        
        # # Plot on the Canvas
        # self.Plot_Layout.addWidget(self.Canvas)

        # # Toolbar
        # Toolbar = NavigationToolbar(self.Canvas, self)
        # self.Plot_Layout.addWidget(Toolbar)

        # # Widget container
        # Param_Container = QWidget()
        # Param_Container.setLayout(self.Param_Layout)
        # self.Splitter.addWidget(Param_Container)

        # Plot_Container = QWidget()
        # Plot_Container.setLayout(self.Plot_Layout)
        # self.Splitter.addWidget(Plot_Container)


        # self.Subplot = self.Canvas.fig.add_subplot(111)



        # Finder's widget inclusion

        self.setCentralWidget(self.Splitter)


if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowWithFinder()
    WindowParam.show()
    app.exec() # Application execution