#! /Users/lacquema/Oracle.env/bin/python3
import sys

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QApplication, QSplitter, QMainWindow, QStatusBar, QSizePolicy

from WidgetParam import WidgetParam
from WidgetPlot import WidgetPlot

class WindowPlot(QMainWindow):

    SignalCloseWindowPlot = pyqtSignal()  # initiation of the closeEvent signal

    def __init__(self, ToolName):
        super().__init__()

        # Window settings
        self.setWindowTitle(ToolName)

        # Splitter widget
        self.Splitter = QSplitter()

        # Parameters widget
        self.WidgetParam = WidgetParam()
        self.Splitter.addWidget(self.WidgetParam)

        # Plotting widgets
        self.WidgetPlots = []

        # Status bar
        self.setStatusBar(QStatusBar(self))

        # Container
        self.setCentralWidget(self.Splitter)


    def add_WidgetPlot(self, plot, layout=None, xlim=False, ylim=False, zlim=False, azim=False, elev=False, xlabel=True, ylabel=True, zlabel=True, title=True, legend=True):
        """
        Creates a new WidgetPlot connected to the WidgetParam.

        The parameters will be save if True.
        Initial values provided during plot creation are overwritten if modified.
        """
        # print(title, 'on add_WidgetPlot')
        widget_plot = WidgetPlot(plot, xlim, ylim, zlim, azim, elev, xlabel, ylabel, zlabel, title, legend)
        self.WidgetPlots.append(widget_plot)
        if layout is None:
            self.Splitter.addWidget(widget_plot)
        else:
            layout.addWidget(widget_plot)
        return widget_plot
    
    
    # Emission of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowPlot.emit()


if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    Window = WindowPlot('test')
    Window.show()
    app.exec() # Application execution