#! /Users/lacquema/Oracle.env/bin/python3
import sys

from PyQt6.QtGui import QFileSystemModel
from PyQt6.QtCore import QDir, Qt
from PyQt6.QtWidgets import QApplication, QColumnView, QSplitter, QTreeView, QMainWindow

class WindowWithFinder(QMainWindow):

    def __init__(self):
        super().__init__()

        # Splitter widget
        self.Splitter = QSplitter()

        # Model of finder
        self.Model = QFileSystemModel()
        self.Model.setRootPath(QDir.rootPath())

        # Finder widget
        self.Finder = QTreeView()
        self.Finder.setModel(self.Model)

        # Just show folders name
        self.Finder.setHeaderHidden(True)
        self.Finder.setColumnHidden(1, True)
        self.Finder.setColumnHidden(2, True)
        self.Finder.setColumnHidden(3, True)

        # Finder's look
        self.Finder.setFixedWidth(300)

        # Set the path to home
        self.Finder.setRootIndex(self.Model.index(QDir.homePath()))

        # Finder's widget inclusion
        self.Splitter.addWidget(self.Finder)

        self.setCentralWidget(self.Splitter)


if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowWithFinder()
    WindowParam.show()
    app.exec() # Application execution

