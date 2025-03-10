#! /Users/lacquema/ByeGildas/bin/python3

### --- Packages --- ###

# Transverse packages
import sys
import numpy as np
from random import random, randint
import corner
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import matplotlib.pyplot as plt

# PyQt packages
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QDateEdit, QGroupBox, QGridLayout
from PyQt6.QtCore import QDateTime, QDate, QSize, QTimer
from Utils import date_to_jd, jd_to_mjd, mjd_to_jd, jd_to_date

# My packages
from Parameters import *
from BestOrbits import BestOrbitsClass
from SelectOrbits import SelectOrbitsClass

from WindowPlot import WindowPlot
import itertools


### --- Tools Generating --- ###

class GeneralToolClass(QWidget):

    def __init__(self, ToolName, ToolStatus, InputData, OutputParams, SelectOrbitsParams, SelectOrbitsEllipses, BestOrbitsParams, BestOrbitsEllipses):
        super().__init__()

        self.colorList = ['black', 'blue', 'red', 'green', 'orange', 'pink']

        # Layout
        Layout = QHBoxLayout()

        # Label
        LblTool = QLabel("{} :".format(ToolName))
        LblTool.setStatusTip(ToolStatus)
        Layout.addWidget(LblTool)

        # Plot button
        self.BtnPlot = QPushButton('Plot')
        self.BtnPlot.clicked.connect(self.Toggle_WindowPlot)
        Layout.addWidget(self.BtnPlot)

        # Initialisation of parameters windows
        self.WindowPlot = WindowPlot(ToolName)
        self.WindowPlot.SignalCloseWindowPlot.connect(lambda: self.BtnPlot.setEnabled(True))  # Enable button when window is closed

        # Connections between parameters and plot
        self.WindowPlot.WidgetParam.BtnReset.clicked.connect(self.ResetParams)
        self.WindowPlot.WidgetParam.BtnRefresh.clicked.connect(self.Refresh_ActivePlots)

        # Initialisation of Data
        if InputData is not None:
            (self.NbInputData_RA, self.I_RA, self.MJD_RA, self.JJ_RA, self.MM_RA, self.YY_RA, self.Ra, self.Dec, self.DRa, self.DDec, self.Corr_DecRa, self.Sep, self.Pa, self.DSep, self.DPa, self.Corr_SepPa, self.Source_RA) = InputData

        if OutputParams is not None:
            (self.NbBodies, self.NbOrbits, self.P, self.a, self.e, self.i, self.w, self.W, self.tp, self.m, self.Mdyn, self.Chi2, self.map) = OutputParams

        if SelectOrbitsParams is not None:
            (self.NbBodies, self.NbSelectOrbits, self.SelectP, self.Selecta, self.Selecte, self.Selecti, self.Selectw, self.SelectW, self.Selecttp, self.Selectm, self.SelectMdyn, self.SelectChi2) = SelectOrbitsParams

        if SelectOrbitsEllipses is not None:
            (self.NbBodies, self.NbSelectOrbits, self.NbPtsEllipse, self.SelectP, self.Selectt, self.SelectRa, self.SelectDec, self.SelectZ, self.SelectSep, self.SelectPa) = SelectOrbitsEllipses

        if BestOrbitsParams is not None:
            (self.NbBodies, self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.BestMdyn, self.BestChi2) = BestOrbitsParams

        if BestOrbitsEllipses is not None:
            (self.NbBodies, self.NbPtsEllipse, self.BestP, self.Bestt, self.BestRa, self.BestDec, self.BestZ, self.BestSep, self.BestPa) = BestOrbitsEllipses

        self.InitParams()

        # Fix the width to avoid resizing of parameters
        left_width = self.WindowPlot.WidgetParam.sizeHint().width()
        self.WindowPlot.WidgetParam.setFixedWidth(left_width)

        # Widget container
        self.setLayout(Layout)  # GeneralToolClass is directly the widget container

    def InitParams(self):
        """Initialize parameters. This method should be overridden by subclasses."""
        return

    def Refresh_ActivePlots(self):
        """Refresh all active plots when the refresh button is clicked."""
        if self.WindowPlot.isVisible():
            self.Plot()

    def closeEvent(self, e):
        """Close program when the main window is closed."""
        app.closeAllWindows()

    def ResetParams(self):
        """Reset all widgets of the parameters window."""
        for i in reversed(range(2, self.WindowPlot.WidgetParam.Layout.count())):
            WidgetToRemove = self.WindowPlot.WidgetParam.Layout.itemAt(i).widget()
            self.WindowPlot.WidgetParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        self.InitParams()

    def LabelOf(self, var=str):
        """Return the label for a given variable."""
        labels = {
            'P': 'Period [yr]',
            'a': 'Semi-major axis [AU]',
            'e': 'Eccentricity',
            'i': 'Inclinaison [°]',
            'w': 'Argument of periastron [°]',
            'W': 'Longitude of ascending node [°]',
            'tp': 'Periastron time passage [MJD]',
            'm': 'Body mass [Mjup]',
            'Mdyn': 'Dynamical mass [Mjup]',
            'Chi2': 'Chi square'
        }
        return labels.get(var, 'Unknown variable')

    def Toggle_WindowPlot(self):
        """Open the plot window when the Plot button is clicked."""
        self.Plot()
        self.WindowPlot.show()
        self.BtnPlot.setEnabled(False)


class SpaceView(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Space view', 'Space view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

    def InitParams(self):
        """Initialize parameters for the SpaceView tool."""

        # Number of studied body
        self.ListBody = [str(k+1) for k in range(self.NbBodies)]
        self.ListBody.append('all')
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Number of shown orbits
        self.NbShownOrbits = 500
        self.NbShownOrbitsWidget = SpinBox('Number of orbits', 'Number of shown orbits', ParamDefault=self.NbShownOrbits, ParamMin=0, ParamMax=self.NbSelectOrbits)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbShownOrbitsWidget)

        # Type of view
        self.indexView = 0
        self.ViewWidget = ComboBox('View', 'Dimension', ['2D XY', '2D XZ', '3D'])
        self.ViewWidget.ComboParam.setCurrentIndex(self.indexView)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ViewWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)

        # Show observations points
        self.CheckObs = CheckBox('Observations', 'Show the observations points with its error bar')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckObs)

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = 'all' if self.nBodyWidget.ComboParam.currentText() == 'all' else int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.NbShownOrbits = self.NbShownOrbitsWidget.SpinParam.value()

    def Plot(self):
        """Plot the orbits based on the selected parameters."""

        # Clear axis
        self.WindowPlot.WidgetPlot.Canvas.fig.clear()

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        if self.indexView == 0:  # 2D (x,y)
            self._plot_2d(self.SelectRa, self.SelectDec, self.BestRa, self.BestDec, r'$\delta Ra$ [mas]', r'$\delta Dec$ [mas]')
        elif self.indexView == 1:  # 2D (x,z)
            self._plot_2d(self.SelectRa, self.SelectZ, self.BestRa, self.BestZ, r'$\delta Ra$ [mas]', 'Depth [mas]')
        elif self.indexView == 2:  # 3D
            self._plot_3d()

        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.fig.tight_layout()
        self.WindowPlot.WidgetPlot.Canvas.draw()

    def _plot_2d(self, SelectX, SelectY, BestX, BestY, xlabel, ylabel):
        """Helper function to plot 2D views."""
        self.Subplot2D = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111, aspect='equal')
        self.Subplot2D.plot(0, 0, marker='*', color='orange', markersize=10)

        if self.nBody == 'all':
            for k in range(self.NbBodies):
                if self.CheckBestFit.CheckParam.isChecked():
                    self.Subplot2D.plot(BestX[k], BestY[k], color='r')
                for n in range(self.NbShownOrbits):
                    self.Subplot2D.plot(SelectX[k][n], SelectY[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
        else:
            if self.CheckBestFit.CheckParam.isChecked():
                self.Subplot2D.plot(BestX[self.nBody], BestY[self.nBody], color='r')
            for n in range(self.NbShownOrbits):
                self.Subplot2D.plot(SelectX[self.nBody][n], SelectY[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)

        if self.CheckObs.CheckParam.isChecked():
            self.Subplot2D.errorbar(self.Ra, self.Dec, self.DRa, self.DDec, linestyle='', color='b')  # Observed data

        self.Subplot2D.set_xlabel(xlabel)
        self.Subplot2D.set_ylabel(ylabel)
        self.Subplot2D.invert_xaxis()
        self.Subplot2D.set_aspect('equal', adjustable='box')
        self.Subplot2D.set_title(' ')

    def _plot_3d(self):
        """Helper function to plot 3D views."""
        self.Subplot3D = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111, projection='3d', aspect='equal')
        self.Subplot3D.plot(0, 0, 0, marker='*', color='orange', markersize=10)

        if self.nBody == 'all':
            for k in range(self.NbBodies):
                if self.CheckBestFit.CheckParam.isChecked():
                    self.Subplot3D.plot(self.BestRa[k], self.BestDec[k], self.BestZ[k], color='r')
                for n in range(self.NbShownOrbits):
                    self.Subplot3D.plot(self.SelectRa[k][n], self.SelectDec[k][n], self.SelectZ[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
        else:
            if self.CheckBestFit.CheckParam.isChecked():
                self.Subplot3D.plot(self.BestRa[self.nBody], self.BestDec[self.nBody], self.BestZ[self.nBody], color='r')
            for n in range(self.NbShownOrbits):
                self.Subplot3D.plot(self.SelectRa[self.nBody][n], self.SelectDec[self.nBody][n], self.SelectZ[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)

        self.Subplot3D.set_xlabel(r'$\delta Ra$ [mas]')
        self.Subplot3D.set_ylabel(r'$\delta Dec$ [mas]')
        self.Subplot3D.set_zlabel('Depth [mas]')
        self.Subplot3D.invert_xaxis()
        self.Subplot3D.set_aspect('equal', adjustable='box')
        self.Subplot3D.set_title(' ')


class TempoView(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Temporal view', 'Temporal view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Initialize subplots
        self.Subplot1 = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(211)
        self.Subplot2 = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(212)

    def InitParams(self):
        """Initialize parameters for the TempoView tool."""

        # Orbit number
        self.ListBody = [str(k+1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Number of shown orbits
        self.NbShownOrbits = 500
        self.NbShownOrbitsWidget = SpinBox('Number of orbits', 'Number of shown orbits', ParamDefault=self.NbShownOrbits, ParamMin=0, ParamMax=self.NbSelectOrbits)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbShownOrbitsWidget)

        # Choice of coordinate
        self.CoordinateWidget = ComboBox('Choice of coordinate', 'Coordinates', ['dRa', 'dDec', 'Sep', 'Pa'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CoordinateWidget)

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.CoordinateIndex = self.CoordinateWidget.ComboParam.currentIndex()
        if self.CoordinateIndex == 0:
            self.Coordinate = r'$\delta Ra$'
        elif self.CoordinateIndex == 1:
            self.Coordinate = r'$\delta Dec$'
        elif self.CoordinateIndex == 2:
            self.Coordinate = 'Sep'
        elif self.CoordinateIndex == 3:
            self.Coordinate = 'Pa'
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.NbShownOrbits = self.NbShownOrbitsWidget.SpinParam.value()

    def Plot(self):
        """Plot the temporal view based on the selected parameters."""

        # Clear the plot
        self.Subplot1.cla()
        self.Subplot2.cla()

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        # Determine the coordinate to plot
        if self.CoordinateIndex == 0:
            YplotOutput = self.SelectRa 
            BestYplotOutput = self.BestRa
            YplotInput = self.Ra
            YplotInputErr = self.DRa
        elif self.CoordinateIndex == 1:
            YplotOutput = self.SelectDec
            BestYplotOutput = self.BestDec
            YplotInput = self.Dec
            YplotInputErr = self.DDec
        elif self.CoordinateIndex == 2:
            YplotOutput = self.SelectSep
            BestYplotOutput = self.BestSep
            YplotInput = self.Sep
            YplotInputErr = self.DSep
        elif self.CoordinateIndex == 3:
            YplotOutput = self.SelectPa
            BestYplotOutput = self.BestPa
            YplotInput = self.Pa
            YplotInputErr = self.DPa

        # Plot output data
        for n in range(self.NbShownOrbits):
            Selectt3P = np.concatenate((self.Selectt[self.nBody][n] - self.SelectP[self.nBody][n] * 365.25, self.Selectt[self.nBody][n], self.Selectt[self.nBody][n] + self.SelectP[self.nBody][n] * 365.25))
            YplotOutput3P = np.concatenate((YplotOutput[self.nBody][n], YplotOutput[self.nBody][n], YplotOutput[self.nBody][n]))
            self.Subplot1.plot(Selectt3P, YplotOutput3P, color=self.colorList[0], linestyle='-', linewidth=0.2, alpha=0.1)

        Bestt3P = np.concatenate((self.Bestt[self.nBody] - self.BestP[self.nBody] * 365.25, self.Bestt[self.nBody], self.Bestt[self.nBody] + self.BestP[self.nBody] * 365.25))
        BestYplotOutput3P = np.concatenate((BestYplotOutput[self.nBody], BestYplotOutput[self.nBody], BestYplotOutput[self.nBody]))
        self.Subplot1.plot(Bestt3P, BestYplotOutput3P, linestyle='-', linewidth=0.5, color='r')

        # Plot input data and residuals
        for k in range(self.NbInputData_RA):
            if self.I_RA[k] == self.nBody + 1:
                self.Subplot1.errorbar(self.MJD_RA[k], YplotInput[k], YplotInputErr[k], linestyle='', color='b')
                indext = np.argmin(np.abs(Bestt3P - self.MJD_RA[k]))  # index of time of output data closer than time of input data
                Res = BestYplotOutput3P[indext] - YplotInput[k]  # Residual
                self.Subplot2.errorbar(Bestt3P[indext], Res, YplotInputErr[k], color='b')

        self.Subplot2.hlines(0, np.min(Bestt3P), np.max(Bestt3P), color='red', linewidth=0.5)

        # Plot features
        if self.CoordinateIndex == 3:
            self.Subplot1.set_ylabel(self.Coordinate + ' [°]')
            self.Subplot2.set_ylabel(self.Coordinate + ' - Bestfit [°]')
        else:
            self.Subplot1.set_ylabel(self.Coordinate + ' [mas]')
            self.Subplot2.set_ylabel(self.Coordinate + ' - Bestfit [mas]')
        self.Subplot2.set_xlabel('Time [MJD]')
        L = np.max(self.MJD_RA) - np.min(self.MJD_RA)
        self.Subplot1.set_xlim(np.min(self.MJD_RA) - 0.1 * L, np.max(self.MJD_RA) + 0.1 * L)
        self.Subplot2.set_xlim(self.Subplot1.get_xlim())
        self.Subplot1.grid()
        self.Subplot2.grid()

        # Update canvas
        self.Subplot1.set_title(' ')
        self.WindowPlot.WidgetPlot.Canvas.fig.tight_layout()
        self.WindowPlot.WidgetPlot.Canvas.draw()


class Conv(GeneralToolClass):
    def __init__(self, OutputParams):
        super().__init__('Convergence', 'Convergence of the fit orbit parameters', None, OutputParams, None, None, None, None)

        # Plot initialization
        self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

    def InitParams(self):
        """Initialize parameters for the Convergence tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Orbital parameter', 'Orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn', 'Chi2'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()

    def Plot(self):
        """Plot the convergence of the fit orbit parameters."""

        # Clear the plot
        self.Subplot.cla()

        # Update parameters
        try:
            self.UpdateParams()
        except Exception as e:
            print(f'Wrong Parameters: {e}')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        # Plot with current parameters
        self.Steps = range(self.NbOrbits)
        self.EvalParamOrbit = eval(f'self.{self.ParamOrbit}')[self.nBody]

        self.Subplot.plot(self.Steps, self.EvalParamOrbit, marker=',', linestyle='')

        # Plot features
        self.Subplot.set_xlabel('Step')
        self.Subplot.set_ylabel(self.LabelOf(self.ParamOrbit))

        # Update canvas
        self.Subplot.set_title(' ')
        self.WindowPlot.WidgetPlot.Canvas.fig.tight_layout()
        self.WindowPlot.WidgetPlot.Canvas.draw()


class Hist(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitsParams):
        super().__init__('Histogram', "Histogram of an orbital parameter", None, OutputParams, None, None, BestOrbitsParams, None)

        # Plot initialization
        self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

    def InitParams(self):
        """Initialize parameters for the Histogram tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Orbital parameter', 'Orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn', 'Chi2'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)

        # Histogram binning
        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)

        # Show confidence interval
        self.CheckMedian = CheckBox('Median:', 'Show the median and the 1 sigma confidence interval')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckMedian)

        # Confidence interval bounds
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.EvalParamOrbit = eval('self.' + self.ParamOrbit)[self.nBody]
        self.leftWidget = DoubleSpinBox(None, 'Left bound of the selected histogram', np.min(self.EvalParamOrbit), np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))
        self.CheckMedian.Layout.addWidget(self.leftWidget)
        self.leftWidget.setEnabled(self.CheckMedian.CheckParam.isChecked())
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.leftWidget.setEnabled(state))

        self.Llbl = QLabel(' < ')
        self.CheckMedian.Layout.addWidget(self.Llbl)
        self.Llbl.setEnabled(False)
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.Llbl.setEnabled(state))

        self.IntConf = 68
        self.IntConfWidget = SpinBox(None, 'Acceptable level of confidence', self.IntConf, 0, 100, 1)
        self.CheckMedian.Layout.addWidget(self.IntConfWidget)
        self.IntConfWidget.setEnabled(self.CheckMedian.CheckParam.isChecked())
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.IntConfWidget.setEnabled(state))

        self.Rlbl = QLabel(' > ')
        self.CheckMedian.Layout.addWidget(self.Rlbl)
        self.Rlbl.setEnabled(False)
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.Rlbl.setEnabled(state))

        self.rightWidget = DoubleSpinBox(None, 'Right bound of the selected histogram', np.max(self.EvalParamOrbit), np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))
        self.CheckMedian.Layout.addWidget(self.rightWidget)
        self.rightWidget.setEnabled(self.CheckMedian.CheckParam.isChecked())
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.rightWidget.setEnabled(state))

        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.ChangeRightandLeftBound)

    def ChangeRightandLeftBound(self):
        """Update the bounds of the histogram when the orbital parameter changes."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.EvalParamOrbit = eval('self.' + self.ParamOrbit)[self.nBody]

        self.leftWidget.SpinParam.setRange(np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))
        self.leftWidget.SpinParam.setValue(np.min(self.EvalParamOrbit))

        self.rightWidget.SpinParam.setRange(np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))
        self.rightWidget.SpinParam.setValue(np.max(self.EvalParamOrbit))

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.EvalParamOrbit = eval('self.' + self.ParamOrbit)[self.nBody]
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.IntConf = self.IntConfWidget.SpinParam.value()
        self.rightBound = self.rightWidget.SpinParam.value()
        self.leftBound = self.leftWidget.SpinParam.value()

    def Plot(self):
        """Plot the histogram based on the selected parameters."""

        # Clear the plot
        self.Subplot.cla()

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        # Plot histogram
        self.Subplot.hist(self.EvalParamOrbit, self.NbBins)

        # Plot best fit
        if self.CheckBestFit.CheckParam.isChecked():
            BestParam = eval('self.Best' + self.ParamOrbit)[self.nBody]
            self.Subplot.axvline(BestParam, color='red')
            self.Subplot.text(BestParam, 0.5 * self.Subplot.get_ylim()[1], s='{}'.format(np.around(BestParam, 3)), color='r', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='r'), fontsize=9, horizontalalignment='center', verticalalignment='center', rotation=0)

        # Plot median and confidence interval
        if self.CheckMedian.CheckParam.isChecked():
            counts, bin_edges = np.histogram([p for p in self.EvalParamOrbit if self.leftBound <= p <= self.rightBound], self.NbBins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            cumsum = np.cumsum(counts)
            total = cumsum[-1]

            seuil_median = 0.50 * total
            seuil_inf = (0.50 - self.IntConf / 100 / 2) * total
            seuil_sup = (0.50 + self.IntConf / 100 / 2) * total

            mediane = bin_centers[np.searchsorted(cumsum, seuil_median)]
            lower_bound = bin_centers[np.searchsorted(cumsum, seuil_inf)]
            upper_bound = bin_centers[np.searchsorted(cumsum, seuil_sup)]

            SigmaMinus = np.abs(mediane - lower_bound)
            SigmaPlus = np.abs(mediane - upper_bound)

            self.Subplot.axvline(mediane, linestyle='-', color='black')
            self.Subplot.text(mediane, 0.83 * self.Subplot.get_ylim()[1], s='{}'.format(np.around(mediane, 3)), bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize=9, horizontalalignment='center', verticalalignment='bottom')

            self.Subplot.axvline(lower_bound, linestyle='--', color='black')
            self.Subplot.text(lower_bound, 0.8 * self.Subplot.get_ylim()[1], s='-{}'.format(np.around(SigmaMinus, 3)), bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize=8, horizontalalignment='right', verticalalignment='top')

            self.Subplot.axvline(upper_bound, linestyle='--', color='black')
            self.Subplot.text(upper_bound, 0.8 * self.Subplot.get_ylim()[1], s='+{}'.format(np.around(SigmaPlus, 3)), bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize=8, horizontalalignment='left', verticalalignment='top')

            self.Subplot.axvspan(np.min(self.EvalParamOrbit), self.leftBound, facecolor='grey', alpha=0.5)
            self.Subplot.axvspan(self.rightBound, np.max(self.EvalParamOrbit), facecolor='grey', alpha=0.5)

        # Plot features
        self.Subplot.set_xlabel(self.LabelOf(self.ParamOrbit))
        self.Subplot.set_ylabel('Count number')
        self.Subplot.set_xlim(np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))

        # Update canvas
        self.Subplot.set_title(' ')
        self.WindowPlot.WidgetPlot.Canvas.fig.tight_layout()
        self.WindowPlot.WidgetPlot.Canvas.draw()
        
    
class Hist2D(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitsParams):
        super().__init__('Histogram 2D', 'Histogram of an orbital parameter as function of another', None, OutputParams, None, None, BestOrbitsParams, None)

    def InitParams(self):
        """Initialize parameters for the 2D Histogram tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Abscissa orbit parameters
        self.XParamOrbitWidget = ComboBox('X Orbital parameter', 'Abscissa orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XParamOrbitWidget)

        # Ordinate orbit parameters
        self.YParamOrbitWidget = ComboBox('Y Orbital parameter', 'Ordinate orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YParamOrbitWidget)

        # Histogram binning
        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.XParamOrbit = self.XParamOrbitWidget.ComboParam.currentText()
        self.YParamOrbit = self.YParamOrbitWidget.ComboParam.currentText()
        self.NbBins = self.NbBinsWidget.SpinParam.value()

    def Plot(self):
        """Plot the 2D histogram based on the selected parameters."""

        # Clear axis
        self.WindowPlot.WidgetPlot.Canvas.fig.clear()

        # Plot initialization
        self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        # Plot with current parameters
        self.EvalXParamOrbit = eval(f'self.{self.XParamOrbit}')[self.nBody]
        self.EvalYParamOrbit = eval(f'self.{self.YParamOrbit}')[self.nBody]

        hist = self.Subplot.hist2d(self.EvalXParamOrbit, self.EvalYParamOrbit, (self.NbBins, self.NbBins))
        ColorbarAx = make_axes_locatable(self.Subplot).append_axes('right', size='5%', pad=0.1)
        self.WindowPlot.WidgetPlot.Canvas.fig.colorbar(hist[3], ColorbarAx, label='Count number')

        # Best fit
        if self.CheckBestFit.CheckParam.isChecked():
            BestXParam = eval(f'self.Best{self.XParamOrbit}')[self.nBody]
            BestYParam = eval(f'self.Best{self.YParamOrbit}')[self.nBody]
            self.Subplot.plot(BestXParam, BestYParam, color='red', marker='x')

        # Plot features
        self.Subplot.set_xlabel(self.LabelOf(self.XParamOrbit))
        self.Subplot.set_ylabel(self.LabelOf(self.YParamOrbit))

        # Update canvas
        self.Subplot.set_title(' ')
        self.WindowPlot.WidgetPlot.Canvas.fig.tight_layout()
        self.WindowPlot.WidgetPlot.Canvas.draw()


class Corner(GeneralToolClass):
    def __init__(self, SelectOrbitsParams, BestOrbitsParams):
        super().__init__('Corner', 'Corner plot of parameters', None, None, SelectOrbitsParams, None, BestOrbitsParams, None)

    def InitParams(self):
        """Initialize parameters for the Corner plot tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Parameters to include in the corner plot
        self.ParamContainer = QWidget()
        self.ParamLayoutV = QVBoxLayout()

        self.ParamLayoutV.addWidget(QLabel('Orbit parameters:'))

        self.ParamCheckGroupBox = QGroupBox()
        self.ParamCheckLayout = QGridLayout()

        self.OrbitParams = [
            ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn'],
            ['Period', 'Semi-major axis', 'Eccentricity', 'Inclination', 'Argument of periastron', 'Longitude of ascending node', 'Periastron time passage', 'Body mass', 'Dynamical mass']
        ]
        self.WidgetOrbitParams = []

        for i in range(len(self.OrbitParams[0])):
            ParamCheckBox = CheckBox(self.OrbitParams[0][i], self.OrbitParams[1][i])

            # Check the first three parameters by default
            if i in (0,1,2):
                ParamCheckBox.CheckParam.setChecked(True)
            else:
                ParamCheckBox.CheckParam.setChecked(False)

            self.WidgetOrbitParams.append(ParamCheckBox)
            self.ParamCheckLayout.addWidget(ParamCheckBox, i // 3, i % 3)

        self.ParamCheckGroupBox.setLayout(self.ParamCheckLayout)

        self.ParamLayoutV.addWidget(self.ParamCheckGroupBox)
        self.ParamContainer.setLayout(self.ParamLayoutV)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamContainer)

        # Check if the parameter has variance
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.CheckParamsVar)
        self.CheckParamsVar()

        # Histogram binning
        self.NbBins = 20
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)


    def CheckParamsVar(self):
        """Check if the parameter has variance."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        for x in self.WidgetOrbitParams:
            param_data = eval(f'self.Select{x.CheckParam.text()}')[self.nBody]
            if np.var(param_data) > 0:  # Check if the parameter has variance
                x.CheckParam.setEnabled(True)
            else:
                x.CheckParam.setEnabled(False)


    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.NbBins = self.NbBinsWidget.SpinParam.value()

    def Plot(self):
        """Plot the corner plot based on the selected parameters."""

        # Clear axis
        self.WindowPlot.WidgetPlot.Canvas.fig.clear()

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        Data = []
        DataLabels = []
        for x in self.WidgetOrbitParams:
            if x.CheckParam.isChecked():
                Data.append(eval(f'self.Select{x.CheckParam.text()}')[self.nBody])
                DataLabels.append(str(x.CheckParam.text()))

        Data = np.array(Data).T
        if len(Data) == 0:
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        grid = corner.corner(Data, labels=DataLabels, bins=self.NbBins, fig=self.WindowPlot.WidgetPlot.Canvas.fig)

        # Adjust the labels to not be slanted
        for k in range(len(grid.get_axes())):
            ax = grid.get_axes()[k]
            ax.tick_params(axis='x', rotation=0)
            ax.tick_params(axis='y', rotation=0)

            # Best fit
            if self.CheckBestFit.CheckParam.isChecked():

                if k % (len(DataLabels) + 1) == 0:
                    param_index = k // (len(DataLabels) + 1)
                    BestParam = eval(f'self.Best{DataLabels[param_index]}')[self.nBody]
                    ax.axvline(BestParam, color='red', linestyle='-')
                else:
                    row = k // len(DataLabels)
                    col = k % len(DataLabels)
                    if row > col:
                        BestXParam = eval(f'self.Best{DataLabels[col]}')[self.nBody]
                        BestYParam = eval(f'self.Best{DataLabels[row]}')[self.nBody]
                        ax.plot(BestXParam, BestYParam, color='red', marker='x')


            # if k in (1,3,6,)

            # BestParam = eval('self.Best' + DataLabels[k])[self.nBody]
            # ax.axvline(BestParam, color='red')


        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.fig.tight_layout()
        self.WindowPlot.WidgetPlot.Canvas.draw()




class PosAtDate(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Position at date', 'Position of bodies at a given date', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

    def InitParams(self):
        """Initialize parameters for the Position at Date tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Date of wanted observation
        self.DateWidget = DateAndMJDEdit('Date', 'Date of wanted observation')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.DateWidget)

        # Histogram binning
        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)

        # Show observations points
        self.CheckObs = CheckBox('Observations', 'Show the observations points with its error bar')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckObs)

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.Date = self.DateWidget.MJDWidget.SpinParam.value()

    def Plot(self):
        """Plot the position at the given date based on the selected parameters."""

        # Clear axis
        self.WindowPlot.WidgetPlot.Canvas.fig.clear()

        # Plot initialization
        self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        # Compute date
        SelectXAtDate, SelectYAtDate = [np.zeros(self.NbSelectOrbits) for _ in range(2)]

        for k in range(self.NbSelectOrbits):
            SelectPeriod = np.max(self.Selectt[self.nBody][k]) - np.min(self.Selectt[self.nBody][k])
            SelectDate = self.Date

            while SelectDate < np.min(self.Selectt[self.nBody][k]):
                SelectDate += SelectPeriod

            while SelectDate > np.max(self.Selectt[self.nBody][k]):
                SelectDate -= SelectPeriod

            indexBestDate = np.argmin(np.abs(self.Selectt[self.nBody][k] - SelectDate))
            SelectXAtDate[k] = self.SelectRa[self.nBody][k][indexBestDate]
            SelectYAtDate[k] = self.SelectDec[self.nBody][k][indexBestDate]

        # Plot with current parameters
        hist = self.Subplot.hist2d(SelectXAtDate, SelectYAtDate, bins=(self.NbBins, self.NbBins), range=((np.min(self.SelectRa), np.max(self.SelectRa)), (np.min(self.SelectDec), np.max(self.SelectDec))))
        ColorbarAx = make_axes_locatable(self.Subplot).append_axes('right', size='5%', pad=0.1)
        self.WindowPlot.WidgetPlot.Canvas.fig.colorbar(hist[3], ColorbarAx, ticks=[], label='Probability')

        self.Subplot.plot(0, 0, marker='*', color='orange', markersize=10)

        if self.CheckBestFit.CheckParam.isChecked():
            self.Subplot.plot(self.BestRa[self.nBody], self.BestDec[self.nBody], color='r', lw=0.5)

            BestPeriod = np.max(self.Bestt[self.nBody]) - np.min(self.Bestt[self.nBody])
            BestDate = self.Date

            while BestDate < np.min(self.Bestt[self.nBody]):
                BestDate += BestPeriod

            while BestDate > np.max(self.Bestt[self.nBody]):
                BestDate -= BestPeriod

            indexBestDate = np.argmin(np.abs(self.Bestt[self.nBody] - BestDate))
            self.Subplot.plot(self.BestRa[self.nBody][indexBestDate], self.BestDec[self.nBody][indexBestDate], marker='x', color='red')

        if self.CheckObs.CheckParam.isChecked():
            self.Subplot.errorbar(self.Ra, self.Dec, self.DRa, self.DDec, linestyle='')  # Observed data

        # Plot features
        self.Subplot.set_xlabel(r'$\delta$ Ra')
        self.Subplot.set_ylabel(r'$\delta$ Dec')
        self.Subplot.invert_xaxis()
        self.Subplot.set_aspect('equal', adjustable='box')

        # Update canvas
        self.Subplot.set_title(' ')
        self.WindowPlot.WidgetPlot.Canvas.fig.tight_layout()
        self.WindowPlot.WidgetPlot.Canvas.draw()

