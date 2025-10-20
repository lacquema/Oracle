#! /Users/lacquema/ByeGildas/bin/python3

### --- Packages --- ###

# Transverse packages
import sys
import numpy as np
from random import random, randint
import corner
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re

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

        self.colorList = ['darkslategrey', 'darkolivegreen', 'indianred', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

        # Layout
        Layout = QHBoxLayout()

        # Plot button
        self.BtnPlot = QPushButton(ToolName)
        self.BtnPlot.clicked.connect(self.Toggle_WindowPlot)
        self.BtnPlot.setStatusTip(ToolStatus)
        Layout.addWidget(self.BtnPlot)

        # Initialisation of plot windows
        self.WindowPlot = WindowPlot(ToolName)
        self.WindowPlot.SignalCloseWindowPlot.connect(lambda: self.BtnPlot.setEnabled(True))  # Enable button when window is closed

        # Connections between parameters and plot
        self.WindowPlot.WidgetParam.BtnReset.clicked.connect(self.ResetParams)
        self.WindowPlot.WidgetParam.BtnRefresh.clicked.connect(self.refresh_plots)

        # Initialisation of Data
        self.InputData = InputData

        if OutputParams is not None:
            (self.NbBodies, self.NbOrbits, self.P, self.a, self.e, self.i, self.w, self.W, self.tp, self.m, self.m0, self.Chi2, self.map) = OutputParams

        if SelectOrbitsParams is not None:
            (self.NbBodies, self.NbSelectOrbits, self.SelectP, self.Selecta, self.Selecte, self.Selecti, self.Selectw, self.SelectW, self.Selecttp, self.Selectm, self.Selectm0, self.SelectChi2) = SelectOrbitsParams

        if SelectOrbitsEllipses is not None:
            (self.NbBodies, self.NbSelectOrbits, self.NbPtsEllipse, self.SelectP, self.Selectt, self.SelectRa, self.SelectDec, self.SelectZ, self.SelectSep, self.SelectPa, self.SelectRV) = SelectOrbitsEllipses

        if BestOrbitsParams is not None:
            (self.NbBodies, self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.Bestm0, self.BestChi2) = BestOrbitsParams

        if BestOrbitsEllipses is not None:
            (self.NbBodies, self.NbPtsEllipse, self.BestP, self.Bestt, self.BestRa, self.BestDec, self.BestZ, self.BestSep, self.BestPa, self.BestRV) = BestOrbitsEllipses

        # # Fix the width to avoid resizing of parameters
        # left_width = self.WindowPlot.WidgetParam.sizeHint().width()
        # self.WindowPlot.WidgetParam.setFixedWidth(left_width)

        # Widget container
        self.setLayout(Layout)  # GeneralToolClass is directly the widget container

    def InitParams(self):
        """Initialize parameters. This method should be overridden by subclasses."""
        return

    def refresh_plots(self):
        """Refresh all active plots when the refresh button is clicked."""
        for WidgetPlot in self.WindowPlot.WidgetPlots:
            if WidgetPlot.isVisible():
                WidgetPlot.refresh_plot()

    def reset_plots(self):
        """Reset all active plots when the reset button is clicked."""
        for WidgetPlot in self.WindowPlot.WidgetPlots:
            WidgetPlot.reset_plot()
        self.refresh_plots()
    
    def Toggle_WindowPlot(self):
        """Open the plot window when the Plot button is clicked."""
        self.WindowPlot.show()
        self.BtnPlot.setEnabled(False)
        self.refresh_plots()
        
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
        self.reset_plots()

    def LabelOf(self, var=str):
        """Return the label for a given variable."""
        labels = {
            'P': 'Period',
            'a': 'Semi-major axis',
            'e': 'Eccentricity',
            'i': 'Inclinaison',
            'w': 'Argument of periastron',
            'W': 'Longitude of ascending node',
            'tp': 'Periastron time passage',
            'm': 'Body mass',
            'm0': 'Central body mass',
            'Chi2': 'Chi square'
        }
        return labels.get(var, 'Unknown variable')
    
    def UnitOf(self, var=str):
        """Return the unit for a given variable."""
        units = {
            'P': '[yr]',
            'a': '[AU]',
            'e': '',
            'i': '[°]',
            'w': '[°]',
            'W': '[°]',
            'tp': '[MJD]',
            'm': '[Mjup]',
            'm0': '[Mjup]',
            'Chi2': ''
        }
        return units.get(var, 'Unknown variable')
        
    def replace_params_in_formula(self, formula, prefixe, nOrbitDefault):
        """Replace parameters and functions in the formula with their corresponding values."""
        # print(formula)
        for num in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
            formula = formula.replace(f'[{num}]', f'[{str(int(num)-1)}]') # Replace [n] by [n-1]
        for param in ['Chi2', 'P', 'a', 'e', 'tp', 'm', 'm0']:
            formula = re.sub(r'\b' + param + r'\b(?!\[)', f'{param}[{nOrbitDefault}]', formula) # Add [nOrbitDefault] to the parameter
            formula = re.sub(r'\b' + param + r'\b', f'{prefixe}{param}', formula) # Replace the parameter by its value
        for param in ['i', 'w', 'W']:
            formula = re.sub(r'\b' + param + r'\b(?!\[)', f'{param}[{nOrbitDefault}]', formula) # add [nOrbitDefault]
            formula = re.sub(rf'\b{param}\[(\d+)\]', rf'np.radians({prefixe}{param}[\1])', formula)
        for fonction in ['sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2', 'hypot', 'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh', 'exp', 'expm1', 'exp2', 'log', 'log10', 'log2', 'log1p', 'sqrt', 'square', 'cbrt', 'power', 'erf', 'erfc', 'gamma', 'lgamma', 'digamma', 'beta']:
            formula = re.sub(r'\b' + fonction + r'\b', f'np.{fonction}', formula) # Replace the function by its numpy equivalent
        # Convert the result to degrees if it is an angle in radians
        if any(f'np.{angle_func}' in formula for angle_func in ['arcsin', 'arccos', 'arctan', 'arctan2']):
            formula = f'np.degrees({formula})'
        # print('Formula is: '+formula)
        return formula
    
    def evaluate_formula(self, formula, prefix, nOrbitDefault):
        """Evaluate a formula."""
        print(formula)
        if not formula:
            return None
        formula = self.replace_params_in_formula(formula, prefix, nOrbitDefault)
        print(formula)
        try:
            return eval(formula)
        except Exception as e:
            print(f'Error evaluating formula: {e}')
            return None


    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        # This method should be overridden by subclasses to update specific parameters
        pass

    def try_UpdateParams(self, WidgetPlot):
        """Try to update parameters and handle exceptions."""
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            WidgetPlot.Canvas.fig.draw()
    

class SpaceView(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Space view', 'Space view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Window plots initialisation
        self.WidgetPlotXY = self.WindowPlot.add_WidgetPlot(self.PlotXY, xlim=True, ylim=True)
        self.WidgetPlotXZ = self.WindowPlot.add_WidgetPlot(self.PlotXZ, xlim=True, ylim=True)
        self.WidgetPlotXYZ = self.WindowPlot.add_WidgetPlot(self.PlotXYZ, xlim=True, ylim=True, zlim=True, azim=True, elev=True)

        # Parameters initialisation
        self.InitParams()
    
    def InitParams(self):
        """Initialize parameters for the SpaceView tool."""

        # Number of studied body
        self.ListBody = [str(k+1) for k in range(self.NbBodies)]
        if self.NbBodies > 1:
            self.ListBody.append('all')
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots) 

        # Number of shown orbits
        self.NbShownOrbits = 500
        self.NbShownOrbitsWidget = SpinBox('Number of orbits', 'Number of shown orbits', ParamDefault=self.NbShownOrbits, ParamMin=0, ParamMax=self.NbSelectOrbits)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbShownOrbitsWidget)

        # Type of view
        self.ViewWidget = ComboBox('View', 'Dimension', ['2D XY', '2D XZ', '3D'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ViewWidget)
        self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.ViewWidget.ComboParam.currentIndexChanged.connect(self.indexViewChanged)
        self.indexViewChanged(self.indexView)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)

        # Show observations points
        self.CheckObs = CheckBox('Observations', 'Show the observations points with its error bar')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckObs)
        if self.InputData is None:
            self.CheckObs.CheckParam.setEnabled(False)
        self.ViewWidget.ComboParam.currentIndexChanged.connect(lambda: self.CheckObs.CheckParam.setEnabled(self.ViewWidget.ComboParam.currentIndex() == 0))  # Enable observations only for 2D view

        # Show date of observations
        self.CheckDateObs = CheckBox('Date of observations', 'Show the date of observations')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckDateObs)
        self.CheckDateObs.CheckParam.setEnabled(False)  # Initially disabled, enabled in Plot method if observations are shown
        self.ViewWidget.ComboParam.currentIndexChanged.connect(lambda: self.CheckDateObs.CheckParam.setEnabled(self.ViewWidget.ComboParam.currentIndex() == 0 and self.CheckObs.CheckParam.isChecked()))
        self.CheckObs.CheckParam.stateChanged.connect(lambda bool: self.CheckDateObs.CheckParam.setEnabled(bool))

    def indexViewChanged(self, value):
        self.indexView = value
        if self.indexView == 0:
            self.WidgetPlotXY.setVisible(True)
            self.WidgetPlotXZ.setVisible(False)
            self.WidgetPlotXYZ.setVisible(False)
        elif self.indexView == 1:
            self.WidgetPlotXY.setVisible(False)
            self.WidgetPlotXZ.setVisible(True)
            self.WidgetPlotXYZ.setVisible(False)
        elif self.indexView == 2:
            self.WidgetPlotXY.setVisible(False)
            self.WidgetPlotXZ.setVisible(False)
            self.WidgetPlotXYZ.setVisible(True)
        self.refresh_plots()

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = 'all' if self.nBodyWidget.ComboParam.currentText() == 'all' else int(self.nBodyWidget.ComboParam.currentText()) - 1
        # self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.NbShownOrbits = self.NbShownOrbitsWidget.SpinParam.value()


    def PlotXY(self):
        """Plot the 2D view of the orbits in the XY plane."""

        # Update parameters
        self.try_UpdateParams(self.WidgetPlotXY)

        # Add subplot
        self.SubplotXY = self.WidgetPlotXY.Canvas.fig.add_subplot(111, aspect='equal')

        # Central star
        self.SubplotXY.plot(0, 0, marker='*', color='orange', markersize=10)

        # X, Y current limits
        xlim_init = (None, None)
        ylim_init = (None, None)
        (Xmin, Xmax) = self.WidgetPlotXY.history[self.WidgetPlotXY.history_index]['xlim'] if len(self.WidgetPlotXY.history)!=0 else xlim_init
        (Ymin, Ymax) = self.WidgetPlotXY.history[self.WidgetPlotXY.history_index]['ylim'] if len(self.WidgetPlotXY.history)!=0 else ylim_init

        if self.nBody == 'all':
            for k in range(self.NbBodies):
                for n in range(self.NbShownOrbits):
                    self.SubplotXY.plot(self.SelectRa[k][n], self.SelectDec[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
                if self.CheckBestFit.CheckParam.isChecked():
                    self.SubplotXY.plot(self.BestRa[k], self.BestDec[k], color='C3', linewidth=0.5)
        else:
            for n in range(self.NbShownOrbits):
                self.SubplotXY.plot(self.SelectRa[self.nBody][n], self.SelectDec[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)
            if self.CheckBestFit.CheckParam.isChecked():
                self.SubplotXY.plot(self.BestRa[self.nBody], self.BestDec[self.nBody], color='C3', linewidth=0.5)

        # Add observations points if available
        if self.CheckObs.CheckParam.isChecked() and self.indexView == 0:
            if self.nBody == 'all':
                for k in range(self.InputData['Planets']['Nb']):
                    ra = self.InputData['Planets']['DataAstrom']['Ra'][k]
                    dec = self.InputData['Planets']['DataAstrom']['Dec'][k]
                    dra = self.InputData['Planets']['DataAstrom']['dRa'][k]
                    ddec = self.InputData['Planets']['DataAstrom']['dDec'][k]
                    dates = self.InputData['Planets']['DataAstrom']['Date'][k]
                    self.SubplotXY.errorbar(ra, dec, ddec, dra, linestyle='', color='black')
                    if self.CheckDateObs.CheckParam.isChecked():
                        if Xmin!=None and Ymin!=None: 
                            self.annotate_dates(dates, ra, dec, Xmin, Xmax, Ymin, Ymax)

            else:
                ra = self.InputData['Planets']['DataAstrom']['Ra'][self.nBody]
                dec = self.InputData['Planets']['DataAstrom']['Dec'][self.nBody]
                dra = self.InputData['Planets']['DataAstrom']['dRa'][self.nBody]
                ddec = self.InputData['Planets']['DataAstrom']['dDec'][self.nBody]
                dates = self.InputData['Planets']['DataAstrom']['Date'][self.nBody]
                self.SubplotXY.errorbar(ra, dec, ddec, dra, linestyle='', color='black')
                if self.CheckDateObs.CheckParam.isChecked():
                    if Xmin!=None and Ymin!=None: 
                        self.annotate_dates(dates, ra, dec, Xmin, Xmax, Ymin, Ymax)
        
        # Set axis
        self.SubplotXY.set_xlabel(r'$\delta Ra$ [mas]')
        self.SubplotXY.set_ylabel(r'$\delta Dec$ [mas]')
        self.SubplotXY.invert_xaxis()
        self.SubplotXY.set_aspect('equal', adjustable='box')


    def PlotXZ(self):
        """Plot the 2D view of the orbits in the XZ plane."""

        # Update parameters
        self.try_UpdateParams(self.WidgetPlotXZ)

        # Add subplot
        self.SubplotXZ = self.WidgetPlotXZ.Canvas.fig.add_subplot(111, aspect='equal')

        # Central star
        self.SubplotXZ.plot(0, 0, marker='*', color='orange', markersize=10)

        if self.nBody == 'all':
            for k in range(self.NbBodies):
                for n in range(self.NbShownOrbits):
                    self.SubplotXZ.plot(self.SelectRa[k][n], self.SelectZ[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
                if self.CheckBestFit.CheckParam.isChecked():
                    self.SubplotXZ.plot(self.BestRa[k], self.BestZ[k], color='C3', linewidth=0.5)
        else:
            for n in range(self.NbShownOrbits):
                self.SubplotXZ.plot(self.SelectRa[self.nBody][n], self.SelectZ[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)
            if self.CheckBestFit.CheckParam.isChecked():
                self.SubplotXZ.plot(self.BestRa[self.nBody], self.BestZ[self.nBody], color='C3', linewidth=0.5)

        # Set axis
        self.SubplotXZ.set_xlabel(r'$\delta Ra$ [mas]')
        self.SubplotXZ.set_ylabel('Depth [mas]')
        self.SubplotXZ.invert_xaxis()
        self.SubplotXZ.set_aspect('equal', adjustable='box')

    def PlotXYZ(self):
        """Helper function to plot 3D views."""

        # Update parameters
        self.try_UpdateParams(self.WidgetPlotXYZ)

        # Add subplot
        self.SubplotXYZ = self.WidgetPlotXYZ.Canvas.fig.add_subplot(111, aspect='equal', projection='3d')

        # Central star
        self.SubplotXYZ.plot(0, 0, 0, marker='*', color='orange', markersize=10)

        if self.nBody == 'all':
            for k in range(self.NbBodies):
                if self.CheckBestFit.CheckParam.isChecked():
                    self.SubplotXYZ.plot(self.BestRa[k], self.BestDec[k], self.BestZ[k], color='r', linewidth=0.5)
                for n in range(self.NbShownOrbits):
                    self.SubplotXYZ.plot(self.SelectRa[k][n], self.SelectDec[k][n], self.SelectZ[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
        else:
            if self.CheckBestFit.CheckParam.isChecked():
                self.SubplotXYZ.plot(self.BestRa[self.nBody], self.BestDec[self.nBody], self.BestZ[self.nBody], color='r', linewidth=0.5)
            for n in range(self.NbShownOrbits):
                self.SubplotXYZ.plot(self.SelectRa[self.nBody][n], self.SelectDec[self.nBody][n], self.SelectZ[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)

        self.SubplotXYZ.set_xlabel(r'$\delta Ra$ [mas]')
        self.SubplotXYZ.set_ylabel(r'$\delta Dec$ [mas]')
        self.SubplotXYZ.set_zlabel('Depth [mas]')
        self.SubplotXYZ.invert_xaxis()
        self.SubplotXYZ.set_aspect('equal', adjustable='box')


    def annotate_dates(self, dates, ra, dec, Xmin, Xmax, Ymin, Ymax):
        """Annotate dates on the plot."""
        ra_range = abs(Xmax - Xmin)
        dec_range = abs(Ymax - Ymin)
        min_dist_x = 0.1 * ra_range
        min_dist_y = 0.05 * dec_range
        min_dist = np.hypot(min_dist_x, min_dist_y)
        annotated_points = []
        for idx in range(len(dates)):
            y = dec[idx] + 0.01 * dec_range
            # Try right, if too close to previous, try left
            x_right = ra[idx] - 0.01 * ra_range 
            overlap = any(np.hypot(x_right - px, y - py) < min_dist for px, py in annotated_points)
            if not overlap:
                self.SubplotXY.annotate(f"{dates[idx]:.0f}", (x_right, y), color='black', fontsize=8, ha='left')
                annotated_points.append((x_right, y))
            else:
                x_left = ra[idx] + 0.01 * ra_range
                overlap = any(np.hypot(x_left - px, y - py) < min_dist for px, py in annotated_points)
                if not overlap:
                    self.SubplotXY.annotate(f"{dates[idx]:.0f}", (x_left, y), color='black', fontsize=8, ha='right')
                    annotated_points.append((x_left, y))
            

class TempoView(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Temporal view', 'Temporal view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Window plots initialisation
        self.WidgetPlot1 = self.WindowPlot.add_WidgetPlot(self.Plot1, xlim=True, ylim=True)
        self.WidgetPlot2 = self.WindowPlot.add_WidgetPlot(self.Plot2, layout=self.WidgetPlot1.Layout)

        # Parameters initialisation
        self.InitParams()

    def InitParams(self):
        """Initialize parameters for the TempoView tool."""

        # Orbit number
        self.ListBody = [str(k+1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)

        # Number of shown orbits
        self.NbShownOrbits = 500
        self.NbShownOrbitsWidget = SpinBox('Number of orbits', 'Number of shown orbits', ParamDefault=self.NbShownOrbits, ParamMin=0, ParamMax=self.NbSelectOrbits)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbShownOrbitsWidget)

        # Choice of coordinate
        self.CoordinateWidget = ComboBox('Choice of coordinate', 'Coordinates', ['dRa', 'dDec', 'Sep', 'Pa', 'RV'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CoordinateWidget)
        self.CoordinateWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)


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
        elif self.CoordinateIndex == 4:
            self.Coordinate = 'RV'
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.NbShownOrbits = self.NbShownOrbitsWidget.SpinParam.value()

    def general_plot(self):
        """Plot the temporal view based on the selected parameters."""
        self.YplotInput = []
        self.YplotInputErr = []
        # Determine the coordinate to plot
        if self.CoordinateIndex == 0:
            self.YplotOutput = self.SelectRa 
            self.BestYplotOutput = self.BestRa
            if self.InputData is not None: 
                self.YplotInput = self.InputData['Planets']['DataAstrom']['Ra'][self.nBody]
                self.YplotInputErr = self.InputData['Planets']['DataAstrom']['dRa'][self.nBody]
        elif self.CoordinateIndex == 1:
            self.YplotOutput = self.SelectDec
            self.BestYplotOutput = self.BestDec
            if self.InputData is not None: 
                self.YplotInput = self.InputData['Planets']['DataAstrom']['Dec'][self.nBody]
                self.YplotInputErr = self.InputData['Planets']['DataAstrom']['dDec'][self.nBody]
        elif self.CoordinateIndex == 2:
            self.YplotOutput = self.SelectSep
            self.BestYplotOutput = self.BestSep
            if self.InputData is not None: 
                self.YplotInput = self.InputData['Planets']['DataAstrom']['Sep'][self.nBody]
                self.YplotInputErr = self.InputData['Planets']['DataAstrom']['dSep'][self.nBody]
        elif self.CoordinateIndex == 3:
            self.YplotOutput = self.SelectPa
            self.BestYplotOutput = self.BestPa
            if self.InputData is not None: 
                self.YplotInput = self.InputData['Planets']['DataAstrom']['Pa'][self.nBody]
                self.YplotInputErr = self.InputData['Planets']['DataAstrom']['dPa'][self.nBody]
        elif self.CoordinateIndex == 4:
            self.YplotOutput = self.SelectRV
            self.BestYplotOutput = self.BestRV
            if self.InputData is not None:
                if self.InputData['Planets']['NbDataRV'][self.nBody] != 0:
                    self.YplotInput = self.InputData['Planets']['DataRV']['RV'][self.nBody]
                    self.YplotInputErr = self.InputData['Planets']['DataRV']['dRV'][self.nBody]

        # 3 periods of best fit
        self.Bestt3P = np.concatenate((self.Bestt[self.nBody] - self.BestP[self.nBody] * 365.25, self.Bestt[self.nBody], self.Bestt[self.nBody] + self.BestP[self.nBody] * 365.25))
        self.BestYplotOutput3P = np.concatenate((self.BestYplotOutput[self.nBody], self.BestYplotOutput[self.nBody], self.BestYplotOutput[self.nBody]))

        # Range of dates
        if len(self.YplotInput) != 0: 
            if self.CoordinateIndex != 4:
                self.DateRange = np.max(self.InputData['Planets']['DataAstrom']['Date'][self.nBody]) - np.min(self.InputData['Planets']['DataAstrom']['Date'][self.nBody])
            else:
                self.DateRange = np.max(self.InputData['Planets']['DataRV']['Date'][self.nBody]) - np.min(self.InputData['Planets']['DataRV']['Date'][self.nBody])

    def Plot1(self):

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot1)

        # Add subplot
        self.Subplot1 = self.WidgetPlot1.Canvas.fig.add_subplot(111, label='Main plot')

        # General plot
        self.general_plot()

        # X current limits
        if len(self.YplotInput) != 0: 
            if self.CoordinateIndex != 4:
                xlim_init = (np.min(self.InputData['Planets']['DataAstrom']['Date'][self.nBody]) - 0.1 * self.DateRange, np.max(self.InputData['Planets']['DataAstrom']['Date'][self.nBody]) + 0.1 * self.DateRange)
            else:
                xlim_init = (np.min(self.InputData['Planets']['DataRV']['Date'][self.nBody]) - 0.1 * self.DateRange, np.max(self.InputData['Planets']['DataRV']['Date'][self.nBody]) + 0.1 * self.DateRange)
        else:
            xlim_init = (None, None)
        xlim = self.WidgetPlot1.history[self.WidgetPlot1.history_index]['xlim'] if len(self.WidgetPlot1.history)!=0 else xlim_init

        # Plot output data
        for n in range(self.NbShownOrbits):
            Selectt3P = np.concatenate((self.Selectt[self.nBody][n] - self.SelectP[self.nBody][n] * 365.25, self.Selectt[self.nBody][n], self.Selectt[self.nBody][n] + self.SelectP[self.nBody][n] * 365.25))
            YplotOutput3P = np.concatenate((self.YplotOutput[self.nBody][n], self.YplotOutput[self.nBody][n], self.YplotOutput[self.nBody][n]))
            self.Subplot1.plot(Selectt3P, YplotOutput3P, color=self.colorList[0], linestyle='-', linewidth=0.2, alpha=0.1)
        
        if len(self.YplotInput) != 0: 
            if self.CoordinateIndex != 4:
                for k in range(self.InputData['Planets']['NbDataAstrom'][self.nBody]):
                    self.Subplot1.errorbar(self.InputData['Planets']['DataAstrom']['Date'][self.nBody][k], self.YplotInput[k], self.YplotInputErr[k], linestyle='', color='b')
            else:
                for k in range(self.InputData['Planets']['NbDataRV'][self.nBody]):
                    self.Subplot1.errorbar(self.InputData['Planets']['DataRV']['Date'][self.nBody][k], self.YplotInput[k], self.YplotInputErr[k], linestyle='', color='b')

        self.Subplot1.plot(self.Bestt3P, self.BestYplotOutput3P, linestyle='-', linewidth=0.5, color='r')

        # Plot features
        if self.CoordinateIndex == 3:
            self.Subplot1.set_ylabel(self.Coordinate + ' [°]')
        elif self.CoordinateIndex == 4:
            self.Subplot1.set_ylabel(self.Coordinate + ' [km/s]')
        else:
            self.Subplot1.set_ylabel(self.Coordinate + ' [mas]')
        self.Subplot1.grid()
        self.Subplot1.set_xlim(xlim_init)


    def Plot2(self):

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot2)

        # Add subplot
        self.Subplot2 = self.WidgetPlot2.Canvas.fig.add_subplot(111, label='Main plot')

        # General plot
        self.general_plot()

        # Plot input data and residuals
        if len(self.YplotInput) != 0: 
            if self.CoordinateIndex != 4:
                for k in range(self.InputData['Planets']['NbDataAstrom'][self.nBody]):
                    indext = np.argmin(np.abs(self.Bestt3P - self.InputData['Planets']['DataAstrom']['Date'][self.nBody][k]))  # index of time of output data closer than time of input data
                    Res = self.BestYplotOutput3P[indext] - self.YplotInput[k]  # Residual
                    self.Subplot2.errorbar(self.Bestt3P[indext], Res, self.YplotInputErr[k], color='b')
            else:
                for k in range(self.InputData['Planets']['NbDataRV'][self.nBody]):
                    indext = np.argmin(np.abs(self.Bestt3P - self.InputData['Planets']['DataRV']['Date'][self.nBody][k]))  # index of time of output data closer than time of input data
                    Res = self.BestYplotOutput3P[indext] - self.YplotInput[k]  # Residual
                    self.Subplot2.errorbar(self.Bestt3P[indext], Res, self.YplotInputErr[k], color='b')

        self.Subplot2.hlines(0, np.min(self.Bestt3P), np.max(self.Bestt3P), color='red', linewidth=0.5)

        # Plot features
        if self.CoordinateIndex == 3:
            self.Subplot2.set_ylabel(self.Coordinate + ' - Bestfit [°]')
        elif self.CoordinateIndex == 4:
            self.Subplot2.set_ylabel(self.Coordinate + ' - Bestfit [km/s]')
        else:
            self.Subplot2.set_ylabel(self.Coordinate + ' - Bestfit [mas]')
        self.Subplot2.set_xlabel('Time [MJD]')
        self.Subplot2.set_xlim(self.Subplot1.get_xlim())
        self.Subplot2.grid()


class Conv(GeneralToolClass):
    def __init__(self, OutputParams):
        super().__init__('Convergence', 'Convergence of the fit orbit parameters', None, OutputParams, None, None, None, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylim=True)

    def InitParams(self):
        """Initialize parameters for the Convergence tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Variable', 'Orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0', 'Chi2'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)


    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()

    def Plot(self):
        """Plot the convergence of the fit orbit parameters."""

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Plot with current parameters
        self.Steps = range(self.NbOrbits)
        self.EvalParamOrbit = eval(f'self.{self.ParamOrbit}')[self.nBody]

        self.Subplot.plot(self.Steps, self.EvalParamOrbit, marker=',', linestyle='')

        # Plot features
        self.Subplot.set_xlabel('Step')
        self.Subplot.set_ylabel(self.LabelOf(self.ParamOrbit)+' '+self.UnitOf(self.ParamOrbit))


class Hist(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitsParams):
        super().__init__('Histogram', "Histogram of orbital parameters", None, OutputParams, None, None, BestOrbitsParams, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True)

    def InitParams(self):
        """Initialize parameters for the Histogram tool."""

        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Variable', 'Variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0', 'Chi2', 'irel', 'other'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # i relatif between 2 bodies
        self.IrelWidget = QWidget()
        self.IrelLayout = QHBoxLayout()
        self.nBodyRelWidget = ComboBox('Relative body', 'Relative body to compare with the reference', self.ListBody)
        self.IrelLayout.addWidget(self.nBodyRelWidget)
        self.nBodyRelWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        self.IrelWidget.setLayout(self.IrelLayout)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.IrelWidget)
        self.IrelWidget.setVisible(False)

        # TextEdit for general formula
        self.FormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, m0, Chi2 with [n] for orbit number and usual mathematical functions', None)
        self.FormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.FormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.FormulaTextEdit)
        self.FormulaTextEdit.EditParam.textChanged.connect(self.reset_plots)

        # Connect ComboBox change event
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(lambda: self.FormulaTextEdit.setVisible(self.ParamOrbitWidget.ComboParam.currentText() == 'other'))
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(lambda: self.IrelWidget.setVisible(self.ParamOrbitWidget.ComboParam.currentText() == 'irel'))

        # Orbit number
        self.nBodyWidget = ComboBox(None, 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.ParamOrbitWidget.Layout.addWidget(self.nBodyWidget)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)

        # Histogram binning
        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)

        # Show confidence interval
        self.CheckMedian = CheckBox('Median :', 'Show the median and the 1 sigma confidence interval')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckMedian)

        # Confidence interval
        self.IntConf = 68
        self.IntConfWidget = SpinBox('Confidence', 'Acceptable level of confidence [%]', self.IntConf, 0, 100, 1)
        self.CheckMedian.Layout.addWidget(self.IntConfWidget)
        self.IntConfWidget.setEnabled(self.CheckMedian.CheckParam.isChecked())
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.IntConfWidget.setEnabled(state))
        self.PercentLbl = QLabel(' %')
        self.CheckMedian.Layout.addWidget(self.PercentLbl)
        self.PercentLbl.setEnabled(self.CheckMedian.CheckParam.isChecked())
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.PercentLbl.setEnabled(state))

    # def ToggleFormulaTextEdit(self):
    #     """Toggle the visibility of the formula text edit based on the ComboBox selection."""
    #     if self.ParamOrbitWidget.ComboParam.currentIndex() == self.ParamOrbitWidget.ComboParam.count() - 1:
    #         self.FormulaTextEdit.setVisible(True)
    #     else:
    #         self.FormulaTextEdit.setVisible(False)
            
    def evaluate_ParamOrbit(self, prefixe):
        """Evaluate the parameter orbit based on the current widget values."""
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        if self.ParamOrbit == 'irel' or self.ParamOrbit == 'other':
            if self.ParamOrbitWidget.ComboParam.currentText() == 'irel':
                if self.nBodyRel == self.nBody:
                    print('nBodyRel can not be the same as nBody')
                    return None
                else:
                    formula = f'arccos(cos(i[{self.nBody+1}])*cos(i[{self.nBodyRel+1}])+cos(W[{self.nBody+1}]-W[{self.nBodyRel+1}])*sin(i[{self.nBody+1}])*sin(i[{self.nBodyRel+1}]))'
            elif self.ParamOrbitWidget.ComboParam.currentText() == 'other':
                formula = self.FormulaTextEdit.EditParam.text()
            return self.evaluate_formula(formula, prefixe, self.nBody)
        else:
            self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
            # print('not irel or other')
            return eval(f'{prefixe}{self.ParamOrbit}')[self.nBody]

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.nBodyRel = int(self.nBodyRelWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.EvalParamOrbit = self.evaluate_ParamOrbit('self.')
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.IntConf = self.IntConfWidget.SpinParam.value()

    def Plot(self):
        """Plot the histogram based on the selected parameters."""

        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # Check if the parameters are valid for plotting
        if self.EvalParamOrbit is None:
            self.Subplot.figure.canvas.draw()
            return
        if np.var(self.EvalParamOrbit) == 0 or self.EvalParamOrbit[0] == float('inf'):
            print('No data to plot or variance is zero.')
            self.Subplot.figure.canvas.draw()
            return
    
        # Actual X limits on the plot or initial limits
        # Overwritten if modified by the user
        xlim_init = (np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))
        xlim = self.WidgetPlot.history[self.WidgetPlot.history_index]['xlim'] if self.WidgetPlot.history else xlim_init
        
        # Plot histogram
        self.Subplot.hist(self.EvalParamOrbit, self.NbBins, xlim)

        # Plot best fit
        if self.CheckBestFit.CheckParam.isChecked():
            BestParam = self.evaluate_ParamOrbit('self.Best')
            self.Subplot.axvline(BestParam, color='red')
            self.Subplot.text(BestParam, 0.5 * self.Subplot.get_ylim()[1], s='{}'.format(np.around(BestParam, 3)), color='r', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='r'), fontsize=9, horizontalalignment='center', verticalalignment='center', rotation=0)

        # Plot median and confidence interval
        if self.CheckMedian.CheckParam.isChecked():
            self.leftBound = xlim[0]
            self.rightBound = xlim[1]
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

            # self.Subplot.axvspan(np.min(self.EvalParamOrbit), self.leftBound, facecolor='grey', alpha=0.5)
            # self.Subplot.axvspan(self.rightBound, np.max(self.EvalParamOrbit), facecolor='grey', alpha=0.5)

        # Initial X and Y labels
        # Overwritten if modified by the user
        if self.ParamOrbitWidget.ComboParam.currentText() == 'irel':
            xlabel_init = r'i$_{rel}$ [°]'
        elif self.ParamOrbitWidget.ComboParam.currentText() == 'other':
            xlabel_init = self.FormulaTextEdit.EditParam.text()
        else:
            xlabel_init = self.LabelOf(self.ParamOrbit)+' '+self.UnitOf(self.ParamOrbit)
        self.Subplot.set_xlabel(xlabel_init)
        self.Subplot.set_ylabel('Count number')

        
class Hist2D(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitsParams):
        super().__init__('Histogram 2D', 'Histogram of an orbital parameter as function of another', None, OutputParams, None, None, BestOrbitsParams, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylim=True)

    def InitParams(self):
        """Initialize parameters for the 2D Histogram tool."""

        # Abscissa orbit parameters
        self.XParamOrbitWidget = ComboBox('X variable', 'Abscissa variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0', 'Chi2', 'other'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XParamOrbitWidget)
        self.XParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # Orbit number for X parameter
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.XnBodyWidget = ComboBox(None, 'Number of the orbit counting from the center of the system outwards', self.ListBody)
        self.XParamOrbitWidget.Layout.addWidget(self.XnBodyWidget)
        self.XnBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        # if self.NbBodies == 1:
        #     self.XnBodyWidget.setEnabled(False)

        # TextEdit for general x formula
        self.XFormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, m0, Chi2 with [n] for orbit number and usual mathematical functions', None)
        self.XFormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.XFormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XFormulaTextEdit)
        self.XFormulaTextEdit.EditParam.textChanged.connect(self.reset_plots)

        # Ordinate orbit parameters
        self.YParamOrbitWidget = ComboBox('Y variable', 'Ordinate variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0', 'Chi2', 'other'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YParamOrbitWidget)
        self.YParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # Orbit number for Y parameter
        self.YnBodyWidget = ComboBox(None, 'Number of the orbit counting from the center of the system outwards', self.ListBody)
        self.YParamOrbitWidget.Layout.addWidget(self.YnBodyWidget)
        self.YnBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        # if self.NbBodies == 1:
        #     self.YnBodyWidget.setEnabled(False)

        # TextEdit for general y formula
        self.YFormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, m0, Chi2 with [n] for orbit number and usual mathematical functions', None)
        self.YFormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.YFormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YFormulaTextEdit)
        self.YFormulaTextEdit.EditParam.textChanged.connect(self.reset_plots)

        # Connect ComboBox change event
        self.XParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.ToggleXFormulaTextEdit)
        self.YParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.ToggleYFormulaTextEdit)

        # Histogram binning
        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)

    def ToggleXFormulaTextEdit(self):
        """Toggle the visibility of the X formula text edit and orbit number ComboBox based on the ComboBox selection."""
        if self.XParamOrbitWidget.ComboParam.currentIndex() == self.XParamOrbitWidget.ComboParam.count() - 1:
            self.XFormulaTextEdit.setVisible(True)
            # self.XnBodyWidget.setVisible(False)
        else:
            self.XFormulaTextEdit.setVisible(False)
            # self.XnBodyWidget.setVisible(True)

    def ToggleYFormulaTextEdit(self):
        """Toggle the visibility of the Y formula text edit and orbit number ComboBox based on the ComboBox selection."""
        if self.YParamOrbitWidget.ComboParam.currentIndex() == self.YParamOrbitWidget.ComboParam.count() - 1:
            self.YFormulaTextEdit.setVisible(True)
            # self.YnBodyWidget.setVisible(False)
        else:
            self.YFormulaTextEdit.setVisible(False)
            # self.YnBodyWidget.setVisible(True)

    def evaluate_ParamOrbit(self, prefixe, param_widget, formula_text_edit, nBodyWidget):
        """Evaluate the parameter orbit based on the current widget values."""
        nBody = int(nBodyWidget.ComboParam.currentText()) - 1
        param_orbit = param_widget.ComboParam.currentText()
        if param_widget.ComboParam.currentIndex() == param_widget.ComboParam.count() - 1:
            formula = formula_text_edit.EditParam.text()
            return self.evaluate_formula(formula, prefixe, nBody)
        return eval(f'{prefixe}{param_orbit}')[nBody]

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.XParamOrbit = self.XParamOrbitWidget.ComboParam.currentText()
        self.YParamOrbit = self.YParamOrbitWidget.ComboParam.currentText()
        self.EvalXParamOrbit = self.evaluate_ParamOrbit('self.', self.XParamOrbitWidget, self.XFormulaTextEdit, self.XnBodyWidget)
        self.EvalYParamOrbit = self.evaluate_ParamOrbit('self.', self.YParamOrbitWidget, self.YFormulaTextEdit, self.YnBodyWidget)
        self.NbBins = self.NbBinsWidget.SpinParam.value()

    def Plot(self):
        """Plot the 2D histogram based on the selected parameters."""

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Check if the parameters are valid for plotting
        if self.EvalXParamOrbit is None or self.EvalYParamOrbit is None: 
            self.Subplot.figure.canvas.draw()
            return
        if np.var(self.EvalXParamOrbit) == 0 or np.var(self.EvalYParamOrbit) == 0 or self.EvalXParamOrbit[0] == float('inf') or self.EvalYParamOrbit[0] == float('inf'):
            self.Subplot.figure.canvas.draw()
            return
        
        # X, Y limits
        xlim_init=(np.min(self.EvalXParamOrbit), np.max(self.EvalXParamOrbit))
        ylim_init=(np.min(self.EvalYParamOrbit), np.max(self.EvalYParamOrbit))
        xlim = self.WidgetPlot.history[self.WidgetPlot.history_index]['xlim'] if len(self.WidgetPlot.history) != 0 else xlim_init
        ylim = self.WidgetPlot.history[self.WidgetPlot.history_index]['ylim'] if len(self.WidgetPlot.history) != 0 else ylim_init

        range = (xlim, ylim)

        hist = self.Subplot.hist2d(self.EvalXParamOrbit, self.EvalYParamOrbit, (self.NbBins, self.NbBins), range)
        ColorbarAx = make_axes_locatable(self.Subplot).append_axes('right', size='5%', pad=0.1)
        self.WidgetPlot.Canvas.fig.colorbar(hist[3], ColorbarAx, ticks=[], label='Count number')

        # Best fit
        if self.CheckBestFit.CheckParam.isChecked():
            BestXParam = self.evaluate_ParamOrbit('self.Best', self.XParamOrbitWidget, self.XFormulaTextEdit, self.XnBodyWidget)
            BestYParam = self.evaluate_ParamOrbit('self.Best', self.YParamOrbitWidget, self.YFormulaTextEdit, self.YnBodyWidget)
            self.Subplot.plot(BestXParam, BestYParam, color='red', marker='x')

        # Plot features
        if self.XParamOrbitWidget.ComboParam.currentIndex() == self.XParamOrbitWidget.ComboParam.count() - 1:
            self.Subplot.set_xlabel(self.XFormulaTextEdit.EditParam.text())
        else:
            self.Subplot.set_xlabel(self.LabelOf(self.XParamOrbit)+' '+self.UnitOf(self.XParamOrbit))
        if self.YParamOrbitWidget.ComboParam.currentIndex() == self.YParamOrbitWidget.ComboParam.count() - 1:
            self.Subplot.set_ylabel(self.YFormulaTextEdit.EditParam.text())
        else:
            self.Subplot.set_ylabel(self.LabelOf(self.YParamOrbit)+' '+self.UnitOf(self.YParamOrbit))


class Corner(GeneralToolClass):
    def __init__(self, SelectOrbitsParams, BestOrbitsParams):
        super().__init__('Corner', 'Corner plot of parameters', None, None, SelectOrbitsParams, None, BestOrbitsParams, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot)

    def InitParams(self):
        """Initialize parameters for the Corner plot tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)

        # Parameters to include in the corner plot
        self.ParamContainer = QWidget()
        self.ParamLayoutV = QVBoxLayout()

        self.ParamLayoutV.addWidget(QLabel('Variables :'))

        self.ParamCheckGroupBox = QGroupBox()
        self.ParamCheckLayout = QGridLayout()

        self.OrbitParams = [
            ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0'],
            ['Period', 'Semi-major axis', 'Eccentricity', 'Inclination', 'Argument of periastron', 'Longitude of ascending node', 'Periastron time passage', 'Body mass', 'Central body mass']
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

        # Short labels
        self.CheckShortLabels = CheckBox('Short labels', 'Show short labels')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckShortLabels)


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

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # No subplot because corner plots are standalone

        # Collect data for the corner plot
        Data = []
        DataNames = []
        DataLabels = []
        for x in self.WidgetOrbitParams:
            if x.CheckParam.isChecked():
                Data.append(eval(f'self.Select{x.CheckParam.text()}')[self.nBody])
                DataNames.append(x.CheckParam.text())
                if self.CheckShortLabels.CheckParam.isChecked():
                    DataLabels.append(x.CheckParam.text()+' '+self.UnitOf(x.CheckParam.text()))
                else:
                    DataLabels.append(self.LabelOf(x.CheckParam.text())+' '+self.UnitOf(x.CheckParam.text()))
        Data = np.array(Data).T

        # Check if there is data to plot
        if len(Data) == 0:
            print('No parameters selected or all parameters have no variance.')
            self.WidgetPlot.Canvas.fig.canvas.draw()
            return

        # Create corner plot
        grid = corner.corner(Data, labels=DataLabels, bins=self.NbBins, fig=self.WidgetPlot.Canvas.fig)

        # Adjust the labels to not be slanted
        for k in range(len(grid.get_axes())):
            ax = grid.get_axes()[k]
            ax.tick_params(axis='x', rotation=0)
            ax.tick_params(axis='y', rotation=0)

            # Best fit
            if self.CheckBestFit.CheckParam.isChecked():

                if k % (len(DataLabels) + 1) == 0:
                    param_index = k // (len(DataLabels) + 1)
                    BestParam = eval(f'self.Best{DataNames[param_index]}')[self.nBody]
                    ax.axvline(BestParam, color='red', linestyle='-', linewidth=0.75)
                else:
                    row = k // len(DataLabels)
                    col = k % len(DataLabels)
                    if row > col:
                        BestXParam = eval(f'self.Best{DataNames[col]}')[self.nBody]
                        BestYParam = eval(f'self.Best{DataNames[row]}')[self.nBody]
                        ax.plot(BestXParam, BestYParam, color='red', marker='x')


class PosAtDate(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Position at date', 'Position of bodies at a given date', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylim=True)

    def InitParams(self):
        """Initialize parameters for the Position at Date tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)

        # Date of wanted observation
        self.DateWidget = DateAndMJDEdit('Date', 'Date of wanted observation')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.DateWidget)
        self.DateWidget.MJDWidget.SpinParam.valueChanged.connect(self.refresh_plots)

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
        if self.InputData is None:
            self.CheckObs.CheckParam.setEnabled(False)

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.Date = self.DateWidget.MJDWidget.SpinParam.value()

    def Plot(self):
        """Plot the position at the given date based on the selected parameters."""

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Compute date
        SelectRaAtDate, SelectDecAtDate = [np.zeros(self.NbSelectOrbits) for _ in range(2)]

        for k in range(self.NbSelectOrbits):
            SelectPeriod = np.max(self.Selectt[self.nBody][k]) - np.min(self.Selectt[self.nBody][k])
            SelectDate = self.Date

            while SelectDate < np.min(self.Selectt[self.nBody][k]):
                SelectDate += SelectPeriod

            while SelectDate > np.max(self.Selectt[self.nBody][k]):
                SelectDate -= SelectPeriod

            indexBestDate = np.argmin(np.abs(self.Selectt[self.nBody][k] - SelectDate))
            SelectRaAtDate[k] = self.SelectRa[self.nBody][k][indexBestDate]
            SelectDecAtDate[k] = self.SelectDec[self.nBody][k][indexBestDate]

        # X and Y limits
        xlim_init = (np.max(self.SelectRa[self.nBody]), np.min(self.SelectRa[self.nBody])) # Inverted X axis for Ra
        ylim_init = (np.min(self.SelectDec[self.nBody]), np.max(self.SelectDec[self.nBody]))
        (Xmin, Xmax) = self.WidgetPlot.history[self.WidgetPlot.history_index]['xlim'] if len(self.WidgetPlot.history) != 0 else xlim_init
        ylim = self.WidgetPlot.history[self.WidgetPlot.history_index]['ylim'] if len(self.WidgetPlot.history) != 0 else ylim_init

        self.range = ((Xmax, Xmin), ylim) # Inverted X axis for Ra

        # Plot with current parameters
        hist = self.Subplot.hist2d(SelectRaAtDate, SelectDecAtDate, bins=(self.NbBins, self.NbBins), range=self.range)
        ColorbarAx = make_axes_locatable(self.Subplot).append_axes('right', size='5%', pad=0.1)
        self.WidgetPlot.Canvas.fig.colorbar(hist[3], ColorbarAx, ticks=[], label='Count number')

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
            ra = self.InputData['Planets']['DataAstrom']['Ra'][self.nBody]
            dec = self.InputData['Planets']['DataAstrom']['Dec'][self.nBody]
            dra = self.InputData['Planets']['DataAstrom']['dRa'][self.nBody]
            ddec = self.InputData['Planets']['DataAstrom']['dDec'][self.nBody]
            dates = self.InputData['Planets']['DataAstrom']['Date'][self.nBody]
            self.Subplot.errorbar(ra, dec, ddec, dra, linestyle='', color='white')

        # Plot features
        self.Subplot.set_xlabel(r'$\delta$ Ra [mas]')
        self.Subplot.set_ylabel(r'$\delta$ Dec [mas]')
        self.Subplot.invert_xaxis()
        self.Subplot.set_aspect('equal', adjustable='box')
        self.Subplot.set_xlim(xlim_init)
        self.Subplot.set_ylim(ylim_init)

