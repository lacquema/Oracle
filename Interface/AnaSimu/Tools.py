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
        self.WindowPlot.WidgetParam.BtnRefresh.clicked.connect(self.Refresh_active_plots)

        # Initialisation of Data
        self.InputData = InputData

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

        # # Fix the width to avoid resizing of parameters
        # left_width = self.WindowPlot.WidgetParam.sizeHint().width()
        # self.WindowPlot.WidgetParam.setFixedWidth(left_width)

        # Widget container
        self.setLayout(Layout)  # GeneralToolClass is directly the widget container

    def InitParams(self):
        """Initialize parameters. This method should be overridden by subclasses."""
        return

    def Refresh_active_plots(self):
        """Refresh all active plots when the refresh button is clicked."""
        for WidgetPlot in self.WindowPlot.WidgetPlots:
            if WidgetPlot.isVisible():
                WidgetPlot.refresh_plot()

    def Toggle_WindowPlot(self):
        """Open the plot window when the Plot button is clicked."""
        self.WindowPlot.show()
        self.BtnPlot.setEnabled(False)
        self.Refresh_active_plots()

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
        for widget_plot in self.WindowPlot.WidgetPlots:
            widget_plot.reset_history()
            # widget_plot.connect_events_to_reset_history()  # Reconnect events to reset history
            # print(widget_plot.events_to_reset_history)
            # self.WindowPlot.connect_events_to_reset_history(widget_plot, widget_plot.events_to_reset_history)
        self.Refresh_active_plots()

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
            'Mdyn': 'Dynamical mass',
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
            'Mdyn': '[Mjup]',
            'Chi2': ''
        }
        return units.get(var, 'Unknown variable')
        
    # Compute of the limits of subplot 2d
    def subplot_lim_2d(self, widget_plot, xlim_init=None, ylim_init=None):
        """
        Determine the axis limits for a subplot 2d based on the widget's history.
        """
        if not widget_plot.history:  # If history is empty
            if xlim_init is None: xlim_init = (None, None)
            if ylim_init is None: ylim_init = (None, None)
            return xlim_init, ylim_init
        else:
            # Retrieve limits from history
            index = widget_plot.history_index
            xlim = widget_plot.history[index].get('xlim', xlim_init)
            ylim = widget_plot.history[index].get('ylim', ylim_init)
            return xlim, ylim
        
    # Compute of the limits of subplot 2d
    def subplot_lim_3d(self, widget_plot, xlim_init=None, ylim_init=None, zlim_init=None):
        """
        Determine the axis limits for a subplot 3d based on the widget's history.
        """
        if not widget_plot.history:  # If history is empty
            if xlim_init is None: xlim_init = (None, None)
            if ylim_init is None: ylim_init = (None, None)
            if zlim_init is None: zlim_init = (None, None)
            return xlim_init, ylim_init, zlim_init
        else:
            # Retrieve limits from history
            index = widget_plot.history_index
            xlim = widget_plot.history[index].get('xlim', xlim_init)
            ylim = widget_plot.history[index].get('ylim', ylim_init)
            zlim = widget_plot.history[index].get('zlim', zlim_init)
            return xlim, ylim, zlim
        
    def replace_params_in_formula(self, formula, prefixe, nOrbitDefault):
        """Replace parameters and functions in the formula with their corresponding values."""
        print(formula)
        for num in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
            formula = formula.replace(f'[{num}]', f'[{str(int(num)-1)}]') # Replace [n] by [n-1]
        for param in ['Chi2', 'P', 'a', 'e', 'tp', 'm', 'Mdyn']:
            formula = re.sub(r'\b' + param + r'\b(?!\[)', f'{param}[{nOrbitDefault}]', formula) # Add [nOrbitDefault] to the parameter
            formula = re.sub(r'\b' + param + r'\b', f'{prefixe}{param}', formula) # Replace the parameter by its value
        for param in ['i', 'w', 'W']:
            formula = re.sub(r'\b' + param + r'\b(?!\[)', f'np.radians({param}[{nOrbitDefault}])', formula) # Convert to radians and add [nOrbitDefault]
            formula = re.sub(r'\b' + param + r'\b', f'np.radians({prefixe}{param})', formula) # Convert to radians and replace the parameter by its value
        for fonction in ['sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2', 'hypot', 'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh', 'exp', 'expm1', 'exp2', 'log', 'log10', 'log2', 'log1p', 'sqrt', 'square', 'cbrt', 'power', 'erf', 'erfc', 'gamma', 'lgamma', 'digamma', 'beta']:
            formula = re.sub(r'\b' + fonction + r'\b', f'np.{fonction}', formula) # Replace the function by its numpy equivalent
        # Convert the result to degrees if it is an angle in radians
        if any(f'np.{angle_func}' in formula for angle_func in ['arcsin', 'arccos', 'arctan', 'arctan2']):
            formula = f'np.degrees({formula})'
        print('Formula is: '+formula)
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



class SpaceView(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Space view', 'Space view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Parameters initialisation
        self.InitParams()

        # Window plots initialisation
        self.WidgetPlotXY = self.WindowPlot.add_WidgetPlot(self.PlotXY)
        self.WidgetPlotXZ = self.WindowPlot.add_WidgetPlot(self.PlotXZ)
        self.WidgetPlotXYZ = self.WindowPlot.add_WidgetPlot(self.PlotXYZ)

        self.indexViewChanged(self.indexView)
    
    def InitParams(self):
        """Initialize parameters for the SpaceView tool."""

        # Number of studied body
        self.ListBody = [str(k+1) for k in range(self.NbBodies)]
        if self.NbBodies > 1:
            self.ListBody.append('all')
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history) 

        # Number of shown orbits
        self.NbShownOrbits = 500
        self.NbShownOrbitsWidget = SpinBox('Number of orbits', 'Number of shown orbits', ParamDefault=self.NbShownOrbits, ParamMin=0, ParamMax=self.NbSelectOrbits)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbShownOrbitsWidget)

        # Type of view
        self.ViewWidget = ComboBox('View', 'Dimension', ['2D XY', '2D XZ', '3D'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ViewWidget)
        self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.ViewWidget.ComboParam.currentIndexChanged.connect(self.indexViewChanged)

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

    def reset_WidgetPlots_history(self):
        """Reset the history of the plot when the parameters are reset."""
        self.WidgetPlotXY.reset_history()
        self.WidgetPlotXZ.reset_history()
        self.WidgetPlotXYZ.reset_history()
        self.Refresh_active_plots()

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
        self.Refresh_active_plots()

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = 'all' if self.nBodyWidget.ComboParam.currentText() == 'all' else int(self.nBodyWidget.ComboParam.currentText()) - 1
        # self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.NbShownOrbits = self.NbShownOrbitsWidget.SpinParam.value()

    def general_plot(self):
        """Plot the orbits based on the selected parameters."""

        # Update parameters
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return

    def PlotXY(self):
        """Plot the 2D view of the orbits in the XY plane."""

        # General plot
        self.general_plot()

        # Add subplot
        self.SubplotXY = self.WidgetPlotXY.Canvas.fig.add_subplot(111, aspect='equal', label='Main plot')

        # X, Y limits
        (Xmin, Xmax), (Ymin, Ymax) = self.subplot_lim_2d(self.WidgetPlotXY)

        # Central star
        self.SubplotXY.plot(0, 0, marker='*', color='orange', markersize=10)

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
                    self.SubplotXY.errorbar(ra, dec, dra, ddec, linestyle='', color='black')
                    if self.CheckDateObs.CheckParam.isChecked():
                        if Xmin!=None: self.annotate_dates(dates, ra, dec, Xmin, Xmax, Ymin, Ymax)

            else:
                ra = self.InputData['Planets']['DataAstrom']['Ra'][self.nBody]
                dec = self.InputData['Planets']['DataAstrom']['Dec'][self.nBody]
                dra = self.InputData['Planets']['DataAstrom']['dRa'][self.nBody]
                ddec = self.InputData['Planets']['DataAstrom']['dDec'][self.nBody]
                dates = self.InputData['Planets']['DataAstrom']['Date'][self.nBody]
                self.SubplotXY.errorbar(ra, dec, dra, ddec, linestyle='', color='black')
                if self.CheckDateObs.CheckParam.isChecked():
                    if Xmin!=None: self.annotate_dates(dates, ra, dec, Xmin, Xmax, Ymin, Ymax)

        # Set axis
        self.SubplotXY.set_xlabel(r'$\delta Ra$ [mas]')
        self.SubplotXY.set_ylabel(r'$\delta Dec$ [mas]')
        self.SubplotXY.invert_xaxis()
        self.SubplotXY.set_aspect('equal', adjustable='box')
        self.SubplotXY.set_title(' ')


    def PlotXZ(self):
        """Plot the 2D view of the orbits in the XZ plane."""

        # General plot
        self.general_plot()

        # Add subplot
        self.SubplotXZ = self.WidgetPlotXZ.Canvas.fig.add_subplot(111, aspect='equal', label='Main plot')

        # # X, Y limits
        # (Xmin, Xmax), (Zmin, Zmax) = self.subplot_lim_2d(self.WidgetPlotXZ)

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
        self.SubplotXZ.set_title(' ')

    def PlotXYZ(self):
        """Helper function to plot 3D views."""

        # General plot
        self.general_plot()

        # Add subplot
        self.SubplotXYZ = self.WidgetPlotXYZ.Canvas.fig.add_subplot(111, aspect='equal', projection='3d', label='Main plot')

        # # X, Y, Z limits
        # (Xmin, Xmax), (Ymin, Ymax), (Zmin, Zmax) = self.subplot_lim_3d(self.WidgetPlotXYZ)

        # Central star
        self.SubplotXYZ.plot(0, 0, 0, marker='*', color='orange', markersize=10)

        if self.nBody == 'all':
            for k in range(self.NbBodies):
                if self.CheckBestFit.CheckParam.isChecked():
                    self.SubplotXYZ.plot(self.BestRa[k], self.BestDec[k], self.BestZ[k], color='r')
                for n in range(self.NbShownOrbits):
                    self.SubplotXYZ.plot(self.SelectRa[k][n], self.SelectDec[k][n], self.SelectZ[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
        else:
            if self.CheckBestFit.CheckParam.isChecked():
                self.SubplotXYZ.plot(self.BestRa[self.nBody], self.BestDec[self.nBody], self.BestZ[self.nBody], color='r')
            for n in range(self.NbShownOrbits):
                self.SubplotXYZ.plot(self.SelectRa[self.nBody][n], self.SelectDec[self.nBody][n], self.SelectZ[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)

        self.SubplotXYZ.set_xlabel(r'$\delta Ra$ [mas]')
        self.SubplotXYZ.set_ylabel(r'$\delta Dec$ [mas]')
        self.SubplotXYZ.set_zlabel('Depth [mas]')
        self.SubplotXYZ.invert_xaxis()
        self.SubplotXYZ.set_aspect('equal', adjustable='box')
        self.SubplotXYZ.set_title(' ')


    def annotate_dates(self, dates, ra, dec, Xmin, Xmax, Ymin, Ymax):
        """Annotate dates on the plot."""
        ra_range = abs(Xmax - Xmin)
        dec_range = abs(Ymax - Ymin)
        min_dist_x = 0.05 * ra_range
        min_dist_y = 0.1 * dec_range
        min_dist = np.hypot(min_dist_x, min_dist_y)
        annotated_points = []
        for idx in range(len(dates)):
            x_right = ra[idx] - 0.01 * ra_range 
            y = dec[idx] + 0.01 * dec_range
            x_left = ra[idx] + 0.01 * ra_range

            # Try right, if too close to previous, try left
            overlap = any(np.hypot(x_right - px, y - py) < min_dist for px, py in annotated_points)
            if not overlap:
                self.SubplotXY.annotate(f"{dates[idx]:.0f}", (x_right, y), color='black', fontsize=8, ha='left')
                annotated_points.append((x_right, y))
            else:
                annotated_points.append((x_left, y))
                if not overlap:
                    self.SubplotXY.annotate(f"{dates[idx]:.0f}", (x_left, y), color='black', fontsize=8, ha='right')
            



class TempoView(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Temporal view', 'Temporal view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Parameters initialisation
        self.InitParams()

        # Window plots initialisation
        self.WidgetPlot1 = self.WindowPlot.add_WidgetPlot(self.Plot1)
        self.WidgetPlot2 = self.WindowPlot.add_WidgetPlot(self.Plot2, layout=self.WidgetPlot1.Layout)

    def InitParams(self):
        """Initialize parameters for the TempoView tool."""

        # Orbit number
        self.ListBody = [str(k+1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)

        # Number of shown orbits
        self.NbShownOrbits = 500
        self.NbShownOrbitsWidget = SpinBox('Number of orbits', 'Number of shown orbits', ParamDefault=self.NbShownOrbits, ParamMin=0, ParamMax=self.NbSelectOrbits)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbShownOrbitsWidget)

        # Choice of coordinate
        self.CoordinateWidget = ComboBox('Choice of coordinate', 'Coordinates', ['dRa', 'dDec', 'Sep', 'Pa'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CoordinateWidget)
        self.CoordinateWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)

    
    def reset_WidgetPlots_history(self):
        """Reset the history of the plot when the parameters are reset."""
        self.WidgetPlot1.reset_history()
        self.WidgetPlot2.reset_history()
        self.Refresh_active_plots()

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

    def general_plot(self):
        """Plot the temporal view based on the selected parameters."""

        # Update parameters
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return

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

        # 3 periods of best fit
        self.Bestt3P = np.concatenate((self.Bestt[self.nBody] - self.BestP[self.nBody] * 365.25, self.Bestt[self.nBody], self.Bestt[self.nBody] + self.BestP[self.nBody] * 365.25))
        self.BestYplotOutput3P = np.concatenate((self.BestYplotOutput[self.nBody], self.BestYplotOutput[self.nBody], self.BestYplotOutput[self.nBody]))

        # Range of dates
        if self.InputData is not None: 
            self.DateRange = np.max(self.InputData['Planets']['DataAstrom']['Date'][self.nBody]) - np.min(self.InputData['Planets']['DataAstrom']['Date'][self.nBody])


    def Plot1(self):

        # General plot
        self.general_plot()

        # Add subplot
        self.Subplot1 = self.WidgetPlot1.Canvas.fig.add_subplot(111, label='Main plot')

        # X, Y limits
        if self.InputData is not None: 
            xlim_init = (np.min(self.InputData['Planets']['DataAstrom']['Date'][self.nBody]) - 0.1 * self.DateRange, np.max(self.InputData['Planets']['DataAstrom']['Date'][self.nBody]) + 0.1 * self.DateRange)
        else:
            xlim_init = None
        (Xmin, Xmax), (Ymin, Ymax) = self.subplot_lim_2d(self.WidgetPlot1, xlim_init)

        # Plot output data
        for n in range(self.NbShownOrbits):
            Selectt3P = np.concatenate((self.Selectt[self.nBody][n] - self.SelectP[self.nBody][n] * 365.25, self.Selectt[self.nBody][n], self.Selectt[self.nBody][n] + self.SelectP[self.nBody][n] * 365.25))
            YplotOutput3P = np.concatenate((self.YplotOutput[self.nBody][n], self.YplotOutput[self.nBody][n], self.YplotOutput[self.nBody][n]))
            self.Subplot1.plot(Selectt3P, YplotOutput3P, color=self.colorList[0], linestyle='-', linewidth=0.2, alpha=0.1)
        
        if self.InputData is not None: 
            for k in range(self.InputData['Planets']['NbDataAstrom'][self.nBody]):
                self.Subplot1.errorbar(self.InputData['Planets']['DataAstrom']['Date'][self.nBody][k], self.YplotInput[k], self.YplotInputErr[k], linestyle='', color='b')

        self.Subplot1.plot(self.Bestt3P, self.BestYplotOutput3P, linestyle='-', linewidth=0.5, color='r')

        # Plot features
        if self.CoordinateIndex == 3:
            self.Subplot1.set_ylabel(self.Coordinate + ' [°]')
        else:
            self.Subplot1.set_ylabel(self.Coordinate + ' [mas]')
        self.Subplot1.set_xlim(Xmin, Xmax)
        self.Subplot1.set_ylim(Ymin, Ymax)
        self.Subplot1.grid()
        self.Subplot1.set_title(' ')


    def Plot2(self):

        # General plot
        self.general_plot()

        # Add subplot
        self.Subplot2 = self.WidgetPlot2.Canvas.fig.add_subplot(111, label='Main plot')

        # X, Y limits
        (Xmin, Xmax), (Ymin, Ymax) = self.subplot_lim_2d(self.WidgetPlot2)

        # Plot input data and residuals
        if self.InputData is not None: 
            for k in range(self.InputData['Planets']['NbDataAstrom'][self.nBody]):
                indext = np.argmin(np.abs(self.Bestt3P - self.InputData['Planets']['DataAstrom']['Date'][self.nBody][k]))  # index of time of output data closer than time of input data
                Res = self.BestYplotOutput3P[indext] - self.YplotInput[k]  # Residual
                self.Subplot2.errorbar(self.Bestt3P[indext], Res, self.YplotInputErr[k], color='b')

        self.Subplot2.hlines(0, np.min(self.Bestt3P), np.max(self.Bestt3P), color='red', linewidth=0.5)

        # Plot features
        if self.CoordinateIndex == 3:
            self.Subplot2.set_ylabel(self.Coordinate + ' - Bestfit [°]')
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
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot)

    def InitParams(self):
        """Initialize parameters for the Convergence tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Variable', 'Orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn', 'Chi2'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)

    def reset_WidgetPlots_history(self):
        """Reset the history of the plot when the parameters are reset."""
        self.WidgetPlot.reset_history()
        self.Refresh_active_plots()

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()

    def Plot(self):
        """Plot the convergence of the fit orbit parameters."""

        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return

        # Plot with current parameters
        self.Steps = range(self.NbOrbits)
        self.EvalParamOrbit = eval(f'self.{self.ParamOrbit}')[self.nBody]

        self.Subplot.plot(self.Steps, self.EvalParamOrbit, marker=',', linestyle='')

        # Plot features
        self.Subplot.set_xlabel('Step')
        self.Subplot.set_ylabel(self.LabelOf(self.ParamOrbit)+' '+self.UnitOf(self.ParamOrbit))

        # Update canvas
        self.Subplot.set_title(' ')


class Hist(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitsParams):
        super().__init__('Histogram', "Histogram of orbital parameters", None, OutputParams, None, None, BestOrbitsParams, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot)

    def InitParams(self):
        """Initialize parameters for the Histogram tool."""

        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Variable', 'Variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn', 'Chi2', 'irel', 'other'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)

        # i relatif between 2 bodies
        self.IrelWidget = QWidget()
        self.IrelLayout = QHBoxLayout()
        self.nBodyRel = ComboBox('Relative body', 'Relative body to compare with the reference', self.ListBody)
        self.IrelLayout.addWidget(self.nBodyRel)
        self.nBodyRel.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)
        self.IrelWidget.setLayout(self.IrelLayout)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.IrelWidget)
        self.IrelWidget.setVisible(False)

        # TextEdit for general formula
        self.FormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, Mdyn, Chi2 with [n] for orbit number and usual mathematical functions', None)
        self.FormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.FormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.FormulaTextEdit)
        self.FormulaTextEdit.EditParam.textChanged.connect(self.reset_WidgetPlots_history)

        # Connect ComboBox change event
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(lambda: self.FormulaTextEdit.setVisible(self.ParamOrbitWidget.ComboParam.currentText() == 'other'))
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(lambda: self.IrelWidget.setVisible(self.ParamOrbitWidget.ComboParam.currentText() == 'irel'))

        # Orbit number
        self.nBodyWidget = ComboBox(None, 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.ParamOrbitWidget.Layout.addWidget(self.nBodyWidget)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)
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

        # Confidence interval bounds
        # self.EvalParamOrbit = self.evaluate_ParamOrbit('self.')
        
        # self.leftWidget = DoubleSpinBox(None, 'Left bound of the selected histogram', np.min(self.EvalParamOrbit), np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))
        # self.CheckMedian.Layout.addWidget(self.leftWidget)
        # self.leftWidget.setEnabled(self.CheckMedian.CheckParam.isChecked())
        # self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.leftWidget.setEnabled(state))

        # self.Llbl = QLabel(' < ')
        # self.CheckMedian.Layout.addWidget(self.Llbl)
        # self.Llbl.setEnabled(False)
        # self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.Llbl.setEnabled(state))

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

        # self.Rlbl = QLabel(' > ')
        # self.CheckMedian.Layout.addWidget(self.Rlbl)
        # self.Rlbl.setEnabled(False)
        # self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.Rlbl.setEnabled(state))

        # self.rightWidget = DoubleSpinBox(None, 'Right bound of the selected histogram', np.max(self.EvalParamOrbit), np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))
        # self.CheckMedian.Layout.addWidget(self.rightWidget)
        # self.rightWidget.setEnabled(self.CheckMedian.CheckParam.isChecked())
        # self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.rightWidget.setEnabled(state))

        # self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.ChangeRightandLeftBound)
        # self.FormulaTextEdit.EditParam.textChanged.connect(self.ChangeRightandLeftBound)

    def ToggleFormulaTextEdit(self):
        """Toggle the visibility of the formula text edit based on the ComboBox selection."""
        if self.ParamOrbitWidget.ComboParam.currentIndex() == self.ParamOrbitWidget.ComboParam.count() - 1:
            self.FormulaTextEdit.setVisible(True)
            # self.nBodyWidget.setVisible(False)
        else:
            self.FormulaTextEdit.setVisible(False)
            # self.nBodyWidget.setVisible(True)

    # def ChangeRightandLeftBound(self):
    #     """Update the bounds of the histogram when the orbital parameter changes."""
    #     try:
    #         self.EvalParamOrbit = self.evaluate_ParamOrbit('self.')
    #         if self.EvalParamOrbit is not None:
    #             min_val, max_val = np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit)
    #             self.leftWidget.SpinParam.setRange(min_val, max_val)
    #             self.leftWidget.SpinParam.setValue(min_val)
    #             self.rightWidget.SpinParam.setRange(min_val, max_val)
    #             self.rightWidget.SpinParam.setValue(max_val)
    #         else:
    #             self.leftWidget.SpinParam.setRange(0, 0)
    #             self.leftWidget.SpinParam.setValue(0)
    #             self.rightWidget.SpinParam.setRange(0, 0)
    #             self.rightWidget.SpinParam.setValue(0)
    #     except Exception as e:
    #         print(f"Error updating bounds: {e}")
            
    def evaluate_ParamOrbit(self, prefixe):
        """Evaluate the parameter orbit based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        print(self.ParamOrbitWidget.ComboParam.currentText())
        if self.ParamOrbit == 'irel' or self.ParamOrbit == 'other':
            if self.ParamOrbitWidget.ComboParam.currentText() == 'irel':
                nBodyRel = int(self.nBodyRel.ComboParam.currentText()) - 1
                if nBodyRel == self.nBody:
                    print('nBodyRel can not be the same as nBody')
                    return None
                else:
                    formula = f'arccos(cos(i[{self.nBody}])*cos(i[{nBodyRel}])+cos(W[{self.nBody}]-W[{nBodyRel}])*sin(i[{self.nBody}])*sin(i[{nBodyRel}]))'
                # return np.arccos(np.cos(self.i[self.nBody])*np.cos(self.i[nBodyRel])+np.cos(self.W[self.nBody]-self.W[nBodyRel])*np.sin(self.i[self.nBody])*np.sin(self.i[nBodyRel]))
            elif self.ParamOrbitWidget.ComboParam.currentText() == 'other':
                formula = self.FormulaTextEdit.EditParam.text()
            return self.evaluate_formula(formula, prefixe, self.nBody)
        else:
            self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
            print('not irel or other')
            return eval(f'{prefixe}{self.ParamOrbit}')[self.nBody]
    
    def reset_WidgetPlots_history(self):
        """Reset the history of the plot when the parameters are reset."""
        self.WidgetPlot.reset_history()
        self.WidgetPlot.plotting()
        # self.Refresh_active_plots()

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.EvalParamOrbit = self.evaluate_ParamOrbit('self.')
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.IntConf = self.IntConfWidget.SpinParam.value()
        # self.rightBound = self.rightWidget.SpinParam.value()
        # self.leftBound = self.leftWidget.SpinParam.value()

    def Plot(self):
        """Plot the histogram based on the selected parameters."""

        # Update parameters
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return
        if self.EvalParamOrbit is None or np.var(self.EvalParamOrbit) == 0 or self.EvalParamOrbit[0] == float('inf'):
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return
        
        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)
        
        # X, Y limits
        (Xmin, Xmax), (Ymin, Ymax) = self.subplot_lim_2d(self.WidgetPlot, xlim_init=(np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit)))
        
        # Range for histogram
        # if Xmin is None and Xmax is None: # I want keep None value at the beggining
        #     range = (np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))
        # else:

        range = (Xmin, Xmax)
        
        # Plot histogram
        self.Subplot.hist(self.EvalParamOrbit, self.NbBins, range)

        # Plot best fit
        if self.CheckBestFit.CheckParam.isChecked():
            BestParam = self.evaluate_ParamOrbit('self.Best')
            self.Subplot.axvline(BestParam, color='red')
            self.Subplot.text(BestParam, 0.5 * self.Subplot.get_ylim()[1], s='{}'.format(np.around(BestParam, 3)), color='r', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='r'), fontsize=9, horizontalalignment='center', verticalalignment='center', rotation=0)

        # Plot median and confidence interval
        if self.CheckMedian.CheckParam.isChecked():
            print(Xmin, Xmax)
            self.leftBound = Xmin
            self.rightBound = Xmax
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

        # Plot features
        if self.ParamOrbitWidget.ComboParam.currentText() == 'irel':
            self.Subplot.set_xlabel('i_rel [°]')
        if self.ParamOrbitWidget.ComboParam.currentText() == 'other':
            self.Subplot.set_xlabel(self.FormulaTextEdit.EditParam.text())
        else:
            self.Subplot.set_xlabel(self.LabelOf(self.ParamOrbit)+' '+self.UnitOf(self.ParamOrbit))
        self.Subplot.set_ylabel('Count number')
        self.Subplot.set_xlim(Xmin, Xmax)
        # self.Subplot.set_xlim(np.min(self.EvalParamOrbit), np.max(self.EvalParamOrbit))

        # Update canvas
        self.Subplot.set_title(' ')
        
    
class Hist2D(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitsParams):
        super().__init__('Histogram 2D', 'Histogram of an orbital parameter as function of another', None, OutputParams, None, None, BestOrbitsParams, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot)

    def InitParams(self):
        """Initialize parameters for the 2D Histogram tool."""

        # Abscissa orbit parameters
        self.XParamOrbitWidget = ComboBox('X variable', 'Abscissa variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn', 'Chi2', 'other'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XParamOrbitWidget)
        self.XParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)

        # Orbit number for X parameter
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.XnBodyWidget = ComboBox(None, 'Number of the orbit counting from the center of the system outwards', self.ListBody)
        self.XParamOrbitWidget.Layout.addWidget(self.XnBodyWidget)
        self.XnBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)
        # if self.NbBodies == 1:
        #     self.XnBodyWidget.setEnabled(False)

        # TextEdit for general x formula
        self.XFormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, Mdyn, Chi2 with [n] for orbit number and usual mathematical functions', None)
        self.XFormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.XFormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XFormulaTextEdit)
        self.XFormulaTextEdit.EditParam.textChanged.connect(self.reset_WidgetPlots_history)

        # Ordinate orbit parameters
        self.YParamOrbitWidget = ComboBox('Y variable', 'Ordinate variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn', 'Chi2', 'other'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YParamOrbitWidget)
        self.YParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)

        # Orbit number for Y parameter
        self.YnBodyWidget = ComboBox(None, 'Number of the orbit counting from the center of the system outwards', self.ListBody)
        self.YParamOrbitWidget.Layout.addWidget(self.YnBodyWidget)
        self.YnBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_WidgetPlots_history)
        # if self.NbBodies == 1:
        #     self.YnBodyWidget.setEnabled(False)

        # TextEdit for general y formula
        self.YFormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, Mdyn, Chi2 with [n] for orbit number and usual mathematical functions', None)
        self.YFormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.YFormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YFormulaTextEdit)
        self.YFormulaTextEdit.EditParam.textChanged.connect(self.reset_WidgetPlots_history)

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

    def reset_WidgetPlots_history(self):
        """Reset the history of the plot when the parameters are reset."""
        self.WidgetPlot.reset_history()
        # self.WidgetPlot.plotting()
        self.Refresh_active_plots()


    def Plot(self):
        """Plot the 2D histogram based on the selected parameters."""

        # Update parameters
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return
        
        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)
        
        # X, Y limits
        xlim_init=(np.min(self.EvalXParamOrbit), np.max(self.EvalXParamOrbit))
        ylim_init=(np.min(self.EvalYParamOrbit), np.max(self.EvalYParamOrbit))
        (Xmin, Xmax), (Ymin, Ymax) = self.subplot_lim_2d(self.WidgetPlot, xlim_init, ylim_init)

        range = ((Xmin, Xmax), (Ymin, Ymax))

        # Plot with current parameters
        if self.EvalXParamOrbit is None or self.EvalYParamOrbit is None:
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return

        hist = self.Subplot.hist2d(self.EvalXParamOrbit, self.EvalYParamOrbit, (self.NbBins, self.NbBins), range)
        ColorbarAx = make_axes_locatable(self.Subplot).append_axes('right', size='5%', pad=0.1)
        self.WidgetPlot.Canvas.fig.colorbar(hist[3], ColorbarAx, label='Count number')

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
        
        # Update canvas
        self.Subplot.set_title(' ')
        self.Subplot.set_xlim(Xmin, Xmax)
        self.Subplot.set_ylim(Ymin, Ymax)


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

        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return

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
        if len(Data) == 0:
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return

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


            # if k in (1,3,6,)

            # BestParam = eval('self.Best' + DataLabels[k])[self.nBody]
            # ax.axvline(BestParam, color='red')


class PosAtDate(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Position at date', 'Position of bodies at a given date', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot)

    def InitParams(self):
        """Initialize parameters for the Position at Date tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)

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
        if self.InputData is None:
            self.CheckObs.CheckParam.setEnabled(False)

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.Date = self.DateWidget.MJDWidget.SpinParam.value()

    def Plot(self):
        """Plot the position at the given date based on the selected parameters."""

        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return

        # Compute date
        SelectRaAtDate, SelectYAtDate = [np.zeros(self.NbSelectOrbits) for _ in range(2)]

        for k in range(self.NbSelectOrbits):
            SelectPeriod = np.max(self.Selectt[self.nBody][k]) - np.min(self.Selectt[self.nBody][k])
            SelectDate = self.Date

            while SelectDate < np.min(self.Selectt[self.nBody][k]):
                SelectDate += SelectPeriod

            while SelectDate > np.max(self.Selectt[self.nBody][k]):
                SelectDate -= SelectPeriod

            indexBestDate = np.argmin(np.abs(self.Selectt[self.nBody][k] - SelectDate))
            SelectRaAtDate[k] = self.SelectRa[self.nBody][k][indexBestDate]
            SelectYAtDate[k] = self.SelectDec[self.nBody][k][indexBestDate]

        # Plot with current parameters
        hist = self.Subplot.hist2d(SelectRaAtDate, SelectYAtDate, bins=(self.NbBins, self.NbBins), range=((np.min(self.SelectRa), np.max(self.SelectRa)), (np.min(self.SelectDec), np.max(self.SelectDec))))
        ColorbarAx = make_axes_locatable(self.Subplot).append_axes('right', size='5%', pad=0.1)
        self.WidgetPlot.Canvas.fig.colorbar(hist[3], ColorbarAx, ticks=[], label='Probability')

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
            self.Subplot.errorbar(self.InputData['Planets']['DataAstrom']['Ra'], self.InputData['Planets']['DataAstrom']['Dec'], self.InputData['Planets']['DataAstrom']['dRa'], self.InputData['Planets']['DataAstrom']['dDec'], linestyle='')  # Observed data

        # Plot features
        self.Subplot.set_xlabel(r'$\delta$ Ra [mas]')
        self.Subplot.set_ylabel(r'$\delta$ Dec [mas]')
        self.Subplot.invert_xaxis()
        self.Subplot.set_aspect('equal', adjustable='box')

        # Update canvas
        self.Subplot.set_title(' ')

