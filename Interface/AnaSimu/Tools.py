#! /Users/lacquema/ByeGildas/bin/python3

### --- Packages --- ###

# Transverse packages
import sys
import logging
import warnings
from matplotlib.pyplot import subplots_adjust
import numpy as np
from random import random, randint
import corner
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re
from matplotlib.transforms import Bbox

# PyQt packages
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QDateEdit, QGroupBox, QGridLayout
from PyQt6.QtCore import QDateTime, QDate, QSize, QTimer
from Utils import date_to_jd, jd_to_mjd, mjd_to_jd, jd_to_date

# My packages
from Parameters import *
from SelectOrbits import SelectOrbitsClass
from TransferData import BuildOrbitFilterMask

from WindowPlot import WindowPlot
import itertools


class _MatplotlibFinitePosFilter(logging.Filter):
    def filter(self, record):
        return 'posx and posy should be finite values' not in record.getMessage()


logging.getLogger('matplotlib.text').addFilter(_MatplotlibFinitePosFilter())
warnings.filterwarnings('ignore', category=RuntimeWarning)


### --- Tools Generating --- ###

class GeneralToolClass(QWidget):

    def __init__(self, ToolName, ToolStatus, InputData, OutputParams, SelectOrbitsParams, SelectOrbitsEllipses, BestOrbitParams, BestOrbitEllipse, LMOrbitParams, LMOrbitEllipse):
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
            (self.NbBodies, self.NbOrbits, self.PlanetsMassUnit, self.P, self.a, self.e, self.i, self.w, self.W, self.tp, self.m, self.m0, self.V0, self.Jitter, self.Chi2, self.map) = OutputParams

        if SelectOrbitsParams is not None:
            (
                self.NbBodies,
                self.NbSelectOrbits,
                self.PlanetsMassUnit,
                self.SelectP,
                self.Selecta,
                self.Selecte,
                self.Selecti,
                self.Selectw,
                self.SelectW,
                self.Selecttp,
                self.Selectm,
                self.Selectm0,
                self.SelectChi2,
                *SelectExtraParams,
            ) = SelectOrbitsParams
            self.SelectV0 = SelectExtraParams[0] if len(SelectExtraParams) > 0 else None
            self.SelectJitter = SelectExtraParams[1] if len(SelectExtraParams) > 1 else None
            self.SelectMap = SelectExtraParams[2] if len(SelectExtraParams) > 2 else None

        if SelectOrbitsEllipses is not None:
            (self.NbBodies, self.NbSelectOrbits, self.NbPtsEllipse, self.SelectP, self.Selectt, self.SelectRa, self.SelectDec, self.SelectZ, self.SelectSep, self.SelectPa, self.SelectRV) = SelectOrbitsEllipses

        if BestOrbitParams is not None:
            (self.NbBodies, self.PlanetsMassUnit, self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.Bestm0, self.BestChi2) = BestOrbitParams

        if BestOrbitEllipse is not None:
            (self.NbBodies, self.NbPtsEllipse, self.BestP, self.Bestt, self.BestRa, self.BestDec, self.BestZ, self.BestSep, self.BestPa, self.BestRV) = BestOrbitEllipse

        if LMOrbitParams is not None:
            (self.NbBodies, self.PlanetsMassUnit, self.LMP, self.LMa, self.LMe, self.LMi, self.LMw, self.LMW, self.LMtp, self.LMm, self.LMm0, self.LMChi2) = LMOrbitParams  

        if LMOrbitEllipse is not None:
            (self.NbBodies, self.NbPtsEllipse, self.LMP, self.LMt, self.LMRa, self.LMDec, self.LMZ, self.LMSep, self.LMPa, self.LMRV) = LMOrbitEllipse

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

    @staticmethod
    def _format_short_label_index(index_label):
        if index_label is None:
            return None

        cleaned_label = str(index_label).strip()
        if len(cleaned_label) == 0:
            return None

        cleaned_label = cleaned_label.replace('\\', r'\backslash ')
        cleaned_label = cleaned_label.replace('{', r'\{').replace('}', r'\}')
        cleaned_label = cleaned_label.replace('_', r'\_').replace(' ', r'\ ')
        return cleaned_label

    def ShortLabelOf(self, var=str, index_label=None):
        """Return the short label for a given variable."""
        formatted_index = self._format_short_label_index(index_label)

        if formatted_index is None:
            labels = {
                'P': r'$P$',
                'a': r'$a$',
                'e': r'$e$',
                'i': r'$i$',
                'w': r'$\omega$',
                'W': r'$\Omega$',
                'tp': r'$T_\mathrm{p}$',
                'm': r'$m$',
                'm0': r'$m_0$',
                'Chi2': r'$\chi^2$'
            }
            return labels.get(var, 'Unknown variable')

        labels = {
            'P': rf'$P_\mathrm{{{formatted_index}}}$',
            'a': rf'$a_\mathrm{{{formatted_index}}}$',
            'e': rf'$e_\mathrm{{{formatted_index}}}$',
            'i': rf'$i_\mathrm{{{formatted_index}}}$',
            'w': rf'$\omega_\mathrm{{{formatted_index}}}$',
            'W': rf'$\Omega_\mathrm{{{formatted_index}}}$',
            'tp': rf'$T_{{p,\mathrm{{{formatted_index}}}}}$',
            'm': rf'$m_\mathrm{{{formatted_index}}}$',
            'm0': rf'$m_{{0,\mathrm{{{formatted_index}}}}}$',
            'Chi2': rf'$\chi^2_\mathrm{{{formatted_index}}}$'
        }
        return labels.get(var, 'Unknown variable')

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
            'a': '[au]',
            'e': '',
            'i': '[deg]',
            'w': '[deg]',
            'W': '[deg]',
            'tp': '[MJD]',
            'm': f'[{self.PlanetsMassUnit}]',
            'm0': '[Msun]',
            'Chi2': ''
        }
        return units.get(var, 'Unknown variable')
        
    @staticmethod
    def _replace_params_in_expression(expression, prefixe, nOrbitDefault, convert_angles_to_radians):
        """Translate user parameters to Python expressions understood by the code."""
        if expression is None:
            return None

        default_index = max(int(nOrbitDefault), 0)
        linear_params = ['Chi2', 'P', 'a', 'e', 'tp', 'm', 'm0']
        angular_params = ['i', 'w', 'W']
        all_params = sorted(linear_params + angular_params, key=len, reverse=True)
        parameter_pattern = '|'.join(map(re.escape, all_params))

        def replace_parameter(match):
            param = match.group('param')
            raw_index = match.group('index')
            if raw_index is None or len(raw_index.strip()) == 0:
                index = default_index
            else:
                index = max(int(raw_index.strip()) - 1, 0)

            replacement = f'{prefixe}{param}[{index}]'
            if convert_angles_to_radians and param in angular_params:
                replacement = f'np.radians({replacement})'
            return replacement

        expression = re.sub(
            rf'\b(?P<param>{parameter_pattern})\b(?:\[\s*(?P<index>\d*)\s*\])?',
            replace_parameter,
            expression,
        )

        for fonction in ['sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2', 'hypot', 'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh', 'exp', 'expm1', 'exp2', 'log', 'log10', 'log2', 'log1p', 'sqrt', 'square', 'cbrt', 'power', 'erf', 'erfc', 'gamma', 'lgamma', 'digamma', 'beta', 'radians', 'degrees']:
            expression = re.sub(rf'(?<![\w.]){re.escape(fonction)}\b', f'np.{fonction}', expression)

        for const, replacement in {
            'G': '6.67430e-11',
            'Msun': '1.98847e30',
            'Mjup': '1.89818e27',
            'pi': 'np.pi',
            'e': 'np.e'
        }.items():
            expression = re.sub(rf'(?<![\w.]){re.escape(const)}\b', replacement, expression)

        if convert_angles_to_radians and any(f'np.{angle_func}' in expression for angle_func in ['arcsin', 'arccos', 'arctan', 'arctan2']):
            expression = f'np.degrees({expression})'

        return expression

    def replace_params_in_formula(self, formula, prefixe, nOrbitDefault):
        """Replace parameters and functions in the formula with their corresponding values."""
        return self._replace_params_in_expression(formula, prefixe, nOrbitDefault, convert_angles_to_radians=True)

    def replace_params_in_condition(self, condition, prefixe, nOrbitDefault):
        """Translate a boolean mask condition to the local data notation."""
        return self._replace_params_in_expression(condition, prefixe, nOrbitDefault, convert_angles_to_radians=False)
    
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
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitEllipse, LMOrbitEllipse):
        super().__init__('Space view', 'Space view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitEllipse, None, LMOrbitEllipse)

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

        # Show LM fit
        self.CheckLMFit = CheckBox('LM fit', 'Show the Levenberg-Marquardt fit')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckLMFit)
        self.CheckLMFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)
        self.CheckBestFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show observations points
        self.CheckObs = CheckBox('Observations', 'Show the observations points with its error bar')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckObs)
        if self.InputData is None:
            self.CheckObs.CheckParam.setEnabled(False)
        self.ViewWidget.ComboParam.currentIndexChanged.connect(lambda: self.CheckObs.CheckParam.setEnabled(self.ViewWidget.ComboParam.currentIndex() == 0))  # Enable observations only for 2D view
        self.CheckObs.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show date of observations
        self.CheckDateObs = CheckBox('Date of observations', 'Show the date of observations')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckDateObs)
        self.CheckDateObs.CheckParam.setEnabled(False)  # Initially disabled, enabled in Plot method if observations are shown
        self.ViewWidget.ComboParam.currentIndexChanged.connect(lambda: self.CheckDateObs.CheckParam.setEnabled(self.ViewWidget.ComboParam.currentIndex() == 0 and self.CheckObs.CheckParam.isChecked()))
        self.CheckObs.CheckParam.stateChanged.connect(lambda bool: self.CheckDateObs.CheckParam.setEnabled(bool))
        self.CheckDateObs.CheckParam.stateChanged.connect(self.refresh_plots)

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
                if self.CheckLMFit.CheckParam.isChecked():
                    self.SubplotXY.plot(self.LMRa[k], self.LMDec[k], color='orange', linewidth=1, label='LM fit' if k==0 else None)
                if self.CheckBestFit.CheckParam.isChecked():
                    self.SubplotXY.plot(self.BestRa[k], self.BestDec[k], color='r', linewidth=1, label='Best fit' if k==0 else None)
        else:
            for n in range(self.NbShownOrbits):
                self.SubplotXY.plot(self.SelectRa[self.nBody][n], self.SelectDec[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)
            if self.CheckLMFit.CheckParam.isChecked():
                self.SubplotXY.plot(self.LMRa[self.nBody], self.LMDec[self.nBody], color='orange', linewidth=1, label='LM fit')
            if self.CheckBestFit.CheckParam.isChecked():
                self.SubplotXY.plot(self.BestRa[self.nBody], self.BestDec[self.nBody], color='r', linewidth=1, label='Best fit')

        # Add observations points if available
        if self.CheckObs.CheckParam.isChecked() and self.indexView == 0:
            if self.nBody == 'all':
                for k in range(self.InputData['Planets']['Nb']):
                    ra = self.InputData['Planets']['DataAstrom']['Ra'][k]
                    dec = self.InputData['Planets']['DataAstrom']['Dec'][k]
                    dra = self.InputData['Planets']['DataAstrom']['dRA'][k]
                    ddec = self.InputData['Planets']['DataAstrom']['dDec'][k]
                    dates = self.InputData['Planets']['DataAstrom']['Date'][k]
                    self.SubplotXY.errorbar(ra, dec, ddec, dra, linestyle='', color='blue', linewidth=0.7)
                    if self.CheckDateObs.CheckParam.isChecked():
                        if Xmin!=None and Ymin!=None: 
                            self.annotate_dates(dates, ra, dec, Xmin, Xmax, Ymin, Ymax)

            else:
                ra = self.InputData['Planets']['DataAstrom']['Ra'][self.nBody]
                dec = self.InputData['Planets']['DataAstrom']['Dec'][self.nBody]
                dra = self.InputData['Planets']['DataAstrom']['dRA'][self.nBody]
                ddec = self.InputData['Planets']['DataAstrom']['dDec'][self.nBody]
                dates = self.InputData['Planets']['DataAstrom']['Date'][self.nBody]
                self.SubplotXY.errorbar(ra, dec, ddec, dra, linestyle='', color='blue', linewidth=1)
                if self.CheckDateObs.CheckParam.isChecked():
                    if Xmin!=None and Ymin!=None: 
                        self.annotate_dates(dates, ra, dec, Xmin, Xmax, Ymin, Ymax)
        
        # Set axis
        self.SubplotXY.set_xlabel(r'$\delta$RA [mas]')
        self.SubplotXY.set_ylabel(r'$\delta$Dec [mas]')
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
                if self.CheckLMFit.CheckParam.isChecked():
                    self.SubplotXZ.plot(self.LMRa[k], self.LMZ[k], color='orange', linewidth=1, label='LM fit' if k==0 else None)
                if self.CheckBestFit.CheckParam.isChecked():
                    self.SubplotXZ.plot(self.BestRa[k], self.BestZ[k], color='r', linewidth=1, label='Best fit' if k==0 else None)
        else:
            for n in range(self.NbShownOrbits):
                self.SubplotXZ.plot(self.SelectRa[self.nBody][n], self.SelectZ[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)
            if self.CheckLMFit.CheckParam.isChecked():
                self.SubplotXZ.plot(self.LMRa[self.nBody], self.LMZ[self.nBody], color='orange', linewidth=1, label='LM fit')
            if self.CheckBestFit.CheckParam.isChecked():
                self.SubplotXZ.plot(self.BestRa[self.nBody], self.BestZ[self.nBody], color='r', linewidth=1, label='Best fit')

        # Set axis
        self.SubplotXZ.set_xlabel(r'$\delta$RA [mas]')
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
                if self.CheckLMFit.CheckParam.isChecked():
                    self.SubplotXYZ.plot(self.LMRa[k], self.LMDec[k], self.LMZ[k], color='orange', linewidth=0.5, label='LM fit' if k==0 else None)
                if self.CheckBestFit.CheckParam.isChecked():
                    self.SubplotXYZ.plot(self.BestRa[k], self.BestDec[k], self.BestZ[k], color='r', linewidth=0.5, label='Best fit' if k==0 else None)
                for n in range(self.NbShownOrbits):
                    self.SubplotXYZ.plot(self.SelectRa[k][n], self.SelectDec[k][n], self.SelectZ[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
        else:
            if self.CheckLMFit.CheckParam.isChecked():
                self.SubplotXYZ.plot(self.LMRa[self.nBody], self.LMDec[self.nBody], self.LMZ[self.nBody], color='orange', linewidth=1, label='LM fit')
            if self.CheckBestFit.CheckParam.isChecked():
                self.SubplotXYZ.plot(self.BestRa[self.nBody], self.BestDec[self.nBody], self.BestZ[self.nBody], color='r', linewidth=1, label='Best fit')
            for n in range(self.NbShownOrbits):
                self.SubplotXYZ.plot(self.SelectRa[self.nBody][n], self.SelectDec[self.nBody][n], self.SelectZ[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)

        self.SubplotXYZ.set_xlabel(r'$\delta$RA [mas]')
        self.SubplotXYZ.set_ylabel(r'$\delta$Dec [mas]')
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
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitEllipse, LMOrbitEllipse):
        super().__init__('Temporal view', 'Temporal view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitEllipse, None, LMOrbitEllipse)

        self._tempo_xlim_sync_in_progress = False

        # Window plots initialisation
        VertLayoutPlots = QVBoxLayout()
        VertLayoutPlots.setContentsMargins(0, 0, 0, 0)
        VertLayoutPlots.setSpacing(0)
        ContainerPlots = QWidget()
        ContainerPlots.setLayout(VertLayoutPlots)
        self.WindowPlot.Splitter.addWidget(ContainerPlots)
        self.WidgetPlot1 = self.WindowPlot.add_WidgetPlot(self.Plot1, xlim=True, ylim=True, layout=VertLayoutPlots)
        self.WidgetPlot2 = self.WindowPlot.add_WidgetPlot(self.Plot2, layout=VertLayoutPlots)

        self.WidgetPlot1.setFixedHeight(430)
        self.WidgetPlot2.setFixedHeight(430)

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
        self.CoordinateWidget = ComboBox('Choice of coordinate', 'Coordinates', ['dRA', 'dDec', 'Sep', 'Pa', 'RV'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CoordinateWidget)
        self.CoordinateWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # Choice of reference solution for temporal model and residuals
        ListReference = ['Best fit', 'LM fit']
        self.ReferenceWidget = ComboBox('Reference orbit', 'Reference solution used for residuals', ListReference)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ReferenceWidget)
        self.ReferenceWidget.ComboParam.currentIndexChanged.connect(self.refresh_plots)


    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.CoordinateIndex = self.CoordinateWidget.ComboParam.currentIndex()
        if self.CoordinateIndex == 0:
            self.Coordinate = r'$\delta$RA'
        elif self.CoordinateIndex == 1:
            self.Coordinate = r'$\delta$Dec'
        elif self.CoordinateIndex == 2:
            self.Coordinate = 'Sep'
        elif self.CoordinateIndex == 3:
            self.Coordinate = 'Pa'
        elif self.CoordinateIndex == 4:
            self.Coordinate = 'RV'
        self.nBody = int(self.nBodyWidget.ComboParam.currentText()) - 1
        self.NbShownOrbits = self.NbShownOrbitsWidget.SpinParam.value()
        self.Reference = self.ReferenceWidget.ComboParam.currentText()

    def _get_input_dates(self):
        if self.InputData is None:
            return np.array([])
        if self.CoordinateIndex != 4:
            return np.asarray(self.InputData['Planets']['DataAstrom']['Date'][self.nBody])
        if self.InputData['Planets']['NbDataRV'][self.nBody] == 0:
            return np.array([])
        return np.asarray(self.InputData['Planets']['DataRV']['Date'][self.nBody])

    def _build_periodic_series(self, times, values, period_years, window_min=None, window_max=None):
        times = np.asarray(times)
        values = np.asarray(values)

        if times.size == 0:
            return np.array([]), np.array([])

        period_days = float(period_years) * 365.25
        if not np.isfinite(period_days) or period_days <= 0:
            order = np.argsort(times)
            return times[order], values[order]

        if window_min is None or window_max is None:
            window_min = np.min(times) - period_days/2
            window_max = np.max(times) + period_days/2

        shift_min = int(np.ceil((window_min - np.max(times)) / period_days))
        shift_max = int(np.floor((window_max - np.min(times)) / period_days))

        if shift_min > shift_max:
            shift_min, shift_max = 0, 0

        time_chunks = []
        value_chunks = []
        for shift in range(shift_min, shift_max + 1):
            shifted_times = times + shift * period_days
            mask = (shifted_times >= window_min) & (shifted_times <= window_max)
            if np.any(mask):
                time_chunks.append(shifted_times[mask])
                value_chunks.append(values[mask])

        if len(time_chunks) == 0:
            return np.array([]), np.array([])

        expanded_times = np.concatenate(time_chunks)
        expanded_values = np.concatenate(value_chunks)
        order = np.argsort(expanded_times)
        expanded_times = expanded_times[order]
        expanded_values = expanded_values[order]
        unique_times, unique_indices = np.unique(expanded_times, return_index=True)
        return unique_times, expanded_values[unique_indices]

    def _reference_values_at(self, times):
        if len(self.RefTimesWindow) == 0:
            return np.array([])
        return np.interp(times, self.RefTimesWindow, self.RefValuesWindow)

    def _sync_widget_history_xlim(self, widget_plot, xlim):
        if len(widget_plot.history) == 0:
            return
        synced_state = dict(widget_plot.history[widget_plot.history_index])
        synced_state['xlim'] = tuple(xlim)
        widget_plot.history[widget_plot.history_index] = synced_state

    def _sync_tempo_xlim(self, source_widget, target_widget):
        if self._tempo_xlim_sync_in_progress:
            return
        if not source_widget.Canvas.fig.axes or not target_widget.Canvas.fig.axes:
            return

        source_xlim = tuple(source_widget.Canvas.fig.axes[0].get_xlim())
        target_xlim = tuple(target_widget.Canvas.fig.axes[0].get_xlim())
        if np.allclose(source_xlim, target_xlim):
            return

        self._tempo_xlim_sync_in_progress = True
        try:
            self._sync_widget_history_xlim(target_widget, source_xlim)
            target_widget.Canvas.fig.axes[0].set_xlim(source_xlim)
            target_widget.Canvas.draw_idle()
        finally:
            self._tempo_xlim_sync_in_progress = False

    def _connect_tempo_xlim_sync(self, subplot, source_widget, target_widget):
        subplot.callbacks.connect('xlim_changed', lambda ax: self._sync_tempo_xlim(source_widget, target_widget))

    def general_plot(self):
        """Plot the temporal view based on the selected parameters."""
        self.YplotInput = []
        self.YplotInputErr = []
        # Determine the coordinate to plot
        if self.CoordinateIndex == 0:
            self.YplotOutput = self.SelectRa 
            self.BestYplotOutput = self.BestRa
            self.LMYplotOutput = self.LMRa
            if self.InputData is not None: 
                self.YplotInput = self.InputData['Planets']['DataAstrom']['Ra'][self.nBody]
                self.YplotInputErr = self.InputData['Planets']['DataAstrom']['dRA'][self.nBody]
        elif self.CoordinateIndex == 1:
            self.YplotOutput = self.SelectDec
            self.BestYplotOutput = self.BestDec
            self.LMYplotOutput = self.LMDec
            if self.InputData is not None: 
                self.YplotInput = self.InputData['Planets']['DataAstrom']['Dec'][self.nBody]
                self.YplotInputErr = self.InputData['Planets']['DataAstrom']['dDec'][self.nBody]
        elif self.CoordinateIndex == 2:
            self.YplotOutput = self.SelectSep
            self.BestYplotOutput = self.BestSep
            self.LMYplotOutput = self.LMSep
            if self.InputData is not None: 
                self.YplotInput = self.InputData['Planets']['DataAstrom']['Sep'][self.nBody]
                self.YplotInputErr = self.InputData['Planets']['DataAstrom']['dSep'][self.nBody]
        elif self.CoordinateIndex == 3:
            self.YplotOutput = self.SelectPa
            self.BestYplotOutput = self.BestPa
            self.LMYplotOutput = self.LMPa
            if self.InputData is not None: 
                self.YplotInput = self.InputData['Planets']['DataAstrom']['Pa'][self.nBody]
                self.YplotInputErr = self.InputData['Planets']['DataAstrom']['dPa'][self.nBody]
        elif self.CoordinateIndex == 4:
            self.YplotOutput = self.SelectRV
            self.BestYplotOutput = self.BestRV
            self.LMYplotOutput = self.LMRV
            if self.InputData is not None:
                if self.InputData['Planets']['NbDataRV'][self.nBody] != 0:
                    self.YplotInput = self.InputData['Planets']['DataRV']['RV'][self.nBody]
                    self.YplotInputErr = self.InputData['Planets']['DataRV']['dRV'][self.nBody]

        # Selected reference orbit (Best or LM)
        self.ReferenceLabel = 'Best fit'
        self.ReferenceColor = 'r'
        self.Reft = self.Bestt[self.nBody]
        self.RefP = self.BestP[self.nBody]
        self.RefYplotOutput = self.BestYplotOutput[self.nBody]

        if self.Reference == 'LM fit':
            self.ReferenceLabel = 'LM fit'
            self.ReferenceColor = 'orange'
            self.Reft = self.LMt[self.nBody]
            self.RefP = self.LMP[self.nBody]
            self.RefYplotOutput = self.LMYplotOutput[self.nBody]

        self.InputDates = self._get_input_dates()
        self.TimeMargin = self.RefP * 365.25
        if len(self.InputDates) != 0:
            self.TimeWindowMin = np.min(self.InputDates) - self.TimeMargin
            self.TimeWindowMax = np.max(self.InputDates) + self.TimeMargin
        else:
            self.TimeWindowMin = np.min(self.Reft) - self.TimeMargin
            self.TimeWindowMax = np.max(self.Reft) + self.TimeMargin

        self.RefTimesWindow, self.RefValuesWindow = self._build_periodic_series(
            self.Reft,
            self.RefYplotOutput,
            self.RefP,
            self.TimeWindowMin,
            self.TimeWindowMax,
        )

    def Plot1(self):

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot1)

        # Add subplot
        self.Subplot1 = self.WidgetPlot1.Canvas.fig.add_subplot(111, label='Main plot')

        # General plot
        self.general_plot()

        # X current limits
        xlim_init = (self.TimeWindowMin, self.TimeWindowMax)
        xlim = self.WidgetPlot1.history[self.WidgetPlot1.history_index]['xlim'] if len(self.WidgetPlot1.history)!=0 else xlim_init

        # Plot output data
        for n in range(self.NbShownOrbits):
            SelectTimesWindow, SelectValuesWindow = self._build_periodic_series(
                self.Selectt[self.nBody][n],
                self.YplotOutput[self.nBody][n],
                self.SelectP[self.nBody][n],
                xlim[0],
                xlim[1],
            )
            if len(SelectTimesWindow) != 0:
                self.Subplot1.plot(SelectTimesWindow, SelectValuesWindow, color=self.colorList[self.nBody], linestyle='-', linewidth=0.2, alpha=0.1)
        
        if len(self.YplotInput) != 0: 
            for k in range(len(self.InputDates)):
                self.Subplot1.errorbar(self.InputDates[k], self.YplotInput[k], self.YplotInputErr[k], linestyle='', color='b')

        self.Subplot1.plot(self.RefTimesWindow, self.RefValuesWindow, linestyle='-', linewidth=1, color=self.ReferenceColor)

        # Plot features
        if self.CoordinateIndex == 3:
            self.Subplot1.set_ylabel(self.Coordinate + ' [deg]')
        elif self.CoordinateIndex == 4:
            self.Subplot1.set_ylabel(self.Coordinate + ' [km/s]')
        else:
            self.Subplot1.set_ylabel(self.Coordinate + ' [mas]')
        self.Subplot1.grid()
        self.Subplot1.set_xlim(xlim)
        self._connect_tempo_xlim_sync(self.Subplot1, self.WidgetPlot1, self.WidgetPlot2)


    def Plot2(self):

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot2)

        # Add subplot
        self.Subplot2 = self.WidgetPlot2.Canvas.fig.add_subplot(111, label='Main plot')

        # General plot
        self.general_plot()

        xlim = self.Subplot1.get_xlim() if hasattr(self, 'Subplot1') else (self.TimeWindowMin, self.TimeWindowMax)

        for n in range(self.NbShownOrbits):
            SelectTimesWindow, SelectValuesWindow = self._build_periodic_series(
                self.Selectt[self.nBody][n],
                self.YplotOutput[self.nBody][n],
                self.SelectP[self.nBody][n],
                xlim[0],
                xlim[1],
            )
            if len(SelectTimesWindow) == 0:
                continue
            ReferenceValues = self._reference_values_at(SelectTimesWindow)
            self.Subplot2.plot(SelectTimesWindow, ReferenceValues - SelectValuesWindow, color=self.colorList[self.nBody], linestyle='-', linewidth=0.2, alpha=0.1)

        # Plot input data and residuals
        if len(self.YplotInput) != 0: 
            ReferenceAtInputDates = self._reference_values_at(self.InputDates)
            for k in range(len(self.InputDates)):
                Res = ReferenceAtInputDates[k] - self.YplotInput[k]
                self.Subplot2.errorbar(self.InputDates[k], Res, self.YplotInputErr[k], color='b', linestyle='')

        self.Subplot2.hlines(0, xlim[0], xlim[1], color=self.ReferenceColor, linewidth=1)

        # Plot features
        if self.CoordinateIndex == 3:
            self.Subplot2.set_ylabel(self.Coordinate + ' - ' + self.ReferenceLabel + ' [deg]')
        elif self.CoordinateIndex == 4:
            self.Subplot2.set_ylabel(self.Coordinate + ' - ' + self.ReferenceLabel + ' [km/s]')
        else:
            self.Subplot2.set_ylabel(self.Coordinate + ' - ' + self.ReferenceLabel + ' [mas]')
        self.Subplot2.set_xlabel('Time [MJD]')
        self.Subplot2.set_xlim(xlim)
        self.Subplot2.grid()
        self._connect_tempo_xlim_sync(self.Subplot2, self.WidgetPlot2, self.WidgetPlot1)


class Conv(GeneralToolClass):
    def __init__(self, OutputParams):
        super().__init__('Convergence', 'Convergence of the fit orbit parameters', None, OutputParams, None, None, None, None, None, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylim=True)

    def InitParams(self):
        """Initialize parameters for the Convergence tool."""
        # Equal aspect ratio
        self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)

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
    def __init__(self, OutputParams, BestOrbitParams, LMOrbitParams):
        super().__init__('Histogram', "Histogram of orbital parameters", None, OutputParams, None, None, BestOrbitParams, None, LMOrbitParams, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True)

    def InitParams(self):
        """Initialize parameters for the Histogram tool."""

        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]

        # Equal aspect ratio
        self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Variable', 'Variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0', 'Chi2', 'irel', 'other'])
        if self.NbBodies == 1:
            self.ParamOrbitWidget.ComboParam.removeItem(self.ParamOrbitWidget.ComboParam.count() - 2)  # Remove 'irel' option if only one body
        # if self.V0 is None:
        #     self.ParamOrbitWidget.ComboParam.removeItem(9)  # Remove 'V0' option if no RV data
        # if self.Jitter is None:
        #     self.ParamOrbitWidget.ComboParam.removeItem(10)  # Remove 'Jitter' option if no RV data
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
        self.FormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, m0, Chi2 with optional [n>0] for body number. Use [] or omit the index to use the current orbit.', None)
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

        # Show LM fit
        self.CheckLMFit = CheckBox('LM fit', 'Show the LM fit')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckLMFit)
        self.CheckLMFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)
        self.CheckBestFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show confidence interval
        self.CheckMedian = CheckBox('Median :', 'Show the median and the 1 sigma confidence interval')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckMedian)
        self.CheckMedian.CheckParam.stateChanged.connect(self.refresh_plots)

        # Confidence interval
        self.IntConf = 68
        self.IntConfWidget = SpinBox('Confidence', 'Acceptable level of confidence [%]', self.IntConf, 0, 100, 1)
        self.CheckMedian.Layout.addWidget(self.IntConfWidget)
        self.IntConfWidget.setEnabled(self.CheckMedian.CheckParam.isChecked())
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.IntConfWidget.setEnabled(state))
        self.PercentLbl = QLabel('%')
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
        if self.CheckEqualAspect.CheckParam.isChecked(): self.Subplot.set_box_aspect(1)

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

        # Plot LM fit
        if self.CheckLMFit.CheckParam.isChecked():
            LMParam = self.evaluate_ParamOrbit('self.LM')
            self.Subplot.axvline(LMParam, color='orange')
            self.Subplot.text(LMParam, 0.25 * self.Subplot.get_ylim()[1], s='{}'.format(np.around(LMParam, 3)), color='orange', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='orange'), fontsize=9, horizontalalignment='center', verticalalignment='center', rotation=0)

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
            xlabel_init = r'i$_{rel}$ between '+self.nBodyWidget.ComboParam.currentText()+' and '+self.nBodyRelWidget.ComboParam.currentText()+' [deg]'
        elif self.ParamOrbitWidget.ComboParam.currentText() == 'other':
            xlabel_init = self.FormulaTextEdit.EditParam.text()
        else:
            xlabel_init = self.LabelOf(self.ParamOrbit)+' of '+self.nBodyWidget.ComboParam.currentText()+' '+self.UnitOf(self.ParamOrbit)
        self.Subplot.set_xlabel(xlabel_init)
        self.Subplot.set_ylabel('Count number')

        
class Hist2D(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitParams, LMOrbitParams):
        super().__init__('Histogram 2D', 'Histogram of an orbital parameter as function of another', None, OutputParams, None, None, BestOrbitParams, None, LMOrbitParams, None)

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylim=True)

    def InitParams(self):
        """Initialize parameters for the 2D Histogram tool."""
        # Equal aspect ratio
        self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)
        self.CheckEqualAspect.CheckParam.stateChanged.connect(self.refresh_plots)

        # Abscissa orbit parameters
        self.XParamOrbitWidget = ComboBox('X variable', 'Abscissa variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0', 'Chi2', 'irel', 'other'])
        if self.NbBodies == 1:
            self.XParamOrbitWidget.ComboParam.removeItem(self.XParamOrbitWidget.ComboParam.count() - 2)  # Remove 'irel' option if only one body
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XParamOrbitWidget)
        self.XParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # Orbit number for X parameter
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.XnBodyWidget = ComboBox(None, 'Number of the orbit counting from the center of the system outwards', self.ListBody)
        self.XParamOrbitWidget.Layout.addWidget(self.XnBodyWidget)
        self.XnBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        # if self.NbBodies == 1:
        #     self.XnBodyWidget.setEnabled(False)

        # i relatif between 2 bodies on X axis
        self.XIrelWidget = QWidget()
        self.XIrelLayout = QHBoxLayout()
        self.XnBodyRelWidget = ComboBox('Relative body', 'Relative body to compare with the reference', self.ListBody)
        self.XIrelLayout.addWidget(self.XnBodyRelWidget)
        self.XnBodyRelWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        self.XIrelWidget.setLayout(self.XIrelLayout)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XIrelWidget)
        self.XIrelWidget.setVisible(False)

        # TextEdit for general x formula
        self.XFormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, m0, Chi2 with optional [n>0] for body number. Use [] or omit the index to use the current orbit.', None)
        self.XFormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.XFormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XFormulaTextEdit)
        self.XFormulaTextEdit.EditParam.textChanged.connect(self.reset_plots)

        # Ordinate orbit parameters
        self.YParamOrbitWidget = ComboBox('Y variable', 'Ordinate variable studied in histogram', ['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0', 'Chi2', 'irel', 'other'])
        if self.NbBodies == 1:
            self.YParamOrbitWidget.ComboParam.removeItem(self.YParamOrbitWidget.ComboParam.count() - 2)  # Remove 'irel' option if only one body
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YParamOrbitWidget)
        self.YParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # Orbit number for Y parameter
        self.YnBodyWidget = ComboBox(None, 'Number of the orbit counting from the center of the system outwards', self.ListBody)
        self.YParamOrbitWidget.Layout.addWidget(self.YnBodyWidget)
        self.YnBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        # if self.NbBodies == 1:
        #     self.YnBodyWidget.setEnabled(False)

        # i relatif between 2 bodies on Y axis
        self.YIrelWidget = QWidget()
        self.YIrelLayout = QHBoxLayout()
        self.YnBodyRelWidget = ComboBox('Relative body', 'Relative body to compare with the reference', self.ListBody)
        self.YIrelLayout.addWidget(self.YnBodyRelWidget)
        self.YnBodyRelWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        self.YIrelWidget.setLayout(self.YIrelLayout)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YIrelWidget)
        self.YIrelWidget.setVisible(False)

        # TextEdit for general y formula
        self.YFormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, tp, m, m0, Chi2 with optional [n>0] for body number. Use [] or omit the index to use the current orbit.', None)
        self.YFormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.YFormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YFormulaTextEdit)
        self.YFormulaTextEdit.EditParam.textChanged.connect(self.reset_plots)

        # Connect ComboBox change event
        self.XParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.ToggleXFormulaTextEdit)
        self.YParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.ToggleYFormulaTextEdit)
        self.XParamOrbitWidget.ComboParam.currentIndexChanged.connect(lambda: self.XIrelWidget.setVisible(self.XParamOrbitWidget.ComboParam.currentText() == 'irel'))
        self.YParamOrbitWidget.ComboParam.currentIndexChanged.connect(lambda: self.YIrelWidget.setVisible(self.YParamOrbitWidget.ComboParam.currentText() == 'irel'))

        # Histogram binning
        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show LM fit
        self.CheckLMFit = CheckBox('LM fit', 'Show the Levenberg-Marquardt fit')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckLMFit)
        self.CheckLMFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)
        self.CheckBestFit.CheckParam.stateChanged.connect(self.refresh_plots)

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
    
    def evaluate_ParamOrbit(self, prefixe, param_widget, formula_text_edit, nbody, nbody_rel):
        """Evaluate the parameter orbit based on the current widget values."""
        parameter = param_widget.ComboParam.currentText()
        if parameter == 'irel' or parameter == 'other':
            if parameter == 'irel':
                if nbody_rel == nbody:
                    print('nBodyRel can not be the same as nBody')
                    return None
                else:
                    formula = f'arccos(cos(i[{nbody+1}])*cos(i[{nbody_rel+1}])+cos(W[{nbody+1}]-W[{nbody_rel+1}])*sin(i[{nbody+1}])*sin(i[{nbody_rel+1}]))'
            elif parameter == 'other':
                formula = formula_text_edit.EditParam.text()
            return self.evaluate_formula(formula, prefixe, nbody)
        else:
            parameter = param_widget.ComboParam.currentText()
            # print('not irel or other')
            return eval(f'{prefixe}{parameter}')[nbody]

    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        self.XnBody = int(self.XnBodyWidget.ComboParam.currentText()) - 1
        self.XnBodyRel = int(self.XnBodyRelWidget.ComboParam.currentText()) - 1
        self.YnBody = int(self.YnBodyWidget.ComboParam.currentText()) - 1
        self.YnBodyRel = int(self.YnBodyRelWidget.ComboParam.currentText()) - 1
        self.XParamOrbit = self.XParamOrbitWidget.ComboParam.currentText()
        self.YParamOrbit = self.YParamOrbitWidget.ComboParam.currentText()
        self.EvalXParamOrbit = self.evaluate_ParamOrbit('self.', self.XParamOrbitWidget, self.XFormulaTextEdit, self.XnBody, self.XnBodyRel)
        self.EvalYParamOrbit = self.evaluate_ParamOrbit('self.', self.YParamOrbitWidget, self.YFormulaTextEdit, self.YnBody, self.YnBodyRel)
        self.NbBins = self.NbBinsWidget.SpinParam.value()

    def Plot(self):
        """Plot the 2D histogram based on the selected parameters."""

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)
        if self.CheckEqualAspect.CheckParam.isChecked(): self.Subplot.set_box_aspect(1)

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

        # LM fit
        if self.CheckLMFit.CheckParam.isChecked():
            LMXParam = self.evaluate_ParamOrbit('self.LM', self.XParamOrbitWidget, self.XFormulaTextEdit, self.XnBody, self.XnBodyRel)
            LMYParam = self.evaluate_ParamOrbit('self.LM', self.YParamOrbitWidget, self.YFormulaTextEdit, self.YnBody, self.YnBodyRel)
            self.Subplot.plot(LMXParam, LMYParam, color='orange', marker='x', markersize=8, markeredgewidth=2)

        # Best fit
        if self.CheckBestFit.CheckParam.isChecked():
            BestXParam = self.evaluate_ParamOrbit('self.Best', self.XParamOrbitWidget, self.XFormulaTextEdit, self.XnBody, self.XnBodyRel)
            BestYParam = self.evaluate_ParamOrbit('self.Best', self.YParamOrbitWidget, self.YFormulaTextEdit, self.YnBody, self.YnBodyRel)
            self.Subplot.plot(BestXParam, BestYParam, color='red', marker='x', markersize=8, markeredgewidth=2)

        # Plot features
        # about X axis
        if self.XParamOrbitWidget.ComboParam.currentIndex() == self.XParamOrbitWidget.ComboParam.count() - 1: # 'other' selected
            self.Subplot.set_xlabel(self.XFormulaTextEdit.EditParam.text())
        elif self.XParamOrbitWidget.ComboParam.currentIndex() == self.XParamOrbitWidget.ComboParam.count() - 2: # 'irel' selected
            self.Subplot.set_xlabel(r'i$_{rel}$ between '+self.XnBodyWidget.ComboParam.currentText()+' and '+self.XnBodyRelWidget.ComboParam.currentText()+' [deg]')
        else:
            self.Subplot.set_xlabel(self.LabelOf(self.XParamOrbit)+' of '+self.XnBodyWidget.ComboParam.currentText()+' '+self.UnitOf(self.XParamOrbit))
        # about Y axis
        if self.YParamOrbitWidget.ComboParam.currentIndex() == self.YParamOrbitWidget.ComboParam.count() - 1: # 'other' selected
            self.Subplot.set_ylabel(self.YFormulaTextEdit.EditParam.text())
        elif self.YParamOrbitWidget.ComboParam.currentIndex() == self.YParamOrbitWidget.ComboParam.count() - 2: # 'irel' selected
            self.Subplot.set_ylabel(r'i$_{rel}$ between '+self.YnBodyWidget.ComboParam.currentText()+' and '+self.YnBodyRelWidget.ComboParam.currentText()+' [deg]')
        else:
            self.Subplot.set_ylabel(self.LabelOf(self.YParamOrbit)+' of '+self.YnBodyWidget.ComboParam.currentText()+' '+self.UnitOf(self.YParamOrbit))


class Corner(GeneralToolClass):
    def __init__(self, SelectOrbitsParams, BestOrbitParams, LMOrbitParams):
        super().__init__('Corner', 'Corner plot of parameters', None, None, SelectOrbitsParams, None, BestOrbitParams, None, LMOrbitParams, None)
        self.MaskCondition = None
        self._corner_sync_in_progress = False
        self._corner_limit_sync_pending = False
        self._corner_axis_sync_ready = False
        self._corner_layout_refresh_pending = False
        self._corner_layout_refresh_in_progress = False
        self._corner_layout_reference = None
        self._corner_param_limits = {}
        self._corner_param_names = []
        self._corner_default_limits = {}
        self.corner_label_fontsize = 11
        self.corner_tick_fontsize = 9
        self.corner_tick_fractions = np.array([0.10, 0.30, 0.50, 0.70, 0.90])
        self.corner_min_significant_digits = 3
        self.corner_bottom_tick_rotation = 35
        self.corner_tick_pad = 2.5
        self.corner_label_min_gap_px = 8.0
        self.corner_left_label_extra_gap_px = 8.0

        # Parameters initialisation
        self.InitParams()

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot)
        self.WidgetPlot.Canvas.mpl_connect('button_release_event', self.on_corner_button_release)
        self.WidgetPlot.Canvas.mpl_connect('resize_event', self.on_corner_resize_event)
        self.WidgetPlot.Toolbar._actions['home'].triggered.connect(self.queue_corner_limit_reset)

    def queue_corner_layout_refresh(self):
        if self._corner_layout_refresh_pending:
            return
        self._corner_layout_refresh_pending = True
        QTimer.singleShot(0, self.refresh_corner_layout)

    def refresh_corner_layout(self):
        self._corner_layout_refresh_pending = False
        if self._corner_layout_refresh_in_progress or not self._corner_axis_sync_ready:
            return
        if len(self._corner_param_names) == 0 or len(self.WidgetPlot.Canvas.fig.axes) == 0:
            return

        self._corner_layout_refresh_in_progress = True
        try:
            self.apply_corner_square_layout()
            self.apply_corner_uniform_label_padding()
        finally:
            self._corner_layout_refresh_in_progress = False

    def on_corner_resize_event(self, event):
        self.queue_corner_layout_refresh()

    @staticmethod
    def normalize_corner_limit(limit):
        if limit is None or len(limit) != 2:
            return None
        lower = float(min(limit))
        upper = float(max(limit))
        if not np.isfinite(lower) or not np.isfinite(upper):
            return None
        if np.isclose(lower, upper, rtol=1e-9, atol=1e-12):
            return None
        return (lower, upper)

    @staticmethod
    def corner_limits_equal(limit1, limit2):
        if limit1 is None and limit2 is None:
            return True
        if limit1 is None or limit2 is None:
            return False
        return np.allclose(limit1, limit2, rtol=1e-9, atol=1e-9)

    @staticmethod
    def finite_range(values, fallback=None):
        values = np.asarray(values)
        values = values[np.isfinite(values)]
        if values.size == 0:
            return fallback
        lower = float(np.min(values))
        upper = float(np.max(values))
        if np.isclose(lower, upper):
            return fallback
        return (lower, upper)

    @staticmethod
    def decimals_for_significant_digits(value, minimum_significant_digits):
        value = abs(float(value))

        if not np.isfinite(value):
            return 0

        if np.isclose(value, 0.0):
            return max(0, minimum_significant_digits - 1)

        magnitude = int(np.floor(np.log10(value)))
        return max(0, minimum_significant_digits - 1 - magnitude)

    def axis_decimals_for_ticks(self, ticks):
        ticks = np.asarray(ticks, dtype=float)
        finite_ticks = ticks[np.isfinite(ticks)]

        if finite_ticks.size == 0:
            return 0

        max_abs_tick = np.max(np.abs(finite_ticks))
        return self.decimals_for_significant_digits(max_abs_tick, self.corner_min_significant_digits)

    @staticmethod
    def corner_tick_labels_are_distinct(labels):
        previous_label = None
        for label in labels:
            if previous_label is not None and label == previous_label:
                return False
            previous_label = label
        return True

    def format_corner_ticks(self, ticks, decimals_override=None):
        axis_decimals = self.axis_decimals_for_ticks(ticks) if decimals_override is None else int(decimals_override)
        max_decimals = max(axis_decimals, 12)

        for decimals in range(axis_decimals, max_decimals + 1):
            labels = [f'{tick:.{decimals}f}' for tick in ticks]
            if self.corner_tick_labels_are_distinct(labels):
                return labels

        return [f'{tick:.{max_decimals}f}' for tick in ticks]

    def set_corner_fraction_ticks(self, axis, axis_name, decimals_override=None):
        fractions = np.asarray(self.corner_tick_fractions, dtype=float)

        if axis_name == 'x':
            axis_min, axis_max = axis.get_xlim()
            set_ticks = axis.set_xticks
            set_labels = axis.set_xticklabels
        elif axis_name == 'y':
            axis_min, axis_max = axis.get_ylim()
            set_ticks = axis.set_yticks
            set_labels = axis.set_yticklabels
        else:
            raise ValueError("axis_name must be 'x' or 'y'.")

        if not np.isfinite(axis_min) or not np.isfinite(axis_max) or np.isclose(axis_min, axis_max):
            return

        ticks = axis_min + (axis_max - axis_min) * fractions
        set_ticks(ticks)
        set_labels(self.format_corner_ticks(ticks, decimals_override=decimals_override))

    @staticmethod
    def normalize_histogram(histogram_values):
        peak = np.max(histogram_values) if histogram_values.size != 0 else 0.0
        return histogram_values / peak if peak > 0 else histogram_values

    def compute_corner_normalized_histogram(self, values, value_range, weights=None):
        lower, upper = value_range
        edges = np.linspace(lower, upper, self.NbBins + 1)
        if not np.isfinite(lower) or not np.isfinite(upper) or np.isclose(lower, upper):
            return np.zeros(self.NbBins), edges

        values = np.asarray(values, dtype=float)
        finite_mask = np.isfinite(values)

        histogram_weights = None
        if weights is not None:
            histogram_weights = np.asarray(weights, dtype=float)
            finite_mask &= np.isfinite(histogram_weights)

        values = values[finite_mask]
        if histogram_weights is not None:
            histogram_weights = histogram_weights[finite_mask]

        if values.size == 0:
            return np.zeros(self.NbBins), edges

        if histogram_weights is not None and np.allclose(histogram_weights, 0.0):
            return np.zeros(self.NbBins), edges

        histogram, edges = np.histogram(
            values,
            bins=self.NbBins,
            range=value_range,
            density=True,
            weights=histogram_weights,
        )
        histogram = np.nan_to_num(histogram, nan=0.0, posinf=0.0, neginf=0.0)
        return self.normalize_histogram(histogram), edges

    def redraw_corner_diagonal_histograms(self, grid, data, corner_ranges, main_weights, overlay_weights=None):
        for axis in grid.get_axes():
            axis_position = self.get_corner_axis_position(axis)
            if axis_position is None:
                continue

            row, col = axis_position
            if row != col:
                continue

            value_range = corner_ranges[col]
            if value_range is None:
                continue

            axis.clear()
            base_histogram, edges = self.compute_corner_normalized_histogram(
                data[:, col],
                value_range,
                weights=main_weights,
            )
            axis.stairs(base_histogram, edges, color='black', linewidth=1.5)
            histogram_peak = float(np.max(base_histogram)) if base_histogram.size != 0 else 0.0

            if overlay_weights is not None and np.any(np.asarray(overlay_weights, dtype=float) > 0):
                overlay_histogram, _ = self.compute_corner_normalized_histogram(
                    data[:, col],
                    value_range,
                    weights=overlay_weights,
                )
                axis.stairs(overlay_histogram, edges, color='C0', linewidth=1.7)
                if overlay_histogram.size != 0:
                    histogram_peak = max(histogram_peak, float(np.max(overlay_histogram)))

            axis.set_xlim(value_range)
            if histogram_peak > 0:
                axis.set_ylim(0.0, 1.05 * histogram_peak)
            else:
                axis.set_ylim(0.0, 1.0)

    def style_corner_axis(self, axis, row, col, nb_params, data_labels):
        axis.set_box_aspect(1)
        axis.tick_params(axis='both', labelsize=self.corner_tick_fontsize, pad=self.corner_tick_pad)
        axis.xaxis.labelpad = 0.0
        axis.yaxis.labelpad = 0.0

        if row < col:
            axis.axis('off')
            return

        self.set_corner_fraction_ticks(axis, 'x')

        if row == col:
            self.set_corner_fraction_ticks(axis, 'y', decimals_override=1 if col == 0 else None)
            if col == 0:
                axis.set_ylabel('Normalized density', fontsize=self.corner_label_fontsize)
                axis.tick_params(axis='y', labelleft=True)
                for label in axis.get_yticklabels():
                    label.set_rotation(0)
                    label.set_ha('right')
                    label.set_rotation_mode('anchor')
            else:
                axis.set_ylabel('')
                axis.tick_params(axis='y', labelleft=False)
        else:
            self.set_corner_fraction_ticks(axis, 'y')
            if col == 0:
                axis.set_ylabel(data_labels[row], fontsize=self.corner_label_fontsize)
                axis.tick_params(axis='y', labelleft=True)
                for label in axis.get_yticklabels():
                    label.set_rotation(0)
                    label.set_ha('right')
                    label.set_rotation_mode('anchor')
            else:
                axis.set_ylabel('')
                axis.tick_params(axis='y', labelleft=False)

        if row == nb_params - 1:
            axis.set_xlabel(data_labels[col], fontsize=self.corner_label_fontsize)
            axis.tick_params(axis='x', labelbottom=True, labelrotation=self.corner_bottom_tick_rotation)
            for label in axis.get_xticklabels():
                label.set_ha('right')
                label.set_rotation_mode('anchor')
        else:
            axis.set_xlabel('')
            axis.tick_params(axis='x', labelbottom=False)

    def active_corner_limits(self, data_names):
        return {name: self._corner_param_limits[name] for name in data_names if name in self._corner_param_limits}

    def build_corner_interactive_mask(self, data, data_names, active_limits):
        interactive_mask = np.ones(len(data), dtype=bool)
        for index, name in enumerate(data_names):
            if name not in active_limits:
                continue
            lower, upper = active_limits[name]
            interactive_mask &= (data[:, index] >= lower) & (data[:, index] <= upper)
        return interactive_mask

    def build_corner_ranges(self, data, data_names, interactive_mask, active_limits):
        if np.any(interactive_mask):
            filtered_data = data[interactive_mask]
        else:
            filtered_data = data

        ranges = []
        for index, name in enumerate(data_names):
            default_range = self._corner_default_limits.get(name)
            if name in active_limits:
                ranges.append(active_limits[name])
                continue
            filtered_range = self.finite_range(filtered_data[:, index], default_range)
            ranges.append(filtered_range)
        return ranges

    def set_corner_param_limit(self, name, limit):
        normalized_limit = self.normalize_corner_limit(limit)
        if normalized_limit is None:
            return False
        default_limit = self._corner_default_limits.get(name)

        if default_limit is not None and self.corner_limits_equal(normalized_limit, default_limit):
            if name in self._corner_param_limits:
                del self._corner_param_limits[name]
                return True
            return False

        current_limit = self._corner_param_limits.get(name)
        if self.corner_limits_equal(current_limit, normalized_limit):
            return False

        self._corner_param_limits[name] = normalized_limit
        return True

    def queue_corner_limit_reset(self):
        QTimer.singleShot(0, self.reset_corner_limits)

    def reset_corner_limits(self):
        if len(self._corner_param_limits) == 0:
            return
        self._corner_param_limits = {}
        self.WidgetPlot.refresh_plot()

    def get_corner_axis_position(self, axis):
        nb_params = len(self._corner_param_names)
        if nb_params == 0 or axis is None:
            return None

        subplot_spec = axis.get_subplotspec() if hasattr(axis, 'get_subplotspec') else None
        if subplot_spec is None:
            return None

        row = subplot_spec.rowspan.start
        col = subplot_spec.colspan.start
        if row >= nb_params or col >= nb_params:
            return None
        return row, col

    def get_corner_axes_by_position(self):
        axes_by_position = {}
        for axis in self.WidgetPlot.Canvas.fig.axes:
            axis_position = self.get_corner_axis_position(axis)
            if axis_position is None:
                continue
            axes_by_position[axis_position] = axis
        return axes_by_position

    def capture_corner_layout_reference(self):
        nb_params = len(self._corner_param_names)
        if nb_params == 0:
            self._corner_layout_reference = None
            return

        axes_by_position = self.get_corner_axes_by_position()
        if len(axes_by_position) != nb_params * nb_params:
            self._corner_layout_reference = None
            return

        axis_boxes = [axis.get_position() for axis in axes_by_position.values()]
        left = min(box.x0 for box in axis_boxes)
        right = max(box.x1 for box in axis_boxes)
        bottom = min(box.y0 for box in axis_boxes)
        top = max(box.y1 for box in axis_boxes)

        fig_width, fig_height = self.WidgetPlot.Canvas.fig.get_size_inches()
        if fig_width <= 0 or fig_height <= 0:
            self._corner_layout_reference = None
            return

        horizontal_gaps_in = []
        for row in range(nb_params):
            for col in range(nb_params - 1):
                left_axis = axes_by_position[(row, col)].get_position()
                right_axis = axes_by_position[(row, col + 1)].get_position()
                gap_in = (right_axis.x0 - left_axis.x1) * fig_width
                if np.isfinite(gap_in):
                    horizontal_gaps_in.append(max(float(gap_in), 0.0))

        vertical_gaps_in = []
        for row in range(nb_params - 1):
            for col in range(nb_params):
                top_axis = axes_by_position[(row, col)].get_position()
                bottom_axis = axes_by_position[(row + 1, col)].get_position()
                gap_in = (top_axis.y0 - bottom_axis.y1) * fig_height
                if np.isfinite(gap_in):
                    vertical_gaps_in.append(max(float(gap_in), 0.0))

        self._corner_layout_reference = {
            'nb_params': nb_params,
            'bounds': (left, right, bottom, top),
            'gap_x_in': float(np.mean(horizontal_gaps_in)) if len(horizontal_gaps_in) != 0 else 0.0,
            'gap_y_in': float(np.mean(vertical_gaps_in)) if len(vertical_gaps_in) != 0 else 0.0,
        }

    def apply_corner_uniform_label_padding(self):
        nb_params = len(self._corner_param_names)
        if nb_params == 0:
            return

        axes_by_position = self.get_corner_axes_by_position()
        if len(axes_by_position) != nb_params * nb_params:
            return

        bottom_row_axes = [axes_by_position[(nb_params - 1, col)] for col in range(nb_params) if (nb_params - 1, col) in axes_by_position]
        left_column_axes = [axes_by_position[(row, 0)] for row in range(nb_params) if (row, 0) in axes_by_position]

        for axis in bottom_row_axes:
            if axis.get_xlabel():
                axis.xaxis.label.set_transform(axis.transAxes)
                axis.xaxis.label.set_horizontalalignment('center')
                axis.xaxis.label.set_verticalalignment('center')
                axis.xaxis.label.set_rotation_mode('anchor')
                axis.xaxis.set_label_coords(0.5, 0.0)

        for axis in left_column_axes:
            if axis.get_ylabel():
                axis.yaxis.label.set_transform(axis.transAxes)
                axis.yaxis.label.set_horizontalalignment('center')
                axis.yaxis.label.set_verticalalignment('center')
                axis.yaxis.label.set_rotation_mode('anchor')
                axis.yaxis.set_label_coords(0.0, 0.5)

        self.WidgetPlot.Canvas.draw()
        renderer = self.WidgetPlot.Canvas.fig.canvas.get_renderer()

        common_x_label_center_px = None
        for axis in bottom_row_axes:
            if not axis.get_xlabel():
                continue
            tick_bboxes = [
                tick.get_window_extent(renderer)
                for tick in axis.get_xticklabels()
                if tick.get_visible() and len(tick.get_text()) != 0
            ]
            if len(tick_bboxes) == 0:
                continue
            tick_bbox = Bbox.union(tick_bboxes)
            label_bbox = axis.xaxis.label.get_window_extent(renderer)
            candidate_center_px = tick_bbox.y0 - self.corner_label_min_gap_px - 0.5 * label_bbox.height
            common_x_label_center_px = candidate_center_px if common_x_label_center_px is None else min(common_x_label_center_px, candidate_center_px)

        common_y_label_center_px = None
        for axis in left_column_axes:
            if not axis.get_ylabel():
                continue
            tick_bboxes = [
                tick.get_window_extent(renderer)
                for tick in axis.get_yticklabels()
                if tick.get_visible() and len(tick.get_text()) != 0
            ]
            if len(tick_bboxes) == 0:
                continue
            tick_bbox = Bbox.union(tick_bboxes)
            label_bbox = axis.yaxis.label.get_window_extent(renderer)
            candidate_center_px = tick_bbox.x0 - self.corner_label_min_gap_px - self.corner_left_label_extra_gap_px - 0.5 * label_bbox.width
            common_y_label_center_px = candidate_center_px if common_y_label_center_px is None else min(common_y_label_center_px, candidate_center_px)

        for axis in bottom_row_axes:
            if axis.get_xlabel():
                axis_bbox = axis.get_window_extent(renderer)
                if common_x_label_center_px is not None and axis_bbox.height > 0:
                    label_y = (common_x_label_center_px - axis_bbox.y0) / axis_bbox.height
                    axis.xaxis.set_label_coords(0.5, label_y)

        for axis in left_column_axes:
            if axis.get_ylabel():
                axis_bbox = axis.get_window_extent(renderer)
                if common_y_label_center_px is not None and axis_bbox.width > 0:
                    label_x = (common_y_label_center_px - axis_bbox.x0) / axis_bbox.width
                    axis.yaxis.set_label_coords(label_x, 0.5)

        self.WidgetPlot.Canvas.draw()

        self.WidgetPlot.Canvas.draw()

    def apply_corner_square_layout(self):
        nb_params = len(self._corner_param_names)
        if nb_params == 0:
            return

        axes_by_position = self.get_corner_axes_by_position()

        if len(axes_by_position) != nb_params * nb_params:
            return

        layout_reference = self._corner_layout_reference
        if layout_reference is None or layout_reference.get('nb_params') != nb_params:
            self.capture_corner_layout_reference()
            layout_reference = self._corner_layout_reference
        if layout_reference is None:
            return

        left, right, bottom, top = layout_reference['bounds']
        available_width = right - left
        available_height = top - bottom
        if available_width <= 0 or available_height <= 0:
            return

        fig_width, fig_height = self.WidgetPlot.Canvas.fig.get_size_inches()
        if fig_width <= 0 or fig_height <= 0:
            return

        available_width_in = available_width * fig_width
        available_height_in = available_height * fig_height
        gap_x_in = layout_reference.get('gap_x_in', 0.0)
        gap_y_in = layout_reference.get('gap_y_in', 0.0)
        box_size_in = min(
            max((available_width_in - (nb_params - 1) * gap_x_in) / nb_params, 0.0),
            max((available_height_in - (nb_params - 1) * gap_y_in) / nb_params, 0.0),
        )
        if box_size_in <= 0:
            gap_x_in = 0.0
            gap_y_in = 0.0
            box_size_in = min(available_width_in / nb_params, available_height_in / nb_params)
        if box_size_in <= 0:
            return

        box_width = box_size_in / fig_width
        box_height = box_size_in / fig_height
        horizontal_gap = gap_x_in / fig_width
        vertical_gap = gap_y_in / fig_height

        total_width = nb_params * box_width + (nb_params - 1) * horizontal_gap
        total_height = nb_params * box_height + (nb_params - 1) * vertical_gap
        x_offset = left + 0.5 * max(available_width - total_width, 0.0)
        y_offset = bottom + 0.5 * max(available_height - total_height, 0.0)

        for (row, col), axis in axes_by_position.items():
            x0 = x_offset + col * (box_width + horizontal_gap)
            y0 = y_offset + (nb_params - 1 - row) * (box_height + vertical_gap)
            axis.set_position([x0, y0, box_width, box_height])
            axis.set_box_aspect(1)

    def connect_corner_axis_limit_sync(self):
        for axis in self.WidgetPlot.Canvas.fig.axes:
            axis.callbacks.connect('xlim_changed', self.on_corner_axis_limit_changed)
            axis.callbacks.connect('ylim_changed', self.on_corner_axis_limit_changed)

    def on_corner_axis_limit_changed(self, axis):
        if not self._corner_axis_sync_ready or self._corner_sync_in_progress:
            return
        self.queue_corner_limit_sync(axis)

    def queue_corner_limit_sync(self, axis):
        if self._corner_limit_sync_pending:
            return
        self._corner_limit_sync_pending = True
        QTimer.singleShot(0, lambda current_axis=axis: self.apply_corner_axis_limits(current_axis))

    def apply_corner_axis_limits(self, axis):
        self._corner_limit_sync_pending = False
        if self._corner_sync_in_progress or not self._corner_axis_sync_ready or axis is None:
            return

        axis_position = self.get_corner_axis_position(axis)
        if axis_position is None:
            return

        row, col = axis_position
        if row < col:
            return

        changed = self.set_corner_param_limit(self._corner_param_names[col], axis.get_xlim())
        if row > col:
            changed = self.set_corner_param_limit(self._corner_param_names[row], axis.get_ylim()) or changed

        if not changed:
            return

        self._corner_sync_in_progress = True
        try:
            self.WidgetPlot.refresh_plot()
        finally:
            self._corner_sync_in_progress = False

    def on_corner_button_release(self, event):
        if self._corner_sync_in_progress or not self._corner_axis_sync_ready or event.inaxes is None:
            return
        if len(self._corner_param_names) == 0:
            return
        if event.inaxes not in self.WidgetPlot.Canvas.fig.axes:
            return
        self.queue_corner_limit_sync(event.inaxes)

    def InitParams(self):
        """Initialize parameters for the Corner plot tool."""

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies)]
        self.nBodyWidget = ComboBox("Orbit number", 'Number of the studied orbit counting from the center of the system outwards', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        # if self.NbBodies == 1:
        #     self.nBodyWidget.setEnabled(False)
        self.nBodyWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # Parameters to include in the corner plot
        self.ParamContainer = QWidget()
        self.ParamLayoutV = QVBoxLayout()

        self.ParamLayoutV.addWidget(QLabel('Variables :'))

        self.ParamCheckGroupBox = QGroupBox()
        self.ParamCheckLayout = QGridLayout()

        self.OrbitParams = [
            ['a', 'e', 'i', 'w', 'W', 'tp', 'm', 'm0', 'P'],
            ['Semi-major axis', 'Eccentricity', 'Inclination', 'Argument of periastron', 'Longitude of ascending node', 'Periastron time passage', 'Body mass', 'Central body mass', 'Period']
        ]
        self.WidgetOrbitParams = []

        for i in range(len(self.OrbitParams[0])):
            ParamCheckBox = CheckBox(self.OrbitParams[0][i], self.OrbitParams[1][i])
            ParamCheckBox.CheckParam.stateChanged.connect(self.refresh_plots)

            # Keep the default selection independent from the display order.
            if self.OrbitParams[0][i] in ('P', 'a', 'e'):
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

        # Show LM fit
        self.CheckLMFit = CheckBox('LM fit', 'Show the Levenberg-Marquardt fit')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckLMFit)
        self.CheckLMFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)
        self.CheckBestFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Mask overlay
        self.CheckMask = CheckBox('Mask', 'Show the masked subset overlay')
        self.CheckMask.CheckParam.setChecked(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckMask)
        self.MaskConditionWidget = LineEdit(
            '   Condition',
            'Only variables P, a, e, i, w, W, tp, m, m0, Chi2 with optional [n>0] for body number. Use [] or omit the index to use the current orbit.',
            '',
        )
        self.MaskConditionWidget.EditParam.setPlaceholderText('Examples: P<3, P[]<3, P[1]<3 and Chi2<2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.MaskConditionWidget)
        self.MaskConditionWidget.EditParam.textChanged.connect(self.refresh_plots)
        self.CheckMask.CheckParam.stateChanged.connect(self.ToggleMaskConditionWidget)
        self.ToggleMaskConditionWidget(self.CheckMask.CheckParam.isChecked())
        self.CheckMask.CheckParam.stateChanged.connect(self.refresh_plots)

        # Cosmetic options
        self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter(Title='Options :'))

        self.CheckShortLabels = CheckBox('Short labels', 'Show short labels')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckShortLabels)
        self.CheckShortLabels.CheckParam.stateChanged.connect(self.refresh_plots)

        self.ShortLabelIndexWidget = LineEdit('Short label index', 'Optional index appended to corner short labels', '')
        self.ShortLabelIndexWidget.EditParam.setPlaceholderText('Examples: in, out, b, 1')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ShortLabelIndexWidget)
        self.ShortLabelIndexWidget.EditParam.textChanged.connect(self.refresh_plots)
        self.ShortLabelIndexWidget.setVisible(self.CheckShortLabels.CheckParam.isChecked())
        self.CheckShortLabels.CheckParam.stateChanged.connect(lambda state: self.ShortLabelIndexWidget.setVisible(bool(state)))

        self.CheckContour = CheckBox('Contours', 'Show contours')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckContour)
        self.CheckContour.CheckParam.stateChanged.connect(self.refresh_plots)

        self.CheckDensity = CheckBox('Density', 'Show density representation')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckDensity)
        self.CheckDensity.CheckParam.stateChanged.connect(self.refresh_plots)


    def ToggleMaskConditionWidget(self, state):
        """Show the mask condition field only when the mask overlay is enabled."""
        self.MaskConditionWidget.setVisible(bool(state))


    def build_selected_mask(self):
        """Build the overlay mask from the general analysis condition."""
        if self.MaskCondition is None:
            return None

        filter_variables = {
            'a': self.Selecta,
            'P': self.SelectP,
            'e': self.Selecte,
            'w': self.Selectw,
            'i': self.Selecti,
            'W': self.SelectW,
            'tp': self.Selecttp,
            'm': self.Selectm,
            'V0': self.SelectV0,
            'Jitter': self.SelectJitter,
            'm0': self.Selectm0,
            'Chi2': self.SelectChi2,
            'Map': self.SelectMap,
        }
        filter_variables = {name: values for name, values in filter_variables.items() if values is not None}

        translated_condition = self.replace_params_in_condition(self.MaskCondition, '', self.nBody)

        try:
            return BuildOrbitFilterMask(translated_condition, filter_variables, self.NbSelectOrbits, normalize_indices=False)
        except Exception as exc:
            print(f'Error evaluating corner mask condition: {exc}')
            return None


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
        mask_condition = self.MaskConditionWidget.EditParam.text().strip()
        self.MaskCondition = mask_condition if len(mask_condition) != 0 else None
        short_label_index = self.ShortLabelIndexWidget.EditParam.text().strip()
        self.ShortLabelIndex = short_label_index if len(short_label_index) != 0 else None

    def Plot(self):
        """Plot the corner plot based on the selected parameters."""

        self._corner_axis_sync_ready = False

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # No subplot because corner plots are standalone

        # Collect data for the corner plot
        Data = []
        DataNames = []
        DataLabels = []
        DataDefaultRanges = []
        for x in self.WidgetOrbitParams:
            if not x.CheckParam.isChecked() or not x.CheckParam.isEnabled():
                continue

            param_name = x.CheckParam.text()
            param_values = np.asarray(eval(f'self.Select{param_name}')[self.nBody])
            default_range = self.finite_range(param_values)
            if default_range is None:
                continue

            Data.append(param_values)
            DataNames.append(param_name)
            DataDefaultRanges.append(default_range)
            if self.CheckShortLabels.CheckParam.isChecked():
                DataLabels.append(f'{self.ShortLabelOf(param_name, self.ShortLabelIndex)} {self.UnitOf(param_name)}')
            else:
                DataLabels.append(f'{self.LabelOf(param_name)} {self.UnitOf(param_name)}')
        Data = np.array(Data).T

        # Check if there is data to plot
        if len(Data) == 0:
            print('No parameters selected or all parameters have no variance.')
            self.WidgetPlot.Canvas.fig.canvas.draw()
            return

        self._corner_param_names = list(DataNames)
        self._corner_default_limits = {
            DataNames[index]: DataDefaultRanges[index]
            for index in range(len(DataNames))
        }
        self._corner_param_limits = {
            name: limit
            for name, limit in self._corner_param_limits.items()
            if name in self._corner_default_limits and not self.corner_limits_equal(limit, self._corner_default_limits[name])
        }

        active_limits = self.active_corner_limits(DataNames)
        interactive_mask = self.build_corner_interactive_mask(Data, DataNames, active_limits)
        corner_ranges = self.build_corner_ranges(Data, DataNames, interactive_mask, active_limits)
        has_interactive_limits = len(active_limits) != 0
        has_interactive_selection = has_interactive_limits and np.any(interactive_mask)

        if has_interactive_limits and not has_interactive_selection:
            self._corner_param_limits = {}
            active_limits = {}
            interactive_mask = np.ones(len(Data), dtype=bool)
            corner_ranges = self.build_corner_ranges(Data, DataNames, interactive_mask, active_limits)
            has_interactive_limits = False
            has_interactive_selection = False

        mask = self.build_selected_mask()
        weights = mask.astype(float) if mask is not None else None
        if weights is not None and has_interactive_limits:
            weights = weights * interactive_mask.astype(float)

        if has_interactive_limits:
            main_weights = interactive_mask.astype(float)
        else:
            main_weights = None

        show_contours = self.CheckContour.CheckParam.isChecked()
        show_density = self.CheckDensity.CheckParam.isChecked()
        show_mask = self.CheckMask.CheckParam.isChecked() and weights is not None and np.any(weights)

        if has_interactive_limits and not np.any(main_weights):
            show_contours = False
            show_density = False

        # Create global corner plot without automatic datapoints. Points are
        # added afterwards so they remain visible above density/contour layers.
        grid = corner.corner(
            Data,
            labels=DataLabels,
            bins=self.NbBins,
            fig=self.WidgetPlot.Canvas.fig,
            range=corner_ranges,
            plot_contours=show_contours,
            fill_contours=show_contours and not show_density,
            plot_density=show_density,
            plot_datapoints=False,
            weights=main_weights,
            contour_kwargs={'linewidths': 1.0},
        )

        # Create weighted corner plot on top of the global one.
        if show_mask:
            grid = corner.corner(
                Data,
                labels=DataLabels,
                bins=self.NbBins,
                fig=self.WidgetPlot.Canvas.fig,
                range=corner_ranges,
                plot_contours=show_contours,
                fill_contours=show_contours and not show_density,
                plot_density=show_density,
                plot_datapoints=False,
                weights=weights,
                color='C0',
                contour_kwargs={'linewidths': 1.3},
            )

        # When contours or density are enabled, keep points very faint so the
        # structure remains readable.
        point_alpha = 0.05 if (show_contours or show_density) else 0.15
        point_size = 0.8 if (show_contours or show_density) else 2
        masked_point_alpha = 0.08 if (show_contours or show_density) else 0.35
        masked_point_size = 1.0 if (show_contours or show_density) else 2

        point_data = Data[interactive_mask] if has_interactive_selection else Data
        corner.overplot_points(self.WidgetPlot.Canvas.fig, point_data, marker='.', color='black', alpha=point_alpha, ms=point_size)
        if show_mask:
            corner.overplot_points(self.WidgetPlot.Canvas.fig, Data[weights > 0], marker='.', color='C0', alpha=masked_point_alpha, ms=masked_point_size)

        self.redraw_corner_diagonal_histograms(
            grid,
            Data,
            corner_ranges,
            main_weights,
            weights if show_mask else None,
        )

        nb_params = len(DataLabels)
        for ax in grid.get_axes():
            axis_position = self.get_corner_axis_position(ax)
            if axis_position is None:
                continue

            row, col = axis_position
            self.style_corner_axis(ax, row, col, nb_params, DataLabels)

            # LM fit
            if self.CheckLMFit.CheckParam.isChecked():
                if row == col:
                    BestParam = eval(f'self.LM{DataNames[col]}')[self.nBody]
                    ax.axvline(BestParam, color='orange', linestyle='-', linewidth=0.75)
                elif row > col:
                    BestXParam = eval(f'self.LM{DataNames[col]}')[self.nBody]
                    BestYParam = eval(f'self.LM{DataNames[row]}')[self.nBody]
                    ax.plot(BestXParam, BestYParam, color='orange', marker='x')

            # Best fit
            if self.CheckBestFit.CheckParam.isChecked():
                if row == col:
                    BestParam = eval(f'self.Best{DataNames[col]}')[self.nBody]
                    ax.axvline(BestParam, color='red', linestyle='-', linewidth=0.75)
                elif row > col:
                    BestXParam = eval(f'self.Best{DataNames[col]}')[self.nBody]
                    BestYParam = eval(f'self.Best{DataNames[row]}')[self.nBody]
                    ax.plot(BestXParam, BestYParam, color='red', marker='x')

        self.WidgetPlot.Canvas.fig.subplots_adjust(left=0.10, bottom=0.13, right=0.94, top=0.95, wspace=0.05, hspace=0.07)
        self.capture_corner_layout_reference()
        self.apply_corner_square_layout()
        self.apply_corner_uniform_label_padding()
        self.connect_corner_axis_limit_sync()
        self._corner_axis_sync_ready = True

class PosAtDate(GeneralToolClass):
    def __init__(self, InputData, OutputParams, SelectOrbitsEllipses, BestOrbitParams, BestOrbitEllipse, LMOrbitParams, LMOrbitEllipse, SystDist):
        super().__init__('Position at date', 'Position of bodies at a given date', InputData, OutputParams, None, SelectOrbitsEllipses, BestOrbitParams, BestOrbitEllipse, LMOrbitParams, LMOrbitEllipse)

        # Parameters initialisation
        self.InitParams()
        self.SystDist = SystDist

        # Plot initialization
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylim=True)

    def InitParams(self):
        """Initialize parameters for the Position at Date tool."""
        # Equal aspect ratio
        # self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)

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

        # Show LM fit
        self.CheckLMFit = CheckBox('LM fit', 'Show the Levenberg-Marquardt fit')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckLMFit)
        self.CheckLMFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)
        self.CheckBestFit.CheckParam.stateChanged.connect(self.refresh_plots)

        # Show observations points
        self.CheckObs = CheckBox('Observations', 'Show the observations points with its error bar')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckObs)
        if self.InputData is None:
            self.CheckObs.CheckParam.setEnabled(False)
        self.CheckObs.CheckParam.stateChanged.connect(self.refresh_plots)

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
        # if self.CheckEqualAspect.CheckParam.isChecked(): self.Subplot.set_box_aspect(1)

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
        # print(self.Date)
        # xR, yR, zR = kepler_position(self.a[self.nBody], self.P[self.nBody], self.e[self.nBody], self.w[self.nBody], self.i[self.nBody], self.W[self.nBody], self.tp[self.nBody], self.Date, self.SystDist)

        # RaAtDate = -xR/self.SystDist*1000
        # DecAtDate = yR/self.SystDist*1000

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

        if self.CheckLMFit.CheckParam.isChecked():
            self.Subplot.plot(self.LMRa[self.nBody], self.LMDec[self.nBody], color='orange', lw=1)

            LMPeriod = np.max(self.LMt[self.nBody]) - np.min(self.LMt[self.nBody])
            LMDate = self.Date

            while LMDate < np.min(self.LMt[self.nBody]):
                LMDate += LMPeriod

            while LMDate > np.max(self.LMt[self.nBody]):
                LMDate -= LMPeriod

            indexLMDate = np.argmin(np.abs(self.LMt[self.nBody] - LMDate))
            self.Subplot.plot(self.LMRa[self.nBody][indexLMDate], self.LMDec[self.nBody][indexLMDate], marker='x', color='orange', markersize=8, markeredgewidth=2)

        if self.CheckBestFit.CheckParam.isChecked():
            self.Subplot.plot(self.BestRa[self.nBody], self.BestDec[self.nBody], color='r', lw=1)

            BestPeriod = np.max(self.Bestt[self.nBody]) - np.min(self.Bestt[self.nBody])
            BestDate = self.Date

            while BestDate < np.min(self.Bestt[self.nBody]):
                BestDate += BestPeriod

            while BestDate > np.max(self.Bestt[self.nBody]):
                BestDate -= BestPeriod

            indexBestDate = np.argmin(np.abs(self.Bestt[self.nBody] - BestDate))
            self.Subplot.plot(self.BestRa[self.nBody][indexBestDate], self.BestDec[self.nBody][indexBestDate], marker='x', color='red', markersize=8, markeredgewidth=2)
            

        if self.CheckObs.CheckParam.isChecked():
            ra = self.InputData['Planets']['DataAstrom']['Ra'][self.nBody]
            dec = self.InputData['Planets']['DataAstrom']['Dec'][self.nBody]
            dra = self.InputData['Planets']['DataAstrom']['dRA'][self.nBody]
            ddec = self.InputData['Planets']['DataAstrom']['dDec'][self.nBody]
            dates = self.InputData['Planets']['DataAstrom']['Date'][self.nBody]
            self.Subplot.errorbar(ra, dec, ddec, dra, linestyle='', color='white', linewidth=1)

        # Plot features
        self.Subplot.set_xlabel(r'$\delta$RA [mas]')
        self.Subplot.set_ylabel(r'$\delta$Dec [mas]')
        self.Subplot.invert_xaxis()
        self.Subplot.set_aspect('equal', adjustable='box')
        self.Subplot.set_xlim(xlim_init)
        self.Subplot.set_ylim(ylim_init)

