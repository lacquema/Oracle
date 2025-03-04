#! /Users/lacquema/ByeGildas/bin/python3

### --- Packages --- ###

# Transverse packages
import sys
import numpy as np
from random import random, randint
import corner
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import gaussian_filter1d

# PyQt packages
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QDateEdit
from PyQt6.QtCore import QDateTime, QDate, QSize
from Utils import date_to_jd, jd_to_mjd, mjd_to_jd, jd_to_date

# My packages
from Parameters import *
from BestOrbits import BestOrbitsClass
from SelectOrbits import SelectOrbitsClass

from WindowPlot import WindowPlot



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
        self.WindowPlot.SignalCloseWindowPlot.connect(lambda: self.BtnPlot.setEnabled(True)) # reception of the closeEvent of the parameter window and set enabled the associed button

        # Connections between parameters and plot
        self.WindowPlot.WidgetParam.BtnReset.clicked.connect(self.ResetParams)
        self.WindowPlot.WidgetParam.BtnRefresh.clicked.connect(self.Refresh_ActivePlots)

        # Initialisation of Data
        if InputData != None: self.NbInputData_RA, self.I_RA, self.MJD_RA, self.JJ_RA, self.MM_RA, self.YY_RA, self.Ra, self.Dec, self.DRa,  self.DDec, self.Corr_DecRa, self.Sep, self.Pa, self.DSep, self.DPa, self.Corr_SepPa, self.Source_RA = InputData
        
        if OutputParams != None: self.NbBodies, self.NbOrbits, self.P, self.a, self.e, self.i, self.w, self.W, self.tp, self.m, self.Mdyn, self.Chi2, self.map = OutputParams
       
        if SelectOrbitsParams != None: self.NbBodies, self.NbSelectOrbits, self.SelectP, self.Selecta, self.Selecte, self.Selecti, self.Selectw, self.SelectW, self.Selecttp, self.Selectm, self.SelectMdyn, self.SelectChi2 = SelectOrbitsParams
        if SelectOrbitsEllipses != None: self.NbBodies, self.NbSelectOrbits, self.NbPtsEllipse, self.SelectP, self.Selectt, self.SelectX, self.SelectY, self.SelectZ = SelectOrbitsEllipses

        if BestOrbitsParams != None: self.NbBodies, self.BestP, self.Besta, self.Beste, self.Besti, self.Bestw, self.BestW, self.Besttp, self.Bestm, self.BestMdyn, self.BestChi2 = BestOrbitsParams
        if BestOrbitsEllipses != None: self.NbBodies, self.NbPtsEllipse, self.BestP, self.Bestt, self.BestX, self.BestY, self.BestZ = BestOrbitsEllipses

        self.InitParams()
        # self.WindowPlot.WidgetParam.setFixedSize(self.WindowPlot.WidgetParam.baseSize())

        # Fixer la taille pour éviter le redimensionnement des parametres.
        left_width = self.WindowPlot.WidgetParam.sizeHint().width()
        self.WindowPlot.WidgetParam.setFixedWidth(left_width)

        # Widget container
        self.setLayout(Layout) # GeneralToolClass is directly the widget container

    # # Open the parameters window when the parameters button is clicked
    # def Toggle_WindowParam(self):
    #         # self.WindowPlot.WidgetParam.move(self.x()-500, self.pos().y())
    #         self.WindowPlot.WidgetParam.show()
    #         self.BtnParam.setEnabled(False)

    # # Open the plot window when the Plot button is clicked
    # def Toggle_WindowPlot(self):
    #         self.Plot()
    #         self.WindowPlot.WidgetPlot.show()
    #         self.BtnPlot.setEnabled(False)

    def InitParams(self):
        return

    # Refresh all active plot when the refresh button is clicked
    def Refresh_ActivePlots(self):
        if self.WindowPlot.isVisible():
            self.Plot()

    # Close programme when the main window are closed
    def closeEvent(self, e):
        app.closeAllWindows()

    # Reset all widgets of the parameters window
    def ResetParams(self):
        for i in reversed(range(2, self.WindowPlot.WidgetParam.Layout.count())): 
            WidgetToRemove = self.WindowPlot.WidgetParam.Layout.itemAt(i).widget()
            self.WindowPlot.WidgetParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        self.InitParams()

    # Giver of label
    def LabelOf(self, var=str):
        if var == 'P':
            return 'Period [yr]'
        elif var == 'a':
            return 'Semi-major axis [AU]'
        elif var == 'e':
            return 'Eccentricity'
        elif var == 'i':
            return 'Inclinaison [°]'
        elif var == 'w':
            return 'Argument of periastron [°]'
        elif var == 'W':
            return 'Longitude of ascending node [°]'
        elif var == 'tp':
            return 'Periastron time passage [MJD]'
        elif var == 'm':
            return 'Body mass [Mjup]'
        elif var == 'Mdyn':
            return 'Dynamical mass [Mjup]'
        elif var == 'Chi2':
            return 'Chi square'
        else:
            return 'Unkown variable'

    def Toggle_WindowPlot(self):
        self.Plot()
        self.WindowPlot.show()
        self.BtnPlot.setEnabled(False)


class SpaceView(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Space view', 'Space view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Parameters initialisation
        # self.InitParams()
        


    # Parameters initialisation
    def InitParams(self):

        # Number of studied body
        self.ListBody = []
        for k in range(self.NbBodies):
            self.ListBody.append(str(k+1))
        self.ListBody.append('all')
        self.nBodyWidget = ComboBox("Orbit number",'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
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

    # Parameters update
    def UpdateParams(self):
        if self.nBodyWidget.ComboParam.currentText() == 'all':
            self.nBody = 'all'
        else:
            self.nBody = int(self.nBodyWidget.ComboParam.currentText())-1
        self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.NbShownOrbits = self.NbShownOrbitsWidget.SpinParam.value()

    # Plot
    def Plot(self):

        # Clear axis
        for i in range(len(self.WindowPlot.WidgetPlot.Canvas.fig.axes)):
            self.WindowPlot.WidgetPlot.Canvas.fig.delaxes(self.WindowPlot.WidgetPlot.Canvas.fig.axes[0])

        # Update of parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        if self.indexView == 0: # 2D (x,y)

            # Initialisation of axis
            self.Subplot2D = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111, aspect='equal')

            # Plot
            self.Subplot2D.plot(0, 0, marker='*', color='orange', markersize=10)

            if self.nBody == 'all': 
                for k in range(self.NbBodies):
                    if self.CheckBestFit.CheckParam.isChecked(): 
                        self.Subplot2D.plot(self.BestX[k], self.BestY[k], color='r')
                    for n in range(self.NbShownOrbits):
                        self.Subplot2D.plot(self.SelectX[k][n], self.SelectY[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
            else: 
                if self.CheckBestFit.CheckParam.isChecked(): 
                    self.Subplot2D.plot(self.BestX[self.nBody], self.BestY[self.nBody], color='r')
                for n in range(self.NbShownOrbits):
                        self.Subplot2D.plot(self.SelectX[self.nBody][n], self.SelectY[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)

            if self.CheckObs.CheckParam.isChecked():
                self.Subplot2D.errorbar(self.Ra, self.Dec, self.DRa, self.DDec, linestyle='', color='b') # Observed data

            # Plot features
            self.Subplot2D.set_xlabel(r'$\delta Ra$ [mas]')
            self.Subplot2D.set_ylabel(r'$\delta Dec$ [mas]')
            self.Subplot2D.invert_xaxis()


        elif self.indexView == 1: # 2D (x,z)

            # Initialisation of axis
            self.Subplot2D = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111, aspect='equal')

            # Plot
            self.Subplot2D.plot(0, 0, marker='*', color='orange', markersize=10)

            if self.nBody == 'all': 
                for k in range(self.NbBodies):
                    if self.CheckBestFit.CheckParam.isChecked(): 
                        self.Subplot2D.plot(self.BestX[k], self.BestZ[k], color='r')
                    for n in range(self.NbShownOrbits):
                        self.Subplot2D.plot(self.SelectX[k][n], self.SelectZ[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
            else: 
                if self.CheckBestFit.CheckParam.isChecked(): 
                    self.Subplot2D.plot(self.BestX[self.nBody], self.BestZ[self.nBody], color='r')
                for n in range(self.NbShownOrbits):
                        self.Subplot2D.plot(self.SelectX[self.nBody][n], self.SelectZ[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)

            # Plot features
            self.Subplot2D.set_xlabel(r'$\delta Ra$ [mas]')
            self.Subplot2D.set_ylabel('Depth [mas]')
            self.Subplot2D.invert_xaxis()
        

        elif self.indexView == 2: # 3D

            # Initialisation of axis
            self.Subplot3D = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111, projection='3d', aspect='equal')

            # Plot
            self.Subplot3D.plot(0, 0, 0, marker='*', color='orange', markersize=10)

            if self.nBody == 'all': 
                for k in range(self.NbBodies):
                    if self.CheckBestFit.CheckParam.isChecked(): 
                        self.Subplot3D.plot(self.BestX[k], self.BestZ[k], self.BestZ[k], color='r')
                    for n in range(self.NbShownOrbits):
                        self.Subplot3D.plot(self.SelectX[k][n], self.SelectY[k][n], self.SelectZ[k][n], color=self.colorList[k], linestyle='-', linewidth=0.3, alpha=0.1)
            else: 
                if self.CheckBestFit.CheckParam.isChecked(): 
                    self.Subplot3D.plot(self.BestX[self.nBody], self.BestY[self.nBody], self.BestZ[self.nBody], color='r')
                for n in range(self.NbShownOrbits):
                        self.Subplot3D.plot(self.SelectX[self.nBody][n], self.SelectY[self.nBody][n], self.SelectZ[self.nBody][n], color=self.colorList[self.nBody], linestyle='-', linewidth=0.3, alpha=0.1)

            # Plot features
            self.Subplot3D.set_xlabel(r'$\delta Ra$ [mas]')
            self.Subplot3D.set_ylabel(r'$\delta Dec$ [mas]')
            self.Subplot3D.set_zlabel('Depth [mas]')
            self.Subplot3D.invert_xaxis()



        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.draw()


class TempoView(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Temporal view', 'Temporal view of fit orbits', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Plots initialisation
        self.Subplot1 = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(211)
        self.Subplot2 = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(212)

        # Parameters initialisation
        # self.InitParams()

    # Parameters initialisation
    def InitParams(self):

        # Orbit number
        self.ListBody = []
        for k in range(self.NbBodies):
            self.ListBody.append(str(k+1))
        self.nBodyWidget = ComboBox("Orbit number",'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Number of shown orbits
        self.NbShownOrbits = 500
        self.NbShownOrbitsWidget = SpinBox('Number of orbits', 'Number of shown orbits', ParamDefault=self.NbShownOrbits, ParamMin=0, ParamMax=self.NbSelectOrbits)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbShownOrbitsWidget)

        # Choice of coordinate
        self.CoordinateWidget = ComboBox('Choice of coordinate', 'Coordinates', ['dRa', 'dDec'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CoordinateWidget)

    # Parameters update
    def UpdateParams(self):
        self.CoordinateIndex = self.CoordinateWidget.ComboParam.currentIndex()
        if self.CoordinateIndex == 0:
            self.Coordinate = r'$\delta Ra$'
        elif self.CoordinateIndex == 1:
            self.Coordinate = r'$\delta Dec$'
        self.nBody = int(self.nBodyWidget.ComboParam.currentText())-1
        self.NbShownOrbits = self.NbShownOrbitsWidget.SpinParam.value()


    # Plot
    def Plot(self):

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

        # Plots with current paramaters
        if self.CoordinateIndex == 0:
            YplotOutput = self.SelectX
            BestYplotOutput = self.BestX
            YplotInput = self.Ra
            YplotInputErr = self.DRa

        elif self.CoordinateIndex == 1:
            YplotOutput = self.SelectY
            BestYplotOutput = self.BestY
            YplotInput = self.Dec
            YplotInputErr = self.DDec

        # Output data
        for n in range(self.NbShownOrbits):
            Selectt3P = np.concatenate((self.Selectt[self.nBody][n]-self.SelectP[self.nBody][n]*365.25, self.Selectt[self.nBody][n], self.Selectt[self.nBody][n]+self.SelectP[self.nBody][n]*365.25))
            YplotOutput3P = np.concatenate((YplotOutput[self.nBody][n], YplotOutput[self.nBody][n], YplotOutput[self.nBody][n]))
            self.Subplot1.plot(Selectt3P, YplotOutput3P, color=self.colorList[0], linestyle='-', linewidth=0.2, alpha=0.1)

        Bestt3P = np.concatenate((self.Bestt[self.nBody]-self.BestP[self.nBody]*365.25, self.Bestt[self.nBody], self.Bestt[self.nBody]+self.BestP[self.nBody]*365.25))
        BestYplotOutput3P = np.concatenate((BestYplotOutput[self.nBody], BestYplotOutput[self.nBody], BestYplotOutput[self.nBody]))
        self.Subplot1.plot(Bestt3P, BestYplotOutput3P, linestyle='-', linewidth=0.5, color='r')

        # Input data
        for k in range(self.NbInputData_RA):
            if self.I_RA[k] == self.nBody+1:
                self.Subplot1.errorbar(self.MJD_RA[k], YplotInput[k], YplotInputErr[k], linestyle='', color='b')
                indext = np.argmin(np.abs(Bestt3P-self.MJD_RA[k])) # index of time of output data closer than time of input data
                Res = BestYplotOutput3P[indext] - YplotInput[k] # Residu
                self.Subplot2.errorbar(Bestt3P[indext], Res, YplotInputErr[k], color='b')

        self.Subplot2.hlines(0, np.min(Bestt3P), np.max(Bestt3P), color='red', linewidth = 0.5)

        # Plot features
        self.Subplot1.set_ylabel(self.Coordinate+' [mas]')
        self.Subplot2.set_ylabel(self.Coordinate+' - Bestfit [mas]')
        self.Subplot2.set_xlabel('Time [MJD]')
        L = np.max(self.MJD_RA)-np.min(self.MJD_RA)
        self.Subplot1.set_xlim(np.min(self.MJD_RA)-0.1*L, np.max(self.MJD_RA)+0.1*L)
        self.Subplot2.set_xlim(self.Subplot1.get_xlim())
        # self.Subplot1.set_xticklabels([])
        self.Subplot1.grid()
        self.Subplot2.grid()


        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.draw()


class Conv(GeneralToolClass):
    def __init__(self, OutputParams):
        super().__init__('Convergence', 'Convergence of the fit orbit parameters', None, OutputParams, None, None, None, None)

        # Plot initialisation
        self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

        # Parameters initialisation
        # self.InitParams()

    # Parameters initialisation
    def InitParams(self):

        # Orbit number
        self.ListBody = []
        for k in range(self.NbBodies):
            self.ListBody.append(str(k+1))
        self.nBodyWidget = ComboBox("Orbit number",'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Orbital parameter', 'Orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp','m', 'Mdyn', 'Chi2'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)

    # Parameters update
    def UpdateParams(self):
        self.nBody = int(self.nBodyWidget.ComboParam.currentText())-1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()

    # Plot
    def Plot(self):

        # Clear the plot
        self.Subplot.cla()

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return
        

        # Plot with current paramaters
        self.Steps = range(self.NbOrbits)
        self.EvalParamOrbit = eval('self.'+self.ParamOrbit)[self.nBody]

        self.Subplot.plot(self.Steps, self.EvalParamOrbit, marker=',', linestyle='')

        # Plot features
        self.Subplot.set_xlabel('Step')
        self.Subplot.set_ylabel(self.LabelOf(self.ParamOrbit))

        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.draw()


class Hist(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitsParams):
        super().__init__('Histogram', "Histogram of an orbital parameter", None, OutputParams, None, None, BestOrbitsParams, None)

        # Plot initialisation
        self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

        # Parameters initialisation
        # self.InitParams()

    # Parameters initialisation
    def InitParams(self):

        # Orbit number
        self.ListBody = []
        for k in range(self.NbBodies):
            self.ListBody.append(str(k+1))
        self.nBodyWidget = ComboBox("Orbit number",'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Orbital parameter', 'Orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp','m', 'Mdyn', 'Chi2'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)

        # Histogram binning
        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        # self.CheckBestFit.CheckParam.setChecked(True)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)

        # Show confidence interval
        self.CheckMedian = CheckBox('Median', 'Show the median and the 1 sigma confidence interval')
        # self.CheckBestFit.CheckParam.setChecked(True)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckMedian)

        self.CheckMedian.Layout.addSpacing(20)

        self.IntConf = 68
        self.IntConfWidget = SpinBox('Confidence', 'Acceptable level of confidence', self.IntConf, 0, 100, 1)
        self.CheckMedian.Layout.addWidget(self.IntConfWidget)
        self.IntConfWidget.setEnabled(self.CheckMedian.CheckParam.isChecked())
        self.CheckMedian.CheckParam.stateChanged.connect(lambda state: self.IntConfWidget.setEnabled(state))

    # Parameters update
    def UpdateParams(self):
        self.nBody = int(self.nBodyWidget.ComboParam.currentText())-1
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.IntConf = self.IntConfWidget.SpinParam.value()/100
        
    # Plot
    def Plot(self):
        
        # Clear the plot
        self.Subplot.cla()

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return
        
        # Plot with current parameters
        self.EvalParamOrbit = eval('self.'+self.ParamOrbit)[self.nBody]
        self.Subplot.hist(self.EvalParamOrbit, self.NbBins)

        # Best fit
        if self.CheckBestFit.CheckParam.isChecked():
            BestParam = eval('self.Best'+self.ParamOrbit)[self.nBody]

            self.Subplot.axvline(BestParam, color='red')
            self.Subplot.text(BestParam, 0.5*self.Subplot.get_ylim()[1], s='{}'.format(np.around(BestParam,3)), color='r', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='r'), fontsize = 9, horizontalalignment='center', verticalalignment='center', rotation=0)

        # Median
        # if self.CheckMedian.CheckParam.isChecked():

        #     counts, bin_edges = np.histogram(self.EvalParamOrbit, self.NbBins)
        #     bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Centres des bins

        #     extended_counts = np.concatenate([counts, counts])  # Copie de l’histogramme
        #     extended_bins = np.concatenate([bin_edges[:-1], bin_edges[:-1] + 360])  # Double les bins

        #     # Détection des pics avec une contrainte sur la proéminence
        #     prom = np.median(counts[counts > 0]) * 0.3
        #     peaks, properties = find_peaks(counts, prominence=prom)

        #     # Ramener les pics dans [0, 360]
        #     peaks = np.sort(peaks % len(counts))

        #     # Fusionner les pics en double (distance proche de 360°)
        #     merged_peaks = []
        #     for i, p in enumerate(peaks):
        #         print(i, p)
        #         if i == 0 or p - peaks[i - 1] > len(counts) * 0.2:  # On fusionne si trop proche
        #             merged_peaks.append(p)

        #     merged_peaks = np.array(merged_peaks)

        #     # Calcul des largeurs des pics à mi-hauteur
        #     rh = 0.5 if prom > np.median(counts) * 0.3 else 0.3
        #     width_results = peak_widths(extended_counts, merged_peaks, rel_height=rh)

        #     # Extraction des valeurs utiles
        #     left_ips = width_results[2] % len(counts)  # Positions des bords gauches en indices flottants
        #     right_ips = width_results[3] % len(counts)  # Positions des bords droits en indices flottants
        #     left_bases = properties["left_bases"]  # Vallée à gauche
        #     right_bases = properties["right_bases"]  # Vallée à droite

        #     # Conversion des indices flottants en entiers pour l'accès aux bins
        #     left_indices = left_ips.astype(int)
        #     right_indices = right_ips.astype(int)

        #     # Calcul des bornes finales en prenant la vallée la plus lointaine
        #     left_bounds = np.minimum(left_indices, left_bases)
        #     right_bounds = np.maximum(right_indices, right_bases)

        #     for i in range(len(peaks)):
        #         left = max(0, left_bounds[i])
        #         right = min(right_bounds[i], len(counts) - 1)
    
        #         # Calcul de la distribution cumulée des comptages
        #         cumsum = np.cumsum(counts[left:right+1])  # Somme cumulée des valeurs de l'histogramme
        #         total = cumsum[-1]  # Nombre total de comptages

        #         # Niveaux de fraction de comptage à atteindre
        #         Ic = self.IntConf/100  # Intervalle de confiance (90%)
        #         seuil_median = 0.50 * total  # Médiane à 50% du total des comptages
        #         seuil_inf = (0.50 - Ic / 2) * total  # Borne inférieure
        #         seuil_sup = (0.50 + Ic / 2) * total  # Borne supérieure

        #         # Trouver les valeurs correspondantes aux seuils
        #         mediane = bin_centers[np.searchsorted(cumsum, seuil_median)+left]
        #         lower_bound = bin_centers[np.searchsorted(cumsum, seuil_inf)+left]
        #         upper_bound = bin_centers[np.searchsorted(cumsum, seuil_sup)+left]

        #         SigmaMinus = np.abs(mediane-lower_bound)
        #         SigmaPlus = np.abs(mediane-upper_bound)

        #         self.Subplot.axvline(mediane, linestyle='-', color='black')
        #         self.Subplot.text(mediane, 0.83*self.Subplot.get_ylim()[1], s='{}'.format(np.around(mediane,3)), bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize = 9, horizontalalignment='center', verticalalignment='bottom')

        #         self.Subplot.axvline(lower_bound, linestyle='--', color='black')
        #         self.Subplot.text(lower_bound, 0.8*self.Subplot.get_ylim()[1], s='-{}'.format(np.around(SigmaMinus,3)), bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize = 8, horizontalalignment='right', verticalalignment='top')

        #         self.Subplot.axvline(upper_bound, linestyle='--', color='black')
        #         self.Subplot.text(upper_bound, 0.8*self.Subplot.get_ylim()[1], s='+{}'.format(np.around(SigmaPlus,3)), bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize = 8, horizontalalignment='left', verticalalignment='top')





            # data_case2 = np.concatenate([
            #     np.random.normal(loc=180, scale=20, size=500),  # Pic à 180°
            #     np.random.normal(loc=0, scale=20, size=500) % 360  # Pic réparti entre 0° et 360°
            #                             ])

            # # Création de l'histogramme
            # counts, bin_edges = np.histogram(data_case2, self.NbBins)
            # bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Centres des bins

            # # Extension de l'histogramme pour une détection cyclique des pics
            # extended_counts = np.concatenate([counts, counts])
            # extended_bins = np.concatenate([bin_edges[:-1], bin_edges[:-1] + 360])

            # # Détection des pics avec une contrainte sur la proéminence
            # prominence_threshold = np.median(counts[counts > 0]) * 0.3
            # peaks, properties = find_peaks(counts, prominence=prominence_threshold)

            # # Ramener les pics dans [0, 360] et fusionner les pics proches
            # peaks = np.sort(peaks % len(counts))
            # merged_peaks = [p for i, p in enumerate(peaks) if i == 0 or p - peaks[i - 1] > len(counts) * 0.2]

            # # Calcul des largeurs des pics à mi-hauteur
            # rel_height = 0.5 if prominence_threshold > np.median(counts) * 0.3 else 0.3
            # width_results = peak_widths(extended_counts, merged_peaks, rel_height=rel_height)

            # # Extraction des bornes des pics
            # left_bounds = width_results[2].astype(int) % len(counts)
            # right_bounds = width_results[3].astype(int) % len(counts)

            # # Intervalle de confiance
            # Ic = self.IntConf / 100  

            # for peak, left, right in zip(merged_peaks, left_bounds, right_bounds):
            #     left = max(0, left)
            #     right = min(right, len(counts) - 1)

            #     # Calcul de la distribution cumulée
            #     cumsum = np.cumsum(counts[left:right + 1])
            #     total = cumsum[-1]

            #     # Calcul des seuils
            #     seuils = {
            #         "median": 0.50 * total,
            #         "inf": (0.50 - Ic / 2) * total,
            #         "sup": (0.50 + Ic / 2) * total
            #     }

            #     # Détermination des valeurs associées
            #     mediane = bin_centers[left + np.searchsorted(cumsum, seuils["median"])]
            #     lower_bound = bin_centers[left + np.searchsorted(cumsum, seuils["inf"])]
            #     upper_bound = bin_centers[left + np.searchsorted(cumsum, seuils["sup"])]

            #     SigmaMinus = abs(mediane - lower_bound)
            #     SigmaPlus = abs(mediane - upper_bound)

            #     # Tracés des lignes et annotations sur le graphique
            #     self.Subplot.axvline(mediane, linestyle='-', color='black')
            #     self.Subplot.text(mediane, 0.83 * self.Subplot.get_ylim()[1], f'{mediane:.3f}', 
            #                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'),
            #                     fontsize=9, ha='center', va='bottom')

            #     for bound, sigma, ha in [(lower_bound, SigmaMinus, 'right'), (upper_bound, SigmaPlus, 'left')]:
            #         self.Subplot.axvline(bound, linestyle='--', color='black')
            #         self.Subplot.text(bound, 0.8 * self.Subplot.get_ylim()[1], f'{"+" if ha == "left" else "-"}{sigma:.3f}', 
            #                         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'),
            #                         fontsize=8, ha=ha, va='top')

        if self.CheckMedian.CheckParam.isChecked():

            # Calcul de l'histogramme
            counts, bin_edges = np.histogram(self.EvalParamOrbit, self.NbBins)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Centres des bins

            # Dupliquer et concaténer les données pour gérer la cyclicité
            extended_counts = np.concatenate([counts, counts])
            extended_bin_centers = np.concatenate([bin_centers, bin_centers + 360])

            # Appliquer un filtre gaussien pour lisser les données
            smoothed_counts = gaussian_filter1d(extended_counts, sigma=2)

            # Détection des pics
            prom = np.median(smoothed_counts[smoothed_counts > 0]) * 0.3
            peaks, properties = find_peaks(smoothed_counts, prominence=prom)

            if len(peaks) == 3:
                unique_peaks = [peaks[1], peaks[2]]
                unique_left_bases = [properties["left_bases"][1], properties["left_bases"][2]]
                unique_right_bases = [properties["right_bases"][1], properties["right_bases"][2]]
            elif len(peaks) == 2:
                unique_peaks = peaks
                unique_left_bases = properties["left_bases"]
                unique_right_bases = properties["right_bases"]

        
            left_ips = unique_left_bases
            right_ips = unique_right_bases

            # Tracé des résultats
            # self.Subplot.plot(bin_centers, smoothed_counts[:len(bin_centers)], label='Données lissées')

            for i in range(len(unique_peaks)):
                left = max(0, left_ips[i])
                right = min(right_ips[i], len(extended_counts) - 1)

                mediane, lower_bound, upper_bound = self.calculate_bounds(smoothed_counts, extended_bin_centers, left, right, Ic=self.IntConf)
                
                self.plot_bounds(self.Subplot, mediane, lower_bound, upper_bound)
                

        # Plot features
        self.Subplot.set_xlabel(self.LabelOf(self.ParamOrbit))
        self.Subplot.set_ylabel('Count number')
        
        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.draw()



    def calculate_bounds(self, counts, bin_centers, left, right, Ic=0.68):

        cumsum = np.cumsum(counts[left:right+1])

        total = cumsum[-1]

        seuil_median = 0.50 * total
        seuil_inf = (0.50 - Ic / 2) * total
        seuil_sup = (0.50 + Ic / 2) * total

        mediane = bin_centers[np.searchsorted(cumsum, seuil_median) + left]
        lower_bound = bin_centers[np.searchsorted(cumsum, seuil_inf) + left]
        upper_bound = bin_centers[np.searchsorted(cumsum, seuil_sup) + left]

        return mediane, lower_bound, upper_bound
    

    def plot_bounds(self, Subplot, mediane, lower_bound, upper_bound):

        SigmaMinus = np.abs(mediane - lower_bound)
        SigmaPlus = np.abs(mediane - upper_bound)

        mediane = mediane % 360
        lower_bound = lower_bound % 360
        upper_bound = upper_bound % 360

        Subplot.axvline(mediane, linestyle='-', color='black')
        Subplot.text(mediane, 0.83 * Subplot.get_ylim()[1], s='{}'.format(np.around(mediane, 3)),
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize=9,
                    horizontalalignment='center', verticalalignment='bottom')

        Subplot.axvline(lower_bound, linestyle='--', color='black')
        Subplot.text(lower_bound, 0.8 * Subplot.get_ylim()[1], s='-{}'.format(np.around(SigmaMinus, 3)),
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize=8,
                    horizontalalignment='right', verticalalignment='top')

        Subplot.axvline(upper_bound, linestyle='--', color='black')
        Subplot.text(upper_bound, 0.8 * Subplot.get_ylim()[1], s='+{}'.format(np.around(SigmaPlus, 3)),
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='black'), fontsize=8,
                    horizontalalignment='left', verticalalignment='top')


class Hist2D(GeneralToolClass):
    def __init__(self, OutputParams, BestOrbitsParams):
        super().__init__('Histogram 2D', 'Histogram of an orbital parameter as fonction of another', None, OutputParams, None, None, BestOrbitsParams, None)

        # Parameters initialisation
        # self.InitParams()

    # Parameters initialisation
    def InitParams(self):

        # Orbit number
        self.ListBody = []
        for k in range(self.NbBodies):
            self.ListBody.append(str(k+1))
        self.nBodyWidget = ComboBox("Orbit number",'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Abscissa orbit parameters
        self.XParamOrbitWidget = ComboBox('X Orbital parameter', 'Abscissa orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp','m', 'Mdyn'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XParamOrbitWidget)

        # Ordinate orbit parameters
        self.YParamOrbitWidget = ComboBox('Y Orbital parameter', 'Ordinate orbit Parameter', ['P', 'a', 'e', 'i', 'w', 'W', 'tp','m', 'Mdyn'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YParamOrbitWidget)

        # Histogram binning
        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Show best fit
        self.CheckBestFit = CheckBox('Best fit', 'Show the fit with the best Chi2')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBestFit)


    # Parameters update
    def UpdateParams(self):
        self.nBody = int(self.nBodyWidget.ComboParam.currentText())-1
        self.XParamOrbit = self.XParamOrbitWidget.ComboParam.currentText()
        self.YParamOrbit = self.YParamOrbitWidget.ComboParam.currentText()
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        
    # Plot
    def Plot(self):
        
        # Clear axis
        for i in range(len(self.WindowPlot.WidgetPlot.Canvas.fig.axes)):
            self.WindowPlot.WidgetPlot.Canvas.fig.delaxes(self.WindowPlot.WidgetPlot.Canvas.fig.axes[0])

        # Plot initialisation
        self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return
        
        # Plot with current parameters
        self.EvalXParamOrbit = eval('self.'+self.XParamOrbit)[self.nBody]
        self.EvalYParamOrbit = eval('self.'+self.YParamOrbit)[self.nBody]

        hist = self.Subplot.hist2d(self.EvalXParamOrbit, self.EvalYParamOrbit, (self.NbBins, self.NbBins))
        ColorbarAx = make_axes_locatable(self.Subplot).append_axes('right', size='5%', pad=0.1)
        self.WindowPlot.WidgetPlot.Canvas.fig.colorbar(hist[3], ColorbarAx, label='Count number')

        # Best fit
        if self.CheckBestFit.CheckParam.isChecked():
            BestXParam = eval('self.Best'+self.XParamOrbit)[self.nBody]
            BestYParam = eval('self.Best'+self.YParamOrbit)[self.nBody]

            self.Subplot.plot(BestXParam, BestYParam, color='red', marker='x')

        # Plot features
        self.Subplot.set_xlabel(self.LabelOf(self.XParamOrbit))
        self.Subplot.set_ylabel(self.LabelOf(self.YParamOrbit))
        
        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.draw()
    

class PosAtDate(GeneralToolClass):
    def __init__(self, InputData, SelectOrbitsEllipses, BestOrbitsEllipses):
        super().__init__('Position at date', 'Position of bodies at a given date', InputData, None, None, SelectOrbitsEllipses, None, BestOrbitsEllipses)

        # Parameters initialisation
        # self.InitParams()

    # Parameters initialisation
    def InitParams(self):

        # Orbit number
        self.ListBody = []
        for k in range(self.NbBodies):
            self.ListBody.append(str(k+1))
        self.nBodyWidget = ComboBox("Orbit number",'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
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

    # Parameters update
    def UpdateParams(self):
        self.nBody = int(self.nBodyWidget.ComboParam.currentText())-1
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.Date = self.DateWidget.MJDWidget.SpinParam.value()
        
    # Plot
    def Plot(self):
        
        # Clear axis
        for i in range(len(self.WindowPlot.WidgetPlot.Canvas.fig.axes)):
            self.WindowPlot.WidgetPlot.Canvas.fig.delaxes(self.WindowPlot.WidgetPlot.Canvas.fig.axes[0])

        # Plot initialisation
        self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlot.Canvas.draw()
            return

        # Compute date
        SelectXAtDate, SelectYAtDate = [np.zeros(self.NbSelectOrbits) for k in range(2)]

        for k in range(self.NbSelectOrbits):
            
            SelectPeriod = np.max(self.Selectt[self.nBody][k]) - np.min(self.Selectt[self.nBody][k])
            SelectDate = self.Date

            while SelectDate < np.min(self.Selectt[self.nBody][k]):
                SelectDate += SelectPeriod

            while SelectDate > np.max(self.Selectt[self.nBody][k]):
                SelectDate -= SelectPeriod
            
            indexBestDate = np.argmin(np.abs(self.Selectt[self.nBody][k] - SelectDate))
            SelectXAtDate[k] = self.SelectX[self.nBody][k][indexBestDate]
            SelectYAtDate[k] = self.SelectY[self.nBody][k][indexBestDate]

        # Plot with current parameters
        hist = self.Subplot.hist2d(SelectXAtDate, SelectYAtDate, bins=(self.NbBins, self.NbBins), range=((np.min(self.SelectX), np.max(self.SelectX)),(np.min(self.SelectY), np.max(self.SelectY))))
        ColorbarAx = make_axes_locatable(self.Subplot).append_axes('right', size='5%', pad=0.1)
        self.WindowPlot.WidgetPlot.Canvas.fig.colorbar(hist[3], ColorbarAx, ticks=[], label='Probability')

        self.Subplot.plot(0, 0, marker='*', color='orange', markersize=10)

        if self.CheckBestFit.CheckParam.isChecked():
            self.Subplot.plot(self.BestX[self.nBody], self.BestY[self.nBody], color='r', lw=0.5)

            BestPeriod = np.max(self.Bestt[self.nBody]) - np.min(self.Bestt[self.nBody])
            BestDate = self.Date

            while BestDate < np.min(self.Bestt[self.nBody]):
                BestDate += BestPeriod

            while BestDate > np.max(self.Bestt[self.nBody]):
                BestDate -= BestPeriod

            indexBestDate = np.argmin(np.abs(self.Bestt[self.nBody] - BestDate))
            self.Subplot.plot(self.BestX[self.nBody][indexBestDate], self.BestY[self.nBody][indexBestDate], marker='x', color='red')

        if self.CheckObs.CheckParam.isChecked():
            self.Subplot.errorbar(self.Ra, self.Dec, self.DRa, self.DDec, linestyle='') # Observed data

        # Plot features
        self.Subplot.set_xlabel(r'$\delta$ Ra')
        self.Subplot.set_ylabel(r'$\delta$ Dec')
        self.Subplot.invert_xaxis()

        
        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.draw()
    

class Corner(GeneralToolClass):
    def __init__(self, SelectOrbitsParams):
        super().__init__('Corner', 'Corner plot of parameters', None, None, SelectOrbitsParams, None, None, None)

        # Parameters initialisation
        # self.InitParams()

    # Parameters initialisation
    def InitParams(self):

        # Orbit number
        self.ListBody = []
        for k in range(self.NbBodies):
            self.ListBody.append(str(k+1))
        self.nBodyWidget = ComboBox("Orbit number",'Number of the orbit which is studied (counting from the center of the system outwards)', self.ListBody)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nBodyWidget)
        if self.NbBodies == 1:
            self.nBodyWidget.setEnabled(False)

        # Parameters include on the corner
        self.CheckLayout = QHBoxLayout()

        self.CheckLabel = QLabel('Orbit parameters :')
        self.CheckLayout.addWidget(self.CheckLabel)

        self.OrbitParams = [['P', 'a', 'e', 'i', 'w', 'W', 'tp', 'm', 'Mdyn'],
                            ['Period', 'Semi-major axis', 'Eccentricity', 'Inclinaison', 'Argument of periastron', 'Longitude of ascending node', 'Periastron time passage', 'Body mass', 'Dynamical mass']]
        self.WidgetOrbitParams = []

        for i in range(len(self.OrbitParams[0])):
            CheckParamWidget = CheckBox(self.OrbitParams[0][i], self.OrbitParams[1][i])
            CheckParamWidget.CheckParam.setChecked(False)
            self.WidgetOrbitParams.append(CheckParamWidget)
            self.CheckLayout.addWidget(CheckParamWidget)

        self.CheckWidget = QWidget() # Container
        self.CheckWidget.setLayout(self.CheckLayout)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckWidget)

        # Histogram binning
        self.NbBins = 20
        self.NbBinsWidget = SpinBox('Number of bins', 'Number of bins', ParamDefault=self.NbBins, ParamMin=1, ParamMax=1000000)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

    # Parameters update
    def UpdateParams(self):
        self.nBody = int(self.nBodyWidget.ComboParam.currentText())-1
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        
    # Plot
    def Plot(self):
        
        # Clear axis
        for i in range(len(self.WindowPlot.WidgetPlot.Canvas.fig.axes)):
            self.WindowPlot.WidgetPlot.Canvas.fig.delaxes(self.WindowPlot.WidgetPlot.Canvas.fig.axes[0])

        # Plot initialisation
        # self.Subplot = self.WindowPlot.WidgetPlot.Canvas.fig.add_subplot(111)

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
                Data.append(eval('self.Select'+x.CheckParam.text())[self.nBody])
                DataLabels.append(str(x.CheckParam.text()))


        # Plot with current parameters
        try:
            CornerFig = corner.corner(np.array(Data).T, labels=DataLabels, fig=self.WindowPlot.WidgetPlot.Canvas.fig, bins=self.NbBins)
        except:
            print('There is a problem with the corner plot: make sure that the orbit parameters selected are variable in this adjustment.')
        
        # Plot features
        
        # Update canvas
        self.WindowPlot.WidgetPlot.Canvas.draw()



### --- Check --- ###
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    n='1'
    System = 'GGTau'
    PathInputData = '/Users/lacquema/Documents/Oracledata/'+System+f'/simu_ggtau_Ab12_'+n+'/ggtau_Ab12_hci.dat'
    PathOutputData = '/Users/lacquema/Documents/Oracledata/'+System+f'/simu_ggtau_Ab12_'+n+'/solggtauAB12.dat'
    StarDist = 145
    NbPtsEllipse = 1000
    NbSelectOrbits = 500
    AngleRange = (-3*np.pi, 3*np.pi)
    InputData = InputDataClass(PathInputData).InputData
    OutputParams = OutputDataClass(PathOutputData).OutputParams
    BestOrbits = BestOrbitsClass(*OutputParams, NbPtsEllipse, StarDist)
    BestOrbitsParams = BestOrbits.BestParams
    BestOrbitsEllipses = BestOrbits.BestEllipses
    SelectOrbits = SelectOrbitsClass(*OutputParams, NbSelectOrbits, NbPtsEllipse, StarDist)
    SelectOrbitsParams = SelectOrbits.SelectParams
    SelectOrbitsEllipses = SelectOrbits.SelectEllipses
    # ToolWidget = TempoView(InputData, SelectOrbitsEllipses, BestOrbitsEllipses)
    PosAtDateWidget = PosAtDate(SelectOrbitsEllipses, BestOrbitsEllipses)
    PosAtDateWidget.show()
    app.exec() # Application execution
