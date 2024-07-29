


### --- Packages --- ###

# Transverse packages
import numpy as np
from UtilsAnaSimu import SEPPAtoRADEC, AllDate2jd
from TimeConvertor import jd_to_mjd
from PyQt6.QtWidgets import QWidget


### --- Plot Window Generating --- ###

# Open output data
class OutputDataClass(QWidget):
    def __init__(self, PathOutputData):
        super().__init__()

        # Header
        Header = np.loadtxt(PathOutputData, max_rows = 1, dtype = int)
        self.NbOrbits, self.NbParams, self.NbBodies = Header

        # Data
        Data = np.loadtxt(PathOutputData, skiprows = 1, max_rows = 2).reshape(self.NbBodies, self.NbParams, self.NbOrbits)

        self.a, self.P, self.e, self.w, self.i, self.W, self.tp, self.m, self.Mdyn, self.Chi2, self.map = [np.zeros((self.NbBodies, self.NbOrbits)) for k in range(11)]

        for j in range(self.NbBodies):
            self.a[j], self.P[j], self.e[j], self.w[j], self.i[j], self.W[j], self.tp[j], self.m[j], self.Mdyn[j], self.Chi2[j], self.map[j] = [Data[j][k][:] for k in range(11)]

        # Conversions
        self.P = self.P/365.25
        self.i = np.rad2deg(self.i)
        self.w = np.rad2deg(self.w+np.pi)
        self.W = np.rad2deg(self.W+np.pi)
        self.tp = jd_to_mjd(self.tp)

        self.OutputParams = [self.NbBodies, self.NbOrbits, self.P, self.a, self.e, self.i, self.w, self.W, self.tp, self.m, self.Mdyn, self.Chi2, self.map]


# Open input data
class InputDataClass(QWidget):
    def __init__(self, PathInputData):
        super().__init__()

        # Header
        Header = np.loadtxt(PathInputData, max_rows = 1, dtype = str)

        # Data
        Data = np.loadtxt(PathInputData, skiprows = 2, dtype = float, usecols=range(len(Header)-1), unpack=True)
        self.Source = np.loadtxt(PathInputData, skiprows = 2, dtype = str, usecols=len(Header)-1)

        self.NbInputData = np.shape(Data)[1]

        if 'RA(mas)' in Header:
            self.I, self.JJ, self.MM, self.YY, self.Dec, self.Ra, self.DDec, self.DRa, self.Corr = Data

        elif 'SEP(mas)' in Header:
            self.I, self.JJ, self.MM, self.YY, self.Sep, self.Pa, self.DSep, self.DPa, self.Corr = Data

            self.Dec, self.Ra, self.DDec, self.DRa = SEPPAtoRADEC(self.Sep, self.Pa, self.DSep, self.DPa) # conversion (Sep, PA) to (Ra, Dec)

        self.MJD = AllDate2jd(self.YY, self.MM, self.JJ)

        self.InputData = [self.NbInputData, self.I, self.MJD, self.JJ, self.MM, self.YY, self.Dec, self.Ra, self.DDec, self.DRa, self.Corr, self.Source]


