#! /Users/lacquema/ByeGildas/bin/python3

import numpy as np
from TimeConvertor import date_to_jd, jd_to_mjd
import TransferData
import BestFit

def OpenOutputData(PathOutputData):

    # Header
    Header = np.loadtxt(PathOutputData, max_rows = 1, dtype = int)
    NbOrbits, NbParams, NbBodies = Header

    # Data
    Data = np.loadtxt(PathOutputData, skiprows = 1, max_rows = 2).reshape(NbBodies, NbParams, NbOrbits)

    a, P, e, w, i, W, tp, m, Mdyn, Chi2, map = [np.zeros((NbBodies, NbOrbits)) for k in range(11)]

    for j in range(NbBodies):
        a[j], P[j], e[j], w[j], i[j], W[j], tp[j], m[j], Mdyn[j], Chi2[j], map[j] = [Data[j][k][:] for k in range(11)]

    return NbBodies, NbOrbits, P/365.25, a, e, np.rad2deg(i), np.rad2deg(w+np.pi), np.rad2deg(W+np.pi), jd_to_mjd(tp), m, Mdyn, Chi2, map