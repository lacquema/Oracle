#! /Users/lacquema/Oracle.env/bin/python3

### --- Packages --- ###

# Transverse packages
import numpy as np
from Utils import *
from PyQt6.QtWidgets import QWidget


### --- Plot Window Generating --- ###

# Open output data
def TransfertSimu(PathOutputData):

    # Header
    Header = np.loadtxt(PathOutputData, max_rows = 1, dtype = int)
    NbOrbits, NbParams, NbBodies = Header

    # Data
    Data = np.loadtxt(PathOutputData, skiprows = 1, max_rows = 2).reshape(NbBodies, NbParams, NbOrbits)

    a, P, e, w, i, W, tp, m, Mdyn, Chi2, map = [np.zeros((NbBodies, NbOrbits)) for k in range(11)]

    for j in range(NbBodies):
        a[j], P[j], e[j], w[j], i[j], W[j], tp[j], m[j], Mdyn[j], Chi2[j], map[j] = [Data[j][k][:] for k in range(11)]

    # Conversions
    P = P/365.25
    i = np.rad2deg(i)
    w = np.rad2deg(w+np.pi)
    W = np.rad2deg(W+np.pi)
    tp = jd_to_mjd(tp)

    return NbBodies, NbOrbits, P, a, e, i, w, W, tp, m, Mdyn, Chi2, map


# Open input data
def TransfertData(PathInputData):

    # Header
    with open(PathInputData, 'r') as Data:
        DataLines = Data.readlines()
        
        JDoffset = float(DataLines[0].split()[0])

        NbDataAbsAstro = 0
        NbDataRelRV = 0
        NbDataAbsRV = 0

        DataAbsAstro = []
        DataRelRV = []
        DataAbsRV = []

        # Relative astrometric variables
        DataRelAstro = []
        NbDataRelAstro = 0
        IRelAstro = []
        JJRelAstro = []
        MMRelAstro = []
        YYRelAstro = []
        MJDRelAstro = []
        Ra = []
        Dec = []
        DRa = []
        DDec = []
        Sep = []
        Pa = []
        DSep = []
        DPa = []
        CorrDecRa = []
        CorrSepPa = []
        SourceRelAstro = []

        RvAbsRv = []
        Drv = []


        for i in range(len(DataLines)):
            if i>2 and DataLines[i-3][0] == '-' and DataLines[i][0] == '-' and len(DataLines[i-2].split()) == 1: # True delimitations
                Header = DataLines[i-1]
                count = 1

                if DataLines[i-2][0] == '1': # Relative astrometry
                    while len(DataLines[i+count].split()) > 1:
                        DataLine = DataLines[i+count].split()
                        IRelAstro.append(int(DataLine[0]))
                        SourceRelAstro.append(DataLine[-1])
                        if 'JJ' in Header:
                            JJRelAstro.append(int(DataLine[1]))
                            MMRelAstro.append(int(DataLine[2]))
                            YYRelAstro.append(int(DataLine[3]))
                            NbVarDate = 3
                            Conv = DatetoMJD(YYRelAstro[-1], MMRelAstro[-1], JJRelAstro[-1])
                            MJDRelAstro.append(Conv)
                        elif 'JD' in Header:
                            MJDRelAstro.append(float(DataLine[1]))
                            NbVarDate = 1
                            Conv = MJDtoDate(MJDRelAstro[-1])
                            JJRelAstro.append(int(Conv[2]))
                            MMRelAstro.append(int(Conv[1]))
                            YYRelAstro.append(int(Conv[0]))
                        if 'RA' in Header:
                            if len(DataLine) != (8 or 10): CorrDecRa.append(0) # no Corr
                            else: CorrDecRa.append(float(DataLine[5+NbVarDate]))
                            Ra.append(float(DataLine[2+NbVarDate]))
                            Dec.append(float(DataLine[1+NbVarDate]))
                            DRa.append(float(DataLine[4+NbVarDate]))
                            DDec.append(float(DataLine[3+NbVarDate]))
                            Conv = cartesian_to_polar_with_errors(Dec[-1], Ra[-1], DDec[-1], DRa[-1], CorrDecRa[-1])
                            Sep.append(Conv[0])
                            Pa.append(Conv[1])
                            DSep.append(Conv[2])
                            DPa.append(Conv[3])
                            CorrSepPa.append(Conv[4])
                        elif 'SEP' in Header:
                            if len(DataLine) != (8 or 10): CorrSepPa.append(0) # no Corr
                            else: CorrSepPa.append(float(DataLine[5+NbVarDate]))
                            Sep.append(float(DataLine[1+NbVarDate]))
                            Pa.append(float(DataLine[2+NbVarDate]))
                            DSep.append(float(DataLine[3+NbVarDate]))
                            DPa.append(float(DataLine[4+NbVarDate]))
                            Conv = polar_to_cartesian_with_errors(Sep[-1], Pa[-1], DSep[-1], DPa[-1], CorrSepPa[-1])
                            Dec.append(Conv[0])
                            Ra.append(Conv[1])
                            DDec.append(Conv[2])
                            DRa.append(Conv[3])
                            CorrDecRa.append(Conv[4])
                        NbDataRelAstro += 1
                        count += 1
                        if i+count == len(DataLines): break
                
                # if DataLines[i-2][0] == '2': # Absolute RV
                #     while len(DataLines[i+count].split()) > 1:
                #         DataAbsRV.append(DataLines[i+count].split())
                #         NbDataAbsRV += 1
                #         count += 1
                #         if i+count == len(DataLines): break
                
                # if DataLines[i-2][0] == '3': # Relative RV
                #     while len(DataLines[i+count].split()) > 1:
                #         DataRelRV.append(DataLines[i+count].split())
                #         NbDataRelRV += 1
                #         count += 1
                #         if i+count == len(DataLines): break
                
                # if DataLines[i-2][0] == '4': # Absolute astrometry (not yet)
                #         return

        DataRelAstro = [NbDataRelAstro, IRelAstro, MJDRelAstro, JJRelAstro, MMRelAstro, YYRelAstro, Ra, Dec, DRa, DDec, CorrDecRa, Sep, Pa, DSep, DPa, CorrSepPa, SourceRelAstro] # format I, MJD, JJ, MM, YY, Dec, Ra, DDec, DRa, Sep, Pa, DSep, DPa, Corr, Source

        # if len(DataAbsRV) !=0 :
        #     DataAbsRV = np.transpose(DataAbsRV)
        #     MJDAbsRV = np.flip(AllDate2jd(*np.float64(np.flip(DataAbsRV[1:4]))))
        #     DataAbsRV = [NbDataAbsRV, MJDAbsRV, *np.float64(DataAbsRV[:-1]), DataAbsRV[-1]] # format I, MJD, JJ, MM, YY, Rv, DRv, Source

        # if len(DataRelRV) !=0 :
        #     DataRelRV = np.transpose(DataRelRV)
        #     MJDRelRV = np.flip(AllDate2jd(*np.float64(np.flip(DataRelRV[1:4]))))
        #     DataRelRV = [NbDataAbsRV, MJDRelRV, *np.float64(DataRelRV[:-1]), DataRelRV[-1]] # format I, MJD, JJ, MM, YY, Rv, DRv, Source

        return JDoffset, DataRelAstro #, DataAbsRV, DataRelRV




# TransfertData('/Users/lacquema/Oracle.env/Oracle/Algorithm/Tests/ggtau_triple.dat')
# TransfertData('/Users/lacquema/Oracle.env/Simulations/BetaPic/simu_bpic_b_1/bpicb_hci.dat')

      