#! /Users/lacquema/Oracle.env/bin/python3

### --- Packages --- ###

# Transverse packages
import numpy as np
from Utils import *
from PyQt6.QtWidgets import QWidget


# Check if the file has a header
def HeaderDataIn(PathOutputData):
    with open(PathOutputData, 'r') as file:
        lines = file.readlines()
        return len(lines) > 2
    

# Open output data
def TransfertSimu(PathOutputData, UnivYN):

    with open(PathOutputData, 'r') as file:
        lines = file.readlines()

        if HeaderDataIn(PathOutputData):

            # Lecture du header
            UnivYN, DateFormat, NbPlanets, Multiplanet, RVYN, JitterYN, NbPriors, NbRanges = map(int, lines[0].split())
            current_line = 1
            IsData = list(map(int, lines[1].split()))
            current_line += 1
            JDOffset, SystDist_pc, SystDist_mas = map(float, lines[2].split())
            current_line += 1
            NbDataStarAstrom, NbDataStarRV = map(int, lines[3].split())
            current_line += 1

            # Lecture des données astrométriques de l'étoile
            if NbDataStarAstrom > 0:
                StarAstromData = np.array([list(map(float, lines[current_line + i].split())) for i in range(NbDataStarAstrom)])
                StarAstromDate, StarAstromX, StarAstromY, StarAstromdX, StarAstromdY, StarAstromCorr = StarAstromData.T
                current_line += NbDataStarAstrom
            else:
                StarAstromDate = StarAstromX = StarAstromY = StarAstromdX = StarAstromdY = StarAstromCorr = np.array([])

            # Lecture des données de vitesse radiale absolue de l'étoile
            if NbDataStarRV > 0:
                DataStarRV = np.array([list(map(float, lines[current_line + i].split())) for i in range(NbDataStarRV)])
                StarDate, StarRV, StardRV = DataStarRV.T
                current_line += NbDataStarRV
            else:
                StarDate = StarRV = StardRV = np.array([])

            # Lecture des données planétaires

            # Lecture du nombre de données astrométriques et de vitesse radiale pour chaque planète
            NbDataPlanetsAstrom = [int(float(x)) for x in lines[current_line].split()]
            current_line += 1
            NbDataPlanetsRV = [int(float(x)) for x in lines[current_line].split()]
            current_line += 1

            # Initialize lists for planetary data
            PlanetsAstromDate = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromRa = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromDec = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromdRa = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromdDec = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromCorrRaDec = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromSep = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromPa = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromdSep = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromdPa = [np.array([]) for _ in range(NbPlanets)]
            PlanetsAstromCorrSepPa = [np.array([]) for _ in range(NbPlanets)]

            PlanetsRVDate = [np.array([]) for _ in range(NbPlanets)]
            PlanetsRV = [np.array([]) for _ in range(NbPlanets)]
            PlanetsdRV = [np.array([]) for _ in range(NbPlanets)]

            for i in range(NbPlanets):
                # Données astrométriques pour la planète
                if NbDataPlanetsAstrom[i] > 0:
                    planet_astrom_data = np.array([list(map(float, lines[current_line + j].split())) for j in range(NbDataPlanetsAstrom[i])])
                    PlanetsAstromDate[i], PlanetsAstromDec[i], PlanetsAstromRa[i], PlanetsAstromdDec[i], PlanetsAstromdRa[i], PlanetsAstromCorrRaDec[i] = planet_astrom_data.T
                    current_line += NbDataPlanetsAstrom[i]

                    # Conversion des coordonnées cartésiennes en coordonnées polaires
                    PlanetsAstromSep[i], PlanetsAstromPa[i], PlanetsAstromdSep[i], PlanetsAstromdPa[i], PlanetsAstromCorrSepPa[i] = map(
                        np.array, zip(*map(cartesian_to_polar_with_errors, 
                                           PlanetsAstromDec[i], PlanetsAstromRa[i], 
                                           PlanetsAstromdDec[i], PlanetsAstromdRa[i], 
                                           PlanetsAstromCorrRaDec[i])))

                # Données de vitesse radiale pour la planète
                if NbDataPlanetsRV[i] > 0:
                    planet_rv_data = np.array([list(map(float, lines[current_line + j].split())) for j in range(NbDataPlanetsRV[i])])
                    PlanetsRVDate[i], PlanetsRV[i], PlanetsdRV[i] = planet_rv_data.T
                    current_line += NbDataPlanetsRV[i]

            # Lecture des unités de masse
            PlanetsMassUnit = [int(float(x)) for x in lines[current_line].split()]
            current_line += 1

            # Lecture des priors
            PriorsType, PriorsCoef_m, PriorsCoef_CM = [], [], []
            PriorsMassUnit, PriorsMean, PriorsSDev, PriorsMin, PriorsMax = [], [], [], [], []
            
            for _ in range(NbPriors - NbPlanets):  # The first NbPlanets Priors are positive mass for the planets
                PriorsType.append(list(map(int, map(float, lines[current_line].split()))))
                current_line += 1
                PriorsCoef_m.append(list(map(float, lines[current_line].split())))
                current_line += 1
                PriorsCoef_CM.append(list(map(float, lines[current_line].split())))
                current_line += 1
                prior = list(map(float, lines[current_line].split()))
                PriorsMassUnit.append(prior[0])
                PriorsMean.append(prior[1])
                PriorsSDev.append(prior[2])
                PriorsMin.append(prior[3])
                PriorsMax.append(prior[4])
                current_line += 1

            # Lecture des limites
            RangesMin = []
            RangesMax = []
            for i in range(NbRanges):
                Ranges = list(map(float, lines[current_line].split()))
                RangesMin.append(Ranges[0])
                RangesMax.append(Ranges[1])
                current_line += 1

            NbData = NbDataStarAstrom + NbDataStarRV
            for i in range(NbPlanets):
                NbData += NbDataPlanetsAstrom[i] + NbDataPlanetsRV[i]

            InputData = {
                "UnivYN": UnivYN,
                "DateFormat": DateFormat,
                "Multiplanet": Multiplanet,
                "RVYN": RVYN,
                "JitterYN": JitterYN,
                "NbData": NbData,
                "JDOffset": JDOffset,
                "SystDist": {
                    "pc": SystDist_pc,
                    "mas": SystDist_mas
                },
                "IsData": {
                    "PlanetsAstrom": IsData[0],
                    "StarRV": IsData[1],
                    "PlanetsRV": IsData[2],
                    "StarAstrom": IsData[3]
                },
                "Star": { # Star mass in Msun
                    "NbDataAstrom": NbDataStarAstrom,
                    "NbDataRV": NbDataStarRV,
                    "DataAstrom": {
                        "Date": StarAstromDate,
                        "X": StarAstromX,
                        "Y": StarAstromY,
                        "dX": StarAstromdX,
                        "dY": StarAstromdY,
                        "Corr": StarAstromCorr
                    },
                    "DataRV": {
                        "Date": StarDate,
                        "RV": StarRV,
                        "dRV": StardRV
                    }
                },
                "Planets": {
                    "Nb": NbPlanets,
                    "MassUnits": PlanetsMassUnit,
                    "NbDataAstrom": NbDataPlanetsAstrom,
                    "NbDataRV": NbDataPlanetsRV,
                    "DataAstrom": {
                        "Date": PlanetsAstromDate,
                        "Ra": [ra*SystDist_mas for ra in PlanetsAstromRa],
                        "Dec": [dec*SystDist_mas for dec in PlanetsAstromDec],
                        "dRa": [dra*SystDist_mas for dra in PlanetsAstromdRa],
                        "dDec": [ddec*SystDist_mas for ddec in PlanetsAstromdDec],
                        "CorrRaDec": PlanetsAstromCorrRaDec,
                        "Sep": [sep*SystDist_mas for sep in PlanetsAstromSep],
                        "Pa": PlanetsAstromPa,
                        "dSep": [dsep*SystDist_mas for dsep in PlanetsAstromdSep],
                        "dPa": PlanetsAstromdPa,
                        "CorrSepPa": PlanetsAstromCorrSepPa
                    },
                    "DataRV": {
                        "Date": PlanetsRVDate,
                        "RV": [RV*149597870.7/86400 for RV in PlanetsRV],
                        "dRV": [dRV*149597870.7/86400 for dRV in PlanetsdRV]
                    }  
                },
                "Priors": {
                    "NbPriors": NbPriors,
                    "Type": PriorsType,
                    "Coef_m": PriorsCoef_m,
                    "Coef_CM": PriorsCoef_CM,
                    "MassUnit": PriorsMassUnit,
                    "Mean": PriorsMean,
                    "SDev": PriorsSDev,
                    "Min": PriorsMin,
                    "Max": PriorsMax
                },
                "Ranges": {
                    "NbRanges": NbRanges,
                    "Min": RangesMin,
                    "Max": RangesMax
                }
            } 
            
            # Data
            NbOrbits, NbParams, NbBodies = [int(float(x)) for x in lines[current_line].split()]
            current_line += 1
            # Data = np.loadtxt(PathOutputData, skiprows = current_line, max_rows = current_line).reshape(NbBodies, NbParams, NbOrbits)
            Data = np.array(list(map(float, lines[current_line].split()))).reshape(NbBodies, NbParams, NbOrbits)

        else:

            InputData = None

            # Data
            NbOrbits, NbParams, NbBodies = [int(float(x)) for x in lines[0].split()]
            Data = np.array(list(map(float, lines[1].split()))).reshape(NbBodies, NbParams, NbOrbits)


        a, P, e, w, i, W, tp, m, V0, Jitter, m0, Chi2, Map = [np.zeros((NbBodies, NbOrbits)) for k in range(13)]

        for j in range(NbBodies):
            if NbParams == 13: # if RV data, V0 and Jitter are included
                a[j], P[j], e[j], w[j], i[j], W[j], tp[j], m[j], V0[j], Jitter[j], m0[j], Chi2[j], Map[j] = [Data[j][k][:] for k in range(13)]
            elif NbParams == 11:
                a[j], P[j], e[j], w[j], i[j], W[j], tp[j], m[j], m0[j], Chi2[j], Map[j] = [Data[j][k][:] for k in range(11)]

        # Conversions
        P = P/365.25
        i = np.rad2deg(i)
        w = np.rad2deg(w+np.pi)
        W = np.rad2deg(W+np.pi)
        W += 90.0 # to have PA=0 at North and increasing Eastward (observer's point of view)
        W = W % 360.0 # keep W between 0 and 360 degrees
        if np.mean(tp) > 2400000.5: tp = jd_to_mjd(tp)

        # print(NbDataPlanetsRV)

        if UnivYN == 2:
            a = a/(1-e) 
            P = P/((1-e)**(3/2))

        OutputParams = [NbBodies, NbOrbits, P, a, e, i, w, W, tp, m, m0, V0, Jitter, Chi2, Map]

        # for m_i in m[0]:
        #     if m_i != np.float('inf'):
        #         print(m_i)
        # print(a[0]**3 / P[0]**2)

    return InputData, OutputParams

# TransfertSimu('/Users/lacquema/Simulations/ggtau/ggtau_Ab2_Ab1_fit/ggtau_Ab2_Ab1_fit_1/solggtauAb12.dat')
# TransfertSimu('/Users/lacquema/Simulations/ggtau/ggtau_Ab2Aa_Ab1_fit/ggtau_Ab2Aa_Ab1_fit_Herve_3f/solggtaua_3f.dat')

# # Open input data
# def TransfertData(PathInputData):

#     # Header
#     with open(PathInputData, 'r') as Data:
#         DataLines = Data.readlines()
        
#         JDoffset = float(DataLines[0].split()[0])

#         NbDataAbsAstro = 0
#         NbDataRelRV = 0
#         NbDataAbsRV = 0

#         DataAbsAstro = []
#         DataRelRV = []
#         DataAbsRV = []

#         # Relative astrometric variables
#         DataRelAstro = []
#         NbDataRelAstro = 0
#         IRelAstro = []
#         JJRelAstro = []
#         MMRelAstro = []
#         YYRelAstro = []
#         MJDRelAstro = []
#         Ra = []
#         Dec = []
#         DRa = []
#         DDec = []
#         Sep = []
#         Pa = []
#         DSep = []
#         DPa = []
#         CorrDecRa = []
#         CorrSepPa = []
#         SourceRelAstro = []

#         RvAbsRv = []
#         Drv = []


#         for i in range(len(DataLines)):
#             if i>2 and DataLines[i-3][0] == '-' and DataLines[i][0] == '-' and len(DataLines[i-2].split()) == 1: # True delimitations
#                 Header = DataLines[i-1]
#                 count = 1

#                 if DataLines[i-2][0] == '1': # Relative astrometry
#                     while len(DataLines[i+count].split()) > 1:
#                         DataLine = DataLines[i+count].split()
#                         IRelAstro.append(int(DataLine[0]))
#                         SourceRelAstro.append(DataLine[-1])
#                         if 'JJ' in Header:
#                             JJRelAstro.append(int(DataLine[1]))
#                             MMRelAstro.append(int(DataLine[2]))
#                             YYRelAstro.append(int(DataLine[3]))
#                             NbVarDate = 3
#                             Conv = DatetoMJD(YYRelAstro[-1], MMRelAstro[-1], JJRelAstro[-1])
#                             MJDRelAstro.append(Conv)
#                         elif 'JD' in Header:
#                             MJDRelAstro.append(float(DataLine[1]))
#                             NbVarDate = 1
#                             Conv = MJDtoDate(MJDRelAstro[-1])
#                             JJRelAstro.append(int(Conv[2]))
#                             MMRelAstro.append(int(Conv[1]))
#                             YYRelAstro.append(int(Conv[0]))
#                         if 'RA' in Header:
#                             if len(DataLine) != (8 or 10): CorrDecRa.append(0) # no Corr
#                             else: CorrDecRa.append(float(DataLine[5+NbVarDate]))
#                             Ra.append(float(DataLine[2+NbVarDate]))
#                             Dec.append(float(DataLine[1+NbVarDate]))
#                             DRa.append(float(DataLine[4+NbVarDate]))
#                             DDec.append(float(DataLine[3+NbVarDate]))
#                             Conv = cartesian_to_polar_with_errors(Dec[-1], Ra[-1], DDec[-1], DRa[-1], CorrDecRa[-1])
#                             Sep.append(Conv[0])
#                             Pa.append(Conv[1])
#                             DSep.append(Conv[2])
#                             DPa.append(Conv[3])
#                             CorrSepPa.append(Conv[4])
#                         elif 'SEP' in Header:
#                             if len(DataLine) != (8 or 10): CorrSepPa.append(0) # no Corr
#                             else: CorrSepPa.append(float(DataLine[5+NbVarDate]))
#                             Sep.append(float(DataLine[1+NbVarDate]))
#                             Pa.append(float(DataLine[2+NbVarDate]))
#                             DSep.append(float(DataLine[3+NbVarDate]))
#                             DPa.append(float(DataLine[4+NbVarDate]))
#                             Conv = polar_to_cartesian_with_errors(Sep[-1], Pa[-1], DSep[-1], DPa[-1], CorrSepPa[-1])
#                             Dec.append(Conv[0])
#                             Ra.append(Conv[1])
#                             DDec.append(Conv[2])
#                             DRa.append(Conv[3])
#                             CorrDecRa.append(Conv[4])
#                         NbDataRelAstro += 1
#                         count += 1
#                         if i+count == len(DataLines): break
                
#                 # if DataLines[i-2][0] == '2': # Absolute RV
#                 #     while len(DataLines[i+count].split()) > 1:
#                 #         DataAbsRV.append(DataLines[i+count].split())
#                 #         NbDataAbsRV += 1
#                 #         count += 1
#                 #         if i+count == len(DataLines): break
                
#                 # if DataLines[i-2][0] == '3': # Relative RV
#                 #     while len(DataLines[i+count].split()) > 1:
#                 #         DataRelRV.append(DataLines[i+count].split())
#                 #         NbDataRelRV += 1
#                 #         count += 1
#                 #         if i+count == len(DataLines): break
                
#                 # if DataLines[i-2][0] == '4': # Absolute astrometry (not yet)
#                 #         return

#         DataRelAstro = [NbDataRelAstro, IRelAstro, MJDRelAstro, JJRelAstro, MMRelAstro, YYRelAstro, Ra, Dec, DRa, DDec, CorrDecRa, Sep, Pa, DSep, DPa, CorrSepPa, SourceRelAstro] # format I, MJD, JJ, MM, YY, Dec, Ra, DDec, DRa, Sep, Pa, DSep, DPa, Corr, Source

#         # if len(DataAbsRV) !=0 :
#         #     DataAbsRV = np.transpose(DataAbsRV)
#         #     MJDAbsRV = np.flip(AllDate2jd(*np.float64(np.flip(DataAbsRV[1:4]))))
#         #     DataAbsRV = [NbDataAbsRV, MJDAbsRV, *np.float64(DataAbsRV[:-1]), DataAbsRV[-1]] # format I, MJD, JJ, MM, YY, Rv, DRv, Source

#         # if len(DataRelRV) !=0 :
#         #     DataRelRV = np.transpose(DataRelRV)
#         #     MJDRelRV = np.flip(AllDate2jd(*np.float64(np.flip(DataRelRV[1:4]))))
#         #     DataRelRV = [NbDataAbsRV, MJDRelRV, *np.float64(DataRelRV[:-1]), DataRelRV[-1]] # format I, MJD, JJ, MM, YY, Rv, DRv, Source

#         return JDoffset, DataRelAstro #, DataAbsRV, DataRelRV



# # TransfertSimu('/Users/lacquema/Simulations/ggtau/ggtau_Ab2Aa_Ab1_fit/ggtau_Ab2Aa_Ab1_fit_Herve_3f/solggtaua_3f.dat')
# # TransfertData('/Users/lacquema/Oracle.env/Oracle/Algorithm/Tests/ggtau_triple.dat')
# # TransfertData('/Users/lacquema/Oracle.env/Simulations/BetaPic/simu_bpic_b_1/bpicb_hci.dat')

