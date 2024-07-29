#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import numpy as np
from TimeConvertor import date_to_jd, jd_to_mjd

# My packages



def Ellipse(P, a, e, i, w, W, tp, NbPtsEllipse=100, Period=False, Time=False, AnomT=False, AnomE=False):
    '''
    Make and 3D ellipse from orbit parameters
    
    input:
    P, period [yr]
    a, semi-major axis [AU]
    e, eccentricity
    i, inclinaison [deg]
    w, argument of periastron [deg]
    W, longitude of ascending node [deg]
    tp, periastron time passage [MJD]
    NbPts=100, number of points
    
    return:
    teta, if AnomT=False, true anomaly along ellipse
    E, if AnomE=True, eccentric anomaly along ellipse
    t [MJD], if Time=True, time along ellipse
    xR [AU], x after rotation
    yR [AU], y after rotation
    zR [AU], z after rotation
    '''

    # Convertion
    i = np.deg2rad(i)
    w = np.deg2rad(w)
    W = np.deg2rad(W)+np.pi/2
    P*=365.25

    Output = []

    if Period==True: Output.append(P)

    # True anomaly
    if AnomT==True:
        teta = np.linspace(-np.pi, np.pi, NbPtsEllipse)
        Output.append(teta)
        E = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(teta/2))
    else: 
        E = np.linspace(-np.pi, np.pi, NbPtsEllipse)

    # Eccentric anomaly
    if AnomE == True: Output.append(E)

    # Time
    if Time==True:
        t = P/(2*np.pi)*(E-e*np.sin(E))+tp
        Output.append(t)

    # Ellipse on R0
    x = a*(np.cos(E)-e)
    y = a*np.sqrt(1-e**2)*np.sin(E)
    z = np.zeros(NbPtsEllipse)

    # Change of referential
    Mp = RotMatrix(w, 'z')@RotMatrix(i, 'x')@RotMatrix(W, 'z')
    
    xR = x*Mp[0,0] + y*Mp[1,0] + z*Mp[2,0]
    Output.append(xR)
    yR = x*Mp[0,1] + y*Mp[1,1] + z*Mp[2,1]
    Output.append(yR)
    zR = x*Mp[0,2] + y*Mp[1,2] + z*Mp[2,2]
    Output.append(zR)

    return Output


 


def RotMatrix(AngleRot, Axis='x'or'y'or'z'):
    '''
    Return rotation matrix associed a rotation of AngleRot around the Axis.
    AngleRot, rotation angle [rad]
    Axis, rotation axis ('x', 'y' or 'z')
    '''
    if Axis == 'x':
        M = np.array([[1, 0, 0],
                        [0, np.cos(AngleRot), np.sin(AngleRot)],
                        [0, -np.sin(AngleRot), np.cos(AngleRot)]])
        
    if Axis == 'y':
        M = np.array([[np.cos(AngleRot), 0, -np.sin(AngleRot)],
                        [0, 1, 0],
                        [np.sin(AngleRot), 0, np.cos(AngleRot)]])
        
    if Axis == 'z':
        M = np.array([[np.cos(AngleRot), np.sin(AngleRot), 0],
                        [-np.sin(AngleRot), np.cos(AngleRot), 0],
                        [0, 0, 1]])
    
    return M
    


# Conversion (Sep, PA) to (Ra, Dec)
def SEPPAtoRADEC(Sep, Pa, DSep, DPa):

    Pa = np.deg2rad(Pa) # Conversion de Pa en rad
    DPa = DPa/3600/1000 # Conversion de DPa en mas

    Ra = Sep*np.sin(Pa)
    DRa = np.sqrt((np.sin(Pa)*DSep)**2 + (Sep*np.cos(Pa)*DPa)**2)

    Dec = Sep*np.cos(Pa)
    DDec = np.sqrt((np.cos(Pa)*DSep)**2 + (Sep*np.sin(Pa)*DPa)**2)

    return Dec, Ra, DDec, DRa


# Convert all date to jd
def AllDate2jd(YY, MM, JJ):
    NbDate = np.shape(YY)[0]
    MJD = np.zeros(NbDate)
    for k in range(NbDate):
        MJD[k] = jd_to_mjd(date_to_jd(YY[k], MM[k], JJ[k]))
    return MJD



def DelAllWidgetsBtw(Layout, indexMin, indexMax):
        for i in reversed(range(indexMin, indexMax)): 
            WidgetToRemove = Layout.itemAt(i).widget()
            Layout.removeWidget(WidgetToRemove)
            # print(self.Layout.count())
            WidgetToRemove.setParent(None)
