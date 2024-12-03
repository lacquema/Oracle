#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import numpy as np
import math
import datetime as dt


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
    


# # Conversion (Sep, PA) to (Ra, Dec) en mas
# def SEPPAtoRADEC(Sep, Pa, DSep, DPa):

#     Pa = np.deg2rad(Pa) # Conversion de Pa en rad
#     # DPa = DPa/3600/1000 # Conversion de DPa en mas
#     DPa = np.deg2rad(DPa)

#     Ra = Sep*np.sin(Pa)
#     DRa = np.sqrt((np.sin(Pa)*DSep)**2 + (Sep*np.cos(Pa)*DPa)**2)

#     Dec = Sep*np.cos(Pa)
#     DDec = np.sqrt((np.cos(Pa)*DSep)**2 + (Sep*np.sin(Pa)*DPa)**2)

#     return Ra, Dec, DRa, DDec



def RADECtoSEPPA(Ra, Dec, DRa, DDec, Corr): # Ra [mas] and Dec [mas]

    if type(Ra)=='int' or 'float':
        NbData = 1
        Ra = [Ra]
        Dec = [Dec]
        DRa = [DRa]
        DDec = [DDec]
    else:
        NbData = np.shape(Ra)[0]

    if Corr==None:
        Corr = 0

    # Initialisation outputs
    Sep, Pa, DSep, DPa = [np.zeros(NbData) for i in range(4)]
    
    # Conversions inputs
    Ra = np.array(Ra)
    Dec = np.array(Dec)
    DRa = np.array(DRa)
    DDec = np.array(DDec)

    for i in range(NbData):

        V = np.array([[DDec[i]**2, DDec[i]*DRa[i]*Corr],
                    [DDec[i]*DRa[i]*Corr, DRa[i]**2]])

        Sep[i] = np.sqrt(Dec[i]**2 + Ra[i]**2)
        gSep = np.array([Dec[i]/np.sqrt(Dec[i]**2 + Ra[i]**2), Ra[i]/np.sqrt(Dec[i]**2 + Ra[i]**2)])
        DSep[i] = np.sqrt(np.matmul(np.transpose(gSep), np.matmul(V, np.transpose(gSep))))
        
        Pa[i] = np.arctan(Ra[i]/Dec[i])
        gPa = np.array([1/(Dec[i]+Ra[i]**2/Dec[i]), -1/(Ra[i]+Dec[i]**2/Ra[i])])
        DPa[i] = np.sqrt(np.matmul(np.transpose(gPa), np.matmul(V, np.transpose(gPa))))

    # Conversions outputs
    Pa = np.rad2deg(Pa)
    DPa = np.rad2deg(DPa)
    
    return Sep, Pa, DSep, DPa # Sep [mas] and Pa [deg]




def SEPPAtoRADEC(Sep, Pa, DSep, DPa, Corr=None): # Sep [mas] and Pa [deg]

    if type(Sep)=='int' or 'float':
        NbData = 1
        Sep = [Sep]
        Pa = [Pa]
        DSep = [DSep]
        DPa = [DPa]
    else:
        NbData = np.shape(Sep)[0]

    if Corr==None:
        Corr = 0

    # Initialisation outputs
    Ra, Dec, DRa, DDec = [np.zeros(NbData) for i in range(4)]
    
    # Conversions inputs
    Sep = np.array(Sep)
    Pa = np.deg2rad(np.array(Pa))
    DSep = np.array(DSep)
    DPa = np.deg2rad(np.array(DPa))

    for i in range(NbData):

        V = np.array([[DSep[i]**2, DSep[i]*DPa[i]*Corr],
                      [DSep[i]*DPa[i]*Corr, DPa[i]**2]])
        
        Ra[i] = Sep[i]*np.sin(Pa[i])
        gRa = np.array([np.sin(Pa[i]), Sep[i]*np.cos(Pa[i])])
        DRa[i] = np.sqrt(np.matmul(np.transpose(gRa), np.matmul(V, np.transpose(gRa))))

        Dec[i] = Sep[i]*np.cos(Pa[i])
        gDec = np.array([np.cos(Pa[i]), -Sep[i]*np.sin(Pa[i])])
        DDec[i] = np.sqrt(np.matmul(np.transpose(gDec), np.matmul(V, np.transpose(gDec))))
    
    return Ra, Dec, DRa, DDec # Dep [mas] and Ra [mas]


    

# Convert all date to jd
def AllDate2jd(YY, MM, JJ):

    if type(YY)=='int' or 'float':
        NbDate = 1
        YY = [YY]
        MM = [MM]
        JJ = [JJ]
    else:
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


def MJDtoDate(MJD):
    JD = mjd_to_jd(MJD)
    YY, MM, JJ = jd_to_date(JD)
    return int(YY), int(MM), int(JJ)


def DatetoMJD(YY, MM, JJ):
    JD = date_to_jd(YY, MM, JJ)
    MJD = jd_to_mjd(JD)
    return int(MJD)


def mjd_to_jd(mjd):
    """
    Convert Modified Julian Day to Julian Day.
        
    Parameters
    ----------
    mjd : float
        Modified Julian Day
        
    Returns
    -------
    jd : float
        Julian Day
    
        
    """
    return mjd + 2400000.5

    
def jd_to_mjd(jd):
    """
    Convert Julian Day to Modified Julian Day
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    mjd : float
        Modified Julian Day
    
    """
    return jd - 2400000.5

    
def date_to_jd(year,month,day):
    """
    Convert a date to Julian Day.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
    
    Returns
    -------
    jd : float
        Julian Day
        
    Examples
    --------
    Convert 6 a.m., February 17, 1985 to Julian Day
    
    >>> date_to_jd(1985,2,17.25)
    2446113.75
    
    """
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month
    
    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
        (year == 1582 and month < 10) or
        (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)
        
    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)

    D = math.trunc(30.6001 * (monthp + 1))
    
    jd = B + C + D + day + 1720994.5
    
    return jd
    
    
def jd_to_date(jd):
    """
    Convert Julian Day to date.
    
    Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
        4th ed., Duffet-Smith and Zwart, 2011.
    
    Parameters
    ----------
    jd : float
        Julian Day
        
    Returns
    -------
    year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
        
    month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
    
    day : float
        Day, may contain fractional part.
        
    Examples
    --------
    Convert Julian Day 2446113.75 to year, month, and day.
    
    >>> jd_to_date(2446113.75)
    (1985, 2, 17.25)
    
    """
    jd = jd + 0.5
    
    F, I = math.modf(jd)
    I = int(I)
    
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    
    D = math.trunc((C - 122.1) / 365.25)
    
    E = math.trunc(365.25 * D)
    
    G = math.trunc((C - E) / 30.6001)
    
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715
        
    return year, month, day
    
    
def hmsm_to_days(hour=0,min=0,sec=0,micro=0):
    """
    Convert hours, minutes, seconds, and microseconds to fractional days.
    
    Parameters
    ----------
    hour : int, optional
        Hour number. Defaults to 0.
    
    min : int, optional
        Minute number. Defaults to 0.
    
    sec : int, optional
        Second number. Defaults to 0.
    
    micro : int, optional
        Microsecond number. Defaults to 0.
        
    Returns
    -------
    days : float
        Fractional days.
        
    Examples
    --------
    >>> hmsm_to_days(hour=6)
    0.25
    
    """
    days = sec + (micro / 1.e6)
    
    days = min + (days / 60.)
    
    days = hour + (days / 60.)
    
    return days / 24.
    
    
def days_to_hmsm(days):
    """
    Convert fractional days to hours, minutes, seconds, and microseconds.
    Precision beyond microseconds is rounded to the nearest microsecond.
    
    Parameters
    ----------
    days : float
        A fractional number of days. Must be less than 1.
        
    Returns
    -------
    hour : int
        Hour number.
    
    min : int
        Minute number.
    
    sec : int
        Second number.
    
    micro : int
        Microsecond number.
        
    Raises
    ------
    ValueError
        If `days` is >= 1.
        
    Examples
    --------
    >>> days_to_hmsm(0.1)
    (2, 24, 0, 0)
    
    """
    hours = days * 24.
    hours, hour = math.modf(hours)
    
    mins = hours * 60.
    mins, min = math.modf(mins)
    
    secs = mins * 60.
    secs, sec = math.modf(secs)
    
    micro = round(secs * 1.e6)
    
    return int(hour), int(min), int(sec), int(micro)
    

def datetime_to_jd(date):
    """
    Convert a `datetime.datetime` object to Julian Day.
    
    Parameters
    ----------
    date : `datetime.datetime` instance
    
    Returns
    -------
    jd : float
        Julian day.
        
    Examples
    --------
    >>> d = datetime.datetime(1985,2,17,6)  
    >>> d
    datetime.datetime(1985, 2, 17, 6, 0)
    >>> jdutil.datetime_to_jd(d)
    2446113.75
    
    """
    days = date.day + hmsm_to_days(date.hour,date.minute,date.second,date.microsecond)
    
    return date_to_jd(date.year,date.month,days)
    
    
def jd_to_datetime(jd):
    """
    Convert a Julian Day to an `jdutil.datetime` object.
    
    Parameters
    ----------
    jd : float
        Julian day.
        
    Returns
    -------
    dt : `jdutil.datetime` object
        `jdutil.datetime` equivalent of Julian day.
    
    Examples
    --------
    >>> jd_to_datetime(2446113.75)
    datetime(1985, 2, 17, 6, 0)
    
    """
    year, month, day = jd_to_date(jd)
    
    frac_days,day = math.modf(day)
    day = int(day)
    
    hour,min,sec,micro = days_to_hmsm(frac_days)
    
    return datetime(year,month,day,hour,min,sec,micro)


def timedelta_to_days(td):
    """
    Convert a `datetime.timedelta` object to a total number of days.
    
    Parameters
    ----------
    td : `datetime.timedelta` instance
    
    Returns
    -------
    days : float
        Total number of days in the `datetime.timedelta` object.
        
    Examples
    --------
    >>> td = datetime.timedelta(4.5)
    >>> td
    datetime.timedelta(4, 43200)
    >>> timedelta_to_days(td)
    4.5
    
    """
    seconds_in_day = 24. * 3600.
    
    days = td.days + (td.seconds + (td.microseconds * 10.e6)) / seconds_in_day
    
    return days
    
    
class datetime(dt.datetime):
    """
    A subclass of `datetime.datetime` that performs math operations by first
    converting to Julian Day, then back to a `jdutil.datetime` object.
    
    Addition works with `datetime.timedelta` objects, subtraction works with
    `datetime.timedelta`, `datetime.datetime`, and `jdutil.datetime` objects.
    Not all combinations work in all directions, e.g.
    `timedelta - datetime` is meaningless.
    
    See Also
    --------
    datetime.datetime : Parent class.
    
    """
    def __add__(self,other):
        if not isinstance(other,dt.timedelta):
            s = "jdutil.datetime supports '+' only with datetime.timedelta"
            raise TypeError(s)
        
        days = timedelta_to_days(other)
        
        combined = datetime_to_jd(self) + days
        
        return jd_to_datetime(combined)
        
    def __radd__(self,other):
        if not isinstance(other,dt.timedelta):
            s = "jdutil.datetime supports '+' only with datetime.timedelta"
            raise TypeError(s)
        
        days = timedelta_to_days(other)
        
        combined = datetime_to_jd(self) + days
        
        return jd_to_datetime(combined)
        
    def __sub__(self,other):
        if isinstance(other,dt.timedelta):
            days = timedelta_to_days(other)
            
            combined = datetime_to_jd(self) - days
            
            return jd_to_datetime(combined)
            
        elif isinstance(other, (datetime,dt.datetime)):
            diff = datetime_to_jd(self) - datetime_to_jd(other)
            
            return dt.timedelta(diff)
            
        else:
            s = "jdutil.datetime supports '-' with: "
            s += "datetime.timedelta, jdutil.datetime and datetime.datetime"
            raise TypeError(s)
            
    def __rsub__(self,other):
        if not isinstance(other, (datetime,dt.datetime)):
            s = "jdutil.datetime supports '-' with: "
            s += "jdutil.datetime and datetime.datetime"
            raise TypeError(s)
            
        diff = datetime_to_jd(other) - datetime_to_jd(self)
            
        return dt.timedelta(diff)
        

    def to_jd(self):
        """
        Return the date converted to Julian Day.
        
        """
        return datetime_to_jd(self)
        
    def to_mjd(self):
        """
        Return the date converted to Modified Julian Day.
        
        """
        return jd_to_mjd(self.to_jd())