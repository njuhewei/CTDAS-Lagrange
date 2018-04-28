#! /usr/bin/env python
import sys
import calendar
import copy
from datetime import datetime, timedelta
from da.tools.general import date2num, num2date
from numpy import array, zeros, newaxis, logical_and, arange
    

def Fromdatetime(date):
    dt = date.timetuple()
    return datetime(*dt[0:6])

def increase_time(dd, **kwargs):
    """ Function increases the time by specified amount"""
    return dd + timedelta(**kwargs)

def chardate(dd, cut=8):
    return dd.strftime('%Y%m%d%H%M')[0:cut]

def timegen(sd, ed, dt):
    dd = []
    while sd <= ed:
        dd.append(Fromdatetime(sd))
        sd = sd + dt
    return dd

def itau2datetime(itau, iyear0):
    """ Function returns a datetime object from TM5s itau times"""
    date0 = datetime(iyear0, 1, 1, 0, 0, 0)
    if len(itau) == 1:
        itau = [itau] 
    for time in itau:
        sec = time % 60
        min = (time / 60) % 60
        hrs = (time / 3600) % 24
        day = (time / 86400) 
        dt = timedelta(days=day, hours=hrs, minutes=min, seconds=sec)
        yield date0 + dt 

def date2dec(date):
    """ Function converts datetime object to a Decimal number time  (e.g., 1991.875 such as in IDL, CCG) """
    if not isinstance(date, list): date = [date]

    newdate = []
    for dd in date:
        Days0 = date2num(datetime(dd.year, 1, 1))
        if calendar.isleap(dd.year):
            DaysPerYear = 366.
        else:
            DaysPerYear = 365.
        DayFrac = date2num(dd)
        newdate.append(dd.year + (DayFrac - Days0) / DaysPerYear)
    if len(newdate) == 1: return newdate[0] 
    return newdate

def dec2date(dectime):
    """ Function converts decimal time from year fraction (e.g., 1991.875 such as in IDL, CCG) to a python datetime object """
    dt = num2date(dec2num(dectime)).timetuple()
    return datetime(*dt[0:7])

def dec2num(dectime):
    """ Function converts decimal time from year fraction (e.g., 1991.875 such as in IDL, CCG) to a python decimal numtime """
    from pylab import floor, drange, num2date, date2num
    if not isinstance(dectime, list): dectime = [dectime]

    newdectime = []
    for dd in dectime:
        yr = floor(dd)
        Days0 = date2num(datetime(int(yr), 1, 1))
        if calendar.isleap(yr):
            DaysPerYear = 366.
        else:
            DaysPerYear = 365.
        DayFrac = (dd - yr) * DaysPerYear
        newdectime.append(Days0 + DayFrac)
    if len(newdectime) == 1: return newdectime[0] 
    return newdectime

def num2dec(numtime):
    """ Function converts python decimal numtime to an IDL decimal time """
    from pylab import floor, drange, num2date, date2num
    res = date2dec(num2mydate(numtime))
    return res

def num2mydate(num):
    """ Function converts decimal time from year fraction (e.g., 1991.875 such as in IDL, CCG) to a python datetime object """
    dt = num2date(num).timetuple()
    return datetime(*dt[0:7])

def monthgen(sd, ed):
    """ Generate sequence of datetime objects spaced by one month"""
    from pylab import arange
    if ed < sd: 
        raise ValueError, 'start date exceeds end date'
        sys.exit(2)
    dates = []
    for year in arange(sd.year, ed.year + 2):
        for month in arange(1, 13):
            date = datetime(year, month, 1)
            if date > ed: return dates
            else: dates.append(date)

def nextmonth(dd):
    """ Find next 1st of the month following the date dd"""

    if dd.month == 12: 
        cc = dd.replace(year=dd.year + 1)
        ee = cc.replace(month=1)
    else:
        ee = dd.replace(month=dd.month + 1)
    ff = ee.replace(day=1)
    return ff

def in_interval(start, stop, times_in):
    """ returns a list of fractions in time interval """
    times = copy.copy(times_in)

    interval = times[1] - times[0]
    times.append(times[-1] + interval)  # extend by one interval
    times_filled = [times[0] + timedelta(days=d) for d in range((times[-1] - times[0]).days)]

    b = []
    in_int = 0.0
    for t in times_filled:   # loop over days
        if t in times[1:]:   # if new interval starts
            b.append(in_int) # add previous aggregate to output
            in_int = 0.0         # reset counter
        in_int += int(logical_and(t >= start, t < stop))  # count if in interval [start,stop >
    b.append(in_int)

    if len(b) != len(times_in) : raise ValueError

    return b

def yearly_avg(time, data, sdev=False):
    """ make monthly average from array using rundat and data"""

    years = array([d.year for d in time])
    
    aa = []
    ss = []
    tt = []
    dd = time[0]
    ed = time[-1]
    while dd <= ed:
        ddnext = datetime(dd.year + 1, 1, 1)
        weights = in_interval(dd, ddnext, time)
        if len(weights) > 1:
            weights = array(weights)
            if weights.sum() > 0.0:
                weights = weights / weights.sum()
            else:
                weights = weights

            if weights.shape[0] != data.shape[0]:
                raise ValueError, 'yearly_avg has wrongly shaped weights (%d) for data of (%d)' % (weights.shape[0], data.shape[0])
                 
            sel = (weights != 0.0).nonzero()[0]     
            #print sel,array(time).take(sel),dd,ddnext
            if data.ndim == 1:
                avg_data = (weights.take(sel) * data.take(sel, axis=0)).sum(axis=0)
                std_data = (weights.take(sel) * data.take(sel, axis=0)).std(axis=0)
            elif data.ndim == 2:
                avg_data = (weights.take(sel)[:, newaxis] * data.take(sel, axis=0)).sum(axis=0).squeeze()
                std_data = (weights.take(sel)[:, newaxis] * data.take(sel, axis=0)).std(axis=0).squeeze()
            elif data.ndim == 3:
                avg_data = (weights.take(sel)[:, newaxis, newaxis] * data.take(sel, axis=0)).sum(axis=0).squeeze()
                std_data = (weights.take(sel)[:, newaxis, newaxis] * data.take(sel, axis=0)).std(axis=0).squeeze()
            else:
                raise ValueError, 'yearly_avg takes 1, 2, or 3d arrays only'
        elif len(weights) == 1:    
            avg_data = data[0]
            std_data = 0.0
        else:
            continue # next year

        aa.append(avg_data)
        ss.append(std_data)
        tt.append(datetime(dd.year, 6, 15))

        dd = ddnext

    aa = array(aa).squeeze()
    ss = array(ss).squeeze()
    time = tt
    if len(tt) == 1: 
        aa = aa.reshape(1, *aa.shape)
        ss = ss.reshape(1, *ss.shape)
    if sdev: return time, aa, ss
    else : return time, aa

def monthly_avg(time, data, sdev=False):
    """ make monthly average from array using rundat and data"""

    years = array([d.year for d in time])
    months = array([d.month for d in time])
    
    mm = []
    ss = []
    tt = []
    dd = time[0]
    ed = time[-1]

    while dd <= ed:
        ddnext = nextmonth(dd)
        weights = in_interval(dd, ddnext, time)
        if len(weights) > 1:
            weights = array(weights)
            if weights.sum() > 0.0:
                weights = weights / weights.sum()
            else:
                weights = weights

            if weights.shape[0] != data.shape[0]:
                raise ValueError, 'yearly_avg has wrongly shaped weights (%d) for data of (%d)' % (weights.shape[0], data.shape[0])
            
            sel = (weights != 0.0).nonzero()[0]     
            #print sel,array(time).take(sel),dd,nextmonth(dd)
            if data.ndim == 1:
                avg_data = (weights.take(sel) * data.take(sel, axis=0)).sum(axis=0)
                std_data = (weights.take(sel) * data.take(sel, axis=0)).std(axis=0)
            elif data.ndim == 2:
                avg_data = (weights.take(sel)[:, newaxis] * data.take(sel, axis=0)).sum(axis=0).squeeze()
                std_data = (weights.take(sel)[:, newaxis] * data.take(sel, axis=0)).std(axis=0).squeeze()
            elif data.ndim == 3:
                avg_data = (weights.take(sel)[:, newaxis, newaxis] * data.take(sel, axis=0)).sum(axis=0).squeeze()
                std_data = (weights.take(sel)[:, newaxis, newaxis] * data.take(sel, axis=0)).std(axis=0).squeeze()
            else:
                raise ValueError, 'monthly_avg takes 1, 2, or 3d arrays only'
        elif len(weights) == 1:    
            avg_data = data[0]
            std_data = 0.0
        else:
            continue # next month

        mm.append(avg_data)
        ss.append(std_data)
        tt.append(datetime(dd.year, dd.month, 15))

        dd = ddnext


    mm = array(mm).squeeze()
    ss = array(ss).squeeze()
    time = tt

    if len(tt) == 1: 
        mm = mm.reshape(-1, *mm.shape)
        ss = ss.reshape(-1, *ss.shape)

    if sdev: return time, mm, ss
    else : return time, mm

def season_avg(time, data, sdev=False):
    """ make season average from array using rundat and data"""

    seas = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]

    mm = []
    ss = []
    tt = []
    dd = time[0]
    ed = time[-1]

    while dd <= ed:
        ddmid = nextmonth(dd)
        ddnext = nextmonth(nextmonth(nextmonth(dd)))
        weights = in_interval(dd, ddnext, time)
        if len(weights) > 1:
            weights = array(weights)
            if weights.sum() > 0.0:
                weights = weights / weights.sum()
            else:
                weights = weights

            if weights.shape[0] != data.shape[0]:
                raise ValueError, 'yearly_avg has wrongly shaped weights (%d) for data of (%d)' % (weights.shape[0], data.shape[0])
            
            sel = (weights != 0.0).nonzero()[0]     
            #print sel,array(time).take(sel),dd,nextmonth(dd)
            if data.ndim == 1:
                avg_data = (weights.take(sel) * data.take(sel, axis=0)).sum(axis=0)
                std_data = (weights.take(sel) * data.take(sel, axis=0)).std(axis=0)
            elif data.ndim == 2:
                avg_data = (weights.take(sel)[:, newaxis] * data.take(sel, axis=0)).sum(axis=0).squeeze()
                std_data = (weights.take(sel)[:, newaxis] * data.take(sel, axis=0)).std(axis=0).squeeze()
            elif data.ndim == 3:
                avg_data = (weights.take(sel)[:, newaxis, newaxis] * data.take(sel, axis=0)).sum(axis=0).squeeze()
                std_data = (weights.take(sel)[:, newaxis, newaxis] * data.take(sel, axis=0)).std(axis=0).squeeze()
            else:
                raise ValueError, 'season_avg takes 1, 2, or 3d arrays only'
        elif len(weights) == 1:    
            avg_data = data[0]
            std_data = 0.0
        else:
            continue # next month

        mm.append(avg_data)
        ss.append(std_data)
        tt.append(datetime(ddmid.year, ddmid.month, 15))

        dd = ddnext


    mm = array(mm).squeeze()
    ss = array(ss).squeeze()
    time = tt

    if len(tt) == 1: 
        mm = mm.reshape(-1, *mm.shape)
        ss = ss.reshape(-1, *ss.shape)

    if sdev: return time, mm, ss
    else : return time, mm

def longterm_avg(time, data):
    """ Create long term mean """

    time_avg = num2date(date2num(time).mean())
    data_avg = data.mean(axis=0)

    return time_avg, data_avg



if __name__ == '__main__':
    #print monthgen(datetime(2000,1,1),datetime(2006,5,1))
    dd = datetime(2002, 3, 1)
    print nextmonth(dd), dd

        
         
