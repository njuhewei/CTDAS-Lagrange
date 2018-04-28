#!/usr/bin/env python
# tools_da.py

"""
Author : peters 

Revision History:
File created on 03 Oct 2008.

Temporary module to hold classes and methods that are in development

"""

import logging
import os
import shutil
import datetime
import re

from dateutil.rrule import rrule, MO, TU, WE, TH, FR, SA, SU, YEARLY, \
     MONTHLY, WEEKLY, DAILY, HOURLY, MINUTELY, SECONDLY

HOURS_PER_DAY = 24.
MINUTES_PER_DAY  = 60.*HOURS_PER_DAY
SECONDS_PER_DAY =  60.*MINUTES_PER_DAY
MUSECONDS_PER_DAY = 1e6*SECONDS_PER_DAY
SEC_PER_MIN = 60
SEC_PER_HOUR = 3600
SEC_PER_DAY = SEC_PER_HOUR * 24
SEC_PER_WEEK = SEC_PER_DAY * 7
MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY = (
    MO, TU, WE, TH, FR, SA, SU)
WEEKDAYS = (MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY)


def validate_rc(rcfile, needed_items):
    """ validate the contents of an rc-file given a dictionary of required keys """

    for k, v in rcfile.iteritems():
        if v == 'True' :
            rcfile[k] = True
        if v == 'False':
            rcfile[k] = False
        if 'date' in k: 
            rcfile[k] = datetime.datetime.strptime(v, '%Y-%m-%d %H:%M:%S')

    for key in needed_items:
        if not rcfile.has_key(key):
            msg = 'Missing a required value in rc-file : %s' % key
            logging.error(msg)
            raise IOError, msg
    logging.debug('rc-file has been validated succesfully')  


def create_dirs(dirname, forceclean=False):
    """ Create a directory and report success, only if non-existent """

    if forceclean:
        try:
            shutil.rmtree(dirname)
        except:
            pass

    if not os.path.exists(dirname):
        os.makedirs(dirname)
        logging.info('Creating new directory %s' % dirname)
    else:
        logging.debug('Using existing directory %s' % dirname)
    return dirname


def advance_time(time_in, interval):
    """ Advance time_in by a specified interval"""

    time_out = time_in

    if interval == 'month':                       # if monthly, this run will go to the first day of the next month
        if time_in.month != 12: 
            time_out = datetime.datetime(time_in.year, time_in.month + 1, 1, time_in.hour, 0, 0)
        else: 
            time_out = datetime.datetime(time_in.year + 1, 1, 1, time_in.hour, 0, 0)  # end of year provision
    elif interval == 'week':
        time_out = time_in + datetime.timedelta(days=7)
    elif isinstance(interval, datetime.timedelta):
        time_out = time_in + interval
    else:                    # assume that the interval specified is the number of days to run forward before resubmitting
        time_out = time_in + datetime.timedelta(days=float(interval))

    return time_out



def to_datetime(datestring, fmt=None):
    """ convert a date string to a datetime object """
  
    if fmt == 'TM5':
        datestring = '%04s-%02s-%02s %02s:%02s:00' % (datestring[0:4], datestring[4:6], datestring[6:8], datestring[8:10], datestring[10:12])
    elif fmt == 'pycasso-TM5':
        pass # Format already compatible
    else:
        pass 

    try:
        return datetime.datetime.strptime(datestring, '%Y-%m-%d %H:%M:%S')
    except:
        date, time = datestring.split(' ')
        year, month, day = map(int, date.split('-'))
        hour, minute, second = map(int, time.split(':'))
        return datetime.datetime(year, month, day, hour, minute, second)


def name_convert(name=None, to=None):
    """ Convert between old GLOBALVIEW-style and new ObsPack-style

      print name_convert(name="lef_surface-pfp_1", to='GV' )
            lef_01P0

      print name_convert(name="hun_35C3", to='OP' )
            hun_tower-insitu_35
   """

    identifier = 'name_convert'

    if name == None or to == None: 
        return ""
    
    platform_dict = { 'surface':0, 'shipboard':1, 'aircraft':2, 'tower':3 }
    strategy_dict = { 'flask':'D', 'pfp':'P', 'insitu':'C' }

    #----------------------------------->>>
    if to.upper() == "GV":
        # alt_surface-flask_1 -> alt_01D0
        try:
            fields = name.split('_')
            code = fields[0]
            lab_num = int(fields[2])
        except:
            return ""

        try:
            fields = fields[1].split('-')
            platform = fields[0]
            strategy = fields[1]
        except:
            return ""

        platform_num = [ v for k, v in platform_dict.iteritems() if platform.lower() == k ]
        if len(platform_num) == 0:
            print '%s: Platform %s not found in platform dict.' % (identifier, platform)
            return ""

        strategy_char = [ v for k, v in strategy_dict.iteritems() if strategy.lower() == k ]
        if len(strategy_char) == 0:
            print '%s: Strategy %s not found in strategy dict.' % (identifier, strategy)
            return "" 
        return "%s_%2.2d%1s%1d" % (code, lab_num, strategy_char[0].upper(), int(platform_num[0]))

    #----------------------------------->>>
    if to.upper() == "OP":

        # hun_35C3 -> hun_tower-insitu_35
        try:
            fields = name.split('_')
  
            code = fields[0]
            lab_num = int(fields[1][:2])
            strategy_char = fields[1][2]
            platform_num = int(fields[1][3])
        except:
            return ""

        platform = [ k for k, v in platform_dict.iteritems() if v == platform_num ]
        if len(platform) == 0:
            print '%s: Platform number %s not found in platform dict.' % (identifier, platform_num)
            return ""
        
        pattern = re.compile(strategy_char, re.IGNORECASE)
        strategy = [ k for k, v in strategy_dict.iteritems() if pattern.search(v) ]
        if len(strategy) == 0:
            print '%s: Strategy character %s not found in strategy list.' % (identifier, strategy_char)
            return ""
        return "%s_%s-%s_%d" % (code, platform[0], strategy[0], lab_num)

def _to_ordinalf(dt):
    """
    Convert :mod:`datetime` to the Gregorian date as UTC float days,
    preserving hours, minutes, seconds and microseconds.  Return value
    is a :func:`float`.
    """

    if hasattr(dt, 'tzinfo') and dt.tzinfo is not None:
        delta = dt.tzinfo.utcoffset(dt)
        if delta is not None:
            dt -= delta

    base =  float(dt.toordinal())
    if hasattr(dt, 'hour'):
        base += (dt.hour/HOURS_PER_DAY + dt.minute/MINUTES_PER_DAY +
                 dt.second/SECONDS_PER_DAY + dt.microsecond/MUSECONDS_PER_DAY
                 )
    return base

def _from_ordinalf(x, tz=None):
    """
    Convert Gregorian float of the date, preserving hours, minutes,
    seconds and microseconds.  Return value is a :class:`datetime`.
    """
    if tz is None: tz = _get_rc_timezone()
    ix = int(x)
    dt = datetime.datetime.fromordinal(ix)
    remainder = float(x) - ix
    hour, remainder = divmod(24*remainder, 1)
    minute, remainder = divmod(60*remainder, 1)
    second, remainder = divmod(60*remainder, 1)
    microsecond = int(1e6*remainder)
    if microsecond<10: microsecond=0 # compensate for rounding errors
    dt = datetime.datetime(
        dt.year, dt.month, dt.day, int(hour), int(minute), int(second),
        microsecond, tzinfo=UTC).astimezone(tz)

    if microsecond>999990:  # compensate for rounding errors
        dt += datetime.timedelta(microseconds=1e6-microsecond)

    return dt

def date2num(d):
    """
    *d* is either a :class:`datetime` instance or a sequence of datetimes.

    Return value is a floating point number (or sequence of floats)
    which gives the number of days (fraction part represents hours,
    minutes, seconds) since 0001-01-01 00:00:00 UTC, *plus* *one*.
    The addition of one here is a historical artifact.  Also, note
    that the Gregorian calendar is assumed; this is not universal
    practice.  For details, see the module docstring.
    """
    try: 
	return np.asarray([_to_ordinalf(val) for val in d])
    except:
        return _to_ordinalf(d)


def num2date(x, tz=None):
    """
    *x* is a float value which gives the number of days
    (fraction part represents hours, minutes, seconds) since
    0001-01-01 00:00:00 UTC *plus* *one*.
    The addition of one here is a historical artifact.  Also, note
    that the Gregorian calendar is assumed; this is not universal
    practice.  For details, see the module docstring.

    Return value is a :class:`datetime` instance in timezone *tz* (default to
    rcparams TZ value).

    If *x* is a sequence, a sequence of :class:`datetime` objects will
    be returned.
    """
    if tz is None: tz = _get_rc_timezone()
    try: 
        return [_from_ordinalf(val, tz) for val in x]
    except:
        return _from_ordinalf(x, tz)


if __name__ == "__main__":
    pass
