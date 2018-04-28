#!/usr/bin/env python
# time_avg_fluxes.py

"""
Author : peters 

Revision History:
File created on 20 Dec 2012.

"""
import sys
sys.path.append('../../')
import os
import sys
import shutil
from dateutil.relativedelta import relativedelta
import datetime
import subprocess

def time_avg(dacycle,avg='transcom'):
    """ Function to create a set of averaged files in a folder, needed to make longer term means """
    
    if avg not in ['transcom','transcom_extended','olson','olson_extended','country','flux1x1']:
        raise IOError,'Choice of averaging invalid'

    analysisdir = dacycle['dir.analysis']

    if not os.path.exists(analysisdir):
        raise IOError,'analysis dir requested (%s) does not exist, exiting...'%analysisdir

    daily_avg(dacycle,avg)

    monthly_avg(dacycle,avg)

    yearly_avg(dacycle,avg)

    longterm_avg(dacycle,avg)

def new_month(dacycle):
    """ check whether we just entered a new month"""

    this_month = dacycle['time.start'].month
    prev_month = (dacycle['time.start']-dacycle['cyclelength']).month

    return (this_month != prev_month)

def new_year(dacycle):
    """ check whether we just entered a new year"""

    this_year = dacycle['time.start'].year
    prev_year = (dacycle['time.start']-dacycle['cyclelength']).year

    return (this_year != prev_year)

def daily_avg(dacycle,avg):
    """ Function to create a set of daily files in a folder, needed to make longer term means """
    
    if avg not in ['transcom','transcom_extended','olson','olson_extended','country','flux1x1']:
        raise IOError,'Choice of averaging invalid'

    analysisdir = dacycle['dir.analysis']
    weekdir = os.path.join(analysisdir , 'data_%s_weekly'%avg)
    daydir = os.path.join(analysisdir , 'data_%s_daily'%avg)

    if not os.path.exists(daydir):
        print "Creating new output directory " + daydir
        os.makedirs(daydir)

    files  = os.listdir(weekdir)
    files = [f for f in files if '-' in f and f.endswith('.nc')]

    fileinfo = {}
    for filename in files:
        date=datetime.datetime.strptime(filename.split('.')[-2],'%Y-%m-%d')
        fileinfo[filename] = date
    
    dt = dacycle['cyclelength']

    for k,v in fileinfo.iteritems():
        cycle_file = os.path.join(weekdir,k)
        for i in range(abs(dt.days)):
            daily_file = os.path.join(daydir,'%s_fluxes.%s.nc'%(avg,(v+datetime.timedelta(days=i)).strftime('%Y-%m-%d')))
            if not os.path.lexists(daily_file):
                os.symlink(cycle_file,daily_file)
                #print daily_file,cycle_file

def monthly_avg(dacycle,avg):
    """ Function to average a set of files in a folder from daily to monthly means """
    
    if avg not in ['transcom','transcom_extended','olson','olson_extended','country','flux1x1']:
        raise IOError,'Choice of averaging invalid'

    analysisdir = dacycle['dir.analysis']

    daydir = os.path.join(analysisdir , 'data_%s_daily'%avg)
    monthdir = os.path.join(analysisdir,'data_%s_monthly'%avg)

    if not os.path.exists(monthdir):
        print "Creating new output directory " + monthdir
        os.makedirs(monthdir)


    files  = os.listdir(daydir)  # get daily files
    files = [f for f in files if '-' in f and f.endswith('.nc')]

    if len(files) < 28:
        print 'No month is yet complete, skipping monthly average'
        return

    fileinfo = {}
    for filename in files:  # parse date from each of them
        date=datetime.datetime.strptime(filename.split('.')[-2],'%Y-%m-%d')
        fileinfo[filename] = date

    years = [d.year for d in fileinfo.values()]   # get actual years
    months = set([d.month for d in fileinfo.values()])  # get actual months
   
    sd = datetime.datetime(min(years),1,1)
    ed = datetime.datetime(max(years)+1,1,1)

    while sd < ed: 

        nd = sd + relativedelta(months=+1)

        ndays_in_month = (nd-sd).days
        
        avg_files = [os.path.join(daydir,k) for k,v in fileinfo.iteritems() if v < nd and v >= sd]
       
        if len(avg_files) != ndays_in_month: # only once month complete 
            #print 'New month (%02d) is not yet complete, skipping monthly average'%(sd.month)
            pass
        else:
            targetfile = os.path.join(monthdir,'%s_fluxes.%s.nc'%(avg,sd.strftime('%Y-%m')))
            if not os.path.exists(targetfile):
                print "New month (%02d) is complete, I have %d days for the next file"%(sd.month,ndays_in_month)
                command = ['ncra','-O']+ avg_files + [targetfile]
                status = subprocess.check_call(command)
            else:
                pass

        sd = nd

def yearly_avg(dacycle,avg):
    """ Function to average a set of files in a folder from monthly to yearly means """

    if avg not in ['transcom','transcom_extended','olson','olson_extended','country','flux1x1']:
        raise IOError,'Choice of averaging invalid'

    analysisdir = dacycle['dir.analysis']
    monthdir = os.path.join(analysisdir , 'data_%s_monthly'%avg )
    yeardir = os.path.join(analysisdir,'data_%s_yearly'%avg)

    if not os.path.exists(yeardir):
        print "Creating new output directory " + yeardir
        os.makedirs(yeardir)

    files  = os.listdir(monthdir)  # get monthly files
    files = [f for f in files if '-' in f and f.endswith('.nc')]

    if not files:
        print "No full year finished yet, skipping yearly average..."
        return

    fileinfo = {}
    for filename in files:
        date=datetime.datetime.strptime(filename.split('.')[-2],'%Y-%m')
        fileinfo[filename] = date

    years = set([d.year for d in fileinfo.values()])

    sd = datetime.datetime(min(years),1,1)
    ed = datetime.datetime(max(years)+1,1,1)

    while sd < ed: 

        nd = sd + relativedelta(years=+1)
        
        avg_files = [os.path.join(monthdir,k) for k,v in fileinfo.iteritems() if v < nd and v >= sd]
       
        if not len(avg_files) == 12 : 
            print "Year %04d not finished yet, skipping yearly average..."%sd.year
        else:
            targetfile = os.path.join(yeardir,'%s_fluxes.%s.nc'%(avg,sd.strftime('%Y')))
        
            if not os.path.exists(targetfile):
                print "Year %04d is complete, I have 12 months for the next file"%sd.year
                command = ['ncra','-O']+ avg_files + [targetfile]
                status = subprocess.check_call(command)

        sd = nd

def longterm_avg(dacycle,avg):
    """ Function to average a set of files in a folder from monthly to yearly means """

    if avg not in ['transcom','transcom_extended','olson','olson_extended','country','flux1x1']:
        raise IOError,'Choice of averaging invalid'

    analysisdir = dacycle['dir.analysis']

    yeardir = os.path.join(analysisdir , 'data_%s_yearly'%avg )
    longtermdir = os.path.join(analysisdir,'data_%s_longterm'%avg)

    if not os.path.exists(longtermdir):
        print "Creating new output directory " + longtermdir
        os.makedirs(longtermdir)

    files  = os.listdir(yeardir)
    files = [f for f in files if '-' in f and f.endswith('.nc')]

    if not files:
        print "No full year finished yet, skipping longterm average..."
        return

    dates = []
    for filename in files:
        date=datetime.datetime.strptime(filename.split('.')[-2],'%Y')
        dates.append( date )

    avg_files = [os.path.join(yeardir,k) for k in files]
   
    if len(avg_files) > 0 : 
        command = ['ncra','-O']+ avg_files + [os.path.join(longtermdir,'%s_fluxes.%04d-%04d.nc'%(avg,min(dates).year, max(dates).year))]
        status = subprocess.check_call(command)

if __name__ == "__main__":

    from da.tools.initexit import CycleControl

    sys.path.append('../../')

    dacycle = CycleControl(args={'rc':'../../ctdas-ei-nobcb-zoom-ecoregions.rc'})
    dacycle.setup()
    dacycle.parse_times()

    while dacycle['time.end'] < dacycle['time.finish']:
        time_avg(dacycle,avg='flux1x1')
        time_avg(dacycle,avg='transcom')
        time_avg(dacycle,avg='transcom_extended')
        time_avg(dacycle,avg='olson')
        time_avg(dacycle,avg='olson_extended')
        time_avg(dacycle,avg='country')
        dacycle.advance_cycle_times()

