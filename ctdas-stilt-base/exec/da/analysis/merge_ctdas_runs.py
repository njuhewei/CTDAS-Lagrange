#!/usr/bin/env python
# merge_ctdas_runs.py

"""
Author : peters 

Revision History:
File created on 14 Jul 2014.

This scrip merges the analysis directory from multiple projects into one new folder.
It steps over existing analysis output files from weekly means, and then averages these to daily/monthy/yearly values.

"""

import datetime as dt
import os
import sys
import shutil
import time_avg_fluxes as tma

basedir = '/Storage/CO2/ingrid/'
basedir2 = '/Storage/CO2/peters/'
targetproject = 'geocarbon-ei-sibcasa-gfed4-zoom-gridded-combined-convec-20011230-20130101'
targetdir = os.path.join(basedir2,targetproject)

sources = {
            '2000-01-01 through 2011-12-31': os.path.join(basedir,'carbontracker','geocarbon-ei-sibcasa-gfed4-zoom-gridded-convec-combined'),
            '2012-01-01 through 2012-12-31': os.path.join(basedir2,'geocarbon-ei-sibcasa-gfed4-zoom-gridded-convec-20111231-20140101'),
            }

dirs = ['flux1x1','transcom','country','olson']

dacycle = {}
dacycle['time.start'] = dt.datetime(2000,12,30)
dacycle['time.end'] = dt.datetime(2013,1,1)
dacycle['cyclelength'] = dt.timedelta(days=7)
dacycle['dir.analysis'] = os.path.join(targetdir,'analysis')


if __name__ == "__main__":

    if not os.path.exists(targetdir):
        os.makedirs(targetdir)
    if not os.path.exists(os.path.join(targetdir,'analysis')):
        os.makedirs(os.path.join(targetdir,'analysis') )
    for nam in dirs:
        if not os.path.exists(os.path.join(targetdir,'analysis','data_%s_weekly'%nam)):
            os.makedirs(os.path.join(targetdir,'analysis','data_%s_weekly'%nam) )

    timedirs=[]
    for ss,vv in sources.iteritems():
        sds,eds = ss.split(' through ')
        sd = dt.datetime.strptime(sds,'%Y-%m-%d')
        ed = dt.datetime.strptime(eds,'%Y-%m-%d')
        timedirs.append([sd,ed,vv])
        print sd,ed, vv

    while dacycle['time.start'] < dacycle['time.end']:

        # copy the weekly flux1x1 file from the original dir to the new project dir

        for td in timedirs:
            if dacycle['time.start'] >= td[0] and dacycle['time.start'] <= td[1]:
                indir=td[2]

        # Now time avg new fluxes

        infile = os.path.join(indir,'analysis','data_flux1x1_weekly','flux_1x1.%s.nc'%(dacycle['time.start'].strftime('%Y-%m-%d') ) ) 
        #print os.path.exists(infile),infile
        shutil.copy(infile,infile.replace(indir,targetdir) )
        tma.time_avg(dacycle,avg='flux1x1')

        infile = os.path.join(indir,'analysis','data_transcom_weekly','transcom_fluxes.%s.nc'%(dacycle['time.start'].strftime('%Y-%m-%d') ) ) 
        #print os.path.exists(infile),infile
        shutil.copy(infile,infile.replace(indir,targetdir) )
        tma.time_avg(dacycle,avg='transcom')

        infile = os.path.join(indir,'analysis','data_olson_weekly','olson_fluxes.%s.nc'%(dacycle['time.start'].strftime('%Y-%m-%d') ) ) 
        #print os.path.exists(infile),infile
        shutil.copy(infile,infile.replace(indir,targetdir) )
        tma.time_avg(dacycle,avg='olson')

        infile = os.path.join(indir,'analysis','data_country_weekly','country_fluxes.%s.nc'%(dacycle['time.start'].strftime('%Y-%m-%d') ) ) 
        #print os.path.exists(infile),infile
        shutil.copy(infile,infile.replace(indir,targetdir) )
        tma.time_avg(dacycle,avg='country')

        dacycle['time.start'] += dacycle['cyclelength']

    

    
