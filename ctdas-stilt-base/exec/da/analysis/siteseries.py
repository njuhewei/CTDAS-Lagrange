#!/usr/bin/env python
# siteseries_c.py

"""
Author : peters

Revision History:
File created on 23 Dec 2008.

"""
import sys
sys.path.append('../../')

import matplotlib
matplotlib.use('pdf')

import sys
import os
from da.tools.general import create_dirs
import matplotlib.dates as pltdt
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf
from matplotlib.font_manager import FontProperties
import numpy as np
import datetime as dt
import da.tools.io4 as io
import logging
import copy
from da.analysis.summarize_obs import nice_lon, nice_lat, nice_alt
##import Image
import urllib2
import StringIO

"""
    General data needed to set up proper aces inside a figure instance
"""

ts_x1 = 0.125
ts_x2 = 0.825
ts_xspan = ts_x2 - ts_x1

ts_y1 = 0.18
ts_y2 = 0.72
ts_yspan = ts_y2 - ts_y1

markersize = 2
fontsize = 16

"""
    Three routines that provide loops around the number of sites. These also create the rundirs and such and make sure the
    figures that were created are stamped and saved properly
"""

def site_timeseries(analysisdir,option='final'):
    """***************************************************************************************
    Call example:

    "***************************************************************************************"""
    #
    # Repeat all options to user
    #
    #
    # Create directories if needed
    #

    mrdir = os.path.join(analysisdir, 'timeseries_molefractions')
    if not os.path.exists(mrdir):
        create_dirs(mrdir)
    #
    # Make a dictionary of available sites and the NetCDF file names associated with them
    #
    sitelist = os.listdir(os.path.join(analysisdir, 'data_molefractions'))
    sitelist = [f for f in sitelist if f.endswith('.nc')]
    #
    # Loop over site and sitefiles
    #
    for sitefile in sitelist:
        #
        # Create filename and extract site codes for the sites that had day/night data separated
        #
        filename = os.path.join(analysisdir, 'data_molefractions', sitefile)
        saveas = os.path.join(mrdir, sitefile[:-3] + '_timeseries')

        if not os.path.exists(saveas+'.pdf'):
            logging.debug('Making timeseries figures for %s: ' % sitefile)
            #
            # Create a figure instance to hold plot
            #
            fig = plt.figure(1, figsize=(15, 7,))#,frameon=False)
            #
            # Make plot
            #
            fig = timevssite_new(fig, filename)
            #
            # Save image
            #
            fig.savefig(saveas+'.pdf', dpi=100)
            fig.savefig(saveas+'.png', dpi=50)
            fig.savefig(saveas+'.large.png', dpi=150)
            #
            plt.close(fig)
        #
        # Residuals next
        #
        saveas = os.path.join(mrdir, sitefile[:-3] + '_residuals')
        if not os.path.exists(saveas+'.pdf'):
            logging.debug('Making residuals figures for %s: ' % sitefile)
            #
            # Create a figure instance to hold plot
            #
            fig = plt.figure(1, figsize=(15, 7,))#,frameon=False)
            #
            # Make plot
            #
            fig = residuals_new(fig, filename, option)
            #
            # Save image
            #
            fig.savefig(saveas+'.pdf', dpi=100)
            fig.savefig(saveas+'.png', dpi=50)
            fig.savefig(saveas+'.large.png', dpi=150)
            #
            # next in loop over sites
            #
            plt.close(fig)
        #
        # histograms next
        #
        saveas = os.path.join(mrdir, sitefile[:-3] + '_histograms')
        if not os.path.exists(saveas+'.pdf'):
            logging.debug('Making histograms figures for %s: ' % sitefile)
            #
            # Create a figure instance to hold plot
            #
            fig = plt.figure(1, figsize=(15, 7,))#,frameon=False)
            #
            # Make plot
            #
            fig = timehistograms_new(fig, filename, option)
            #
            # Save image
            #
            fig.savefig(saveas+'.pdf', dpi=100)
            fig.savefig(saveas+'.png', dpi=50)
            fig.savefig(saveas+'.large.png', dpi=150)
            #
            # next in loop over sites
            #
            plt.close(fig)
        #
        # html table next
        #
        """
        saveas = os.path.join(mrdir, sitefile[:-3] + '.html')
        if not os.path.exists(saveas):
            logging.debug('Making html info table for %s' % sitefile)
            f = io.CT_CDF(filename, 'read')
            with open(saveas, "wt") as fout:
                with open("sitetable.html", "rt") as fin:
                    for lineout in fin:
                        lineout = lineout.replace('site_name',f.site_name)
                        lineout = lineout.replace('site_country',f.site_country)
                        lineout = lineout.replace('site_latitude',nice_lat(f.site_latitude,'html'))
                        lineout = lineout.replace('site_longitude',nice_lon(f.site_longitude,'html'))
                        lineout = lineout.replace('site_elevation',str(f.site_elevation))
                        lineout = lineout.replace('intake_height',str(f.variables['intake_height'][:].max()))
                        if 'co2_hip' in sitefile:
                            lineout = lineout.replace('site_map','http://carbontracker.eu/development/data/ct_europe/timeseries_molefractions/co2_hip_map.jpg')
                        else: lineout = lineout.replace('site_map',str(f.site_map))
                        lineout = lineout.replace('lab_abbr',f.lab_abbr)
                        lineout = lineout.replace('lab_1_name',f.lab_name)
                        lineout = lineout.replace('lab_1_country',f.lab_country)
                        if 'lab_provider' in lineout:
                            if 'Morgui' in f.provider_1_name:
                                lineout = lineout.replace('lab_provider','J.A. Morgui / M. Ramonet')
                            else: lineout = lineout.replace('lab_provider',f.provider_1_name)
                        if 'lab_url' in lineout:
                            if f.lab_url <> '':
                                lineout = lineout.replace('lab_url',f.lab_url)
                            else: lineout = ''
                        lineout = lineout.replace('lab_logo',f.lab_logo)
                        lineout = lineout.replace('dataset_selection',f.dataset_selection_tag)
                        lineout = lineout.replace('dataset_project',f.dataset_project)
                        lineout = lineout.replace('dataset_calibration_scale',f.dataset_calibration_scale)
                        if f.variables['modeldatamismatch'][:].max() > 0.0009:
                            lineout = lineout.replace('assimilated','No')
                        else: lineout = lineout.replace('assimilated','Yes')
                        fout.write(lineout)
            f.close()
        """
def timehistograms_new(fig, infile, option='final'):
    """
    This routine makes two side-by-side histograms representing summer and winter PDFs of the residuals. It uses the special
    x-axis and y-axis definitions from above. Note that currently, the PDFs are based on forecast-observed CO2, and not on
    optimized-observed CO2.
    """

    fontsize = 17
    #
    # Get data
    #
    f = io.CT_CDF(infile, 'read')
    species = f.dataset_parameter
    if species == 'co2':
        molefac=1e6
        units = '$\mu$mol mol$^{-1}$'
        species = "CO$_2$"
    if species == 'co2c13':
        molefac=1.0
        units = 'permil'
        species = "$\delta^{13}$C"

    date = f.get_variable('time')
    obs = f.get_variable('value') * molefac
    mdm = f.get_variable('modeldatamismatch') * molefac
    hphtr = f.get_variable('totalmolefractionvariance_forecast') * molefac * molefac
    if option == 'final':
        simulated = f.get_variable('modelsamplesmean') * molefac
    if option == 'forecast':
        simulated = f.get_variable('modelsamplesmean_forecast') * molefac
    flags = f.get_variable('flag_forecast')

    longsitestring = f.site_name + ', ' + f.site_country
    location = nice_lat(f.site_latitude,'python') + ', ' + nice_lon(f.site_longitude,'python') + ', ' + nice_alt(f.site_elevation)

    SDSInfo = {}
    for k in f.ncattrs():
        SDSInfo[k] = f.getncattr(k)

    f.close()

    pydates = np.array([dt.datetime(1970, 1, 1) + dt.timedelta(seconds=int(d)) for d in date])

    sampled = (np.ma.getmaskarray(simulated) == False)

    if len(sampled.nonzero()[0]) < 2:
        logging.warning("Too few simulated values found, continuing...")
        return fig

    simulated = simulated.compress(sampled)
    obs = obs.compress(sampled)
    pydates = pydates.compress(sampled)
    mdm = mdm.compress(sampled)
    hphtr = hphtr.compress(sampled)
    flags = flags.compress(sampled)
    #mdm=ma.masked_invalid(mdm)

    residual = simulated - obs
    if option == 'final':
        chisquared = (residual ** 2) / mdm
    elif option == 'forecast':
        chisquared = (residual ** 2) / hphtr

    rejected = (flags == 2.0)
    notused = (flags == 99.0)

    #if notused.all():
    #    return fig
    #else:
    obslabel = 'Residual'

    sd = pydates[0]
    ed = pydates[-1]

    summer = [i for i, d in enumerate(pydates) if d.month in [6, 7, 8, 9] ]   # JJAS
    winter = [i for i, d in enumerate(pydates) if d.month in [11, 12, 1, 2, 3, 4] ]  # NDJFMA

    # Create two side-by-side axes, turn off their frame

    ax1 = fig.add_axes([0.05, 0.18, 0.4, 0.7])
    ax2 = fig.add_axes([0.55, 0.18, 0.4, 0.7])


    # Loop simultaneously over ax1/ax2 and summer/winter values

    for ax, sel in zip([ax1, ax2], [summer, winter]):

        if not np.array(sel).any(): continue

        # Subselect data for winter/summer

        sel_obs = obs.take(sel)
        sel_fc = simulated.take(sel)
        sel_hqhr = hphtr.take(sel)
        sel_mdm = mdm.take(sel)
        sel_flags = flags.take(sel)
        sel_rej = rejected.take(sel)


        # Calculate residual and chi squared values

        #res = sel_fc - sel_obs
        if option == 'final':
            res = sel_fc.compress(sel_flags != 2) - sel_obs.compress(sel_flags != 2)
            chi = res / np.sqrt(sel_mdm.compress(sel_flags != 2))
        elif option == 'forecast':
            res = sel_fc - sel_obs
            chi = res / np.sqrt(sel_hqhr)
        #res=ma.masked_invalid(res)

        # Get a scaling factor for the x-axis range. Now we will include 5 standard deviations

        sc = res.std()

        # If there is too little data for a reasonable PDF, skip to the next value in the loop

        if res.shape[0] < 10: continue

        # make a histogram plot of the residuals with a minimum of 10 bins, and maximum of N/10 bins, normalize the PDF to an area of 1.0

        n, bins, patches = ax.hist(res, max(res.shape[0] / 10, 10), normed=1)

        #print res.mean(), res.sum(),n

        # Change the colors on the bars

        p = plt.setp(patches, 'facecolor', 'tan', 'edgecolor', 'tan', label='None')

        # Create two normal distributions for the line plots over the interval of the x-axis

        bins = np.arange(-5 * sc, 5 * sc, 0.1)
        n = normpdf(bins, res.mean(), res.std())
        l = ax.plot(bins, n, 'b-', linewidth=2) # plot the PDF of the histogram in blue
        n = normpdf(bins, 0.0 , sel_mdm.mean())
        l = ax.plot(bins, n, 'g-', linewidth=2) # plot the PDF of the model-data-mismatch in green
        #
        # Add a legend, not as a legend object but simply as text labels
        #
        if option == 'final':
            strX = ''
        elif option == 'forecast':
            strX = 'Inn. '
        if chi.mean() != chi.mean() or mdm.mean() < 900:
            labs = [
                '%.2f $\pm$ %.2f' % (res.mean(), res.std()) , \
                'N = %d' % sel_obs.shape[0], \
                '%s$\chi^2$= %.2f'%(strX, (chi**2).mean())
                ]
        else:
            labs = [
                '%.2f $\pm$ %.2f' % (res.mean(), res.std()) , \
                'N = %d' % sel_obs.shape[0]
                ]

        # print the above labels onto the figure. Note that I use relative coordinates for their position by specifying the transform=ax.transAxes

        for i, l in enumerate(labs):
            ax.text(0.75, 0.9 - 0.07 * i, l, transform=ax.transAxes, fontsize=fontsize, horizontalalignment='center', color='blue')
        #
        # Set Tick Font Size on x and y labels
        #
        #dummy = [lab.set_fontsize(20) for lab in ax.get_xticklabels()]
        #dummy = [lab.set_fontsize(20) for lab in ax.get_yticklabels()]

        # set limits on x-axis and get limits on y-axis to determine the position of the x-axis labels (offset keyword to make_yaxis)

        ax.set_xlim(-5 * sc, 5 * sc)

        ax.spines['left'].set_position(('data', 0))
        ax.spines['right'].set_color('none')
        ax.spines['right'].axis.set_ticks([])
        ax.spines['bottom'].set_position(('data', 0))
        ax.spines['top'].set_color('none')
        #ax.spines['left'].set_smart_bounds(True)
        #ax.spines['bottom'].set_smart_bounds(True)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['bottom'].set_position(('outward', 10))

        matplotlib.rcParams.update({'font.size': 18})
        ax.xaxis.set_ticks_position('bottom')

        ax.set_xlabel('[%s]'%units,size=16)
    #
    # All custom titles and auxiliary info are placed onto the figure directly (fig.text) in relative coordinates
    #
    fig.text(0.5, 0.02, 'Simulated - Observed %s [%s]\nData from %s to %s' %(species,units,pydates[0].strftime('%d-%b-%Y'), pydates[-1].strftime('%d-%b-%Y')), horizontalalignment='center', fontsize=fontsize)
    #fig.text(0.75,0.02,'Simulated - Observed\n CO$_2$ ($\mu$mol/mol)',horizontalalignment='center',fontsize=fontsize)
    fig.text(0.5, 0.35, 'model-data\nmismatch:\n%.2f %s' % (sel_mdm.mean(), units), horizontalalignment='center', fontsize=fontsize, color='green')
    fig.text(0.12, 0.75, 'NH Summer\n(Jun-Sep)', horizontalalignment='center', fontsize=fontsize)
    fig.text(0.62, 0.75, 'NH Winter\n(Nov-Apr)', horizontalalignment='center', fontsize=fontsize)
    #
    # Title
    #

    plt.suptitle('%s [%s]\n%s, %s, %s ' % (longsitestring, location , SDSInfo['dataset_project'], SDSInfo['lab_1_name'], SDSInfo['lab_1_country'],), fontsize=fontsize + 4)

    #
    # Add info to plot
    #
    font0= FontProperties(size=15,style='italic',weight='bold')
    txt='' #'CTDAS-WRF-STILT\n $\copyright$ University of Groningen'
    clr='red'
    #fig.text(0.82,0.01,txt,ha='left',font_properties = font0, color=clr )

    #now = dt.datetime.today()
    #str1 = 'CTDAS2012\n' + now.strftime('%d/%m/%y')
    #fig.text(0.93, 0.95, str1, fontsize=0.75 * fontsize, color='0.5')
    #str1 = 'data provided by %s'%SDSInfo['provider_1_name']
    #fig.text(0.12,0.16,str1,fontsize=0.8*fontsize,color='0.75')

    try:
        img = urllib2.urlopen(SDSInfo['lab_logo']).read()
    except:
        logging.warning("No logo found for this program, continuing...")
        return fig

    ##im = Image.open(StringIO.StringIO(img))
    ##height = im.size[1]
    ##width = im.size[0]

    # We need a float array between 0-1, rather than
    # a uint8 array between 0-255
    ##im = np.array(im).astype(np.float)[::-1, :] / 255

    # With newer (1.0) versions of matplotlib, you can
    # use the "zorder" kwarg to make the image overlay
    # the plot, rather than hide behind it... (e.g. zorder=10)
    ax3 = fig.add_axes([0.425, 0.65, 0.15, 0.15 * height / width])
    ax3.axis('off')
    ##ax3.imshow(im, interpolation='None')

    return fig

def timevssite_new(fig, infile):
    fontsize = 17
    #
    # Get data
    #
    f = io.CT_CDF(infile, 'read')
    species = f.dataset_parameter
    if species == 'co2':
        molefac=1e6
        units = '$\mu$mol mol$^{-1}$'
        species = "CO$_2$"
    if species == 'co2c13':
        molefac=1.0
        units = 'permil'
        species = "$\delta^{13}$C"
    date = f.get_variable('time')
    obs = f.get_variable('value') * molefac
    mdm = f.get_variable('modeldatamismatch') * molefac
    simulated = f.get_variable('modelsamplesmean') * molefac
    flags = f.get_variable('flag_forecast')

    longsitestring = f.site_name + ', ' + f.site_country
    location = nice_lat(f.site_latitude,'python') + ', ' + nice_lon(f.site_longitude,'python') + ', ' + nice_alt(f.site_elevation)

    SDSInfo = {}
    for k in f.ncattrs():
        SDSInfo[k] = f.getncattr(k)

    f.close()

    pydates = np.array([dt.datetime(1970, 1, 1) + dt.timedelta(seconds=int(d)) for d in date])
    select = [i for i,d in enumerate(pydates) if d.year == 2010]

    sampled = (np.ma.getmaskarray(simulated) == False)

    if len(sampled.nonzero()[0]) < 2:
        logging.warning("Too few simulated values found, continuing...")
        return fig

    simulated = simulated[select].compress(sampled[select])
    obs = obs[select].compress(sampled[select])
    pydates = pydates[select].compress(sampled[select])
    mdm = mdm[select].compress(sampled[select])
    flags = flags[select].compress(sampled[select])
    #mdm=ma.masked_invalid(mdm[select])

    residual = simulated - obs

    #print 'sampled',sampled

    rejected = (flags == 2.0)
    notused = (flags == 99.0)

    if notused.any():
        obslabel = 'Observed (not assimilated)'
    else:
        obslabel = 'Observed (assimilated)'

    sd = pydates[0]
    ed = pydates[-1]

    ax1 = fig.add_axes([0.1, 0.12, 0.7, 0.75])
    ax2 = fig.add_axes([0.85, 0.12, 0.12, 0.75])

    ax1.spines['right'].set_color('none')
    ax1.spines['top'].set_color('none')
    ax1.spines['left'].set_linewidth(1.5)
    ax1.spines['bottom'].set_linewidth(1.5)
    ax1.spines['left'].set_position(('outward', 10))
    ax1.spines['bottom'].set_position(('outward', 10))

    ax2.spines['right'].set_color('none')
    ax2.spines['top'].set_color('none')
    ax2.spines['left'].set_linewidth(1.5)
    ax2.spines['bottom'].set_linewidth(1.5)
    ax2.spines['left'].set_position(('outward', 10))
    ax2.spines['bottom'].set_position(('outward', 10))

    markersize = 8
    fontsize = 16

    #print 'flags rejected',flags.compress(rejected),date.compress(rejected)


    p = ax1.plot(pydates, obs, marker='o', markeredgewidth=1, linestyle='None', markerfacecolor='None', \
           markeredgecolor='k', label=obslabel , markersize=markersize)
    #
    # Add the simulated values
    #
    q = ax1.plot(pydates, simulated, marker='o', markeredgewidth=1, linestyle='None', markerfacecolor='None', \
           markeredgecolor='lightblue', label='Simulated', markersize=markersize)
    #
    # Add the rejected values if available
    #
    if rejected.any():
        r = ax1.plot(pydates.compress(rejected), simulated.compress(rejected), marker='s', markeredgewidth=1, markeredgecolor='r', \
                                  linestyle='None', label='Model Rejected (N=%d)' % len(pydates.compress(rejected)), markersize=markersize)

    #
    # Set up x axis labels
    #
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax1.get_xticklabels()]
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax1.get_yticklabels()]
    #
    # Location and format of xticks
    #
    ax1.xaxis.set_major_locator(pltdt.MonthLocator([2,4,6,8,10,12]))#[7],bymonthday=7))
    ax1.xaxis.set_major_formatter(pltdt.DateFormatter('%Y-%b'))
    #
    # Legend
    #
    leg = ax1.legend(prop=FontProperties(size=(0.75 * fontsize)), borderpad=0.1, loc='upper left')
    #leg.get_frame().set_visible(False)
    leg.set_zorder(20)
    leg.get_frame().set_color('1.0')
    dummy = [lab.set_fontsize(16) for lab in leg.get_texts()]
    #
    # include grid
    #
    ax1.grid(True, ls='-', color='0.75', axis='y')
    ax1.autoscale(enable=True, axis='y', tight=False)
    #ax1.set_ylim(obs.min()-3*residual.std(),obs.max()+5*residual.std())
    #ax1.set_xlim(pltdt.date2num(dt.datetime(sd.year, 1, 1)), pltdt.date2num(dt.datetime(ed.year + 1, 1, 1)))
    ax1.set_xlim(pltdt.date2num(dt.datetime(sd.year, 1, 1)), pltdt.date2num(dt.datetime(ed.year + 1, 1, 1)))
    #ax1.set_ylim(360,430) #LUT

    ym = ax1.get_ylim()
    ymin=ym[0] ; ymax =ym[1]
    for yr in range(sd.year,ed.year+1,2):
        x1=dt.datetime(yr,1,1)
        x2=dt.datetime(yr+1,1,1)
        ax1.fill([x1,x2,x2,x1],[ymin,ymin,ymax,ymax],color='0.9',zorder=1)

    ax1.set_ylim(ymin,ymax)
    #
    #
    # Set Tick Font Size
    #
    #matplotlib.rcParams.update({'font.size': 30})
    ax1.xaxis.set_ticks_position('bottom')
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax1.get_xticklabels()]
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax1.get_yticklabels()]

    #xtitle='Time'
    #ax1.set_xlabel(xtitle, fontsize=fontsize)   # label x axis
    ax1.set_ylabel(r"%s [%s]"% (species,units), fontsize=fontsize + 5)  # label y-axis

    #
    # Axes 2
    #
    residual = residual.compress(flags != 2)
    offset = 0.0
    n, bins, patches = ax2.hist(residual, max(residual.shape[0] / 15, 15), normed=1, orientation='horizontal')
    p = plt.setp(patches, 'facecolor', 'tan' , 'edgecolor', 'tan', label='None', alpha=0.25)

    # Create normal distributions for the line plots over the interval of the x-axis
    sc = residual.std()
    bins = np.arange(-4 * sc, 4 * sc, 0.1)
    n = normpdf(bins, residual.mean(), residual.std())
    l = ax2.plot(n, bins, linestyle='-', color='lightblue', linewidth=1) # plot the PDF of the histogram in blue

    dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax2.get_xticklabels()]
    dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax2.get_yticklabels()]
    labs = [
        '%+.2f $\pm$ %.2f\nN=%d' % (residual.mean(), residual.std(), residual.shape[0],)
        ]
    # print the above labels onto the figure. Note that I use relative coordinates for their position by specifying the transform=ax.transAxes

    ax2.text(0.6, 0.01 + offset, labs[0], transform=ax2.transAxes, fontsize=1.1 * fontsize, horizontalalignment='center', color='k')
    offset += -0.05

    ax2.set_ylim(-4 * sc, 4 * sc)

    ax2.spines['left'].set_position(('axes', 0.0))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].axis.set_ticks([])
    ax2.spines['bottom'].set_position(('axes', 0.5))

    ax2.spines['top'].set_color('none')
    ax2.spines['left'].set_smart_bounds(True)
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.spines['left'].set_linewidth(1.5)
    ax2.spines['bottom'].set_linewidth(1.5)
    ax2.spines['bottom'].set_position(('outward', 10))

    matplotlib.rcParams.update({'font.size': 18})
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticklabels([])

    #ax2.set_ylabel(r"CO$_2$ [ppm]", fontsize=fontsize)  # label y-axis
    #ax2.set_xlabel("frequency", fontsize=fontsize)  # label x-axis
    #ax2.grid(True, axis='y')
    ax2.grid(True, ls='-', color='0.75', axis='y')

    #
    # Title
    #

    plt.suptitle('%s [%s]\n%s, %s, %s ' % (longsitestring, location , SDSInfo['dataset_project'], SDSInfo['lab_1_name'], SDSInfo['lab_1_country'],), fontsize=fontsize + 5)

    #
    # Add info to plot
    #
    font0= FontProperties(size=15,style='italic',weight='bold')
    txt='' #'CTDAS-WRF-STILT\n $\copyright$ University of Groningen'
    clr='red'
    #fig.text(0.82,0.01,txt,ha='left',font_properties = font0, color=clr )

    #now = dt.datetime.today()
    #str1 = 'CTDAS2012\n' + now.strftime('%d/%m/%y')
    #fig.text(0.93, 0.95, str1, fontsize=0.75 * fontsize, color='0.5')
    #str1 = 'data provided by %s' % SDSInfo['provider_1_name']
    #fig.text(0.12, 0.16, str1, fontsize=0.8 * fontsize, color='0.75')

    try:
        img = urllib2.urlopen(SDSInfo['lab_logo']).read()
    except:
        logging.warning("No logo found for this program, continuing...")
        return fig

    ##im = Image.open(StringIO.StringIO(img))
    ##height = im.size[1]
    ##width = im.size[0]

    # We need a float array between 0-1, rather than
    # a uint8 array between 0-255
    ##im = np.array(im).astype(np.float)[::-1, :] / 255

    # With newer (1.0) versions of matplotlib, you can
    # use the "zorder" kwarg to make the image overlay
    # the plot, rather than hide behind it... (e.g. zorder=10)
    ax3 = fig.add_axes([0.7, 0.16, 0.15, 0.15 * height / width])
    ax3.axis('off')
    ##ax3.imshow(im, interpolation='None')

    return fig

def residuals_new(fig, infile, option):

    fontsize = 17
    #
    # Get data
    #
    f = io.CT_CDF(infile, 'read')
    species = f.dataset_parameter
    if species == 'co2':
        molefac=1e6
        units = '$\mu$mol mol$^{-1}$'
        species = "CO$_2$"
    if species == 'co2c13':
        molefac=1.0
        units = 'permil'
        species = "$\delta^{13}$C"
    date = f.get_variable('time')
    obs = f.get_variable('value') * molefac
    mdm = f.get_variable('modeldatamismatch') * molefac
    if option == 'final':
        simulated = f.get_variable('modelsamplesmean') * molefac
    if option == 'forecast':
        simulated = f.get_variable('modelsamplesmean_forecast') * molefac
    hphtr = f.get_variable('totalmolefractionvariance_forecast') * molefac * molefac
    flags = f.get_variable('flag_forecast')

    longsitestring = f.site_name + ', ' + f.site_country
    location = nice_lat(f.site_latitude,'python') + ', ' + nice_lon(f.site_longitude,'python') + ', ' + nice_alt(f.site_elevation)

    SDSInfo = {}
    for k in f.ncattrs():
        SDSInfo[k] = f.getncattr(k)

    f.close()

    pydates = np.array([dt.datetime(1970, 1, 1) + dt.timedelta(seconds=int(d)) for d in date])

    select = [i for i,d in enumerate(pydates) if d.year == 2010]

    sampled = (np.ma.getmaskarray(simulated) == False)

    if len(sampled.nonzero()[0]) < 2:
        logging.warning("Too few simulated values found, continuing...")
        return fig

    simulated = simulated[select].compress(sampled[select])
    obs = obs[select].compress(sampled[select])
    pydates = pydates[select].compress(sampled[select])
    mdm = mdm[select].compress(sampled[select])
    hphtr = hphtr[select].compress(sampled[select])
    flags = flags[select].compress(sampled[select])
    #mdm=ma.masked_invalid(mdm)

    rejected = (flags == 2.0)
    notused = (flags == 99.0)

    residual = simulated - obs
    #if notused.all():
    #    return fig
    #else:
    obslabel = 'Residual'

    sd = pydates[0]
    ed = pydates[-1]

    ax1 = fig.add_axes([0.1, 0.12, 0.7, 0.75])
    ax2 = fig.add_axes([0.85, 0.12, 0.12, 0.75])

    ax1.spines['right'].set_color('none')
    ax1.spines['top'].set_color('none')
    ax1.spines['left'].set_linewidth(1.5)
    ax1.spines['bottom'].set_linewidth(1.5)
    ax1.spines['left'].set_position(('outward', 10))
    ax1.spines['bottom'].set_position(('outward', 10))

    ax2.spines['right'].set_color('none')
    ax2.spines['top'].set_color('none')
    ax2.spines['left'].set_linewidth(1.5)
    ax2.spines['bottom'].set_linewidth(1.5)
    ax2.spines['left'].set_position(('outward', 10))
    ax2.spines['bottom'].set_position(('outward', 10))

    markersize = 8
    fontsize = 16

    p = ax1.plot(pydates, residual, marker='o', markeredgewidth=1, linestyle='None', markerfacecolor='None', \
           markeredgecolor='k', label=obslabel , markersize=markersize)
    #
    # Add the model-data mismatch
    #
    q = ax1.fill_between(pydates, mdm, -1.0 * mdm, label='model-data mismatch', color='tan', alpha=0.25, zorder=5)

    #
    # Add the rejected values if available
    #

    #for i in range(len(residual.compress(rejected))):

        #print "rejected",residual.compress(rejected)[i],pydates.compress(rejected)[i]

    if rejected.any():
       r = ax1.plot(pydates.compress(rejected), residual.compress(rejected), marker='s', markeredgewidth=1, markeredgecolor='r', markerfacecolor='red', \
                                  linestyle='None', label='Model Rejected (N=%d)' % len(pydates.compress(rejected)), markersize=markersize)

    #
    # Axes 2
    #
    if option == 'final':
        residual = simulated.compress(flags != 2) - obs.compress(flags != 2)
        pydates = pydates.compress(flags != 2)
        mdm = mdm.compress(flags != 2)
        chisquared = (residual ** 2) / mdm
    elif option == 'forecast':
        chisquared = (residual ** 2) / hphtr
    offset = 0.0

    n, bins, patches = ax2.hist(residual, max(residual.shape[0] / 15, 15), normed=1, orientation='horizontal')
    p = plt.setp(patches, 'facecolor', 'tan' , 'edgecolor', 'tan', label='None', alpha=0.25)

    # Create normal distributions for the line plots over the interval of the x-axis

    sc = residual.std()
    bins = np.arange(-4 * sc, 4 * sc, 0.1)
    n = normpdf(bins, residual.mean(), residual.std())
    l = ax2.plot(n, bins, linestyle='-', color='lightblue', linewidth=1) # plot the PDF of the histogram in blue

    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax2.get_xticklabels()]
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax2.get_yticklabels()]
    if option == 'final':
        strX = ''
    elif option == 'forecast':
        strX = 'Inn. '
    if chisquared.mean() != chisquared.mean() or mdm.mean() < 900:
        labs = [
            '%+.2f $\pm$ %.2f\nN=%d\n%s $\chi^2$ = %5.2f'%(residual.mean(), residual.std(), residual.shape[0], strX, chisquared.mean(),)
            ]
    else:
        labs = [
            '%+.2f $\pm$ %.2f\nN=%d'%(residual.mean(), residual.std(), residual.shape[0],)
            ]

    # print the above labels onto the figure. Note that I use relative coordinates for their position by specifying the transform=ax.transAxes

    ax2.text(0.6, 0.01 + offset, labs[0], transform=ax2.transAxes, fontsize=1.1 * fontsize, horizontalalignment='center', color='k')
    offset += -0.05

    ax2.set_ylim(-4 * sc, 4 * sc)

    ax2.spines['left'].set_position(('axes', 0.0))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].axis.set_ticks([])
    ax2.spines['bottom'].set_position(('axes', 0.5))

    ax2.spines['top'].set_color('none')
    ax2.spines['left'].set_smart_bounds(True)
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.spines['left'].set_linewidth(1.5)
    ax2.spines['bottom'].set_linewidth(1.5)
    ax2.spines['bottom'].set_position(('outward', 10))

    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticklabels([])

    #ax2.set_ylabel(r"CO$_2$ [ppm]", fontsize=fontsize)  # label y-axis
    #ax2.set_xlabel("frequency", fontsize=fontsize)  # label x-axis
    ax2.grid(True, ls='-', color='0.75', axis='y')

    #
    # Set up x axis labels
    #
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax1.get_xticklabels()]
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax1.get_yticklabels()]
    #
    # Location and format of xticks
    #
    ax1.xaxis.set_major_locator(pltdt.MonthLocator([2,4,6,8,10,12]))#[7],bymonthday=7))
    ax1.xaxis.set_major_formatter(pltdt.DateFormatter('%Y-%b'))
    #
    # Legend
    #
    leg = ax1.legend(prop=FontProperties(size=(0.75 * fontsize)), borderpad=0.1, loc='upper left')
    #leg.get_frame().set_visible(False)
    leg.set_zorder(20)
    leg.get_frame().set_color('1.0')
    dummy = [lab.set_fontsize(16) for lab in leg.get_texts()]
    #
    # include grid
    #
    ax1.grid(True, ls='-', color='0.75', axis='y')
    ax1.set_xlim(pltdt.date2num(dt.datetime(sd.year, 1, 1)), pltdt.date2num(dt.datetime(ed.year + 1, 1, 1)))

    ym = ax1.get_ylim()
    ymin=ym[0] ; ymax =ym[1]
    for yr in range(sd.year,ed.year+1,2):
        x1=dt.datetime(yr,1,1)
        x2=dt.datetime(yr+1,1,1)
        ax1.fill([x1,x2,x2,x1],[ymin,ymin,ymax,ymax],color='0.9',zorder=1)

    #ax1.set_ylim(ymin,ymax)
    ax1.set_ylim(-4 * sc, 4 * sc)
    #
    #
    # Set Tick Font Size
    #
    matplotlib.rcParams.update({'font.size': 18})
    ax1.xaxis.set_ticks_position('bottom')
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax1.get_xticklabels()]
    #dummy = [lab.set_fontsize(0.9 * fontsize) for lab in ax1.get_yticklabels()]

    #xtitle='Time'
    #ax1.set_xlabel(xtitle, fontsize=fontsize)   # label x axis
    ax1.set_ylabel(r"%s [%s]"%(species,units), fontsize=fontsize+5)  # label y-axis
    #
    # Title
    #

    plt.suptitle('%s [%s]\n%s, %s, %s ' % (longsitestring, location , SDSInfo['dataset_project'], SDSInfo['lab_1_name'], SDSInfo['lab_1_country'],), fontsize=fontsize + 5)

    #
    # Add info to plot
    #
    font0= FontProperties(size=15,style='italic',weight='bold')
    txt=''  #'CTDAS-WRF-STILT\n $\copyright$ University of Groningen'
    clr='red'
    #fig.text(0.82,0.01,txt,ha='left',font_properties = font0, color=clr )

    #now = dt.datetime.today()
    #str1 = 'CTDAS2012\n' + now.strftime('%d/%m/%y')
    #fig.text(0.93, 0.95, str1, fontsize=0.75 * fontsize, color='0.5')
    #str1 = 'data provided by %s' % SDSInfo['provider_1_name']
    #fig.text(0.12, 0.16, str1, fontsize=0.8 * fontsize, color='0.75')

    try:
        img = urllib2.urlopen(SDSInfo['lab_logo']).read()
    except:
        logging.warning("No logo found for this program, continuing...")
        return fig

    ##im = Image.open(StringIO.StringIO(img))
    ##height = im.size[1]
    ##width = im.size[0]

    # We need a float array between 0-1, rather than
    # a uint8 array between 0-255
    ##im = np.array(im).astype(np.float)[::-1, :] / 255

    # With newer (1.0) versions of matplotlib, you can
    # use the "zorder" kwarg to make the image overlay
    # the plot, rather than hide behind it... (e.g. zorder=10)
    ax3 = fig.add_axes([0.7, 0.16, 0.15, 0.15 * height / width])
    ax3.axis('off')
    ##ax3.imshow(im, interpolation='None')

    return fig


# main body if called as script

if __name__ == '__main__':    # started as script

    sys.path.append('../../')

    logging.root.setLevel(logging.DEBUG)

    analysisdir = "../../../analysis/"
    site_timeseries(analysisdir,option='forecast')

    sys.exit(0)


