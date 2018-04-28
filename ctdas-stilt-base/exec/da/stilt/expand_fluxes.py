#!/usr/bin/env python
# expand_fluxes.py
import sys
sys.path.append('../../')
import os
from datetime import datetime, timedelta

import logging
import numpy as np
from da.tools.general import date2num, num2date
import da.tools.io4 as io
from da.analysis.tools_regions import globarea, state_to_grid
from da.tools.general import create_dirs
from da.analysis.tools_country import countryinfo  # needed here
from da.analysis.tools_transcom import transcommask, ExtendedTCRegions


import da.analysis.tools_transcom as tc
import da.analysis.tools_country as ct
import da.analysis.tools_time as timetools



"""
Author: Wouter Peters (Wouter.Peters@noaa.gov)

Revision History:
File created on 21 Ocotber 2008.

"""

def proceed_dialog(txt, yes=['y', 'yes'], all=['a', 'all', 'yes-to-all']):
    """ function to ask whether to proceed or not """
    response = raw_input(txt)
    if response.lower() in yes:
        return 1
    if response.lower() in all:
        return 2
    return 0

def save_weekly_avg_1x1_data(dacycle, statevector):
    """
        Function creates a NetCDF file with output on 1x1 degree grid. It uses the flux data written by the
        :class:`~da.baseclasses.obsoperator.ObsOperator.py`, and multiplies these with the mapped parameters and
        variance (not covariance!) from the :class:`~da.baseclasses.statevector.StateVector`.

           :param dacycle: a :class:`~da.tools.initexit.CycleControl` object
           :param statevector: a :class:`~da.baseclasses.statevector.StateVector`
           :rtype: None
    """
#
    dirname = create_dirs(os.path.join(dacycle['dir.analysis'], 'data_flux1x1_weekly'))
#
# Some help variables
#
    dectime0 = date2num(datetime(2000, 1, 1))
    dt = dacycle['cyclelength']
    startdate = dacycle['time.start']
    enddate = dacycle['time.end']
    nlag = statevector.nlag

    logging.debug("DA Cycle start date is %s" % startdate.strftime('%Y-%m-%d %H:%M'))
    logging.debug("DA Cycle end   date is %s" % enddate.strftime('%Y-%m-%d %H:%M'))

#
# Create or open NetCDF output file
#
    saveas = os.path.join(dirname, 'flux_1x1.%s.nc' % startdate.strftime('%Y-%m-%d'))
    ncf = io.CT_CDF(saveas, 'write')

#
# Create dimensions and lat/lon grid
#
    dimgrid = ncf.add_latlon_dim()
    dimensemble = ncf.add_dim('members', statevector.nmembers)
    dimdate = ncf.add_date_dim()
#
# set title and tell GMT that we are using "pixel registration"
#
    setattr(ncf, 'Title', 'CarbonTracker fluxes')
    setattr(ncf, 'node_offset', 1)
#
# skip dataset if already in file
#
    ncfdate = date2num(startdate) - dectime0 + dt.days / 2.0
    skip = ncf.has_date(ncfdate)
    if skip:
        logging.warning('Skipping writing of data for date %s : already present in file %s' % (startdate.strftime('%Y-%m-%d'), saveas))
    else:

#
# if not, process this cycle. Start by getting flux input data from CTDAS
#
        filename = os.path.join(dacycle['dir.output'], 'flux1x1_%s_%s.nc' % (startdate.strftime('%Y%m%d%H'), enddate.strftime('%Y%m%d%H')))
        file = io.ct_read(filename, 'read')
        #bio = np.array(file.get_variable(dacycle.dasystem['background.co2.bio.flux']))
        gpp = np.array(file.get_variable(dacycle.dasystem['background.co2.gpp.flux']))
        res = np.array(file.get_variable(dacycle.dasystem['background.co2.res.flux']))
        bio = res - gpp
        ocean = np.array(file.get_variable(dacycle.dasystem['background.co2.ocean.flux']))
        fire = np.array(file.get_variable(dacycle.dasystem['background.co2.fires.flux']))
        fossil = np.array(file.get_variable(dacycle.dasystem['background.co2.fossil.flux']))
        #mapped_parameters   = np.array(file.get_variable(dacycle.dasystem['final.param.mean.1x1']))
        if dacycle.dasystem['background.co2.biosam.flux'] in file.variables.keys():
            sam = True
            biosam = np.array(file.get_variable(dacycle.dasystem['background.co2.biosam.flux']))
            firesam = np.array(file.get_variable(dacycle.dasystem['background.co2.firesam.flux']))
        else: sam = False
        file.close()

        if sam:
            bio = bio + biosam
            fire = fire + firesam

        next = ncf.inq_unlimlen()[0]


# Start adding datasets from here on, both prior and posterior datasets for bio and ocn

        for prior in [True, False]:
#
# Now fill the statevector with the prior values for this time step. Note that the prior value for this time step
# occurred nlag time steps ago, so we make a shift in the output directory, but only if we are more than nlag cycle away from the start date..
#

            if prior:
                qual_short = 'prior'
                for n in range(nlag, 0, -1):
                    priordate = startdate + n*dt - timedelta(dt.days * n)
                    savedir = dacycle['dir.output'].replace(startdate.strftime('%Y%m%d'), priordate.strftime('%Y%m%d'))
                    filename = os.path.join(savedir, 'savestate_%s.nc' % priordate.strftime('%Y%m%d'))
                    if os.path.exists(filename):
                        statevector.read_from_file(filename, qual=qual_short)
                        gridmean, gridensemble = statevector.state_to_grid(lag=n)

# Replace the mean statevector by all ones (assumed priors)

                        #gridmean = statevector.vector2grid(vectordata=np.ones(statevector.nparams,))
                        gridmean = statevector.vector2grid(vectordata=np.zeros(statevector.nparams,))

                        logging.debug('Read prior dataset from file %s, sds %d: ' % (filename, n))
                        break
            else:
                qual_short = 'opt'
                savedir = dacycle['dir.output']
                filename = os.path.join(savedir, 'savestate_%s.nc' % startdate.strftime('%Y%m%d'))
                statevector.read_from_file(filename, qual=qual_short)
                gridmean, gridensemble = statevector.state_to_grid(lag=1)

                logging.debug('Read posterior dataset from file %s, sds %d: ' % (filename, 1))
#
# if prior, do not multiply fluxes with parameters, otherwise do
#
            #print gridensemble.shape, bio.shape, gridmean.shape
            #biomapped = bio * gridmean
            #oceanmapped = ocean * gridmean
            #biovarmapped = bio * gridensemble
            #oceanvarmapped = ocean * gridensemble
            biomapped = bio + gridmean
            oceanmapped = ocean + gridmean
            biovarmapped = gridensemble
            oceanvarmapped = gridensemble

#
#
#  For each dataset, get the standard definitions from the module mysettings, add values, dimensions, and unlimited count, then write
#
            savedict = ncf.standard_var(varname='bio_flux_' + qual_short)
            savedict['values'] = biomapped.tolist()
            savedict['dims'] = dimdate + dimgrid
            savedict['count'] = next
            ncf.add_data(savedict)
#
            savedict = ncf.standard_var(varname='ocn_flux_' + qual_short)
            savedict['values'] = oceanmapped.tolist()
            savedict['dims'] = dimdate + dimgrid
            savedict['count'] = next
            ncf.add_data(savedict)

            #print biovarmapped.shape
            savedict = ncf.standard_var(varname='bio_flux_%s_ensemble' % qual_short)
            savedict['values'] = biovarmapped.tolist()
            savedict['dims'] = dimdate + dimensemble + dimgrid
            savedict['count'] = next
            ncf.add_data(savedict)
#
            savedict = ncf.standard_var(varname='ocn_flux_%s_ensemble' % qual_short)
            savedict['values'] = oceanvarmapped.tolist()
            savedict['dims'] = dimdate + dimensemble + dimgrid
            savedict['count'] = next
            ncf.add_data(savedict)

        # End prior/posterior block

        savedict = ncf.standard_var(varname='fire_flux_imp')
        savedict['values'] = fire.tolist()
        savedict['dims'] = dimdate + dimgrid
        savedict['count'] = next
        ncf.add_data(savedict)
#
        savedict = ncf.standard_var(varname='fossil_flux_imp')
        savedict['values'] = fossil.tolist()
        savedict['dims'] = dimdate + dimgrid
        savedict['count'] = next
        ncf.add_data(savedict)

        area = globarea()
        savedict = ncf.standard_var(varname='cell_area')
        savedict['values'] = area.tolist()
        savedict['dims'] = dimgrid
        ncf.add_data(savedict)
#
        savedict = ncf.standard_var(varname='date')
        savedict['values'] = date2num(startdate) - dectime0 + dt.days / 2.0
        savedict['dims'] = dimdate
        savedict['count'] = next
        ncf.add_data(savedict)

        sys.stdout.write('.')
        sys.stdout.flush()
#
#   Done, close the new NetCDF file
#
    ncf.close()
#
#   Return the full name of the NetCDF file so it can be processed by the next routine
#
    logging.info("Gridded weekly average fluxes now written")

    return saveas

def save_weekly_avg_state_data(dacycle, statevector):
    """
        Function creates a NetCDF file with output for all parameters. It uses the flux data written by the
        :class:`~da.baseclasses.obsoperator.ObsOperator.py`, and multiplies these with the mapped parameters and
        variance (not covariance!) from the :class:`~da.baseclasses.statevector.StateVector`.

           :param dacycle: a :class:`~da.tools.initexit.CycleControl` object
           :param statevector: a :class:`~da.baseclasses.statevector.StateVector`
           :rtype: None
    """

    dirname = create_dirs(os.path.join(dacycle['dir.analysis'], 'data_state_weekly'))
#
# Some help variables
#
    dectime0 = date2num(datetime(2000, 1, 1))
    dt = dacycle['cyclelength']
    startdate = dacycle['time.start']
    enddate = dacycle['time.end']
    nlag = statevector.nlag

    area = globarea()
    vectorarea = statevector.grid2vector(griddata=area, method='sum')

    logging.debug("DA Cycle start date is %s" % startdate.strftime('%Y-%m-%d %H:%M'))
    logging.debug("DA Cycle end   date is %s" % enddate.strftime('%Y-%m-%d %H:%M'))

#
# Create or open NetCDF output file
#
    saveas = os.path.join(dirname, 'statefluxes.nc')
    ncf = io.CT_CDF(saveas, 'write')

#
# Create dimensions and lat/lon grid
#
    dimregs = ncf.add_dim('nparameters', statevector.nparams)
    dimmembers = ncf.add_dim('nmembers', statevector.nmembers)
    dimdate = ncf.add_date_dim()
#
# set title and tell GMT that we are using "pixel registration"
#
    setattr(ncf, 'Title', 'CarbonTracker fluxes')
    setattr(ncf, 'node_offset', 1)
#
# skip dataset if already in file
#
    ncfdate = date2num(startdate) - dectime0 + dt.days / 2.0
    skip = ncf.has_date(ncfdate)
    if skip:
        logging.warning('Skipping writing of data for date %s : already present in file %s' % (startdate.strftime('%Y-%m-%d'), saveas))
    else:
        next = ncf.inq_unlimlen()[0]

#
# if not, process this cycle. Start by getting flux input data from CTDAS
#
        filename = os.path.join(dacycle['dir.output'], 'flux1x1_%s_%s.nc' % (startdate.strftime('%Y%m%d%H'), enddate.strftime('%Y%m%d%H')))

        file = io.ct_read(filename, 'read')
        #bio = np.array(file.get_variable(dacycle.dasystem['background.co2.bio.flux']))
        gpp = np.array(file.get_variable(dacycle.dasystem['background.co2.gpp.flux']))
        res = np.array(file.get_variable(dacycle.dasystem['background.co2.res.flux']))
        bio = res - gpp
        ocean = np.array(file.get_variable(dacycle.dasystem['background.co2.ocean.flux']))
        fire = np.array(file.get_variable(dacycle.dasystem['background.co2.fires.flux']))
        fossil = np.array(file.get_variable(dacycle.dasystem['background.co2.fossil.flux']))
        #mapped_parameters   = np.array(file.get_variable(dacycle.dasystem['final.param.mean.1x1']))
        if dacycle.dasystem['background.co2.biosam.flux'] in file.variables.keys():
            sam = True
            biosam = np.array(file.get_variable(dacycle.dasystem['background.co2.biosam.flux']))
            firesam = np.array(file.get_variable(dacycle.dasystem['background.co2.firesam.flux']))
        else: sam = False
        file.close()

        if sam:
            bio = bio + biosam
            fire = fire + firesam

        next = ncf.inq_unlimlen()[0]

        vectorbio = statevector.grid2vector(griddata=bio * area, method='sum')
        vectorocn = statevector.grid2vector(griddata=ocean * area, method='sum')
        vectorfire = statevector.grid2vector(griddata=fire * area, method='sum')
        vectorfossil = statevector.grid2vector(griddata=fossil * area, method='sum')


# Start adding datasets from here on, both prior and posterior datasets for bio and ocn

        for prior in [True, False]:
#
# Now fill the statevector with the prior values for this time step. Note that the prior value for this time step
# occurred nlag time steps ago, so we make a shift in the output directory, but only if we are more than nlag cycle away from the start date..
#
            if prior:
                qual_short = 'prior'
                for n in range(nlag, 0, -1):
                    priordate = enddate - timedelta(dt.days * n)
                    priordate = startdate + n*dt - timedelta(dt.days * n)
                    savedir = dacycle['dir.output'].replace(startdate.strftime('%Y%m%d'), priordate.strftime('%Y%m%d'))
                    filename = os.path.join(savedir,'savestate_%s.nc' % priordate.strftime('%Y%m%d'))
                    if os.path.exists(filename):
                        statevector.read_from_file(filename, qual=qual_short)
# Replace the mean statevector by all ones (assumed priors)
                        #statemean = np.ones((statevector.nparams,))
                        statemean = np.zeros((statevector.nparams,))

                        choicelag = n
                        logging.debug('Read prior dataset from file %s, lag %d: ' % (filename, choicelag))
                        break
            else:
                qual_short = 'opt'
                savedir = dacycle['dir.output']
                filename = os.path.join(savedir, 'savestate_%s.nc' % startdate.strftime('%Y%m%d'))
                statevector.read_from_file(filename)
                choicelag = 1
                statemean = statevector.ensemble_members[choicelag - 1][0].param_values
                logging.debug('Read posterior dataset from file %s, lag %d: ' % (filename, choicelag))
#
# if prior, do not multiply fluxes with parameters, otherwise do
#
            #data = statemean * vectorbio # units of mole region-1 s-1
            data = statemean + vectorbio
            print statemean.mean(),vectorbio.mean()

            savedict = ncf.standard_var(varname='bio_flux_%s' % qual_short)
            savedict['values'] = data
            savedict['dims'] = dimdate + dimregs
            savedict['count'] = next
            ncf.add_data(savedict)

#
# Here comes a special provision for the posterior flux covariances: these are calculated relative to the prior flux covariance to
# ensure they are indeed smaller due to the data assimilation. If they would be calculated relative to the mean posterior flux, the
# uncertainties would shift just because the mean flux had increased or decreased, which is not what we want.
#
# The implementation is done by multiplying the ensemble with the vectorbio only, and not with the statemean values
# which are assumed 1.0 in the prior always.
#

            members = statevector.ensemble_members[choicelag - 1]
            #deviations = np.array([mem.param_values * vectorbio for mem in members])
            deviations = np.array([mem.param_values[:statevector.nparams] for mem in members])
            deviations = deviations - deviations[0, :]

            savedict = ncf.standard_var(varname='bio_flux_%s_ensemble' % qual_short)

            savedict['values'] = deviations.tolist()
            savedict['dims'] = dimdate + dimmembers + dimregs
            savedict['comment'] = "This is the matrix square root, use (M x M^T)/(nmembers-1) to make covariance"
            savedict['units'] = "mol region-1 s-1"
            savedict['count'] = next
            ncf.add_data(savedict)

            savedict = ncf.standard_var('unknown')
            savedict['name'] = 'bio_flux_%s_std' % qual_short
            savedict['long_name'] = 'Biosphere flux standard deviation, %s' % qual_short
            savedict['values'] = deviations.std(axis=0)
            savedict['dims'] = dimdate + dimregs
            savedict['comment'] = "This is the standard deviation on each parameter"
            savedict['units'] = "mol region-1 s-1"
            savedict['count'] = next
            ncf.add_data(savedict)

            #data = statemean * vectorocn # units of mole region-1 s-1
            data = statemean + vectorocn

            savedict = ncf.standard_var(varname='ocn_flux_%s' % qual_short)
            savedict['values'] = data
            savedict['dims'] = dimdate + dimregs
            savedict['count'] = next
            ncf.add_data(savedict)


#
# Here comes a special provision for the posterior flux covariances: these are calculated relative to the prior flux covariance to
# ensure they are indeed smaller due to the data assimilation. If they would be calculated relative to the mean posterior flux, the
# uncertainties would shift just because the mean flux had increased or decreased, which is not what we want.
#
# The implementation is done by multiplying the ensemble with the vectorocn only, and not with the statemean values
# which are assumed 1.0 in the prior always.
#

            #deviations = np.array([mem.param_values * vectorocn for mem in members])
            deviations = np.array([mem.param_values[:statevector.nparams] for mem in members])
            deviations = deviations - deviations[0, :]

            savedict = ncf.standard_var(varname='ocn_flux_%s_ensemble' % qual_short)
            savedict['values'] = deviations.tolist()
            savedict['dims'] = dimdate + dimmembers + dimregs
            savedict['comment'] = "This is the matrix square root, use (M x M^T)/(nmembers-1) to make covariance"
            savedict['units'] = "mol region-1 s-1"
            savedict['count'] = next
            ncf.add_data(savedict)

            savedict = ncf.standard_var('unknown')
            savedict['name'] = 'ocn_flux_%s_std' % qual_short
            savedict['long_name'] = 'Ocean flux standard deviation, %s' % qual_short
            savedict['values'] = deviations.std(axis=0)
            savedict['dims'] = dimdate + dimregs
            savedict['comment'] = "This is the standard deviation on each parameter"
            savedict['units'] = "mol region-1 s-1"
            savedict['count'] = next
            ncf.add_data(savedict)

        data = vectorfire

        savedict = ncf.standard_var(varname='fire_flux_imp')
        savedict['values'] = data
        savedict['dims'] = dimdate + dimregs
        savedict['count'] = next
        ncf.add_data(savedict)

        data = vectorfossil

        savedict = ncf.standard_var(varname='fossil_flux_imp')
        savedict['values'] = data
        savedict['dims'] = dimdate + dimregs
        savedict['count'] = next
        ncf.add_data(savedict)

        savedict = ncf.standard_var(varname='date')
        savedict['values'] = ncfdate
        savedict['dims'] = dimdate
        savedict['count'] = next
        ncf.add_data(savedict)

        sys.stdout.write('.')
        sys.stdout.flush()
#
#   Done, close the new NetCDF file
#
    ncf.close()
#
#   Return the full name of the NetCDF file so it can be processed by the next routine
#
    logging.info("Vector weekly average fluxes now written")

    return saveas


def save_weekly_avg_tc_data(dacycle, statevector):
    """
        Function creates a NetCDF file with output on TransCom regions. It uses the flux input from the
        function `save_weekly_avg_1x1_data` to create fluxes of length `nparameters`, which are then projected
        onto TC regions using the internal methods from :class:`~da.baseclasses.statevector.StateVector`.

           :param dacycle: a :class:`~da.tools.initexit.CycleControl` object
           :param statevector: a :class:`~da.baseclasses.statevector.StateVector`
           :rtype: None

        This function only read the prior fluxes from the flux_1x1.nc files created before, because we want to convolve
        these with the parameters in the statevector. This creates posterior fluxes, and the posterior covariance for the complete
        statevector in units of mol/box/s which we then turn into TC fluxes and covariances.
    """

#
    dirname = create_dirs(os.path.join(dacycle['dir.analysis'], 'data_tc_weekly'))
#
# Some help variables
#
    dectime0 = date2num(datetime(2000, 1, 1))
    dt = dacycle['cyclelength']
    startdate = dacycle['time.start']
    enddate = dacycle['time.end']
    ncfdate = date2num(startdate) - dectime0 + dt.days / 2.0

    logging.debug("DA Cycle start date is %s" % startdate.strftime('%Y-%m-%d %H:%M'))
    logging.debug("DA Cycle end   date is %s" % enddate.strftime('%Y-%m-%d %H:%M'))

    # Write/Create NetCDF output file
    #
    saveas = os.path.join(dirname, 'tcfluxes.nc')
    ncf = io.CT_CDF(saveas, 'write')
    dimdate = ncf.add_date_dim()
    dimidateformat = ncf.add_date_dim_format()
    dimregs = ncf.add_region_dim(type='tc')
#
# set title and tell GMT that we are using "pixel registration"
#
    setattr(ncf, 'Title', 'CarbonTracker TransCom fluxes')
    setattr(ncf, 'node_offset', 1)
    #

    skip = ncf.has_date(ncfdate)
    if skip:
        logging.warning('Skipping writing of data for date %s : already present in file %s' % (startdate.strftime('%Y-%m-%d'), saveas))
    else:

        # Get input data

        area = globarea()

        infile = os.path.join(dacycle['dir.analysis'], 'data_state_weekly', 'statefluxes.nc')
        if not os.path.exists(infile):
            logging.error("Needed input file (%s) does not exist yet, please create file first, returning..." % infile)
            return None

        ncf_in = io.ct_read(infile, 'read')

        # Transform data one by one

        # Get the date variable, and find index corresponding to the dacycle date

        try:
            dates = ncf_in.variables['date'][:]
        except KeyError:
            logging.error("The variable date cannot be found in the requested input file (%s) " % infile)
            logging.error("Please make sure you create gridded fluxes before making TC fluxes ")
            raise KeyError

        try:
            index = dates.tolist().index(ncfdate)
        except ValueError:
            logging.error("The requested cycle date is not yet available in file %s " % infile)
            logging.error("Please make sure you create state based fluxes before making TC fluxes")
            raise ValueError

        # First add the date for this cycle to the file, this grows the unlimited dimension

        savedict = ncf.standard_var(varname='date')
        savedict['values'] = ncfdate
        savedict['dims'] = dimdate
        savedict['count'] = index
        ncf.add_data(savedict)

        # Now convert other variables that were inside the flux_1x1 file

        vardict = ncf_in.variables
        for vname, vprop in vardict.iteritems():

            data = ncf_in.get_variable(vname)[index]

            if vname in ['latitude','longitude', 'date', 'idate'] or 'std' in vname:
                continue
            elif 'ensemble' in vname:
                tcdata = []
                for member in data:
                    tcdata.append(statevector.vector2tc(vectordata=member))

                tcdata = np.array(tcdata)
                try:
                    cov = tcdata.transpose().dot(tcdata) / (statevector.nmembers - 1)
                except:
                    cov = np.dot(tcdata.transpose(), tcdata) / (statevector.nmembers - 1) # Huygens fix

                #print vname,cov.sum()

                tcdata = cov

                savedict = ncf.standard_var(varname=vname.replace('ensemble', 'cov'))
                savedict['units'] = '[mol/region/s]**2'
                savedict['dims'] = dimdate + dimregs + dimregs

            else:

                tcdata = statevector.vector2tc(vectordata=data) # vector to TC

                savedict = ncf.standard_var(varname=vname)
                savedict['dims'] = dimdate + dimregs
                savedict['units'] = 'mol/region/s'

            savedict['count'] = index
            savedict['values'] = tcdata
            ncf.add_data(savedict)

        ncf_in.close()
    ncf.close()

    logging.info("TransCom weekly average fluxes now written")

    return saveas

def save_weekly_avg_ext_tc_data(dacycle):
    """ Function SaveTCDataExt saves surface flux data to NetCDF files for extended TransCom regions

        *** Inputs ***
        rundat : a RunInfo object

        *** Outputs ***
        NetCDF file containing n-hourly global surface fluxes per TransCom region

        *** Example ***
        ./expand_savestate project=enkf_release sd=20000101 ed=20010101 """


#
    dirname = create_dirs(os.path.join(dacycle['dir.analysis'], 'data_tc_weekly'))
#
# Some help variables
#
    dectime0 = date2num(datetime(2000, 1, 1))
    dt = dacycle['cyclelength']
    startdate = dacycle['time.start']
    enddate = dacycle['time.end']
    ncfdate = date2num(startdate) - dectime0 + dt.days / 2.0

    logging.debug("DA Cycle start date is %s" % startdate.strftime('%Y-%m-%d %H:%M'))
    logging.debug("DA Cycle end   date is %s" % enddate.strftime('%Y-%m-%d %H:%M'))

    # Write/Create NetCDF output file
    #
    saveas = os.path.join(dirname, 'tc_extfluxes.nc')
    ncf = io.CT_CDF(saveas, 'write')
    dimdate = ncf.add_date_dim()
    dimidateformat = ncf.add_date_dim_format()
    dimregs = ncf.add_region_dim(type='tc_ext')
#
# set title and tell GMT that we are using "pixel registration"
#
    setattr(ncf, 'Title', 'CarbonTracker TransCom fluxes')
    setattr(ncf, 'node_offset', 1)
    #

    skip = ncf.has_date(ncfdate)
    if skip:
        logging.warning('Skipping writing of data for date %s : already present in file %s' % (startdate.strftime('%Y-%m-%d'), saveas))
    else:
        infile = os.path.join(dacycle['dir.analysis'], 'data_tc_weekly', 'tcfluxes.nc')
        if not os.path.exists(infile):
            logging.error("Needed input file (%s) does not exist yet, please create file first, returning..." % infile)
            return None

        ncf_in = io.ct_read(infile, 'read')

        # Transform data one by one

        # Get the date variable, and find index corresponding to the dacycle date

        try:
            dates = ncf_in.variables['date'][:]
        except KeyError:
            logging.error("The variable date cannot be found in the requested input file (%s) " % infile)
            logging.error("Please make sure you create gridded fluxes before making extended TC fluxes")
            raise KeyError

        try:
            index = dates.tolist().index(ncfdate)
        except ValueError:
            logging.error("The requested cycle date is not yet available in file %s " % infile)
            logging.error("Please make sure you create state based fluxes before making extended TC fluxes ")
            raise ValueError

        # First add the date for this cycle to the file, this grows the unlimited dimension

        savedict = ncf.standard_var(varname='date')
        savedict['values'] = ncfdate
        savedict['dims'] = dimdate
        savedict['count'] = index
        ncf.add_data(savedict)

        # Now convert other variables that were inside the tcfluxes.nc file

        vardict = ncf_in.variables
        for vname, vprop in vardict.iteritems():

            data = ncf_in.get_variable(vname)[index]

            if vname == 'latitude': continue
            elif vname == 'longitude': continue
            elif vname == 'date': continue
            elif vname == 'idate': continue
            elif 'cov' in vname:

                tcdata = ExtendedTCRegions(data, cov=True)

                savedict = ncf.standard_var(varname=vname)
                savedict['units'] = '[mol/region/s]**2'
                savedict['dims'] = dimdate + dimregs + dimregs

            else:

                tcdata = ExtendedTCRegions(data, cov=False)

                savedict = ncf.standard_var(varname=vname)
                savedict['dims'] = dimdate + dimregs
                savedict['units'] = 'mol/region/s'

            savedict['count'] = index
            savedict['values'] = tcdata
            ncf.add_data(savedict)

        ncf_in.close()
    ncf.close()

    logging.info("TransCom weekly average extended fluxes now written")

    return saveas

def save_weekly_avg_agg_data(dacycle, region_aggregate='olson'):
    """
        Function creates a NetCDF file with output on TransCom regions. It uses the flux input from the
        function `save_weekly_avg_1x1_data` to create fluxes of length `nparameters`, which are then projected
        onto TC regions using the internal methods from :class:`~da.baseclasses.statevector.StateVector`.

           :param dacycle: a :class:`~da.tools.initexit.CycleControl` object
           :param StateVector: a :class:`~da.baseclasses.statevector.StateVector`
           :rtype: None

        This function only read the prior fluxes from the flux_1x1.nc files created before, because we want to convolve
        these with the parameters in the statevector. This creates posterior fluxes, and the posterior covariance for the complete
        statevector in units of mol/box/s which we then turn into TC fluxes and covariances.
    """

#
    dirname = create_dirs(os.path.join(dacycle['dir.analysis'], 'data_%s_weekly' % region_aggregate))
#
# Some help variables
#
    dectime0 = date2num(datetime(2000, 1, 1))
    dt = dacycle['cyclelength']
    startdate = dacycle['time.start']
    enddate = dacycle['time.end']
    ncfdate = date2num(startdate) - dectime0 + dt.days / 2.0

    logging.debug("DA Cycle start date is %s" % startdate.strftime('%Y-%m-%d %H:%M'))
    logging.debug("DA Cycle end   date is %s" % enddate.strftime('%Y-%m-%d %H:%M'))

    logging.debug("Aggregating 1x1 fluxes to %s totals" % region_aggregate)


    # Write/Create NetCDF output file
    #
    saveas = os.path.join(dirname, '%s_fluxes.%s.nc' % (region_aggregate, startdate.strftime('%Y-%m-%d')))
    ncf = io.CT_CDF(saveas, 'write')
    dimdate = ncf.add_date_dim()
    dimidateformat = ncf.add_date_dim_format()
    dimgrid = ncf.add_latlon_dim()  # for mask
#
#   Select regions to aggregate to
#

    if region_aggregate == "olson":
        regionmask = tc.olson240mask
        dimname = 'olson'
        dimregs = ncf.add_dim(dimname, regionmask.max())

        regionnames = []
        for i in range(11):
            for j in range(19):
                regionnames.append("%s_%s" % (tc.transnams[i], tc.olsonnams[j],))
        regionnames.extend(tc.oifnams)
        xform = False

        for i, name in enumerate(regionnames):
            lab = 'Aggregate_Region_%03d' % (i + 1,)
            setattr(ncf, lab, name)

    elif region_aggregate == "olson_extended":
        regionmask = tc.olson_ext_mask
        dimname = 'olson_ext'
        dimregs = ncf.add_dim(dimname, regionmask.max())
        xform = False

        for i, name in enumerate(tc.olsonextnams):
            lab = 'Aggreate_Region_%03d'%(i+1)
            setattr(ncf, lab, name)

    elif region_aggregate == "transcom":
        regionmask = tc.transcommask
        dimname = 'tc'
        dimregs = ncf.add_region_dim(type='tc')
        xform = False

    elif region_aggregate == "transcom_extended":
        regionmask = tc.transcommask
        dimname = 'tc_ext'
        dimregs = ncf.add_region_dim(type='tc_ext')
        xform = True


    elif region_aggregate == "country":

        xform = False
        countrydict = ct.get_countrydict()
        selected = ['Russia', 'Canada', 'China', 'United States', 'EU27', 'Brazil', 'Australia', 'India'] #,'G8','UNFCCC_annex1','UNFCCC_annex2']
        regionmask = np.zeros((180, 360,), 'float')

        for i, name in enumerate(selected):
            lab = 'Country_%03d' % (i + 1,)
            setattr(ncf, lab, name)

            if name == 'EU27':
                namelist = ct.EU27
            elif name == 'EU25':
                namelist = ct.EU25
            elif name == 'G8':
                namelist = ct.G8
            elif name == 'UNFCCC_annex1':
                namelist = ct.annex1
            elif name == 'UNFCCC_annex2':
                namelist = ct.annex2
            else:
                namelist = [name]

            for countryname in namelist:
                try:
                    country = countrydict[countryname]
                    regionmask.put(country.gridnr, i + 1)
                except:
                    continue

        dimname = 'country'
        dimregs = ncf.add_dim(dimname, regionmask.max())

    #

    skip = ncf.has_date(ncfdate)
    if skip:
        logging.warning('Skipping writing of data for date %s : already present in file %s' % (startdate.strftime('%Y-%m-%d'), saveas))
    else:
        #
        # set title and tell GMT that we are using "pixel registration"
        #
        setattr(ncf, 'Title', 'CTDAS Aggregated fluxes')
        setattr(ncf, 'node_offset', 1)

        savedict = ncf.standard_var('unknown')
        savedict['name'] = 'regionmask'
        savedict['comment'] = 'numerical mask used to aggregate 1x1 flux fields, each integer 0,...,N is one region aggregated'
        savedict['values'] = regionmask.tolist()
        savedict['units'] = '-'
        savedict['dims'] = dimgrid
        savedict['count'] = 0
        ncf.add_data(savedict)

        # Get input data from 1x1 degree flux files

        area = globarea()

        infile = os.path.join(dacycle['dir.analysis'], 'data_flux1x1_weekly', 'flux_1x1.%s.nc' % startdate.strftime('%Y-%m-%d'))
        if not os.path.exists(infile):
            logging.error("Needed input file (%s) does not exist yet, please create file first, returning..." % infile)
            return None

        ncf_in = io.ct_read(infile, 'read')

        # Transform data one by one

        # Get the date variable, and find index corresponding to the dacycle date

        try:
            dates = ncf_in.variables['date'][:]
        except KeyError:
            logging.error("The variable date cannot be found in the requested input file (%s) " % infile)
            logging.error("Please make sure you create gridded fluxes before making TC fluxes ")
            raise KeyError

        try:
            index = dates.tolist().index(ncfdate)
        except ValueError:
            logging.error("The requested cycle date is not yet available in file %s " % infile)
            logging.error("Please make sure you create state based fluxes before making TC fluxes ")
            raise ValueError

        # First add the date for this cycle to the file, this grows the unlimited dimension

        savedict = ncf.standard_var(varname='date')
        savedict['values'] = ncfdate
        savedict['dims'] = dimdate
        savedict['count'] = index
        ncf.add_data(savedict)

        # Now convert other variables that were inside the statevector file

        vardict = ncf_in.variables
        for vname, vprop in vardict.iteritems():
            if vname == 'latitude': continue
            elif vname == 'longitude': continue
            elif vname == 'date': continue
            elif vname == 'idate': continue
            elif 'std' in vname: continue
            elif 'ensemble' in vname:

                data = ncf_in.get_variable(vname)[index]

                dimensemble = ncf.add_dim('members', data.shape[0])

                regiondata = []
                for member in data:
                    aggdata = state_to_grid(member * area, regionmask, reverse=True, mapname=region_aggregate)
                    regiondata.append(aggdata)

                regiondata = np.array(regiondata)
                try:
                    regioncov = regiondata.transpose().dot(regiondata) / (data.shape[0] - 1)
                except:
                    regioncov = np.dot(regiondata.transpose(), regiondata) / (data.shape[0] - 1) # Huygens fix

                if xform:
                    regiondata = ExtendedTCRegions(regiondata,cov=False)
                    regioncov = ExtendedTCRegions(regioncov,cov=True)

                savedict = ncf.standard_var(varname=vname)
                savedict['name'] = vname.replace('ensemble','covariance')
                savedict['units'] = '[mol/region/s]^2'
                savedict['dims'] = dimdate + dimregs + dimregs
                savedict['count'] = index
                savedict['values'] = regioncov
                ncf.add_data(savedict)

                savedict = ncf.standard_var(varname=vname)
                savedict['name'] = vname
                savedict['units'] = 'mol/region/s'
                savedict['dims'] = dimdate + dimensemble + dimregs


            elif 'flux' in vname:

                data = ncf_in.get_variable(vname)[index]

                regiondata = state_to_grid(data * area, regionmask, reverse=True, mapname=region_aggregate)

                if xform:
                    regiondata = ExtendedTCRegions(regiondata)

                savedict = ncf.standard_var(varname=vname)
                savedict['dims'] = dimdate + dimregs
                savedict['units'] = 'mol/region/s'

            else:

                data = ncf_in.get_variable(vname)[:]
                regiondata = state_to_grid(data, regionmask, reverse=True, mapname=region_aggregate)
                if xform:
                    regiondata = ExtendedTCRegions(regiondata)

                savedict = ncf.standard_var(varname=vname)
                savedict['dims'] = dimdate + dimregs

            savedict['count'] = index
            savedict['values'] = regiondata
            ncf.add_data(savedict)

        ncf_in.close()
    ncf.close()

    logging.info("%s aggregated weekly average fluxes now written" % dimname)

    return saveas

def save_time_avg_data(dacycle, infile, avg='monthly'):
    """ Function saves time mean surface flux data to NetCDF files

        *** Inputs ***
        rundat : a RunInfo object

        *** Outputs ***
        daily NetCDF file containing 1-hourly global surface fluxes at 1x1 degree

        *** Example ***
        ./expand_savestate project=enkf_release sd=20000101 ed=20010101 """

    if 'weekly' in infile:
        intime = 'weekly'
    if 'monthly' in infile:
        intime = 'monthly'
    if 'yearly' in infile:
        intime = 'yearly'

    dirname, filename = os.path.split(infile)
    outdir = create_dirs(os.path.join(dacycle['dir.analysis'], dirname.replace(intime, avg)))

    dectime0 = date2num(datetime(2000, 1, 1))

# Create NetCDF output file
#
    saveas = os.path.join(outdir, filename)
    ncf = io.CT_CDF(saveas, 'create')
    dimdate = ncf.add_date_dim()
#
# Open input file specified from the command line
#
    if not os.path.exists(infile):
        logging.error("Needed input file (%s) not found. Please create this first:" % infile)
        logging.error("returning...")
        return None
    else:
        pass

    file = io.ct_read(infile, 'read')
    datasets = file.variables.keys()
    date = file.get_variable('date')
    globatts = file.ncattrs()

    for att in globatts:
        attval = file.getncattr(att)
        if not att in ncf.ncattrs():
            ncf.setncattr(att, attval)


    time = [datetime(2000, 1, 1) + timedelta(days=d) for d in date]

# loop over datasets in infile, skip idate and date as we will make new time axis for the averaged data

    for sds in ['date'] + datasets:

# get original data

        data = file.get_variable(sds)
        varatts = file.variables[sds].ncattrs()
        vardims = file.variables[sds].dimensions
#
# Depending on dims of input dataset, create dims for output dataset. Note that we add the new dimdate now.
#

        for d in vardims:
            if 'date' in d:
                continue
            if d in ncf.dimensions.keys():
                pass
            else:
                dim = ncf.createDimension(d, size=len(file.dimensions[d]))

        savedict = ncf.standard_var(sds)
        savedict['name'] = sds
        savedict['dims'] = vardims
        savedict['units'] = file.variables[sds].units
        savedict['long_name'] = file.variables[sds].long_name
        savedict['comment'] = file.variables[sds].comment
        savedict['standard_name'] = file.variables[sds].standard_name
        savedict['count'] = 0

        if not 'date' in vardims:
            savedict['values'] = data
            ncf.add_data(savedict)
        else:

            if avg == 'monthly':
                time_avg, data_avg = timetools.monthly_avg(time, data)
            elif avg == 'seasonal':
                time_avg, data_avg = timetools.season_avg(time, data)
            elif avg == 'yearly':
                time_avg, data_avg = timetools.yearly_avg(time, data)
            elif avg == 'longterm':
                time_avg, data_avg = timetools.longterm_avg(time, data)
                time_avg = [time_avg]
                data_avg = [data_avg]
            else:
                raise ValueError, 'Averaging (%s) does not exist' % avg

            count = -1
            for dd, data in zip(time_avg, data_avg):
                count = count + 1
                if sds == 'date':
                    savedict['values'] = date2num(dd) - dectime0
                else:
                    savedict['values'] = data
                savedict['count'] = count
                ncf.add_data(savedict, silent=True)

                sys.stdout.write('.')

            sys.stdout.write('\n')
            sys.stdout.flush()

# end NetCDF file access
    file.close()
    ncf.close()

    logging.info("------------------- Finished time averaging---------------------------------")

    return saveas

if __name__ == "__main__":
    from da.tools.initexit import CycleControl
    from da.carbondioxide.dasystem import CO2DaSystem
    from da.carbondioxide.statevector import CO2StateVector

    sys.path.append('../../')

    logging.root.setLevel(logging.DEBUG)

    dacycle = CycleControl(args={'rc':'../../ctdas-stilt-proto.rc'})
    dasystem = CO2DaSystem('carbontracker_griddedNA.rc')
    dacycle.dasystem = dasystem
    dacycle.setup()
    dacycle.parse_times()



    statevector = CO2StateVector()
    statevector.setup(dacycle)

    while dacycle['time.start'] < dacycle['time.finish']:
        save_weekly_avg_1x1_data(dacycle, statevector)
        save_weekly_avg_state_data(dacycle, statevector)
        save_weekly_avg_tc_data(dacycle, statevector)
        save_weekly_avg_ext_tc_data(dacycle)
#        save_weekly_avg_agg_data(dacycle, region_aggregate='olson')
#        save_weekly_avg_agg_data(dacycle, region_aggregate='olson_extended')
        save_weekly_avg_agg_data(dacycle, region_aggregate='transcom')
        save_weekly_avg_agg_data(dacycle, region_aggregate='transcom_extended')
#        save_weekly_avg_agg_data(dacycle, region_aggregate='country')

        dacycle.advance_cycle_times()

    statevector = None # free memory

    sys.exit(0)

