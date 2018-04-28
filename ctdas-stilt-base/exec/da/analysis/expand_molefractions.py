#!/usr/bin/env python
# expand_fluxes.py
import sys
import os
import getopt
import shutil
import logging
import netCDF4
import numpy as np
from string import join
from datetime import datetime, timedelta
sys.path.append('../../')
from da.tools.general import date2num, num2date
from da.tools.general import create_dirs
import da.tools.io4 as io


"""
Author: Wouter Peters (Wouter.Peters@wur.nl)

Revision History:
File created on 11 May 2012.

"""

def write_mole_fractions(dacycle):
    """ 
    
    Write Sample information to NetCDF files. These files are organized by site and 
    have an unlimited time axis to which data is appended each cycle.

    The needed information is obtained from the sample_auxiliary.nc files and the original input data files from ObsPack. 

    The steps are:

    (1) Create a directory to hold timeseries output files
    (2) Read the sample_auxiliary.nc file for this cycle and get a list of original files they were obtained from
    (3) For each file, copy the original data file from ObsPack (if not yet present)
    (4) Open the copied file, find the index of each observation, fill in the simulated data
    
    """

    dirname = create_dirs(os.path.join(dacycle['dir.analysis'], 'data_molefractions'))
#
# Some help variables
#
    dectime0 = date2num(datetime(2000, 1, 1))
    dt = dacycle['cyclelength']
    startdate = dacycle['time.start'] 
    enddate = dacycle['time.end'] 

    logging.debug("DA Cycle start date is %s" % startdate.strftime('%Y-%m-%d %H:%M'))
    logging.debug("DA Cycle end   date is %s" % enddate.strftime('%Y-%m-%d %H:%M'))

    dacycle['time.sample.stamp'] = "%s_%s" % (startdate.strftime("%Y%m%d%H"), enddate.strftime("%Y%m%d%H"),)

    # Step (1): Get the posterior sample output data file for this cycle

    infile = os.path.join(dacycle['dir.output'], 'sample_auxiliary_%s.nc' % dacycle['time.sample.stamp'])

    ncf_in = io.ct_read(infile, 'read')

    obs_num = ncf_in.get_variable('obs_num')
    obs_val = ncf_in.get_variable('observed')
    simulated = ncf_in.get_variable('modelsamples')
    infilename = ncf_in.get_variable('inputfilename')
    infiles1 = netCDF4.chartostring(infilename).tolist()
    # In case of reanalysis on different platform, obspack-input-directory might have a different name. 
    # This is checked here, and the filenames are corrected
    dir_from_rc = dacycle.dasystem['obspack.input.dir']
    dir_from_output = infiles1[0]
    d1 = dir_from_rc[:dir_from_rc.find('obspacks')]
    d2 = dir_from_output[:dir_from_output.find('obspacks')]
    if d1 == d2:
        infiles = infiles1
    else:
        infiles = []
        for ff in infiles1:
            infiles.append(ff.replace(d2,d1))

    #infiles   = [join(s.compressed(),'') for s in infilename]

    ncf_in.close()

    # Step (2): Get the prior sample output data file for this cycle

    infile = os.path.join(dacycle['dir.output'], 'optimizer.%s.nc' % startdate.strftime('%Y%m%d'))

    if os.path.exists(infile): 
        optimized_present = True
    else:
        optimized_present = False

    if optimized_present:

        ncf_fc_in = io.ct_read(infile, 'read')

        fc_obs_num = ncf_fc_in.get_variable('obspack_num')
        fc_obs_val = ncf_fc_in.get_variable('observed')
        fc_simulated = ncf_fc_in.get_variable('modelsamplesmean_prior')
        fc_simulated_ens = ncf_fc_in.get_variable('modelsamplesdeviations_prior')
        fc_flag      = ncf_fc_in.get_variable('flag')
        if not dacycle.dasystem.has_key('opt.algorithm'):
            fc_r         = ncf_fc_in.get_variable('modeldatamismatchvariance')
            fc_hphtr     = ncf_fc_in.get_variable('totalmolefractionvariance')
        elif dacycle.dasystem['opt.algorithm'] == 'serial':
            fc_r         = ncf_fc_in.get_variable('modeldatamismatchvariance')
            fc_hphtr     = ncf_fc_in.get_variable('totalmolefractionvariance')
        elif dacycle.dasystem['opt.algorithm'] == 'bulk':
            fc_r         = ncf_fc_in.get_variable('modeldatamismatchvariance').diagonal()
            fc_hphtr     = ncf_fc_in.get_variable('totalmolefractionvariance').diagonal()
        filesitecode = ncf_fc_in.get_variable('sitecode')

        fc_sitecodes = netCDF4.chartostring(filesitecode).tolist()
        #fc_sitecodes = [join(s.compressed(),'') for s in filesitecode]

        ncf_fc_in.close()

        # Expand the list of input files with those available from the forecast list

        infiles_rootdir = os.path.split(infiles[0])[0]
        infiles.extend(os.path.join(infiles_rootdir, f + '.nc') for f in fc_sitecodes)


    #Step (2): For each observation timeseries we now have data for, open it and fill with data

    for orig_file in set(infiles):

        if not os.path.exists(orig_file):
            logging.error("The original input file (%s) could not be found, continuing to next file..." % orig_file)
            continue


        copy_file = os.path.join(dirname, os.path.split(orig_file)[-1])
        if not os.path.exists(copy_file):
            shutil.copy(orig_file, copy_file)
            logging.debug("Copied a new original file (%s) to the analysis directory" % orig_file)

            ncf_out = io.CT_CDF(copy_file, 'write')

            # Modify the attributes of the file to reflect added data from CTDAS properly

	    try:
	       host=os.environ['HOSTNAME']
	    except:
	       host='unknown'
     

            ncf_out.Caution = '==================================================================================='
            try:
                ncf_out.History += '\nOriginal observation file modified by user %s on %s\n' % (os.environ['USER'], datetime.today().strftime('%F'),)
            except:
                ncf_out.History = '\nOriginal observation file modified by user %s on %s\n' % (os.environ['USER'], datetime.today().strftime('%F'),)
            ncf_out.CTDAS_info = 'Simulated values added from a CTDAS run by %s on %s\n' % (os.environ['USER'], datetime.today().strftime('%F'),)\
                               + '\nCTDAS was run on platform %s' % (host,)\
                               + '\nCTDAS job directory was %s' % (dacycle['dir.da_run'],)\
                               + '\nCTDAS Da System was %s' % (dacycle['da.system'],)\
                               + '\nCTDAS Da ObsOperator was %s' % (dacycle['da.obsoperator'],)
            ncf_out.CTDAS_startdate = dacycle['time.start'].strftime('%F')
            ncf_out.CTDAS_enddate = dacycle['time.finish'].strftime("%F")
            ncf_out.original_file = orig_file


            # get nobs dimension

            if ncf_out.dimensions.has_key('id'): 
                dimidob = ncf_out.dimensions['id']
                dimid = ('id',)
            elif ncf_out.dimensions.has_key('obs'): 
                dimidob = ncf_out.dimensions['obs']
                dimid = ('obs',)

            if dimidob.isunlimited:
                nobs = ncf_out.inq_unlimlen()
            else:
                nobs = len(dimid)


            # add nmembers dimension

            dimmembersob = ncf_out.createDimension('nmembers', size=simulated.shape[1])
            dimmembers = ('nmembers',)
            nmembers = len(dimmembers)

            # Create empty arrays for posterior samples, as well as for forecast sample statistics

            savedict = io.std_savedict.copy()
            savedict['name'] = "flag_forecast"
            savedict['long_name'] = "flag_for_obs_model in forecast"
            savedict['units'] = "None"
            savedict['dims'] = dimid
            savedict['comment'] = 'Flag (0/1/2/99) for observation value, 0 means okay, 1 means QC error, 2 means rejected, 99 means not sampled'
            ncf_out.add_variable(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "modeldatamismatch"
            savedict['long_name'] = "modeldatamismatch"
            savedict['units'] = "[mol mol-1]^2"
            savedict['dims'] = dimid
            savedict['comment'] = 'Variance of mole fractions resulting from model-data mismatch'
            ncf_out.add_variable(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "totalmolefractionvariance_forecast"
            savedict['long_name'] = "totalmolefractionvariance of forecast"
            savedict['units'] = "[mol mol-1]^2"
            savedict['dims'] = dimid
            savedict['comment'] = 'Variance of mole fractions resulting from prior state and model-data mismatch'
            ncf_out.add_variable(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "modelsamplesmean"
            savedict['long_name'] = "mean modelsamples"
            savedict['units'] = "mol mol-1"
            savedict['dims'] = dimid
            savedict['comment'] = 'simulated mole fractions based on optimized state vector'
            ncf_out.add_variable(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "modelsamplesmean_forecast"
            savedict['long_name'] = "mean modelsamples from forecast"
            savedict['units'] = "mol mol-1"
            savedict['dims'] = dimid
            savedict['comment'] = 'simulated mole fractions based on prior state vector'
            ncf_out.add_variable(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "modelsamplesstandarddeviation"
            savedict['long_name'] = "standard deviaton of modelsamples over all ensemble members"
            savedict['units'] = "mol mol-1"
            savedict['dims'] = dimid
            savedict['comment'] = 'std dev of simulated mole fractions based on optimized state vector'
            ncf_out.add_variable(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "modelsamplesstandarddeviation_forecast"
            savedict['long_name'] = "standard deviaton of modelsamples from forecast over all ensemble members"
            savedict['units'] = "mol mol-1"
            savedict['dims'] = dimid
            savedict['comment'] = 'std dev of simulated mole fractions based on prior state vector'
            ncf_out.add_variable(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "modelsamplesensemble"
            savedict['long_name'] = "modelsamples over all ensemble members"
            savedict['units'] = "mol mol-1"
            savedict['dims'] = dimid + dimmembers
            savedict['comment'] = 'ensemble of simulated mole fractions based on optimized state vector'
            ncf_out.add_variable(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "modelsamplesensemble_forecast"
            savedict['long_name'] = "modelsamples from forecast over all ensemble members"
            savedict['units'] = "mol mol-1"
            savedict['dims'] = dimid + dimmembers
            savedict['comment'] = 'ensemble of simulated mole fractions based on prior state vector'
            ncf_out.add_variable(savedict)

        else:
            logging.debug("Modifying existing file (%s) in the analysis directory" % copy_file)

            ncf_out = io.CT_CDF(copy_file, 'write')


        # Get existing file obs_nums to determine match to local obs_nums

        if ncf_out.variables.has_key('id'):
            file_obs_nums = ncf_out.get_variable('id')
        elif ncf_out.variables.has_key('obspack_num'):
            file_obs_nums = ncf_out.get_variable('obspack_num')

        # Get all obs_nums related to this file, determine their indices in the local arrays

        selected_obs_nums = [num for infile, num in zip(infiles, obs_num) if infile == orig_file]


        # Optimized data 1st: For each index, get the data and add to the file in the proper file index location

        for num in selected_obs_nums:

            model_index = obs_num.tolist().index(num)
            file_index = file_obs_nums.tolist().index(num)

            #var = ncf_out.variables['modeldatamismatch']   # Take from optimizer.yyyymmdd.nc file instead
            #var[file_index] = mdm[model_index]

            var = ncf_out.variables['modelsamplesmean']
            var[file_index] = simulated[model_index, 0]

            var = ncf_out.variables['modelsamplesstandarddeviation']
            var[file_index] = simulated[model_index, 1:].std()

            var = ncf_out.variables['modelsamplesensemble']
            var[file_index] = simulated[model_index, :]

        # Now forecast data too: For each index, get the data and add to the file in the proper file index location

        if optimized_present:

            selected_fc_obs_nums = [num for sitecode, num in zip(fc_sitecodes, fc_obs_num) if sitecode in orig_file]

            for num in selected_fc_obs_nums:

                model_index = fc_obs_num.tolist().index(num)
                file_index = file_obs_nums.tolist().index(num)

                var = ncf_out.variables['modeldatamismatch']
                var[file_index] = np.sqrt(fc_r[model_index])

                var = ncf_out.variables['modelsamplesmean_forecast']
                var[file_index] = fc_simulated[model_index]

                var = ncf_out.variables['modelsamplesstandarddeviation_forecast']
                var[file_index] = fc_simulated_ens[model_index, 1:].std()

                var = ncf_out.variables['modelsamplesensemble_forecast']
                var[file_index] = fc_simulated_ens[model_index, :]

                var = ncf_out.variables['totalmolefractionvariance_forecast']
                var[file_index] = fc_hphtr[model_index]

                var = ncf_out.variables['flag_forecast']
                var[file_index] = fc_flag[model_index]


        # close the file

        status = ncf_out.close()

    return None



if __name__ == "__main__":
    import logging
    from da.tools.initexit import CycleControl
    from da.carbondioxide.dasystem import CO2DaSystem

    sys.path.append('../../')

    logging.root.setLevel(logging.DEBUG)

    dacycle = CycleControl(args={'rc':'../../ctdas-od-gfed2-glb6x4-obspack-full.rc'})
    dasystem = CO2DaSystem('../rc/carbontracker_ct09_opfnew.rc')
    dacycle.dasystem = dasystem
    dacycle.setup()
    dacycle.parse_times()



    while dacycle['time.start'] < dacycle['time.finish']:

        outfile = write_mole_fractions(dacycle)

        dacycle.advance_cycle_times()

    sys.exit(0)

