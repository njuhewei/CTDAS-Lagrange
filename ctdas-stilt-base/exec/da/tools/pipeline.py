#!/usr/bin/env python
# pipeline.py

"""
.. module:: pipeline
.. moduleauthor:: Wouter Peters 

Revision History:
File created on 06 Sep 2010.

The pipeline module holds methods that execute consecutive tasks with each of the objects of the DA system. 

"""
import logging
import os
import sys
import datetime
import copy

header = '\n\n    ***************************************   '
footer = '    *************************************** \n  '

def ensemble_smoother_pipeline(dacycle, platform, dasystem, samples, statevector, obsoperator, optimizer):
    """ The main point of entry for the pipeline """
    sys.path.append(os.getcwd())

    logging.info(header + "Initializing current cycle" + footer)
    start_job(dacycle, dasystem, platform, statevector, samples, obsoperator)

    prepare_state(dacycle, statevector)
    sample_state(dacycle, samples, statevector, obsoperator)
    
    invert(dacycle, statevector, optimizer)
    
    advance(dacycle, samples, statevector, obsoperator)
        
    save_and_submit(dacycle, statevector)    
    logging.info("Cycle finished...exiting pipeline")

def forward_pipeline(dacycle, platform, dasystem, samples, statevector, obsoperator):
    """ The main point of entry for the pipeline """
    sys.path.append(os.getcwd())

    logging.info(header + "Initializing current cycle" + footer)
    start_job(dacycle, dasystem, platform, statevector, samples, obsoperator)                   

    if dacycle.has_key('forward.savestate.exceptsam'):
        sam = (dacycle['forward.savestate.exceptsam'].upper() in ["TRUE","T","YES","Y"])
    else:
        sam = False

    if dacycle.has_key('forward.savestate.dir'):
        fwddir = dacycle['forward.savestate.dir']
    else:
        logging.debug("No forward.savestate.dir key found in rc-file, proceeding with self-constructed prior parameters")
        fwddir = False

    if dacycle.has_key('forward.savestate.legacy'):
        legacy = (dacycle['forward.savestate.legacy'].upper() in ["TRUE","T","YES","Y"])
    else:
        legacy = False
        logging.debug("No forward.savestate.legacy key found in rc-file")
    
    if not fwddir:
        # Simply make a prior statevector using the normal method


        prepare_state(dacycle, statevector)#LU tutaj zamiast tego raczej to stworzenie nowej kowariancji i ensembli bo pozostale rzeczy sa na gorze i na doel.

    else: 
        # Read prior information from another simulation into the statevector.
        # This loads the results from another assimilation experiment into the current statevector

        if sam:
            filename = os.path.join(fwddir, dacycle['time.start'].strftime('%Y%m%d'), 'savestate.nc')
            statevector.read_from_file_exceptsam(filename, 'prior') 
        elif not legacy:
            filename = os.path.join(fwddir, dacycle['time.start'].strftime('%Y%m%d'), 'savestate_%s.nc'%dacycle['time.start'].strftime('%Y%m%d'))        
            statevector.read_from_file(filename, 'prior') 
        else:
            filename = os.path.join(fwddir, dacycle['time.start'].strftime('%Y%m%d'), 'savestate.nc')
            statevector.read_from_legacy_file(filename, 'prior') 


    # We write this "prior" statevector to the restart directory, so we can later also populate it with the posterior statevector
    # Note that we could achieve the same by just copying the wanted forward savestate.nc file to the restart folder of our current
    # experiment, but then it would already contain a posterior field as well which we will try to write in save_and_submit. 
    # This could cause problems. Moreover, this method allows us to read older formatted savestate.nc files (legacy) and write them into
    # the current format through the "write_to_file" method.

    savefilename = os.path.join(dacycle['dir.restart'], 'savestate_%s.nc' % dacycle['time.start'].strftime('%Y%m%d'))
    statevector.write_to_file(savefilename, 'prior')

    # Now read optimized fluxes which we will actually use to propagate through the system
    
    if not fwddir:

        # if there is no forward dir specified, we simply run forward with unoptimized prior fluxes in the statevector
        logging.info("Running forward with prior savestate from: %s"%savefilename)

    else: 

        # Read posterior information from another simulation into the statevector.
        # This loads the results from another assimilation experiment into the current statevector

        if sam:
            statevector.read_from_file_exceptsam(filename, 'opt') 
        elif not legacy:
            statevector.read_from_file(filename, 'opt') 
        else:
            statevector.read_from_legacy_file(filename, 'opt') 

        logging.info("Running forward with optimized savestate from: %s"%filename)

    # Finally, we run forward with these parameters

    advance(dacycle, samples, statevector, obsoperator)

    # In save_and_submit, the posterior statevector will be added to the savestate.nc file, and it is added to the copy list.
    # This way, we have both the prior and posterior data from another run copied into this assimilation, for later analysis.

    save_and_submit(dacycle, statevector) 

    logging.info("Cycle finished...exiting pipeline")
####################################################################################################

def analysis_pipeline(dacycle, platform, dasystem, samples, statevector):
    """ Main entry point for analysis of ctdas results """

    from da.analysis.expand_fluxes import save_weekly_avg_1x1_data, save_weekly_avg_state_data, save_weekly_avg_tc_data, save_weekly_avg_ext_tc_data, save_weekly_avg_agg_data
    from da.analysis.expand_molefractions import write_mole_fractions
    from da.analysis.summarize_obs import summarize_obs
    from da.analysis.time_avg_fluxes import time_avg

    logging.info(header + "Starting analysis" + footer) 

    dasystem.validate()                               
    dacycle.dasystem = dasystem                       
    dacycle.daplatform = platform                    
    dacycle.setup()                              
    statevector.setup(dacycle)   

    logging.info(header + "Starting mole fractions" + footer) 

    write_mole_fractions(dacycle)
    summarize_obs(dacycle['dir.analysis'])

    logging.info(header + "Starting weekly averages" + footer) 

    save_weekly_avg_1x1_data(dacycle, statevector)
    save_weekly_avg_state_data(dacycle, statevector)
    save_weekly_avg_tc_data(dacycle, statevector)
    save_weekly_avg_ext_tc_data(dacycle)
    save_weekly_avg_agg_data(dacycle,region_aggregate='transcom')
    save_weekly_avg_agg_data(dacycle,region_aggregate='transcom_extended')
    save_weekly_avg_agg_data(dacycle,region_aggregate='olson')
    save_weekly_avg_agg_data(dacycle,region_aggregate='olson_extended')
    save_weekly_avg_agg_data(dacycle,region_aggregate='country')

    logging.info(header + "Starting monthly and yearly averages" + footer) 

    time_avg(dacycle,'flux1x1')
    time_avg(dacycle,'transcom')
    time_avg(dacycle,'transcom_extended')
    time_avg(dacycle,'olson')
    time_avg(dacycle,'olson_extended')
    time_avg(dacycle,'country')

    logging.info(header + "Finished analysis" + footer) 

def archive_pipeline(dacycle, platform, dasystem):
    """ Main entry point for archiving of output from one disk/system to another """

    if not dacycle.has_key('task.rsync'):
        logging.info('rsync task not found, not starting automatic backup...')
        return
    else:
        logging.info('rsync task found, starting automatic backup...')

    for task in dacycle['task.rsync'].split():
        sourcedirs = dacycle['task.rsync.%s.sourcedirs'%task]
        destdir = dacycle['task.rsync.%s.destinationdir'%task]

        rsyncflags = dacycle['task.rsync.%s.flags'%task]

        # file ID and names
        jobid = dacycle['time.end'].strftime('%Y%m%d') 
        targetdir = os.path.join(dacycle['dir.exec'])
        jobfile = os.path.join(targetdir, 'jb.rsync.%s.%s.jb' % (task,jobid) )
        logfile = os.path.join(targetdir, 'jb.rsync.%s.%s.log' % (task,jobid) )
        # Template and commands for job
        jobparams = {'jobname':"r.%s" % jobid, 'jobnodes': '1', 'jobtime': '1:00:00', 'joblog': logfile, 'errfile': logfile}

        if platform.ID == 'cartesius':
            jobparams['jobqueue'] = 'staging'

        template = platform.get_job_template(jobparams)
        for sourcedir in sourcedirs.split():
            execcommand = "\nrsync %s %s %s\n" % (rsyncflags, sourcedir,destdir,)
            template += execcommand

        # write and submit 
        platform.write_job(jobfile, template, jobid)
        jobid = platform.submit_job(jobfile, joblog=logfile)


def start_job(dacycle, dasystem, platform, statevector, samples, obsoperator):
    """ Set up the job specific directory structure and create an expanded rc-file """

    dasystem.validate()                               
    dacycle.dasystem = dasystem                       
    dacycle.daplatform = platform                    
    dacycle.setup()                              
    #statevector.dacycle = dacycle # also embed object in statevector so it can access cycle information for I/O etc   
    #samples.dacycle = dacycle # also embed object in samples object so it can access cycle information for I/O etc    
    #obsoperator.dacycle = dacycle # also embed object in obsoperator object so it can access cycle information for I/O etc
    obsoperator.setup(dacycle)  # Setup Observation Operator                                                          
    statevector.setup(dacycle)   

def prepare_state(dacycle, statevector):
    """ Set up the input data for the forward model: obs and parameters/fluxes"""

    # We now have an empty statevector object that we need to populate with data. If this is a continuation from a previous cycle, we can read
    # the previous statevector values from a NetCDF file in the restart directory. If this is the first cycle, we need to populate the statevector
    # with new values for each week. After we have constructed the statevector, it will be propagated by one cycle length so it is ready to be used
    # in the current cycle

    logging.info(header + "starting prepare_state" + footer)

    if not dacycle['time.restart']:                    

        # Fill each week from n=1 to n=nlag with a new ensemble

        for n in range(statevector.nlag):
            date = dacycle['time.start'] + datetime.timedelta(days=(n + 0.5) * int(dacycle['time.cycle']))            
            cov = statevector.get_covariance(date, dacycle)
            statevector.make_new_ensemble(n, cov)

    else:

        # Read the statevector data from file

        #saved_sv = os.path.join(dacycle['dir.restart.current'], 'savestate.nc')               
        
        
        saved_sv = os.path.join(dacycle['dir.restart'], 'savestate_%s.nc' % dacycle['da.restart.tstamp'].strftime('%Y%m%d'))
        statevector.read_from_file(saved_sv) # by default will read "opt"(imized) variables, and then propagate

        # Now propagate the ensemble by one cycle to prepare for the current cycle
        statevector.propagate(dacycle)

    # Finally, also write the statevector to a file so that we can always access the a-priori information

    current_sv = os.path.join(dacycle['dir.restart'], 'savestate_%s.nc' % dacycle['time.start'].strftime('%Y%m%d'))
    statevector.write_to_file(current_sv, 'prior')  # write prior info 

def sample_state(dacycle, samples, statevector, obsoperator):
    """ Sample the filter state for the inversion """


    # Before a forecast step, save all the data to a save/tmp directory so we can later recover it before the propagation step.
    # This is especially important for:
    #  (i) The transport model restart data which holds the background mole fractions. This is needed to run the model one cycle forward
    #  (ii) The random numbers (or the seed for the random number generator) so we can recreate the ensembles if needed

    #status   = dacycle.MoveSaveData(io_option='store',save_option='partial',filter=[])
    #msg     = "All restart data have been copied to the save/tmp directory for future use"    ; logging.debug(msg)
    logging.info(header + "starting sample_state" + footer) 
    nlag = int(dacycle['time.nlag'])
    logging.info("Sampling model will be run over %d cycles" % nlag)           

    obsoperator.get_initial_data()

    for lag in range(nlag):                                                               
        logging.info(header + ".....Ensemble Kalman Filter at lag %d" % (lag + 1)) 

        ############# Perform the actual sampling loop #####################

        sample_step(dacycle, samples, statevector, obsoperator, lag)

        logging.debug("statevector now carries %d samples" % statevector.nobs)


def sample_step(dacycle, samples, statevector, obsoperator, lag, advance=False):
    """ Perform all actions needed to sample one cycle """
    

    # First set up the information for time start and time end of this sample
    dacycle.set_sample_times(lag)

    startdate = dacycle['time.sample.start'] 
    enddate = dacycle['time.sample.end'] 
    dacycle['time.sample.window'] = lag
    dacycle['time.sample.stamp'] = "%s_%s" % (startdate.strftime("%Y%m%d%H"), enddate.strftime("%Y%m%d%H"))

    logging.info("New simulation interval set : ")
    logging.info("                  start date : %s " % startdate.strftime('%F %H:%M'))
    logging.info("                  end   date : %s " % enddate.strftime('%F %H:%M'))
    logging.info("                  file  stamp: %s " % dacycle['time.sample.stamp'])


    # Implement something that writes the ensemble member parameter info to file, or manipulates them further into the 
    # type of info needed in your transport model

    statevector.write_members_to_file(lag, dacycle['dir.input']) 

    samples.setup(dacycle)       
    samples.add_observations() 

    # Add model-data mismatch to all samples, this *might* use output from the ensemble in the future??

    samples.add_model_data_mismatch(dacycle.dasystem['obs.sites.rc']) 
    
    sampling_coords_file = os.path.join(dacycle['dir.input'], 'sample_coordinates_%s.nc' % dacycle['time.sample.stamp'])
    samples.write_sample_coords(sampling_coords_file) 
    # Write filename to dacycle, and to output collection list
    dacycle['ObsOperator.inputfile'] = sampling_coords_file

    # Run the observation operator

    obsoperator.run_forecast_model()


    # Read forecast model samples that were written to NetCDF files by each member. Add them to the exisiting
    # Observation object for each sample loop. This data fill be written to file in the output folder for each sample cycle. 
    

    # We retrieve all model samples from one output file written by the ObsOperator. If the ObsOperator creates
    # one file per member, some logic needs to be included to merge all files!!!

    if os.path.exists(sampling_coords_file):
        samples.add_simulations(obsoperator.simulated_file)
    else: logging.warning("No simulations added, because input file does not exist (no samples found in obspack)")

    # Now write a small file that holds for each observation a number of values that we need in postprocessing

    #filename = samples.write_sample_coords() 

    # Now add the observations that need to be assimilated to the statevector. 
    # Note that obs will only be added to the statevector if either this is the first step (restart=False), or lag==nlag
    # This is to make sure that the first step of the system uses all observations available, while the subsequent
    # steps only optimize against the data at the front (lag==nlag) of the filter. This way, each observation is used only 
    # (and at least) once # in the assimilation


    if not advance:
        if dacycle['time.restart'] == False or lag == int(dacycle['time.nlag']) - 1:
            statevector.obs_to_assimilate += (copy.deepcopy(samples),)
            statevector.nobs += samples.getlength()
            logging.debug("Added samples from the observation operator to the assimilated obs list in the statevector")

        else:
            statevector.obs_to_assimilate += (None,)


def invert(dacycle, statevector, optimizer):
    """ Perform the inverse calculation """
    logging.info(header + "starting invert" + footer)
    dims = (int(dacycle['time.nlag']),
                  int(dacycle['da.optimizer.nmembers']),
                  int(dacycle.dasystem['nparameters']),
                  statevector.nobs)

    if not dacycle.dasystem.has_key('opt.algorithm'):
        logging.info("There was no minimum least squares algorithm specified in the DA System rc file (key : opt.algorithm)") 
        logging.info("...using serial algorithm as default...")
        optimizer.set_algorithm('Serial')
    elif dacycle.dasystem['opt.algorithm'] == 'serial':
        logging.info("Using the serial minimum least squares algorithm to solve ENKF equations")
        optimizer.set_algorithm('Serial')
    elif dacycle.dasystem['opt.algorithm'] == 'bulk':
        logging.info("Using the bulk minimum least squares algorithm to solve ENKF equations")
        optimizer.set_algorithm('Bulk')
    
    optimizer.setup(dims)
    optimizer.state_to_matrix(statevector)    
    
    diagnostics_file = os.path.join(dacycle['dir.output'], 'optimizer.%s.nc' % dacycle['time.start'].strftime('%Y%m%d'))
        
    optimizer.write_diagnostics(diagnostics_file, 'prior')
    optimizer.set_localization(dacycle['da.system.localization'])
    
    if optimizer.algorithm == 'Serial':
        optimizer.serial_minimum_least_squares()
    else: 
        optimizer.bulk_minimum_least_squares()
    
    optimizer.matrix_to_state(statevector)
    optimizer.write_diagnostics(diagnostics_file, 'optimized')
    
    

def advance(dacycle, samples, statevector, obsoperator):
    """ Advance the filter state to the next step """

    # This is the advance of the modeled CO2 state. Optionally, routines can be added to advance the state vector (mean+covariance)

    # Then, restore model state from the start of the filter
    logging.info(header + "starting advance" + footer)
    logging.info("Sampling model will be run over 1 cycle")

    obsoperator.get_initial_data()

    sample_step(dacycle, samples, statevector, obsoperator, 0, True) 

    dacycle.restart_filelist.extend(obsoperator.restart_filelist)
    dacycle.output_filelist.extend(obsoperator.output_filelist)
    logging.debug("Appended ObsOperator restart and output file lists to dacycle for collection ")
    
    dacycle.output_filelist.append(dacycle['ObsOperator.inputfile'])
    logging.debug("Appended Observation filename to dacycle for collection ")

    sampling_coords_file = os.path.join(dacycle['dir.input'], 'sample_coordinates_%s.nc' % dacycle['time.sample.stamp'])
    if os.path.exists(sampling_coords_file):
        outfile = os.path.join(dacycle['dir.output'], 'sample_auxiliary_%s.nc' % dacycle['time.sample.stamp'])
        samples.write_sample_auxiliary(outfile)
    else: logging.warning("Sample auxiliary output not written, because input file does not exist (no samples found in obspack)")

def save_and_submit(dacycle, statevector):
    """ Save the model state and submit the next job """
    logging.info(header + "starting save_and_submit" + footer)
    
    filename = os.path.join(dacycle['dir.restart'], 'savestate_%s.nc' % dacycle['time.start'].strftime('%Y%m%d'))
    statevector.write_to_file(filename, 'opt')
    
    dacycle.output_filelist.append(filename)
    dacycle.finalize()




