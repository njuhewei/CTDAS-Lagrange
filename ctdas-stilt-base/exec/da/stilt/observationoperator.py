#!/usr/bin/env python
# stilt_tools.py

"""
Author : W. He

Revision History:
Basing on Wouter's codes, replace the TM5 model with the STILT model, April 2015.

This module holds specific functions needed to use the STILT model within the data assimilation shell. It uses the information
from the DA system in combination with the generic stilt.rc files.

The STILT model is now controlled by a python subprocess. This subprocess consists of an MPI wrapper (written in C) that spawns
a large number ( N= nmembers) of STILT model instances under mpirun, and waits for them all to finish.

#The design of the system assumes that the stilt function are executed using R codes, and is residing in a
#directory specified by the ${RUNDIR} of a stilt rc-file. This stilt rc-file name is taken from the data assimilation rc-file. Thus,

"""
import shutil
import os
import sys
import logging
import shutil
import datetime
import subprocess    #use to get info from running status
import numpy as np
from string import join
import glob
#sys.path.append(os.getcwd())
#sys.path.append("../../")

#import da.tools.rc as rc
#from da.tools.general import create_dirs, to_datetime
#from da.baseclasses.observationoperator import ObservationOperator

# global constants, which will be used in the following classes
identifier = 'WRF-STILT'
version = '1.0'
#mpi_shell_filename = 'STILT_mpi_wrapper'   #for STILT, not use MPI
#mpi_shell_location = 'da/bin/'


################### Begin Class STILT ###################

class STILTObservationOperator(object):

    def __init__(self, dacycle=None):   #only the filename used to specify the location of the stavector file for wrf-stilt runs
        """ The instance of an STILTObservationOperator is application dependent """
        self.ID = identifier    # the identifier gives the model name
        self.version = version       # the model version used
        self.restart_filelist = []
        self.output_filelist = []
        self.outputdir = None # Needed for opening the samples.nc files created

        #self.simulated_file = None
        #self.forecast_nmembers = None

        # we need the folder of the input statevector and the path of the excecutable STILT procedure
        #self.rc=filename

        logging.info('Observation Operator initialized: %s (%s)' % (self.ID, self.version))

        if dacycle != None:
            self.dacycle = dacycle
        else:
            self.dacycle = {}

        self.startdate=None

    def get_initial_data(self):
        """ This method places all initial data needed by an ObservationOperator in the proper folder for the model """

    def setup(self, dacycle):
        """ Execute all steps needed to prepare the ObsOperator for use inside CTDAS, only done at the very first cycle normally """
        self.dacycle = dacycle
        self.outputdir = dacycle['dir.output']
        self.rc = dacycle['da.obsoperator.rc']

    def prepare_run(self,proc):
        """ Prepare the running of the actual forecast model, for example compile code """
        import os

        # Define the name of the file that will contain the modeled output of each observation
        self.simulated_file = os.path.join(self.outputdir, 'samples_simulated.%s-%s.nc' % (self.dacycle['time.sample.stamp'],proc))
        self.forecast_nmembers = int(self.dacycle['da.optimizer.nmembers'])

        self.update_rc(self.rc,proc)

    def update_rc(self,name,proc):
        if 'stilt.rc' in name:
            path,dummy= name.split('stilt.rc')
            shutil.copyfile(name,path+'stilt_%s.rc'%proc)
            name=os.path.join(path,'stilt_%s.rc'%proc)
        self.rc_filename = name


        logging.info('RC name %s'%name)
        with open(name) as f:
            data = f.read().split()

        f.close()

        data=np.array(data)
        starttime = self.dacycle['time.sample.start']
        starttime = datetime.date(starttime.year,starttime.month,starttime.day)  #YYYY-MM-DD

        data[17] = self.outputdir  #simulated_file
        data[20] = self.forecast_nmembers
        data[29] = starttime
        self.startdate=starttime
        t=data.reshape((10, 3))

        output = open(name, 'w')
        for i in range(0,t.shape[0],1) :
           output.write('%12s %s %s\n' % (t[i][0],t[i][1],t[i][2]))

        logging.debug('STILT rc-file updated successfully')


    def run_forecast_model(self,proc,out_q):   #extral interface
        self.prepare_run(proc)
        self.run(proc)

        outdict = {}
        outdict[proc] = self.simulated_file
        out_q.put(outdict)

    def run(self,proc):
        """
         Start the STILT executable. A new log file is started for the STILT model IO, and then a subprocess is
         spawned with the STILT_mpi_wrapper and the STILT.x executable. The exit code of the model is caught and
         only if successfull on all processors will execution of the shell continue.

        """
        cwd = os.getcwd()

        # (1) Where an mpi process is forked to do a STILT instance with N tracers, each an ensemble member

        # parallel running for the statevectors to get stilt-based concentration
        code = self.runstilt(proc)

        if code == 0:
            logging.info('Finished model executable succesfully (%s)' % code)
            self.Status = 'Success'
        else:
            logging.error('Error in model executable return code: %s ' % code)
            self.Status = 'Failed'
            raise OSError

        # Return to working directory

        os.chdir(cwd)

        return code

    def runstilt(self,proc):
        """
        Call stilt R executable : //give input file lists
        """
        okfile = 'stilt_%d.ok'%proc
        if os.path.exists(okfile):
            os.remove(okfile)

        Rdir=self.dacycle['da.obsoperator.home']
        args='--args %d'%proc
        submitcommand=['R','CMD','BATCH','--args %s'%proc,Rdir+'/stilt.co2.simu.2015.Sep.fpbc.r']
        logging.info('Submit command %s' %submitcommand)
        logging.info('Submitting job at %s' % datetime.datetime.now())
        code = subprocess.call(submitcommand)
        logging.info('Resuming job at %s' % datetime.datetime.now())

        if not os.path.exists(okfile):
            code = -1
        else:
            code = 0

        return code


    def save_data(self):
        """ Copy the STILT recovery data from the outputdir to the STILT savedir, also add the restart files to a list of names
            that is used by the dacycle object to collect restart data for the filter.

            WP Note: with the new pycasso restart files we no longer need to copy save files from outdir to savedir

            Note 2: also adding the weekly mean flux output to the output_filelist for later collection
         """

        sourcedir = os.path.join(self.STILT_settings[self.savedirkey])
        filterlist = ['%s' % self.STILT_settings[self.timefinalkey].strftime('%Y%m%d')]

        logging.debug("Creating a new list of STILT restart data")
        logging.debug("           from directory: %s " % sourcedir)
        logging.debug("           with filter: %s " % filterlist)


        # Start from empty lists for each STILT run. Note that these "private" lists from the obs operator are later on appended to the system
        # lists

        self.restart_filelist = []

        for fil in os.listdir(sourcedir):
            fil = os.path.join(sourcedir, fil)
            if os.path.isdir(fil): # skip dirs
                skip = True
            elif filterlist == []:      # copy all
                skip = False
            else:                   # check filter
                skip = True         # default skip
                for f in filterlist:
                    if f in fil:
                        skip = False # unless in filterlist
                        break

            if skip:
                logging.debug("           [skip] .... %s " % fil)
                continue

            self.restart_filelist.append(fil)
            logging.debug("           [added to restart list] .... %s " % fil)

        sourcedir = os.path.join(self.STILT_settings[self.outputdirkey])
        sd_ed = self.dacycle['time.sample.stamp']
        filterlist = ['flask_output.%s' % sd_ed, 'flux1x1_%s' % sd_ed]

        logging.debug("Creating a new list of STILT output data to collect")
        logging.debug("           from directory: %s " % sourcedir)
        logging.debug("           with filter: %s " % filterlist)


        # Start from empty lists for each STILT run. Note that these "private" lists from the obs operator are later on appended to the system
        # lists

        self.output_filelist = []

        for fil in os.listdir(sourcedir):
            fil = os.path.join(sourcedir, fil)

            if os.path.isdir(fil): # skip dirs
                skip = True
            elif filterlist == []:      # copy all
                skip = False
            else:                   # check filterlist
                skip = True         # default skip
                for f in filterlist:
                    if f in fil:
                        skip = False # unless in filterlist
                        break

            if skip:
                logging.debug("           [skip] .... %s " % fil)
                continue

            self.output_filelist.append(fil)
            logging.debug("           [added to output list] .... %s " % fil)


################### End Class STILT ###################


if __name__ == "__main__":
    pass


