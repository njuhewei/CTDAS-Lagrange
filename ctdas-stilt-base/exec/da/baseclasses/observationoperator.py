#!/usr/bin/env python
# model.py

"""
.. module:: observationoperator
.. moduleauthor:: Wouter Peters 

Revision History:
File created on 30 Aug 2010.

"""

import logging

identifier = 'RandomizerObservationOperator'
version = '1.0'

################### Begin Class ObservationOperator ###################
class ObservationOperator(object):
    """
    Testing
    =======
    This is a class that defines an ObervationOperator. This object is used to control the sampling of
    a statevector in the ensemble Kalman filter framework. The methods of this class specify which (external) code
    is called to perform the sampling, and which files should be read for input and are written for output.

    The baseclasses consist mainly of empty methods that require an application specific application. The baseclass will take observed values, and perturb them with a random number chosen from the model-data mismatch distribution. This means no real operator will be at work, but random normally distributed residuals will come out of y-H(x) and thus the inverse model can proceed. This is mainly for testing the code...

    """

    def __init__(self, dacycle=None):
        """ The instance of an ObservationOperator is application dependent """
        self.ID = identifier
        self.version = version
        self.restart_filelist = []
        self.output_filelist = []
        self.outputdir = None # Needed for opening the samples.nc files created 

        logging.info('Observation Operator object initialized: %s' % self.ID)

        # The following code allows the object to be initialized with a dacycle object already present. Otherwise, it can
        # be added at a later moment.

        if dacycle != None:
            self.dacycle = dacycle
        else:
            self.dacycle = {}

    
    def get_initial_data(self):
        """ This method places all initial data needed by an ObservationOperator in the proper folder for the model """

    def setup(self,dacycle):
        """ Perform all steps necessary to start the observation operator through a simple Run() call """

        self.dacycle = dacycle
	self.outputdir = dacycle['dir.output']

    def prepare_run(self):
        """ Prepare the running of the actual forecast model, for example compile code """

	import os

	# Define the name of the file that will contain the modeled output of each observation

    	self.simulated_file = os.path.join(self.outputdir, 'samples_simulated.%s.nc' % self.dacycle['time.sample.stamp'])
    	self.forecast_nmembers = int(self.dacycle['da.optimizer.nmembers'])

    def validate_input(self):
        """ Make sure that data needed for the ObservationOperator (such as observation input lists, or parameter files)
            are present.
        """
    def save_data(self):
        """ Write the data that is needed for a restart or recovery of the Observation Operator to the save directory """

    def run(self):
	"""
	 This Randomizer will take the original observation data in the Obs object, and simply copy each mean value. Next, the mean 
	 value will be perturbed by a random normal number drawn from a specified uncertainty of +/- 2 ppm
	"""

	import da.tools.io4 as io
	import numpy as np

	# Create a flask output file in TM5-style (to be updated later?) to hold simulated values for later reading

    	f = io.CT_CDF(self.simulated_file, method='create')
    	logging.debug('Creating new simulated observation file in ObservationOperator (%s)' % self.simulated_file)
	
        dimid = f.createDimension('obs_num', size=None)
	dimid = ('obs_num',)
        savedict = io.std_savedict.copy() 
        savedict['name'] = "obs_num"
        savedict['dtype'] = "int"
        savedict['long_name'] = "Unique_Dataset_observation_index_number"
        savedict['units'] = ""
        savedict['dims'] = dimid
        savedict['comment'] = "Unique index number within this dataset ranging from 0 to UNLIMITED."
        f.add_data(savedict,nsets=0)

        dimmember = f.createDimension('nmembers', size=self.forecast_nmembers)
	dimmember = ('nmembers',)
        savedict = io.std_savedict.copy() 
        savedict['name'] = "flask"
        savedict['dtype'] = "float"
        savedict['long_name'] = "mole_fraction_of_trace_gas_in_air"
        savedict['units'] = "mol tracer (mol air)^-1"
        savedict['dims'] = dimid + dimmember
        savedict['comment'] = "Simulated model value created by RandomizerObservationOperator"
        f.add_data(savedict,nsets=0)

	# Open file with x,y,z,t of model samples that need to be sampled

        f_in = io.ct_read(self.dacycle['ObsOperator.inputfile'],method='read') 

	# Get simulated values and ID

        ids = f_in.get_variable('obs_num')
        obs = f_in.get_variable('observed')
        mdm = f_in.get_variable('modeldatamismatch')

	# Loop over observations, add random white noise, and write to file

	for i,data in enumerate(zip(ids,obs,mdm)):
	    f.variables['obs_num'][i] = data[0]		
	    f.variables['flask'][i,:] = data[1]+np.random.randn(self.forecast_nmembers)*data[2]

	f.close()
	f_in.close()

	# Report success and exit

    	logging.info('ObservationOperator finished successfully, output file written (%s)' % self.simulated_file)

    def run_forecast_model(self):
        self.prepare_run()
        self.validate_input()
        self.run()
        self.save_data()


################### End Class ObservationOperator ###################

class RandomizerObservationOperator(ObservationOperator):
    """ This class holds methods and variables that are needed to use a random number generated as substitute
        for a true observation operator. It takes observations and returns values for each obs, with a specified 
        amount of white noise added 
    """



if __name__ == "__main__":
    pass
