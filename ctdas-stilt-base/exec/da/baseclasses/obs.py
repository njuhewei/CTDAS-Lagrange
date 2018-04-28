#!/usr/bin/env python
# obs.py

"""
.. module:: obs
.. moduleauthor:: Wouter Peters 

Revision History:
File created on 28 Jul 2010.

.. autoclass:: da.baseclasses.obs.Observations 
   :members: setup, Validate, add_observations, add_simulations, add_model_data_mismatch, write_sample_coords  

.. autoclass:: da.baseclasses.obs.ObservationList 
   :members: __init__

"""

import logging
from numpy import array, ndarray

identifier = 'Observations baseclass'
version = '0.0'

################### Begin Class Observations ###################

class Observations(object):
    """ 
    The baseclass Observations is a generic object that provides a number of methods required for any type of observations used in 
    a data assimilation system. These methods are called from the CarbonTracker pipeline. 

    .. note:: Most of the actual functionality will need to be provided through a derived Observations class with the methods 
              below overwritten. Writing your own derived class for Observations is one of the first tasks you'll likely 
              perform when extending or modifying the CarbonTracker Data Assimilation Shell.

    Upon initialization of the class, an object is created that holds no actual data, but has a placeholder attribute `self.Data` 
    which is an empty list of type :class:`~da.baseclasses.obs.ObservationList`. An ObservationList object is created when the 
    method :meth:`~da.baseclasses.obs.Observations.add_observations` is invoked in the pipeline. 

    From the list of observations, a file is written  by method 
    :meth:`~da.baseclasses.obs.Observations.write_sample_info`
    with the sample info needed by the 
    :class:`~da.baseclasses.observationoperator.ObservationOperator` object. The values returned after sampling 
    are finally added by :meth:`~da.baseclasses.obs.Observations.add_simulations`

    """ 

    def __init__(self):
        """
        create an object with an identifier, version, and an empty ObservationList
        """
        self.ID = identifier
        self.version = version
        self.datalist = []  # initialize with an empty list of obs

        # The following code allows the object to be initialized with a dacycle object already present. Otherwise, it can
        # be added at a later moment.

        logging.info('Observations object initialized: %s' % self.ID)

    def getlength(self):
        return len(self.datalist)

    def setup(self, cycleparams):
        """ Perform all steps needed to start working with observational data, this can include moving data, concatenating files,
            selecting datasets, etc.
        """

    def add_observations(self):
        """ 
        Add actual observation data to the Observations object. This is in a form of an 
        :class:`~da.baseclasses.obs.ObservationList` that is contained in self.Data. The 
        list has as only requirement that it can return the observed+simulated values 
        through the method :meth:`~da.baseclasses.obs.ObservationList.getvalues`

        """

    def add_simulations(self):
        """ Add the simulation data to the Observations object. 
        """

    def add_model_data_mismatch(self):
        """ 
            Get the model-data mismatch values for this cycle.
        """

    def write_sample_coords(self,obsinputfile):
        """ 
            Write the information needed by the observation operator to a file. Return the filename that was written for later use
        """

    def write_sample_auxiliary(self, auxoutputfile):
        """ 
            Write selected additional information contained in the Observations object to a file for later processing. 

        """

    def getvalues(self, name, constructor=array):

        result = constructor([getattr(o, name) for o in self.datalist])
        if isinstance(result, ndarray): 
            return result.squeeze()
        else:
            return result


################### End Class Observations ###################


