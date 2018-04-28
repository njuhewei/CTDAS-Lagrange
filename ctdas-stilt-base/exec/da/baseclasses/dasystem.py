#!/usr/bin/env python
# control.py

"""
.. module:: dasystem
.. moduleauthor:: Wouter Peters 

Revision History:
File created on 26 Aug 2010.

The DaSystem class is found in the module :mod:`dasystem`, or in a specific implementation under the da/ source tree. It is derived from the standard python :class:`dictionary` object. 

It describes the details of the data assimilation system used (i.e., CarbonTracker, or CT Methane, or ....) ::

    datadir         : /Volumes/Storage/CO2/carbontracker/input/ct08/   ! The directory where input data is found
    obs.input.dir   : ${datadir}/obsnc/with_fillvalue                  ! the observation input dir
    obs.input.fname : obs_forecast.nc                                  ! the observation input file
    ocn.covariance  : ${datadir}/oif_p3_era40.dpco2.2000.01.hdf        ! the ocean flux covariance file
    bio.covariance  : ${datadir}/covariance_bio_olson19.nc             ! the biosphere flux covariance file
    deltaco2.prefix : oif_p3_era40.dpco2                               ! the type of ocean product used
    regtype         : olson19_oif30                                    ! the ecoregion definitions
    nparameters     : 240                                              ! the number of parameters to solve for
    random.seed     : 4385                                             ! the random seed for the first cycle
    regionsfile     : transcom_olson19_oif30.hdf                       ! the ecoregion defintion mask file

    ! Info on the sites file used

    obs.sites.rc        : ${datadir}/sites_and_weights_co2.ct10.rc     ! the weights in the covariance matric of each obs

The full baseclass description:

.. autoclass:: da.baseclasses.dasystem.DaSystem
   :members:

"""

import logging
import da.tools.rc as rc 
################### Begin Class DaSystem ###################

class DaSystem(dict):
    """ 
    Information on the data assimilation system used. This is normally an rc-file with settings.
    """

    def __init__(self, rcfilename):
        """
        Initialization occurs from passed rc-file name, items in the rc-file will be added
        to the dictionary
        """

        self.ID = 'CarbonTracker CO2'    # the identifier gives the platform name
        self.load_rc(rcfilename)

        logging.debug("Data Assimilation System initialized: %s" % self.ID)

    def load_rc(self, rcfilename):
        """ 
        This method loads a DA System Info rc-file with settings for this simulation 
        """
        for k, v in rc.read(rcfilename).iteritems():
            self[k] = v
        
        logging.debug("DA System Info rc-file (%s) loaded successfully" % rcfilename)


    def validate(self):
        """ 
        validate the contents of the rc-file given a dictionary of required keys
        """
        needed_rc_items = {}

        for k, v in self.iteritems():
            if v == 'True' : 
                self[k] = True
            if v == 'False': 
                self[k] = False

        for key in needed_rc_items:
            if not self.has_key(key):
                msg = 'Missing a required value in rc-file : %s' % key
                logging.error(msg)
                raise IOError, msg
        logging.debug('DA System Info settings have been validated succesfully')

################### End Class DaSystem ###################


if __name__ == "__main__":
    pass
