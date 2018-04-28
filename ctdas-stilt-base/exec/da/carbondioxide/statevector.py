#!/usr/bin/env python
# ct_statevector_tools.py

"""
Author : peters 

Revision History:
File created on 28 Jul 2010.

"""

import os
import sys
sys.path.append(os.getcwd())

import logging
import numpy as np
from da.baseclasses.statevector import StateVector, EnsembleMember

import da.tools.io4 as io

identifier = 'CarbonTracker Statevector '
version = '0.0'

################### Begin Class CO2StateVector ###################

class CO2StateVector(StateVector):
    """ This is a StateVector object for CarbonTracker. It has a private method to make new ensemble members """

    def get_covariance(self, date, dacycle):
        """ Make a new ensemble from specified matrices, the attribute lag refers to the position in the state vector. 
            Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
            The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]
        """    
        try:
            import matplotlib.pyplot as plt
        except:
            pass

        # Get the needed matrices from the specified covariance files

        file_ocn_cov = dacycle.dasystem['ocn.covariance'] 
        file_bio_cov = dacycle.dasystem['bio.covariance'] 

        # replace YYYY.MM in the ocean covariance file string

        file_ocn_cov = file_ocn_cov.replace('2000.01', date.strftime('%Y.%m'))

        for fil in [file_ocn_cov, file_bio_cov]:
            if not os.path.exists(fil):
                msg = "Cannot find the specified file %s" % fil
                logging.error(msg)
                raise IOError, msg
            else:
                logging.info("Using covariance file: %s" % fil)

        f_ocn = io.ct_read(file_ocn_cov, 'read')
        f_bio = io.ct_read(file_bio_cov, 'read')

        cov_ocn = f_ocn.get_variable('CORMAT')
        if f_bio.variables.has_key('covariance'):
            cov_bio = f_bio.get_variable('covariance')  # newly created CTDAS covariance files
        else:
            cov_bio = f_bio.get_variable('qprior')  # old CarbonTracker covariance files

        f_ocn.close()
        f_bio.close()

        logging.debug("Succesfully closed files after retrieving prior covariance matrices")

        # Once we have the matrices, we can start to make the full covariance matrix, and then decompose it

        fullcov = np.zeros((self.nparams, self.nparams), float)

        nocn = cov_ocn.shape[0]
        nbio = cov_bio.shape[0]

        fullcov[0:nbio, 0:nbio] = cov_bio
        fullcov[nbio:nbio + nocn, nbio:nbio + nocn] = cov_ocn
        fullcov[nocn + nbio, nocn + nbio] = 1.e-10


        try:
            plt.imshow(fullcov)
            plt.colorbar()
            plt.savefig('fullcovariancematrix.png')
            plt.close('all')
            logging.debug("Covariance matrix visualized for inspection")
        except:
            pass

        return fullcov

    def read_from_legacy_file(self, filename, qual='opt'):
        """ 
        :param filename: the full filename for the input NetCDF file
        :param qual: a string indicating whether to read the 'prior' or 'opt'(imized) StateVector from file
        :rtype: None

        Read the StateVector information from a NetCDF file and put in a StateVector object
        In principle the input file will have only one four datasets inside 
        called:
            * `meanstate_prior`, dimensions [nlag, nparamaters]
            * `ensemblestate_prior`, dimensions [nlag,nmembers, nparameters]
            * `meanstate_opt`, dimensions [nlag, nparamaters]
            * `ensemblestate_opt`, dimensions [nlag,nmembers, nparameters]

        This NetCDF information can be written to file using 
        :meth:`~da.baseclasses.statevector.StateVector.write_to_file`

        """
        

        f = io.ct_read(filename, 'read')

        for n in range(self.nlag):
            if qual == 'opt':
                meanstate = f.get_variable('xac_%02d' % (n + 1))
                EnsembleMembers = f.get_variable('adX_%02d' % (n + 1))

            elif qual == 'prior':
                meanstate = f.get_variable('xpc_%02d' % (n + 1))
                EnsembleMembers = f.get_variable('pdX_%02d' % (n + 1))

            if not self.ensemble_members[n] == []:
                self.ensemble_members[n] = []
                logging.warning('Existing ensemble for lag=%d was removed to make place for newly read data' % (n + 1))

            for m in range(self.nmembers):
                newmember = EnsembleMember(m)
                newmember.param_values = EnsembleMembers[m, :].flatten() + meanstate  # add the mean to the deviations to hold the full parameter values
                self.ensemble_members[n].append(newmember)

        f.close()

        logging.info('Successfully read the State Vector from file (%s) ' % filename)
    
    def read_from_file_exceptsam(self, filename, qual='opt'):
        """ 
        :param filename: the full filename for the input NetCDF file
        :param qual: a string indicating whether to read the 'prior' or 'opt'(imized) StateVector from file
        :rtype: None

        Read the StateVector information from a NetCDF file and put in a StateVector object
        In principle the input file will have only one four datasets inside 
        called:
            * `meanstate_prior`, dimensions [nlag, nparamaters]
            * `ensemblestate_prior`, dimensions [nlag,nmembers, nparameters]
            * `meanstate_opt`, dimensions [nlag, nparamaters]
            * `ensemblestate_opt`, dimensions [nlag,nmembers, nparameters]

        This NetCDF information can be written to file using 
        :meth:`~da.baseclasses.statevector.StateVector.write_to_file`

        """

        f = io.ct_read(filename, 'read')
        
        meanstate = f.get_variable('statevectormean_' + qual)
        meanstate[:,39:77] = 1
        ensmembers = f.get_variable('statevectorensemble_' + qual)
        f.close()

        for n in range(self.nlag):
            if not self.ensemble_members[n] == []:
                self.ensemble_members[n] = []
                logging.warning('Existing ensemble for lag=%d was removed to make place for newly read data' % (n + 1))

            for m in range(self.nmembers):
                newmember = EnsembleMember(m)
                newmember.param_values = ensmembers[n, m, :].flatten() + meanstate[n]  # add the mean to the deviations to hold the full parameter values
                self.ensemble_members[n].append(newmember)

        logging.info('Successfully read the State Vector from file (%s) ' % filename)
        logging.info('State Vector set to 1 for South American regions')

################### End Class CO2StateVector ###################


if __name__ == "__main__":
    pass
