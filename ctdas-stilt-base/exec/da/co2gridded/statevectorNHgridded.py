#!/usr/bin/env python
# ct_statevector_tools.py

"""
Author : peters 

Revision History:
File created on 28 Jul 2010.

"""

import os
import sys
import logging
import numpy as np

sys.path.append(os.getcwd())
sys.path.append('../../')

import da.tools.io4 as io
from da.baseclasses.statevector import StateVector, EnsembleMember
identifier = 'CarbonTracker Gridded Statevector '
version = '0.0'

################### Begin Class CtStateVector ###################

class CO2GriddedStateVector(StateVector):
    """ This is a StateVector object for CarbonTracker. It has a private method to make new ensemble members """

    def __init__(self):
        self.ID = identifier
        self.version = version

        # The following code allows the object to be initialized with a dacycle object already present. Otherwise, it can
        # be added at a later moment.

        logging.info('Statevector object initialized: %s' % self.ID)


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

        cov_files = os.listdir(dacycle.dasystem['bio.cov.dir'])
        cov_files = [os.path.join(dacycle.dasystem['bio.cov.dir'], f) for f in cov_files if dacycle.dasystem['bio.cov.prefix'] in f]

        logging.debug("Found %d covariances to use for biosphere" % len(cov_files))

        # replace YYYY.MM in the ocean covariance file string

        file_ocn_cov = file_ocn_cov.replace('2000.01', date.strftime('%Y.%m'))

        cov_files.append(file_ocn_cov)

        covariancematrixlist = []
        for file in cov_files:
            if not os.path.exists(file):
                msg = "Cannot find the specified file %s" % file 
                logging.error(msg)
                raise IOError, msg
            else:
                logging.debug("Using covariance file: %s" % file)

            f = io.ct_read(file, 'read')

            if 'pco2' in file or 'cov_ocean' in file: 
                cov_ocn = f.get_variable('CORMAT')
                parnr = range(9805,9835)
                cov = cov_ocn
            elif 'tropic' in file:
                cov = f.get_variable('covariance')
                parnr = f.get_variable('parameternumber')
            else: 
                cov = f.get_variable('covariance')
                parnr = f.get_variable('parameternumber')
                cov_sf = 90.0 / np.sqrt(cov.diagonal().sum())  # this scaling factor makes the total variance close to the value of a single ecoregion
                cov = cov * cov_sf

            f.close()
            covariancematrixlist.append([cov,parnr,file])

        logging.debug("Succesfully closed files after retrieving prior covariance matrices")

        # Once we have the matrices, we can start to make the full covariance matrix, and then decompose it

        return covariancematrixlist

    def make_new_ensemble(self, lag, covariancematrixlist=[None]):
        """ 
        :param lag: an integer indicating the time step in the lag order
        :param covariancematrix: a list of matrices specifying the covariance distribution to draw from
        :rtype: None
    
        Make a new ensemble, the attribute lag refers to the position in the state vector. 
        Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
        The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]

        The covariance list object to be passed holds a list of matrices with a total number of dimensions [nparams, nparams], which is
        used to draw ensemblemembers from. Each draw is done on a matrix from the list, to make the computational burden smaller when
        the StateVector nparams becomes very large.

        """    
        try:
            import matplotlib.pyplot as plt
        except:
            pass

        if not isinstance(covariancematrixlist, list):
            logging.error("The covariance matrix or matrices must be passed as a list of array objects, exiting..." )
            raise ValueError

        # Check dimensions of covariance matrix list, must add up to nparams

        dims = 1  # start from 1.0 to account for the last parameter that scales Ice+Non-optimized, we have no covariance matrix for this though

        for matrix,parnr,file in covariancematrixlist: 
            dims += matrix.shape[0]

        if dims != self.nparams:
            logging.error("The total dimension of the covariance matrices passed (%d) does not add up to the prescribed nparams (%d), exiting..." % (dims, self.nparams))
            raise ValueError

        # Loop over list if identity matrices and create a matrix of (nparams,nmembers) with the deviations

        dof = 0.0
        dev_matrix = np.zeros((self.nparams, self.nmembers,), 'float')
        randstate = np.random.get_state()

        for matrix,parnr,file in covariancematrixlist:
            # Make a cholesky decomposition of the covariance matrix

            _, s, _ = np.linalg.svd(matrix)
            dof += np.sum(s) ** 2 / sum(s ** 2)
            try:
                C = np.linalg.cholesky(matrix)
            except np.linalg.linalg.LinAlgError, err:
                logging.error('Cholesky decomposition has failed ')
                logging.error('For a matrix of dimensions: %d' % matrix.shape[0])
                logging.debug(err)
                raise np.linalg.linalg.LinAlgError


            # Draw nmembers instances of this distribution

            npoints = matrix.shape[0]

            istart = parnr[0]-1
            istop = parnr[-1]

            for member in range(1, self.nmembers):
                rands = np.random.randn(npoints)
                deviations = np.dot(C, rands)
                dev_matrix[istart:istop, member - 1] = deviations

                dev_matrix[-1, member - 1] = 1.e-10 * np.random.randn()

            logging.debug("Inserted prior covariance matrix from file %s of %d elements into statevector location %d through %d"%(file,len(parnr),istart,istop) )

            #cov2 = np.dot(dev_matrix[istart:istop,:],np.transpose(dev_matrix[istart:istop,:])) / (self.nmembers-1)
            #print matrix.sum(),cov2.sum(),abs(matrix.diagonal()-cov2.diagonal()).max(), matrix.shape,cov2.shape


        logging.debug('Successfully constructed a deviation matrix from covariance structure')
        logging.info('Appr. degrees of freedom in full covariance matrix is %s' % (int(dof)))

        # Now fill the ensemble members with the deviations we have just created


        # Create mean values 

        new_mean = np.ones(self.nparams, float) # standard value for a new time step is 1.0

        # If this is not the start of the filter, average previous two optimized steps into the mix

        if lag == self.nlag - 1 and self.nlag >= 3:
            new_mean += self.ensemble_members[lag - 1][0].param_values + \
                                           self.ensemble_members[lag - 2][0].param_values 
            new_mean = new_mean / 3.0

        # Create the first ensemble member with a deviation of 0.0 and add to list

        new_member = EnsembleMember(0)
        new_member.param_values = new_mean.flatten()  # no deviations
        self.ensemble_members[lag].append(new_member)

        # Create members 1:nmembers and add to ensemble_members list

        for member in range(1, self.nmembers):
            new_member = EnsembleMember(member)
            new_member.param_values = dev_matrix[:, member - 1] + new_mean
            self.ensemble_members[lag].append(new_member)

        logging.debug('%d new ensemble members were added to the state vector # %d' % (self.nmembers, (lag + 1)))





################### End Class CtStateVector ###################


if __name__ == "__main__":
    pass
