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

        #file_ocn_cov = dacycle.dasystem['ocn.covariance']

        cov_files = os.listdir(dacycle.dasystem['bio.cov.dir'])
        cov_files = [os.path.join(dacycle.dasystem['bio.cov.dir'], f) for f in cov_files if dacycle.dasystem['bio.cov.prefix'] in f]

        logging.debug("Found %d covariances to use for biosphere" % len(cov_files))

        # replace YYYY.MM in the ocean covariance file string

        #file_ocn_cov = file_ocn_cov.replace('2000.01', date.strftime('%Y.%m'))

        #cov_files.append(file_ocn_cov)

        covariancematrixlist = []
        for file in cov_files:
            if not os.path.exists(file):
                msg = "Cannot find the specified file %s" % file
                logging.error(msg)
                raise IOError, msg
            else:
                logging.debug("Using covariance file: %s" % file)

            f = io.ct_read(file, 'read')

            if 'pco2' in file:
                cov_ocn = f.get_variable('CORMAT')
                cov = cov_ocn
            else:
                cov = f.get_variable('covariance')
                #cov_sf      = 10.0/np.sqrt(cov.diagonal().sum())  # this scaling factor makes the total variance close to the value of a single ecoregion
                cov_sf = 360. / np.sqrt(cov.diagonal().sum())  # this scaling factor makes the total variance close to the value of a single ecoregion #I use 360 to boost up the P matrix uncertainty
                cov1 = cov * cov_sf * (1.e-6)**2  # here you assume that your P matrix has units of mol m-2 s-1 squared.

            f.close()
            covariancematrixlist.append(cov1)

        # Boundary conditions covariance

        cov = np.array([[2*2]])
        covariancematrixlist.append(cov)
        covariancematrixlist.append(cov)
        covariancematrixlist.append(cov)
        covariancematrixlist.append(cov)

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

        #dims = 1  # start from 1.0 to account for the last parameter that scales Ice+Non-optimized, we have no covariance matrix for this though

        dims = 0

        for matrix in covariancematrixlist:
            dims += matrix.shape[0]

        if dims != self.nparams:
            logging.error("The total dimension of the covariance matrices passed (%d) does not add up to the prescribed nparams (%d), exiting..." % (dims, self.nparams))
            raise ValueError

        # Loop over list if identity matrices and create a matrix of (nparams,nmembers) with the deviations

        istart = 0
        istop = 0
        dof = 0.0
        dev_matrix = np.zeros((self.nparams, self.nmembers,), 'float')
        randstate = np.random.get_state()

        for i,matrix in enumerate(covariancematrixlist):
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

            istop = istop + npoints

            #if i < len(covariancematrixlist)-1:
            if i < len(covariancematrixlist)-4:
                for member in range(1, self.nmembers):
                    rands = np.random.randn(npoints)
                    deviations = np.dot(C, rands)
                    dev_matrix[istart:istop, member - 1] = deviations
                    #dev_matrix[istop, member - 1] = 1.e-10 * np.random.randn()



            #if i == len(covariancematrixlist)-1:
            if i >= len(covariancematrixlist)-4 and i <= len(covariancematrixlist)-1:
                for member in range(1, self.nmembers):
                    rands = np.random.randn(npoints)
                    deviations = np.dot(C, rands)
                    #logging.info('deviation boundary conditions %s,%s,%s'%(deviations,istart,istop))
                    dev_matrix[istart:istop, member - 1] = deviations

            istart = istart + npoints


        #for i in range(self.nparams):
            #logging.info('i=%s,deviation=%s'%(i,dev_matrix[i,:]))


        logging.debug('Successfully constructed a deviation matrix from covariance structure')
        logging.info('Appr. degrees of freedom in full covariance matrix is %s' % (int(dof)))

        # Now fill the ensemble members with the deviations we have just created


        # Create mean values

        #new_mean = np.ones(self.nparams, float) # standard value for a new time step is 1.0
        new_mean = np.zeros(self.nparams, float)
        new_mean[-4] = np.zeros((1),float) # standard value for boundary conditions for new time step is 0.0
        new_mean[-3] = np.zeros((1),float)
        new_mean[-2] = np.zeros((1),float)
        new_mean[-1] = np.zeros((1),float)

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


    def write_members_to_file(self, lag, outdir,endswith='.nc'):
        """
           :param: lag: Which lag step of the filter to write, must lie in range [1,...,nlag]
           :param: outdir: Directory where to write files
           :param: endswith: Optional label to add to the filename, default is simply .nc
           :rtype: None

           Write ensemble member information to a NetCDF file for later use. The standard output filename is
           *parameters.DDD.nc* where *DDD* is the number of the ensemble member. Standard output file location
           is the `dir.input` of the dacycle object. In principle the output file will have only two datasets inside
           called `parametervalues` which is of dimensions `nparameters` and `parametermap` which is of dimensions (180,360).
           This dataset can be read and used by a :class:`~da.baseclasses.observationoperator.ObservationOperator` object.

           .. note:: if more, or other information is needed to complete the sampling of the ObservationOperator you
                     can simply inherit from the StateVector baseclass and overwrite this write_members_to_file function.

        """

        # These import statements caused a crash in netCDF4 on MacOSX. No problems on Jet though. Solution was
        # to do the import already at the start of the module, not just in this method.

        #import da.tools.io as io
        #import da.tools.io4 as io

        members = self.ensemble_members[lag]

        for mem in members:
            filename = os.path.join(outdir, 'parameters.%03d%s' % (mem.membernumber, endswith))
            ncf = io.CT_CDF(filename, method='create')
            #dimparams = ncf.add_params_dim(self.nparams-1)
            dimparams = ncf.add_params_dim(self.nparams-4)
            #dimparams_bc  = ncf.add_dim('nparameters_bc',1)
            dimparams_bc  = ncf.add_dim('nparameters_bc',4)
            dimgrid = ncf.add_latlon_dim()

            data = mem.param_values

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametervalues"
            savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimparams
            #savedict['values'] = data[:self.nparams-1]
            savedict['values'] = data[:self.nparams-4]
            savedict['comment'] = 'These are parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            griddata = self.vector2grid(vectordata=data)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap"
            savedict['long_name'] = "parametermap_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid
            savedict['values'] = griddata.tolist()
            savedict['comment'] = 'These are gridded parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametervalues_bc"
            savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimparams_bc
            #savedict['values'] = data[-1]
            savedict['values'] = data[-4:]   #data[-4:-1]
            savedict['comment'] = 'These are parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)



            ncf.close()

            logging.debug('Successfully wrote data from ensemble member %d to file (%s) ' % (mem.membernumber, filename))



################### End Class CtStateVector ###################


if __name__ == "__main__":
    pass
