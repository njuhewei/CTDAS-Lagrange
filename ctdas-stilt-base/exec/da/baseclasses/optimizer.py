#!/usr/bin/env python
# optimizer.py

"""
.. module:: optimizer
.. moduleauthor:: Wouter Peters 

Revision History:
File created on 28 Jul 2010.

"""

import logging
import numpy as np
import numpy.linalg as la
import da.tools.io4 as io

identifier = 'Optimizer baseclass'
version = '0.0'

################### Begin Class Optimizer ###################

class Optimizer(object):
    """
        This creates an instance of an optimization object. It handles the minimum least squares optimization
        of the state vector given a set of sample objects. Two routines will be implemented: one where the optimization
        is sequential and one where it is the equivalent matrix solution. The choice can be made based on considerations of speed
        and efficiency.
    """

    def __init__(self):
        self.ID = identifier
        self.version = version

        logging.info('Optimizer object initialized: %s' % self.ID)

    def setup(self, dims):
        self.nlag = dims[0]
        self.nmembers = dims[1]
        self.nparams = dims[2]
        self.nobs = dims[3]
        self.create_matrices()

    def create_matrices(self):
        """ Create Matrix space needed in optimization routine """

        # mean state  [X]
        self.x = np.zeros((self.nlag * self.nparams,), float)
        # deviations from mean state  [X']
        self.X_prime = np.zeros((self.nlag * self.nparams, self.nmembers,), float)
        # mean state, transported to observation space [ H(X) ]
        self.Hx = np.zeros((self.nobs,), float)
        # deviations from mean state, transported to observation space [ H(X') ]
        self.HX_prime = np.zeros((self.nobs, self.nmembers), float)
        # observations
        self.obs = np.zeros((self.nobs,), float)
        # observation ids
        self.obs_ids = np.zeros((self.nobs,), float)
        # covariance of observations
        # Total covariance of fluxes and obs in units of obs [H P H^t + R]
        if self.algorithm == 'Serial':
            self.R = np.zeros((self.nobs,), float)
            self.HPHR = np.zeros((self.nobs,), float)
        else:
            self.R = np.zeros((self.nobs, self.nobs,), float)
            self.HPHR = np.zeros((self.nobs, self.nobs,), float)
        # localization of obs
        self.may_localize = np.zeros(self.nobs, bool)
        # rejection of obs
        self.may_reject = np.zeros(self.nobs, bool)
        # flags of obs
        self.flags = np.zeros(self.nobs, int)
        # species type
        self.species = np.zeros(self.nobs, str)
        # species type
        self.sitecode = np.zeros(self.nobs, str)

        # species mask
        self.speciesmask = {}

        # Kalman Gain matrix
        self.KG = np.zeros((self.nlag * self.nparams, self.nobs,), float)
        #self.KG = np.zeros((self.nlag * self.nparams,), float)

    def state_to_matrix(self, statevector):
        allsites = []      # collect all obs for n=1,..,nlag
        allobs = []      # collect all obs for n=1,..,nlag
        allmdm = []      # collect all mdm for n=1,..,nlag
        allids = []  # collect all model samples for n=1,..,nlag
        allreject = []  # collect all model samples for n=1,..,nlag
        alllocalize = []  # collect all model samples for n=1,..,nlag
        allflags = []  # collect all model samples for n=1,..,nlag
        allspecies = []  # collect all model samples for n=1,..,nlag
        allsimulated = []  # collect all members model samples for n=1,..,nlag

        for n in range(self.nlag):
            samples = statevector.obs_to_assimilate[n]
            members = statevector.ensemble_members[n]
            self.x[n * self.nparams:(n + 1) * self.nparams] = members[0].param_values
            self.X_prime[n * self.nparams:(n + 1) * self.nparams, :] = np.transpose(np.array([m.param_values for m in members]))

            if samples != None:        
                self.rejection_threshold = samples.rejection_threshold

                allreject.extend(samples.getvalues('may_reject'))
                alllocalize.extend(samples.getvalues('may_localize'))
                allflags.extend(samples.getvalues('flag'))
                allspecies.extend(samples.getvalues('species'))
                allobs.extend(samples.getvalues('obs'))
                allsites.extend(samples.getvalues('code'))
                allmdm.extend(samples.getvalues('mdm'))
                allids.extend(samples.getvalues('id'))

                simulatedensemble = samples.getvalues('simulated')
		for s in range(simulatedensemble.shape[0]):
                	allsimulated.append(simulatedensemble[s])

        self.obs[:] = np.array(allobs)
        self.obs_ids[:] = np.array(allids)
        self.HX_prime[:, :] = np.array(allsimulated)
        self.Hx[:] = self.HX_prime[:, 0]

        self.may_reject[:] = np.array(allreject)
        self.may_localize[:] = np.array(alllocalize)
        self.flags[:] = np.array(allflags)
        self.species[:] = np.array(allspecies)
        self.sitecode = allsites

        self.X_prime = self.X_prime - self.x[:, np.newaxis] # make into a deviation matrix
        self.HX_prime = self.HX_prime - self.Hx[:, np.newaxis] # make a deviation matrix

        if self.algorithm == 'Serial':
            for i, mdm in enumerate(allmdm):
                self.R[i] = mdm ** 2
        else:
            for i, mdm in enumerate(allmdm):
                self.R[i, i] = mdm ** 2
                    
    def matrix_to_state(self, statevector):
        for n in range(self.nlag):
            members = statevector.ensemble_members[n]
            for m, mem in enumerate(members):
                members[m].param_values[:] = self.X_prime[n * self.nparams:(n + 1) * self.nparams, m] + self.x[n * self.nparams:(n + 1) * self.nparams]     

        logging.debug('Returning optimized data to the StateVector, setting "StateVector.isOptimized = True" ')

    def write_diagnostics(self, filename, type):
        """
            Open a NetCDF file and write diagnostic output from optimization process:

                - calculated residuals
                - model-data mismatches
                - HPH^T
                - prior ensemble of samples
                - posterior ensemble of samples
                - prior ensemble of fluxes
                - posterior ensemble of fluxes

            The type designation refers to the writing of prior or posterior data and is used in naming the variables"
        """

        # Open or create file

        if type == 'prior':
            f = io.CT_CDF(filename, method='create')
            logging.debug('Creating new diagnostics file for optimizer (%s)' % filename)
        elif type == 'optimized':
            f = io.CT_CDF(filename, method='write')
            logging.debug('Opening existing diagnostics file for optimizer (%s)' % filename)

        # Add dimensions 

        dimparams = f.add_params_dim(self.nparams)
        dimmembers = f.add_members_dim(self.nmembers)
        dimlag = f.add_lag_dim(self.nlag, unlimited=False)
        dimobs = f.add_obs_dim(self.nobs)
        dimstate = f.add_dim('nstate', self.nparams * self.nlag)
        dim200char = f.add_dim('string_of200chars', 200)

        # Add data, first the ones that are written both before and after the optimization

        savedict = io.std_savedict.copy() 
        savedict['name'] = "statevectormean_%s" % type
        savedict['long_name'] = "full_statevector_mean_%s" % type
        savedict['units'] = "unitless"
        savedict['dims'] = dimstate
        savedict['values'] = self.x.tolist()
        savedict['comment'] = 'Full %s state vector mean ' % type
        f.add_data(savedict)

        savedict = io.std_savedict.copy()
        savedict['name'] = "statevectordeviations_%s" % type
        savedict['long_name'] = "full_statevector_deviations_%s" % type
        savedict['units'] = "unitless"
        savedict['dims'] = dimstate + dimmembers
        savedict['values'] = self.X_prime.tolist()
        savedict['comment'] = 'Full state vector %s deviations as resulting from the optimizer' % type
        f.add_data(savedict)

        savedict = io.std_savedict.copy()
        savedict['name'] = "modelsamplesmean_%s" % type
        savedict['long_name'] = "modelsamplesforecastmean_%s" % type
        savedict['units'] = "mol mol-1"
        savedict['dims'] = dimobs
        savedict['values'] = self.Hx.tolist()
        savedict['comment'] = '%s mean mole fractions based on %s state vector' % (type, type)
        f.add_data(savedict)

        savedict = io.std_savedict.copy()
        savedict['name'] = "modelsamplesdeviations_%s" % type
        savedict['long_name'] = "modelsamplesforecastdeviations_%s" % type
        savedict['units'] = "mol mol-1"
        savedict['dims'] = dimobs + dimmembers
        savedict['values'] = self.HX_prime.tolist()
        savedict['comment'] = '%s mole fraction deviations based on %s state vector' % (type, type)
        f.add_data(savedict)

        # Continue with prior only data

        if type == 'prior':

            savedict = io.std_savedict.copy()
            savedict['name'] = "sitecode"
            savedict['long_name'] = "site code propagated from observation file"
            savedict['dtype'] = "char"
            savedict['dims'] = dimobs + dim200char
            savedict['values'] = self.sitecode
            savedict['missing_value'] = '!'
            f.add_data(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "observed"
            savedict['long_name'] = "observedvalues"
            savedict['units'] = "mol mol-1"
            savedict['dims'] = dimobs
            savedict['values'] = self.obs.tolist()
            savedict['comment'] = 'Observations used in optimization'
            f.add_data(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "obspack_num"
            savedict['dtype'] = "int64"
            savedict['long_name'] = "Unique_ObsPack_observation_number"
            savedict['units'] = ""
            savedict['dims'] = dimobs
            savedict['values'] = self.obs_ids.tolist()
            savedict['comment'] = 'Unique observation number across the entire ObsPack distribution'
            f.add_data(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "modeldatamismatchvariance"
            savedict['long_name'] = "modeldatamismatch variance"
            savedict['units'] = "[mol mol-1]^2"
            if self.algorithm == 'Serial':
                savedict['dims'] = dimobs
            else: savedict['dims'] = dimobs + dimobs
            savedict['values'] = self.R.tolist()
            savedict['comment'] = 'Variance of mole fractions resulting from model-data mismatch'
            f.add_data(savedict)

        # Continue with posterior only data

        elif type == 'optimized':
            
            savedict = io.std_savedict.copy()
            savedict['name'] = "totalmolefractionvariance"
            savedict['long_name'] = "totalmolefractionvariance"
            savedict['units'] = "[mol mol-1]^2"
            if self.algorithm == 'Serial':
                savedict['dims'] = dimobs
            else: savedict['dims'] = dimobs + dimobs
            savedict['values'] = self.HPHR.tolist()
            savedict['comment'] = 'Variance of mole fractions resulting from prior state and model-data mismatch'
            f.add_data(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "flag"
            savedict['long_name'] = "flag_for_obs_model"
            savedict['units'] = "None"
            savedict['dims'] = dimobs
            savedict['values'] = self.flags.tolist()
            savedict['comment'] = 'Flag (0/1/2/99) for observation value, 0 means okay, 1 means QC error, 2 means rejected, 99 means not sampled'
            f.add_data(savedict)

            savedict = io.std_savedict.copy()
            savedict['name'] = "kalmangainmatrix"
            savedict['long_name'] = "kalmangainmatrix"
            savedict['units'] = "unitless molefraction-1"
            savedict['dims'] = dimstate + dimobs
            savedict['values'] = self.KG.tolist()
            savedict['comment'] = 'Kalman gain matrix of all obs and state vector elements'
            f.add_data(savedict)

        f.close()
        logging.debug('Diagnostics file closed')


    def serial_minimum_least_squares(self):
        """ Make minimum least squares solution by looping over obs"""
        for n in range(self.nobs):

            # Screen for flagged observations (for instance site not found, or no sample written from model)

            if self.flags[n] != 0:
                logging.debug('Skipping observation (%s,%i) because of flag value %d' % (self.sitecode[n], self.obs_ids[n], self.flags[n]))
                continue

            # Screen for outliers greather than 3x model-data mismatch, only apply if obs may be rejected

            res = self.obs[n] - self.Hx[n]

            if self.may_reject[n]:
                threshold = self.rejection_threshold * np.sqrt(self.R[n])
                if np.abs(res) > threshold:
                    logging.debug('Rejecting observation (%s,%i) because residual (%f) exceeds threshold (%f)' % (self.sitecode[n], self.obs_ids[n], res, threshold))
                    self.flags[n] = 2
                    continue

            logging.debug('Proceeding to assimilate observation %s, %i' % (self.sitecode[n], self.obs_ids[n]))

            PHt = 1. / (self.nmembers - 1) * np.dot(self.X_prime, self.HX_prime[n, :])
            self.HPHR[n] = 1. / (self.nmembers - 1) * (self.HX_prime[n, :] * self.HX_prime[n, :]).sum() + self.R[n]

            self.KG[:] = PHt / self.HPHR[n]

            if self.may_localize[n]:
                logging.debug('Trying to localize observation %s, %i' % (self.sitecode[n], self.obs_ids[n]))
                self.localize(n)
            else:
                logging.debug('Not allowed to localize observation %s, %i' % (self.sitecode[n], self.obs_ids[n]))

            alpha = np.double(1.0) / (np.double(1.0) + np.sqrt((self.R[n]) / self.HPHR[n]))

            self.x[:] = self.x + self.KG[:] * res

            for r in range(self.nmembers):
                self.X_prime[:, r] = self.X_prime[:, r] - alpha * self.KG[:] * (self.HX_prime[n, r])

#WP !!!! Very important to first do all obervations from n=1 through the end, and only then update 1,...,n. The current observation
#WP      should always be updated last because it features in the loop of the adjustments !!!!

            for m in range(n + 1, self.nobs):
                res = self.obs[n] - self.Hx[n]
                fac = 1.0 / (self.nmembers - 1) * (self.HX_prime[n, :] * self.HX_prime[m, :]).sum() / self.HPHR[n]
                self.Hx[m] = self.Hx[m] + fac * res
                self.HX_prime[m, :] = self.HX_prime[m, :] - alpha * fac * self.HX_prime[n, :]

            for m in range(1, n + 1):
                res = self.obs[n] - self.Hx[n]
                fac = 1.0 / (self.nmembers - 1) * (self.HX_prime[n, :] * self.HX_prime[m, :]).sum() / self.HPHR[n]
                self.Hx[m] = self.Hx[m] + fac * res
                self.HX_prime[m, :] = self.HX_prime[m, :] - alpha * fac * self.HX_prime[n, :]

            
    def bulk_minimum_least_squares(self):
        """ Make minimum least squares solution by solving matrix equations"""
        

        # Create full solution, first calculate the mean of the posterior analysis

        HPH = np.dot(self.HX_prime, np.transpose(self.HX_prime)) / (self.nmembers - 1)   # HPH = 1/N * HX' * (HX')^T
        self.HPHR[:, :] = HPH + self.R                                                            # HPHR = HPH + R
        HPb = np.dot(self.X_prime, np.transpose(self.HX_prime)) / (self.nmembers - 1)    # HP = 1/N X' * (HX')^T
        self.KG[:, :] = np.dot(HPb, la.inv(self.HPHR))                                         # K = HP/(HPH+R)

        for n in range(self.nobs):
            self.localize(n)

        self.x[:] = self.x + np.dot(self.KG, self.obs - self.Hx)                             # xa = xp + K (y-Hx)

        # And next make the updated ensemble deviations. Note that we calculate P by using the full equation (10) at once, and 
        # not in a serial update fashion as described in Whitaker and Hamill. 
        # For the current problem with limited N_obs this is easier, or at least more straightforward to do.

        I = np.identity(self.nlag * self.nparams)
        sHPHR = la.cholesky(self.HPHR)                                  # square root of HPH+R
        part1 = np.dot(HPb, np.transpose(la.inv(sHPHR)))                 # HP(sqrt(HPH+R))^-1
        part2 = la.inv(sHPHR + np.sqrt(self.R))                           # (sqrt(HPH+R)+sqrt(R))^-1
        Kw = np.dot(part1, part2)                                     # K~
        self.X_prime[:, :] = np.dot(I, self.X_prime) - np.dot(Kw, self.HX_prime)         # HX' = I - K~ * HX'


        # Now do the adjustments of the modeled mole fractions using the linearized ensemble. These are not strictly needed but can be used
        # for diagnosis.

        part3 = np.dot(HPH, np.transpose(la.inv(sHPHR)))                           # HPH(sqrt(HPH+R))^-1
        Kw = np.dot(part3, part2)                                               # K~
        self.Hx[:] = self.Hx + np.dot(np.dot(HPH, la.inv(self.HPHR)), self.obs - self.Hx)  # Hx  = Hx+ HPH/HPH+R (y-Hx)
        self.HX_prime[:, :] = self.HX_prime - np.dot(Kw, self.HX_prime)                            # HX' = HX'- K~ * HX'

        logging.info('Minimum Least Squares solution was calculated, returning')

    def set_localization(self):
        """ determine which localization to use """
        self.localization = True
        self.localizetype = "None"
        logging.info("Current localization option is set to %s" % self.localizetype)

    def localize(self, n):
        """ localize the Kalman Gain matrix """
        logging.debug('Not localized observation %s, %i' % (self.sitecode[n], self.obs_ids[n]))
    
    def set_algorithm(self):
        self.algorithm = 'Serial'
        logging.info("Current algorithm is set to %s in baseclass optimizer" % self.algorithm)

################### End Class Optimizer ###################



if __name__ == "__main__":
    pass
