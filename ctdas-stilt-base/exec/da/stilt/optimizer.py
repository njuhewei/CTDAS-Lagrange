#!/usr/bin/env python
# optimizer.py

"""
Author : peters

Revision History:
File created on 28 Jul 2010.

"""

import os
import sys
import logging
sys.path.append(os.getcwd())
from da.baseclasses.optimizer import Optimizer


identifier = 'Ensemble Square Root Filter'
version = '0.0'

################### Begin Class CO2Optimizer ###################

class CO2Optimizer(Optimizer):
    """
        This creates an instance of a CarbonTracker optimization object. The base class it derives from is the optimizer object.
        Additionally, this CO2Optimizer implements a special localization option following the CT2007 method.

        All other methods are inherited from the base class Optimizer.
    """

    def set_localization(self, loctype='None'):
        """ determine which localization to use """

        if loctype == 'CT2007':
            self.localization = True
            self.localizetype = 'CT2007'
            #T-test values for two-tailed student's T-test using 95% confidence interval for some options of nmembers
            if self.nmembers == 50:
                self.tvalue = 2.0086
            elif self.nmembers == 100:
                self.tvalue = 1.9840
            elif self.nmembers == 150:
                self.tvalue = 1.97591
            elif self.nmembers == 200:
                self.tvalue = 1.9719
            else: self.tvalue = 0
        else:
            self.localization = False
            self.localizetype = 'None'

        logging.info("Current localization option is set to %s" % self.localizetype)
        if self.localization == True:
            if self.tvalue == 0:
                logging.error("Critical tvalue for localization not set for %i ensemble members"%(self.nmembers))
                sys.exit(2)
            else: logging.info("Used critical tvalue %0.05f is based on 95%% probability and %i ensemble members in a two-tailed student's T-test"%(self.tvalue,self.nmembers))

    def localize(self, n):
        """ localize the Kalman Gain matrix """
        import numpy as np

        if not self.localization:
            logging.debug('Not localized observation %i' % self.obs_ids[n])
            return
        if self.localizetype == 'CT2007':
            count_localized = 0
            for r in range(self.nlag * self.nparams):
                corr = np.corrcoef(self.HX_prime[n, :], self.X_prime[r, :].squeeze())[0, 1]
                prob = corr / np.sqrt((1.000000001 - corr ** 2) / (self.nmembers - 2))
                if abs(prob) < self.tvalue:
                    self.KG[r] = 0.0
                    count_localized = count_localized + 1
            logging.debug('Localized observation %i, %i%% of values set to 0' % (self.obs_ids[n],count_localized*100/(self.nlag * self.nparams)))

    def set_algorithm(self, algorithm='Serial'):
        """ determine which minimum least squares algorithm to use """

        if algorithm == 'Serial':
            self.algorithm = 'Serial'
        else:
            self.algorithm = 'Bulk'

        logging.info("Current minimum least squares algorithm is set to %s" % self.algorithm)




    def serial_minimum_least_squares(self):
        """ Make minimum least squares solution by looping over obs"""
        import numpy as np
        #for n in range(self.nobs):

         #   res = self.obs[n] - self.Hx[n]

          #  if self.may_reject[n]:
           #     threshold = self.rejection_threshold * np.sqrt(self.R[n])
            #    if np.abs(res) > threshold:
             #       logging.debug('Rejecting observation (%s,%i) because residual (%f) exceeds threshold (%f)' % (self.sitecode[n], self.obs_ids[n], res, threshold))
              #      self.flags[n] = 2
               #     continue

        #test=np.zeros((self.Hx.shape[0]),)
        #test[:]=self.Hx[:]
        #logging.info("Total HX array: %s"%test)

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

            self.KG[:,n] = PHt / self.HPHR[n]


            #if 'surface' in self.sitecode[n]:
            #    self.KG[-4:,n]=0.
            #    self.KG[3078:3082,n]=0.
            #    logging.debug('BC KG value set to zero for %s' %(self.sitecode[n]))
            #if 'aircraft' in self.sitecode[n]:
            #    self.KG[0:3078,n]=0.
            #    self.KG[3082:-4,n]=0.
            #    logging.debug('Flux KG values set to zero for %s' %(self.sitecode[n]))


            if self.may_localize[n]:
                logging.debug('Trying to localize observation %s, %i' % (self.sitecode[n], self.obs_ids[n]))
                self.localize(n)
            else:
                logging.debug('Not allowed to localize observation %s, %i' % (self.sitecode[n], self.obs_ids[n]))

            alpha = np.double(1.0) / (np.double(1.0) + np.sqrt((self.R[n]) / self.HPHR[n]))


            self.x[:] = self.x + self.KG[:,n] * res
            #logging.debug('Residual %s'%res)
            #logging.debug('obs = %s, Hx = %s,Hx(CAR) = %s, Hxp(CAR) = %s, Hx_prime_std(CAR) = %s'%(self.obs[n],self.Hx[n],self.Hx[19],test[19],self.HX_prime[19,:].std()))
            #logging.debug('New self.KG BC1 %s' %self.KG[3078,n])
            #logging.debug('New self.KG BC2 %s' %self.KG[-1,n])
            #logging.debug('New self.x BC1 %s' %self.x[3078])
            #logging.debug('New self.x BC2 %s' %self.x[-1])

            for r in range(self.nmembers):
                self.X_prime[:, r] = self.X_prime[:, r] - alpha * self.KG[:,n] * (self.HX_prime[n, r])



#WP !!!! Very important to first do all obervations from n=1 through the end, and only then update 1,...,n. The current observation
#WP      should always be updated last because it features in the loop of the adjustments !!!!

            for m in range(n + 1, self.nobs):
                #if 'aircraft' in self.sitecode[m] and not 'aircraft' in self.sitecode[n]:
                #    continue

                res = self.obs[n] - self.Hx[n]
                fac = 1.0 / (self.nmembers - 1) * (self.HX_prime[n, :] * self.HX_prime[m, :]).sum() / self.HPHR[n]
                self.Hx[m] = self.Hx[m] + fac * res
                #if n==0 and m==19:
                #    logging.debug('self.HX_prime[n, :]= %s'%self.HX_prime[n, :])
                #    logging.debug('self.HX_prime[m, :]= %s'%self.HX_prime[m, :])

                self.HX_prime[m, :] = self.HX_prime[m, :] - alpha * fac * self.HX_prime[n, :]

                #logging.debug('m = %s, Corrcoef = %s, fac = %s'%(m, np.corrcoef(self.HX_prime[n, :],self.HX_prime[m, :])[0,1],fac))



            for m in range(1, n + 1):
                #if 'aircraft' in self.sitecode[m] and not 'aircraft' in self.sitecode[n]:
                #    continue
                res = self.obs[n] - self.Hx[n]
                fac = 1.0 / (self.nmembers - 1) * (self.HX_prime[n, :] * self.HX_prime[m, :]).sum() / self.HPHR[n]
                self.Hx[m] = self.Hx[m] + fac * res
                self.HX_prime[m, :] = self.HX_prime[m, :] - alpha * fac * self.HX_prime[n, :]

            #    logging.debug('m = %s, Corrcoef = %s, fac = %s'%(m, np.corrcoef(self.HX_prime[n, :],self.HX_prime[m, :])[0,1],fac))


################### End Class CO2Optimizer ###################

if __name__ == "__main__":
    pass
