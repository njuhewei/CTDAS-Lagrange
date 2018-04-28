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

################### End Class CO2Optimizer ###################

if __name__ == "__main__":
    pass
