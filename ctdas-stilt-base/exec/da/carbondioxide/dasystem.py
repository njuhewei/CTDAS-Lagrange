#!/usr/bin/env python
# control.py

"""
Author : peters 

Revision History:
File created on 26 Aug 2010.

"""

import logging

################### Begin Class CO2DaSystem ###################

from da.baseclasses.dasystem import DaSystem

class CO2DaSystem(DaSystem):
    """ Information on the data assimilation system used. This is normally an rc-file with settings.
    """
    def validate(self):
        """ 
        validate the contents of the rc-file given a dictionary of required keys
        """

        needed_rc_items = ['obs.input.dir',
                           'obs.input.fname',
                           'obspack.input.id',
                           'obspack.input.dir',
                           'ocn.covariance',
                           'nparameters',
                           'bio.covariance',
                           'deltaco2.prefix',
                           'regtype']


        for k, v in self.iteritems():
            if v == 'True' : 
                self[k] = True
            if v == 'False': 
                self[k] = False

        for key in needed_rc_items:
            if not self.has_key(key):
                logging.warning('Missing a required value in rc-file : %s' % key)
        logging.debug('DA System Info settings have been validated succesfully')

################### End Class CO2DaSystem ###################


if __name__ == "__main__":
    pass
