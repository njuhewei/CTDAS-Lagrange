#!/usr/bin/env python
# MolefractionToPermil.py

"""
Author : ivar 

Revision History:
File created on 11 May 2012.

"""


def molefraction_to_permil(filename,simulated):
    """ Converts 13C mole fractions to permil values"""
    import da.tools.io4 as io
    import logging
    import numpy as np

    pdb         = 0.011112

    trlength=len(simulated[0,:])
    memlength=trlength/2
    simulated=simulated #*1.e6 #convert to ppm
 #   np.set_printoptions(threshold=np.nan)
 #   msg='simulated shape',simulated.shape;logging.info(msg)
 #   msg='simulated',simulated[0,0],simulated[0,40];logging.info(msg)
 #   msg='simulated',simulated ;logging.info(msg)
    simulated=np.float64(simulated)
    simulated[:,memlength:trlength]=1000.*((simulated[:,memlength:trlength]/simulated[:,0:memlength]/pdb)-1.)
 #   msg='simulated',simulated[0,0],simulated[0,40],memlength,trlength,pdb;logging.info(msg)
 #   msg='simulated',simulated ;logging.info(msg)




    return simulated


if __name__ == "__main__":
    pass
