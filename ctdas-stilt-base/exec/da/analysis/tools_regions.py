#!/usr/bin/env python

import numpy as np
import cPickle
from da.analysis.tools_transcom import *

# Aggregated olson ecosystem regions for CT Europe

aggregates = {
            "Forest" : (1, 2, 3, 5, 8, 10,) , \
            "Grass" : (4, 6, 13) , \
            "Crops": (14,) , \
            "Tundra" : (7, 9, 16) , \
            "Coastal" : (11, 15, 17) , \
            "IceWaterDeserts" : (12, 18, 19) \
}

ext_econams = [a for a, b in aggregates.iteritems()]
ext_ecocomps = [b for a, b in aggregates.iteritems()]

eco19_to_ecosums = zeros((19, 6), float)
for i, k in enumerate(ext_ecocomps):
    indices = [x - 1 for x in k]
    eco19_to_ecosums[:, i].put(indices, 1.0)

##### END OF REGION DEFINITIONS

def state_to_grid(values, regionmap, reverse=False, avg=False, mapname=None):
    """ 
    This method converts parameters from a CarbonTracker StateVector object to a gridded map of linear multiplication values. These
    can subsequently be used in the transport model code to multiply/manipulate fluxes

    """
    nregions = regionmap.max()
    try:
        if not mapname: 
            raise Exception 

        regionselect = cPickle.load(open('%s_regiondict.pickle' % mapname, 'rb'))
    except:

        # dictionary for region <-> map conversions
        regs = {}
        for r in np.arange(1, nregions + 1):
            sel = (regionmap.flat == r).nonzero()
            if len(sel[0]) > 0: 
                regs[r] = sel

        regionselect = regs
        
        cPickle.dump(regionselect, open('%s_regiondict.pickle' % mapname, 'wb'), -1)
        print 'Pickling region map'

    if reverse:
        """ project 1x1 degree map onto ecoregions """

        result = np.zeros(nregions, float)
        for k, v in regionselect.iteritems():
            if avg: 
                result[k - 1] = values.ravel().take(v).mean()
            else : 
                result[k - 1] = values.ravel().take(v).sum()
        return result

    else:
        """ project ecoregion properties onto 1x1 degree map """

        result = np.zeros((180, 360,), float)
        for k, v in regionselect.iteritems():
            result.put(v, values[k - 1])

        return result

def globarea(im=360, jm=180, silent=True):
    """ Function calculates the surface area according to TM5 definitions"""

    radius = 6.371e6  # the earth radius in meters
    deg2rad = np.pi / 180.
    g = 9.80665 

    dxx = 360.0 / im * deg2rad 
    dyy = 180.0 / jm * deg2rad 
    lat = np.arange(-90 * deg2rad, 90 * deg2rad, dyy)
    dxy = dxx * (np.sin(lat + dyy) - np.sin(lat)) * radius ** 2
    area = np.resize(np.repeat(dxy, im, axis=0) , [jm, im])
    if not silent:
        print 'total area of field = ', np.sum(area.flat)
        print 'total earth area    = ', 4 * np.pi * radius ** 2
    return area

