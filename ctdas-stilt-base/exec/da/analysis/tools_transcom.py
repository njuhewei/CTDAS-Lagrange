#!/usr/bin/env python

import os
import sys
sys.path.append('../../')
rootdir = os.getcwd().split('da/')[0]
analysisdir = os.path.join(rootdir, 'da/analysis')

from string import join, split
from numpy import array, identity, zeros, arange, dot
import da.tools.io4 as io

# Get masks of different region definitions

matrix_file = os.path.join(analysisdir, 'copied_regions.nc')
cdf_temp = io.CT_CDF(matrix_file, 'read')
transcommask = cdf_temp.get_variable('transcom_regions')
if transcommask.max() < 23:
    if 'transcom_regions_original' in cdf_temp.variables:
        transcommask = cdf_temp.get_variable('transcom_regions_original')
olson240mask = cdf_temp.get_variable('regions')
olsonmask = cdf_temp.get_variable('land_ecosystems')
#oifmask = cdf_temp.get_variable('ocean_regions')
dummy = cdf_temp.close()

matrix_file = os.path.join(analysisdir, 'copied_regions_extended.nc')
cdf_temp = io.CT_CDF(matrix_file, 'read')
olson_ext_mask = cdf_temp.get_variable('regions')
dummy = cdf_temp.close()

# Names and short names of TransCom regions

transshort = []
transnams = []
transland = []
temp = open(os.path.join(analysisdir, 't3_region_names'), 'r').readlines()
for line in temp:
    items = line.split()
    if items:
        num, abbr, name = (items[0], items[1], join(items[2:]),)
        transnams.append(name.strip('"'))
        transshort.append(abbr)
        if abbr.startswith('T'):
            transland.append(name.strip('"'))

transnams.append("Non-optimized")
transshort.append("I-NNOP")

# Names and short names of Olson regions

olsonnams = []
olsonshort = []
temp = open(os.path.join(analysisdir, 'olson19_region_names'), 'r').readlines()
for line in temp:
    items = line.split()
    if items:
        num, abbr, name = (items[0], items[1], join(items[2:]),)
        olsonnams.append(name.strip('"'))
        olsonshort.append(abbr)

olsonextnams = []
matrix_file = os.path.join(analysisdir, 'copied_regions_extended.nc')
cdf_temp = io.CT_CDF(matrix_file, 'read')
keys = cdf_temp.ncattrs()
keys.sort()
for k in keys:
    if 'Region' in k:
        olsonextnams.append(getattr(cdf_temp, k))
cdf_temp.close()

ext_transnams = []
ext_transshort = []
ext_transcomps = []

# Get names of aggregated regions for post aggregation

matrix_file = os.path.join(analysisdir, 'postagg_definitions.nc')
cdf_temp = io.CT_CDF(matrix_file, 'read')
xform = cdf_temp.get_variable('xform')
keys = cdf_temp.ncattrs()

keys.sort()
for k in keys:
    if 'longname' in k:
        ext_transnams.append(getattr(cdf_temp, k))
    if 'shortname' in k:
        ext_transshort.append(getattr(cdf_temp, k))
    if 'component' in k:
        ext_transcomps.append(map(int, getattr(cdf_temp, k).split(',')))

cdf_temp.close()


# Names of the ocean inversion flux regions, to go along with oifmask

oifnams = ['(1) NOCN Arctic Ocean', \
         '(2) NAH  North Atlantic (49 - 76N)', \
         '(3) NAM  North Atlantic (36 - 49N)', \
         '(4) NAL  North Atlantic (18 - 36N)', \
         '(5) NAT  North Atlantic ( 0 - 18N)', \
         '(6) SAT  South Atlantic ( 0 - 18S)', \
         '(7) SAL  South Atlantic (18 - 31S)', \
         '(8) SAM  South Atlantic (31 - 44S)', \
         '(9) SAH  South Atlantic (44 - 58S)', \
         '(10) SOCN Southern Ocean (S of 58S)', \
         '(11) NPHW North Pacific (N of 49N, W of 195E)', \
         '(12) NPHE North Pacific (N of 36N, E of 195E)', \
         '(13) NPK  North Pacific (Kuroshio Extension)', \
         '(14) NPLW North Pacific (18N - K.Ext, W of 195E)', \
         '(15) NPLE North Pacific (18 - 36N, E of 195E)', \
         '(16) NPTW North Pacific ( 0 - 18N, W of 199E)', \
         '(17) NPTE North Pacific ( 0 - 18N, E of 199E)', \
         '(18) SPTW South Pacific ( 0 - 18S, W of 199E)', \
         '(19) SPTE South Pacific ( 0 - 18S, E of 199E)', \
         '(20) SPLW South Pacific (18 - 31S, W of 233E)', \
         '(21) SPLE South Pacific (18 - 31S, E of 233E)', \
         '(22) SPMW South Pacific (31 - 44S, W of 248E)', \
         '(23) SPME South Pacific (31 - 44S, E of 248E, W of 278E)', \
         '(24) SPMC South Pacific (31 - 44S, coastal E of 278E)', \
         '(25) SPH  South Pacific (44 - 58S)              ', \
         '(26) NI   North Indian', \
         '(27) SIT  South Indian (0 - 18S)', \
         '(28) SIL  South Indian (18 - 31S)', \
         '(29) SIM  South Indian (31 - 44S)', \
         '(30) SIH  South Indian (44 - 58S)']

oiflocs = [ (200, 80,), \
          (330, 55,), \
          (330, 40,), \
          (330, 22,), \
          (330, 8,), \
          (350, -12,), \
          (350, -27,), \
          (350, -40,), \
          (350, -53,), \
          (200, -70,), \
          (178, 54,), \
          (210, 40,), \
          (165, 38,), \
          (178, 25,), \
          (215, 25,), \
          (170, 8,), \
          (230, 8,), \
          (175, -10,), \
          (240, -10,), \
          (195, -27,), \
          (265, -27,), \
          (195, -40,), \
          (262, -40,), \
          (283, -40,), \
          (220, -53,), \
          (68, 8,), \
          (75, -10,), \
          (75, -27,), \
          (75, -40,), \
          (75, -53,)]


translocs = [ (-177, 0), \
            (-92, 53,), \
            (-108, 34,), \
            (-66, 4,), \
            (-50, -17,), \
            (15, 17,), \
            (26, -12,), \
            (84, 63,), \
            (103, 30,), \
            (115, 0,), \
            (132, -25,), \
            (9, 50,), \
            (-174, 46,), \
            (136, 6,), \
            (-108, 6,), \
            (-123, -15,), \
            (-32, 58,), \
            (-32, 38,), \
            (-32, 0,), \
            (-32, -38,), \
            (-14, -65,), \
            (68, 2,)]

#olsonshort=[str(name.split()[1:2]).join('  ') for name in olsonnams]
old_olsonshort = [join(split(name, ' ')[1:2], ' ') for name in olsonnams]

olsonlabs = ['Conifer Forest', 'Broadleaf Forest', 'Mixed Forest', 'Grass/Shrub', 'Tropical Forest', 'Scrub/Woods', 'Semitundra', 'Fields/Woods/\nSavanna', \
     'Northern Taiga', 'Forest/Field', 'Wetland', 'Deserts', 'Shrub/Tree/\nSuc ', 'Crops', 'Conifer\n Snowy/Coastal', \
     'Wooded tundra', 'Mangrove', 'Ice and \nPolar desert', 'Water'] 

ecmwfnams = [ ' 1 CRPSMF Crops, mixed farming', \
            ' 2 SHGRSS Short Grass', \
            ' 3 EVNDLF Evergreen Needleleaf', \
            ' 4 DECNDLF Deciduous Needleleaf', \
            ' 5 EVBRDLF Evergreen Broadleaf', \
            ' 6 DECBRLF Deciduous Broadleaf', \
            ' 7 TLGRSS Tall Grass', \
            ' 8 DES Desert', \
            ' 9 TDR Tundra', \
            '10 IRRCR Irrigated Crops', \
            '11 SMDES Semidesert', \
            '12 ICE Ice Caps', \
            '13 BGM Bogs and Marches', \
            '14 INW Inland Water', \
            '15 OCE Ocean', \
            '16 EVSHRB Evergreen Shrubs', \
            '17 DECSHR Deciduous shrubs', \
            '18 MXFRST Mixed Forest', \
            '19 INTFRST Interrupted Forest'] 

ecmwfshort = [str(name.split()[1:2]).join('  ') for name in ecmwfnams]

ecmwflabs = ['Crops, mixed farming', 'Short Grass', 'Evergreen Needleleaf', 'Deciduous Needleleaf', 'Evergreen Broadleaf', \
      'Deciduous Broadleaf', 'Tall Grass', 'Desert', \
     'Tundra', 'Irrigated Crops', 'Semidesert', 'Ice Caps', 'Bogs and Marches', 'Inland Water', 'Ocean', \
     'Evergreen Shrubs', 'Deciduous shrubs', 'Mixed Forest', 'Interrupted Forest'] 

a = array([\
  0 , 0 , 0 , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 1  , 1  , 1  , 1  , 1  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0 , \
  0 , 0 , 0 , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 1  , 0  , 1  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0 , \
  0 , 0 , 0 , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 1  , 0  , 1  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0 , \
  0 , 0 , 0 , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 1  , 1  , 1  , 1  , 1  , 0  , 0  , 0  , 0  , 0  , 0 , \
  1 , 1 , 0 , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0 , \
  0 , 0 , 1 , 1  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0 , \
  0 , 0 , 0 , 0  , 1  , 1  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0 , \
  0 , 0 , 0 , 0  , 0  , 0  , 1  , 1  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0 , \
  0 , 0 , 0 , 0  , 0  , 0  , 0  , 0  , 1  , 1  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 1  , 0  , 0  , 0  , 0  , 1 , \
  0 , 0 , 0 , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 1  , 1  , 0  , 0  , 0 , \
  0 , 0 , 0 , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 0  , 1  , 1  , 0])
O30to11 = a.reshape(11, 30).transpose()

O11to11 = identity(11)

ntcland = 11 # TC standard
ntcocean = 11 # TC standard

Ols259_to_TC23 = zeros((259, 23), float)
Ols240_to_TC23 = zeros((240, 23), float)
Ols221_to_TC23 = zeros((221, 23), float)
for i in arange(ntcland):
    Ols259_to_TC23[i * 19:(i + 1) * 19, i] = 1.0

Ols259_to_TC23[190:228, 10] = 1.0  # Europe
Ols259_to_TC23[228:258, 11:22] = O30to11
for i in arange(ntcland):
    Ols240_to_TC23[i * 19:(i + 1) * 19, i] = 1.0
Ols240_to_TC23[209:239, 11:22] = O30to11
for i in arange(ntcland):
    Ols221_to_TC23[i * 19:(i + 1) * 19, i] = 1.0
Ols221_to_TC23[209:220:, 11:22] = O11to11

Ols221_to_TC23[220, 22] = 1.0
Ols240_to_TC23[239, 22] = 1.0
Ols259_to_TC23[258, 22] = 1.0


ntcland = 11 # TC standard
ntcocean = 11 # TC standard

ExtendedTCRegionsFile = 'postagg_definitions.nc'

def ExtendedTCRegions(data, cov=False):
    """ convert to extended transcom shaped regions"""

    nparams = data.shape[-1]
    if nparams != 23:
        raise ValueError('Do not know how to convert %s regions to 37 extended transcom regions' % (nparams,))
    M = xform

    if not cov:
        return dot(array(data).squeeze(), M)
    else:
        try:
            return M.transpose().dot(data).dot(M)
        except:
            return dot(dot(M.transpose(), data), M) #Huygens fix

def cov2corr(A):
    b = 1. / sqrt(A.diagonal())
    return A * dot(b[:, newaxis], b[newaxis, :])

    """ function projects 1x1 degree map onto TransCom regions by adding gridboxes over larger areas """
    from hdf2field import Sds2field 
    import cPickle
    import os
    from plottools import rebin

    transcommapfile = 'tc_land11_oif30.hdf'
    transcomconversionfile = 'map_to_tc.pickle'
    try:
        regionselect = cPickle.load(open(transcomconversionfile, 'rb'))
    except:
        # read map from NetCDF
        print '[in map_to_tc() in tctools.py:] ' + \
              'Creating conversion map and pickle file for future quick use, patience please...'
        map = Sds2field(transcommapfile, 'tc_region')

        # create dictionary for region <-> map conversions based on 1x1 map
        
        regs = {}
        nregions = map.max()
        for r in arange(1, nregions + 1):
            sel = (map.flat == r).nonzero()
            if len(sel[0]) > 0:
                regs[r] = sel
        regionselect = regs
        dummy = cPickle.dump(regionselect, open(transcomconversionfile, 'wb'), -1)
    
    result = zeros(len(regionselect.keys()), float)
    for k, v in regionselect.iteritems():
        result[k - 1] = data.ravel().take(v).sum()
    return result

    """ return name of region number reg """

    if longnames:
        econames = olsonnams
    else :
        econames = olsonshort

    if tc:
        return (transnams[reg - 1],)
    elif eco:
        if reg > rundat.npparameters:
            raise IOError, 'Region number exceeds definitions'
        elif reg > rundat.n_land and reg != rundat.nparameters:
            ret = ('Ocean', oifnams[reg - rundat.n_land - 1])
        elif reg > 209 and reg <= rundat.n_land:
            ret = ('Europe', econames[(reg - 1) % 19] + "_East")
        elif reg == rundat.nparameters:
            ret = (transnams[-1])
        else:
            ret = (transnams[(reg - 1) / 19], econames[(reg - 1) % 19])
        return ret
    elif olson:
        return (econames[(reg - 1) % 19],)

if __name__ == '__main__':
    print transnams
    print transshort
    print ext_transnams
    print ext_transshort
    print olsonnams
    print olsonshort
    print ext_transcomps
    print olsonextnams    

