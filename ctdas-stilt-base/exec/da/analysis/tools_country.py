#!/usr/bin/env python
# tools_country.py

"""
Author : peters 

Revision History:
File created on 25 Jul 2008.

This module provides an interface to go from arrays of gridded data to aggregates for countries.
It uses the country information from the UNEP database, and the database created by them for this purpose. See:

http://na.unep.net/metadata/unep/GRID/GRIDCTRY.html

The module sets up a class 'countryinfo' that holds a number of attributes. These attributes can be used to 
extract gridded information about the country. A build-in method agg_1x1 is provided for convencience. 
The routine get_countrydict() creates a dictionary with this information for a large number of countries. 

CAUTION:

The country data only covers the land areas of a nation and aggregation will exclude the fraction of a land covered
by oceans, sea, or open water. The aggregation will thus work best on arrays that are *not* defined by unit area!

"""
import sys
import cPickle
import os
sys.path.append('../../')
rootdir = os.getcwd().split('da/')[0]
analysisdir = os.path.join(rootdir, 'da/analysis')

from numpy import sum, array
try:
    from dbfpy import dbf
except:
    print "the python DBF lib might be needed, please install from:"
    print "http://dbfpy.sourceforge.net/"
    print " Trying to complete anyway..."

sys.path.append('../../')
from da.analysis.tools_regions import globarea

EU25 = ["Austria", "Belgium", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland", "Portugal", "Slovakia", "Slovenia", "Spain", "Sweden", "United Kingdom"]
EU27 = ["Austria", "Belgium", "Bulgaria", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "United Kingdom"]
G8 = [ "France", "Germany", "Italy", "Japan", "United Kingdom", "United States"]
annex1 = [ "Australia", "Austria", "Belarus", "Belgium", "Bulgaria", "Canada", "Croatia", "Czech Republic", "Denmark", "European Union", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Japan", "Latvia", "Liechtenstein", "Lithuania", "Luxembourg", "Monaco", "Netherlands", "New Zealand", "Norway", "Poland", "Portugal", "Romania", "Russia", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland", "Turkey", "Ukraine", "United Kingdom", "United States"]
annex2 = [ "Australia", "Austria", "Belgium", "Canada", "Denmark", "Finland", "France", "Germany", "Greece", "Iceland", "Ireland", "Italy", "Japan", "Luxembourg", "Netherlands", "New Zealand", "Norway", "Portugal", "Spain", "Sweden", "Switzerland", "Turkey", "United Kingdom", "United States"]

class countryinfo(object):
    """ Defines the information about 1x1 gridboxes covering each country """

    def __init__(self, code):
        self.country = code
        self.ngrids = 0
        self.gridij = []
        self.gridnr = []
        self.gridcoords = []

        self.gridlandmask = []

        self.gridsharedocean = []
        self.gridsharedborder = []

    def add_gridinfo(self, grid_i, grid_j, gridlandmask, shared_border, shared_water):
        """ add information from one gridbox to the object """ 

        self.gridij.append((grid_j, grid_i,))  # add tuple with j,i coordinates
        lon = -180 + grid_i + 0.5
        lat = -90 + grid_j + 0.5
        self.gridcoords.append((lat, lon,))  # add tuple with lon,lat coordinates
        self.gridnr.append(grid_j * 360 + grid_i)          # add grid number for take() function

        self.gridlandmask.append(gridlandmask)  # this gives the total fraction of land

        self.gridsharedocean.append(shared_water)
        self.gridsharedborder.append(shared_border)
        self.ngrids = self.ngrids + 1

    def agg_1x1(self, field):
        """ aggregate a 1x1 input field to country total """

        #print field.take(self.gridnr)
        #print self.gridlandmask

        return (field.take(self.gridnr) * array(self.gridlandmask)).sum()

    def __str__(self):
        return '  Country name         : %s\n' % self.country + \
               '  Number of gridpoints : %d ' % self.ngrids 
        #'  Gridpoint indices : %s ' % self.gridnr

def fix_eu(rec):
    """ fix Czech Republic and Slovakia manually """

    alternative_slov = {
    '140202': (2.0, 28.1), \
    '140201': (2.0, 34.0), \
    '140200': (2.0, 44.0), \
    '140199': (3.0, 25.5), \
    '139203': (3.0, 18.9), \
    '139202': (2.0, 57.5), \
    '139201': (2.0, 59.7), \
    '139200': (2.0, 87.2), \
    '139199': (1.0, 100.0), \
    '139198': (2.0, 72.8), \
    '138198': (3.0, 7.7), \
    '138199':(2.0, 10.0)  }
    
    alternative_czech = {
    '141193': (2.0, 23.0), \
    '141194': (2.0, 62.1), \
    '141195': (3.0, 89.5), \
    '141196': (2.0, 79.4), \
    '141197': (2.0, 42.3), \
    '141198': (2.0, 24.5), \
    '141199': (2.0, 0.1), \
    '140193': (2.0, 20.6), \
    '140194': (2.0, 88.9), \
    '140195': (1.0, 100.0), \
    '140196': (1.0, 100.0), \
    '140197': (1.0, 100.0), \
    '140198': (1.0, 100.0), \
    '140199': (3.0, 50.0), \
    '139195': (2.0, 70.6), \
    '139196': (2.0, 12.4), \
    '139197': (2.0, 30.9), \
    '139198': (2.0, 25.0) }

    id = str(int(rec['GRID']))
    for dict in [alternative_slov, alternative_czech]:
        if id in dict:
            rec['COVER_ID'] = dict[id][0]
            rec['RATE_IN_GR'] = dict[id][1]

    return rec

def get_countrydict():
    """ Create a dictionary with grid-to-country information from a dbf file"""

    countrydict = countryinfo('Test')

    file = os.path.join(analysisdir,'country_dictionary.dat')
    
    try:
        countrydict = cPickle.load(open(file, 'rb'))
    except:     
        db = dbf.Dbf(os.path.join(analysisdir,'GRIDCTRY.DBF'))

        countrydict = {}
        for n, rec in enumerate(db):
            code = rec['COUNTRY']
            gridid = str(int(rec['GRID']))

            if code in ['Czech Republic', 'Slovakia']:
                rec = fix_eu(rec)

            rate_in_gr = rec['RATE_IN_GR'] * 1.e-2

            i = int(gridid[-3::])
            j = int(gridid[0:-3])
            lat = -91 + j + 0.5
            lon = -181 + i + 0.5
            if code in countrydict: 
                a = countrydict[code]
            else:
                a = countryinfo(code)


            shared_border = False
            shared_water = False
            if rec['COVER_ID'] == 0.0:
                shared_border = False
                shared_water = True
            if rec['COVER_ID'] >= 2.0:
                shared_border = True
            if rec['COVER_ID'] >= 10.0:
                shared_water = True

            a.add_gridinfo(i - 1, j - 1, rate_in_gr, shared_border, shared_water)

            countrydict[code] = a

        db.close()

        cPickle.dump(countrydict, open(file, 'wb'), -1)

    return countrydict

if __name__ == "__main__":

    countrydict = get_countrydict()

    area = globarea()

    areas = []
    for k, v in countrydict.items():
        ar = v.agg_1x1(area) / 1.e6
        areas.append((ar, k))

    areas.sort()
    areas.reverse()
    for a in areas: print a

    v = countrydict['Ocean']
    print v.agg_1x1(area)
    v = countrydict['Netherlands']
    print v.agg_1x1(area)
    v = countrydict['Slovakia']
    print v.agg_1x1(area)
    v = countrydict['Czech Republic']
    print v.agg_1x1(area)
    v = countrydict['Czechoslovakia']
    print v.agg_1x1(area)








