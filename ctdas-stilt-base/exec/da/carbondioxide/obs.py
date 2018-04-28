#!/usr/bin/env python
# obs.py

"""
Author : peters 

Revision History:
File created on 28 Jul 2010.

"""
import os
import sys
import logging
#from da.baseclasses.statevector import filename
import datetime as dtm
from string import strip
from numpy import array, logical_and

sys.path.append(os.getcwd())
sys.path.append('../../')

identifier = 'CarbonTracker CO2 mole fractions'
version = '0.0'

from da.baseclasses.obs import Observations
import da.tools.io4 as io
import da.tools.rc as rc

################### Begin Class CO2Observations ###################

class CO2Observations(Observations):
    """ an object that holds data + methods and attributes needed to manipulate mole fraction values """

    def setup(self, dacycle):
        self.startdate = dacycle['time.sample.start']
        self.enddate = dacycle['time.sample.end']
        
        sfname = dacycle.dasystem['obs.input.fname']
        if sfname.endswith('.nc'):
            filename = os.path.join(dacycle.dasystem['obs.input.dir'], sfname)
        else:
            filename = os.path.join(dacycle.dasystem['obs.input.dir'], sfname + '.' + self.startdate.strftime('%Y%m%d') + '.nc')

        if not os.path.exists(filename):
            msg = 'Could not find  the required observation input file (%s) ' % filename
            logging.error(msg)
            raise IOError, msg
        else:
            self.obs_filename = filename
        self.datalist = []

    def add_observations(self):
        """ Returns a MoleFractionList holding individual MoleFractionSample objects for all obs in a file
      
            The CarbonTracker mole fraction files are provided as one long list of obs for all possible dates. So we can 
            either:
            
            (1) read all, and the subselect the data we will use in the rest of this cycle
            (2) Use nco to make a subset of the data
            
            For now, we will stick with option (1) 
        
        """
        ncf = io.ct_read(self.obs_filename, 'read')
        idates = ncf.get_variable('date_components')
        dates = array([dtm.datetime(*d) for d in idates])

        subselect = logical_and(dates >= self.startdate, dates <= self.enddate).nonzero()[0]

        dates = dates.take(subselect, axis=0)
        
        ids = ncf.get_variable('id').take(subselect, axis=0)
        evn = ncf.get_variable('eventnumber').take(subselect, axis=0)
        evn = [s.tostring().lower() for s in evn]
        evn = map(strip, evn)
        sites = ncf.get_variable('site').take(subselect, axis=0)
        sites = [s.tostring().lower() for s in sites]
        sites = map(strip, sites)
        lats = ncf.get_variable('lat').take(subselect, axis=0)
        lons = ncf.get_variable('lon').take(subselect, axis=0)
        alts = ncf.get_variable('alt').take(subselect, axis=0)
        obs = ncf.get_variable('obs').take(subselect, axis=0) * 1.e-6
        logging.info("Converting observed values from ppm to mol/mol!!!!")
        species = ncf.get_variable('species').take(subselect, axis=0)
        species = [s.tostring().lower() for s in species]
        species = map(strip, species)
        strategy = ncf.get_variable('sampling_strategy').take(subselect, axis=0)
        flags = ncf.get_variable('NOAA_QC_flags').take(subselect, axis=0)
        flags = [s.tostring().lower() for s in flags]
        flags = map(strip, flags)
        flags = [int(f == '...') for f in flags]
        ncf.close()

        logging.debug("Successfully read data from obs file (%s)" % self.obs_filename)

        for n in range(len(dates)): 
            self.datalist.append(MoleFractionSample(ids[n], dates[n], sites[n], obs[n], 0.0, 0.0, 0.0, 0.0, flags[n], alts[n], lats[n], lons[n], evn[n], species[n], strategy[n], 0.0, self.obs_filename))
        logging.debug("Added %d observations to the Data list" % len(dates))

    def add_simulations(self, filename, silent=True):
        """ Adds model simulated values to the mole fraction objects """


        if not os.path.exists(filename):
            msg = "Sample output filename for observations could not be found : %s" % filename
            logging.error(msg)
            logging.error("Did the sampling step succeed?")
            logging.error("...exiting")
            raise IOError, msg

        ncf = io.ct_read(filename, method='read')
        ids = ncf.get_variable('obs_num')
        simulated = ncf.get_variable('flask')
        ncf.close()
        logging.info("Successfully read data from model sample file (%s)" % filename)

        obs_ids = self.getvalues('id')

        obs_ids = obs_ids.tolist()
        ids = map(int, ids)

        missing_samples = []

        for idx, val in zip(ids, simulated): 
            if idx in obs_ids:
                index = obs_ids.index(idx)
                #print id,val,val.shape
                self.datalist[index].simulated = val
            else:     
                missing_samples.append(idx)

        if not silent and missing_samples != []:
            logging.warning('Model samples were found that did not match any ID in the observation list. Skipping them...')
            #msg = '%s'%missing_samples ; logging.warning(msg)

        logging.debug("Added %d simulated values to the Data list" % (len(ids) - len(missing_samples)))

    def write_sample_coords(self, obsinputfile):
        """ 
            Write the information needed by the observation operator to a file. Return the filename that was written for later use

        """
        f = io.CT_CDF(obsinputfile, method='create')
        logging.debug('Creating new observations file for ObservationOperator (%s)' % obsinputfile)

        dimid = f.add_dim('obs', len(self.datalist))
        dim200char = f.add_dim('string_of200chars', 200)
        dimcalcomp = f.add_dim('calendar_components', 6)

        if len(self.datalist) == 0:
            f.close()
            #return obsinputfile

        data = self.getvalues('id')

        savedict = io.std_savedict.copy() 
        savedict['name'] = "obs_num"
        savedict['dtype'] = "int"
        savedict['long_name'] = "Unique_Dataset_observation_index_number"
        savedict['units'] = ""
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['comment'] = "Unique index number within this dataset ranging from 0 to UNLIMITED."
        f.add_data(savedict)

        data = [[d.year, d.month, d.day, d.hour, d.minute, d.second] for d in self.getvalues('xdate') ]

        savedict = io.std_savedict.copy() 
        savedict['dtype'] = "int"
        savedict['name'] = "date_components"
        savedict['units'] = "integer components of UTC date/time"
        savedict['dims'] = dimid + dimcalcomp
        savedict['values'] = data
        savedict['missing_value'] = -9
        savedict['comment'] = "Calendar date components as integers. Times and dates are UTC." 
        savedict['order'] = "year, month, day, hour, minute, second"
        f.add_data(savedict)

        data = self.getvalues('lat')

        savedict = io.std_savedict.copy() 
        savedict['name'] = "latitude"
        savedict['units'] = "degrees_north"
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['missing_value'] = -999.9
        f.add_data(savedict)

        data = self.getvalues('lon')

        savedict = io.std_savedict.copy() 
        savedict['name'] = "longitude"
        savedict['units'] = "degrees_east"
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['missing_value'] = -999.9
        f.add_data(savedict)

        data = self.getvalues('height')

        savedict = io.std_savedict.copy() 
        savedict['name'] = "altitude"
        savedict['units'] = "meters_above_sea_level"
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['missing_value'] = -999.9
        f.add_data(savedict)

        data = self.getvalues('samplingstrategy')

        savedict = io.std_savedict.copy() 
        savedict['dtype'] = "int"
        savedict['name'] = "sampling_strategy"
        savedict['units'] = "NA"
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['missing_value'] = -9
        f.add_data(savedict)

        data = self.getvalues('evn')

        savedict = io.std_savedict.copy() 
        savedict['dtype'] = "char"
        savedict['name'] = "obs_id"
        savedict['units'] = "NOAA database identifier"
        savedict['dims'] = dimid + dim200char
        savedict['values'] = data
        savedict['missing_value'] = '!'
        f.add_data(savedict)

        f.close()

        logging.debug("Successfully wrote data to obs file")
        logging.info("Sample input file for obs operator now in place [%s]" % obsinputfile)

        


    def add_model_data_mismatch(self, filename):
        """ 
            Get the model-data mismatch values for this cycle.

                (1) Open a sites_weights file
                (2) Parse the data
                (3) Compare site list against data
                (4) Take care of double sites, etc

        """
        

        

        if not os.path.exists(filename):
            msg = 'Could not find  the required sites.rc input file (%s)' % filename
            logging.error(msg)
            raise IOError, msg
        else:
            self.sites_file = filename

        sites_weights = rc.read(self.sites_file)

        self.rejection_threshold = int(sites_weights['obs.rejection.threshold'])
        self.global_R_scaling = float(sites_weights['global.R.scaling'])
        self.n_site_categories = int(sites_weights['n.site.categories'])
        self.n_sites_active = int(sites_weights['n.sites.active'])
        self.n_sites_moved = int(sites_weights['n.sites.moved'])

        logging.debug('Model-data mismatch rejection threshold: %d ' % self.rejection_threshold)
        logging.debug('Model-data mismatch scaling factor     : %f ' % self.global_R_scaling)
        logging.debug('Model-data mismatch site categories    : %d ' % self.n_site_categories)
        logging.debug('Model-data mismatch active sites       : %d ' % self.n_sites_active)
        logging.debug('Model-data mismatch moved sites        : %d ' % self.n_sites_moved)
   
        cats = [k for k in sites_weights.keys() if 'site.category' in k] 

        SiteCategories = {}
        for key in cats:
            name, error, may_localize, may_reject = sites_weights[key].split(';')
            name = name.strip().lower()
            error = float(error)
            may_reject = ("TRUE" in may_reject.upper())
            may_localize = ("TRUE" in may_localize.upper())
            SiteCategories[name] = {'error':error, 'may_localize':may_localize, 'may_reject':may_reject}
            #print name,SiteCategories[name]


        active = [k for k in sites_weights.keys() if 'site.active' in k] 

        site_info = {}
        for key in active:
            sitename, sitecategory = sites_weights[key].split(';')
            sitename = sitename.strip().lower()
            sitecategory = sitecategory.strip().lower()
            site_info[sitename] = SiteCategories[sitecategory]
            #print sitename,site_info[sitename]

        for obs in self.datalist:
            obs.mdm = 1000.0  # default is very high model-data-mismatch, until explicitly set by script
            if site_info.has_key(obs.code): 
                logging.debug("Observation found (%s)" % obs.code)
                obs.mdm = site_info[obs.code]['error'] * self.global_R_scaling
                obs.may_localize = site_info[obs.code]['may_localize']
                obs.may_reject = site_info[obs.code]['may_reject']
            else:
                logging.warning("Observation NOT found (%s, %s), please check sites.rc file  (%s)  !!!" % (obs.code, identifier, self.sites_file))
                obs.flag = 99

            # Add site_info dictionary to the Observations object for future use

            self.site_info = site_info

    def write_sample_auxiliary(self, auxoutputfile):
        """ 
            Write selected information contained in the Observations object to a file. 

        """
        
        f = io.CT_CDF(auxoutputfile, method='create')
        logging.debug('Creating new auxiliary sample output file for postprocessing (%s)' % auxoutputfile)

        dimid = f.add_dim('obs', len(self.datalist))
        dim200char = f.add_dim('string_of200chars', 200)
        dimcalcomp = f.add_dim('calendar_components', 6)

        if len(self.datalist) == 0:
            f.close()
            #return outfile

        data = self.getvalues('id')

        savedict = io.std_savedict.copy() 
        savedict['name'] = "obs_num"
        savedict['dtype'] = "int"
        savedict['long_name'] = "Unique_Dataset_observation_index_number"
        savedict['units'] = ""
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['comment'] = "Unique index number within this dataset ranging from 0 to UNLIMITED."
        f.add_data(savedict)

        data = [[d.year, d.month, d.day, d.hour, d.minute, d.second] for d in self.getvalues('xdate')]

        savedict = io.std_savedict.copy() 
        savedict['dtype'] = "int"
        savedict['name'] = "date_components"
        savedict['units'] = "integer components of UTC date/time"
        savedict['dims'] = dimid + dimcalcomp
        savedict['values'] = data
        savedict['missing_value'] = -9
        savedict['comment'] = "Calendar date components as integers. Times and dates are UTC." 
        savedict['order'] = "year, month, day, hour, minute, second"
        f.add_data(savedict)

        data = self.getvalues('obs')

        savedict = io.std_savedict.copy()
        savedict['name'] = "observed"
        savedict['long_name'] = "observedvalues"
        savedict['units'] = "mol mol-1"
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['comment'] = 'Observations used in optimization'
        f.add_data(savedict)

        data = self.getvalues('mdm')

        savedict = io.std_savedict.copy()
        savedict['name'] = "modeldatamismatch"
        savedict['long_name'] = "modeldatamismatch"
        savedict['units'] = "[mol mol-1]"
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['comment'] = 'Standard deviation of mole fractions resulting from model-data mismatch'
        f.add_data(savedict)

        data = self.getvalues('simulated') 

        dimmembers = f.add_dim('members', data.shape[1])

        savedict = io.std_savedict.copy()
        savedict['name'] = "modelsamples"
        savedict['long_name'] = "modelsamples for all ensemble members"
        savedict['units'] = "mol mol-1"
        savedict['dims'] = dimid + dimmembers
        savedict['values'] = data.tolist()
        savedict['comment'] = 'simulated mole fractions based on optimized state vector'
        f.add_data(savedict)

        data = self.getvalues('fromfile') 

        savedict = io.std_savedict.copy()
        savedict['name'] = "inputfilename"
        savedict['long_name'] = "name of file where original obs data was taken from"
        savedict['dtype'] = "char"
        savedict['dims'] = dimid + dim200char
        savedict['values'] = data
        savedict['missing_value'] = '!'
        f.add_data(savedict)

        data = self.getvalues('code') 

        savedict = io.std_savedict.copy()
        savedict['name'] = "sitecode"
        savedict['long_name'] = "site code propagated from observation file"
        savedict['dtype'] = "char"
        savedict['dims'] = dimid + dim200char
        savedict['values'] = data
        savedict['missing_value'] = '!'
        f.add_data(savedict)

        f.close()

        logging.debug("Successfully wrote data to auxiliary sample output file (%s)" % auxoutputfile)

        #return outfile


################### End Class CO2Observations ###################



################### Begin Class MoleFractionSample ###################

class MoleFractionSample(object):
    """ 
        Holds the data that defines a mole fraction Sample in the data assimilation framework. Sor far, this includes all
        attributes listed below in the __init__ method. One can additionally make more types of data, or make new
        objects for specific projects.

    """
    
    def __init__(self, idx, xdate, code='XXX', obs=0.0, simulated=0.0, resid=0.0, hphr=0.0, mdm=0.0, flag=0, height=0.0, lat= -999., lon= -999., evn='0000', species='co2', samplingstrategy=1, sdev=0.0, fromfile='none.nc'):
        self.code = code.strip()      # Site code
        self.xdate = xdate             # Date of obs
        self.obs = obs               # Value observed
        self.simulated = simulated         # Value simulated by model
        self.resid = resid             # Mole fraction residuals
        self.hphr = hphr              # Mole fraction prior uncertainty from fluxes and (HPH) and model data mismatch (R)
        self.mdm = mdm               # Model data mismatch
        self.may_localize = True           # Whether sample may be localized in optimizer
        self.may_reject = True              # Whether sample may be rejected if outside threshold
        self.flag = flag              # Flag
        self.height = height            # Sample height
        self.lat = lat               # Sample lat
        self.lon = lon               # Sample lon
        self.id = idx               # ID number
        self.evn = evn               # Event number
        self.sdev = sdev              # standard deviation of ensemble
        self.masl = True              # Sample is in Meters Above Sea Level
        self.mag = not self.masl     # Sample is in Meters Above Ground
        self.species = species.strip()
        self.samplingstrategy = samplingstrategy
        self.fromfile = fromfile          # netcdf filename inside observation distribution, to write back later

################### End Class MoleFractionSample ###################


if __name__ == "__main__":
    pass
