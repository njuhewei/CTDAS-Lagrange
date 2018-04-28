#!/usr/bin/env python
# io.py

"""
Author : peters 

Revision History:
File created on 15 Oct 2008.
File modified for CT data assimilation system in July 2010, Wouter Peters

"""
import standardvariables
import netCDF4
#import pyhdf.SD as hdf
import datetime as dt
from numpy import array, arange
import os
import logging
import sys

disclaimer = "This data belongs to the CarbonTracker project"
email = "wouter.peters@wur.nl"
url = "http://carbontracker.wur.nl"
institution = "Wageningen University and Research Center"
source 	= "CarbonTracker release 2.0" 
conventions = "CF-1.1"
historytext	= 'created on '+dt.datetime.now().strftime('%B %d, %Y')+' by %s'%os.environ['USER']

std_savedict={'name':'unknown','values':[],'dims':(0,0,),'units':'','long_name':'','comment':''}

def ct_read(filename='',method=''):
    """ read from an HDF or NetCDF file. Function choses itself which type is needed """

    if 'hdf' in filename.split('.'):
        return CT_HDF(filename,method)
    elif 'nc' in filename.split('.'):
        return CT_CDF(filename,method)
    else:
        msg = 'Could not determine whether input file was NetCDF or HDF trying both: ' ; logging.warning(msg)
        try:
            return CT_CDF(filename,method)
        except:
            return CT_HDF(filename,method)

class CT_CDF(netCDF4.Dataset):
    """ function opens a NetCDF file for writing of output"""

    def __init__(self,filename, method='read'):

        if method not in ['read','write','create']:
            raise ValueError, 'Method %s is not defined for a CarbonTracker NetCDF file object' % method

        if method == 'read':
            try:
                super(CT_CDF,self).__init__(filename, 'r') 
            except RuntimeError: 
                msg = 'Requested file not found for opening: %s'%filename ; logging.error(msg)
                msg = "Exiting" ; logging.info(msg)
                sys.exit(2)
        elif method == 'write':
            try:
                super(CT_CDF,self).__init__(filename, 'a')
            except:
                super(CT_CDF,self).__init__(filename, 'w',format='NETCDF4') 

            #self.AddCTHeader()
        elif method == 'create':
            if os.path.exists(filename): os.remove(filename)
            super(CT_CDF,self).__init__(filename, 'w',format='NETCDF4') 
            self.add_tc_header()


    def add_tc_header(self):

        #
        self.setncattr('Institution',institution)
        self.setncattr('Contact',email)
        self.setncattr('URL',url)
        self.setncattr('Source',source)
        self.setncattr('Convention',conventions)
        self.setncattr('Disclaimer',disclaimer)
        self.setncattr('History',historytext)

    def add_params_dim(self,nparams):

        if 'nparameters' in self.dimensions.keys():
            pass
        else:
            dimparams=self.createDimension('nparameters',size=nparams)

        return ('nparameters',)

    def add_members_dim(self,nmembers):

        if 'nmembers' in self.dimensions.keys():
            pass
        else:
            dimmembers=self.createDimension('nmembers',size=nmembers)

        return ('nmembers',)

    def add_lag_dim(self,nlag,unlimited=True):

        if 'nlag' in self.dimensions.keys():
            pass
        else:
            if unlimited:
                dimlag = self.createDimension('nlag',size=None)
            else:
                dimlag = self.createDimension('nlag',size=nlag)

        return ('nlag',)

    def add_obs_dim(self,nobs):

        if 'nobs' in self.dimensions.keys():
            pass
        else:
            dimobs = self.createDimension('nobs',size=nobs)

        return ('nobs',)

    def add_latlon_dim(self,istart=0,iend=360,jstart=0,jend=180):

        from numpy import arange, float64

        if 'latitude'  in self.dimensions.keys(): return ('latitude','longitude',)

        lons=-180+arange(360)*1.0+0.5
        lats=-90+arange(180)*1.0+0.5
        #
        lats=lats[jstart:jend]
        lons=lons[istart:iend]
        #
        dimlon = self.createDimension('longitude',size=lons.shape[0])
        dimlat = self.createDimension('latitude',size=lats.shape[0])

        savedict=self.standard_var(varname='latitude')
        savedict['values']=lats.tolist()
        savedict['actual_range']=(float(lats[0]),float(lats[-1]))
        savedict['dims']=('latitude',)
        self.add_data(savedict)

        savedict=self.standard_var(varname='longitude')
        savedict['values']=lons.tolist()
        savedict['actual_range']=(float(lons[0]),float(lons[-1]))
        savedict['dims']=('longitude',)
        self.add_data(savedict)

        return ('latitude','longitude',)

    def add_region_dim(self,type='eco',dimsize=None):

        from da.analysis.tools_transcom import olsonlabs, transnams, ext_transnams, ext_transcomps, olsonnams
        from da.analysis.tools_regions import ext_econams, ext_ecocomps

        if type not in ['eco','eco_ext','tc','tc_ext','olson']: 
            raise ValueError,'Type of dimension for regions requested (%s) is not possible' %type

        dimname='regions_%s' % type

        if dimname in self.dimensions.keys(): 
            return (dimname,)

        if type == 'olson':

            dim = self.createDimension(dimname,size=len(olsonlabs))

            for i,name in enumerate(olsonnams):
                att = setattr(self, 'OlsonEcosystem_%03d'%(i+1,), name  )

        elif type == 'tc':

            dim = self.createDimension(dimname,size=len(transnams))
            for i,name in enumerate(transnams):
                att = setattr(self, 'TransComRegion_%03d'%(i+1,), name  )

        elif type == 'tc_ext':

            dim = self.createDimension(dimname,size=len(ext_transnams))

            for i,name in enumerate(ext_transnams):
                    lab='Aggregate_Region_%03d'%(i+1,)
                    setattr(self,lab,name)
            for i,name in enumerate(ext_transcomps):
                    lab='Aggregate_Components_%03d'%(i+1,)
                    setattr(self,lab,name)

        elif type == 'eco':

            dim = self.createDimension(dimname,size=dimsize)

        return (dimname,)


    def add_date_dim(self,unlimited=False):

        if 'date' in self.dimensions.keys():
            pass
        else:
            dimdate = self.createDimension('date',size=None)

        return ('date',)

    def add_date_dim_format(self):

        if 'yyyymmddhhmmss'  in self.dimensions.keys(): 
            pass
        else:
            dimdateformat = self.createDimension('yyyymmddhhmmss',size=6)
            return ('yyyyymmddhhmmss',)

    def add_dim(self,dimname,dimsize):

        if dimname  in self.dimensions.keys(): 
            pass
        else:
            newdim = self.createDimension(dimname,dimsize)
        return (dimname,)

    def has_date(self,dd):

        if not self.dimensions.has_key('date'): 
            return False
        if not self.variables.has_key('date'): 
            return False
        if self.dimensions['date'].isunlimited:
            if dd in self.get_variable('date').tolist():
                return True
            else:
                return False
        else:
            return False
            
    def get_variable(self,varname):
        """ get variable from ncf file"""
        return self.variables[varname][:]

    def get_attribute(self,attname):
        """ get attribute from ncf file"""
        return getattr(self,attname)

    def add_attribute(self,attname,attvalue):
        """ set attribute in ncf file"""
        self.setncattr(attname,attvalue)

    def standard_var(self,varname):
        """ return properties of standard variables """
        import standardvariables

        if varname in standardvariables.standard_variables.keys():
            return standardvariables.standard_variables[varname]
        else:
            return standardvariables.standard_variables['unknown']

    def inq_unlimlen(self):
        """ return lenght of unlimited dimenion(s) """

        unlims=()
        for dimname, dimobj in self.dimensions.iteritems():
            if dimobj.isunlimited() : unlims += (len(dimobj),)

        return unlims

    def has_unlimlen(self,dims):
        """ return T/F whether dimensions include an unlimited dimenion(s) """

        for dimname, dimobj in self.dimensions.iteritems():
            if dimname in dims:
                if dimobj.isunlimited() : return True

        return False

    def add_variable(self,datadict,silent=True):
        """ add variables to file, but no data"""
        import numpy as np

        existing_vars=self.variables

        if existing_vars.has_key(datadict['name']):
            return
        else:
            if not silent: print 'Creating new dataset: '+datadict['name']

            if datadict.has_key('dtype'):
                if datadict['dtype'] == 'int':
                    var = self.createVariable(datadict['name'],'i4',datadict['dims'])
                elif datadict['dtype'] == 'int64':
                    var = self.createVariable(datadict['name'],'i8',datadict['dims'])
                elif datadict['dtype'] == 'char':
                    var = self.createVariable(datadict['name'],'S1',datadict['dims'],fill_value='!')
                elif datadict['dtype'] == 'float':
                    var = self.createVariable(datadict['name'],'f4',datadict['dims'])
                elif datadict['dtype'] == 'double':
                    var = self.createVariable(datadict['name'],'f8',datadict['dims'])
                else:
                    var = self.createVariable(datadict['name'],'f8',datadict['dims'])
            else:
                var = self.createVariable(datadict['name'],'f4',datadict['dims'])

            for k,v in datadict.iteritems(): 
                if k not in ['name','dims','values','_FillValue','count']: 
                    var.setncattr(k,v)


    def add_data(self,datadict,nsets=1,silent=True):
        """ add fields to file, at end of unlimited dimension"""
        import numpy as np

        existing_vars=self.variables

        try: 
            next = datadict['count']
        except:
            next=0


        if existing_vars.has_key(datadict['name']):
            var = self.variables[datadict['name']] 
            ndims = var.ndim

            datadict = ConvertCharDims(var,datadict)

            if ndims == 1:
                var[next:next+nsets]=datadict['values']
            elif ndims == 2:
                var[next:next+nsets,:]=datadict['values']
            elif ndims == 3:
                var[next:next+nsets,:,:]=datadict['values']
            elif ndims == 4:
                var[next:next+nsets,:,:,:]=datadict['values']
            elif ndims == 5:
                var[next:next+nsets,:,:,:,:]=datadict['values']
            else:
                print 'More than 5 dimensions in array not implemented yet'
                raise ValueError

        else:
            if not silent: print 'Creating new dataset: '+datadict['name']

            if datadict.has_key('dtype'):
                if datadict['dtype'] == 'int':
                    var = self.createVariable(datadict['name'],'i4',datadict['dims'])#,fill_value=datadict['_FillValue'])
                elif datadict['dtype'] == 'int64':
                    var = self.createVariable(datadict['name'],'i8',datadict['dims'])#,fill_value=datadict['_FillValue'])
                elif datadict['dtype'] == 'char':
                    var = self.createVariable(datadict['name'],'S1',datadict['dims'],fill_value='!')
                elif datadict['dtype'] == 'float':
                    var = self.createVariable(datadict['name'],'f4',datadict['dims'])#,fill_value=datadict['_FillValue'])
                elif datadict['dtype'] == 'double':
                    var = self.createVariable(datadict['name'],'f8',datadict['dims'])#,fill_value=datadict['_FillValue'])
                else:
                    var = self.createVariable(datadict['name'],'f8',datadict['dims'])#,fill_value=datadict['_FillValue'])
            else:
                var = self.createVariable(datadict['name'],'f4',datadict['dims'])#,fill_value=datadict['_FillValue'])

            for k,v in datadict.iteritems(): 
                if k not in ['name','dims','values','_FillValue','count']: 
                    var.setncattr(k,v)

            #if nsets > 1 or self.has_unlimlen(datadict['dims']) == True:
            if nsets > 1 or (nsets > 0 and self.has_unlimlen(datadict['dims']) ) == True:
                ndims = var.ndim

                datadict = ConvertCharDims(var,datadict)
                if ndims == 1:
                    var[next:next+nsets]=datadict['values']
                elif ndims == 2:
                    var[next:next+nsets,:]=datadict['values']
                elif ndims == 3:
                    var[next:next+nsets,:,:]=datadict['values']
                elif ndims == 4:
                    var[next:next+nsets,:,:,:]=datadict['values']
                elif ndims == 5:
                    var[next:next+nsets,:,:,:,:]=datadict['values']
                else:
                    print 'More than 5 dimensions in array not implemented yet'
                    raise ValueError
            else:
                ndims = var.ndim

                datadict = ConvertCharDims(var,datadict)

                var[:] = datadict['values']

try:
	import pyhdf.SD as hdf
	class CT_HDF(hdf.SD):
	    """ function opens a HDF file for reading """

	    def __init__(self,filename, method='read'):

		if method in ['write','create']:
		    raise ValueError, 'Method %s is not defined for a CarbonTracker HDF file object' % method

		if method == 'read':
		    #print 'Reading from file'
		    try:
			super(CT_HDF,self).__init__(filename) 
		    except hdf.HDF4Error: 
			msg = 'Requested file not found for opening: %s'%filename ; logging.error(msg)
			msg = "Exiting" ; logging.info(msg)
			sys.exit(2)

	    def get_variable(self,varname):
		""" get variable from ncf file"""
		return self.select(varname).get()

	    def get_attribute(self,attname):
		""" get attribute from ncf file"""
		return getattr(self,attname)

	    def standard_var(self,varname):
		""" return properties of standard variables """
		import standardvariables

		if varname in standardvariables.standard_variables.keys():
		    return standardvariables.standard_variables[varname]
		else:
		    return standardvariables.standard_variables['unknown']

	    def close(self):
		""" close file"""

		return self.end()
except:
	print('IO Class CT_HDF not compiled, no HDF support!!!')




def ConvertCharDims(var,datadict):

    if not var.dtype == 'S1': 
        pass
    else:
        datalen = len(datadict['values'])
        dimlen = list(var.shape)

        dimlen.remove(datalen) # string length remaining 

        slen=dimlen[0]

        #print [d for d in datadict['values'] ]
        values = [netCDF4.stringtoarr(d,slen) for d in datadict['values'] ] 
        datadict['values']  = values

    return datadict

def get_variable(file,varname):
    """ get variable from HDF file"""
    return array(file.select(varname).get())


if __name__ == '__main__':
    import numpy as np

    ncf=CT_CDF('test.nc','create')
    print ncf.file_format
    dimmembers=ncf.add_members_dim(200)
    dimparams=ncf.add_params_dim(200)

    dimdate=ncf.add_date_dim()
    dimidate=ncf.add_date_dim_format()
    dimlon,dimlat=ncf.add_latlon_dim()

    savedict=std_savedict.copy()
    savedict['name']='testvar'
    savedict['values']=np.zeros((2,200,))+2.0
    savedict['dims']=('date','nparameters',)
    ncf.add_data(savedict,nsets=2)

    savedict=std_savedict.copy()
    savedict['name']='testvar'
    savedict['values']=np.zeros((3,200,))+3.0
    savedict['dims']=('date','nparameters',)
    savedict['count']=2
    ncf.add_data(savedict,nsets=3)

    ncf.close()


