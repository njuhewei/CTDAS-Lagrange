#!/usr/bin/env python
# ct_statevector_tools.py

"""
.. module:: statevector
.. moduleauthor:: Wouter Peters 

Revision History:
File created on 28 Jul 2010.

The module statevector implements the data structure and methods needed to work with state vectors (a set of unknown parameters to be optimized by a DA system) of different lengths, types, and configurations. Two baseclasses together form a generic framework:
    * :class:`~da.baseclasses.statevector.StateVector`
    * :class:`~da.baseclasses.statevector.EnsembleMember`

As usual, specific implementations of StateVector objects are done through inheritance form these baseclasses. An example of designing 
your own baseclass StateVector we refer to :ref:`tut_chapter5`.

.. autoclass:: da.baseclasses.statevector.StateVector 

.. autoclass:: da.baseclasses.statevector.EnsembleMember 

"""

import os
import logging
import numpy as np
from datetime import timedelta
import da.tools.io4 as io

identifier = 'Baseclass Statevector '
version = '0.0'

################### Begin Class EnsembleMember ###################

class EnsembleMember(object):
    """ 
        An ensemble member object consists of:
           * a member number
           * parameter values
           * an observation object to hold sampled values for this member

        Ensemble members are initialized by passing only an ensemble member number, all data is added by methods 
        from the :class:`~da.baseclasses.statevector.StateVector`. Ensemble member objects have almost no functionality 
        except to write their data to file using method :meth:`~da.baseclasses.statevector.EnsembleMember.write_to_file`

        .. automethod:: da.baseclasses.statevector.EnsembleMember.__init__ 
        .. automethod:: da.baseclasses.statevector.EnsembleMember.write_to_file 
        .. automethod:: da.baseclasses.statevector.EnsembleMember.AddCustomFields 

    """

    def __init__(self, membernumber):
        """
           :param memberno: integer ensemble number
           :rtype: None

           An EnsembleMember object is initialized with only a number, and holds two attributes as containter for later
           data:
                * param_values, will hold the actual values of the parameters for this data
                * ModelSample, will hold an :class:`~da.baseclasses.obs.Observation` object and the model samples resulting from this members' data

        """
        self.membernumber = membernumber   # the member number
        self.param_values = None           # Parameter values of this member

################### End Class EnsembleMember ###################

################### Begin Class StateVector ###################


class StateVector(object):
    """ 
    The StateVector object first of all contains the data structure of a statevector, defined by 3 attributes that define the 
    dimensions of the problem in parameter space:
        * nlag
        * nparameters
        * nmembers

    The fourth important dimension `nobs` is not related to the StateVector directly but is initialized to 0, and later on 
    modified to be used in other parts of the pipeline:
        * nobs

    These values are set as soon as the :meth:`~da.baseclasses.statevector.StateVector.setup` is called from the :ref:`pipeline`. 
    Additionally, the value of attribute `isOptimized` is set to `False` indicating that the StateVector holds a-priori values 
    and has not been modified by the :ref:`optimizer`.

    StateVector objects can be filled with data in two ways
        1. By reading the data from file
        2. By creating the data through a set of method calls

    Option (1) is invoked using method :meth:`~da.baseclasses.statevector.StateVector.read_from_file`. 
    Option (2) consists of a call to method :meth:`~da.baseclasses.statevector.StateVector.make_new_ensemble`

    Once the StateVector object has been filled with data, it is used in the pipeline and a few more methods are
    invoked from there:
        * :meth:`~da.baseclasses.statevector.StateVector.propagate`, to advance the StateVector from t=t to t=t+1
        * :meth:`~da.baseclasses.statevector.StateVector.write_to_file`, to write the StateVector to a NetCDF file for later use

    The methods are described below:

    .. automethod:: da.baseclasses.statevector.StateVector.setup 
    .. automethod:: da.baseclasses.statevector.StateVector.read_from_file
    .. automethod:: da.baseclasses.statevector.StateVector.write_to_file
    .. automethod:: da.baseclasses.statevector.StateVector.make_new_ensemble
    .. automethod:: da.baseclasses.statevector.StateVector.propagate
    .. automethod:: da.baseclasses.statevector.StateVector.write_members_to_file

    Finally, the StateVector can be mapped to a gridded array, or to a vector of TransCom regions, using:

    .. automethod:: da.baseclasses.statevector.StateVector.grid2vector
    .. automethod:: da.baseclasses.statevector.StateVector.vector2grid
    .. automethod:: da.baseclasses.statevector.StateVector.vector2tc
    .. automethod:: da.baseclasses.statevector.StateVector.state2tc

    """

    def __init__(self):
        self.ID = identifier
        self.version = version

        # The following code allows the object to be initialized with a dacycle object already present. Otherwise, it can
        # be added at a later moment.

        logging.info('Statevector object initialized: %s' % self.ID)

    def setup(self, dacycle):
        """
        setup the object by specifying the dimensions. 
        There are two major requirements for each statvector that you want to build:
        
            (1) is that the statevector can map itself onto a regular grid
            (2) is that the statevector can map itself (mean+covariance) onto TransCom regions

        An example is given below.
        """

        self.nlag = int(dacycle['time.nlag'])
        self.nmembers = int(dacycle['da.optimizer.nmembers'])
        self.nparams = int(dacycle.dasystem['nparameters'])
        self.nobs = 0
        
        self.obs_to_assimilate = ()  # empty containter to hold observations to assimilate later on

        # These list objects hold the data for each time step of lag in the system. Note that the ensembles for each time step consist 
        # of lists of EnsembleMember objects, we define member 0 as the mean of the distribution and n=1,...,nmembers as the spread.

        self.ensemble_members = range(self.nlag)

        for n in range(self.nlag):
            self.ensemble_members[n] = []


        # This specifies the file to read with the gridded mask at 1x1 degrees. Each gridbox holds a number that specifies the parametermember
        #  that maps onto it. From this map, a dictionary is created that allows a reverse look-up so that we can map parameters to a grid.

        mapfile = os.path.join(dacycle.dasystem['regionsfile'])
        ncf = io.ct_read(mapfile, 'read')
        self.gridmap = ncf.get_variable('regions')
        self.tcmap = ncf.get_variable('transcom_regions')
        ncf.close()

        logging.debug("A TransCom  map on 1x1 degree was read from file %s" % dacycle.dasystem['regionsfile'])
        logging.debug("A parameter map on 1x1 degree was read from file %s" % dacycle.dasystem['regionsfile'])

        # Create a dictionary for state <-> gridded map conversions

        nparams = self.gridmap.max()
        self.griddict = {}
        for r in range(1, int(nparams) + 1):
            sel = (self.gridmap.flat == r).nonzero()
            if len(sel[0]) > 0: 
                self.griddict[r] = sel

        logging.debug("A dictionary to map grids to states and vice versa was created")

        # Create a matrix for state <-> TransCom conversions

        self.tcmatrix = np.zeros((self.nparams, 23), 'float') 

        for r in range(1, self.nparams + 1):
            sel = (self.gridmap.flat == r).nonzero()
            if len(sel[0]) < 1: 
                continue
            else:
                n_tc = set(self.tcmap.flatten().take(sel[0]))
                if len(n_tc) > 1: 
                    logging.error("Parameter %d seems to map to multiple TransCom regions (%s), I do not know how to handle this" % (r, n_tc))
                    raise ValueError
                self.tcmatrix[r - 1, n_tc.pop() - 1] = 1.0

        logging.debug("A matrix to map states to TransCom regions and vice versa was created")

        # Create a mask for species/unknowns

        self.make_species_mask()

    def make_species_mask(self):

        """

        This method creates a dictionary with as key the name of a tracer, and as values an array of 0.0/1.0 values 
        specifying which StateVector elements are constrained by this tracer. This mask can be used in 
        the optimization to ensure that certain types of osbervations only update certain unknowns.

        An example would be that the tracer '14CO2' can be allowed to only map onto fossil fuel emissions in the state

        The form of the mask is:

        {'co2': np.ones(self.nparams), 'co2c14', np.zeros(self.nparams)  }

        so that 'co2' maps onto all parameters, and 'co2c14' on none at all. These arrays are used in the Class 
        optimizer when state updates are actually performed

        """
        self.speciesdict = {'co2': np.ones(self.nparams)}
        logging.debug("A species mask was created, only the following species are recognized in this system:")
        for k in self.speciesdict.keys(): 
            logging.debug("   ->    %s" % k)


    def make_new_ensemble(self, lag, covariancematrix=None):
        """ 
        :param lag: an integer indicating the time step in the lag order
        :param covariancematrix: a matrix to draw random values from
        :rtype: None
    
        Make a new ensemble, the attribute lag refers to the position in the state vector. 
        Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
        The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]

        The optional covariance object to be passed holds a matrix of dimensions [nparams, nparams] which is
        used to draw ensemblemembers from. If this argument is not passed it will ne substituted with an 
        identity matrix of the same dimensions.

        """    

        if covariancematrix == None: 
            covariancematrix = np.identity(self.nparams)

        # Make a cholesky decomposition of the covariance matrix


        try:
            _, s, _ = np.linalg.svd(covariancematrix)
        except:
            s = np.linalg.svd(covariancematrix, full_matrices=1, compute_uv=0) #Cartesius fix
        dof = np.sum(s) ** 2 / sum(s ** 2)
        C = np.linalg.cholesky(covariancematrix)

        logging.debug('Cholesky decomposition has succeeded ')
        logging.info('Appr. degrees of freedom in covariance matrix is %s' % (int(dof)))


        # Create mean values 

        newmean = np.ones(self.nparams, float) # standard value for a new time step is 1.0

        # If this is not the start of the filter, average previous two optimized steps into the mix

        if lag == self.nlag - 1 and self.nlag >= 3:
            newmean += self.ensemble_members[lag - 1][0].param_values + \
                                           self.ensemble_members[lag - 2][0].param_values 
            newmean = newmean / 3.0

        # Create the first ensemble member with a deviation of 0.0 and add to list

        newmember = EnsembleMember(0)
        newmember.param_values = newmean.flatten()  # no deviations
        self.ensemble_members[lag].append(newmember)

        # Create members 1:nmembers and add to ensemble_members list

        for member in range(1, self.nmembers):
            rands = np.random.randn(self.nparams)

            newmember = EnsembleMember(member)
            newmember.param_values = np.dot(C, rands) + newmean
            self.ensemble_members[lag].append(newmember)

        logging.debug('%d new ensemble members were added to the state vector # %d' % (self.nmembers, (lag + 1)))


    def propagate(self, dacycle):
        """ 
        :rtype: None

        Propagate the parameter values in the StateVector to the next cycle. This means a shift by one cycle 
        step for all states that will
        be optimized once more, and the creation of a new ensemble for the time step that just 
        comes in for the first time (step=nlag). 
        In the future, this routine can incorporate a formal propagation of the statevector.

        """
        
        # Remove State Vector n=1 by simply "popping" it from the list and appending a new empty list at the front. This empty list will
        # hold the new ensemble for the new cycle 

        self.ensemble_members.pop(0)
        self.ensemble_members.append([])

        # And now create a new time step of mean + members for n=nlag
        date = dacycle['time.start'] + timedelta(days=(self.nlag - 0.5) * int(dacycle['time.cycle']))
        cov = self.get_covariance(date, dacycle)
        self.make_new_ensemble(self.nlag - 1, cov)

        logging.info('The state vector has been propagated by one cycle')


    def write_to_file(self, filename, qual):
        """
        :param filename: the full filename for the output NetCDF file
        :rtype: None

        Write the StateVector information to a NetCDF file for later use. 
        In principle the output file will have only one two datasets inside 
        called:
            * `meanstate`, dimensions [nlag, nparamaters]
            * `ensemblestate`, dimensions [nlag,nmembers, nparameters]

        This NetCDF information can be read back into a StateVector object using 
        :meth:`~da.baseclasses.statevector.StateVector.read_from_file`

        """
        #import da.tools.io4 as io
        #import da.tools.io as io

        if qual == 'prior':
            f = io.CT_CDF(filename, method='create')
            logging.debug('Creating new StateVector output file (%s)' % filename)
            #qual = 'prior'
        else:
            f = io.CT_CDF(filename, method='write')
            logging.debug('Opening existing StateVector output file (%s)' % filename)
            #qual = 'opt'

        dimparams = f.add_params_dim(self.nparams)
        dimmembers = f.add_members_dim(self.nmembers)
        dimlag = f.add_lag_dim(self.nlag, unlimited=True)

        for n in range(self.nlag):
            members = self.ensemble_members[n]
            mean_state = members[0].param_values

            savedict = f.standard_var(varname='meanstate_%s' % qual)
            savedict['dims'] = dimlag + dimparams 
            savedict['values'] = mean_state
            savedict['count'] = n
            savedict['comment'] = 'this represents the mean of the ensemble'
            f.add_data(savedict)

            members = self.ensemble_members[n]
            devs = np.asarray([m.param_values.flatten() for m in members])
            data = devs - np.asarray(mean_state)

            savedict = f.standard_var(varname='ensemblestate_%s' % qual)
            savedict['dims'] = dimlag + dimmembers + dimparams 
            savedict['values'] = data
            savedict['count'] = n
            savedict['comment'] = 'this represents deviations from the mean of the ensemble'
            f.add_data(savedict)
        f.close()

        logging.info('Successfully wrote the State Vector to file (%s) ' % filename)

    def read_from_file(self, filename, qual='opt'):
        """ 
        :param filename: the full filename for the input NetCDF file
        :param qual: a string indicating whether to read the 'prior' or 'opt'(imized) StateVector from file
        :rtype: None

        Read the StateVector information from a NetCDF file and put in a StateVector object
        In principle the input file will have only one four datasets inside 
        called:
            * `meanstate_prior`, dimensions [nlag, nparamaters]
            * `ensemblestate_prior`, dimensions [nlag,nmembers, nparameters]
            * `meanstate_opt`, dimensions [nlag, nparamaters]
            * `ensemblestate_opt`, dimensions [nlag,nmembers, nparameters]

        This NetCDF information can be written to file using 
        :meth:`~da.baseclasses.statevector.StateVector.write_to_file`

        """

        #import da.tools.io as io
        f = io.ct_read(filename, 'read')
        meanstate = f.get_variable('statevectormean_' + qual)
        ensmembers = f.get_variable('statevectorensemble_' + qual)
        f.close()

        for n in range(self.nlag):
            if not self.ensemble_members[n] == []:
                self.ensemble_members[n] = []
                logging.warning('Existing ensemble for lag=%d was removed to make place for newly read data' % (n + 1))

            for m in range(self.nmembers):
                newmember = EnsembleMember(m)
                newmember.param_values = ensmembers[n, m, :].flatten() + meanstate[n]  # add the mean to the deviations to hold the full parameter values
                self.ensemble_members[n].append(newmember)

        logging.info('Successfully read the State Vector from file (%s) ' % filename)

    def write_members_to_file(self, lag, outdir,endswith='.nc'):
        """ 
           :param: lag: Which lag step of the filter to write, must lie in range [1,...,nlag]
           :param: outdir: Directory where to write files
           :param: endswith: Optional label to add to the filename, default is simply .nc
           :rtype: None

           Write ensemble member information to a NetCDF file for later use. The standard output filename is 
           *parameters.DDD.nc* where *DDD* is the number of the ensemble member. Standard output file location 
           is the `dir.input` of the dacycle object. In principle the output file will have only two datasets inside 
           called `parametervalues` which is of dimensions `nparameters` and `parametermap` which is of dimensions (180,360). 
           This dataset can be read and used by a :class:`~da.baseclasses.observationoperator.ObservationOperator` object. 

           .. note:: if more, or other information is needed to complete the sampling of the ObservationOperator you
                     can simply inherit from the StateVector baseclass and overwrite this write_members_to_file function.

        """

        # These import statements caused a crash in netCDF4 on MacOSX. No problems on Jet though. Solution was
        # to do the import already at the start of the module, not just in this method.
           
        #import da.tools.io as io
        #import da.tools.io4 as io

        members = self.ensemble_members[lag]

        for mem in members:
            filename = os.path.join(outdir, 'parameters.%03d%s' % (mem.membernumber, endswith))
            ncf = io.CT_CDF(filename, method='create')
            dimparams = ncf.add_params_dim(self.nparams)
            dimgrid = ncf.add_latlon_dim()

            data = mem.param_values

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametervalues"
            savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimparams 
            savedict['values'] = data
            savedict['comment'] = 'These are parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            griddata = self.vector2grid(vectordata=data)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap"
            savedict['long_name'] = "parametermap_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid 
            savedict['values'] = griddata.tolist()
            savedict['comment'] = 'These are gridded parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            ncf.close()

            logging.debug('Successfully wrote data from ensemble member %d to file (%s) ' % (mem.membernumber, filename))

    def grid2vector(self, griddata=None, method='avg'):
        """ 
            Map gridded data onto a vector of length (nparams,)

           :param griddata: a gridded dataset to use. This dataset is mapped onto a vector of length `nparams`
           :param method: a string that specifies the method to combine grid boxes in case reverse=True. Must be either ['avg','sum','minval']
           :rtype: ndarray: size (nparameters,)

           This method makes use of a dictionary that links every parameter number [1,...,nparams] to a series of gridindices. These 
           indices specify a location on a 360x180 array, stretched into a vector using `array.flat`. There are multiple ways of calling 
           this method::

               values       = self.grid2vector(griddata=mygriddeddata,method='minval') # 
                                                                    using the minimum value of all datapoints covered by that parameter index

               values       = self.grid2vector(griddata=mygriddeddata,method='avg') # 
                                                                    using the average value of all datapoints covered by that parameter index

               values       = self.grid2vector(griddata=mygriddeddata,method='sum') # 
                                                                    using the sum of values of all datapoints covered by that parameter index

           .. note:: This method uses a DaSystem object that must be initialized with a proper parameter map. See :class:`~da.baseclasses.dasystem` for details

        """

        methods = ['avg', 'sum', 'minval']
        if method not in methods:
            logging.error("To put data from a map into the statevector, please specify the method to use (%s)" % methods)
            raise ValueError

        result = np.zeros((self.nparams,), float)
        for k, v in self.griddict.iteritems():
            #print k,k-1,result.shape, v
            if method == "avg": 
                result[k - 1] = griddata.take(v).mean()
            elif method == "sum" : 
                result[k - 1] = griddata.take(v).sum()
            elif method == "minval" : 
                result[k - 1] = griddata.take(v).min()
        return result # Note that the result is returned, but not yet placed in the member.param_values attrtibute!


    def vector2grid(self, vectordata=None):
        """ 
            Map vector elements to a map or vice cersa

           :param vectordata: a vector dataset to use in case `reverse = False`. This dataset is mapped onto a 1x1 grid and must be of length `nparams`
           :rtype: ndarray: an array of size (360,180,) 

           This method makes use of a dictionary that links every parameter number [1,...,nparams] to a series of gridindices. These 
           indices specify a location on a 360x180 array, stretched into a vector using `array.flat`. There are multiple ways of calling 
           this method::

               griddedarray = self.vector2grid(vectordata=param_values) # simply puts the param_values onto a (180,360,) array

           .. note:: This method uses a DaSystem object that must be initialzied with a proper parameter map. See :class:`~da.baseclasses.dasystem` for details

        """
        result = np.zeros(self.gridmap.shape, float)
        for k, v in self.griddict.iteritems():
            #print k,v
            result.put(v, vectordata[k - 1])
        return result         

    def vector2tc(self, vectordata, cov=False):
        """ 
            project Vector onto TransCom regions 

           :param vectordata: a vector dataset to use, must be of length `nparams`
           :param cov: a Boolean to specify whether the input dataset is a vector (mean), or a matrix (covariance)
           :rtype: ndarray: an array of size (23,) (cov:F) or of size (23,23,) (cov:T)
        """

        M = self.tcmatrix
        if cov:
            return np.dot(np.transpose(M), np.dot(vectordata, M))
        else:
            return np.dot(vectordata.squeeze(), M)

    def state_to_grid(self, fluxvector=None, lag=1):
        """ 
            Transforms the StateVector information (mean + covariance) to a 1x1 degree grid.

            :param: fluxvector: a vector of length (nparams,) that holds the fluxes associated with each parameter in the StateVector
            :param: lag: the lag at which to evaluate the StateVector
            :rtype: a tuple of two arrays (gridmean,gridvariance) with dimensions (180,360,)

            If the attribute `fluxvector` is not passed, the function will return the mean parameter value and its variance on a 1x1 map.
            
            ..note:: Although we can return the variance information for each gridbox, the covariance information contained in the original ensemble is lost when mapping to 1x1 degree!

        """

        if fluxvector == None:
            fluxvector = np.ones(self.nparams)

        ensemble = self.ensemble_members[lag - 1]
        ensemblemean = ensemble[0].param_values

        # First transform the mean
        gridmean = self.vector2grid(vectordata=ensemblemean * fluxvector)

        # And now the covariance, first create covariance matrix (!), and then multiply
        deviations = np.array([mem.param_values * fluxvector - ensemblemean for mem in ensemble])
        ensemble = []
        for mem in deviations:
            ensemble.append(self.vector2grid(mem))

        return (gridmean, np.array(ensemble))

    def state2tc(self, fluxvector=None, lag=1):
        """ 
            Transforms the StateVector information (mean + covariance) to the TransCom regions.

            :param: fluxvector: a vector of length (nparams,) that holds the fluxes associated with each parameter in the StateVector
            :param: lag: the lag at which to evaluate the StateVector
            :rtype: a tuple of two arrays (mean,covariance) with dimensions ((23,), (23,23,) )

        """
        ensemble = self.ensemble_members[lag - 1]
        ensemblemean = ensemble[0].param_values

        # First transform the mean

        mean = self.vector2tc(vectordata=ensemble[0].param_values * fluxvector)

        # And now the covariance, first create covariance matrix (!), and then multiply

        deviations = np.array([mem.param_values * fluxvector - ensemblemean for mem in ensemble])
        covariance = np.dot(np.transpose(deviations), deviations) / (self.nmembers - 1)
        cov = self.vector2tc(covariance, cov=True)

        return (mean, cov)

    def get_covariance(self, date, cycleparams):
        pass
    
################### End Class StateVector ###################

if __name__ == "__main__":
    pass

