#!/usr/bin/env python

#################################################################################################
# First order of business is always to make all other python modules accessible through the path
#################################################################################################

import sys
import os
import logging
sys.path.append(os.getcwd())

#################################################################################################
# Next, import the tools needed to initialize a data assimilation cycle
#################################################################################################

from da.tools.initexit import start_logger, validate_opts_args, parse_options, CycleControl
from da.stilt.pipeline import ensemble_smoother_pipeline, forward_pipeline,header, footer, analysis_pipeline, archive_pipeline
from da.platform.kermadec import KermadecPlatform
from da.carbondioxide.dasystem import CO2DaSystem
from da.stilt.optimizer import CO2Optimizer
from da.stilt.obspack import ObsPackObservations
from da.stilt.statevector import CO2GriddedStateVector
from da.stilt.observationoperator import STILTObservationOperator


#################################################################################################
# Parse and validate the command line options, start logging
#################################################################################################

start_logger()
opts, args = parse_options()
opts, args = validate_opts_args(opts, args)

#################################################################################################
# Create the Cycle Control object for this job
#################################################################################################

dacycle = CycleControl(opts, args)

###########################################################################################
### IMPORT THE APPLICATION SPECIFIC MODULES HERE, TO BE PASSED INTO THE MAIN PIPELINE!!! ##
###########################################################################################


platform = KermadecPlatform()
dasystem = CO2DaSystem(dacycle['da.system.rc'])
obsoperator = STILTObservationOperator()
samples = ObsPackObservations()
statevector = CO2GriddedStateVector()
optimizer = CO2Optimizer()
##########################################################################################
################### ENTER THE PIPELINE WITH THE OBJECTS PASSED BY THE USER ###############
##########################################################################################



logging.info(header + "Entering Pipeline " + footer)

ensemble_smoother_pipeline(dacycle, platform, dasystem, samples,statevector, obsoperator, optimizer)
#forward_pipeline(dacycle, platform, dasystem, samples, statevector, obsoperator)



##########################################################################################
################### All done, extra stuff can be added next, such as analysis
##########################################################################################

analysis_pipeline(dacycle, platform, dasystem, samples, statevector )

archive_pipeline(dacycle, platform, dasystem)

sys.exit(0)


