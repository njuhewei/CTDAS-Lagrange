#!/usr/bin/env python
# jobcontrol.py

"""
.. module:: platform
.. moduleauthor:: Wouter Peters 

Revision History:
File created on 06 Sep 2010.

The Platform class is found in the module :mod:`platform`, or in a specific implementation under the da/source tree. 

The platform object holds attributes and methods that allow job control on each specific platform. This includes methods to create and submit jobs, but also to obtain process and/or job ID's. These are needed to control the flow of 
the system on each platform.

Typically, every platform needs specific implementations of this object (through inheritance), and you should refer to your specific Platform object documentation for details (see *da/platform/*).

.. autoclass:: da.baseclasses.platform.Platform
   :members:
   :inherited-members:

"""

import os
import logging
import subprocess

std_joboptions = {'jobname':'test', 'jobaccount':'co2', 'jobnodes':'nserial 1', 'jobshell':'/bin/sh', 'depends':'', 'jobtime':'01:00:00'}

class Platform(object):
    """ 
    This specifies platform dependent options under generic object calls. A platform object is used to control and submit jobs
    """

    def __init__(self):
        """
        The init function reports the hard-coded ``Identifier`` and ``Version`` of the Platform. Since each new
        computer/user requires their own Platform object modifications, the init function is usually overwritten
        in the specific implementation of this class
        """
        self.ID = 'iPad'    # the identifier gives the plaform name
        self.version = '1.0'     # the platform version used

        logging.debug('%s object initialized' % self.ID)
        logging.debug('%s version: %s' % (self.ID, self.version))

    def give_blocking_flag(self):
        return ""

    def give_queue_type(self):
        return "foreground"

    def get_job_template(self, joboptions={}, block=False):
        """ 
        Returns the job template for a given computing system, and fill it with options from the dictionary provided as argument.
        The job template should return the preamble of a job that can be submitted to a queue on your platform, 
        examples of popular queuing systems are:
            - SGE
            - MOAB
            - XGrid
            -

        A list of job options can be passed through a dictionary, which are then filled in on the proper line,
        an example is for instance passing the dictionary {'account':'co2'} which will be placed 
        after the ``-A`` flag in a ``qsub`` environment.

        An extra option ``block`` has been added that allows the job template to be configured to block the current
        job until the submitted job in this template has been completed fully.
        """

        template = """## \n""" + \
                   """## This is a set of dummy names, to be replaced by values from the dictionary \n""" + \
                   """## Please make your own platform specific template with your own keys and place it in a subfolder of the da package.\n """ + \
                   """## \n""" + \
                   """ \n""" + \
                   """#$ jobname \n""" + \
                   """#$ jobaccount \n""" + \
                   """#$ jobnodes \n""" + \
                   """#$ jobtime \n""" + \
                   """#$ jobshell \n """

        if 'depends' in joboptions:
            template += """#$ -hold_jid depends \n"""

        # First replace from passed dictionary
        for k, v in joboptions.iteritems():
            while k in template:
                template = template.replace(k, v)
        # Fill remaining values with std_options
        for k, v in std_joboptions.iteritems():
            while k in template:
                template = template.replace(k, v)
        return template

    def get_my_id(self):
        """ Return the process ID, or job ID of the current process or job"""
        return os.getpid()

    def write_job(self, jobfile, template, jobid):                   
        """ 
        This method writes a jobfile to the exec dir and makes it executable (mod 477)
        """
        #
        # Done, write jobfile
        # 
        f = open(jobfile, 'w')
        f.write(template)
        f.close()
        os.chmod(jobfile, 477)                  
        logging.debug("A job file was created (%s)" % jobfile)

    def submit_job(self, jobfile, joblog=None, block=False): 
        """ 
           :param jobfile: a string with the filename of a jobfile to run
           :param joblog:  a string with the filename of a logfile to write run output to
           :param block:  Boolean specifying whether to submit and continue (F), or submit and wait (T)
           :rtype: integer

        This method submits a jobfile to the queue, and returns the job ID 
        """
        cmd = ["sh", jobfile]
        logging.info("A new task will be started (%s)" % cmd)
        if block:
            jobid = subprocess.call(cmd)
        else:
            jobid = subprocess.Popen(cmd).pid
        
        logging.info('Summary:')
        logging.info('job script      : %s' % jobfile)
        logging.info('job log         : %s' % joblog)
        logging.info('To manage this process:')
        logging.info('  # kill process:')
        logging.info('  kill %i\n' % jobid)


    def kill_job(self, jobid):                   
        """ This method kills a running job """

    def job_stat(self, jobid):                   
        """ This method gets the status of a running job """
        output = subprocess.Popen(['qstat', jobid], stdout=subprocess.PIPE).communicate()[0]  
        logging.info(output)
        return output


if __name__ == "__main__":
    pass
