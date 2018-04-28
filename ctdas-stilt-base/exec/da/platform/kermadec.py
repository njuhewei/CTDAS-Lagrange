#!/usr/bin/env python
# kermadec.py

"""
Author : peters 

Revision History:
File created on 05 Jan 2015.

"""

import logging

from da.baseclasses.platform import Platform, std_joboptions

class KermadecPlatform(Platform):
    def __init__(self):

        self.ID = 'RUG kermadec'    # the identifier gives the platform name
        self.version = '1.0'     # the platform version used

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

        template = """#!/bin/sh \n""" + \
                   """## \n""" + \
                   """## This is a set of dummy names, to be replaced by values from the dictionary \n""" + \
                   """## Please make your own platform specific template with your own keys and place it in a subfolder of the da package.\n """ + \
                   """## \n""" + \
                   """ \n""" + \
                   """#$ jobname \n""" + \
                   """#$ jobaccount \n""" + \
                   """#$ jobnodes \n""" + \
                   """#$ jobtime \n""" + \
                   """#$ jobshell \n""" + \
                   """\n""" + \
                   """export HOST='kermadec'\n""" + \
                   """\n"""

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



if __name__ == "__main__":
    pass
