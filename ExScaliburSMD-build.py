#!/usr/bin/env python
"""The python wrapper script for creating a Tumor/Normal Paired Exome
   somatic mutation pipeline project.
"""

__author__='Kyle Hernandez'
__email__='khernandez@bsd.uchicago.edu'
__version__='0.1.0'
__license__='LGPLv3'
__url__='https://bitbucket.org/kmhernan/bds-exscalibursmd'

import argparse
import datetime
import os
import logging
import sys
import time

from project.project import Project

def build_project():
    '''Builds the project directory, YAML files, and job files.'''
    # Check if the config file exists
    assert os.access(args.config_file, os.R_OK), \
                     logger.error("Could not find config file '{0}'. Please provide the correct path".format(
                                   args.config_file))

    # Get global application path
    PATH = os.path.dirname(os.path.realpath(__file__))
    project = Project(args, PATH)
    project.initialize_project()
    project.configure_pipeline()

def get_parser():
    '''Set up the main parsing commands.'''
    p = argparse.ArgumentParser(prog='ExScaliburSMD-build',
            epilog="Copyright 2014 - Center for Research Informatics - University of Chicago")

    p_main = p.add_argument_group(title='Main Commands')
    p_main.add_argument('-m', '--metadata',
        help='Tab-delimited file with all sample information.', required=True)
    p_main.add_argument('-o', '--output-directory', required=True,
        help='Top-level directory for all project outputs.')
    p_main.add_argument('-c', '--config-file', required=True,
        help='The configuration file containing relevent options for the pipeline.')
    p_main.add_argument('-t', '--threads', type=int, default=1,
        help='Number of threads to use for all multi-threaded tools; ' + \
             'This can be overridden in the system config. [1]')
    p_main.add_argument('-p', '--project-id', required=True,
        help='The ID for this SMD project')
    p_main.add_argument('--target-bed',
        help='A file of target intervals if desired to limit pipeline to specific regions. ' +
             'The file must be in the bed format.')
    p_main.add_argument('--max-splits', type=int, default=50, 
        help='The maximum number of ways to split the callable regions for parallel processing. ' + \
             'This does not apply to Virmid or Strelka. [50]')

    p_clip = p.add_argument_group(title='Fastq Clipping Options')
    p_clip.add_argument('--no-adapter-clipping', action='store_true',
        help="Don't clip adapters from reads or, if paired, merge overlapping pairs via SeqPrep")
    p_clip.add_argument('--fastq-min-length', type=int, default = 30,
        help="Minimum length of a trimmed or merged read for clipping. " + \
             "If using BWA mem, should be >= 70. [30]")

    p_align = p.add_argument_group(title='Alignment Options')
    p_align.add_argument('--bwa-aln', action='store_true',
        help='Run BWA aln. If no aligner is specified, this will be the one and only alignment tool used.')
    p_align.add_argument('--bwa-mem', action='store_true',
        help='Run BWA mem. This algorithm suggests that reads are >= 70 bp.')
    p_align.add_argument('--novoalign', action='store_true',
        help='Run novoalign alignment software.')

    p_smd = p.add_argument_group(title='Somatic Mutation Detection Options')
    p_smd.add_argument('--mutect', action='store_true',
        help='Run MuTect. If no other SMD callers are specified, this will be the one and ' + \
             'only SMD caller used.')
    p_smd.add_argument('--varscan', action='store_true', help='Run VarScan2.')
    p_smd.add_argument('--somatic-sniper', action='store_true', help='Run SomaticSniper')
    p_smd.add_argument('--strelka', action='store_true', help='Run Strelka')
    p_smd.add_argument('--virmid', action='store_true', help='Run Virmid')
    p_smd.add_argument('--shimmer', action='store_true', help='Run SHIMMER')
    p_smd.add_argument('--min-map-q', type=int, default=20,
        help='Minimum MAPQ score to consider for variant calling. [20]')
    return p

if __name__ == '__main__':
    start = time.time()

    # Set up logger
    logger    = logging.getLogger('SMD')
    logger.setLevel(logging.INFO)
    ch        = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(levelname)s] [%(asctime)s] [%(name)s] - %(message)s',
                                  datefmt='%H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Print header
    logger.info('-'*75)
    logger.info('ExScaliburSMD-build v{0}'.format(__version__))
    logger.info('Copyright (c) 2014, {0}'.format(__license__))
    logger.info('The Center for Research Informatics, The University of Chicago')
    logger.info('For support and documentation go to:')
    logger.info('{0}'.format(__url__))
    logger.info("Program Args: ExScaliburSMD-build.py " + " ".join(sys.argv[1::]))
    logger.info('Date/time: {0}'.format(datetime.datetime.now()))
    logger.info('-'*75)
    logger.info('-'*75)

    # Set up parser
    parser = get_parser()
    args   = parser.parse_args()

    # Send to build_project
    build_project()

    # Done
    logger.info("Finished, took {0} seconds.".format(time.time() - start))
