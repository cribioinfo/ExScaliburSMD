#!/usr/bin/env python
"""The python wrapper script for filtering / converting SMD output files 
   from the ExScalibur somatic mutation pipeline.
"""

__author__='Kyle Hernandez'
__email__='khernandez@bsd.uchicago.edu'
__version__='0.5.0'
__license__='LGPLv3'
__url__='https://bitbucket.org/cribioinformatics/exscalibursmd'

import argparse
import datetime
import os
import logging
import sys
import time

from filters import mutect, shimmer, sniper, strelka, varscan, virmid
 
def get_parser():
    '''Set up the main parsing commands.'''
    p = argparse.ArgumentParser(prog='ExScaliburSMD-filter',
            epilog="Copyright 2014 - Center for Research Informatics - University of Chicago")

    subparsers = p.add_subparsers(help='Choose the Somatic Mutation Detection tool that was used',
                                  dest='choice')

    ### Mutect 
    p_mutect = subparsers.add_parser('mutect', help='Options for processing MuTect output files.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_mutect.add_argument('--min-normal-depth', type=int, default=8,
        help='Minimum acceptable read depth for normal sample')
    p_mutect.add_argument('--min-tumor-depth', type=int, default=8,
        help='Minimum acceptable read depth for tumor sample')
    p_mutect.add_argument('--max-alt-freq-normal', type=float, default=0.05,
        help='Maximum alternate allele frequency allowed in the normal sample')
    p_mutect.add_argument('--min-alt-freq-tumor', type=float, default=0.20,
        help='Minimum alternate allele frequency allowed in the tumor sample')
    p_mutect.add_argument('--min-base-quality', type=int, default=30,
        help='Minimum average base quality for reads supporting alleles for any sample')
    p_mutect.add_argument('--tumor-name', required=True, help='The name of the tumor sample') 
    p_mutect.add_argument('--normal-name', required=True, help='The name of the normal sample') 
    p_mutect.add_argument('input_vcf', help='Input VCF file from MuTect to process')
    p_mutect.add_argument('output_vcf', help='Output filtered VCF file')

    ### Shimmer 
    p_shimmer = subparsers.add_parser('shimmer', help='Options for processing Shimmer output files.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_shimmer.add_argument('--min-normal-depth', type=int, default=8,
        help='Minimum acceptable read depth for normal sample')
    p_shimmer.add_argument('--min-tumor-depth', type=int, default=8,
        help='Minimum acceptable read depth for tumor sample')
    p_shimmer.add_argument('--max-alt-freq-normal', type=float, default=0.05,
        help='Maximum alternate allele frequency allowed in the normal sample')
    p_shimmer.add_argument('--min-alt-freq-tumor', type=float, default=0.20,
        help='Minimum alternate allele frequency allowed in the tumor sample')
    p_shimmer.add_argument('--min-qual', type=int, default=30,
        help='Minimum QUAL score')
    p_shimmer.add_argument('--tumor-name', required=True, help='The name of the tumor sample') 
    p_shimmer.add_argument('--normal-name', required=True, help='The name of the normal sample') 
    p_shimmer.add_argument('--reference', required=True, help='The reference file used')
    p_shimmer.add_argument('input_vcf', help='Input VCF file from Shimmer to process')
    p_shimmer.add_argument('input_vs', help='Input varsifter file from Shimmer to process')
    p_shimmer.add_argument('output_vcf', help='Output filtered VCF file')

    ## SomaticSniper
    p_sniper = subparsers.add_parser('sniper', help='Options for processing SomaticSniper output files.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_sniper.add_argument('--min-normal-depth', type=int, default=8,
        help='Minimum acceptable read depth for normal sample')
    p_sniper.add_argument('--min-tumor-depth', type=int, default=8,
        help='Minimum acceptable read depth for tumor sample')
    p_sniper.add_argument('--min-mapq-normal', type=int, default=20,
        help='Minimum mapping quality for normal sample')
    p_sniper.add_argument('--min-mapq-tumor', type=int, default=30,
        help='Minimum mapping quality for tumor sample')
    p_sniper.add_argument('--min-gq-normal', type=int, default=20,
        help='Minimum genotype quality for normal sample')
    p_sniper.add_argument('--min-gq-tumor', type=int, default=20,
        help='Minimum genotype quality for tumor sample')
    p_sniper.add_argument('--min-somatic-score', type=int, default=40,
        help='Minimum somatic score')
    p_sniper.add_argument('--tumor-name', required=True, help='The name of the tumor sample')
    p_sniper.add_argument('--normal-name', required=True, help='The name of the normal sample')
    p_sniper.add_argument('--reference', required=True, help='The reference file used')
    p_sniper.add_argument('input_vcf', help='Input VCF file from SomaticSniper to process')
    p_sniper.add_argument('output_vcf', help='Output filtered VCF file')

    ## Strelka 
    p_strelka = subparsers.add_parser('strelka', help='Options for processing Strelka output files.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_strelka.add_argument('--min-normal-depth', type=int, default=8,
        help='Minimum acceptable read depth for normal sample')
    p_strelka.add_argument('--min-tumor-depth', type=int, default=8,
        help='Minimum acceptable read depth for tumor sample')
    p_strelka.add_argument('--tumor-name', required=True, help='The name of the tumor sample')
    p_strelka.add_argument('--normal-name', required=True, help='The name of the normal sample')
    p_strelka.add_argument('input_all_snv', help='Input all snv VCF calls file from Strelka to process')
    p_strelka.add_argument('input_all_indel', help='Input all InDel VCF calls file from Strelka to process')
    p_strelka.add_argument('output_vcf', help='Output filtered VCF file')

    ## VarScan2 
    p_varscan = subparsers.add_parser('varscan', help='Options for processing VarScan2 output files.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_varscan.add_argument('--min-normal-depth', type=int, default=8,
        help='Minimum acceptable read depth for normal sample')
    p_varscan.add_argument('--min-tumor-depth', type=int, default=8,
        help='Minimum acceptable read depth for tumor sample')
    p_varscan.add_argument('--max-alt-freq-normal', type=float, default=0.05,
        help='Maximum alternate allele frequency allowed in the normal sample')
    p_varscan.add_argument('--min-alt-freq-tumor', type=float, default=0.20,
        help='Minimum alternate allele frequency allowed in the tumor sample')
    p_varscan.add_argument('--pval-cutoff', type=float, default=0.05,
        help='Somatic Pvalue cutoff')
    p_varscan.add_argument('--tumor-name', required=True, help='The name of the tumor sample')
    p_varscan.add_argument('--normal-name', required=True, help='The name of the normal sample')
    p_varscan.add_argument('--reference', required=True, help='The reference file used')
    p_varscan.add_argument('input_snv', help='Input snv file from VarScan2 to process')
    p_varscan.add_argument('input_indel', help='Input InDel file from VarScan2 to process')
    p_varscan.add_argument('output_vcf', help='Output filtered VCF file')

    ## Virmid 
    p_virmid = subparsers.add_parser('virmid', help='Options for processing Virmid output files.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_virmid.add_argument('--min-normal-depth', type=int, default=8,
        help='Minimum acceptable read depth for normal sample')
    p_virmid.add_argument('--min-tumor-depth', type=int, default=8,
        help='Minimum acceptable read depth for tumor sample')
    p_virmid.add_argument('--max-alt-freq-normal', type=float, default=0.05,
        help='Maximum alternate allele frequency allowed in the normal sample')
    p_virmid.add_argument('--min-alt-freq-tumor', type=float, default=0.20,
        help='Minimum alternate allele frequency allowed in the tumor sample')
    p_virmid.add_argument('--min-qual', type=int, default=20,
        help='QUAL score cutoff')
    p_virmid.add_argument('--tumor-name', required=True, help='The name of the tumor sample')
    p_virmid.add_argument('--normal-name', required=True, help='The name of the normal sample')
    p_virmid.add_argument('--reference', required=True, help='The reference file used')
    p_virmid.add_argument('input_all_somatic', help='Input all somatic VCF calls file from Virmid to process')
    p_virmid.add_argument('input_all_germline', help='Input all germline VCF calls file from Virmid to process') 
    p_virmid.add_argument('input_all_loh', 
                          help='Input all loss-of-heterozygosity (loh) VCF calls file from Virmid to process') 
    p_virmid.add_argument('output_vcf', help='Output filtered VCF file')
    return p

if __name__ == '__main__':
    start = time.time()

    # Set up logger
    logger    = logging.getLogger('SomaticFilter')
    logger.setLevel(logging.INFO)
    ch        = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(levelname)s] [%(asctime)s] [%(name)s] - %(message)s', 
                                  datefmt='%H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Print header
    logger.info('-'*75)
    logger.info('ExScaliburSMD-filter.py v{0}'.format(__version__))
    logger.info('Part of the ExScalibur-SMD Pipeline')
    logger.info('Copyright (c) 2014, {0}'.format(__license__))
    logger.info('The Center for Research Informatics, The University of Chicago')
    logger.info('For support and documentation go to:')
    logger.info('{0}'.format(__url__))
    logger.info("Program Args: ExScaliburSMD-filter.py " + " ".join(sys.argv[1::]))
    logger.info('Date/time: {0}'.format(datetime.datetime.now()))
    logger.info('-'*75)
    logger.info('-'*75)

    # Set up parser
    parser = get_parser()
    args   = parser.parse_args()

    # Send to appropriate function 
    if args.choice == 'mutect':
        mutect.run(args)
    elif args.choice == 'shimmer':
        shimmer.run(args)
    elif args.choice == 'sniper':
        sniper.run(args)
    elif args.choice == 'strelka':
        strelka.run(args)
    elif args.choice == 'varscan':
        varscan.run(args)
    elif args.choice == 'virmid':
        virmid.run(args)

    # Done
    logger.info("Finished, took {0} seconds.".format(time.time() - start))
