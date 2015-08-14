#!/usr/bin/env python
"""Python script for creating the aggregated and summarized files for the
   ExScaliburViz reporter tool. 
"""

__author__='Kyle Hernandez'
__email__='khernandez@bsd.uchicago.edu'
__version__='0.5.0'
__license__='LGPLv3'
__url__='https://bitbucket.org/kmhernan/bds-exscalibursmd'

import argparse
import datetime
from operator import itemgetter
import os
import logging
import sys
import time
import yaml

def check_dir(path):
    '''Checks for/Creates a path'''
    if not os.path.isdir(path):
        os.makedirs(path)

class FastQCReport(object):
    '''Parses the output from FastQC into the appropriate files for the Reporter'''
    def __init__(self, config):
        # I/O Files 
        self.input_normal  = config['reads']['inputs']['normal']
        self.input_tumor   = config['reads']['inputs']['tumor']
        self.output_normal = config['reads']['outputs']['normal']
        self.output_tumor  = config['reads']['outputs']['tumor']

        self.final_files   = []

    def build(self):
        '''Parses the fastqc files to create the appropriate output files'''

        for n in self.generate_io(True):
            self.parse_fastqc_data(n['left_odir'], n['left_in_data'], n['left_o_data'],
                                   n['left_o_gcbd'], n['left_o_pbnc'], n['left_o_qbd'])
            if n['paired']:
                self.parse_fastqc_data(n['right_odir'], n['right_in_data'], n['right_o_data'],
                                       n['right_o_gcbd'], n['right_o_pbnc'], n['right_o_qbd'])

        for t in self.generate_io(False):
            self.parse_fastqc_data(t['left_odir'], t['left_in_data'], t['left_o_data'],
                                   t['left_o_gcbd'], t['left_o_pbnc'], t['left_o_qbd'])
            if t['paired']:
                self.parse_fastqc_data(t['right_odir'], t['right_in_data'], t['right_o_data'],
                                       t['right_o_gcbd'], t['right_o_pbnc'], t['right_o_qbd'])

    def generate_io(self, is_normal):
        '''Generates a dictionary of IO files'''
        curr = zip(self.input_normal, self.output_normal) if is_normal else \
               zip(self.input_tumor, self.output_tumor)

        for ip,op in curr:
            assert ip['readgroup'] == op['readgroup'], 'Inputs and outputs are not in the correct order!!'
            paired = ip['paired']

            # Leftseq
            left_odir    = op['leftseq']
            check_dir(left_odir)
            left_in_data = ip['leftseq']['data']
            left_o_data  = os.path.join(left_odir, 'fastqc_data.txt')
            left_o_gcbd  = os.path.join(left_odir, 'fastqc_gcbd.txt')
            left_o_pbnc  = os.path.join(left_odir, 'fastqc_pbnc.txt')
            left_o_qbd   = os.path.join(left_odir, 'fastqc_qbd.txt')

            # Rightseq
            right_odir    = op['rightseq'] if paired else None
            if paired: check_dir(right_odir)
            right_in_data = ip['rightseq']['data'] if paired else None
            right_o_data  = os.path.join(right_odir, 'fastqc_data.txt') if paired else None
            right_o_gcbd  = os.path.join(right_odir, 'fastqc_gcbd.txt') if paired else None
            right_o_pbnc  = os.path.join(right_odir, 'fastqc_pbnc.txt') if paired else None
            right_o_qbd   = os.path.join(right_odir, 'fastqc_qbd.txt') if paired else None

            yield {'paired': paired, 'left_odir': left_odir, 'left_in_data': left_in_data, 
                    'left_o_data': left_o_data, 'left_o_gcbd': left_o_gcbd, 'left_o_pbnc': left_o_pbnc, 
                    'left_o_qbd': left_o_qbd, 'right_odir': right_odir, 'right_in_data': right_in_data, 
                    'right_o_data': right_o_data, 'right_o_gcbd': right_o_gcbd, 'right_o_pbnc': right_o_pbnc, 
                    'right_o_qbd': right_o_qbd}

    def parse_fastqc_data(self, odir, indata, odata, ogcbd, opbnc, oqbd):
        '''Parses the fastqc_data file'''
        section = None 
        curr    = []

        with open(odata, 'wb') as o:
            for line in open(indata, 'rU'):
                o.write(line)
                if line.startswith('#'): continue 
                elif line.startswith('>>') and not section:
                    section = line.rstrip().split("\t")[0].replace('>>', '')
                elif line.startswith('>>END_MODULE'):
                    self.__process_section(section, curr, ogcbd, opbnc, oqbd)
                    section = None
                    curr = []
                else:
                    curr.append(line.rstrip()) 

    def __process_section(self, section, data, gcbd, pbnc, qbd):
        '''Processes a section of the FastQC data file and writes to appropriate output file'''
        if section == "Per base sequence quality":
            with open(qbd, 'wb') as o:
                o.write("Base\tMean\tMedian\tLQ\tUQ\tP10\tP90\n")
                for d in data:
                    o.write(d + "\n")
        elif section == "Per sequence GC content":
            with open(gcbd, 'wb') as o:
                o.write("PGC\tCount\n")
                for d in data:
                    o.write(d + "\n")
        elif section == "Per base N content":
            with open(pbnc, 'wb') as o:
                o.write("Base\tCount\n")
                for d in data:
                    o.write(d + "\n")
        else: pass

class AlignmentReport(object):
    '''Parses the output from Picard metrics and bedtools into the appropriate files for the Reporter'''
    def __init__(self, config):
        # I/O Files 
        self.inputs      = config['alignments']['inputs']
        self.outputs     = config['alignments']['outputs']
        self.sample      = config['parameters']['sample']
        self.final_files = []

    def build(self):
        '''Main function for creating the output files'''
        self.write_alignment_summary_metrics()
        self.write_insert_size_summary_metrics()
        self.write_coverage_summary()

    def write_alignment_summary_metrics(self):
        '''Summarizes the alignment_summary_metrics files'''
        header  = 'Sample\tType\tAln\tMetric\tCategory\tMeasure\tValue\n'
        keep    = ['TOTAL_READS', 'PF_READS_ALIGNED', 'PCT_PF_READS_ALIGNED',
                   'PF_HQ_ALIGNED_READS', 'PF_MISMATCH_RATE', 'PF_HQ_ERROR_RATE',
                   'PF_INDEL_RATE', 'MEAN_READ_LENGTH', 'PCT_READS_ALIGNED_IN_PAIRS',
                   'STRAND_BALANCE', 'PCT_CHIMERAS']
        lookup  = {'FIRST_OF_PAIR': 'R1', 'SECOND_OF_PAIR': 'R2', 'PAIR': 'Paired',
                   'UNPAIRED': 'Unpaired'} 
        head    = []
        section = "alignment_summary_metrics"
        out     = self.outputs['alignment_summary_metrics']

        with open(out, 'wb') as o:
            o.write(header)
            for aln in self.inputs:
                n_aln   = self.inputs[aln]['normal']['alignment_summary_metrics']
                t_aln   = self.inputs[aln]['tumor']['alignment_summary_metrics']

                for fil,pheno in zip([n_aln, t_aln], ['NORMAL', 'TUMOR']):
                    for line in open(fil, 'rU'):
                        if line.startswith('#'): continue
                        elif line.startswith('\n'): continue
                        else:
                            if line.startswith('CATEGORY'): head = line.rstrip().split("\t")
                            else:
                                cols = line.rstrip("\n").split("\t")
                                dat  = dict(zip(head, cols))
                                for k in keep:
                                    row = "\t".join([
                                        self.sample, pheno, aln, section, lookup[dat['CATEGORY']],
                                        k, dat[k]
                                    ])
                                    o.write(row + "\n") 

    def write_insert_size_summary_metrics(self):
        '''Summarizes the insert size metrics'''
        header  = 'Sample\tType\tAln\tMetric\tCategory\tMeasure\tValue\n'
        section = 'insert_size_metrics'
        out     = self.outputs['insert_size_metrics']
        head    = []

        with open(out, 'wb') as o:
            o.write(header)
            header  = 'Sample\tType\tAln\tMetric\tCategory\tMeasure\tValue\n'
            for aln in self.inputs:
                n_ins   = self.inputs[aln]['normal']['insert_size_metrics']
                t_ins   = self.inputs[aln]['tumor']['insert_size_metrics']
                
                for fil,pheno in zip([n_ins, t_ins], ['NORMAL', 'TUMOR']):
                    metflag  = False
                    histflag = False
                    try:
                        for line in open(fil, 'rU'):
                            if line.startswith('#'): continue
                            elif line.startswith('\n'): continue
                            elif line.startswith('insert_size'): histflag = True
                            elif line.startswith('MEDIAN'):
                                head    = line.rstrip().split('\t')
                                metflag = True
                            elif histflag:
                                row = "\t".join([
                                    self.sample, pheno, aln, section, "HIST", line.rstrip()
                                ])
                                o.write(row + "\n")
                            elif metflag:
                                dat = dict(zip(head, line.rstrip('\n').split('\t')))
                                for c in head:
                                    if c not in ['SAMPLE', 'LIBRARY', 'READ_GROUP']:
                                        row = "\t".join([
                                            self.sample, pheno, aln, section, "TABLE", c, dat[c]
                                        ])
                                        o.write(row + "\n")

                    except IOError:
                        row = "\t".join([self.sample, pheno, aln, section, "TABLE", 
                                         'NA', 'NA', 'NA'])
                        o.write(row + '\n') 

                        row = "\t".join([self.sample, pheno, aln, section, "HIST", 
                                         'NA', 'NA', 'NA'])
                        o.write(row + '\n') 

    def write_coverage_summary(self):
        '''Summarizes the coverage'''
        header  = 'Sample\tType\tAln\tCoverage\n'
        out     = self.outputs['total_coverage']
        head    = []

        with open(out, 'wb') as o:
            o.write(header)
            for aln in self.inputs:
                n_cov = self.inputs[aln]['normal']['total_coverage']
                t_cov = self.inputs[aln]['tumor']['total_coverage']
                for fil,pheno in zip([n_cov, t_cov], ['NORMAL', 'TUMOR']):
                    curr_sum  = 0
                    curr_size = None 
                    for line in open(fil, 'rU'):
                        cols = line.rstrip().split("\t")
                        if cols[0] == "all":
                            binsize = int(cols[1])
                            bincts  = int(cols[2])
                            if curr_size is None:
                                curr_size = int(cols[3])
                            curr_sum += binsize * bincts
                    row = "{0}\t{1}\t{2}\t{3:.2f}".format(
                          self.sample, pheno, aln, curr_sum / float(curr_size))
                    o.write(row + "\n") 

class SomaticReport(object):
    '''Parses the output from Annovar into the appropriate files for the Reporter'''
    def __init__(self, config):
        # I/O Files 
        self.inputs      = config['somatic']['inputs']
        self.outputs     = config['somatic']['outputs']['smd_snp_table']
        self.sample      = config['parameters']['sample']
        self.keep        = config['options']['reporter_colnames'] 
        self.smd_dic     = {}
        self.combos      = []
        self.final_files = []

    def build(self):
        '''Parses the Annovar outputs and creates the appropriately formatted
           somatic SNP table for the reporter'''
        # First, load the SMD dictionary
        self.__load_smd_dic()

        # Then, write the new file
        self.__write_concordance_file()

    def __load_smd_dic(self):
        '''Loops over all Annovar files and builds a dictionary for estimating score'''
        for obj in self.inputs:
            self.combos.append('{smd}.{aln}'.format(**obj))
            vcf_head  = self.__get_vcf_header(obj['vcf'])
            if obj['smd'] == 'mutect': self.__process_mutect(**obj)
            elif obj['smd'] == 'shimmer': self.__process_shimmer(**obj)
            elif obj['smd'] == 'sniper': self.__process_sniper(**obj)
            elif obj['smd'] == 'strelka': self.__process_strelka(**obj)
            elif obj['smd'] == 'varscan': self.__process_varscan(**obj)
            elif obj['smd'] == 'virmid': self.__process_virmid(**obj)
            else:
                raise RuntimeError("Unknown SMD '{0}'".format(obj['smd']))

    def __write_concordance_file(self):
        '''Writes the reporter-formatted variant file'''
        with open(self.outputs, 'wb') as o:
            o.write('sample\t' + '\t'.join(self.keep) + '\tTotalScore\tAligner\tCaller\t' + \
                    '\t'.join(self.combos) + "\n") 

            for k in sorted(self.smd_dic, key=itemgetter(0,1)):
                for a in self.smd_dic[k]:
                    row    = [self.sample, self.smd_dic[k][a]['meta']]
                    combo  = self.smd_dic[k][a]['combo']
                    aln_ct = len(set([i.split('.')[1] for i in combo]))
                    smd_ct = len(set([i.split('.')[0] for i in combo]))
                    score  = aln_ct * smd_ct
                    row    = row + map(str, [score, aln_ct, smd_ct]) + \
                             ["1" if i in combo else "0" for i in self.combos] 
                    o.write("\t".join(row) + "\n")

    def __process_mutect(self, smd, aln, vcf, annovar):
        '''Parses the mutect vcf'''
        vcf_head  = self.__get_vcf_header(vcf)
        anno_head = []
        combo     = 'mutect.' + aln
        for line in open(annovar, 'rU'):
            if not anno_head:
                anno_head = line.rstrip().split("\t")[:-1] + vcf_head
            else:
                cols = line.rstrip("\n").split("\t")
                dat  = dict(zip(anno_head, cols))
                if dat['FILTER'] == "PASS":
                    ttype = [i.split('=')[1] for i in dat['INFO'].split(';') if i.startswith('TTYPE')][0]
                    ntype = [i.split('=')[1] for i in dat['INFO'].split(';') if i.startswith('NTYPE')][0]
                    vtype = [i.split('=')[1] for i in dat['INFO'].split(';') if i.startswith('VT')][0]
                    if ttype == 'SOMATIC' and vtype == 'SNP':
                        self.__update_dic(dat, combo)

    def __process_shimmer(self, smd, aln, vcf, annovar):
        '''Parses the shimmer vcf'''
        vcf_head  = self.__get_vcf_header(vcf)
        anno_head = []
        combo     = 'shimmer.' + aln
        for line in open(annovar, 'rU'):
            if not anno_head:
                anno_head = line.rstrip().split("\t")[:-1] + vcf_head
            else:
                cols = line.rstrip("\n").split("\t")
                dat  = dict(zip(anno_head, cols))
                if dat['FILTER'] == "PASS":
                    self.__update_dic(dat, combo)

    def __process_sniper(self, smd, aln, vcf, annovar):
        '''Parses the somatic sniper vcf'''
        vcf_head  = self.__get_vcf_header(vcf)
        anno_head = []
        combo     = 'sniper.' + aln
        for line in open(annovar, 'rU'):
            if not anno_head:
                anno_head = line.rstrip().split("\t")[:-1] + vcf_head
            else:
                cols = line.rstrip("\n").split("\t")
                dat  = dict(zip(anno_head, cols))
                if dat['FILTER'] == "PASS":
                    ttype = [i.split('=')[1] for i in dat['INFO'].split(';') if i.startswith('TTYPE')][0]
                    ntype = [i.split('=')[1] for i in dat['INFO'].split(';') if i.startswith('NTYPE')][0]
                    if ttype == 'SOMATIC' and ntype == 'REF':
                        self.__update_dic(dat, combo)

    def __process_strelka(self, smd, aln, vcf, annovar):
        '''Parses the strelka vcf'''
        vcf_head  = self.__get_vcf_header(vcf)
        anno_head = []
        combo     = 'strelka.' + aln
        for line in open(annovar, 'rU'):
            if not anno_head:
                anno_head = line.rstrip().split("\t")[:-1] + vcf_head
            else:
                cols = line.rstrip("\n").split("\t")
                dat  = dict(zip(anno_head, cols))
                if dat['FILTER'] == "PASS":
                    if 'SOMATIC' in dat['INFO'] and 'INDEL' not in dat['INFO']:
                        self.__update_dic(dat, combo)

    def __process_varscan(self, smd, aln, vcf, annovar):
        '''Parses the varscan vcf'''
        vcf_head  = self.__get_vcf_header(vcf)
        anno_head = []
        combo     = 'varscan.' + aln
        for line in open(annovar, 'rU'):
            if not anno_head:
                anno_head = line.rstrip().split("\t")[:-1] + vcf_head
            else:
                cols = line.rstrip("\n").split("\t")
                dat  = dict(zip(anno_head, cols))
                if dat['FILTER'] == "PASS":
                    ttype = [i.split('=')[1] for i in dat['INFO'].split(';') if i.startswith('SS')][0]
                    if ttype == 'SOMATIC' and 'INDEL' not in dat['INFO']: 
                        self.__update_dic(dat, combo)

    def __process_virmid(self, smd, aln, vcf, annovar):
        '''Parses the virmid vcf'''
        vcf_head  = self.__get_vcf_header(vcf)
        anno_head = []
        combo     = 'virmid.' + aln
        for line in open(annovar, 'rU'):
            if not anno_head:
                anno_head = line.rstrip().split("\t")[:-1] + vcf_head
            else:
                cols = line.rstrip("\n").split("\t")
                dat  = dict(zip(anno_head, cols))
                if dat['FILTER'] == "PASS":
                    ttype = [i.split('=')[1] for i in dat['INFO'].split(';') if i.startswith('TYPE')][0]
                    if ttype == 'SOMATIC': 
                        self.__update_dic(dat, combo)

    def __update_dic(self, dat, combo):
        '''Updates the smd_dic'''
        key     = (dat['Chr'], int(dat['Start']))
        alleles = tuple(dat['Ref'].upper().split(',') + dat['Alt'].upper().split(','))

        if key not in self.smd_dic: self.smd_dic[key] = {} 
        if alleles not in self.smd_dic[key]: 
            self.smd_dic[key][alleles] = {
                'meta': "\t".join([dat[i] if dat[i] else "NA" for i in self.keep]), 'combo': []}
        self.smd_dic[key][alleles]['combo'].append(combo) 

    def __get_vcf_header(self, vcf):
        '''Parses out the #CHROM line and returns it as a list'''
        head = []
        for line in open(vcf, 'rU'):
            if line.startswith('##'): continue
            elif line.startswith('#'): head = line.rstrip()[1::].split("\t")
            else:
                break
        return head

if __name__ == '__main__':
    start = time.time()

    # Set up logger
    logger    = logging.getLogger('BuildReporter')
    logger.setLevel(logging.INFO)
    ch        = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(levelname)s] [%(asctime)s] [%(name)s] - %(message)s',
                                  datefmt='%H:%M:%S')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Print header
    logger.info('-'*75)
    logger.info('BuildReporterFiles.py v{0}'.format(__version__))
    logger.info('Copyright (c) 2014, {0}'.format(__license__))
    logger.info('The Center for Research Informatics, The University of Chicago')
    logger.info('For support and documentation go to:')
    logger.info('{0}'.format(__url__))
    logger.info("Program Args: BuildReporterFiles.py " + " ".join(sys.argv[1::]))
    logger.info('Date/time: {0}'.format(datetime.datetime.now()))
    logger.info('-'*75)
    logger.info('-'*75)

    # Set up parser
    parser = argparse.ArgumentParser(prog='BuildReporterFiles.py',
               epilog='Copyright 2014, Center for Research Informatics, University of Chicago')
    parser.add_argument('config_file', help='the YAML configuration file')
    args   = parser.parse_args()
    config_file = None

    with open(args.config_file, 'rU') as fh:
        config_file = yaml.safe_load(fh)

    logger.info('Creating FastQC files...')
    fastqc_creator = FastQCReport(config_file)
    fastqc_creator.build()

    logger.info('Creating Alignment files...')
    align_creator = AlignmentReport(config_file)
    align_creator.build()

    logger.info('Creating Somatic Mutation files...')
    somatic_creator = SomaticReport(config_file)
    somatic_creator.build()

    # Done
    logger.info("Finished, took {0} seconds.".format(time.time() - start))
