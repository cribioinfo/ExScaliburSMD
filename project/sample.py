'''Contains a factory object for creating Tumor/Normal Paired samples'''
import os
import re
import logging
import yaml

import utils 
from read import IlluminaSeq 

class SampleFactory(object):
    '''Factory object takes the metadata file and parses it. Validates the table
       and creates a dictionary of each sample.
    '''
    possible_phenotypes = set(['TUMOR', 'NORMAL']) 
    sample_obj  = {}
    logger      = logging.getLogger('SMD.SampleFactory')

    def __iter__(self):
        '''Iterates over each tumor/normal pair in the sample_obj dictionary'''
        for sample in self.sample_obj:
            yield self.sample_obj[sample]

    def __getitem__(self, key):
        return self.sample_obj[key]

    def load_samples(self, metadata, PATH, project_id):
        """Parses the metadata file and validates it. Creates the sample_obj dictionary.

           -- metadata - the metadata file
           -- PATH - the ExScaliburSMD package path
           -- project_id - the user provided project name
        """
        header      = []
        required    = ['sample', 'phenotype', 'lane', 'readgroup',
                       'chemistry', 'path', 'fastq_r1']
        total_lines = 0

        for line in open(metadata, 'rU'):
            # Skip over commented-out lines
            if line.startswith('#'): continue
            # Get header row
            elif not header:
                header = [i.lower() for i in line.rstrip().split("\t")]
            else:
                # Create a dictionary for each element in the row
                cols            = [i.strip() for i in line.rstrip('\n').split("\t")]
                curr            = dict(zip(header, cols))
                curr["project"] = project_id
                total_lines    += 1
                
                # Check required elements
                for r in required:
                    assert r in curr and curr[r], \
                        self.logger.error(
                            'Required element {0} is missing from metadata file line {0}'.format(
                                     r, total_lines))

                # Now, fill in all required elements
                curr['phenotype'] = self._phenotype(curr, total_lines)
                curr_read         = IlluminaSeq(curr, total_lines)
 
                # Create/add to the Sample object and add to the sample dictionary
                if curr['sample'] not in self.sample_obj:
                    self.sample_obj[curr['sample']] = TumorNormalPair(curr['sample'])
                self.sample_obj[curr['sample']].add_read(curr, curr_read)
        self._validate(total_lines)
 
    def _validate(self, total_lines):
        '''Validates the samples'''
        # Each row must have a unique readgroup
        readgroup_list = []
        for s in self.sample_obj:
            [readgroup_list.append(j.readgroup) for j in self.sample_obj[s].reads()]
        readgroup_set  = set(readgroup_list)
        assert len(readgroup_set) == total_lines, \
            self.logger.error('There should be {0} unique read IDs, but there was only {1}!!'.format(
                         total_lines, len(readgroup_set)))
        self.logger.info('There are {0} Tumor/Normal Pairs'.format(len(self.sample_obj)))
        self.logger.info('There are {0} unique readgroups'.format(total_lines))

    def _phenotype(self, dic, line_number):
        '''Sets the current sample read type. This must be either tumor or normal
        
        -- dic - A row of the metadata table as a dictionary
        -- line_number - The current metadata line number
        '''
        chk = dic['phenotype'].upper()
        assert chk in self.possible_phenotypes, \
            self.logger.error('Unknown sample phenotype {0} in row {1} -- Must be either TUMOR or NORMAL'.format(
                          chk, line_number))
        return chk

class TumorNormalPair(object):
    '''Represents a tumor/normal pair.

    -- name - the sample id
    -- project - the project id
    -- normal_reads - list of all the normal IlluminaSeq objects
    -- tumor_reads - list of all the tumor IlluminaSeq objects
    -- reads - a generator for all normal and tumor IlluminaSeq objects
    -- preprocessing - container for the Preprocessing module object for this sample
    -- alignments - container for the Alignments module object for this sample
    -- postprocessing - container for the Alignments Posptrocessing sub-module object for this sample
    -- gatk - container for the GATK Processing module object for this sample
    -- somatic - container for the Somatic Variant Calling module object for this sample
    '''
    def __init__(self, sample_id):
        self.logger         = logging.getLogger('SMD.Pair')
        self.name           = sample_id
        self.project        = None
        self.normal_reads   = []
        self.tumor_reads    = []
        self.reads          = self.__read_generator
        self.reporter_obj   = {} 
        self.sample_cfg     = None
        self.strelka_cfg    = None
        self.report_yaml    = None
        self.reporter_files = []

    ########################################################
    ## Public member functions accessed by the Project object
    def add_read(self, meta, read_obj):
        '''Adds a IlluminaSeq object to the current sample.

        meta - The current row from the metadata file
        read_obj - The current IlluminaSeq instance
        '''
        if not self.project: self.project = meta['project']
        assert meta['project'] == self.project

        if meta['phenotype'].upper() == 'TUMOR':
            self.tumor_reads.append(read_obj)
        elif meta['phenotype'].upper() == 'NORMAL':
            self.normal_reads.append(read_obj)

    def write_sample_bds(self, pobj, settings):
        '''Write the config files for this sample for the BDS
           scripts to use.

           pobj - the project object instance
           settings - the project settings object instance
        '''
        # Set config file
        self.__set_cfg(pobj.files['config_path'])
        
        # Print logs
        self.logger.info('Sample={0.name} | Phenotype=NORMAL | Readgroups={1}'.format(
          self, len(self.normal_reads)))
        self.logger.info('Sample={0.name} | Phenotype=TUMOR | Readgroups={1}'.format(
          self, len(self.tumor_reads)))

        # Set project reporter obj
        self.__set_reporter_obj(pobj, settings)

        # Write the main BDS config
        self.__write_main_bds_config(pobj, settings)

        # Write BDS configs for each readgroup
        self.__write_readgroup_bds_configs(pobj, settings)

    def write_sample_report_yaml(self):
        '''Writes the reporter_obj to a YAML file for the reporter processing'''
        with open(self.report_yaml, 'wb') as o:
            yaml.safe_dump(self.reporter_obj, o)

    ########################################################
    ## Private member functions
    def __read_generator(self):
        '''Create a generator of IlluminaSeq instances that are associated
        with this sample.

        For example:
        >>> for read in self.read_generator(): print read
        '''
        for r in self.normal_reads + self.tumor_reads:
            yield r

    def __set_cfg(self, cfgdir):
        '''Sets the config files and checks the directory'''
        cfgdir = os.path.join(cfgdir, self.name)
        utils.check_dir(cfgdir)
        self.sample_cfg = os.path.join(cfgdir, '{0.name}.cfg'.format(self))
        self.report_yaml = os.path.join(cfgdir, '{0.name}.build_reporter.yaml'.format(self))

    def __set_reporter_obj(self, pobj, settings):
        '''Builds the basics of the reporter_obj'''
        self.reporter_obj = {
          'module': 'PipelineReport',
          'software': {
            'python': settings.system['python']
          },
          'parameters': {'project': pobj.name, 'sample': self.name},
          'options': settings.options['annovar'],
          'reads': {
            'inputs': {'normal': [], 'tumor': []},
            'outputs': {'normal': [], 'tumor': []} 
          },
          'alignments': {'inputs': {}, 'outputs': {}},
          'somatic': {'inputs': [], 'outputs': {}}
        }
       
        # INPUTS
        for a in ['bwa_aln', 'bwa_mem', 'novoalign']:
            if a in settings.aln_list:
                # ALIGNMENTS
                self.reporter_obj['alignments']['inputs'][a] = {
                  'normal': {
                    'alignment_summary_metrics': os.path.join(pobj.files['results_path'], \
                      '{0.name}/01_Alignments/{1}/{0.name}-NORMAL.{1}.metrics.alignment_summary_metrics'.format(
                      self, a)),
                    'insert_size_metrics': os.path.join(pobj.files['results_path'], \
                      '{0.name}/01_Alignments/{1}/{0.name}-NORMAL.{1}.metrics.insert_size_metrics'.format(
                      self, a)),
                    'total_coverage': os.path.join(pobj.files['results_path'], \
                      '{0.name}/01_Alignments/{1}/{0.name}-NORMAL.{1}.{2}.exons.bed'.format(
                      self, a, pobj.assembly['refname']))
                  },
                  'tumor': {
                    'alignment_summary_metrics': os.path.join(pobj.files['results_path'], \
                      '{0.name}/01_Alignments/{1}/{0.name}-TUMOR.{1}.metrics.alignment_summary_metrics'.format(
                      self, a)),
                    'insert_size_metrics': os.path.join(pobj.files['results_path'], \
                      '{0.name}/01_Alignments/{1}/{0.name}-TUMOR.{1}.metrics.insert_size_metrics'.format(
                      self, a)),
                    'total_coverage': os.path.join(pobj.files['results_path'], \
                      '{0.name}/01_Alignments/{1}/{0.name}-TUMOR.{1}.{2}.exons.bed'.format(
                      self, a, pobj.assembly['refname']))
                  },
                }
               
                # SMD
                for s in ['mutect', 'shimmer', 'sniper', 'strelka', 'varscan', 'virmid']:
                    if s in settings.smd_list: 
                        self.reporter_obj['somatic']['inputs'].append({
                          'aln': a, 'smd': s,
                          'annovar': os.path.join(pobj.files['results_path'], \
                            '{0.name}/04_VariantAnnotation/'.format(self) + \
                            '{0.name}.{1}.{2}.vcf.annovar.anno.{3}_multianno.txt'.format(
                            self, s, a, pobj.assembly['refname'])),
                          'vcf': os.path.join(pobj.files['results_path'], \
                            '{0.name}/03_SomaticMutations/'.format(self) + \
                            '{0.name}.{1}.{2}.final.vcf'.format(self, s, a))
                        })

        # OUTPUTS
        self.reporter_obj['alignments']['outputs'] = {
          'alignment_summary_metrics': os.path.join(pobj.files['report_path'], \
            'project_alignments/{0.name}.alignment_summary_metrics.txt'.format(self)),
          'insert_size_metrics': os.path.join(pobj.files['report_path'], \
            'project_alignments/{0.name}.insert_size_metrics.txt'.format(self)),
          'total_coverage': os.path.join(pobj.files['report_path'], \
            'project_alignments/{0.name}.total_coverage.txt'.format(self)),
          'path': os.path.join(pobj.files['report_path'], 'project_alignments')
        }
        self.reporter_obj['somatic']['outputs'] = {
          'smd_snp_table': os.path.join(pobj.files['report_path'], \
            'project_somatic/{0.name}.somatic.snps.tsv'.format(self))
        }

        # Append to reporter_files
        anames = ['alignment_summary_metrics', 'insert_size_metrics', 'total_coverage']
        [self.reporter_files.append(self.reporter_obj['alignments']['outputs'][i]) for i in anames]
        self.reporter_files.append(self.reporter_obj['somatic']['outputs']['smd_snp_table'])

        # Check output paths
        utils.check_dir(self.reporter_obj['alignments']['outputs']['path']) 
        utils.check_dir(os.path.abspath(os.path.dirname(self.reporter_obj['somatic']['outputs']['smd_snp_table']))) 

    def __write_main_bds_config(self, pobj, settings):
        '''Creates the main BDS config for this sample

           pobj - the project object instance
           settings - the project settings object instance
        '''
        # Normal readgroup configs
        normal_configs = [os.path.join(pobj.files['config_path'], '{0.name}/{1.readgroup_cfg}'.format(self, i)) \
                          for i in self.normal_reads]
        tumor_configs = [os.path.join(pobj.files['config_path'], '{0.name}/{1.readgroup_cfg}'.format(self, i)) \
                          for i in self.tumor_reads]

        # Write the main sample configs
        with open(self.sample_cfg, 'wb') as o:
            o.write('# Sample\n')
            o.write('sample = {0.name}\n'.format(self)) 
            o.write('normalConfigs = {0}\n'.format(",".join(normal_configs)))
            o.write('tumorConfigs = {0}\n\n'.format(",".join(tumor_configs)))
            
            o.write('# Project\n')
            o.write('projectDir = {0}\n'.format(pobj.files['project_path']))
            o.write('resultsDir = {0}\n'.format(os.path.join(pobj.files['results_path'], self.name)))
            o.write('logsDir = {0}\n'.format(pobj.files['log_path']))
            o.write('split = {0}\n'.format('true' if settings.is_split else 'false'))
            o.write('max_split = {0}\n'.format(settings.nsplits))
            if settings.args.strelka: 
                self.strelka_cfg = os.path.join(pobj.files['config_path'], \
                  '{0.name}/{0.name}.strelka.cfg'.format(self))
                o.write('strelka_config = {0.strelka_cfg}\n\n'.format(self))
                self.__make_strelka_config(settings)
            else:
                o.write('\n') 

            o.write('# Which preprocessing\n')
            o.write('runFastqc = true\n')
            o.write('runClipping = {0}\n\n'.format('true' if settings.clipping else 'false'))
           
            o.write('# Which aligners\n')
            o.write('runBwaAln = {0}\n'.format('true' if 'bwa_aln' in settings.aln_list else 'false')) 
            o.write('runBwaMem = {0}\n'.format('true' if 'bwa_mem' in settings.aln_list else 'false')) 
            o.write('runNovoalign = {0}\n\n'.format('true' if 'novoalign' in settings.aln_list else 'false')) 

            o.write('# Which SMDs\n')
            o.write('runMutect = {0}\n'.format('true' if 'mutect' in settings.smd_list else 'false')) 
            o.write('runShimmer = {0}\n'.format('true' if 'shimmer' in settings.smd_list else 'false')) 
            o.write('runSomaticSniper = {0}\n'.format('true' if 'sniper' in settings.smd_list else 'false')) 
            o.write('runStrelka = {0}\n'.format('true' if 'strelka' in settings.smd_list else 'false')) 
            o.write('runVarScan = {0}\n'.format('true' if 'varscan' in settings.smd_list else 'false')) 
            o.write('runVirmid = {0}\n\n'.format('true' if 'virmid' in settings.smd_list else 'false'))

            o.write('# Reporter\n')
            o.write('reporterConfig = {0.report_yaml}\n'.format(self))
            o.write('reporterOutputs = {0}\n\n'.format(",".join(self.reporter_files)))

    def __write_readgroup_bds_configs(self, pobj, settings):
        '''Creates the BDS configs for each readgroup from this sample. 

           pobj - the project object instance
           settings - the project settings object instance
        '''
        # Now we need to write the readgroup configs
        for r in self.__read_generator():
            curr_cfg = os.path.join(pobj.files['config_path'], '{0.name}/{1.readgroup_cfg}'.format(self, r))
            with open(curr_cfg, 'wb') as o:
                o.write('# Sample info\n')
                o.write('sample = {0.name}\n'.format(self))
                o.write('phenotype = {0.phenotype}\n'.format(r))
                o.write('offset = {0.offset}\n\n'.format(r))

                o.write('# Readgroup info\n')
                o.write('paired = {0}\n'.format('true' if r.paired else 'false'))
                o.write('readgroup = {0.readgroup}\n'.format(r))
                o.write('rgstring = {0.rgstring}\n'.format(r))
                o.write('fragment = {0.fragment}\n'.format(r))
                o.write('fragment_sd = {0.sdfragment}\n\n'.format(r))

                o.write('# Input fastq\n')
                o.write('left = {0}\n'.format(os.path.join(r.path, r.leftseq)))
                if r.paired: o.write('right = {0}\n\n'.format(os.path.join(r.path, r.rightseq)))
                else: o.write('\n')

                self.__write_fastqc(o, pobj, r)

                if settings.clipping:
                    self.__write_clipping(o, pobj, r)
                else:
                    o.write("# Adapter clipping\n")
                    o.write("clip = false\n\n")

    def __write_fastqc(self, o, pobj, r):
        # Dict for reporter 
        incurr = {
          'leftseq': {'data': None,'summary': None}, 
          'rightseq': {'data': None,'summary': None}, 
          'paired': r.paired, 
          'readgroup': r.readgroup
        }
        ocurr = {
          'leftseq': None, 
          'rightseq': None,
          'paired': r.paired,
          'readgroup': r.readgroup
        }

        # Set files
        fqdir = os.path.join(pobj.files['results_path'], \
                  '{0.name}/00_Preprocessing/fastqc/{1.readgroup}'.format(self, r))
        o.write('# Fastqc\n')
        o.write('fastqcResultsDir = {0}\n'.format(fqdir))

        raw1 = r.leftseq 
        out1 = [i for i in [re.search('([\w\_\-\.\d]+)\.txt\.gz$', raw1),
                    re.search('([\w\_\-\.\d]+)\.txt$', raw1),
                    re.search('([\w\_\-\.\d]+)\.fastq\.gz$', raw1),
                    re.search('([\w\_\-\.\d]+)\.fastq$', raw1)] if i][0].group(1)
        o.write('leftFastqcResultsDir = {0}_fastqc\n'.format(out1)) 
        incurr['leftseq']['data'] = os.path.join(fqdir, '{0}_fastqc/fastqc_data.txt'.format(out1))
        incurr['leftseq']['summary'] = os.path.join(fqdir, '{0}_fastqc/summary.txt'.format(out1))
        ocurr['leftseq'] = os.path.join(pobj.files['report_path'], \
          'project_reads/{0.name}/{1.readgroup}/leftseq'.format(self, r))
        utils.check_dir(ocurr['leftseq'])
        self.reporter_files.append(os.path.join(ocurr['leftseq'], 'fastqc_data.txt')) 
        self.reporter_files.append(os.path.join(ocurr['leftseq'], 'fastqc_gcbd.txt')) 
        self.reporter_files.append(os.path.join(ocurr['leftseq'], 'fastqc_pbnc.txt')) 
        self.reporter_files.append(os.path.join(ocurr['leftseq'], 'fastqc_qbd.txt')) 

        if r.paired:
            raw2 = r.rightseq 
            out2 = [i for i in [re.search('([\w\_\-\.\d]+)\.txt\.gz$', raw2),
                        re.search('([\w\_\-\.\d]+)\.txt$', raw2),
                        re.search('([\w\_\-\.\d]+)\.fastq\.gz$', raw2),
                        re.search('([\w\_\-\.\d]+)\.fastq$', raw2)] if i][0].group(1)
            o.write('rightFastqcResultsDir = {0}_fastqc\n'.format(out2))
            incurr['rightseq']['data'] = os.path.join(fqdir, '{0}_fastqc/fastqc_data.txt'.format(out2))
            incurr['rightseq']['summary'] = os.path.join(fqdir, '{0}_fastqc/fastqc_summary.txt'.format(out2))
            ocurr['rightseq'] = os.path.join(pobj.files['report_path'], \
              'project_reads/{0.name}/{1.readgroup}/rightseq'.format(self, r))
            utils.check_dir(ocurr['rightseq'])
            self.reporter_files.append(os.path.join(ocurr['rightseq'], 'fastqc_data.txt')) 
            self.reporter_files.append(os.path.join(ocurr['rightseq'], 'fastqc_gcbd.txt')) 
            self.reporter_files.append(os.path.join(ocurr['rightseq'], 'fastqc_pbnc.txt')) 
            self.reporter_files.append(os.path.join(ocurr['rightseq'], 'fastqc_qbd.txt')) 
        o.write('\n')

        # Add to dict
        self.reporter_obj['reads']['inputs'][r.phenotype.lower()].append(incurr)
        self.reporter_obj['reads']['outputs'][r.phenotype.lower()].append(ocurr)

    def __write_clipping(self, o, pobj, r):
        o.write('# Adapter trimming\n')
        if not r.paired:
            if pobj.settings.options['cutadapt']['adapter1']:
                o.write('clip = true\n')
                o.write('clipResultsDir = {0}\n'.format(
                  os.path.join(pobj.files['results_path'], '{0.name}/00_Preprocessing/clipping'.format(self))))
                o.write('leftClip = {0.readgroup}_R1.clip.fastq.gz\n'.format(r))
            else:
                o.write('clip = false\n')
        else:
            o.write('clip = true\n')
            o.write('clipResultsDir = {0}\n'.format(
              os.path.join(pobj.files['results_path'], '{0.name}/00_Preprocessing/clipping'.format(self))))
            o.write('leftClip = {0.readgroup}_R1.clip.fastq.gz\n'.format(r)) 
            o.write('rightClip = {0.readgroup}_R2.clip.fastq.gz\n'.format(r)) 
            o.write('mergeClip = {0.readgroup}_ME.clip.fastq.gz\n'.format(r))

        o.write('\n')

    def __make_strelka_config(self, settings):
        '''Creates the config file for Strelka runs, if requested by user'''
        class Strelka(object):
            def __init__(self, dat):
                self.minq          = dat.min_map_q
                self.ignore_rnames = ' --ignore-conflicting-read-names' if \
                                     dat.options['strelka']['ignore_conflicting_read_names'] else ''
                self.indel_repeat  = dat.options['strelka']['indel_max_ref_repeat']
                self.indel_prior   = dat.options['strelka']['indel_prior']
                self.max_depth     = dat.options['strelka']['max_input_depth']
                self.max_spanning  = dat.options['strelka']['max_spanning_deletion_fraction']
                self.min_indel_q   = dat.options['strelka']['min_indel_quality']
                self.min_snv_q     = dat.options['strelka']['min_snv_quality']
                self.min_tier_2    = 0 if self.min_snv_q == 0 else 5 if 5 < self.min_snv_q else 5 / self.min_snv_q
                self.sindel_noise  = dat.options['strelka']['sindel_noise']
                self.skip_dp       = dat.options['strelka']['skip_depth_filters']
                self.snv_prior     = dat.options['strelka']['snv_prior']
                self.ssnv_noise    = dat.options['strelka']['ssnv_noise']

        sobj = Strelka(settings)
        with open(settings.options['strelka']['base_config'], 'rU') as fh:
            template = fh.read()
            config_string = template.format(sobj)
            with open(self.strelka_cfg, 'wb') as oc:
                oc.write(config_string) 
