'''Main class for handling ExScaliburSMD projects. This class contains all samples, settings
   files, etc. for the entire project.
'''
import logging
import os
import stat
import yaml

from sample import SampleFactory
from settings import PipelineSettings

from utils import check_dir

class Project(object):
    '''Object representing a Tumor/Normal Exome Project

       -- args - the argparser object
       -- PATH - the path to the ExScaliburSMD application directory
       -- assembly - the assembly dictionary from the assembly yaml file
       -- settings - the PipelineSettings object
       -- files    - dictionary of main project paths
       -- samples  - the SampleFactory object
       -- name - the project ID
    '''
    def __init__(self, args, PATH):
        self.logger              = logging.getLogger('SMD.Project')
        self.args                = args
        self.PATH                = PATH
        self.assembly            = None
        self.settings            = None
        self.files               = None
        self.name                = args.project_id 
        self.samples             = SampleFactory() 
        self.final_reporter_yaml = None
        self.bds_file            = None

        # Lists holding the job submission commands to write to the
        # submission scripts
        self.preprocessing_jobs = []
        self.alignments_jobs    = []
        self.gatk_jobs          = []
        self.somatic_jobs       = []
        self.annotation_jobs    = []
        self.report_jobs        = []
        self.config_list        = []

    def initialize_project(self):
        '''Wrapper for building the project.
          
           1. Load the settings
           2. Load the assembly
           3. Load the samples
           4. Set the directory structure
        '''
        self.logger.info('Project Name: {0.name}'.format(self))

        # Settings 
        self.settings = PipelineSettings(self.args, self.PATH)
        self.logger.info('Baseline Threads: {0.threads}'.format(self.settings))
        self.logger.info('Clip Adapters: {0.clipping}'.format(self.settings))
        if self.settings.clipping:
            self.logger.info('Minimum read length: {0.min_length}'.format(self.settings))
        self.logger.info('Minimum MAPQ: {0.min_map_q}'.format(self.settings))
        self.logger.info('Aligners: {0}'.format(self.settings.aln_list))
        self.logger.info('Somatic Mutation Detectors: {0}'.format(self.settings.smd_list))

        # Assembly info
        self._set_assembly()
        self.logger.info('Reference Name: {0}'.format(self.assembly['refname']))

        # Samples
        self.samples.load_samples(os.path.abspath(self.args.metadata), self.PATH, self.name)

        # File structure
        self._set_project_files()

    def configure_pipeline(self):
        '''Wrapper for creating the configuration files for each sample'''
        self.settings.write_settings_cfg(self)
        for s in self.samples:
            s.write_sample_bds(self, self.settings)
            s.write_sample_report_yaml()
            self.config_list.append(s.sample_cfg)
        self.__write_reporter_yaml()
        self.__write_bds_config()
        self.__write_project_run_script()

    def _set_assembly(self):
        '''Grabs the yaml file, validates it, and returns the dictionary'''
        # Load
        in_fh   = open(self.args.config_file, 'rU')
        chk     = yaml.safe_load(in_fh)['reference']
        in_fh.close()

        # Check the file
        refname = chk.keys()[0]
        required = ['cosmic', 'exons_bed', 'knowndb', 'knownindel', 
                    'referenceseq', 'sequence_dictionary']
        if 'novoalign' in self.settings.aln_list: required.append('novoalign_index') 
        for i in required:
            assert i in chk[refname], self.logger.error("Required refrence section '%s' missing" % i)
            if i in ['cosmic', 'exons_bed', 'knowndb', 'knownindel', 'referenceseq', 'sequence_dictionary']:
                try:
                    fp = open(chk[refname][i])
                except IOError as e:
                    raise IOError('Unable to open {0} file: {1}'.format(i, chk[refname][i]))
                else:
                    fp.close()

        # Set refname and assign the dictionary to the member 
        chk[refname]['refname'] = refname
        self.assembly  = chk[refname]

    def _set_project_files(self):
        '''Sets the project path, log, jobs, config, and results paths for the project.

           project_path - parent directory for this project
           config_path - parent directory for all config files within this project
           job_path - parent directory for all job files within this project
           log_path - parent directory for all log files within this project
           results_path - parent directory for all results files within this project
           report_path - parent directory for all report files within this project
        '''
        self.files = {
            'project_path': os.path.abspath(self.args.output_directory),
            'config_path': os.path.join(os.path.abspath(self.args.output_directory), 'config'),
            'log_path': os.path.join(os.path.abspath(self.args.output_directory), 'logs'),
            'results_path': os.path.join(os.path.abspath(self.args.output_directory), 'results'),
            'report_path': os.path.join(os.path.abspath(self.args.output_directory), 'report')
        }
        [check_dir(i) for i in self.files.values()]

    def __write_bds_config(self):
        '''Write the default BDS config file for an amazon cloud run'''
        self.bds_file = os.path.join(os.path.abspath(self.args.output_directory), 'bds.config')
        with open(self.bds_file, 'wb') as o:
            o.write('# BDS configuration file\n')
            o.write('# This is set up to work easily with STAR cluster on EC2\n')
            o.write('# Please see http://pcingola.github.io/BigDataScript/bigDataScript_manual.html#config\n')
            o.write('# for more information.\n\n')

            o.write('# Easy way to change mem usage for all jobs\n')
            o.write('#mem = -1\n\n')

            o.write('# Easy way to declare a particular node for all jobs\n')
            o.write('#node = my_node\n\n')

            o.write('# Easy way to declare a particular queue for all jobs\n')
            o.write('#queue = my_queue\n\n')

            o.write('# Easy way to declare the timeout for all jobs (seconds). It is 24 hours by default\n')
            o.write('#timeout = 86400\n\n')

            o.write('# Number of times a failed job can be retried. We find 0 is the safest.\n')
            o.write('retry = 0\n\n')

            o.write('# Sometimes many qsub commands is not able to be handled on some systems.\n')
            o.write('# Use this to set the number of milliseconds to wait after qsub.\n') 
            o.write('waitAfterTaskRun = 200\n\n')

            o.write('# Use bash shell for tasks\n') 
            o.write('taskShell = /bin/bash -e\n\n')

            o.write('# SGE parallel environment\n') 
            o.write('sge.pe = smp\n\n')

            o.write('# SGE mem argument. We set mem_free as consumable, jobs will wait, but it does not\n') 
            o.write('# guarantee that memory will be available (e.g., it will not kill jobs that use more\n')
            o.write('# memory than requested. You can change memory usage for each tool in the ExScalibur configs\n') 
            o.write('sge.mem = mem_free\n\n')

            o.write('# SGE timeout variable\n') 
            o.write('sge.timeout = h_rt\n\n')

    def __write_project_run_script(self):
        '''The script for running ExScalibur'''
        bds_script  = os.path.join(self.PATH, 'ExScaliburSMD-run.bds')
        main_job    = os.path.join(self.files['project_path'], '{0.name}.exscalibur-smd.sh'.format(self))
        bds_log     = os.path.join(self.files['project_path'], 'BDS-System.logs')
        cri_log     = os.path.join(self.files['project_path'], 'CRI-Info.logs')

        # Create the job script
        with open(main_job, 'wb') as o:
            o.write('#!/bin/bash\n')
            o.write('\n')
            o.write('cd {0}\n\n'.format(self.files['project_path']))

            if self.settings.system['module_source']:
                o.write('. {0}\n'.format(self.settings.system['module_source']))
            if self.settings.system['java']['use_module']:
                o.write('module load {0}\n\n'.format(self.settings.system['java']['module']['name']))
            o.write('# Command for running the ExScalibur-SMD pipeline\n\n')
            o.write('bds -c ' + self.bds_file + ' -reportHtml -reportYaml -v -log -s sge ' + \
                    '{0} -in {1} -options {2} > {3} 2> {4}\n'.format(
                bds_script, " ".join(self.config_list), self.settings.settings_cfg, 
                cri_log, bds_log))

        # Make user executable
        import stat
        st = os.stat(main_job)
        os.chmod(main_job, st.st_mode | stat.S_IEXEC)

    ###################################
    ## Reporter functions
    def __write_reporter_yaml(self):
        '''Writes the pipeline reporter yaml file'''
        self.final_reporter_yaml = os.path.join(self.files['report_path'], '{0.name}.exscalibur.smd.yaml'.format(self))
        dic = {'data': [], 'project': self.name}
        for s in self.samples:
            curr = s.reporter_obj
            dic['data'].append({
              'fastqc': self.__reporter_fastqc(curr),
              'alignments': self.__reporter_alignments(curr),
              'somatic': self.__reporter_somatic(curr),
              'sample': s.name
            })
        with open(self.final_reporter_yaml, 'wb') as o:
            yaml.safe_dump(dic, o)

    def __reporter_fastqc(self, dic):
        '''Creates the relative paths for fastqc files'''
        normal = dic['reads']['outputs']['normal']
        curr_normal = [{
            'readgroup': i['readgroup'],
            'leftseq': os.path.relpath(i['leftseq'], start=self.files['report_path']),
            'rightseq': os.path.relpath(i['rightseq'], start=self.files['report_path']) if i['paired'] else None,
            'paired': i['paired']
        } for i in normal]
        tumor = dic['reads']['outputs']['tumor']
        curr_tumor = [{
            'readgroup': i['readgroup'],
            'leftseq': os.path.relpath(i['leftseq'], start=self.files['report_path']),
            'rightseq': os.path.relpath(i['rightseq'], start=self.files['report_path']) if i['paired'] else None,
            'paired': i['paired']
        } for i in tumor]
        return {'inputs': {'normal': curr_normal, 'tumor': curr_tumor}}

    def __reporter_alignments(self, dic):
        '''Creates the relative paths for the alignment files'''
        curr = dic['alignments']['outputs']
        return {
            'inputs': {
                'alignment_summary_metrics': os.path.relpath(curr['alignment_summary_metrics'], 
                                                            start=self.files['report_path']),
                'insert_size_metrics': os.path.relpath(curr['insert_size_metrics'],
                                                            start=self.files['report_path']),
                'total_coverage': os.path.relpath(curr['total_coverage'], start=self.files['report_path'])
            }
        } 

    def __reporter_somatic(self, dic):
        '''Creates the relative paths for the somatic files'''
        return {
            'inputs': {
                'somatic_table': os.path.relpath(dic['somatic']['outputs']['smd_snp_table'], 
                                                 start=self.files['report_path'])
            }
        }


