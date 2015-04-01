'''Contains the object that parses all the settings files'''
import logging
import os
import yaml
import utils

class PipelineSettings(object):
    '''Holder and parser of all settings.

       -- args       - the argparse object
       -- PATH       - the top-level ExScaliburSMD path
       -- job_system - the pbs or sge queue system
       -- script     - full path to the ExScaliburSMD.py script
       -- threads    - global threads as provided by user
       -- is_split   - if this is a split project
       -- start_from - if the user needs to rerun a module from a sub module
       -- clipping   - are we clipping the reads?
       -- min_length - minimum read length when clipping
       -- min_map_q  - the minimum mapq for filtering
       -- aln_list   - list of user-chosen aligners
       -- smd_list   - list of user-chosen SMD tools
       -- system     - dictionary of the system config file
       -- options    - dictionary of the settings config file  
    '''
    def __init__(self, args, PATH):
        self.logger        = logging.getLogger('SMD.Settings')
        self.args          = args
        self.cfgfil        = yaml.safe_load(file(args.config_file))
        self.PATH          = PATH
        self.threads       = args.threads
        self.is_split      = False
        self.nsplits       = args.max_splits
        self.clipping      = args.no_adapter_clipping == False
        self.min_length    = args.fastq_min_length
        self.min_map_q     = args.min_map_q
        self.aln_list      = self._set_aln_list() 
        self.smd_list      = self._set_smd_list() 
        self.target_bed    = args.target_bed
        self._system       = None
        self._options      = None
        self.settings_cfg  = None

    @property
    def options(self):
        '''settings config file'''
        if not self._options:
            self._load_options_file()
        return self._options

    def _load_options_file(self):
        self._options = self.cfgfil['settings']

    @property
    def system(self):
        '''The dictionary of the system configuration yaml file'''
        if not self._system:
            self._set_system()
        return self._system

    def _set_system(self):
        '''Parses the system settings yaml file''' 
        self._system = self.cfgfil['system']['software'] 

    def _set_aln_list(self):
        '''Builds the list of requested aligners'''
        all_aln        = [self.args.bwa_aln, self.args.bwa_mem, self.args.novoalign]
        lst = []
        if all_aln == [0,0,0] or all_aln == [1,0,0]:
            lst.append('bwa_aln')
        else:
            if self.args.bwa_aln: lst.append('bwa_aln')
            if self.args.bwa_mem: lst.append('bwa_mem')
            if self.args.novoalign: lst.append('novoalign') 
        return lst

    def _set_smd_list(self):
        '''Builds the list of requested SMDs'''
        all_smd        = [self.args.mutect, self.args.shimmer, self.args.somatic_sniper,
                          self.args.strelka, self.args.varscan, self.args.virmid]
        lst = [] 
        if all_smd == [0, 0, 0, 0, 0, 0] or all_smd == [1, 0, 0, 0, 0, 0]:
            lst.append('mutect')
            if self.nsplits > 1: self.is_split = True
        else:
            if self.args.mutect: lst.append('mutect')
            if self.args.shimmer: lst.append('shimmer')
            if self.args.somatic_sniper: lst.append('sniper')
            if self.args.strelka: lst.append('strelka')
            if self.args.varscan: lst.append('varscan')
            if self.args.virmid: lst.append('virmid')

            if any([self.args.mutect, self.args.shimmer, self.args.somatic_sniper, self.args.varscan]) \
              and self.nsplits > 1:
                self.is_split = True
        return lst

    def write_settings_cfg(self, pobj):
        '''Creates the main system settings configuration file for BDS'''
        self.settings_cfg = os.path.join(pobj.files['config_path'], '{0.name}.settings.cfg'.format(pobj))
        
        with open(self.settings_cfg, 'wb') as o:
            # General
            o.write('# GENERAL SETTINGS\n')
            o.write('min_map_q = {0}\n'.format(self.min_map_q))
            o.write('min_clip_length = {0.min_length}\n'.format(self))
            if self.system['module_source']: o.write('module_source = {0}\n'.format(self.system['module_source']))
            o.write('\n')

            # Reference
            o.write('# REFERENCE SETTINGS\n')
            o.write('reference = {0}\n'.format(pobj.assembly['referenceseq']))
            o.write('reference_name = {0}\n'.format(pobj.assembly['refname'])) 
            o.write('exons_bed = {0}\n'.format(pobj.assembly['exons_bed'])) 
            o.write('reference_dict = {0}\n'.format(pobj.assembly['sequence_dictionary'])) 
            o.write('known_indel = {0}\n'.format(pobj.assembly['knownindel'])) 
            o.write('knowndb = {0}\n'.format(pobj.assembly['knowndb'])) 
            o.write('cosmic = {0}\n'.format(pobj.assembly['cosmic'])) 
            if pobj.assembly['novoalign_index']:
                utils.check_file(pobj.assembly['novoalign_index'], "Novoalign Index")
                o.write('novoalign_index = {0}\n'.format(pobj.assembly['novoalign_index']))
            if self.target_bed:
                utils.check_file(self.target_bed, "target bed file")
                o.write('target_bed = {0.target_bed}\n'.format(self))
            o.write('\n')

            # Preprocessing 
            self.__fastqc_settings(o)
            self.__cutadapt_settings(o)
            self.__seqprep_settings(o) 

            # Sam and bedtools
            self.__sam_bed_tools_settings(o)

            # Alignments
            if 'bwa_aln' in self.aln_list: self.__bwa_aln_settings(o)
            if 'bwa_mem' in self.aln_list: self.__bwa_mem_settings(o)
            if 'novoalign' in self.aln_list: self.__novoalign_settings(o)
            self.__picard_settings(o)

            # GATK
            self.__gatk_settings(o)

            # Mutation calling
            if 'mutect' in self.smd_list: self.__mutect_settings(o)
            if 'shimmer' in self.smd_list: self.__shimmer_settings(o)
            if 'sniper' in self.smd_list: self.__sniper_settings(o)
            if 'strelka' in self.smd_list: self.__strelka_settings(o)
            if 'varscan' in self.smd_list: self.__varscan_settings(o)
            if 'virmid' in self.smd_list: self.__virmid_settings(o)

            # Annovar
            self.__annovar_settings(o)

            # Others
            o.write('# OTHER SETTINGS\n')
            if self.system['python']['use_module']:
                o.write('python_modname = {0}\n'.format(self.system['python']['module']['name'])) 
                o.write('python_exe = {0}\n'.format(self.system['python']['module']['exe']))
            else:
                o.write('python_exe = {0}\n'.format(self.system['python']['exe']))
           
            if self.system['java']['use_module']:
                o.write('java_modname = {0}\n'.format(self.system['java']['module']['name'])) 
            o.write('somatic_filter_exe = {0}\n'.format(os.path.join(self.PATH, 'ExScaliburSMD-filter.py')))
            utils.check_file(os.path.join(self.PATH, 'ExScaliburSMD-filter.py'), 'ExScaliburSMD-filter.py')

            o.write('build_reporter_exe = {0}\n'.format(os.path.join(self.PATH, 'scripts/BuildReporterFiles.py')))
            utils.check_file(os.path.join(self.PATH, 'scripts/BuildReporterFiles.py'), 'BuildReporterFiles.py')

            o.write('fix_intersect_exe = {0}\n'.format(os.path.join(self.PATH, 'scripts/FixBedOutput.sh')))
            utils.check_file(os.path.join(self.PATH, 'scripts/FixBedOutput.sh'), 'FixBedOutput.sh')

            o.write('split_bed_exe = {0}\n\n'.format(os.path.join(self.PATH, 'scripts/CollapseSplitBed.py')))
            utils.check_file(os.path.join(self.PATH, 'scripts/CollapseSplitBed.py'), 'CollapseSplitBed.py')

    def __fastqc_settings(self, o):
        o.write('# FASTQC SETTINGS\n')
        if self.system['fastqc']['use_module']: 
            o.write('fastqc_modname = {0}\n'.format(self.system['fastqc']['module']['name'])) 
            o.write('fastqc_exe = {0}\n'.format(self.system['fastqc']['module']['exe']))
        else: 
            o.write('fastqc_exe = {0}\n'.format(self.system['fastqc']['exe']))
        o.write('fastqc_threads = {0}\n'.format(self.threads if \
          self.system['fastqc']['threads'] is None else self.system['fastqc']['threads']))
        o.write('fastqc_mem = {0}\n\n'.format(self.system['fastqc']['max_mem']))

    def __cutadapt_settings(self, o): 
        o.write('# CUTADAPT SETTINGS\n')
        if self.system['cutadapt']['use_module']: 
            o.write('cutadapt_modname = {0}\n'.format(self.system['cutadapt']['module']['name'])) 
            o.write('cutadapt_exe = {0}\n'.format(self.system['cutadapt']['module']['exe']))
        else: 
            o.write('cutadapt_exe = {0}\n'.format(self.system['cutadapt']['exe']))

        if self.options['cutadapt']['adapter1']:
            o.write('cutadapt_adapter1 = {0}\n'.format(self.options['cutadapt']['adapter1']))
        if self.options['cutadapt']['adapter2']:
            o.write('cutadapt_adapter2 = {0}\n'.format(self.options['cutadapt']['adapter2']))

        o.write('cutadapt_mem = {0}\n\n'.format(self.system['cutadapt']['max_mem']))
        
    def __seqprep_settings(self, o): 
        o.write('# SEQPREP SETTINGS\n')
        if self.system['seqprep']['use_module']: 
            o.write('seqprep_modname = {0}\n'.format(self.system['seqprep']['module']['name'])) 
            o.write('seqprep_exe = {0}\n'.format(self.system['seqprep']['module']['exe']))
        else: 
            o.write('seqprep_exe = {0}\n'.format(self.system['seqprep']['exe']))

        if self.options['seqprep']['forward_adapter']:
            o.write('seqprep_forward_adapter = {0}\n'.format(self.options['seqprep']['forward_adapter']))
        if self.options['seqprep']['reverse_adapter']:
            o.write('seqprep_reverse_adapter = {0}\n'.format(self.options['seqprep']['reverse_adapter']))

        o.write('seqprep_mem = {0}\n\n'.format(self.system['seqprep']['max_mem']))

    def __sam_bed_tools_settings(self, o):
        o.write('# SAMTOOLS SETTINGS\n')
        if self.system['samtools']['use_module']: 
            o.write('samtools_modname = {0}\n'.format(self.system['samtools']['module']['name'])) 
            o.write('samtools_exe = {0}\n\n'.format(self.system['samtools']['module']['exe']))
        else: 
            o.write('samtools_exe = {0}\n\n'.format(self.system['samtools']['exe']))

        o.write('# BEDTOOLS SETTINGS\n')
        if self.system['bedtools']['use_module']: 
            o.write('bedtools_modname = {0}\n'.format(self.system['bedtools']['module']['name'])) 
            o.write('bedtools_exe = bedtools\n\n')
        else: 
            o.write('bedtools_exe = {0}\n\n'.format(self.system['bedtools']['exe']))

    def __bwa_aln_settings(self, o):
        o.write('# BWA ALN SETTINGS\n')
        if self.system['bwa']['use_module']: 
            o.write('bwa_aln_modname = {0}\n'.format(self.system['bwa']['module']['name'])) 
            o.write('bwa_aln_exe = {0}\n'.format(self.system['bwa']['module']['exe']))
        else: 
            o.write('bwa_aln_exe = {0}\n'.format(self.system['bwa']['exe']))
        o.write('bwa_aln_threads = {0}\n'.format(self.threads if \
          self.system['bwa']['threads'] is None else self.system['bwa']['threads']))
        o.write('bwa_aln_mem = {0}\n'.format(self.system['bwa']['max_mem']))
        o.write('bwa_aln_barcode = {0}\n'.format(self.options['bwa_aln']['barcode_length']))
        o.write('bwa_aln_read_trim = {0}\n'.format(self.options['bwa_aln']['read_trimming_parameter']))
        o.write('bwa_aln_max_occur = {0}\n\n'.format(self.options['bwa_aln']['max_occur']))

    def __bwa_mem_settings(self, o):
        o.write('# BWA MEM SETTINGS\n')
        if self.system['bwa']['use_module']: 
            o.write('bwa_mem_modname = {0}\n'.format(self.system['bwa']['module']['name'])) 
            o.write('bwa_mem_exe = {0}\n'.format(self.system['bwa']['module']['exe']))
        else: 
            o.write('bwa_mem_exe = {0}\n'.format(self.system['bwa']['exe']))
        o.write('bwa_mem_threads = {0}\n'.format(self.threads if \
          self.system['bwa']['threads'] is None else self.system['bwa']['threads']))
        o.write('bwa_mem_mem = {0}\n'.format(self.system['bwa']['max_mem']))
        o.write('bwa_mem_max_occur = {0}\n\n'.format(self.options['bwa_mem']['max_occur']))

    def __novoalign_settings(self, o):
        o.write('# NOVOALIGN SETTINGS\n')
        if self.system['novoalign']['use_module']: 
            o.write('novoalign_modname = {0}\n'.format(self.system['novoalign']['module']['name'])) 
            o.write('novoalign_exe = {0}\n'.format(self.system['novoalign']['module']['exe']))
        else: 
            o.write('novoalign_exe = {0}\n'.format(self.system['novoalign']['exe']))
        o.write('novoalign_threads = {0}\n'.format(self.threads if \
          self.system['novoalign']['threads'] is None else self.system['novoalign']['threads']))
        o.write('novoalign_mem = {0}\n'.format(self.system['novoalign']['max_mem']))
        o.write('novoalign_clip_trailing = {0}\n'.format('true' \
          if self.options['novoalign']['clip_trailing'] else 'false'))
        o.write('novoalign_repeats = {0}\n'.format(self.options['novoalign']['repeats']))
        if self.options['novoalign']['max_multi']: 
            o.write('novoalign_max_multi = {0}\n\n'.format(self.options['novoalign']['max_multi']))
        else: o.write('\n')

    def __picard_settings(self, o):
        o.write('# PICARD SETTINGS\n')
        if self.system['picard']['use_module']: 
            o.write('picard_modname = {0}\n'.format(self.system['picard']['module']['name'])) 
            o.write('picard_exe = {0}\n'.format(self.system['picard']['module']['exe']))
        else: 
            o.write('picard_exe = {0}\n'.format(self.system['picard']['dir']))
        o.write('picard_merge_alignments_mem = {0}\n'.format(self.system['picard']['max_merge_mem']))
        o.write('picard_dedup_mem = {0}\n'.format(self.system['picard']['max_dedup_mem']))
        o.write('picard_metrics_mem = {0}\n'.format(self.system['picard']['max_metrics_mem']))
        o.write('picard_vcf_mem = {0}\n\n'.format(self.system['picard']['max_vcf_mem']))

    def __gatk_settings(self, o):
        o.write('# GATK SETTINGS\n')
        if self.system['gatk']['use_module']: 
            o.write('gatk_modname = {0}\n'.format(self.system['gatk']['module']['name'])) 
            o.write('gatk_exe = {0}\n'.format(self.system['gatk']['module']['exe']))
        else: 
            o.write('gatk_exe = {0}\n'.format(self.system['gatk']['exe']))
        o.write('gatk_threads = {0}\n'.format(self.threads if \
          self.system['gatk']['threads'] is None else self.system['gatk']['threads']))
        o.write('gatk_mem = {0}\n\n'.format(self.system['gatk']['max_mem']))

    def __mutect_settings(self, o):
        o.write('# MUTECT SETTINGS\n')
        if self.system['mutect']['use_module']: 
            o.write('mutect_modname = {0}\n'.format(self.system['mutect']['module']['name'])) 
            o.write('mutect_exe = {0}\n'.format(self.system['mutect']['module']['exe']))
        else: 
            o.write('mutect_exe = {0}\n'.format(self.system['mutect']['exe']))
        o.write('mutect_mem = {0}\n'.format(self.system['mutect']['max_mem']))
        o.write('mutect_extended_output = {0}\n'.format('true' if \
          self.options['mutect']['extended_output'] else 'false'))
        o.write('mutect_only_passing_calls = {0}\n'.format('true' if \
          self.options['mutect']['only_passing_calls'] else 'false'))
        if self.options['mutect']['tumor_lod_threshold'] is not None:
            o.write('mutect_tumor_lod_threshold = {0:.2f}\n'.format(\
              self.options['mutect']['tumor_lod_threshold']))
        if self.options['mutect']['normal_lod_threshold'] is not None:
            o.write('mutect_normal_lod_threshold = {0:.2f}\n'.format(\
              self.options['mutect']['normal_lod_threshold']))
        for flt in ['min_normal_depth', 'min_tumor_depth', 'max_alt_freq_normal', \
                    'min_alt_freq_tumor', 'min_base_quality']:
            if self.options['mutect']['somatic_filtration'][flt]:
                curr = str(self.options['mutect']['somatic_filtration'][flt])
                o.write('mutect_filter_{0} = {1}\n'.format(flt, curr))
        o.write('\n')

    def __shimmer_settings(self, o):
        o.write('# SHIMMER SETTINGS\n')
        if self.system['shimmer']['use_module']: 
            o.write('shimmer_modname = {0}\n'.format(self.system['shimmer']['module']['name'])) 
            o.write('shimmer_exe = {0}\n'.format(self.system['shimmer']['module']['exe']))
        else: 
            o.write('shimmer_exe = {0}\n'.format(os.path.basename(self.system['shimmer']['exe'])))
            o.write('shimmer_path = {0}\n'.format(os.path.abspath(os.path.dirname(self.system['shimmer']['exe']))))
        o.write('shimmer_mem = {0}\n'.format(self.system['shimmer']['max_mem']))
        if self.options['shimmer']['min_base_qual']: 
            o.write('shimmer_min_bq = {0}\n'.format(self.options['shimmer']['min_base_qual']))
        if self.options['shimmer']['max_qscore']: 
            o.write('shimmer_max_qscore = {0:.2f}\n'.format(self.options['shimmer']['max_qscore']))

        for flt in ['min_normal_depth', 'min_tumor_depth', 'max_alt_freq_normal', \
                    'min_alt_freq_tumor', 'min_qual']:
            if self.options['shimmer']['somatic_filtration'][flt]:
                curr = str(self.options['shimmer']['somatic_filtration'][flt])
                o.write('shimmer_filter_{0} = {1}\n'.format(flt, curr))
        o.write('\n')

    def __sniper_settings(self, o):
        o.write('# SOMATIC SNIPER SETTINGS\n')
        if self.system['sniper']['use_module']: 
            o.write('sniper_modname = {0}\n'.format(self.system['sniper']['module']['name'])) 
            o.write('sniper_exe = {0}\n'.format(self.system['sniper']['module']['exe']))
        else: 
            o.write('sniper_exe = {0}\n'.format(self.system['sniper']['exe']))
        o.write('sniper_mem = {0}\n'.format(self.system['sniper']['max_mem']))
        if self.options['sniper']['min_somatic_quality']: 
            o.write('sniper_min_somatic_quality = {0}\n'.format(self.options['sniper']['min_somatic_quality']))
        if self.options['sniper']['joint_genotyping_mode']: 
            o.write('sniper_use_priors = true\n')
        if self.options['sniper']['prior_probability']: 
            o.write('sniper_prior_prob = {0}\n'.format(str(self.options['sniper']['prior_probability'])))
        if self.options['sniper']['n_haplotypes']: 
            o.write('sniper_n_hap = {0}\n'.format(self.options['sniper']['n_haplotypes']))
        if self.options['sniper']['maq_theta']: 
            o.write('sniper_mac_theta = {0}\n'.format(str(self.options['sniper']['maq_theta'])))
        if self.options['sniper']['prior_haplotype_difference']: 
            o.write('sniper_hap_prior = {0}\n'.format(str(self.options['sniper']['prior_haplotype_difference'])))

        for flt in ['min_normal_depth', 'min_tumor_depth', 'min_mapq_normal', \
                    'min_mapq_tumor', 'min_gq_normal', 'min_gq_tumor', 'min_somatic_score']:
            if self.options['sniper']['somatic_filtration'][flt]:
                curr = str(self.options['sniper']['somatic_filtration'][flt])
                o.write('sniper_filter_{0} = {1}\n'.format(flt, curr))
        o.write('\n')

    def __strelka_settings(self, o):
        o.write('# STRELKA SETTINGS\n')
        if self.system['strelka']['use_module']: 
            o.write('strelka_modname = {0}\n'.format(self.system['strelka']['module']['name'])) 
            o.write('strelka_exe = {0}\n'.format(self.system['strelka']['module']['exe']))
        else: 
            o.write('strelka_exe = {0}\n'.format(self.system['strelka']['exe']))
        o.write('strelka_mem = {0}\n'.format(self.system['strelka']['max_mem']))
        o.write('strelka_threads = {0}\n'.format(self.threads if \
          self.system['strelka']['threads'] is None else self.system['strelka']['threads']))

        for flt in ['min_normal_depth', 'min_tumor_depth']: 
            if self.options['strelka']['somatic_filtration'][flt]:
                curr = str(self.options['strelka']['somatic_filtration'][flt])
                o.write('strelka_filter_{0} = {1}\n'.format(flt, curr))
        o.write('\n')

    def __varscan_settings(self, o):
        o.write('# VARSCAN2 SETTINGS\n')
        if self.system['varscan']['use_module']: 
            o.write('varscan_modname = {0}\n'.format(self.system['varscan']['module']['name'])) 
            o.write('varscan_exe = {0}\n'.format(self.system['varscan']['module']['exe']))
        else: 
            o.write('varscan_exe = {0}\n'.format(self.system['varscan']['exe']))
        o.write('varscan_mem = {0}\n'.format(self.system['varscan']['max_mem']))
        if self.options['varscan']['min_freq_for_het']: 
            o.write('varscan_min_freq_het = {0}\n'.format(str(self.options['varscan']['min_freq_for_het'])))
        if self.options['varscan']['min_freq_for_hom']: 
            o.write('varscan_min_freq_hom = {0}\n'.format(str(self.options['varscan']['min_freq_for_hom'])))
        if self.options['varscan']['p_value_het']: 
            o.write('varscan_p_value_het = {0}\n'.format(str(self.options['varscan']['p_value_het'])))
        if self.options['varscan']['p_value_somatic']: 
            o.write('varscan_p_value_hom = {0}\n'.format(str(self.options['varscan']['p_value_somatic'])))
        o.write('varscan_strand_filter = {0}\n'.format(
          'true' if self.options['varscan']['strand_filter'] else 'false'))

        for flt in ['min_normal_depth', 'min_tumor_depth', 'max_alt_freq_normal','min_alt_freq_tumor', 'pval_cutoff']: 
            if self.options['varscan']['somatic_filtration'][flt]:
                curr = str(self.options['varscan']['somatic_filtration'][flt])
                o.write('varscan_filter_{0} = {1}\n'.format(flt, curr))
        o.write('\n')

    def __virmid_settings(self, o):
        o.write('# VIRMID SETTINGS\n')
        if self.system['virmid']['use_module']: 
            o.write('virmid_modname = {0}\n'.format(self.system['virmid']['module']['name']))
            o.write('virmid_exe = {0}\n'.format(self.system['virmid']['module']['exe']))
        else: 
            o.write('virmid_exe = {0}\n'.format(self.system['virmid']['exe']))
        o.write('virmid_mem = {0}\n'.format(self.system['virmid']['max_mem']))
        if self.options['virmid']['max_training_sampling']: 
            o.write('virmid_max_training_sampling = {0}\n'.format(self.options['virmid']['max_training_sampling']))
        o.write('virmid_max_dp = {0}\n'.format(self.options['virmid']['max_read_depth']))

        for flt in ['min_normal_depth', 'min_tumor_depth', 'max_alt_freq_normal', \
                    'min_alt_freq_tumor', 'min_qual']: 
            if self.options['virmid']['somatic_filtration'][flt]:
                curr = str(self.options['virmid']['somatic_filtration'][flt])
                o.write('virmid_filter_{0} = {1}\n'.format(flt, curr))
        o.write('\n')

    def __annovar_settings(self, o):
        o.write('# ANNOVAR SETTINGS\n')
        dbname = self.options['annovar']['annovardb']
        if self.system['annovar']['use_module']: 
            o.write('annovar_modname = {0}\n'.format(self.system['annovar']['module']['name']))
            o.write('annovar_convert = {0}\n'.format(
              os.path.join(self.system['annovar']['module']['exe'], 'convert2annovar.pl')))
            o.write('annovar_table = {0}\n'.format(
              os.path.join(self.system['annovar']['module']['exe'], 'table_annovar.pl')))
            o.write('annovar_db_path = {0}\n'.format(os.path.join(self.system['annovar']['module']['exe'], dbname)))
        else: 
            o.write('annovar_convert = {0}\n'.format(
              os.path.join(self.system['annovar']['dir'], 'convert2annovar.pl')))
            o.write('annovar_table = {0}\n'.format(
              os.path.join(self.system['annovar']['dir'], 'table_annovar.pl')))
            o.write('annovar_db_path = {0}\n'.format(
              os.path.join(self.system['annovar']['dir'], dbname)))
        o.write('annovar_db_name = {0}\n'.format(dbname))

        anno_ops = []
        anno_pro = []
        for a in self.options['annovar']['standard']:
            anno_pro.append(a['protocol'])
            anno_ops.append(a['operation'])

        if self.options['annovar']['other']:
            for b in self.options['annovar']['standard']:
                anno_pro.append(b['protocol'])
                anno_ops.append(b['operation'])
        
        o.write('annovar_protocols = {0}\n'.format(",".join(anno_pro))) 
        o.write('annovar_operations = {0}\n'.format(",".join(anno_ops)))
        o.write('annovar_mem = {0}\n\n'.format(self.system['annovar']['max_mem']))
