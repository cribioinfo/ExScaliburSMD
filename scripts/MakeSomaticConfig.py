#!/usr/bin/env python
"""Makes the default configuration YAML file"""
import ConfigParser
import argparse
import yaml
import os

# File to write
def load_system(config):
    '''Creates the system config object'''
    somatic_filter = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "SomaticFilter.py")
    return {
        'software': {
            'module_source' : '/etc/profile.d/modules.sh',
            'annovar'       : {'use_module': False, 
                               'module': {'name': 'annovar-bio', 'dir': '${ANNOVAR}'}, 
                               'dir': config['annovar_path'],
                               'db': config['annovar_db'], 'max_mem': 4},
            'bedtools'      : {'use_module': True, 'module': {'name': 'bedtools/2.21.0', 'exe': 'bedtools'}, 
                               'dir': None, 'max_mem': 1},
            'bwa'           : {'use_module': True, 'module': {'name': 'bwa/0.7.10', 'exe': 'bwa'}, 
                               'exe': None, 'threads': None, 'max_mem': 6},
            'cutadapt'      : {'use_module': False, 'module': {'name': None, 'exe': None}, 
                               'exe': 'cutadapt', 'max_mem': 1},
            'fastqc'        : {'use_module': True, 'module': {'name': 'fastqc/0.11.2', 'exe': 'fastqc'}, 
                               'exe': None, 'threads': None, 'max_mem': 2},
            'gatk'          : {'use_module': False, 'module': {'name': None, 'exe': '${GATK}'}, 
                               'exe': config['gatk_jar'], 
                               'max_mem': 8, 'threads': None},
            'java'          : {'use_module': True, 'module': {'name': 'java/1.7.0', 'exe': None}, 
                               'exe': None},
            'mutect'        : {'use_module': False, 'module': {'name': None, 'exe': '${MUTECT}'}, 
                               'exe': config['mutect_jar'], 
                               'max_mem': 1},
            'novoalign'     : {'use_module': False, 
                               'module': {'name': None, 'exe':'novoalign'}, 
                               'exe': config['novoalign'], 'threads': None, 'max_mem': 4},
            'picard'        : {'use_module': True, 
                               'module': {'name': 'picard/1.123', 'exe': '${PICARD}'},
                               'dir': None, 'max_merge_mem': 4, 'max_dedup_mem': 4, 'max_metrics_mem': 4,
                               'max_vcf_mem': 4},
            'perl'          : {'use_module': False, 'module': {'name': None, 'exe': 'perl'}, 
                               'dir': None},
            'python'        : {'use_module': False, 'module': {'name': None, 'exe': 'python'}, 
                               'exe': 'python', 'lib': None},
            'r'             : {'use_module': False, 'module': {'name': None, 'exe': None}, 'dir': None},
            'samtools'      : {'use_module': True, 'module': {'name': 'samtools/0.1.19', 'exe': 'samtools'}, 
                               'exe': None},
            'seqprep'       : {'use_module': True, 'module': {'name': 'SeqPrep/b5efabc5f7', 'exe': 'SeqPrep'}, 
                               'exe': None, 'max_mem': 1},
            'shimmer'       : {'use_module': True, 'module': {'name': 'shimmer/b62f433', 'exe': 'shimmer.pl'}, 
                               'exe': None, 'max_mem': 1}, 
            'sniper'        : {'use_module': True, 
                               'module': {'name': 'somatic-sniper/1.0.4', 'exe': 'bam-somaticsniper'}, 
                               'exe': None, 'max_mem': 1}, 
            'strelka'       : {'use_module': True, 'module': {'name': 'strelka/1.0.14', 
                               'exe': 'configureStrelkaWorkflow.pl'}, 'exe': None, 'threads': None, 'max_mem': 4},
            'varscan'       : {'use_module': True, 'module': {'name': 'varscan/2.3.6', 'exe': '${VARSCAN}'}, 
                               'exe': None, 'max_mem': 1, 'nscatter': None, 'max_mem': 1},
            'virmid'        : {'use_module': True, 'module': {'name': 'virmid/1.1.1', 'exe': '${VIRMID}'}, 
                               'exe': None, 'max_mem': 4}
            }
    }


## Program settings
def load_settings():
    '''Creates the software settings config object'''
    return {
        'annovar': {
            'refname': 'hg19',
            'annovardb': 'humandb',
            'standard': [
                {'protocol': 'refGene', 'operation': 'g'},
                {'protocol': 'snp138', 'operation': 'f'}, 
                {'protocol': 'esp6500si_all', 'operation': 'f'}, 
                {'protocol': 'ALL.sites.2014_10', 'operation': 'f'}, 
                {'protocol': 'cosmic70', 'operation': 'f'},
                #{'protocol': 'cadd', 'operation': 'f'}, 
                {'protocol': 'clinvar_20140929', 'operation': 'f'}, 
                {'protocol': 'ljb26_all', 'operation': 'f'}
            ],
            'reporter_colnames': [
                'Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
                'ExonicFunc.refGene','AAChange.refGene', 'snp138', 'esp6500si_all',
                'ALL.sites.2014_10', 'cosmic70', 'clinvar_20140929', 'CADD_raw'
            ],
            'other': None
        },
        'bwa_aln': {
            'barcode_length': 0,
            'read_trimming_parameter': 10,
            'max_occur': 100000
        },
        'bwa_mem': {
            'max_occur': 10000,
            'min_score': 30
        },
        'cutadapt': {
            'adapter1': None, 
            'adapter2': None 
        },
        'novoalign': {
            'clip_trailing': False,
            'repeats': 'Random',
            'max_multi': None
        },
        'varscan': {
            'min_freq_for_het': None, 
            'min_freq_for_hom': None, 
            'p_value_het': None, 
            'p_value_somatic': None, 
            'strand_filter': False,
            'somatic_filtration': {
                'min_normal_depth': None,
                'min_tumor_depth': None,
                'max_alt_freq_normal': None,
                'min_alt_freq_tumor': None,
                'pval_cutoff': None
            }
        },
        'virmid': {
            'max_training_sampling': None, 
            'max_read_depth': -1,
            'somatic_filtration': {
                'min_normal_depth': None,
                'min_tumor_depth': None,
                'max_alt_freq_normal': None,
                'min_alt_freq_tumor': None,
                'min_qual': None
            }
        },
        'strelka': {
            'max_spanning_deletion_fraction': 0.75,
            'indel_max_ref_repeat': 8,
            'snv_prior': 0.000001,
            'indel_prior': 0.000001,
            'min_snv_quality': 15,
            'min_indel_quality': 30,
            'ignore_conflicting_read_names': True,
            'skip_depth_filters': 1,
            'max_input_depth': -1,
            'ssnv_noise': 0.0000005,
            'sindel_noise': 0.000001,
            'base_config': os.path.join(PATH, 'strelka_config_bwa_default_template.ini'), 
            'somatic_filtration': {
                'min_normal_depth': None,
                'min_tumor_depth': None
            }
        },
        'seqprep': {
            'forward_adapter': None,
            'reverse_adapter': None
        },
        'sniper': {
            'min_somatic_quality': 15,
            'joint_genotyping_mode': False,
            'prior_probability': None,
            'n_haplotypes': None,
            'maq_theta': None, 
            'prior_haplotype_difference': None,
            'somatic_filtration': {
                'min_normal_depth': None,
                'min_tumor_depth': None,
                'min_mapq_normal': None,
                'min_mapq_tumor': None,
                'min_gq_normal': None,
                'min_gq_tumor': None,
                'min_somatic_score': None
            }
        },
        'shimmer': {
            'min_base_qual': 10,
            'max_qscore': 0.05,
            'somatic_filtration': {
                'min_normal_depth': None,
                'min_tumor_depth': None,
                'max_alt_freq_normal': None,
                'min_alt_freq_tumor': None,
                'min_qual': None
            }
        },
        'mutect': {
            'extended_output': True,
            'only_passing_calls': False,
            'tumor_lod_threshold': 6.3,
            'normal_lod_threshold': 2.2,
            'somatic_filtration': {
                'min_normal_depth': None,
                'min_tumor_depth': None,
                'max_alt_freq_normal': None,
                'min_alt_freq_tumor': None,
                'min_base_quality': None
            }
        } 
    }

def load_reference(config):
    '''Creates the reference config object'''
    return {'hg19': {
        'referenceseq': config['referenceseq'], 
        'sequence_dictionary': config['sequence_dictionary'], 
        'novoalign_index': config['novoalign_index'],
        'exons_bed': os.path.join(PATH, 'hg19.exons.refGene.bed'), 
        'knowndb': config['knowndb'], 
        'knownindel': config['knownindel'], 
        'cosmic': config['cosmic'] 
        }
    }

#def load_bds_settings():
#    return {
#        'system': 'sge',
#        'mem': 'free_mem'
#    }

def load_user_config(user_config):
    '''Loads the user config into a dictionary'''
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(user_config)

    cfg = {}

    # REFERENCE
    # REQUIRED
    assert config.has_option('REFERENCE', 'referenceseq') and config.get('REFERENCE', 'referenceseq'), \
        'Missing referenceseq' 
    cfg['referenceseq'] = os.path.abspath(config.get('REFERENCE', 'referenceseq'))

    assert config.has_option('REFERENCE', 'sequence_dictionary') and config.get('REFERENCE', 'sequence_dictionary'), \
        'Missing sequence_dictionary' 
    cfg['sequence_dictionary'] = os.path.abspath(config.get('REFERENCE', 'sequence_dictionary'))
 
    assert config.has_option('REFERENCE', 'knowndb') and config.get('REFERENCE', 'knowndb'), \
        'Missing knowndb' 
    cfg['knowndb'] = os.path.abspath(config.get('REFERENCE', 'knowndb'))

    assert config.has_option('REFERENCE', 'knownindel') and config.get('REFERENCE', 'knownindel'), \
        'Missing knownindel' 
    cfg['knownindel'] = os.path.abspath(config.get('REFERENCE', 'knownindel'))

    assert config.has_option('REFERENCE', 'cosmic') and config.get('REFERENCE', 'cosmic'), \
        'Missing cosmic' 
    cfg['cosmic'] = os.path.abspath(config.get('REFERENCE', 'cosmic'))

    # OPTIONAL
    if config.has_option('REFERENCE', 'novoalign_index') and config.get('REFERENCE', 'novoalign_index'):
        cfg['novoalign_index'] = os.path.abspath(config.get('REFERENCE', 'novoalign_index'))
    else: cfg['novoalign_index'] = None

    # SOFTWARE
    assert config.has_option('SOFTWARE', 'annovar_db') and config.get('SOFTWARE', 'annovar_db'), \
        'Missing annovar_db'
    cfg['annovar_db'] = os.path.abspath(config.get('SOFTWARE', 'annovar_db'))

    assert config.has_option('SOFTWARE', 'annovar_path') and config.get('SOFTWARE', 'annovar_path'), \
        'Missing annovar_path'
    cfg['annovar_path'] = os.path.abspath(config.get('SOFTWARE', 'annovar_path'))

    assert config.has_option('SOFTWARE', 'gatk_jar') and config.get('SOFTWARE', 'gatk_jar'), \
        'Missing gatk_jar'
    cfg['gatk_jar'] = os.path.abspath(config.get('SOFTWARE', 'gatk_jar'))

    assert config.has_option('SOFTWARE', 'mutect_jar') and config.get('SOFTWARE', 'mutect_jar'), \
        'Missing mutect_jar'
    cfg['mutect_jar'] = os.path.abspath(config.get('SOFTWARE', 'mutect_jar'))

    # OPTIONAL
    if config.has_option('SOFTWARE', 'novoalign') and config.get('SOFTWARE', 'novoalign'):
        cfg['novoalign'] = os.path.abspath(config.get('SOFTWARE', 'novoalign'))
    else: cfg['novoalign'] = None

    return cfg

if __name__ == '__main__':
    PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "data")
    p = argparse.ArgumentParser(prog='MakeSMDDefaultConfig', 
        epilog="Copyright 2014 - Center for Research Informatics - University of Chicago")

    p.add_argument('-o', '--output-file', required=True, 
        help='The path to the configuration file you want to create. Usually ending with .yml')

    p.add_argument('-c', '--user-config', required=True,
        help='The config file with the appropriate values for user-installed databases and software')
    args = p.parse_args()

    # Load config
    config_dic = load_user_config(args.user_config)

    dic = {}
    dic['reference'] = load_reference(config_dic)
    dic['system']    = load_system(config_dic)
    dic['settings']  = load_settings()
    #dic['bds']       = load_bds_settings()
    with open(args.output_file, 'wb') as o:
        yaml.safe_dump(dic, o, default_flow_style=False)
