'''Functions for filtering Mutect output VCF files'''
import logging

from smd_vcf import VcfReader, VcfRecord

logger = logging.getLogger('SomaticFilter.MuTect')

class MutectRecord(VcfRecord):
    '''Represents a VCF row from MuTect'''
    lookup = {0: 'WT', 1: 'GERM', 2: 'SOMATIC', 3: 'LOH', 4: 'PTM', 5: 'UK'}

    def apply_filters(self, args):
        '''Applies various filters to the VCF row with side-effects on the info
           and filt attributes -- see filters.smd_vcf.VcfRecord 
        '''
        # Initialize the filter list
        filter_list = [self.filt] if self.filt != 'PASS' else [] 

        # Create a list of the INFO column which will likely change
        info_dat    = [self.info] if self.info != "." else []

        # Create dicts of the tumor and normal columns
        tumor_dat   = dict(zip(self.fmt.split(":"), self.tumor.split(":")))
        normal_dat  = dict(zip(self.fmt.split(":"), self.normal.split(":")))
      
        ## Actual filters
        # Check for depth
        if int(normal_dat["DP"]) < args.min_normal_depth: filter_list.append("LowDPN")
        if int(tumor_dat["DP"]) < args.min_tumor_depth: filter_list.append("LowDPT")
        # Check for frequencies
        if float(tumor_dat["FA"]) < args.min_alt_freq_tumor: filter_list.append("TAF")
        if float(normal_dat["FA"]) >= args.max_alt_freq_normal: filter_list.append("NAF")
        # Check for base quality
        try:
            if int(tumor_dat["BQ"]) < args.min_base_quality: filter_list.append("LowBQ")
        except: pass
        # Check for the type of variant 
        if "SS" in tumor_dat:
            if int(tumor_dat["SS"]) == 5: filter_list.append("UK")
            curr = "TTYPE=" + self.lookup[int(tumor_dat["SS"])] 
            info_dat.append(curr)
        if "SS" in normal_dat:
            if int(normal_dat["SS"]) == 5: filter_list.append("UK")
            curr = "NTYPE=" + self.lookup[int(normal_dat["SS"])] 
            info_dat.append(curr)
        
        # Side effects
        self.info = "." if not info_dat else ";".join(info_dat)
        self.filt = "PASS" if not filter_list else ";".join(filter_list)

def run(args):
    '''Main wrapper function for filtering MuTect VCF files'''
    logger.info('<FILTER>LowDPN=Normal depth < {0.min_normal_depth}'.format(args))
    logger.info('<FILTER>LowDPT=Tumor depth < {0.min_tumor_depth}'.format(args))
    logger.info('<FILTER>TAF=(SAMPLE == TUMOR) && (TYPE==SOM) && AF < {0.min_alt_freq_tumor:.2f}'.format(args))
    logger.info('<FILTER>NAF=(SAMPLE == NORMAL) && (TYPE==SOM) && AF >= {0.max_alt_freq_normal:.2f}'.format(args))
    logger.info('<FILTER>LowBQ=BQ < {0.min_base_quality}'.format(args))

    # New filter and info lines to add to the output vcf
    filters = [
        '##FILTER=<ID=LowDPN,Description="Normal depth < {0.min_normal_depth}">'.format(args),
        '##FILTER=<ID=LowDPT,Description="Tumor depth < {0.min_tumor_depth}">'.format(args),
        '##FILTER=<ID=TAF,Description="Tumor: (TYPE == SOM) && AF < {0.min_alt_freq_tumor:.2f}">'.format(args),
        '##FILTER=<ID=NAF,Description="Normal: (TYPE == SOM) && AF >= {0.max_alt_freq_normal:.2f}">'.format(args),
        '##FILTER=<ID=LowBQ,Description="BQ < {0.min_base_quality}">'.format(args),
        '##FILTER=<ID=UK,Description="SS==5">'
    ]
    info = [
        '##INFO=<ID=NTYPE,Number=1,Type=String,Description="Normal type, can be WT, GERM, SOMATIC, LOH, PTM, UK">',
        '##INFO=<ID=TTYPE,Number=1,Type=String,Description="Tumor type, can be WT, GERM, SOMATIC, LOH, PTM, UK">'
    ]

    # Process the file
    with open(args.output_vcf, 'wb') as o:
        mutect_reader = VcfReader(args.choice, args.tumor_name, args.normal_name, args.input_vcf)
        mutect_reader.Open()
        mutect_reader.get_header()
        mutect_reader.write_new_header(o, filters=filters, info=info)
        mutect_reader.apply_filters(o, args, MutectRecord)
        mutect_reader.Close() 
