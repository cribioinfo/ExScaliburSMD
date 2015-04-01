'''Functions for filtering Shimmer output VCF files'''
import logging
import os

from smd_vcf import VcfReader, VcfRecord, load_contigs

logger = logging.getLogger('SomaticFilter.Shimmer')

class ShimmerRecord(VcfRecord):
    '''Represents a VCF row from Shimmer'''
    def apply_filters(self, args, vsdict):
        '''Applies various filters to the VCF row with side-effects on the filter 
           and format attributes
        '''
        # Initialize the filter list
        filter_list = [] 

        # Create a list of the FORMAT column which will likely change
        fmt_dat     = self.fmt.split(':') 

        # Create dicts of the tumor and normal columns
        tumor_dat   = dict(zip(fmt_dat, self.tumor.split(":")))
        normal_dat  = dict(zip(fmt_dat, self.normal.split(":")))

        # Get the data from the VarSifter dictionary
        key     = (self.chrom, self.pos)
        results = vsdict[key] 

        ## Filters
        # Depth
        if int(normal_dat['DP']) < args.min_normal_depth: filter_list.append("LowDPN") 
        if int(tumor_dat['DP']) < args.min_tumor_depth: filter_list.append("LowDPT") 

        # Frequency
        if float(results['tumor']) < args.min_alt_freq_tumor: filter_list.append("TAF")
        if float(results['normal']) >= args.max_alt_freq_normal: filter_list.append("NAF")

        # Quality
        if self.qual < args.min_qual: filter_list.append('LowQual')

        # Side effects
        fmt_dat.append('AF')
        self.fmt    = ':'.join(fmt_dat)
        self.filt   = "PASS" if not filter_list else ";".join(filter_list)
        self.normal = self.normal + ':' + '{0:.3f}'.format(results['normal'])
        self.tumor  = self.tumor + ':' + '{0:.3f}'.format(results['tumor'])

def load_varsifter(args):
    '''Loads the VarSifter file into a dictionary'''
    logger.info('Parsing the VarSifter file...')
    dic = {}
    for line in open(args.input_vs, 'rU'):
        if not line.startswith('Index'):
            cols = line.rstrip().split("\t")
            key  = (cols[1], int(cols[2]) + 1)
            nratio = float(cols[9])
            tratio = float(cols[10])
            dic[key] = {
                'normal': nratio,
                'tumor': tratio
            }
    return dic

def run(args):
    '''Main wrapper function for filtering Shimmer VCF files'''
    # Print the filters to the log
    logger.info('<FILTER>LowDPN=Normal DP < {0.min_normal_depth}'.format(args))
    logger.info('<FILTER>LowDPT=Tumor DP < {0.min_tumor_depth}'.format(args))
    logger.info('<FILTER>TAF=Tumor AF < {0.min_alt_freq_tumor:.3f}>'.format(args))
    logger.info('<FILTER>NAF=Normal AF >= {0.max_alt_freq_normal:.3f}>'.format(args))
    logger.info('<FILTER>LowQual=QUAL < {0.min_qual}>'.format(args))

    # New filter and format lines
    formats = ['##FORMAT=<ID=AF,Number=1,Type=Float,Description="Ratio of reads with alternate base">']
    filters = [
        '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">',
        '##FILTER=<ID=LowDPN,Description="Normal DP < {0.min_normal_depth}">'.format(args),
        '##FILTER=<ID=LowDPT,Description="Tumor DP < {0.min_tumor_depth}">'.format(args),
        '##FILTER=<ID=TAF,Description="Tumor AF < {0.min_alt_freq_tumor:.3f}">'.format(args),
        '##FILTER=<ID=NAF,Description="Normal AF >= {0.max_alt_freq_normal:.3f}">'.format(args),
        '##FILTER=<ID=LowQual,Description="QUAL < {0.min_qual}">'.format(args)
    ]

    # Load the VarSifter dictionary
    varsifter = load_varsifter(args)

    # Load the contigs
    reffile = '##reference=file://' + os.path.abspath(args.reference)
    contigs = load_contigs(args)

    # Process the file
    with open(args.output_vcf, 'wb') as o:
        shimmer_reader = VcfReader(args.choice, args.tumor_name, args.normal_name, args.input_vcf)
        shimmer_reader.Open()
        shimmer_reader.get_header()
        shimmer_reader.write_new_header(o, filters=filters, formats=formats,
                                        refpath=reffile, contigs=contigs)
        shimmer_reader.apply_filters(o, args, ShimmerRecord, vsdict=varsifter)
        shimmer_reader.Close() 
