'''Functions for filtering virmid output VCF files'''
import logging
import os
import re

from smd_vcf import VcfReader, VcfRecord, load_contigs

logger = logging.getLogger('SomaticFilter.Virmid')

class VirmidReader(object):
    '''Object that allows for the simultaneous processing of germ, loh, and som VCFs'''
    pattern = re.compile('##contig=<ID=(.+),.+>')

    def __init__(self, args):
        self.germline    = VcfReader(args.choice, args.tumor_name, args.normal_name, args.input_all_germline)
        self.loh         = VcfReader(args.choice, args.tumor_name, args.normal_name, args.input_all_loh)
        self.somatic     = VcfReader(args.choice, args.tumor_name, args.normal_name, args.input_all_somatic)
        self.info        = []
        self.flt         = []
        self.contigs     = load_contigs(args) 
        self.other       = []
        self.chrom_order = []
        self.tumor_name  = args.tumor_name
        self.normal_name = args.normal_name

    def Open(self):
        '''Opens germline, loh, and somatic files'''
        self.germline.Open()
        self.loh.Open()
        self.somatic.Open()

    def Close(self):
        '''Closes germline, loh, and somatic files'''
        self.germline.Close()
        self.loh.Close()
        self.somatic.Close()

    def get_headers(self):
        '''Loads all headers'''
        self.germline.get_header()
        self.loh.get_header()
        self.somatic.get_header()

        # Combine info
        dic = {}
        for i in self.germline.header.info + self.loh.header.info + self.somatic.header.info:
            k = i.split(',')[0]
            dic[k] = i
        [self.info.append(dic[i]) for i in dic]

        # Combine filter 
        dic = {}
        for i in self.germline.header.flt + self.loh.header.flt + self.somatic.header.flt:
            k = i.split(',')[0]
            dic[k] = i
        [self.flt.append(dic[i]) for i in dic]

        # Get other info
        self.other = list(set(self.germline.header.other + self.loh.header.other + self.somatic.header.other))
        self.other.append('##normalSampleName={0.normal_name}'.format(self))
        self.other.append('##tumorSampleName={0.tumor_name}'.format(self))

    def write_new_header(self, o, filters, info):
        '''Writes the new header lines to the output VCF'''
        # Write fmt, date, reference, and other
        if self.germline.header.vcffmt: o.write(self.germline.header.vcffmt + '\n')
        if self.germline.header.vcfdate: o.write(self.germline.header.vcfdate + '\n')
        if self.other: o.write('\n'.join(self.other) + '\n')
        if self.germline.header.src: o.write(self.germline.header.src + '\n')
        if self.germline.header.ref: o.write(self.germline.header.ref + '\n') 
        
        # Write contigs
        if self.contigs: o.write('\n'.join(self.contigs) + '\n')

        # Write old and new info
        if self.info: o.write('\n'.join(self.info) + '\n')
        if info: o.write('\n'.join(info) + '\n')

        # Write old and new filters
        if self.flt: o.write('\n'.join(self.flt) + '\n')
        if filters: o.write('\n'.join(filters) + '\n')

        # Write new header
        new_header = "\t".join(self.somatic.header.header[:9])
        o.write(new_header + '\n')

    def apply_filters(self, args, o): 
        '''Wrapper function for applying filters to germline, LOH, and somatic VCFs'''
        self.__get_chromosomes()
        
        # First, load germline into a dic
        dic = {}
        for line in self.germline.fh:
            germline_record = VirmidGermlineRecord(line.rstrip().split('\t'))
            germline_record.apply_filters(args)
            if germline_record.chrom not in dic: dic[germline_record.chrom] = {}
            dic[germline_record.chrom][germline_record.pos] = germline_record

        # Next, load LOH into the dic
        for line in self.loh.fh:
            loh_record = VirmidLohRecord(line.rstrip().split('\t')) 
            loh_record.apply_filters(args)
            if loh_record.chrom not in dic: dic[loh_record.chrom] = {}
            dic[loh_record.chrom][loh_record.pos] = loh_record

        # Finally, load somatic into the dic
        for line in self.somatic.fh:
            somatic_record = VirmidSomaticRecord(line.rstrip().split('\t')) 
            somatic_record.apply_filters(args)
            if somatic_record.chrom not in dic: dic[somatic_record.chrom] = {}
            dic[somatic_record.chrom][somatic_record.pos] = somatic_record

        # Now, write out the new VCF file with the correct order
        for c in self.chrom_order:
            if c in dic:
                for p in sorted(dic[c]):
                    dic[c][p].write_record(o)

    def __get_chromosomes(self):
        '''Parse out the chromosome order from the contig lines'''
        for i in self.contigs:
            self.chrom_order.append(self.pattern.match(i).groups()[0])

class VirmidGermlineRecord(VcfRecord):
    '''Represents a VCF row from Virmid germline'''
    def apply_filters(self, args):
        '''Applies various filters to the VCF row with side-effects on the filter 
           and format attributes
        '''
        # Initialize the filter list
        filter_list = [] if self.filt == 'PASS' else [self.filt] 

        # Parse info column
        info = {j[0] : j[1] for j in [i.split('=') for i in self.info.split(';')]}

        ## Filters
        # Depth
        if int(info['NDP']) < args.min_normal_depth: filter_list.append("LowDPN")
        if int(info['DDP']) < args.min_tumor_depth: filter_list.append("LowDPT")

        # Side effects
        self.info   = 'TYPE=GERM;' + self.info
        self.filt   = "PASS" if not filter_list else ";".join(filter_list)

class VirmidLohRecord(VcfRecord):
    '''Represents a VCF row from Virmid loh'''
    def apply_filters(self, args):
        '''Applies various filters to the VCF row with side-effects on the filter 
           and format attributes
        '''
        # Initialize the filter list
        filter_list = [] if self.filt == 'PASS' else [self.filt] 

        # Parse info column
        info = {j[0] : j[1] for j in [i.split('=') for i in self.info.split(';')]}

        ## Filters
        # Depth
        if int(info['NDP']) < args.min_normal_depth: filter_list.append("LowDPN")
        if int(info['DDP']) < args.min_tumor_depth: filter_list.append("LowDPT")

        # QUAL
        if self.qual < args.min_qual: filter_list.append("LowQual")
        
        # Side effects
        self.info   = 'TYPE=LOH;' + self.info
        self.filt   = "PASS" if not filter_list else ";".join(filter_list)

class VirmidSomaticRecord(VcfRecord):
    '''Represents a VCF row from Virmid somatic'''
    def apply_filters(self, args):
        '''Applies various filters to the VCF row with side-effects on the filter 
           and format attributes
        '''
        # Initialize the filter list
        filter_list = [] if self.filt == 'PASS' else [self.filt] 

        # Parse info column
        info = {j[0] : j[1] for j in [i.split('=') for i in self.info.split(';')]}
        naf, daf = None, None

        ## Filters
        # Depth
        if int(info['NDP']) < args.min_normal_depth: filter_list.append("LowDPN")
        if int(info['DDP']) < args.min_tumor_depth: filter_list.append("LowDPT")

        # Frequency
        try: 
            if float(info['naf']) >= args.max_alt_freq_normal: filter_list.append('NAF')
            if float(info['daf']) < args.min_alt_freq_tumor: filter_list.append('TAF')
        except KeyError:
            naf = int(info['NAC']) / float(info['NDP'])
            daf = int(info['DAC']) / float(info['DDP'])
            if naf >= args.max_alt_freq_normal: filter_list.append('NAF')
            if daf < args.min_alt_freq_tumor: filter_list.append('TAF')

        # QUAL
        if self.qual < args.min_qual: filter_list.append("LowQual")
        
        # Side effects
        self.info   = 'TYPE=SOMATIC;' + self.info if naf is None else \
                      'TYPE=SOMATIC;naf={0:.2f};daf={1:.2f};'.format(naf,daf) + self.info 
        self.filt   = "PASS" if not filter_list else ";".join(filter_list)

def run(args):
    '''Main wrapper function for filtering Virmid VCF files'''

    # Print the filters to the log
    logger.info('<FILTER>LowDPN=NDP < {0.min_normal_depth}'.format(args))
    logger.info('<FILTER>LowDPT=DDP < {0.min_tumor_depth}'.format(args))
    logger.info('<FILTER>NAF=(TYPE == SOMATIC) && (naf > {0.max_alt_freq_normal:.2f})'.format(args))
    logger.info('<FILTER>TAF=(TYPE == SOMATIC) && (daf < {0.min_alt_freq_tumor:.2f}'.format(args))
    logger.info('<FILTER>LowQual=QUAL < {0.min_qual}'.format(args))

    ## New filter and info lines
    filters = [
        '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">',
        '##FILTER=<ID=LowDPN,Description="NDP < {0.min_normal_depth}">'.format(args),
        '##FILTER=<ID=LowDPT,Description="DDP < {0.min_tumor_depth}">'.format(args),
        '##FILTER=<ID=NAF,Description="(TYPE == SOMATIC) && (naf > {0.max_alt_freq_normal:.2f})">'.format(args),
        '##FILTER=<ID=TAF,Description="(TYPE == SOMATIC) && (daf < {0.min_alt_freq_tumor:.2f})">'.format(args),
        '##FILTER=<ID=LowQual,Description="QUAL < {0.min_qual}">'.format(args)
    ]
    info   = [
        '##INFO=<ID=TYPE,Number=A,Type=String,Description="GERM - A potential germline mutation; ' + \
        'LOH - A potential loss of heterozygosity event; SOMATIC - potential somatic mutation">',
        '##INFO=<ID=naf,Number=1,Type=Float,Description="Normal alt allele frequency.">',
        '##INFO=<ID=daf,Number=1,Type=Float,Description="Disease alt allele frequency.">'
    ]

    # Process the files
    with open(args.output_vcf, 'wb') as o:
        virmid_reader = VirmidReader(args)
        virmid_reader.Open()
        virmid_reader.get_headers()
        virmid_reader.write_new_header(o, filters, info)
        virmid_reader.apply_filters(args, o)
        virmid_reader.Close()
