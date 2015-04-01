'''Functions for filtering Strelka output VCF files'''
import logging
import os
import re

from smd_vcf import VcfReader, VcfRecord, load_contigs

logger = logging.getLogger('SomaticFilter.Strelka')

class StrelkaReader(object):
    '''Object that allows for the simultaneous processing of SNV and InDel VCFs'''
    pattern = re.compile('##contig=<ID=(.+),.+>')

    def __init__(self, args):
        self.indels      = VcfReader(args.choice, args.tumor_name, args.normal_name, args.input_all_indel)
        self.snvs        = VcfReader(args.choice, args.tumor_name, args.normal_name, args.input_all_snv)
        self.info        = []
        self.fmt         = []
        self.flt         = []
        self.contigs     = []
        self.other       = []
        self.chrom_order = []

    def Open(self):
        '''Opens both indel and snv files'''
        self.indels.Open()
        self.snvs.Open()

    def Close(self):
        '''Closes both indel and snv files'''
        self.indels.Close()
        self.snvs.Close()

    def get_headers(self):
        '''Loads both headers'''
        self.indels.get_header()
        self.snvs.get_header()

        # Combine info
        dic = {}
        for i in self.indels.header.info + self.snvs.header.info:
            k = i.split(',')[0]
            dic[k] = i
        [self.info.append(dic[i]) for i in dic]

        # Combine format 
        dic = {}
        for i in self.indels.header.fmt + self.snvs.header.fmt:
            k = i.split(',')[0]
            dic[k] = i
        [self.fmt.append(dic[i]) for i in dic]

        # Combine filter 
        dic = {}
        for i in self.indels.header.flt + self.snvs.header.flt:
            k = i.split(',')[0]
            dic[k] = i
        [self.flt.append(dic[i]) for i in dic]

        # Get contigs
        self.contigs = self.indels.header.contig

        # Get other info
        self.other = list(set(self.indels.header.other + self.snvs.header.other))

    def write_new_header(self, o, filters, info):
        # Write fmt, date, reference, and other
        if self.indels.header.vcffmt: o.write(self.indels.header.vcffmt + '\n')
        if self.indels.header.vcfdate: o.write(self.indels.header.vcfdate + '\n')
        if self.other: o.write('\n'.join(self.other) + '\n')
        if self.indels.header.ref: o.write(self.indels.header.ref + '\n') 
        
        # Write contigs
        if self.contigs: o.write('\n'.join(self.contigs) + '\n')

        # Write 
        # Write old and new info
        if self.info: o.write('\n'.join(self.info) + '\n')
        if info: o.write('\n'.join(info) + '\n')

        # Write old and new filters
        if self.flt: o.write('\n'.join(self.flt) + '\n')
        if filters: o.write('\n'.join(filters) + '\n')

        # Write formats
        if self.fmt: o.write('\n'.join(self.fmt) + '\n')

        # Write new header
        new_header = "\t".join(self.snvs.header.header[:9]) + '\t' + self.snvs.nname + '\t' + self.snvs.tname
        o.write(new_header + '\n')

    def apply_filters(self, args, o): 
        self.__get_chromosomes()
        
        # First, load indels into a dic
        dic = {}
        for line in self.indels.fh:
            indel_record = StrelkaIndelRecord(line.rstrip().split('\t'), 
                                              normal_idx=self.indels.nidx, tumor_idx=self.indels.tidx)
            indel_record.apply_filters(args)
            if indel_record.chrom not in dic: dic[indel_record.chrom] = {}
            dic[indel_record.chrom][indel_record.pos] = indel_record

        # Next, load snvs into the dic
        for line in self.snvs.fh:
            snv_record = StrelkaSnvRecord(line.rstrip().split('\t'), 
                                            normal_idx=self.snvs.nidx, tumor_idx=self.snvs.tidx)
            snv_record.apply_filters(args)
            if snv_record.chrom not in dic: dic[snv_record.chrom] = {}
            dic[snv_record.chrom][snv_record.pos] = snv_record

        # Now, write out the new VCF file with the correct order
        for c in self.chrom_order:
            if c in dic:
                for p in sorted(dic[c]):
                    dic[c][p].write_record(o)

    def __get_chromosomes(self):
        '''Parse out the chromosome order from the contig lines'''
        for i in self.contigs:
            self.chrom_order.append(self.pattern.match(i).groups()[0])

class StrelkaIndelRecord(VcfRecord):
    '''Represents a VCF row from Strelka'''
    def apply_filters(self, args):
        '''Applies various filters to the VCF row with side-effects on the filter 
           and format attributes
        '''
        # Initialize the filter list
        filter_list = [] if self.filt == 'PASS' else [self.filt] 

        # Create a list of the FORMAT column which will likely change
        fmt_dat     = self.fmt.split(':') 

        # Create dicts of the tumor and normal columns
        tumor_dat   = dict(zip(fmt_dat, self.tumor.split(":")))
        normal_dat  = dict(zip(fmt_dat, self.normal.split(":")))

        ## Filters
        # Depth
        if int(normal_dat['DP']) < args.min_normal_depth: filter_list.append("LowDPN")
        if int(tumor_dat['DP']) < args.min_tumor_depth: filter_list.append("LowDPT")

        # Side effects
        self.info   = self.info + ';INDEL'
        self.filt   = "PASS" if not filter_list else ";".join(filter_list)

class StrelkaSnvRecord(VcfRecord):
    '''Represents a VCF row from Strelka'''
    def apply_filters(self, args):
        '''Applies various filters to the VCF row with side-effects on the filter 
           and format attributes
        '''
        # Initialize the filter list
        filter_list = [] if self.filt == 'PASS' else [self.filt] 

        # Create a list of the FORMAT column which will likely change
        fmt_dat     = self.fmt.split(':') 

        # Create dicts of the tumor and normal columns
        tumor_dat   = dict(zip(fmt_dat, self.tumor.split(":")))
        normal_dat  = dict(zip(fmt_dat, self.normal.split(":")))

        ## Filters
        # Depth
        if int(normal_dat['DP']) < args.min_normal_depth: filter_list.append("LowDPN")
        if int(tumor_dat['DP']) < args.min_tumor_depth: filter_list.append("LowDPT")

        # Side effects
        self.filt   = "PASS" if not filter_list else ";".join(filter_list)

def run(args):
    '''Main wrapper function for filtering Strelka VCF files'''
    # Print the filters to the log
    logger.info('<FILTER>LowDPN=Normal DP < {0.min_normal_depth}'.format(args))
    logger.info('<FILTER>LowDPT=Tumor DP < {0.min_tumor_depth}'.format(args))

    ## New filter and info lines
    filters = [
        '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">',
        '##FILTER=<ID=LowDPN,Description="Normal DP < {0.min_normal_depth}">'.format(args),
        '##FILTER=<ID=LowDPT,Description="Tumor DP < {0.min_tumor_depth}">'.format(args)
    ]
    info   = [
        '##INFO=<ID=INDEL,Number=0,Type=Flag,Description="This is an InDel.">'
    ]

    # Process the files
    with open(args.output_vcf, 'wb') as o:
        strelka_reader = StrelkaReader(args)
        strelka_reader.Open()
        strelka_reader.get_headers()
        strelka_reader.write_new_header(o, filters, info)
        strelka_reader.apply_filters(args, o)
        strelka_reader.Close() 
