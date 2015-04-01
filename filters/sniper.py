'''Functions for filtering Shimmer output VCF files'''
import logging
import os

from smd_vcf import VcfReader, VcfRecord, load_contigs

logger = logging.getLogger('SomaticFilter.SomaticSniper')

class SniperRecord(VcfRecord):
    '''Represents a VCF row from SomaticSniper'''
    lookup = {0: "REF", 1: "GERM", 2: "SOMATIC", 3: "LOH", 4: "UK"}

    def apply_filters(self, args):
        '''Applies various filters to the VCF row with side-effects on the filter 
           and format attributes
        '''
        # Initialize the filter list
        filter_list = []

        # Get info column which will likey change
        info_dat    = [self.info] if self.info != '.' else [] 

        # Create a list of the FORMAT column
        fmt_dat     = self.fmt.split(':')

        # Create dicts of the tumor and normal columns
        tumor_dat   = dict(zip(fmt_dat, self.tumor.split(":")))
        normal_dat  = dict(zip(fmt_dat, self.normal.split(":")))

        ## Filters
        # Depth
        if int(normal_dat['DP']) < args.min_normal_depth: filter_list.append("LowDPN")
        if int(tumor_dat['DP']) < args.min_tumor_depth: filter_list.append("LowDPT")

        # Genotype Quality 
        if int(normal_dat['GQ']) < args.min_gq_normal: filter_list.append("LowGQN")
        if int(tumor_dat['GQ']) < args.min_gq_tumor: filter_list.append("LowGQT")

        # Mapping Quality 
        if int(normal_dat['MQ']) < args.min_mapq_normal: filter_list.append("LowMQN")
        if int(tumor_dat['MQ']) < args.min_mapq_tumor: filter_list.append("LowMQT")

        # Somatic Score 
        if int(tumor_dat['SSC']) < args.min_somatic_score: filter_list.append("LowScore")

        # Check for the type of variant 
        if "SS" in tumor_dat:
            if int(tumor_dat["SS"]) == 4: filter_list.append("UK")
            curr = "TTYPE=" + self.lookup[int(tumor_dat["SS"])]
            info_dat.append(curr)
        if "SS" in normal_dat:
            if int(normal_dat["SS"]) == 4: filter_list.append("UK")
            curr = "NTYPE=" + self.lookup[int(normal_dat["SS"])]
            info_dat.append(curr)

        # Side effects
        self.filt   = "PASS" if not filter_list else ";".join(filter_list)
        self.info   = ";".join(info_dat)

    def write_record(self, o):
        o.write("{0.chrom}\t{0.pos}\t{0.dbid}\t{0.ref}\t{0.alt}\t{0.qual}\t{0.filt}\t".format(self) + \
                "{0.info}\t{0.fmt}\t{0.normal}\t{0.tumor}\n".format(self))

def run(args):
    '''Main wrapper function for filtering SomaticSniper VCF files'''
    # Print the filters to the log
    logger.info('<FILTER>LowDPN="Normal DP < {0.min_normal_depth}"'.format(args))
    logger.info('<FILTER>LowDPT="Tumor DP < {0.min_tumor_depth}"'.format(args))
    logger.info('<FILTER>LowMQT="Tumor MQ < {0.min_mapq_tumor}"'.format(args))
    logger.info('<FILTER>LowMQN="Normal MQ < {0.min_mapq_normal}"'.format(args))
    logger.info('<FILTER>LowGQT="Tumor GQ < {0.min_gq_tumor}"'.format(args))
    logger.info('<FILTER>LowGQN="Normal GQ < {0.min_gq_normal}"'.format(args))
    logger.info('<FILTER>LowScore="Somatic score < {0.min_somatic_score}"'.format(args))

    # New info and filter lines
    info    = [
        '##INFO=<ID=NTYPE,Number=1,Type=String,Description="Normal type, can be REF,GERM,SOMATIC,LOH,UK">',
        '##INFO=<ID=TTYPE,Number=1,Type=String,Description="Tumor type REF,GERM,SOMATIC,LOH,UK">'
    ]
    filters = [
        '##FILTER=<ID=PASS,Description="Accept as a confident somatic mutation">',
        '##FILTER=<ID=LowDPN,Description="Normal DP < {0.min_normal_depth}">'.format(args),
        '##FILTER=<ID=LowDPT,Description="Tumor DP < {0.min_tumor_depth}">'.format(args),
        '##FILTER=<ID=LowMQT,Description="Tumor MQ < {0.min_mapq_tumor}">'.format(args),
        '##FILTER=<ID=LowMQN,Description="Normal MQ < {0.min_mapq_normal}">'.format(args),
        '##FILTER=<ID=LowGQT,Description="Tumor GQ < {0.min_gq_tumor}">'.format(args),
        '##FILTER=<ID=LowGQN,Description="Normal GQ < {0.min_gq_normal}">'.format(args),
        '##FILTER=<ID=LowScore,Description="Somatic score < {0.min_somatic_score}">'.format(args),
        '##FILTER=<ID=UK,Description="Unknown variant type">'
    ]

    # Load the contigs
    contigs = load_contigs(args)

    # Process the file
    with open(args.output_vcf, 'wb') as o:
        sniper_reader = VcfReader(args.choice, args.tumor_name, args.normal_name, args.input_vcf)
        sniper_reader.Open()
        sniper_reader.get_header()
        sniper_reader.write_new_header(o, filters=filters, info=info, contigs=contigs)
        sniper_reader.apply_filters(o, args, SniperRecord)
        sniper_reader.Close()
    logger.info('Filtered and formatted VCF file: {0}'.format(os.path.abspath(args.output_vcf)))
