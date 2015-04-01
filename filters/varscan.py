'''Functions for filtering VarScan2 output VCF files'''
import logging
import os
import re
import datetime

from smd_vcf import load_contigs

logger = logging.getLogger('SomaticFilter.VarScan2')

# IUPAC lookup dictionary
iupac  = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': ['A','G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
    'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G']
}

# Lookup dict
lookup = {"Germline": "GERM", "Somatic": "SOMATIC", "LOH": "LOH"}

# Parse out contigs
pattern = re.compile('##contig=<ID=(.+),.+>')

def load_snv(args, dic):
    '''Loads the SNP VarScan format into a dictionary, applies filters,
       and converts to VCF format'''
    head = []
    linect = 0
    for line in open(args.input_snv, 'rU'):
        linect += 1
        if not head: head = line.rstrip().split("\t")
        else:
            try:
                dat = dict(zip(head, line.rstrip('\n').split("\t")))
                ss  = dat["somatic_status"]
                if ss not in lookup: continue
                ss = lookup[ss]
                c = dat["chrom"]
                p = int(dat["position"])
                if c not in dic: dic[c] = {}
                vcf_row = [dat["chrom"], dat['position'], ".", dat['ref'], dat['var'], "."]
                filter_list = []

                nr1 = int(dat["normal_reads1"])
                nr2 = int(dat["normal_reads2"])
                tr1 = int(dat["tumor_reads1"])
                tr2 = int(dat["tumor_reads2"])
                sspv = float(dat["somatic_p_value"])

                ncov = nr1 + nr2
                tcov = tr1 + tr2
                nfreq = nr2 / float(ncov)
                tfreq = tr2 / float(tcov)
                spv   = float(dat["somatic_p_value"])
                gpv   = float(dat["variant_p_value"])

                if ncov < args.min_normal_depth: filter_list.append("LowDPN")
                if tcov < args.min_tumor_depth: filter_list.append("LowDPT")
                if ss == "GERM" and gpv >= args.pval_cutoff: filter_list.append("GPV")
                if ss == "SOM":
                    if spv >= args.pval_cutoff: filter_list.append("SPV")
                    if nfreq >= args.max_alt_freq_normal: filter_list.append("NAF")
                    if tfreq < args.min_alt_freq_tumor: filter_list.append("TAF")

                if not filter_list: vcf_row.append("PASS")
                else: vcf_row.append(";".join(filter_list))

                info = "SS={0};GPV={1:.2e};SPV={2:.2e}".format(ss, gpv, spv)
                vcf_row.append(info)
                fmt = "GT:RD:AD:FREQ:DP"
                vcf_row.append(fmt)

                alleles = [dat['ref'], dat['var']]
                newngt  = sorted([alleles.index(i) for i in iupac[dat['normal_gt']]])
                newtgt  = sorted([alleles.index(i) for i in iupac[dat['tumor_gt']]])
                ngt     = "{0}/{0}".format(*newngt) if len(newngt) == 1 else "{0}/{1}".format(*newngt)
                tgt     = "{0}/{0}".format(*newtgt) if len(newtgt) == 1 else "{0}/{1}".format(*newtgt)
                nvcf = "{0}:{1}:{2}:{3:.4f}:{4}".format(ngt, nr1, nr2, nfreq, ncov)
                tvcf = "{0}:{1}:{2}:{3:.4f}:{4}".format(tgt, tr1, tr2, tfreq, tcov)
                vcf_row.append(nvcf)
                vcf_row.append(tvcf)
                dic[c][p] = "\t".join(vcf_row)
            except:
                logger.warning('Unable to parse VarScan2 file "{0}" - line {1}: "{2}"... skipping'.format(
                               os.path.basename(args.input_snv), linect, line.rstrip('\n')))
                pass
    return dic

def load_indel(args, dic):
    '''Loads the InDel VarScan format into a dictionary, applies filters,
       and converts to VCF format'''
    head   = []
    linect = 0
    for line in open(args.input_indel, 'rU'):
        linect += 1
        if not head: head = line.rstrip().split('\t')
        else:
            try:
                dat = dict(zip(head, line.rstrip('\n').split('\t')))
                ss  = dat['somatic_status']
                if ss not in lookup: continue
                ss  = lookup[ss]
                c   = dat['chrom']
                p   = int(dat['position'])
                if c not in dic: dic[c] = {}

                vcf_row = [dat['chrom'], dat['position'], '.']

                filter_list = []

                nr1 = int(dat["normal_reads1"])
                nr2 = int(dat["normal_reads2"])
                tr1 = int(dat["tumor_reads1"])
                tr2 = int(dat["tumor_reads2"])
                sspv = float(dat["somatic_p_value"])

                ncov = nr1 + nr2
                tcov = tr1 + tr2
                nfreq = nr2 / float(ncov)
                tfreq = tr2 / float(tcov)
                spv   = float(dat["somatic_p_value"])
                gpv   = float(dat["variant_p_value"])

                # Depth filters
                if ncov < args.min_normal_depth: filter_list.append("LowDPN")
                if tcov < args.min_tumor_depth: filter_list.append("LowDPT")

                # Pvalue filters
                if ss == "GERM" and gpv >= args.pval_cutoff: filter_list.append("GPV")
                if ss == "SOM":
                    if spv >= args.pval_cutoff: filter_list.append("SPV")
                    if nfreq >= args.max_alt_freq_normal: filter_list.append("NAF")
                    if tfreq < args.min_alt_freq_tumor: filter_list.append("TAF")

                info = "SS={0};GPV={1:.2e};SPV={2:.2e};INDEL".format(ss, gpv, spv)
                fmt = "GT:RD:AD:FREQ:DP"

                if dat["var"].startswith('+'):
                    alleles = dat["ref"].split(',') + dat["var"].split(',')
                    var     = dat["var"].replace('+', dat["ref"])
                    normal_gt  = dat['normal_gt']
                    tumor_gt   = dat["tumor_gt"]
                    newngt, newtgt = "", ""

                    vcf_row = vcf_row + [dat["ref"], var, '.']

                    if normal_gt in iupac:
                        newngt = sorted([alleles.index(i) for i in iupac[normal_gt]])
                    else:
                        newngt = sorted([0 if i == '*' else alleles.index(i) for i in normal_gt.split('/')])

                    if tumor_gt in iupac:
                        newtgt = sorted([alleles.index(i) for i in iupac[tumor_gt]])
                    else:
                        newtgt = sorted([0 if i == '*' else alleles.index(i) for i in tumor_gt.split('/')])

                    if not filter_list: vcf_row.append("PASS")
                    else: vcf_row.append(";".join(filter_list))
                    vcf_row = vcf_row + [info, fmt]

                    ngt = "{0}/{0}".format(*newngt) if len(newngt) == 1 else "{0}/{1}".format(*newngt)
                    tgt = "{0}/{0}".format(*newtgt) if len(newtgt) == 1 else "{0}/{1}".format(*newtgt)
                    nvcf = "{0}:{1}:{2}:{3:.4f}:{4}".format(ngt, nr1, nr2, nfreq, ncov)
                    tvcf = "{0}:{1}:{2}:{3:.4f}:{4}".format(tgt, tr1, tr2, tfreq, tcov)
                    vcf_row.append(nvcf)
                    vcf_row.append(tvcf)
                    dic[c][p] = "\t".join(vcf_row)

                elif dat["var"].startswith('-'):
                    alleles = dat["ref"].split(',') + dat["var"].split(',')
                    varallele = dat["var"].replace('-', '')
                    tmp     = dat["ref"]
                    ref     = dat["ref"] + varallele
                    var     = tmp
                    normal_gt  = dat['normal_gt']
                    tumor_gt   = dat["tumor_gt"]
                    newngt, newtgt = "", ""

                    vcf_row = vcf_row + [ref, var, '.']

                    if normal_gt in iupac:
                        newngt = sorted([alleles.index(i) for i in iupac[normal_gt]])
                    else:
                        newngt = sorted([0 if i == '*' else alleles.index(i) for i in normal_gt.split('/')])

                    if tumor_gt in iupac:
                        newtgt = sorted([alleles.index(i) for i in iupac[tumor_gt]])
                    else:
                        newtgt = sorted([0 if i == '*' else alleles.index(i) for i in tumor_gt.split('/')])

                    if not filter_list: vcf_row.append("PASS")
                    else: vcf_row.append(";".join(filter_list))
                    vcf_row = vcf_row + [info, fmt]
                    ngt = "{0}/{0}".format(*newngt) if len(newngt) == 1 else "{0}/{1}".format(*newngt)
                    tgt = "{0}/{0}".format(*newtgt) if len(newtgt) == 1 else "{0}/{1}".format(*newtgt)
                    nvcf = "{0}:{1}:{2}:{3:.4f}:{4}".format(ngt, nr1, nr2, nfreq, ncov)
                    tvcf = "{0}:{1}:{2}:{3:.4f}:{4}".format(tgt, tr1, tr2, tfreq, tcov)
                    vcf_row.append(nvcf)
                    vcf_row.append(tvcf)
                    dic[c][p] = "\t".join(vcf_row)
            except:
                logger.warning('Unable to parse VarScan2 file "{0}" - line {1}: "{2}"... skipping'.format(
                               os.path.basename(args.input_indel), linect, line.rstrip('\n')))
                pass
    return dic

def run(args):
    '''Main wrapper function for filtering VarScan2 files'''
    # Print the filters to the log
    logger.info('<FILTER>LowDPN=Normal DP < {0.min_normal_depth}'.format(args))
    logger.info('<FILTER>LowDPT=Tumor DP < {0.min_tumor_depth}'.format(args))
    logger.info('<FILTER>NAF=Normal FREQ >= {0.max_alt_freq_normal:.2f}'.format(args))
    logger.info('<FILTER>TAF=Tumor FREQ < {0.min_alt_freq_tumor:.2f}'.format(args))
    logger.info('<FILTER>SPV=SPV >= {0.pval_cutoff:.4e}'.format(args))
    logger.info('<FILTER>GPV=GPV >= {0.pval_cutoff:.4e}'.format(args))

    # Chromosome order list
    contigs    = load_contigs(args)
    chrm_order = []
    for i in contigs:
        chrm_order.append(pattern.match(i).groups()[0])

    ## New filter and info lines
    new_header = """##fileformat=VCFv4.1
##fileDate={fdate}
##source=VarScan2
##reference=file://{ref}
{contigs}
##INFO=<ID=SS,Number=1,Type=String,Description="Somatic status call SOMATIC=somatic, LOH=loss of het, GERM=germline">
##INFO=<ID=GPV,Number=1,Type=Float,Description="Variant p-value for germline events">
##INFO=<ID=SPV,Number=1,Type=Float,Description="Somatic p-value for Somatic/LOH events">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="If this is an InDel">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reads supporting reference in sample">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Reads supporting variant in sample">
##FORMAT=<ID=FREQ,Number=1,Type=Float,Description="Variant allele frequency in sample">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Concensus genotype call">
##FORMAT=<ID=DP,Number=4,Type=Integer,Description="REF reads +,REF read -, ALT reads +, ALT reads -">
##FILTER=<ID=LowDPN,Description="Normal DP < {0.min_normal_depth}">
##FILTER=<ID=LowDPT,Description="Tumor DP < {0.min_tumor_depth}">
##FILTER=<ID=NAF,Description="(SS==SOM) && (Normal FREQ >= {0.max_alt_freq_normal:.2f})">
##FILTER=<ID=TAF,Description="(SS==SOM) && (Tumor FREQ < {0.min_alt_freq_tumor:.2f})">
##FILTER=<ID=SPV,Description="(SS==SOM) && (SPV >= {0.pval_cutoff:.4e}">
##FILTER=<ID=GPV,Description="(SS==GERM) && (GPV >= {0.pval_cutoff:.4e}">
""".format(args,fdate=datetime.date.today(),ref=os.path.abspath(args.reference),
           contigs="\n".join(contigs))

    # CHROM header line
    chromline = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}\t{1}\n'.format(
                 args.normal_name, args.tumor_name)

    # Load InDels and SNVs 
    var_dict = {}
    var_dict = load_indel(args, var_dict)
    var_dict = load_snv(args, var_dict)

    # Write the new VCF file 
    with open(args.output_vcf, 'wb') as o:
        # Write header
        o.write(new_header)
        # Write CHROM line
        o.write(chromline)
        # Write results
        for c in chrm_order:
            if c in var_dict:
                for p in sorted(var_dict[c]):
                    o.write(var_dict[c][p] + '\n')
