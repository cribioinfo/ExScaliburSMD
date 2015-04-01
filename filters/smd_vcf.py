'''Functions for filtering Mutect output VCF files'''
import logging

class VcfReader(object):
    '''Object for parsing a VCF file'''
    def __init__(self, smd, tname, nname, infile):
        self.logger = logging.getLogger('SomaticFilter.VcfReader')
        self.smd    = smd 
        self.tname  = tname 
        self.nname  = nname 
        self.infile = infile 
        self.fh     = None
        self.header = VcfHeader()
        self.tidx   = None
        self.nidx   = None

    def Open(self):
        '''Simply opens the file'''
        self.logger.info('Opening VCF file {0.infile}'.format(self))
        self.fh = open(self.infile, 'rU') 

    def Close(self):
        '''Simply closes the file'''
        if self.fh.closed: pass
        else:
            self.logger.info('Closing VCF file {0.infile}'.format(self))
            self.fh.close()

    def get_header(self):
        '''Loads the header information into the VcfHeader instance'''
        self.logger.info('Loading VCF header...')
        self.header.load_header(self.fh)
        if self.smd == 'mutect':
            assert self.tname in self.header.header, \
                self.logger.error('Tumor sample {0.tname} does not match VCF header'.format(self))
            assert self.nname in self.header.header, \
                self.logger.error('Normal sample {0.nname} does not match VCF header'.format(self))
            self.tidx = self.header.header.index(self.tname)
            self.nidx = self.header.header.index(self.nname)
        elif self.smd == 'virmid':
            pass 
        else:
            self.tidx = self.header.header.index('TUMOR')
            self.nidx = self.header.header.index('NORMAL')

    def write_new_header(self, o, filters=None, info=None, contigs=None, refpath=None, formats=None):
        '''Writes the new header to the file given the existing headers and 
           the new filters/info and sample names.
        '''
        self.logger.info('Writing new header...')

        # If available, write out the vcffmt row, vcfdate row, source row and reference path row
        if self.header.vcffmt: o.write(self.header.vcffmt + '\n')
        if self.header.vcfdate: o.write(self.header.vcfdate + '\n')
        if self.header.src: o.write(self.header.src + '\n')
        else: o.write('##source={0.smd}\n'.format(self))
        if self.header.ref: o.write(self.header.ref + '\n')
        if refpath: o.write(refpath + '\n')

        # Write contigs
        if self.header.contig: o.write('\n'.join(self.header.contig) + '\n')
        if contigs: o.write('\n'.join(contigs) + '\n')

        # Write other information 
        if self.header.other: o.write('\n'.join(self.header.other) + '\n')

        # Write old and new info
        if self.header.info: o.write('\n'.join(self.header.info) + '\n')
        if info: o.write('\n'.join(info) + '\n')

        # Write old and new filters
        if self.header.flt: o.write('\n'.join(self.header.flt) + '\n')
        if filters: o.write('\n'.join(filters) + '\n')

        # Write old and new formats 
        if self.header.fmt: o.write('\n'.join(self.header.fmt) + '\n')
        if formats: o.write('\n'.join(formats) + '\n')

        # Write new CHROM header line. We always keep NORMAL column first, then TUMOR column.
        # Strelka does not output these columns, so we skip it.
        if self.smd == 'mutect':
            new_header = "\t".join(self.header.header[:9]) + '\t' + \
                         self.header.header[self.nidx] + '\t' + \
                         self.header.header[self.tidx]
            o.write(new_header + '\n')
        else:
            new_header = "\t".join(self.header.header[:9]) + '\t' + \
                         self.nname + '\t' + self.tname
            o.write(new_header + '\n')

    def apply_filters(self, o, args, RecordClass, **kwargs):
        '''Given the instance of a VcfRecord object, apply filters and write new records''' 
        self.logger.info('Applying {0.smd} filters...'.format(self))
        for line in self.fh:
            record = RecordClass(line.rstrip().split('\t'), normal_idx=self.nidx, tumor_idx=self.tidx)
            record.apply_filters(args, **kwargs)
            record.write_record(o)

    def __del__(self):
        '''Another safeguard for making sure the file stream is closed.'''
        self.Close()
        del self.fh

class VcfHeader(object):
    '''Object for parsing the VCF metadata'''
    def __init__(self):
        self.vcffmt  = None
        self.vcfdate = None
        self.src     = None
        self.flt     = [] 
        self.fmt     = [] 
        self.cmd     = None 
        self.info    = [] 
        self.contig  = [] 
        self.ref     = None
        self.other   = []
        self.header  = None

    def load_header(self, fh):
        for line in fh:
            if line.startswith('##'):
                if line.startswith('##fileformat'): self.vcffmt = line.rstrip()
                elif line.startswith('##fileDate'): self.vcfdate = line.rstrip()
                elif line.startswith('##source'): self.src = line.rstrip()
                elif line.startswith('##FILTER'): self.flt.append(line.rstrip())
                elif line.startswith('##FORMAT'): self.fmt.append(line.rstrip())
                elif line.startswith('##INFO'): self.info.append(line.rstrip())
                elif line.startswith('##contig'): self.contig.append(line.rstrip())
                elif line.startswith('##reference'): self.ref = line.rstrip()
                else: self.other.append(line.rstrip())
            elif line.startswith('#CHROM'):
                self.header = line.rstrip().split('\t')
                break 

class VcfRecord(object):
    '''Abstract VCF row'''
    def __init__(self, row, normal_idx=None, tumor_idx=None):
        self.chrom  = row[0]
        self.pos    = int(row[1])
        self.dbid   = row[2]
        self.ref    = row[3]
        self.alt    = row[4]
        self.qual   = row[5]
        self.filt   = row[6]
        self.info   = row[7]
        self.fmt    = None if normal_idx is None and tumor_idx is None else row[8] 
        self.normal = None if normal_idx is None else row[normal_idx]
        self.tumor  = None if tumor_idx is None else row[tumor_idx]

    def apply_filters(self, args):
        '''Applies various filters to the VCF row with side-effects on the info
           and filt attributes. This should be implemented by the child class.
        '''
        pass
       
    def write_record(self, o):
        '''Writes the VcfRecord to the output stream o'''
        if self.normal is None and self.tumor is None:
            o.write("{0.chrom}\t{0.pos}\t{0.dbid}\t{0.ref}\t{0.alt}\t{0.qual}\t{0.filt}\t".format(self) + \
                "{0.info}\n".format(self))
        else:
            o.write("{0.chrom}\t{0.pos}\t{0.dbid}\t{0.ref}\t{0.alt}\t{0.qual}\t{0.filt}\t".format(self) + \
                "{0.info}\t{0.fmt}\t{0.normal}\t{0.tumor}\n".format(self))

# Function for building the contig lines
def load_contigs(args):
    '''Loads the contig name and size information from the sequence dictionary
       into a list'''
    lst   = []
    infil = args.reference + '.dict'
    try:
        fh = open(infil, 'rU')
    except:
        try:
            infil = args.reference.replace('.fasta', '.dict') if args.reference.endswith('.fasta') else \
                    args.reference.replace('.fa', '.dict')
            fh = open(infil, 'rU')
        except:
            raise IOError("No such file or directory: '{0}'".format(infil)) 
    finally:
        fh.close()

    for line in open(infil, 'rU'):
        if line.startswith('@SQ'):
            cols = line.rstrip().split("\t")
            contig = ""
            length = ""
            for i in cols:
                if i.startswith('SN'): contig = i.split(':')[1]
                elif i.startswith('LN'): length = i.split(':')[1]
            lst.append('##contig=<ID={0},length={1}>'.format(contig, length))
    return lst
