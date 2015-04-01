'''Read-level information for Samples. Contains metadata for each row in the metadata table.'''

from datetime import datetime
import logging
import os
import re

#from .preprocessing.fastqc import FastQC
#from .preprocessing.seqprep import SeqPrep
#from .preprocessing.cutadapt import Cutadapt 
#
#from .alignment.bwa import BwaAln
#from .alignment.bwa import BwaMem
#from .alignment.novoalign import Novoalign

class IlluminaSeq(object):
    '''Class representing metadata for a SE or PE Illumina Fastq file(s).
      
       -- row - the current row in the metadata table
       -- number - the line number
       -- sample - sample name
       -- readgroup - readgroup id
       -- phenotype - either tumor or normal
       -- path - the path to the parent directory containing the fastq files
       -- leftseq - file name for read 1
       -- rightseq - file name for read 2
       -- flowcell - the flowcell id
       -- lane - the lane number
       -- library - the library id (set to sample if not provided)
       -- center - the sequencing center shortname
       -- unit - the platform unit
       -- date - sequencing run date (MM/DD/YYYY format)
       -- barcode - the barcode sequence
       -- fragment - the mean fragment size (250 default)
       -- sdfragment - the sd fragment size (50 default)
       -- offset - Phred-score offset (33 by default; 64 possible)
       -- paired - Boolean whether it's paired or not
       -- readlength - the length of the reads
    '''
    def __init__(self, row, number):
        self.logger     = logging.getLogger('SMD.IlluminaSeq')
        self.row        = row
        self.number     = number
        self.sample     = row['sample'].replace(' ','_')
        self.readgroup  = row['readgroup'].replace(' ', '_') 
        self.phenotype  = row['phenotype'].upper()
        self.path       = row['path']
        self.leftseq    = row['fastq_r1']
        self.rightseq   = row.get('fastq_r2')
        self.flowcell   = row.get('flowcell')
        self.lane       = int(row['lane'])
        self.library    = row.get('library', self.sample) 
        self.center     = row.get('center')
        self.unit       = row.get('unit')
        self.date       = row.get('date')
        self.barcode    = row.get('barcode')
        self.fragment   = int(row.get('mean_fragment_size', 250))
        self.sdfragment = int(row.get('sd_fragment_size', 50))
        self.offset     = self._set_offset() 
        self.readgroup_cfg = self.readgroup + '.cfg' 
        self.paired     = None
        self.readlength = None
        self._get_chemistry()

        # Just another check to be safe
        self.__validate()
     
        # Set the rg string
        self.rgstring   = self._set_rgstring() 

        # Containers for read-level analyses
        self.fastqc     = None
        self.seqprep    = None
        self.cutadapt   = None
        self.bwa_aln    = None
        self.bwa_mem    = None
        self.novoalign  = None

    def __str__(self):
        return "\t".join([self.sample, self.phenotype, self.readgroup, 
                          str(self.offset), str(self.paired), str(self.readlength)])

    def _set_rgstring(self):
        '''Sets the RG string based on the metadata provided into a format that can be used
        by BWA.
        '''
        rg = ['@RG', 'ID:'+self.readgroup,
              'PL:illumina', 'LB:'+self.library,
              'SM:'+self.sample+'-'+self.phenotype]
        if self.center: rg.append('CN:'+self.center)
        if self.unit: rg.append('PU:'+self.unit)
        if self.date: rg.append('DT:'+self.date)
        return "\"" + '\\t'.join(rg) + "\""

    def __validate(self):
        '''Runs checks to make sure everything is ok'''
        assert self.sample, 'Must provide sample ID. Row number: %r' % self.number
        assert self.readgroup, 'Must provide readgroup ID. Row number: %r' % self.number
        assert os.path.isdir(self.path), \
          "Path on row %r is not an existing directory.  '%r'" % (self.number, self.path)

        # Check if fastq files exist
        try:
            fp = open(os.path.join(self.path, self.leftseq))
        except IOError as e:
            raise IOError('ERROR!! Unable to access fastq file: {0}'.format(os.path.join(self.path, self.leftseq)))
        else:
            fp.close()

        if self.paired: 
            try:
                fp = open(os.path.join(self.path, self.rightseq))
            except IOError as e:
                raise IOError('ERROR!! Unable to access fastq file: {0}'.format(
                              os.path.join(self.path, self.rightseq)))
            else:
                fp.close()

        # If the library column was in the table, but was empty row.get() would still get it
        if not self.library: self.library = self.sample
        else: self.library = self.library.replace(' ', '_') 

        # If date, let's format it
        if self.date:
            self.date = datetime.strptime(self.date, '%m/%d/%Y').strftime('%m/%d/%Y')

        # Validate other data
        if self.center: self.center = self.center.replace(' ', '_') 
        if self.unit: self.unit = self.unit.replace(' ', '_') 
        if self.flowcell: self.flowcell = self.flowcell.replace(' ', '_') 

    def _get_chemistry(self):
        '''Sets the chemistry information by using regex. For example, an entry of 2x100 or
        100x2 would be parsed as a paired-end read with 100bp length.

        :raises: KeyError
        '''
        if re.match(r'^(\d+)x2$', self.row['chemistry']):
            self.paired = True
            self.readlength = int(re.match(r'^(\d+)x2$', self.row['chemistry']).group(0))
            assert self.rightseq.strip(), 'You must provide Read 2 file name if paired!! Line %r' % (self.number)
        elif re.match(r'^(\d+)x1$', self.row['chemistry']):
            self.paired = False
            self.readlength = int(re.match(r'^(\d+)x1$', self.row['chemistry']).group(0))
        elif re.match(r'^2x(\d+)$', self.row['chemistry']):
            self.paired = True
            self.readlength = int(re.match(r'^2x(\d+)$', self.row['chemistry']).group(1))
            assert self.rightseq.strip(), 'You must provide Read 2 file name if paired!! Line %r' % (self.number)
        elif re.match(r'^1x(\d+)$', self.row['chemistry']):
            self.paired = False
            self.readlength = int(re.match(r'^1x(\d+)$', self.row['chemistry']).group(1))
        else:
            self.logger.error("Unrecognized chemistry column {0} for fastq_r1 {1}".format(
                self.row['chemistry'], self.row['fastq_r1']))
            raise KeyError("Unrecognized chemistry '%r' for fastq_r1 '%r'" % \
                  (self.row['chemistry'], self.row['fastq_r1']))

    def _set_offset(self):
        '''Sets the Phred offset as an integer. If not provided, a warning is printed adn it is
        set to 33

        :raises: KeyError
        '''
        if not self.row.get('offset'):
            self.logger.warn('Offset not defined for {0}. Assuming Phred33'.format(self.row['fastq_r1']))
            return 33
        else:
            try:
                return int(self.row['offset'])
            except:
                if re.match(r'^Phred33', self.row['offset']):
                    return 33
                elif re.match(r'^Phred64', self.row['offset']):
                    return 64
                else:
                    raise KeyError("Offset {0} not recognized for fastq_r1 {1}".format(
                        self.row['offset'], self.row['fastq_r1']))

    def __getitem__(self, key):
        if key == 'bwa_aln': return self.bwa_aln
        elif key == 'bwa_mem': return self.bwa_mem
        elif key == 'novoalign': return self.novoalign
        else:
            raise KeyError('Unknown aligner {0}'.format(key))


