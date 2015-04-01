#!/usr/bin/env python
# Kyle Hernandez
# Center for Research Informatics
# University of Chicago
#
# This is adapted from P. Cingolani's script to evenly split callable regions

import sys
import os

# Max distance to merge
maxDistance = 1000

outLines = []
sizes = {}

def main(bed, cfg, msplits, pfx):
    '''
    ------------------------------------------------------------------
    Merges and splits callable regions from a bed file

    Usage: SplitBed.py <in.bed> <out.cfg> <max splits> <output prefix>
    ------------------------------------------------------------------
    ''' 
    totalBases = processBed(bed)
    avgBases   = totalBases / msplits + 1
    avgLines   = len(outLines) / msplits + 1
    ffiles = writeBedFiles( pfx, avgBases, msplits )
    opath  = os.path.abspath(pfx) 
    with open(cfg, 'wb') as o:
        o.write('totalBases = {0}\n'.format(totalBases))
        o.write('avgBases = {0}\n'.format(avgBases))
        o.write('avgLines = {0}\n'.format(avgLines))
        o.write('totalSplits = {0}\n'.format(len(ffiles)))
        o.write('outputPrefix = {0}\n'.format(opath))
 
def processBed(bed):
    '''
    Reads the BED file, and calculates the size and distances.
    '''
    # Initialize
    tBases    = 0
    pscaff, pstart, pstop = '', -1, -1
    callStart = -1
 
    # Read bed file
    for line in open(bed, 'rU'):
        cols   = line.rstrip().split("\t")
        cscaff, cstart, cstop = cols[0], int(cols[1]), int(cols[2])
        
        # Calculate size and distance
        size = cstop - cstart
        dist = cstart - pstop
 
        # Need to split?
        if cscaff != pscaff or dist > maxDistance:
            tBases   += addBed(pscaff, callStart, pstop)
            callStart = cstart 
        
        pscaff, pstart, pstop = cscaff, cstart, cstop
 
    if callStart != cstart:
        tBases += addBed(pscaff, callStart, pstop)
    return tBases
 
def addBed(chrm, s, e):
    if s > 0:
        size    = e - s
        outLine = chrm + '\t' + str(s) + '\t' + str(e) 
        outLines.append( outLine )
        sizes[outLine] = size
        return size
    return 0

def writeBedFiles( pfx, avgBases, msplits ):
    bases      = 0
    outFileIdx = 0
    outFiles   = [os.path.abspath('{0}.{1}.callable.bed'.format(pfx, i)) for i in range(msplits)]
    finalFiles = []
    hold       = []

    for outLine in outLines:
        if not hold: 
            hold.append(outLine)
            bases += sizes[ outLine ]

        elif bases > avgBases and outFileIdx < msplits-1 :
            finalFiles.append(outFiles[outFileIdx])
            with open(outFiles[outFileIdx], 'wb') as o: 
                o.write("\n".join(hold) + "\n")

            hold        = []
            outFileIdx += 1
            bases       = 0
            hold.append(outLine)
            bases += sizes[ outLine ]
        else:
            hold.append(outLine)
            bases += sizes[ outLine ]

    if hold:
        finalFiles.append(outFiles[outFileIdx])
        with open(outFiles[outFileIdx], 'wb') as o: 
            o.write("\n".join(hold) + "\n")

    return finalFiles

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print main.__doc__
        sys.exit(1)

    inbed     = sys.argv[1]
    ocfg      = sys.argv[2]
    maxSplits = int(sys.argv[3])
    odir      = sys.argv[4]

    main(inbed, ocfg, maxSplits, odir) 
