import sys
import celltools

from optparse import OptionParser

###############################################################################
USAGE = """
tabulate-supporting-reads.py.py 
                         --bam <10X bam file to use> 
                         --candidate <table of candidate sites, Yifan csv format> 
                         --out <file for output>

Assumes/only deails with biallelic SNVs
Has some built in  filters for quality.
"""
parser = OptionParser(USAGE)
parser.add_option('--out',dest='outFileName', help = 'name of output table file')
parser.add_option('--bam',dest='bamFileName', help = 'name of 10X  single cell BAM file')
parser.add_option('--candidate',dest='candidateFileName', help = 'candidateFileName not given')


(options, args) = parser.parse_args()


if options.outFileName is None:
    parser.error('outFileName not given')
if options.bamFileName is None:
    parser.error('bamFileName not given')
if options.candidateFileName is None:
    parser.error('candidateFileName not given')
###############################################################################


outFile = open(options.outFileName,'w')
inFile = open(options.candidateFileName,'r')

header = ['chrom','pos','ref','alt','totReads','readFractionAlt','refReads','altReads','otherReads',
              'totCells','cellFractionAlt','refCellCount','altCellCount','otherCellCount','refCells','altCells','otherCells']
              
header = '\t'.join(header) + '\n'
outFile.write(header)              

for line in inFile:
    line = line.rstrip()
    line = line.split(',')
    if line[0] == 'chrm':
        continue
    c = line[0]
    pos = int(line[1])
    ref = line[2]
    alt = line[3]
    
    print line[0:7]
    
    counts = celltools.count_reads(options.bamFileName,c,pos,ref,alt)
    
    nl = [c,str(pos),ref,alt]
    
    refReads = counts[3]
    altReads = counts[4]
    totReads = refReads + altReads
    fReads = float(altReads)/refReads
    nl.append(totReads)
    nl.append(fReads)
    nl.append(refReads)
    nl.append(altReads)
    nl.append(counts[5])
    
    
    totCellsRefAlt = len(counts[0]) + len(counts[1])
    nl.append(totCellsRefAlt)
    fAlt = float(len(counts[1]))/totCellsRefAlt
    nl.append(fAlt)
    nl.append(len(counts[0]))
    nl.append(len(counts[1]))
    nl.append(len(counts[2]))
    
    k = counts[0].keys()
    k.sort()
    if len(k) == 0:
        k = ['.']

    k = ':'.join(k)
    nl.append(k)

    k = counts[1].keys()
    k.sort()
    if len(k) == 0:
        k = ['.']
    k = ':'.join(k)
    nl.append(k)

    k = counts[2].keys()
    k.sort()
    if len(k) == 0:
        k = ['.']
    k = ':'.join(k)
    nl.append(k)
    
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)
    
inFile.close()
outFile.close()

###############################################################################
