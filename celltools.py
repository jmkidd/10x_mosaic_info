#celltools.py
import sys
import pysam

#######################################################################
def count_reads(bamFileName,chrom,pos,ref,alt):
    #  counts will be list of cell barcodes supporting ref, alt, and other alleles
    # impose some quality cutoffs
    # cells ref,alt,none, reads ref,alt,none
    counts = [{},{},{},0,0,0]   
    minMapQual = 10
   
   
    samfile = pysam.AlignmentFile(bamFileName, 'rb')
    for read in samfile.fetch(chrom,pos-1,pos):
        # some simple filters, not being very strict at all...
        if read.is_duplicate is True:
            continue
        if read.is_qcfail is True:
            continue    
        if read.mapping_quality < minMapQual:
            continue    
        
        aligned_pairs = read.get_aligned_pairs(with_seq=False)         
        queryAlignedBase = '?' # incase we don't have an aligned base, probably not needed
        numAligned = 0
        for i in aligned_pairs:
            if i[1] == pos-1:
                if i[0] == None: # deletion
                    continue
                queryAlignedBase = read.query_sequence[i[0]]
                numAligned += 1
        if numAligned != 1:  #might see multiple aligned in case of indel
            continue
        if queryAlignedBase == '?' or queryAlignedBase == 'N':  # no alignment...
            continue
        if queryAlignedBase == ref:
            allele_i = 0
        elif queryAlignedBase == alt:
            allele_i = 1
        else:
            allele_i = 2
            
        counts[allele_i + 3] += 1    

        cellTag = 'NoTag'
        if read.has_tag('CB'):
            cellTag = read.get_tag('CB')
        
        if cellTag in counts[allele_i]:
            counts[allele_i][cellTag] += 1
        else:
            counts[allele_i][cellTag] = 1        
    samfile.close()
    return counts
#######################################################################

