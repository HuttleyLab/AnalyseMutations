from cogent import LoadTable, DNA

__author__ = "Yicheng Zhu, Gavin Huttley"
__copyright__ = "Copyright 2015, Gavin Huttley"
__credits__ = ["Yicheng Zhu", "Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

def reverse_complement_record(gene_strand, snp_strand):
    """returns True if the Ensembl 71 SNP record needs to be reverse complemented to put on transcribed strand"""
    # rc if the strands are different
    if gene_strand not in [-1,1] or snp_strand not in [-1, 1]:
        raise ValueError("strand must be -1 or 1")
    
    return snp_strand != gene_strand

def get_rc_record(alleles, ancestor, allele_freqs, flank_5, flank_3):
    """reverse complements the alleles, ancestror, flanking seqs, and allele freqs
    """
    complement = DNA.complement
    alleles_rc = set([complement(b) for b in alleles])
    ancestor_rc = complement(ancestor)
    
    allele_freqs_rc = {}
    for allele, freq in allele_freqs.items():
        allele_freqs_rc[complement(allele)] = freq
    
    flank_5_rc = str(DNA.makeSequence(flank_5).rc())
    flank_3_rc = str(DNA.makeSequence(flank_3).rc())
    return alleles_rc, ancestor_rc, allele_freqs_rc, flank_3_rc, flank_5_rc

