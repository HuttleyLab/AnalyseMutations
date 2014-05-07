import sys, os

from cogent.util.unit_test import TestCase, main

from strand import reverse_complement_record, get_rc_record
from snptables_to_aln import filtered_records

class TestStrand(TestCase):
    def test_reverse_record(self):
        """correctly identify when a record needs to be reversed"""
        # gene is on plus strand and snp is on plus strand
        self.assertFalse(reverse_complement_record(gene_strand=1, snp_strand=1))
        # gene is on plus strand and snp is on minus strand
        self.assertTrue(reverse_complement_record(gene_strand=1, snp_strand=-1))
        # gene is on minus strand and snp is on plus strand
        self.assertTrue(reverse_complement_record(gene_strand=-1, snp_strand=1))
        # gene is on minus strand and snp is on minus strand
        self.assertFalse(reverse_complement_record(gene_strand=-1, snp_strand=-1))
        # invalid data
        self.assertRaises(ValueError, reverse_complement_record, 0, 0)
    
    def test_rc_record(self):
        """reverse complement a record"""
        alleles = set('TC')
        ancestor = 'T'
        af = {'C': 0.9885, 'T': 0.0115}
        f_seq_5 = 'ATAAGATCTTATGTTCTT'
        f_seq_3 = 'GACATTGCAGGAAAATCT'
        alleles_rc, ancestor_rc, allele_freqs_rc, flank_5_rc, flank_3_rc = \
                    get_rc_record(alleles, ancestor, af, f_seq_5, f_seq_3)
        self.assertEqual(alleles_rc, set('AG'))
        self.assertEqual(ancestor_rc, 'A')
        expect_freqs = {'A': 0.0115, 'G': 0.9885}
        self.assertEqual(allele_freqs_rc, expect_freqs)
    
    def test_gen_stranded_fasta(self):
        """correctly generate a stranded fasta formatted sequence file"""
        data = [line.strip() for line in open('test_data/strand.txt')]
        seen = set()
        chroms = set()
        expect_records = [
        '>rs141018741\nATACTGTTGCCTCAATATCAAGGAACTGGAAGCATGTAAAA\n',
        '>rs114239297\nTGTGATGTTGGAGAACTACAAGAACCTGATCTCCCTGGGTG\n',
        '>rs4983150\nTATCATGACTGGGCATCTCTATGCCAAACTGGTCATACTGT\n',
        '>rs4572095\nCTGGATTCTGGGAGTCTGTCAAGATTCCAGGACTGCAGATG\n'
        ]
        
        for record in filtered_records(data[1:], ('A', 'G'), seen, chroms, verbose=False):
            self.assertTrue(record in expect_records)
    
    

if __name__ == "__main__":
    main()
