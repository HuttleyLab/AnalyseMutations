"""sampling SNP data from Ensembl"""
import os, sys, cPickle, gzip
import collections
import csv
import re

from optparse import make_option

from cogent import DNA, Sequence, LoadTable
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.option_parsing import parse_command_line_parameters
from cogent.db.ensembl import Species, HostAccount, Genome
from cogent.db.ensembl.region import Variation

from scitrack import CachingLogger

from mutation_motif.util import abspath, create_path

__author__ = "Yicheng Zhu, Gavin Huttley"
__copyright__ = "Copyright 2015, Gavin Huttley"
__credits__ = ["Yicheng Zhu", "Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

LOGGER = CachingLogger(create_dir=True)
_nucs = set(DNA)

def get_snp_dump_data(record, allele_freqs, alleles, ancestral, flank5, flank3, genic=False):
    """returns all SNP data for dumping."""
    loc_template =\
         "%(Species)s:%(CoordType)s:%(Chrom)s:%(Start)s-%(End)s:%(Strand)s"
    snp_loc = loc_template % record.Location._asdict()
    
    pep_alleles = gene_loc = gene_id = 'None'
    if genic:
        gene_loc = loc_template % record.GeneLocation._asdict()
        gene_id = record.GeneId
        
    output = [record.Symbol, snp_loc, record.Strand,
            record.Effect, allele_freqs, alleles, ancestral,
            flank5, flank3, pep_alleles, gene_id, gene_loc]
    
    return output

def get_location(loc_string):
    """returns genomic location"""
    # Homo sapiens:chromosome:19:107460-111696:1
    header = ["Species", "CoordType", "Chrom", "Start", "End", "Strand"]
    mk_loc = collections.namedtuple("Location", header)
    loc = loc_string.split(':')
    loc = loc[:-2] + map(int, loc[-2].split('-')) + loc[-1:]
    loc = mk_loc(*loc)
    return loc

def read_snp_dump(path, genic=False, limit=None):
    """returns saved SNP record"""
    infile = gzip.open(path)
    record_reader = csv.reader(infile, delimiter='\t')
    header = ["Symbol", "Location", "Strand", "Effect", "AlleleFreqs", "Alleles",
            "Ancestral", "GeneId", "GeneLocation"]
    mk_record = collections.namedtuple("Record", header)
    n = 0
    for record in record_reader:
        n += 1
        if not len(record):
            continue
        
        record[1] = get_location(record[1])
        if genic:
            record[-1] = get_location(record[-1])
            
        record = mk_record(*record)
        yield record
        if limit and n >= limit:
            break
    
    infile.close()

def get_flanks(chrom, location, flank_size):
    """returns sequence of flank_size * 2 + 1"""
    start = location - flank_size
    end = location + flank_size + 1
    segment = chrom[start: end]
    assert len(segment) == (2 * flank_size + 1)
    return segment

assert get_flanks('AAA--X--TTT', 5, 2) == '--X--'

def main(script_info):
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    chroms = {}
    flank_size = 300
    genic = 'missense' in opts.snp_path or 'syn' in opts.snp_path
    record_reader = read_snp_dump(opts.snp_path, genic=genic, limit=opts.limit)
    
    LOGGER.input_file(opts.snp_path)
    
    ancestral_mismatch = 0
    ref_mismatch = 0
    written = 0
    
    out_dir = abspath(opts.output)
    out_path = os.path.join(out_dir, os.path.basename(opts.snp_path))
    runlog_path = "%s.log" % out_path.split('.')[0]
    LOGGER.log_file_path = runlog_path
    
    if not opts.dry_run:
        create_path(out_dir)
        outfile = gzip.open(out_path, 'w')
    
    for record in record_reader:
        if record.Ancestral not in record.Alleles:
            ancestral_mismatch += 1
            continue
        
        snp_strand = record.Strand
        if genic:
            gene_strand = record.GeneLocation.Strand
        else:
            gene_strand = snp_strand
        
        if record.Location.Chrom not in chroms:
            chrom_path = os.path.join(opts.chroms_path,
                            "Chr%s.fa.gz" % record.Location.Chrom)
            
            LOGGER.input_file(chrom_path)
            infile = gzip.open(chrom_path)
            seq = [(l, s) for l, s in MinimalFastaParser(infile)][0][1]
            infile.close()
            chroms[record.Location.Chrom] = seq
        
        chrom = chroms[record.Location.Chrom]
        location = record.Location.Start
        segment = get_flanks(chrom, location, flank_size)
        base = str(segment[flank_size])
        alleles = set(record.Alleles.split('/'))
        ancestral = record.Ancestral
        allele_freqs = eval(record.AlleleFreqs)
        if len(alleles) != 2 or allele_freqs is None:
            continue
        
        # if the Ancestral base does not match an allele, count and continue
        if ancestral not in alleles:
            ancestral_mismatch += 1
            continue
        
        # if genic & the ref-base at location does not match an allele,
        # reverse complement if it's on negative strand
        if genic and base not in alleles and record.GeneLocation.Strand in ('-1', -1):
            alleles = set(DNA.complement(a) for a in alleles)
            allele_freqs = [(DNA.complement(b), f) for b, f in allele_freqs]
        
        if base not in alleles:
            ref_mismatch += 1
            continue
        
        # replace the ref-base with Ancestral (simpler than checking then
        # changing etc..)
        segment = list(segment)
        segment[flank_size] = ancestral
        segment = ''.join(segment)
        
        # slice segment
        flank5 = segment[:flank_size]
        flank3 = segment[flank_size+1:]
        # dump full record
        output = get_snp_dump_data(record, allele_freqs, alleles, ancestral,
                                   flank5, flank3, genic=genic)
        output = map(str, output)
        output = '\t'.join(output) + '\n'
        if not opts.dry_run:
            outfile.write(output)
            written += 1
        
        if written % 1000 == 0:
            sys.stdout.write("Written=%s\r" % (' ' * 20))
            sys.stdout.write("Written=%d\r" % written)
            sys.stdout.flush()
    
    LOGGER.output_file(out_path)
    
    table = LoadTable(header=['category', 'counts'],
         rows=[['ancestral_mismatch', ancestral_mismatch],
               ['ref_mismatch', ref_mismatch],
           ['written', written]])
    
    print table
    
    LOGGER.write("\n" + str(table), label="summary statistics")


script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = "associate the flanks with SNP records."

script_info['required_options'] = [
     make_option('--snp_path', help='SNP dump path.'),
     make_option('--chroms_path', help='Chrom path.'),
     make_option('-o','--output', help='Path to write data.'),
    # see http://asia.ensembl.org/info/docs/variation/predicted_data.html
    ]

script_info['optional_options'] = [
    make_option('-l','--limit', default=None, type=int,
        help='Number of results to return.'),
    make_option('-D','--dry_run', action='store_true', default=False,
        help='Do a dry run of the analysis without writing output.'),
    ]

script_info['version'] = '0.1'
script_info['authors'] = 'Gavin Huttley'


if __name__ == "__main__":
    main(script_info)
    