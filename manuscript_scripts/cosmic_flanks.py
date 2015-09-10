import os
import gzip
import sys
import click

from cogent import Sequence, DNA
from cogent.parse.fasta import MinimalFastaParser
from scitrack import CachingLogger

from read_cosmic import reader

__author__ = "Bob Buckley, Gavin Huttley"
__copyright__ = "Copyright 2015, Bob Buckley"
__credits__ = ["Bob Buckley", "Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

effects = {'Substitution - Missense': 'missense_variant',
           'Substitution - coding silent': 'synonymous_variant',
           'Substitution - Nonsense': 'nonsense_variant'}

LOGGER = CachingLogger(create_dir=False)

class CosmicRef():
    """Get a reference for work with a COSMIC SNPs file (format specific to this project)
    """
    # COSMIC uses numbers for chromosomes ... these tuple maps COSMIC number to names
    chrom_names = tuple(map(str, range(1,23))) + ('X', 'Y', 'MT')
    nums_names = dict(zip(tuple(str(i+1) for i in range(len(chrom_names))), chrom_names))
    # transstrand = { '+':'1', '-':'-1' }
    
    def __init__(self, chrdir):
        """Class to "encapsulate" most of this application ....
        
        Read the chromosome files from the chromosome directory, store them for easy/faster
        access to flanking regions.
        """
        # we know what files are expected ... check they are all there
        if not os.path.isdir(chrdir):
            print "Didn't find directory", chrdir
            sys.exit(1)
            
        chrom_files = tuple(os.path.join(chrdir, 'Chr%02d_%s.fa.gz'%(i+1, n)) 
                            for i, n in enumerate(self.chrom_names) if n)
        missing = [fn for fn in chrom_files if not os.path.exists(fn)]
        
        if missing:
            print "Missing reference files:", missing
            print "Did you run dump_genome.py?"
            sys.exit(1)

        self.chrom_files = dict(zip(self.chrom_names, chrom_files))
        self._chrom_data = {}
        
        return
    
    def chrom_data(self, chrom_name):
        chrom_name = self.nums_names[chrom_name]
        if chrom_name not in self._chrom_data:
            fn = self.chrom_files[chrom_name]
            if os.path.islink(fn):
                fn = os.readlink(fn)
            
            print "\r%s" % (" " * 80),
            print ("\rloading %s" % fn).ljust(60),
            sys.stdout.flush()
            with gzip.open(fn) as infile:
                data = [r for r in MinimalFastaParser(infile)]
                assert len(data) == 1, len(data)
                name, seq = data[0]
                if not data:
                    print fn, data
                chrom = Sequence(moltype=DNA, name=name, seq=seq)
            
            LOGGER.input_file(fn, label='chromosome_seq')
            self._chrom_data[chrom_name] = chrom
        
        return self._chrom_data[chrom_name]
    
    def get_flanks(self, chrno, loc, flanksize):
        """return the flanking regions for a location in a chromosome"""
        rawflank5 = self.chrom_data(chrno)[loc-flanksize: loc]
        rawflank3 = self.chrom_data(chrno)[loc+1: loc+1+flanksize]
        ref = self.chrom_data(chrno)[loc]
        return rawflank5, ref, rawflank3
    
    def get_snp_dump_data(self, snp, flanksize=300):
        snpchr, locx = snp.Location.split(':', 1) 
        locstart, locend = locx.split('-', 1)
        if locstart!=locend:
            return None
        
        loc = int(locstart) - 1   # convert from 1 based location to zero-based
        try:
            chrom_seq = self.chrom_data(snpchr)
        except KeyError:
            print snp
            raise
        
        if flanksize > loc or loc + 1 + flanksize >= len(chrom_seq):
           # too close to terminii of chromosome
           return None
        
        f5, ref, f3 = self.get_flanks(snpchr, loc, flanksize=flanksize)
        
        # note: Stand in the COSMIC data is '+' or '-'
        strand = '1' if snp.Strand=='+' else '-1'
        
        # the alleles
        orig, dest = snp.Change[-3], snp.Change[-1]
        if strand == '-1':
            # since the snptables to aln script corrects for strand, we need
            # to leave things on the original strand and just reverse the
            # alleles to make them consistent with the genome seq strand
            orig = DNA.complement(orig)
            dest = DNA.complement(dest)
        
        if ref != orig:
            return None
        
        location = 'Homo sapiens:chromosome:%s:%s-%s:%s' % (snpchr, loc, loc+1, strand)
        
        # some data validity checks!
        alleles = str(set((snp.Change[-3], snp.Change[-1])))
        pep_alleles = gene_loc = gene_id = allele_freqs = 'None'
        effect = effects.get(snp.Type, None)
        if effect is None:
            raise ValueError(str(snp))

        return (snp.Symbol, location, strand, effect, allele_freqs, alleles,
                str(ref), str(f5), str(f3), pep_alleles, gene_loc, gene_id)

def cosmic_reader(cosmic_snps_path):
    """make a COSMIC SNP reader"""
    LOGGER.input_file(cosmic_snps_path, label='cosmic_snps_path')
    return reader(cosmic_snps_path, name="COSMICrec", rf=['Symbol', 'Location', 'Strand', 'Change', 'Histology'])

# for Click documentation see http://click.pocoo.org/
@click.command()
@click.option('--output_path', required=True, default=None, help='Path to write output file containing flanking seqs.')
@click.option('--chrom_dir', required=True, default=None, help='dir containing chrom fasta files.')
@click.option('--cosmic_snps_path', default="~gavin/COSMIC/derived_data/COSMIC_SNPs.txt", help='Path to file listing COSMIC SNP records.')
@click.option('--limit', required=False, type=int, default=None, help='Number of records to output.')
@click.option('--verbose', required=False, is_flag=True, default=False, help='prints volumes.')
def main(cosmic_snps_path, output_path, chrom_dir, limit, verbose):
    
    print "COSMIC processing - adding flanking data to COSMIC SNPs."
    app = CosmicRef(chrom_dir)
    flanksize = 300
    fcnt = 0
    nonsnp, noflank, snpmm = 0, 0, 0    # counts for problem SNP records
    
    if not os.path.isdir(output_path):
        raise ValueError("Create the output dir first")
    
    print
    print "processing COSMIC data from file:", cosmic_snps_path
    print "results stored in: %s/<primary histology>.txt.gz" % output_path
    
    runlog_path = os.path.join(output_path, 'COSMIC_flanks.log')
    LOGGER.log_file_path = runlog_path
    
    outfiles = {}
    outfile_paths = []
    rcnt = 0
    for snp in cosmic_reader(cosmic_snps_path):
        if snp.Change[-2]=='>': # check that it really is a SNP
            res = app.get_snp_dump_data(snp)
            if res:
                if snp.Histology not in outfiles:
                    outfile_path = os.path.join(output_path, "COSMIC_flanks_%s.txt.gz" % snp.Histology)
                    outfile_paths.append(outfile_path)
                    outfiles[snp.Histology] = gzip.open(outfile_path, "w")
                
                outfile = outfiles[snp.Histology]
                outfile.write('\t'.join(res)+'\n')
                fcnt += 1
            else:
                snpmm += 1
        else:
            nonsnp += 1
        
        rcnt += 1
        if limit and rcnt >= limit:
            break
        
        if rcnt % 100000 == 0:
            print "\r%s" % (" " * 80),
            print "\rread", rcnt, "records. non-SNPs =", nonsnp, "SNP errors =", snpmm,
            sys.stdout.flush()
    
    for cancer in outfiles:
        outfiles[cancer].close()
    
    for outfile_path in outfile_paths:
        LOGGER.output_file(outfile_path)
    
    print
    print "# flanks output:", fcnt
    if nonsnp:
        print "# non-SNP records in data file:", nonsnp
    if noflank:
        print "# SNPs without flank data:", noflank
    if snpmm:
        print "# SNP mismatches from reference:", snpmm
    
    print
    print"Done."
    return

if __name__=="__main__":
    main()
