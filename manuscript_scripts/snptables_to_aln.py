"""export seq files for different mutation types"""
from __future__ import division

import os, sys, time
from itertools import permutations

import click

from mutation_motif.util import open_, create_path, abspath
from scitrack import CachingLogger

import strand

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2015, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

LOGGER = CachingLogger(create_dir=True)

MAF = 0.05

def get_gc_freq(seq):
    l = len(seq)
    GC = sum([seq.count(base) for base in 'GC'])
    freq = GC / l
    return freq

def MakeFreqCompare(cuttoff, ge, get_freq=None, verbose=False):
    """returns a function that takes both flanks, returns True if base freq
    satisfies the specified inequality
    
    Arguments:
        - cuttoff: a proportion of G+C
        - ge: if ge=True, returns result of seq_gc >= cuttoff, otherwise
          returns seq_gc < cuttoff
        - get_freq: a callback function for computing the freq to be used"""
    if ge:
        func = lambda x: x >= cuttoff
    else:
        func = lambda x: x < cuttoff
    
    def call(data):
        if get_freq:
            freq = get_freq(data)
        else:
            freq = data
        result = func(freq)
        if verbose:
            print cuttoff, freq, result
        
        return result
    
    return call

def MakeStranded(stranded):
    """returns a function that reverse complements if seq_type is not-intergenic"""
    
    def unmodified(gene_strand, snp_strand, alleles, ancestor, freqs, flank_5, flank_3):
        """makes the seq"""
        return alleles, ancestor, flank_5 + ancestor + flank_3
    
    def reverse_complementable(gene_strand, snp_strand, alleles, ancestor, freqs, flank_5, flank_3):
        """makes the seq"""
        if gene_strand and \
         strand.reverse_complement_record(gene_strand, snp_strand):
            alleles, ancestor, freqs, flank_5, flank_3 = strand.get_rc_record(alleles, ancestor, freqs, flank_5, flank_3)
        
        return alleles, ancestor, flank_5 + ancestor + flank_3
    
    if stranded:
        func = reverse_complementable
    else:
        func = unmodified
    
    return func

def is_autosome(chrom):
    """docstring for is_autosome"""
    chrom = chrom.upper()
    return 'X' not in chrom and 'Y' not in chrom

def is_xchrom(chrom):
    chrom = chrom.upper()
    return 'X' in chrom

everything = lambda x: True

def filtered_records(records, direction, seen, chroms, correct_chrom=everything, correct_freq=everything, correct_comp=everything, stranded=False, verbose=True):
    """parses records, yielding when condition met"""
    adjust_strand = MakeStranded(stranded)
    progress = ""
    for record_num, line in enumerate(records):
        if record_num % 1000 == 0:
            progress = ' ' * len(progress)
            sys.stdout.write(progress + "\r")
            progress = 'records read = %d\r' % (record_num)
            sys.stdout.write(progress)
            sys.stdout.flush()
        
        
        line = line.strip().split('\t')
        label = line[0].strip()
        if label in seen:
            continue

        coord = line[1]

        ancestor = line[6]
        if len(ancestor) > 1: # it's the string "None"
            continue

        alleles = eval(line[5])
        if ancestor not in alleles or set(direction) != alleles:
            continue

        if len(alleles) != 2:
            continue

        if not correct_chrom(coord):
            continue

        freqs = eval(line[4])
        if freqs is not None:
            freqs = dict(freqs)
        else:
            freqs = {}

        if not correct_freq(freqs.values()):
            continue
        
        flank_5, flank_3 = line[7], line[8]
        
        if 'N' in flank_5 or 'N' in flank_3:
            continue
        
        # handle the strandedness
        snp_strand = int(line[2])
        try:
            # e.g. Homo sapiens:chromosome:20:68350-77174:1
            gene_strand = int(line[11].split(':')[-1])
        except (IndexError, ValueError):
            gene_strand = None
        
        alleles, ancestor, seq = adjust_strand(gene_strand, snp_strand, alleles,
                                    ancestor, freqs, flank_5, flank_3)
        
        got = alleles.difference(ancestor).pop()
        if (ancestor, got) != direction:
            continue
        
        seen.update([label])

        if not correct_comp(seq):
            continue

        
        record = '\n'.join(['>%s' % label, seq, ''])
        yield record
        
        chroms.update([coord.split(':')[2]])
    

def run(input_path, output_path, direction, prefix, chrom_class, gc_class, freq_class, adjust_strand, limit, force_overwrite, dry_run, verbose):
    if not dry_run:
        create_path(output_path)
    
    correct_freq = {'All': everything,
        'Common': MakeFreqCompare(MAF, ge=True, get_freq=min,
                                    verbose=verbose),
        'Rare': MakeFreqCompare(MAF, ge=False, get_freq=min,
                                    verbose=verbose)}[freq_class]
    
    correct_comp = {'All': everything,
        'Hi': MakeFreqCompare(0.5, ge=True, get_freq=get_gc_freq,
                                    verbose=verbose),
        'Lo': MakeFreqCompare(0.4, ge=False, get_freq=get_gc_freq,
                                    verbose=verbose)}[gc_class]
    
    correct_chrom = {'All': everything, 'A': is_autosome,
                     'X': is_xchrom}[chrom_class]
    
    seen = set()
    chroms = set()
    if not os.path.exists(input_path):
        raise IOError("no files matching %s" % input_path)
    
    name_components = dict(freq_class='freq_'+freq_class,
            chrom_class='chrom_'+chrom_class,
            gc_class='GC_'+gc_class, direction=direction,
            prefix = prefix or '')
    
    outfilename = os.path.join(output_path,
    '%(prefix)s%(freq_class)s-%(chrom_class)s-%(gc_class)s-%(direction)s.fasta.gz' % name_components)
    
    runlog_path = os.path.join(output_path,
    '%(prefix)s%(freq_class)s-%(chrom_class)s-%(gc_class)s-%(direction)s.log' % name_components)
    LOGGER.log_file_path = runlog_path
    
    if not force_overwrite and (os.path.exists(outfilename) or os.path.exists(runlog_path)):
        msg = "Either %s or %s already exist. Force overwrite of existing files with -F."
        raise ValueError(msg % (outfilename, runlog_path))
 
    LOGGER.input_file(input_path)
    
    direction = tuple(direction.split('to'))
    
    with open_(input_path) as infile:
        with open_(outfilename, 'w') as outfile:
            num = 0
            for record in filtered_records(infile, direction, seen, chroms,
             correct_chrom=correct_chrom, correct_freq=correct_freq,
             correct_comp=correct_comp, stranded=adjust_strand, verbose=False):
                outfile.write(record)
                num += 1
                if limit and num >= limit:
                    break
        
        LOGGER.output_file(outfilename)
    msg = "Wrote %d records to %s" % (num, outfilename)
    print msg
    LOGGER.log_message(msg + "\n", label="completed")

@click.command()
@click.option('-i','--input_path', required=True, help='glob pattern to data files.')
@click.option('-o','--output_path', required=True, help='Path to write data.')
@click.option('--direction', default=None, required=True,
        type=click.Choice(['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT',
                'GtoA', 'GtoC', 'GtoT', 'TtoA', 'TtoC', 'TtoG']),
                help='Mutation direction.')
@click.option('-p','--prefix', help='Prefix for output file, e.g. intronic-')
@click.option('-c','--chrom_class', default='All',
        type=click.Choice(['All', 'X', 'A', 'Y']), help='Chrom class.')
@click.option('--GC', default='All',
        type=click.Choice(['All', 'Hi', 'Lo']),
        help='GC proportion. Hi is >0.5, Lo is < 0.4.')
@click.option('--freq_class', default='All',
        type=click.Choice(['All', 'Common', 'Rare']),
        help='Frequency class. Common has MAF >0.05, rare <= 0.05.')
@click.option('-a','--adjust_strand', is_flag=True,
        help='Reverse complements records whose gene and snp strands differ.')
@click.option('-l','--limit', default=None, type=int,
        help='Number of results to return.')
@click.option('-F', '--force_overwrite', is_flag=True, help='Overwrite existing files.')
@click.option('-D', '--dry_run', is_flag=True, help='Do a dry run of the analysis without writing output.')
@click.option('--verbose', is_flag=True, help='Verbose output.')
def main(input_path, output_path, direction, prefix, chrom_class, gc, freq_class, adjust_strand, limit, force_overwrite, dry_run, verbose):
    """export fasta formatted seqs matching specified conditions."""
    if not dry_run:
        LOGGER.log_message("%s" % locals(), label="vars")
    
    input_path = abspath(input_path)
    output_path = abspath(output_path)
    
    start_time = time.time()
    
    run(input_path, output_path, direction, prefix, chrom_class, gc, freq_class, adjust_strand, limit, force_overwrite, dry_run, verbose)
    
    # determine runtime
    duration = time.time() - start_time
    if not dry_run:
        LOGGER.log_message("%.2f" % (duration/60.), label="run duration (minutes)")

if __name__ == "__main__":
    main()