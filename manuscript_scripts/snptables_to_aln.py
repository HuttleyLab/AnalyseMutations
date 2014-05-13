"""export seq files for different mutation types"""
from __future__ import division

import os
import glob
from itertools import permutations
from optparse import make_option
from cogent.util.option_parsing import parse_command_line_parameters
from mutation_motif.util import open_, create_path
import strand


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
    for line in records:
        line = line.strip().split('\t')
        label = line[0].strip()
        if label in seen:
            continue

        coord = line[1]

        ancestor = line[6]
        if len(ancestor) > 1: # it's the string "None"
            continue

        alleles = line[5]
        if ancestor not in alleles:
            continue

        alleles = set(alleles.split('/'))
        if len(alleles) != 2:
            continue

        if not correct_chrom(coord):
            continue

        freqs = dict(eval(line[4]))
        if not freqs:
            continue

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
        except IndexError:
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
    

def main(script_info):
    option_parser, opts, args =\
       parse_command_line_parameters(disallow_positional_arguments=False, **script_info)
    
    if len(args) > 1:
        raise RuntimeError("too many positional args")
    
    if args:
        # have a command line label from sumatra which requires we put all
        # results in a sub-directory of the same name
        outpath = os.path.join(opts.outpath, args[0])
    else:
        outpath = opts.outpath
    
    if not opts.dry_run:
        create_path(outpath)
    
    chrom_class = opts.chrom_class
    freq_class = opts.freq_class
    correct_freq = {'All': everything,
        'Common': MakeFreqCompare(MAF, ge=True, get_freq=min,
                                    verbose=opts.verbose),
        'Rare': MakeFreqCompare(MAF, ge=False, get_freq=min,
                                    verbose=opts.verbose)}[freq_class]
    
    gc_class = opts.GC
    correct_comp = {'All': everything,
        'Hi': MakeFreqCompare(0.5, ge=True, get_freq=get_gc_freq,
                                    verbose=opts.verbose),
        'Lo': MakeFreqCompare(0.4, ge=False, get_freq=get_gc_freq,
                                    verbose=opts.verbose)}[gc_class]
    
    chrom_class = opts.chrom_class
    correct_chrom = {'All': everything, 'A': is_autosome,
                     'X': is_xchrom}[chrom_class]
    
    direction = tuple(opts.direction.split('to'))
    seen = set()
    chroms = set()
    if not os.path.exists(opts.input_path):
        raise IOError("no files matching %s" % opts.input_path)
    
    name_components = dict(freq_class='freq_'+freq_class,
            chrom_class='chrom_'+chrom_class,
            gc_class='GC_'+gc_class, direction=opts.direction,
            prefix = opts.prefix or '')
    
    outfilename = os.path.join(outpath,
    '%(prefix)s%(freq_class)s-%(chrom_class)s-%(gc_class)s-%(direction)s.fasta.gz' % name_components)
    with open_(opts.input_path) as infile:
        with open_(outfilename, 'w') as outfile:
            num = 0
            for record in filtered_records(infile, direction, seen, chroms,
             correct_chrom=correct_chrom, correct_freq=correct_freq,
             correct_comp=correct_comp, stranded=opts.adjust_strand, verbose=False):
                outfile.write(record)
                num += 1
                if opts.limit and num >= opts.limit:
                    break

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = "export fasta formatted seqs matching specified conditions."


script_info['required_options'] = [
     make_option('-i','--input_path', help='glob pattern to data files.'),
     make_option('-o','--outpath', help='Path to write data.'),
     make_option('-p','--prefix', help='Prefix for output figure, e.g. intronic-'),
    make_option('--direction', default=None,
     choices=['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT', 'GtoA', 'GtoC',
             'GtoT', 'TtoA', 'TtoC', 'TtoG'], help='Mutation direction.'),
    ]

script_info['optional_options'] = [
    make_option('-c','--chrom_class', type='choice', default='All',
        choices=['All', 'X', 'A', 'Y'], help='Chrom class [default:%default].'),
    make_option('--GC', type='choice', default='All',
        choices=['All', 'Hi', 'Lo'], help='GC proportion. Hi is >0.5, Lo is < 0.4.'),
    make_option('--freq_class', type='choice', default='All',
        choices=['All', 'Common', 'Rare'],
        help='Frequency class. Common has MAF >0.05, rare <= 0.05.'),
    make_option('-l','--limit', default=None, type=int,
        help='Number of results to return.'),
    make_option('-a','--adjust_strand', action='store_true', default=False,
        help='Reverse complements records whose gene and snp strands differ.'),
    make_option('-D','--dry_run', action='store_true', default=False,
        help='Do a dry run of the analysis without writing output.'),
    ]

script_info['version'] = '0.1'
script_info['authors'] = 'Gavin Huttley'

if __name__ == "__main__":
    main(script_info)
    