"""export seq files for different mutation types"""
from __future__ import division

import os, sys, time, re
from itertools import permutations
from optparse import make_option
from cogent.util.option_parsing import parse_command_line_parameters

from sumatra.projects import load_project
from sumatra.programs import get_executable

from mutation_motif.util import open_, create_path, abspath, just_nucs, load_from_fasta
from mutation_motif import profile, motif_count

fn_suffixes = re.compile(r"\.(fa|fasta)\.(gz|gzip|bz2)$")

def get_counts_filename(align_path, output_dir):
    """returns counts output path
    
    Arguments:
        - align_path: path to the alignment file. The basename will be modified to use a .txt suffix
        - output_dir: directory where the counts file is to be written
    """
    
    fn = os.path.basename(align_path)
    fn = fn_suffixes.sub(".txt", fn)
    counts_filename = os.path.join(output_dir, fn)
    return counts_filename

def align_to_counts(opts):
    '''returns counts table from alignment of sequences centred on a SNP'''
    
    if len(args) > 1:
        raise RuntimeError("too many positional args")
    
    if not opts.dry_run:
        create_path(opts.output_path)
    
    print "Deriving counts from sequence file"
    direction = tuple(opts.direction.split('to'))
    chosen_base = direction[0]
    orig_seqs = load_from_fasta(os.path.abspath(opts.alignfile))
    seqs = orig_seqs.ArraySeqs
    seqs = just_nucs(seqs)
    orig, ctl = profile.get_profiles(seqs, chosen_base=chosen_base, step=1,
                                     flank_size=opts.flank_size, seed=opts.seed)
    
    # convert profiles to a motif count table
    orig_counts = motif_count.profile_to_seq_counts(orig, flank_size=opts.flank_size)
    ctl_counts = motif_count.profile_to_seq_counts(ctl, flank_size=opts.flank_size)
    counts_table = motif_count.get_count_table(orig_counts, ctl_counts, opts.flank_size*2)
    counts_table = counts_table.sorted(columns='mut')
    return counts_table


script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = "export tab delimited counts table from "\
"alignment centred on SNP position. Output file is written to the same "\
"path with just the file suffix changed from fasta to txt."

script_info['required_options'] = [
     make_option('-a','--alignfile', help='fasta aligned file centred on mutated position.'),
     make_option('-o','--output_path', help='Path to write data.'),
     make_option('-f','--flank_size', help='Number of bases per side to include.'),
     make_option('--direction', default=None,
                 choices=['AtoC', 'AtoG', 'AtoT', 'CtoA', 'CtoG', 'CtoT', 'GtoA', 'GtoC',
                          'GtoT', 'TtoA', 'TtoC', 'TtoG'], help='Mutation direction.'),
    ]

script_info['optional_options'] = [
    make_option('-S', '--seed',
        help='Seed for random number generator (e.g. 17, or 2015-02-13). Defaults to system time.'),
    make_option('-D','--dry_run', action='store_true', default=False,
        help='Do a dry run of the analysis without writing output.'),
    make_option('-r','--reason', help='Reason for running analysis (for Sumatra log).'),
    ]

script_info['version'] = '0.1'
script_info['authors'] = 'Gavin Huttley'

if __name__ == "__main__":
    option_parser, opts, args =\
       parse_command_line_parameters(disallow_positional_arguments=False, **script_info)
    
    if not opts.seed:
        opts.seed = str(time.time())
        print "NOTE: set random number seed to '%s'" % (opts.seed)
    
    # record run using sumatra
    project = load_project()
    
    # determine the path to input data relative to sumatra's input data store
    input_path = abspath(opts.alignfile)
    if project.input_datastore.root not in input_path:
        raise ValueError("input path not nested under sumatra input path")
    
    input_path = input_path.replace(project.input_datastore.root, '')
    if input_path.startswith('/'):
        input_path = input_path[1:]
    
    # generate unique hash key for input files, then create sumatra record
    input_data = project.input_datastore.generate_keys(input_path)
    record = project.new_record(parameters=opts, input_data=input_data,
                                executable=get_executable(script_file=__file__),
                                main_file=__file__,
                                reason=opts.reason)
    
    # the sumatra_label will be used to create a subdirectory within the user
    # specified output path so that output from simultaneous runs are
    # associated with the correct sumatra record
    opts.sumatra_label = record.label
    opts.output_path = os.path.join(opts.output_path, record.label)
    start_time = time.time()
    
    # run the program
    counts_table = align_to_counts(opts)
    counts_filename = get_counts_filename(opts.alignfile, opts.output_path)
    if not opts.dry_run:
        counts_table.writeToFile(counts_filename, sep='\t')
    
    
    # determine runtime, output file identifiers and save the sumatra record
    record.datastore.root = opts.output_path
    record.duration = time.time() - start_time
    record.output_data = record.datastore.find_new_data(record.timestamp)
    project.add_record(record)
    project.save()
