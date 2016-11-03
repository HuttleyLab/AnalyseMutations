#!/usr/bin/env python

"""creates counts tavbles, using a genomic control"""

import os
from glob import glob
import time
from random import randint, seed as set_seed, choice

import click
from tqdm import tqdm
from cogent3 import DNA
from cogent3.core.alignment import ArrayAlignment
from cogent3.util.misc import open_
from scitrack import CachingLogger

from mutation_motif.util import makedirs, abspath, just_nucs,\
    load_from_fasta
from mutation_motif import profile, motif_count, aln_to_counts

LOGGER = CachingLogger(create_dir=True)

_bases = set(DNA)


def just_nucs_str(seq):
    return set(seq) <= _bases


def load_seq(path):
    LOGGER.input_file(path)
    with open_(path, 'rt') as infile:
        next(infile).strip()
        frags = []
        for i, line in enumerate(infile):
            if line.startswith('>'):
                raise ValueError("should be only 1 sequence per file", line)
            line = line.strip()
            frags.append(line)

    if len(frags) > 1:
        raise ValueError(len(frags))

    return frags[0]


def GenomeSampler(seq, flank_size):
    """factory function for returning random sequences centered on a specified
    base

    This is compute inefficient, making multiple random selections until one
    matches the chosen base.
    """
    seq_length = len(seq)

    def call(base):
        while True:
            index = randint(flank_size, seq_length - flank_size - 1)
            if seq[index] != base:
                continue

            subseq = seq[index - flank_size: index + flank_size + 1]
            if not just_nucs_str(subseq):
                continue

            break

        return subseq

    return call


# inline test that random sampling
__expect = {'A': ["AAA", "AAG"], 'C': ["CCC"],
            'G': ["AGG", "GGG", "GGT"], 'T': ["GTT", "TTT", "TTC"]}
__seq = "AAAAGGGGTTTNCCCC"
__genome_sampler = GenomeSampler(__seq, 1)
for i in range(100):  # just repeat a few times
    for base in __expect:
        got = __genome_sampler(base)
        assert got in __expect[base], \
            "base=%s; got=%s; possibles=%s" % (base, got, __expect[base])


def get_genome_control(num, genome_sampler, chosen_base, step, flank_size,
                       limit, sample_indices=None, circle_range=None,
                       seed=None):
    """returns control profile with sampling from whole genome

    num: number of controls to produce"""
    assert seed is not None, "Must provide a random number seed"
    set_seed(seed)
    seqs = []
    if limit is None:
        # nested in a conditional so we don't evaluate the contional millions
        # of times when limit is None
        for i in range(num):
            seq = genome_sampler(chosen_base)
            seqs.append(['s%s' % i, seq])
    else:
        for i in range(limit):
            seq = genome_sampler(chosen_base)
            seqs.append(['s%s' % i, seq])

    data = ArrayAlignment(data=seqs, moltype=DNA)
    if limit is not None:
        print(data)

    return data.array_seqs


def counts_with_genomic_control(align_path, genome_genome_sampler,
                                output_path, flank_size, direction,
                                step, limit, seed, dry_run):
    '''returns counts table from alignment of sequences centred on a SNP'''
    # this differs from the mutation_motif implementation by being based on
    # strings, rather than array
    if not dry_run:
        makedirs(output_path)

    step = int(step)

    direction = tuple(direction.split('to'))
    chosen_base = direction[0]
    orig_seqs = load_from_fasta(os.path.abspath(align_path))
    seqs = orig_seqs.array_seqs
    seqs = just_nucs(seqs)
    orig = profile.get_observed(seqs, flank_size=flank_size)
    del(seqs)  # for memory issues
    ctl = get_genome_control(orig.shape[0], genome_genome_sampler,
                             chosen_base=chosen_base, step=step,
                             flank_size=flank_size, limit=limit, seed=seed)
    if limit:
        return ctl

    # convert profiles to a motif count table
    orig_counts = motif_count.profile_to_seq_counts(orig,
                                                    flank_size=flank_size)
    ctl_counts = motif_count.profile_to_seq_counts(ctl,
                                                   flank_size=flank_size)
    counts_table = motif_count.get_count_table(
        orig_counts, ctl_counts, flank_size * 2)
    counts_table = counts_table.sorted(columns='mut')
    return counts_table


@click.command()
@click.option("-a", "--align_pattern", required=True,
              help="glob pattern for fasta aligned files centred on mutated "
              "position.")
@click.option("-g", "--genome_path", required=True,
              help="fasta file of whole genome.")
@click.option('-o', '--output_path', required=True, help='Path to write data.')
@click.option('-f', '--flank_size', required=True, type=int,
              help='Number of bases per side to include.')
@click.option('-S', '--seed',
              help='Seed for random number generator (e.g. 17, or 2015-02-13).'
              ' Defaults to system time.')
@click.option('--step', default='1', type=click.Choice(['1', '2', '3']),
              help='Specifies a "frame" for selecting the random base.')
@click.option('-l', '--limit', type=int, default=None,
              help='Limit the number of records processed per file, causes no data to be written.')
@click.option('-D', '--dry_run', is_flag=True,
              help='Do a dry run of the analysis without writing output.')
@click.option('-F', '--force_overwrite', is_flag=True,
              help='Overwrite output and run.log files.')
def main(align_pattern, genome_path, output_path, flank_size, seed,
         step, limit, dry_run, force_overwrite):
    """Export tab delimited counts table from alignment centred on SNP position.

    Output file is written to the same path with just the file suffix changed
    from fasta to txt."""
    args = locals()
    if not seed:
        seed = str(time.time())

    align_pattern = abspath(align_pattern)
    align_paths = glob(align_pattern)
    genome_path = abspath(genome_path)

    click.secho("Loading genome", fg='green')
    genomeseq = load_seq(genome_path)
    genome_genome_sampler = GenomeSampler(genomeseq, flank_size)

    output_path = abspath(output_path)
    start_time = time.time()

    if not dry_run:
        makedirs(output_path)

        LOGGER.log_message(str(args), label='vars')
        LOGGER.log_message(str(seed), label="random_number_seed")
        runlog_path = os.path.join(output_path,
                                   "counts_with_genomic_control.log")
        LOGGER.log_file_path = runlog_path

    if limit is None:
        align_paths = tqdm(align_paths, ncols=80)

    for align_path in align_paths:
        if limit is not None:
            print('\n', align_path)

        x = align_path.find('to')
        assert x > 0, x
        direction = align_path[x - 1: x + 3]

        counts_filename = aln_to_counts.get_counts_filename(
            align_path, output_path)
        if not force_overwrite and os.path.exists(counts_filename):
            msg = "%s already exist. Force overwrite of existing"\
                  " files with -F."
            raise ValueError(msg % (counts_filename))

        LOGGER.input_file(align_path, label="align_path")

        counts_table = counts_with_genomic_control(align_path,
                                                   genome_genome_sampler,
                                                   output_path,
                                                   flank_size,
                                                   direction, step, limit,
                                                   seed, dry_run)

        if not dry_run and limit is None:
            counts_table.write(counts_filename, sep='\t')
            LOGGER.output_file(counts_filename)

    # determine runtime
    duration = time.time() - start_time
    click.secho("run duration (minutes)=%.2f" % (duration / 60.))


if __name__ == "__main__":
    main()
