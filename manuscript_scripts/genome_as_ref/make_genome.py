#!/usr/bin/env python
from subprocess import call
from cogent3.util.misc import open_
from mutation_motif.util import abspath

from tqdm import tqdm

def load_seq(path):
    with open_(path, 'rt') as infile:
        label = next(infile).strip()
        frags = []
        for i, line in enumerate(infile):
            if line.startswith('>'):
                raise ValueError(line)
            line = line.strip()
            frags.append(line)

    if len(frags) > 1:
        raise ValueError(len(frags))

    return frags[0]

fns = \
['/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr1.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr2.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr3.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr4.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr5.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr6.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr7.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr8.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr9.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr10.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr11.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr12.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr13.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr14.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr15.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr16.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr17.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr18.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr19.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr20.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr21.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/Chr22.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/ChrX.fa.gz',
 '/Users/gavin/DevRepos/AnalyseMutationsBloated/data/chroms79/ChrY.fa.gz']


outfile = "../../data/chroms79/genome.fa"
with open(outfile, "wt") as out:
    out.write(">human\n")
    for fn in tqdm(fns, desc="Loading chromosomes", ncols=80):
        seq = load_seq(fn)
        out.write(seq)

call(["gzip", outfile])

print("Done")
