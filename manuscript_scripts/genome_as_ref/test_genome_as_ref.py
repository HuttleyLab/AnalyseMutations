#!/usr/bin/env python
import os
from collections import defaultdict, Counter
from glob import glob

import click

from cogent3 import LoadTable


def get_matched_paths(local_dir, global_dir):
    local_fns = glob(os.path.join(local_dir, "*to*.txt"))
    global_fns = glob(os.path.join(global_dir, "*to*.txt"))
    matched = defaultdict(list)
    for fn in local_fns:
        base = os.path.basename(fn)
        matched[base].append(fn)

    for fn in global_fns:
        base = os.path.basename(fn)
        matched[base].append(fn)

    return list(matched.values())


def _load_mutated(path):
    t = LoadTable(path, sep="\t")
    t = t.filtered("mut=='M'")
    t = t.get_columns(t.header[:-1])
    counts = Counter()
    for r in t.tolist():
        c = r.pop(0)
        counts[tuple(r)] = c
    return counts


def compare_local_global(local_dir, global_dir):
    local_global = get_matched_paths(local_dir, global_dir)
    print("Number of paired files:", len(local_global))
    for l, g in local_global:
        loc = _load_mutated(l)
        glo = _load_mutated(g)
        motifs = set(loc) | set(glo)
        for motif in motifs:
            assert glo[motif] == loc[motif]

    return True


@click.command()
@click.argument("original", type=click.Path(exists=True))
@click.argument("new", type=click.Path(exists=True))
def main(original, new):
    result = compare_local_global(original, new)
    if result:
        print("Observed counts match perfectly between %s and %s" %
              (original, new))


if __name__ == "__main__":
    main()

