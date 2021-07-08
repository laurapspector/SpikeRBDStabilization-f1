import gzip
import itertools
import multiprocessing
from operator import itemgetter
import os

import numpy as np
import pandas as pd

from dms.arguments import parse_args_and_read_config
from dms.merge import merge_all_reads
from dms.mutation import AminoAcidMutation, Mutation, WildType, is_wt
from dms.tile import mutations_in_seq


def open_by_extension(path, mode):
    """Open a file using gzip.open if its name ends with '.gz', otherwise
    use open."""
    return (gzip.open if path.endswith('gz') else open)(path, mode)

def mutation_counts(seqs, tile):
    """Count up mutation sets in an iterable of sequences.

    seqs: An iterable of DNA sequences.
    tile: Tile object describing the DNA sequences.

    Returns: A dict mapping a tuple of Mutations and
    NontargetMutations to count.
    """
    seq_counts = {}
    for seq in seqs:
        seq_counts[seq] = seq_counts.get(seq, 0) + 1
    counts = {}
    for seq, n in seq_counts.items():
        muts = mutations_in_seq(tile, seq)
        counts[muts] = counts.get(muts, 0) + n
    return counts

def library_statistics(tile, counts):
    n_muts = {}
    others = 0
    total = sum(counts.values())
    j = 0
    for muts, n in counts.items():
        # Are all mutations Mutation (not NontargetMutation) objects
        # and are all at targeted tile positions?
        if all(isinstance(mut, Mutation) for mut in muts) and \
           all(mut.pos in tile.positions for mut in muts):
            n_muts[len(muts)] = n_muts.get(len(muts), 0) + n
        # Otherwise, consider them 'other' mutations.
        else:
            others += n
    return total, n_muts, others

def collapsed_and_filtered_counts(tile, counts):
    total = sum(counts.values())
    new_counts = {}
    for muts, n in counts.items():
        if len(muts) == 0:
            mut = WildType
        elif len(muts) == 1 and isinstance(muts[0], Mutation):
            # This will either create an AminoAcidMutation or WildType
            # depending on whether the codon change is synonymous or
            # not.
            mut = AminoAcidMutation.from_mutation(muts[0])
            # We don't want to record this mutation if it is at a
            # position not specified by the tile.
            if (not is_wt(mut)) and (mut.pos not in tile.positions):
                continue
        else:
            continue
        new_counts[mut] = new_counts.get(mut, 0) + n
    return total, new_counts

def get_stats_and_counts(path1, path2, tile, params):
    """Get statistics and mutation counts from two paired-end read FASTQ
    files.

    path1: path to forward read FASTQ file.
    path2: path to reverse read FASTQ file.
    tile: a Tile object describing the amplicon.
    params: parameter dict.

    Returns a tuple (stats, total, counts) where stats is a tuple from
    library_statistics, total is the total number of reads in the
    sample, and counts is a dict mapping mutation => count.

    If a filename ends in '.gz' it will be assumed to be gzipped,
    otherwise it will be assumed to be plain text.
    """
    f1 = open_by_extension(path1, 'rt')
    f2 = open_by_extension(path2, 'rt')
    reads = merge_all_reads(f1, f2,
                            tile.length,
                            params.max_mismatches,
                            params.min_quality)
    #reads = itertools.islice(reads, 10000)
    raw_counts = mutation_counts(reads, tile)
    stats = library_statistics(tile, raw_counts)
    total, counts = collapsed_and_filtered_counts(tile, raw_counts)
    f1.close()
    f2.close()
    return stats, total, counts

def process_all_samples(params, tiles, samples):
    inputs = []
    for sample, (tile_name, filenames) in samples.items():
        tile = tiles[tile_name]
        path1, path2 = [os.path.join(params.fastq_file_dir, f)
                        for f in filenames]
        inputs.append((path1, path2, tile, params))
    if params.use_multiprocessing:
        with multiprocessing.Pool() as pool:
            results = pool.starmap(get_stats_and_counts, inputs)
    else:
        results = list(itertools.starmap(get_stats_and_counts, inputs))
    stats = dict(zip(samples, map(itemgetter(0), results)))
    counts = dict(zip(samples, map(itemgetter(1, 2), results)))
    return stats, counts

def process_all_experiments(params, tiles, samples, experiments, counts):
    results = {}
    for experiment, (ref_sample, sel_sample) in experiments.items():
        tile_name = samples[ref_sample][0]
        tile = tiles[tile_name]

        ref_total, ref_counts = counts[ref_sample]
        sel_total, sel_counts = counts[sel_sample]

        # Remove mutations that don't have enough reference counts.
        muts = sorted([m for (m, n) in ref_counts.items()
                       if n >= params.min_ref_counts],
                      key=AminoAcidMutation.sort_key)

        if not any(is_wt(m) for m in muts):
            print('WARNING: The wild-type sequence will not appear in'
                  f' experiment {experiment}.')

        d = pd.DataFrame({
            'experiment' : experiment,
            'variant' : muts,
            'sel_counts' : [sel_counts.get(m, params.pseudocount)
                            for m in muts],
            'sel_total' : sel_total,
            'ref_counts' : [ref_counts[m] for m in muts],
            'ref_total' : ref_total})
        d['ER'] = np.log2((d['sel_counts'] / d['sel_total'])
                          / (d['ref_counts'] / d['ref_total']))
        results[experiment] = d
    return results

def all_reference_samples(experiments):
    return sorted({ref_sample for (_, (ref_sample, _)) in experiments.items()})

def write_stats(stats, path):
    total, n_muts, others = stats
    with open(path, 'wt') as f:
        for n in sorted(n_muts):
            print(f'{n}\t{n_muts[n]}\t{100*n_muts[n]/total:5.2f}', file=f)
        print(f'others\t{others}\t{100*others/total:5.2f}', file=f)
        print(f'total\t{total}', file=f)


def main(argv):
    if not os.path.exists('Output'):
        os.makedirs('Output')

    params, tiles, samples, experiments, proteins = \
        parse_args_and_read_config(argv)

    stats, counts = process_all_samples(params, tiles, samples)
    # stats and counts both have total reads as their first elements,
    # this is just a sanity check.
    assert set(stats.keys()) == set(counts.keys())
    for sample in stats:
        assert stats[sample][0] == counts[sample][0]

    for sample in all_reference_samples(experiments):
        out_path = os.path.join(params.output_dir, 'Output', f'{sample}_stats.tsv')
        write_stats(stats[sample], out_path)

    data = process_all_experiments(params, tiles, samples, experiments, counts)
    for protein, exps in proteins.items():
        out_path = os.path.join(params.output_dir, 'Output', f'{protein}_counts.csv')
        d = pd.concat([data[exp] for exp in exps])\
              .reset_index(drop=True)\
              .to_csv(out_path, index=False)
