#!/usr/bin/env python2.7
import sys
import os
import os.path
import gzip
from collections import defaultdict
from multiprocessing import Pool
import functools
import cPickle as pik
from Bio import Seq
from Bio import SeqIO
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import seaborn
from scipy import stats
import pandas as pd
from pandas.io import pytables


def main(args):
    if ".fastq" in args.files[0]:
        process_fastq_files(args)
    if ".hf5" in args.files[0]:
        data = pytables.HDFStore(args.files[0])
        ac = data['allele_counts']
        oc = data['other_counts']
        meta = data['metadata']
        wt = ac[meta['Pos'] == 0]
        vswt = np.log2(ac / np.sum(wt))
        r1 = ['R1T0', 'R1T1', 'R1T2']
        pool = Pool(processes=args.numproc)
        result = pool.map(linregress_wrapper, vswt[r1].values)
        fitness = pd.DataFrame(data=result, index=vswt.index, columns=['slope', 'intercept', 'r', 'p', 'stderr'])
        pool.close()
        if len(args.files) == 2:
            output = pytables.HDFStore(args.files[1])
            output['fitness_r1'] = fitness
            output.close()
        data.close()
        if args.c:
            draw_hist(fitness)
        if args.m:
            idx = meta.reset_index().set_index(['Pos', 'AA', 'Codon'])
            aamap = compute_map(fitness, meta, idx, 'AA')
            codmap = compute_map(fitness, meta, idx, 'Codon')
            draw_maps(aamap, codmap, meta)
    return 0


def draw_hist(fitness):
    # cnt, bns = np.histogram(fitness['slope'], bins=20, normed=True)
    fig = plt.figure()
    # plt.bar(bns, cnt)
    plt.hist(fitness['slope'][np.isfinite(fitness['slope'])], bins=50, normed=True)
    fig.set_facecolor('white')
    ax = plt.gca()
    ax.set_xlabel('Fitness Score')
    ax.set_ylabel('Frequency')
    plt.savefig('fitnesshisto.png')
    # plt.show()


def draw_maps(aamap, codmap, meta):
    fig = plt.figure()
    plt.pcolor(aamap, cmap='RdBu')
    plt.xlim(0, aamap.shape[1])
    plt.ylim(0, aamap.shape[0])
    ax = plt.gca()
    fig.set_facecolor('white')
    ax.set_yticks([x+0.6 for x in xrange(0, aamap.shape[0])])
    aa = list(set(meta['AA']) - {None, np.nan})
    ax.set_yticklabels(aa)
    ax.set_ylabel('Residue')
    ax.set_xlabel('Ub Sequence Position')
    cb = plt.colorbar()
    cb.set_label('Relative Fitness')
    plt.savefig('aamap.png', bbox_inches='tight')

    fig = plt.figure()
    plt.pcolor(codmap, cmap='RdBu')
    plt.xlim(0, codmap.shape[1])
    plt.ylim(0, codmap.shape[0])
    ax = plt.gca()
    fig.set_facecolor('white')
    ax.set_yticks([x+0.6 for x in xrange(0, codmap.shape[0])])
    cod = list(set(meta['Codon']) - {None, np.nan})
    ax.set_yticklabels(cod)
    ax.set_ylabel('Codon')
    ax.set_xlabel('Ub Sequence Position')
    cb = plt.colorbar()
    cb.set_label('Relative Fitness')
    plt.savefig('codmap.png', bbox_inches='tight')


def compute_map(fitness, meta, idx, col, avgfun=np.nanmean):
    aa = list(set(meta[col]) - {None, np.nan})
    pos = list(set(meta['Pos']) - {None, np.nan, 0})
    aamap = np.zeros((len(aa), len(pos)))
    for i in xrange(0, len(aa)):
        for p in xrange(0, len(pos)):
            aamap[i, p] = avgfun(fitness.loc[idx.xs([pos[p], aa[i]], level=['Pos', col])['Seq']]['slope'])
    return aamap


def linregress_wrapper(y):
    slope, intercept, rval, pval, err = stats.linregress(x=np.array([0, 1.87, 3.82]), y=y)
    return slope, intercept, rval, pval, err


def process_fastq_files(args):
    store = pytables.HDFStore(args.out)
    metadata = parse_library(args.dbfile)
    store['metadata'] = metadata
    pool = Pool(processes=args.numproc)
    result = pool.map(functools.partial(process_fastq, prepop=metadata.index), args.files)
    pool.close()
    allele_temp = {file2col(f): r[0] for f, r in zip(args.files, result)}
    other_temp = {file2col(f): r[1] for f, r in zip(args.files, result)}
    allele_counts = pd.DataFrame.from_dict(allele_temp, orient="columns")
    other_counts = pd.DataFrame.from_dict(other_temp, orient="columns")
    allele_counts.index.rename('Seq', inplace=True)
    other_counts.index.rename('Seq', inplace=True)
    store['allele_counts'] = allele_counts
    store['other_counts'] = other_counts
    store.close()


def parse_library(dbfile):
    library = pik.load(open(dbfile, 'rb'))
    allele = {}
    for k in library:
        s = str(Seq.Seq(k, Seq.Alphabet.SingleLetterAlphabet()).reverse_complement())
        pos = int(library[k][0])
        if pos != 0:
            cod = Seq.Seq(library[k][1], Seq.Alphabet.SingleLetterAlphabet())
            aa = cod.transcribe().translate()
            allele[s] = [pos, str(cod), str(aa)]
        else:
            allele[s] = [pos, None, None]
    metadata = pd.DataFrame(data=([k]+v for k, v in allele.items()),
                            columns=['Seq', 'Pos', 'Codon', 'AA']).set_index('Seq')
    return metadata


def file2col(x):
    return os.path.basename(x).replace('.fastq', '').replace('.gz', '')


def process_fastq(fastq, prepop):
    alleles = {k: 0 for k in prepop}
    others = defaultdict(int)
    parse_file(fastq, alleles, others)
    return alleles, dict(others)


def parse_file(fastq, counts, others):
    if fastq.endswith('.gz'):
        f = gzip.open(fastq, 'rb')
    else:
        f = open(fastq, 'rU')
    records = SeqIO.parse(f, 'fastq', Seq.Alphabet.DNAAlphabet())
    for r in records:
        putative = str(r.seq[0:18])
        if putative in counts:
            counts[putative] += 1
            # r.letter_annotations['phred_quality']
        else:
            others[putative] += 1
            # r.letter_annotations['phred_quality']
    f.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Map read sequences to known barcodes.")
    parser.add_argument("-p", "--processes", action="store", type=int, default=1,
                        dest="numproc", help="Number of parallel processes.")
    parser.add_argument("-d", "--database", action="store", type=str,
                        default=None, dest="dbfile", help="Barcode definitions or HDF5 file.")
    parser.add_argument("-m", "--map", action="store_true", dest="m", default=False,
                        help="Compute and display fitness heatmaps.")
    parser.add_argument("-c", "--counts", action="store_true", dest="c", default=False,
                        help="Compute and display fitness histogram.")
    parser.add_argument("-o", "--output", action="store", type=str,
                        dest="out", help="Output Pickle with barcode definitions and frequencies.")
    parser.add_argument("files", nargs='+')
    sys.exit(main(parser.parse_args()))
