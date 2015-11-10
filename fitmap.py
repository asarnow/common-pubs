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
    if args.out is None:
        args.out = "."
    if len(args.files) > 0 and ".fastq" in args.files[0]:
        if len(args.files) > 1:
            meta, ac, oc = process_fastq_files(args)
            expdefs = { None: [(file2col(f),) for f in args.files]}
        else:
            if args.expdefs is None:
                print "Experiment definitions required."
                return 1
            with open(args.expdefs, 'rb') as f:
                expdefs = pik.load(f)
            meta, ac, oc = process_fastq_files(args, [v[0] for v in [expdefs[k] for k in expdefs]])
        if args.cache is not None:
            store = pytables.HDFStore(args.out)
            store['metadata'] = meta
            store['allele_counts'] = ac
            store['other_counts'] = oc
            store.close()
    elif args.cache is not None:
        if args.expdefs is None:
            print "Experiment definitions required."
            return 1
        with open(args.expdefs, 'rb') as f:
            expdefs = pik.load(f)
        cache = pytables.HDFStore(args.cache)
        ac = cache['allele_counts']
        oc = cache['other_counts']
        meta = cache['metadata']
        cache.close()
    else:
        print "Specify at least one data source."
        return 1

    wt = ac[meta['Pos'] == 0]
    vswt = np.log2(ac / np.sum(wt))

    for k in expdefs:
        cols = [t[0] for t in expdefs[k]]
        pool = Pool(processes=args.numproc)
        result = pool.map(linregress_wrapper, vswt[cols].values)
        fitness = pd.DataFrame(data=result, index=vswt.index, columns=['slope', 'intercept', 'r', 'p', 'stderr'])
        pool.close()

    if args.cache is not None:
        cache = pytables.HDFStore(args.cache)
        cache['fitness'] = fitness
        cache.close()

    draw_hist(fitness, os.path.join(args.out, 'fitmap_hist.png'))
    idx = meta.reset_index().set_index(['Pos', 'AA', 'Codon'])
    aa = list(set(meta['AA']) - {None, np.nan})
    cod = list(set(meta['Codon']) - {None, np.nan})
    pos = list(set(meta['Pos']) - {None, np.nan, 0})
    aamaps = compute_hmap(fitness['slope'], pos, aa, 'Pos', 'AA', idx, [np.nanmedian, np.nanstd])
    codmaps = compute_hmap(fitness['slope'], pos, cod, 'Pos', 'Codon', idx, [np.nanmedian, np.nanstd])
    draw_hmap(aamaps[0], aa, os.path.join(args.out, 'fitmap_aa.png'))
    draw_hmap(codmaps[0], cod, os.path.join(args.out, 'fitmap_cod.png'))
    return 0


def draw_hist(fitness, fname):
    # cnt, bns = np.histogram(fitness['slope'], bins=20, normed=True)
    fig = plt.figure()
    # plt.bar(bns, cnt)
    plt.hist(fitness['slope'][np.isfinite(fitness['slope'])], bins=50, normed=True)
    fig.set_facecolor('white')
    ax = plt.gca()
    ax.set_xlabel('Fitness Score')
    ax.set_ylabel('Frequency')
    plt.savefig(fname)
    # plt.show()


def draw_hmap(hmap, yvals, fname):
    fig = plt.figure()
    plt.pcolor(hmap, cmap='RdBu')
    plt.xlim(0, hmap.shape[1])
    plt.ylim(0, hmap.shape[0])
    ax = plt.gca()
    fig.set_facecolor('white')
    ax.set_yticks([x+0.6 for x in xrange(0, hmap.shape[0])])
    ax.set_yticklabels(yvals)
    ax.set_ylabel('Residue')
    ax.set_xlabel('Ub Sequence Position')
    cb = plt.colorbar()
    cb.set_label('Relative Fitness')
    plt.savefig(fname, bbox_inches='tight')


def compute_hmap(fitness, xvals, yvals, xcol, ycol, idx, avgfuns=[np.nanmean]):
    maps = [np.zeros((len(yvals), len(xvals))) for f in avgfuns]
    for i in xrange(0, len(yvals)):
        for p in xrange(0, len(xvals)):
            for m, f in zip(maps, avgfuns):
                m[i, p] = f(fitness.loc[idx.xs([xvals[p], yvals[i]], level=[xcol, ycol])['Seq']])
    return maps


def linregress_wrapper(y):
    slope, intercept, rval, pval, err = stats.linregress(x=np.array([0, 1.87, 3.82]), y=y)
    return slope, intercept, rval, pval, err


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


def process_fastq_files(args, indices=None):
    metadata = parse_library(args.dbfile)
    if len(args.files) > 1:
        allele_temp, other_temp = process_split_fastq_files(args, metadata)
    else:
        allele_temp, other_temp = process_fastq_file(args, indices, metadata.index)
    allele_counts = pd.DataFrame.from_dict(allele_temp, orient="columns")
    other_counts = pd.DataFrame.from_dict(other_temp, orient="columns")
    allele_counts.index.rename('Seq', inplace=True)
    other_counts.index.rename('Seq', inplace=True)
    return allele_counts, other_counts, metadata


def process_fastq_file(args, indices, barcodes):
    alleles = {i: {k: 0 for k in barcodes} for i in indices}
    others = {i: defaultdict(int) for i in indices}
    parse_fastq_file(args.files[0], alleles, others)
    return alleles, others


def parse_fastq_file(fastq, counts, others):
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


def process_split_fastq_files(args, metadata):
    pool = Pool(processes=args.numproc)
    result = pool.map(functools.partial(process_split_fastq, prepop=metadata.index), args.files)
    pool.close()
    allele_temp = {file2col(f): r[0] for f, r in zip(args.files, result)}
    other_temp = {file2col(f): r[1] for f, r in zip(args.files, result)}
    return allele_temp, other_temp


def process_split_fastq(fastq, prepop):
    alleles = {k: 0 for k in prepop}
    others = defaultdict(int)
    parse_split_fastq(fastq, alleles, others)
    return alleles, dict(others)


def parse_split_fastq(fastq, counts, others):
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


def file2col(x):
    return os.path.basename(x).replace('.fastq', '').replace('.gz', '')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Map read sequences to known barcodes.")
    parser.add_argument("-p", "--processes", action="store", type=int, default=1,
                        dest="numproc", help="Number of parallel processes.")
    parser.add_argument("-d", "--database", action="store", type=str,
                        default=None, dest="dbfile", help="Barcode definitions (allele dictionary).")
    parser.add_argument("-e", "--expdefs", action="store", type=str, default=None,
                        dest="expdefs", help="Experiment definitions.")
    parser.add_argument("-c", "--cache", action="store", dest="cache", default=None,
                        help="Cache file path.")
    parser.add_argument("-o", "--outdir", action="store", type=str,
                        dest="out", help="Output figure directory.")
    parser.add_argument("files", nargs='*')
    sys.exit(main(parser.parse_args()))
