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
    if args.counts and args.scores is not None:
        print "Scores and counts are mutually exclusive!"
        return 1
    if args.counts:
        print "Not implemented, my bad."
        return 1
    if len(args.files) < 1 or len(args.files) > 2:
        print "Usage: fitmap.py <input> <output_dir>."
        return 1
    if args.dbfile is None:
        print "Must specify allele dictionary!"
        return 1

    meta = parse_library(args.dbfile)
    idx = meta.reset_index().set_index(['Pos', 'AA', 'Codon'])

    data = pik.load(open(args.files[0], 'rb'))
    if type(data) is dict:
        df = pd.DataFrame.from_dict(data, orient="columns")
    elif type(data) is pd.DataFrame:
        df = data
    else:
        print "Input data must be dict of dict or DataFrame!"
        return 1

    df.index.rename('Seq', inplace=True)
    aa = list(set(meta['AA']) - {None, np.nan})
    cod = list(set(meta['Codon']) - {None, np.nan})
    pos = list(set(meta['Pos']) - {None, np.nan, 0})
    aamap = compute_hmap(df[args.score], pos, aa, 'Pos', 'AA', idx, np.nanmedian)
    codmap = compute_hmap(df[args.score], pos, cod, 'Pos', 'Codon', idx, np.nanmedian)
    draw_hmap(aamap, meta, os.path.join(args.files[1], 'aamap.png'))
    draw_hmap(codmap, meta, os.path.join(args.files[1], 'codonmap.png'))
    return 0


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


def compute_hmap(fitness, xvals, yvals, xcol, ycol, idx, avgfun=np.nanmean):
    aamap = np.zeros((len(yvals), len(xvals)))
    for i in xrange(0, len(yvals)):
        for p in xrange(0, len(xvals)):
            aamap[i, p] = avgfun(fitness.loc[idx.xs([xvals[p], yvals[i]], level=[xcol, ycol])['Seq']])
    return aamap


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


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Map read sequences to known barcodes.")
    parser.add_argument("-p", "--processes", action="store", type=int, default=1,
                        dest="numproc", help="Number of parallel processes.")
    parser.add_argument("-d", "--database", action="store", type=str,
                        default=None, dest="dbfile", help="Barcode definitions (allele dictionary).")
    parser.add_argument("-c", "--counts", action="store_true", dest="counts", default=False,
                        help="Input data contains sequence counts.")
    parser.add_argument("-s", "--score", action="store", dest="score", type=str, default='slope',
                        metavar="SCORE", help="Fitness score label (e.g. 'slope').")
    parser.add_argument("files", nargs="+")
    sys.exit(main(parser.parse_args()))
