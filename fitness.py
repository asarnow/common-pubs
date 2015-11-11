#!/usr/bin/env python2.7
import sys
from multiprocessing import Pool
import cPickle as pik
from Bio import Seq
import numpy as np
from scipy import stats
import pandas as pd


def main(args):
    if args.bardefs is None:
        print "Must specify allele dictionary!"
        return 1
    if args.expdefs is None:
        print "Must specify experiment definitions!"
        return 1

    meta = parse_library(args.bardefs)
    # idx = meta.reset_index().set_index(['Pos', 'AA', 'Codon'])

    # Read Tamas-brand counts data and pivot.
    grouped_counts = pd.read_hdf(args.input, key='grouped_data')
    grouped_counts.reset_index(inplace=True)
    counts = grouped_counts.pivot(index='barcodes', columns='index', values='counts')
    # ac = counts.loc[meta.index]
    # oc = counts.drop(meta.index, errors='ignore')  # Unexpected barcodes.
    vswt = cnt2fit(counts, meta)

    with open(args.expdefs, 'rb') as f:
        expdefs = pik.load(f)

    fitness = compute_fitness(vswt, expdefs, args.numproc)
    fitness.to_hdf(args.ouput, key="fitness", mode='w', complevel=9)
    return 0


def cnt2fit(counts, meta):
    ac = counts.loc[meta.index]
    ac[ac == 0] = 0.5
    ac[np.isnan(ac)] = 0.5
    wt = ac[meta['Pos'] == 0]
    vswt = np.log2(ac / np.sum(wt))
    return vswt


def compute_fitness(data, expdefs, numproc=1):
    fitness = {k: [] for k in expdefs}
    for k in expdefs:
        cols = [t[0] for t in expdefs[k]]
        if len(expdefs[k][0]) > 2:
            xvals = [t[1] for t in expdefs[k]]
            poolargs = [data[cols].values, xvals]
        else:
            poolargs = data[cols].values
        pool = Pool(processes=numproc)
        result = pool.map(linregress_wrapper, poolargs)
        fitness[k] = pd.DataFrame(data=result, index=data.index, columns=['slope', 'intercept', 'r', 'p', 'stderr'])
        pool.close()
    panel = pd.Panel(fitness)
    return panel


def linregress_wrapper(y, x=np.array([0, 1.87, 3.82])):
    slope, intercept, rval, pval, err = stats.linregress(x=x, y=y)
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


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Map read sequences to known barcodes.")
    parser.add_argument("-p", "--processes", action="store", type=int, default=1,
                        dest="numproc", help="Number of parallel processes.")
    parser.add_argument("-b", "--bardefs", action="store", type=str, default=None,
                        dest="bardefs", help="Barcode definitions (allele dictionary).")
    parser.add_argument("-d", "--expdefs", action="store", type=str, default=None,
                        dest="expdefs", help="Experiment definitions.")
    parser.add_argument("-e", "--experiment", action="store", type=str, default=None,
                        dest="exp", metavar="EXPNAME", help="Experiment label.")
    parser.add_argument("input", nargs=1, metavar="COUNTS")
    parser.add_argument("output", nargs=1, metavar="FITNESS")
    sys.exit(main(parser.parse_args()))
