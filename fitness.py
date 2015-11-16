#!/usr/bin/env python2.7
from __future__ import print_function
import sys
import os.path
from multiprocessing import Pool
import cPickle as pik
from Bio import Seq
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt


def main(args):
    if args.bardefs is None:
        print("Must specify allele dictionary!")
        return 1

    meta = parse_library(args.bardefs)

    try:  # Read Tamas-brand counts data and pivot.
        grouped_counts = pd.read_hdf(args.input[0], key='grouped_data')
        grouped_counts.reset_index(inplace=True)
        counts = grouped_counts.pivot(index='barcodes', columns='index', values='counts')
        # ac = counts.loc[meta.index]
        # oc = counts.drop(meta.index, errors='ignore')  # Unexpected barcodes.
        vswt = cnt2pop(counts, meta)  # Convert counts to relative population measure.

        if args.expdefs is None:
            print("Must specify experiment definitions!")
            return 1
        with open(args.expdefs, 'rb') as f:
            expdefs = pik.load(f)

        # Compute fitness scores for all experiments.
        fitness = compute_fitness(vswt, expdefs, args.numproc)

        if not os.path.isdir(args.output[0]):  # No figure directory, save fitness scores.
            fitness.to_hdf(args.output[0], key="fitness", mode='w', complevel=9)
            return 0
    except KeyError:
        try:  # Read previously computed fitness data.
            fitness = pd.read_hdf(args.input[0], key='fitness')
        except KeyError:
            print("No known HDF5 key found.")
            return 1

    if not os.path.isdir(args.output[0]):  # Make sure we have a figure output directory.
        print("Output path is not a directory.")
        return 1

    # Variables needed for plotting.
    idx = meta.reset_index().set_index(['Pos', 'AA', 'Codon'])  # Metadata index.
    aa = list(set(meta['AA']) - {None, np.nan})  # Unique amino acids.
    cod = list(set(meta['Codon']) - {None, np.nan})  # Unique codons.
    pos = list(set(meta['Pos']) - {None, np.nan, 0})  # Unique Ub positions.
    # Computing the heat maps using the metadata index.
    # TODO Loop over experiments.
    # TODO Compute std heat maps.
    # TODO Sorted list of positions by: mean, max, std, min, isnan of mean/med/etc. maps.
    # TODO Where isnan set in all experiments.
    # TODO Subtract control values and repeat analysis.
    # TODO Exclude non-surface residues and repeat.
    # TODO Group residues by type and repeat.
    aamap = compute_hmap(fitness[args.exp][args.score], pos, aa, 'Pos', 'AA', idx, np.nanmedian)
    codmap = compute_hmap(fitness[args.exp][args.score], pos, cod, 'Pos', 'Codon', idx, np.nanmedian)
    # Plot heat maps and write image files.
    draw_hmap(aamap, aa, os.path.join(args.output[0], args.prefix + 'heatmap_aa.png'))
    draw_hmap(codmap, cod, os.path.join(args.output[0], args.prefix + 'heatmap_codon.png'))
    # TODO process all exps with replicates.
    # TODO Compute corr coef, difference heat maps.
    return 0


def draw_hmap(hmap, yvals, fname):
    """
    Plot a matrix as a heat map and write an image file.
    :param hmap: Heat map matrix.
    :param yvals: Heat map Y labels (e.g. amino acid names).
    :param fname: Destination image file.
    """
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
    """
    Use the metadata index to efficiently populate fitness heat map.
    :param fitness: DataFrame containing fitness scores (per barcode).
    :param xvals: List of heat map X values / labels (e.g. Ub positions).
    :param yvals: List of heat map Y values / labels (e.g. amino acids).
    :param xcol: Name of X value column in fitness DataFrame.
    :param ycol: Name of Y value column in fitness DataFrame.
    :param idx: DataFrame containing indexed metadata.
    :param avgfun: Aggregation function for grouped barcode data (e.g. with same mutation).
    :return: Heat map matrix.
    """
    aamap = np.zeros((len(yvals), len(xvals)))
    for i in xrange(0, len(yvals)):
        for p in xrange(0, len(xvals)):
            aamap[i, p] = avgfun(fitness.loc[idx.xs([xvals[p], yvals[i]], level=[xcol, ycol])['Seq']])
    return aamap


def cnt2pop(counts, meta):
    """
    Convert from raw barcode sequence counts to relative population measure.
    :param counts: DataFrame containing counts data.
    :param meta: DataFrame containing barcode metadata.
    :return: Relative population measure for non-wild-type barcodes.
    """
    ac = counts.loc[meta.index]
    ac[ac == 0] = 0.5
    ac[np.isnan(ac)] = 0.5
    wt = ac[meta['Pos'] == 0]  # All true wild-type Ub nucleotide sequences.
    vswt = np.log2(ac / np.sum(wt))
    return vswt


def compute_fitness(data, expdefs, numproc=1):
    """
    Compute fitness scores with barcode-level parallelism.
    Currently uses linear regression of relative population values.
    :param data: DataFrame containing relative population data.
    :param expdefs: Dictionary defining experiment names, TrueSeq indices and time points.
    :param numproc: Number of parallel processes in process pool.
    :return: Panel containing fitness computation results for each experiment.
    """
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
    """
    Wrapper function allowing stats.linregress to be used with multiprocessing.Pool.
    :param y: Vector of relative population measure (dependent variable).
    :param x: Vector of time points (independent variable).
    :return: Tuple containing fields of LinearRegressionResult.
    """
    slope, intercept, rval, pval, err = stats.linregress(x=x, y=y)
    return slope, intercept, rval, pval, err


def parse_library(dbfile):
    """
    Read barcode allele dictionary and convert to DataFrame containing barcode metadata.
    :param dbfile: Pickle file holding barcode allele dictionary.
    :return: DataFrame containing barcode metadata.
    """
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
    parser.add_argument("-s", "--score", action="store", dest="score", type=str,
                        default='slope', metavar="SCORE", help="Fitness score label (e.g. 'slope').")
    parser.add_argument("--prefix", action="store", type=str, default="",
                        metavar="PREFIX", help="Prefix for output figure files.")
    parser.add_argument("input", nargs=1, metavar="COUNTS")
    parser.add_argument("output", nargs=1, metavar="FITNESS")
    sys.exit(main(parser.parse_args()))
