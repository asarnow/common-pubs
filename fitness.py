#!/usr/bin/env python2.7
from __future__ import print_function
import sys
import os.path
import itertools
from multiprocessing import Pool
import cPickle as pik
from Bio import Seq
from Bio.Data import CodonTable
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import seaborn
seaborn.set()
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def main(args):
    if args.bardefs is None:
        print("Must specify allele dictionary!")
        return 1

    meta = parse_library(args.bardefs)

    try:  # Read Tamas-brand counts data and pivot.
        grouped_counts = pd.read_hdf(args.input[0], key='grouped_data')
        grouped_counts.reset_index(inplace=True)
        counts = grouped_counts.pivot(index='barcodes', columns='index', values='counts')
        counts.index.rename('Seq', inplace=True)

        # oc = counts.drop(meta.index, errors='ignore')  # Unexpected barcodes.
        ac = counts.loc[meta.index]
        cdata = meta.merge(ac, left_index=True, right_index=True)

        if args.expdefs is None:
            print("Must specify experiment definitions!")
            return 1
        with open(args.expdefs, 'rb') as f:
            expdefs = pik.load(f)

        if args.groupby is not None:
            if args.groupby.lower().startswith('c'):  # E.g. 'residue'.
                cntsum = cdata.groupby(['Pos', 'Codon', 'AA']).sum()
            elif args.groupby.lower().startswith('r'):  # E.g. 'residue'.
                cntsum = cdata.groupby(['Pos', 'AA']).sum()
            else:  # E.g. 'barcode'.
                cntsum = cdata

            vswt = cnt2pop(cntsum)  # Convert counts to relative population measure.
            # Compute fitness scores for all experiments.
            fitdat = compute_fitness(vswt, cntsum, expdefs, args.numproc)
        else:
            cntsums = {'AA': cdata.groupby(['Pos', 'AA']).sum(),
                       'Codon': cdata.groupby(['Pos', 'Codon', 'AA']).sum(),
                       'Barcode': cdata.set_index(['Pos', 'Codon', 'AA'])}
            vswt = {k: cnt2pop(cntsums[k]) for k in cntsums}
            fitdat = {k: compute_fitness(vswt[k], cntsums[k], expdefs, args.numproc) for k in cntsums}

        if not os.path.isdir(args.output[0]):  # No figure directory, save fitness scores.
            if type(fitdat) is dict:
                for k in fitdat:
                    fitdat[k].to_hdf(args.output[0], key=k, mode='w', complevel=9)
            else:
                fitdat.to_hdf(args.output[0], key="fitness", mode='w', complevel=9)
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
    # mask = fitness[args.exp]['stderr'] < 0.05
    aamap = compute_hmap(fitness[args.exp][args.score], pos, aa, 'Pos', 'AA', idx, np.nanmedian)
    codmap = compute_hmap(fitness[args.exp][args.score], pos, cod, 'Pos', 'Codon', idx, np.nanmedian)
    # Plot heat maps and write image files.
    draw_hmap(aamap, aa, os.path.join(args.output[0], args.prefix + 'heatmap_aa.png'))
    draw_hmap(codmap, cod, os.path.join(args.output[0], args.prefix + 'heatmap_codon.png'))
    # TODO process all exps with replicates.
    # TODO Compute corr coef, difference heat maps.
    return 0


def draw_hmap_old(hmap, yvals, fname=None):
    """
    Plot a matrix as a heat map and write an image file.
    :param hmap: Heat map matrix.
    :param yvals: Heat map Y labels (e.g. amino acid names).
    :param fname: Destination image file.
    """
    if np.nanmax(hmap) > abs(np.nanmin(hmap)):
        vmax = np.nanmax(hmap)
        vmin = -np.nanmax(hmap)
    else:
        vmax = abs(np.nanmin(hmap))
        vmin = np.nanmin(hmap)
    fig = plt.figure()
    plt.figure(figsize=(20,10))
    plt.imshow(hmap, cmap='RdBu', interpolation = 'nearest',aspect='auto',vmin = vmin ,vmax = vmax )
    plt.xlim(0, hmap.shape[1])
    plt.ylim(0, hmap.shape[0])
    ax = plt.gca()
    fig.set_facecolor('white')
    ax.set_xlim((-0.5, hmap.shape[1] -0.5))
    ax.set_ylim((-0.5, hmap.shape[0] -0.5))
    ax.set_yticks([x for x in xrange(0, hmap.shape[0])])
    ax.set_yticklabels(yvals)
    ax.set_xticks(range(0,76,5))
    ax.set_xticklabels(range(2,76,5)+['STOP'])
    ax.set_ylabel('Residue')
    ax.set_xlabel('Ub Sequence Position')
    cb = plt.colorbar()
    cb.set_clim(vmin=vmin, vmax=vmax)
    cb.set_label('Relative Fitness')
    if fname is not None:
        plt.savefig(fname, bbox_inches='tight')
    return fig


def draw_hmap(hmap, yvals, fname=None):
    """
    Plot a matrix as a heat map and write an image file.
    :param hmap: Heat map matrix.
    :param yvals: Heat map Y labels (e.g. amino acid names).
    :param fname: Destination image file.
    """
    if np.nanmax(hmap) > abs(np.nanmin(hmap)):
        vmax = np.nanmax(hmap)
        vmin = -np.nanmax(hmap)
    else:
        vmax = abs(np.nanmin(hmap))
        vmin = np.nanmin(hmap)
    if len(yvals[0]) == 1:
        fig = plt.figure(figsize=(20, 10))
    else:
        fig = plt.figure(figsize=(20, 20))
    cm = plt.get_cmap('RdBu')
    cm.set_bad(color='k', alpha=1.0)
    plt.imshow(hmap, cmap=cm, interpolation='nearest', aspect='auto', vmin=vmin, vmax=vmax)
    plt.xlim(0, hmap.shape[1])
    plt.ylim(0, hmap.shape[0])
    ax = plt.gca()
    fig.set_facecolor('white')
    plt.grid('off')
    ax.set_xlim((-0.5, hmap.shape[1] - 0.5))
    ax.set_ylim((-0.5, hmap.shape[0] - 0.5))
    ax.set_yticks([x for x in xrange(0, hmap.shape[0])])
    if len(yvals[0]) == 1:
        ax.set_yticklabels(yvals, size=18, weight='bold')
    else:
        ax.set_yticklabels([gclut[y] + ' ' + y for y in yvals], size=18, weight='bold')
    ax.set_xticks(range(0, 76, 5))
    ax.set_xticklabels(range(2, 76, 5) + ['STOP'], size=24)
    if len(yvals[0]) == 1:
        ax.set_ylabel('Residue', size=32)
    else:
        ax.set_ylabel('Codon', size=32)
    ax.set_xlabel('Ub Sequence Position', size=32)
    cb = plt.colorbar()
    cb.set_clim(vmin=-0.7, vmax=0.7)
    cb.set_label('Relative Fitness', size=32)
    for p in xrange(0, hmap.shape[1]):
        if len(yvals[0]) == 1:
            ax.add_patch(
                Rectangle((p - 0.5, yvals.index(wt_residues()[p + 1]) - 0.5), 1, 1, fill=False, edgecolor='lime', lw=3))
        else:
            ax.add_patch(
                Rectangle((p - 0.5, yvals.index(wt_codons()[p + 1]) - 0.5), 1, 1, fill=False, edgecolor='lime', lw=3))
    if fname is not None:
        plt.savefig(fname, bbox_inches='tight')
        plt.close(fig)


def draw_pos(pos, vals, fname=None, ylabel="Fitness"):
    fig = plt.figure(figsize=(20, 10))
    fig.set_facecolor('white')
    plt.bar(range(0, len(pos)), vals)
    ax = plt.gca()
    ax.set_xticks(np.array(range(0, len(pos), 5)) + 0.25)
    ax.set_xticklabels(pos[0:-2:5] + ["STOP"])
    ax.set_xlabel('Ub Sequence Position')
    ax.set_ylabel(ylabel)
    if fname is not None:
        plt.savefig(fname, bbox_inches='tight')
        plt.close(fig)


def compute_hmap(fitness, xvals, yvals, xcol, ycol, avgfun=np.nanmean):
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
            # aamap[i, p] = avgfun(fitness.loc[idx.xs([xvals[p], yvals[i]], level=[xcol, ycol])['Seq']])
            aamap[i, p] = avgfun(fitness.xs([xvals[p], yvals[i]], level=[xcol, ycol]))
    return aamap


def cnt2pop(counts):
    """
    Convert from raw barcode sequence counts to relative population measure.
    :param counts: DataFrame containing indexed counts data.
    :return: Relative population measure for non-wild-type barcodes.
    """
    cdata = counts.copy(deep=True)
    cdata[np.isnan(cdata)] = 0.5
    wt = cdata.loc[0]  # All true wild-type Ub nucleotide sequences.
    vswt = np.log2(cdata / np.sum(wt))
    return vswt


def compute_weights(cnts, pos, cod, idx):
    w = cnts.copy(deep=True)
    for p, c in itertools.product(pos, cod):
        vals = cnts.loc[idx.xs([p, c], level=['Pos', 'Codon'])['Seq']]
        w.loc[idx.xs([p, c], level=['Pos', 'Codon'])['Seq']] = vals / np.sum(vals, axis=0)
    return w


def compute_weighting(vswt, counts, pos, cod, idx):
    prod = list(itertools.product(pos, cod))
    new = pd.DataFrame(data=None, index=[[p[0] for p in prod], [p[1] for p in prod]], columns=vswt.columns)
    for p, c in prod:
        ind = idx.xs([p, c], level=['Pos', 'Codon'])['Seq']
        w = counts.loc[ind] / np.sum(counts.loc[ind], axis=0)
        new.loc[p, c] = np.sum(vswt.loc[ind] * w, axis=0)
    return new


def compute_fitness(data, counts, expdefs, numproc=1):
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
            poolargs = [data[cols].values, counts[cols].values, xvals]
        else:
            poolargs = data[cols].values
        pool = Pool(processes=numproc)
        result = pool.map(linregress_wrapper, poolargs)
        fitness[k] = pd.DataFrame(data=result, index=data.index, columns=['slope', 'intercept', 'r', 'p', 'stderr'])
        pool.close()
    panel = pd.Panel(fitness)
    return panel


def linregress_wrapper(y, counts=None, x=np.array([0, 1.87, 3.82])):
    """
    Wrapper function allowing stats.linregress to be used with multiprocessing.Pool.
    :param y: Vector of relative population measure (dependent variable).
    :param x: Vector of time points (independent variable).
    :return: Tuple containing fields of LinearRegressionResult.
    """
    if counts is not None:
        if np.isnan(counts[0]):
            return (np.nan,) * 5
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
            allele[s] = [pos, "", ""]
    metadata = pd.DataFrame(data=([k] + v for k, v in allele.items()),
                            columns=['Seq', 'Pos', 'Codon', 'AA']).set_index('Seq')
    # metadata.loc[metadata[metadata.Pos == 0].index] = [0, "", ""]
    return metadata


def wt_codons():
    seq = 'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTTGA'
    n = 3
    codons = [seq[i:i + n] for i in xrange(0, len(seq), n)]
    return codons


def wt_residues():
    seq = Seq.Seq(
        'ATGCAGATTTTCGTCAAGACTTTGACCGGTAAAACCATAACATTGGAAGTTGAATCTTCCGATACCATCGACAACGTTAAGTCGAAAATTCAAGACAAGGAAGGTATCCCTCCAGATCAACAAAGATTGATCTTTGCCGGTAAGCAGCTAGAAGACGGTAGAACGCTGTCTGATTACAACATTCAGAAGGAGTCCACCTTACATCTTGTGCTAAGGCTAAGAGGTGGTTGA',
        Seq.Alphabet.DNAAlphabet())
    aa = seq.translate()
    return str(aa)


def hydrophobicity(aa=None):
    hphob = {'F': 100, 'I': 99, 'W': 97, 'L': 97, 'V': 76, 'M': 74, 'Y': 63, 'C': 49, 'A': 41, 'T': 13, 'H': 8,
             'G': 0, 'S': -5, 'Q': -10, 'R': -14, 'K': -23, 'N': -28, 'E': -31, 'P': -46, 'D': -55}
    if aa is None:
        return hphob
    return hphob[aa]


def _gclut():
    gencode = CodonTable.unambiguous_dna_by_name["Standard"]
    gclut = gencode.forward_table
    for c in gencode.stop_codons:
        gclut[c] = '*'
    return gclut


gclut = _gclut()

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
    parser.add_argument("-g", "--groupby", action="store", type=str, dest="groupby",
                        default="None", help="Group counts by barcode, codon or residue.")
    parser.add_argument("input", nargs=1, metavar="COUNTS")
    parser.add_argument("output", nargs=1, metavar="FITNESS")
    sys.exit(main(parser.parse_args()))
