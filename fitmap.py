#!/usr/bin/env python2.7
import sys
import os
import gzip
from collections import defaultdict
import cPickle as pik
from Bio import Seq
from Bio import SeqIO
import numpy as np
import scipy as sp
from scipy import stats
import pandas as pd
from pandas.io import pytables


def main(args):
    if args.parse:
        process_fastq(args)
    return 0


def process_fastq(args):
    store = pytables.HDFStore(args.out)
    metadata = parse_library(args.dbfile)
    store['metadata'] = metadata

    cols = [fn.replace('.fastq', '').replace('.gz', '') for fn in args.files]

    allele_counts = pd.DataFrame(data=None, columns=cols, index=metadata.index)
    other_counts = pd.DataFrame(data=None, columns=cols)

    for fastq in args.files:
        col = fastq.replace('.fastq', '').replace('.gz', '')
        others = defaultdict(int)
        parse_file(fastq, allele_counts[col], others)
        others = pd.Series(others)
        for i in others.index:
            other_counts.loc[i] = np.nan
            other_counts[col][i] = others[i]
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
    parser.add_argument("-p", "--parse", action="store_true", dest="parse",
                        help="Process FASTQ files")
    parser.add_argument("-d", "--database", action="store", type=str,
                        default=None, dest="dbfile", help="Barcode definitions or HDF5 file.")
    parser.add_argument("-o", "--output", action="store", type=str,
                        dest="out", help="Output Pickle with barcode definitions and frequencies.")
    parser.add_argument("files", nargs='+')
    sys.exit(main(parser.parse_args()))
