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
from scipy import stats
import pandas as pd
from pandas.io import pytables


def main(args):
    if ".fastq" in args.files[0]:
        process_fastq_files(args)
    return 0


def process_fastq_files(args):
    store = pytables.HDFStore(args.out)
    metadata = parse_library(args.dbfile)
    store['metadata'] = metadata
    # pool = Pool(processes=args.numproc)
    result = map(functools.partial(process_fastq, prepop=metadata.index), args.files)
    allele_temp = {file2col(f): r[0] for f, r in zip(args.files, result)}
    other_temp = {file2col(f): r[1] for f, r in zip(args.files, result)}
    allele_counts = pd.DataFrame.from_dict(allele_temp, orient="columns")
    other_counts = pd.DataFrame.from_dict(other_temp, orient="columns")
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
    parser.add_argument("-o", "--output", action="store", type=str,
                        dest="out", help="Output Pickle with barcode definitions and frequencies.")
    parser.add_argument("files", nargs='+')
    sys.exit(main(parser.parse_args()))
