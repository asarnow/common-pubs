#!/usr/bin/env python2.7
import sys
import cPickle as pik
from Bio import Seq


def main(args):
    if args.dbfile is None:
        args.dbfile = 'data/allele_dic_with_WT.pkl'
    library = pik.load(open(args.dbfile, 'rb'))
    allele = {}
    for k in library:
        s = Seq.Seq(k, Seq.Alphabet.SingleLetterAlphabet()).reverse_complement()
        pos = int(library[k][0])
        if pos != 0:
            cod = Seq.Seq(library[k][1], Seq.Alphabet.SingleLetterAlphabet())
            aa = cod.transcribe().translate()
            allele[s] = [pos, cod, aa, 0]
        else:
            allele[s] = [pos, None, None, 0]

    for fname in args.files:
        counts = pik.load(open(fname, 'rb'))
        for k in allele:
            if k in counts:
                allele[k][3] = counts[k]
    pik.dump(allele, open(args.out, 'wb'))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Map read sequences to known barcodes.")
    parser.add_argument("-d", "--database", action="store", type=str,
                        default=None, dest="dbfile", help="Pickle with barcode definitions.")
    parser.add_argument("-o", "--output", action="store", type=str,
                        dest="out", help="Output Pickle with barcode definitions and frequencies.")
    parser.add_argument("files", nargs='+')
    sys.exit(main(parser.parse_args()))
